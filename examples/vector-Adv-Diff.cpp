//                                MFEM Example 9
//
// Compile with: make ex9
//
// Sample runs:
//    ex9 -m ../data/periodic-segment.mesh -p 0 -r 2 -dt 0.005
//    ex9 -m ../data/periodic-square.mesh -p 0 -r 2 -dt 0.01 -tf 10
//    ex9 -m ../data/periodic-hexagon.mesh -p 0 -r 2 -dt 0.01 -tf 10
//    ex9 -m ../data/periodic-square.mesh -p 1 -r 2 -dt 0.005 -tf 9
//    ex9 -m ../data/periodic-hexagon.mesh -p 1 -r 2 -dt 0.005 -tf 9
//    ex9 -m ../data/amr-quad.mesh -p 1 -r 2 -dt 0.002 -tf 9
//    ex9 -m ../data/amr-quad.mesh -p 1 -r 2 -dt 0.02 -s 23 -tf 9
//    ex9 -m ../data/star-q3.mesh -p 1 -r 2 -dt 0.005 -tf 9
//    ex9 -m ../data/star-mixed.mesh -p 1 -r 2 -dt 0.005 -tf 9
//    ex9 -m ../data/disc-nurbs.mesh -p 1 -r 3 -dt 0.005 -tf 9
//    ex9 -m ../data/disc-nurbs.mesh -p 2 -r 3 -dt 0.005 -tf 9
//    ex9 -m ../data/periodic-square.mesh -p 3 -r 4 -dt 0.0025 -tf 9 -vs 20
//    ex9 -m ../data/periodic-cube.mesh -p 0 -r 2 -o 2 -dt 0.02 -tf 8
//    ex9 -m ../data/periodic-square.msh -p 0 -r 2 -dt 0.005 -tf 2
//    ex9 -m ../data/periodic-cube.msh -p 0 -r 1 -o 2 -tf 2
//
// Device sample runs:
//    ex9 -pa
//    ex9 -ea
//    ex9 -fa
//    ex9 -pa -m ../data/periodic-cube.mesh
//    ex9 -pa -m ../data/periodic-cube.mesh -d cuda
//    ex9 -ea -m ../data/periodic-cube.mesh -d cuda
//    ex9 -fa -m ../data/periodic-cube.mesh -d cuda
//    ex9 -pa -m ../data/amr-quad.mesh -p 1 -r 2 -dt 0.002 -tf 9 -d cuda
//
// Description:  This example code solves the time-dependent advection equation
//               du/dt + v.grad(u) = 0, where v is a given fluid velocity, and
//               u0(x)=u(0,x) is a given initial condition.
//
//               The example demonstrates the use of Discontinuous Galerkin (DG)
//               bilinear forms in MFEM (face integrators), the use of implicit
//               and explicit ODE time integrators, the definition of periodic
//               boundary conditions through periodic meshes, as well as the use
//               of GLVis for persistent visualization of a time-evolving
//               solution. The saving of time-dependent data files for external
//               visualization with VisIt (visit.llnl.gov) and ParaView
//               (paraview.org) is also illustrated.

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional>


using namespace std;
using namespace mfem;

// Choice for the problem setup. The fluid velocity, initial condition and
// inflow boundary condition are chosen based on this parameter.
int problem;

// Velocity coefficient
void velocity_function(const Vector &x, Vector &v);

// Initial condition
void u0_function(const Vector &x, Vector &u_func);

// Inflow boundary condition
real_t inflow_function(const Vector &x);

// Vector Forcing Function
void f_exact(const Vector &x, Vector &f_vec);

// Mesh bounding box
Vector bb_min, bb_max;
//

/// @brief Time dependent DG operator for hyperbolic conservation laws
class DGVectorEvolution : public TimeDependentOperator
{
private:
   BilinearForm &K;
   const Vector &b;
   //const int num_equations; // the number of equations
   const int dim;
   FiniteElementSpace &vfes; // vector finite element space
   // Element integration form. Should contain ComputeFlux
   //std::unique_ptr<HyperbolicFormIntegrator> formIntegrator;
   // Base Nonlinear Form
   //std::unique_ptr<NonlinearForm> nonlinearForm;
   // element-wise inverse mass matrix
   std::vector<DenseMatrix> invmass; // local scalar inverse mass

   // auxiliary variable used in Mult
   mutable Vector z;

   // Compute element-wise inverse mass matrix
   void ComputeInvMass();
   // Compute element-wise weak-divergence matrix
   //void ComputeWeakDivergence();

public:
   /**
    * @brief Construct a new DGHyperbolicConservationLaws object
    *
    * @param vfes_ vector finite element space. Only tested for DG [Pₚ]ⁿ
    * @param formIntegrator_ integrator (F(u,x), grad v)
    * @param preassembleWeakDivergence preassemble weak divergence for faster
    *                                  assembly
    */
   DGVectorEvolution(BilinearForm &K_, const Vector &b_, FiniteElementSpace &vfes_);
   /**
    * @brief Apply nonlinear form to obtain M⁻¹(DIVF + JUMP HAT(F))
    *
    * @param x current solution vector
    * @param y resulting dual vector to be used in an EXPLICIT solver
    */
   void Mult(const Vector &x, Vector &y) const override;

   void Update();

};

//////////////////////////////////////////////////////////////////
///        HYPERBOLIC CONSERVATION LAWS IMPLEMENTATION         ///
//////////////////////////////////////////////////////////////////

// Implementation of class DGHyperbolicConservationLaws
DGVectorEvolution::DGVectorEvolution(BilinearForm &K_, const Vector &b_,
   FiniteElementSpace &vfes_)
   : TimeDependentOperator(vfes_.GetTrueVSize()), K(K_), b(b_),
     dim(vfes_.GetMesh()->SpaceDimension()),
     vfes(vfes_),
     z(vfes_.GetTrueVSize())
{
   // Standard local assembly and inversion for energy mass matrices.
   ComputeInvMass();
#ifndef MFEM_USE_MPI
   //nonlinearForm.reset(new NonlinearForm(&vfes));
#else
   ParFiniteElementSpace *pvfes = dynamic_cast<ParFiniteElementSpace *>(&vfes);
   if (pvfes)
   {
      nonlinearForm.reset(new ParNonlinearForm(pvfes));
   }
   else
   {
      nonlinearForm.reset(new NonlinearForm(&vfes));
   }
#endif
//   if (preassembleWeakDivergence)
//   {
//      ComputeWeakDivergence();
//   }
//   else
//   {
//      nonlinearForm->AddDomainIntegrator(formIntegrator.get());
//   }
//   nonlinearForm->AddInteriorFaceIntegrator(formIntegrator.get());
//   nonlinearForm->UseExternalIntegrators();

}

void DGVectorEvolution::ComputeInvMass()
{
   InverseIntegrator inv_mass(new MassIntegrator());

   invmass.resize(vfes.GetNE());
   for (int i=0; i<vfes.GetNE(); i++)
   {
      int dof = vfes.GetFE(i)->GetDof();
      invmass[i].SetSize(dof);
      inv_mass.AssembleElementMatrix(*vfes.GetFE(i),
                                     *vfes.GetElementTransformation(i),
                                     invmass[i]);
   }
    
//    DenseMatrix Block_dense_mat;
//    // Print each block of inverse mass
//    for (int block = 0; block<vfes.GetNE(); block++)
//    {
//        Block_dense_mat = invmass[block];
//        Block_dense_mat.Print();
//        Block_dense_mat = 0.0;
//    }
    
}

void DGVectorEvolution::Mult(const Vector &x, Vector &y) const
{
    /*
     // y = M^{-1} (K x + b)
     K.Mult(x, z);
     //z.Neg(); // z = -z //New change
     z += b;
     M_solver.Mult(z, y);*/
    
    // New implementation
    K.Mult(x, z); // Multiply x by stiffness matrix K
    z += b; // Add linear form RHS
    
    // Apply block inverser mass M^(-1) * (z)
    // Apply block inverse mass
    Vector zval; // z_loc, dof*num_eq
    int num_equations = 2;
    DenseMatrix current_zmat; // view of element auxiliary result, dof x num_eq
    DenseMatrix current_ymat; // view of element result, dof x num_eq
    Array<int> vdofs;
    for (int i=0; i<vfes.GetNE(); i++)
    {
       int dof = vfes.GetFE(i)->GetDof();
       vfes.GetElementVDofs(i, vdofs);
       z.GetSubVector(vdofs, zval);
       current_zmat.UseExternalData(zval.GetData(), dof, num_equations);
       current_ymat.SetSize(dof, num_equations);
       mfem::Mult(invmass[i], current_zmat, current_ymat);
       y.SetSubVector(vdofs, current_ymat.GetData());
        //y.Print();
        //std::cout << "hellp";
    }
    
}

//void DGHyperbolicConservationLaws::ComputeWeakDivergence()
//{
//   TransposeIntegrator weak_div(new GradientIntegrator());
//   DenseMatrix weakdiv_bynodes;
//
//   weakdiv.resize(vfes.GetNE());
//   for (int i=0; i<vfes.GetNE(); i++)
//   {
//      int dof = vfes.GetFE(i)->GetDof();
//      weakdiv_bynodes.SetSize(dof, dof*dim);
//      weak_div.AssembleElementMatrix2(*vfes.GetFE(i), *vfes.GetFE(i),
//                                      *vfes.GetElementTransformation(i),
//                                      weakdiv_bynodes);
//      weakdiv[i].SetSize(dof, dof*dim);
//      // Reorder so that trial space is ByDim.
//      // This makes applying weak divergence to flux value simpler.
//      for (int j=0; j<dof; j++)
//      {
//         for (int d=0; d<dim; d++)
//         {
//            weakdiv[i].SetCol(j*dim + d, weakdiv_bynodes.GetColumn(d*dof + j));
//         }
//      }
//
//   }
//}


void DGVectorEvolution::Update()
{
    mfem_error("Inside nonlinear form update funtion");
}
/** A time-dependent operator for the right-hand side of the ODE. The DG weak
    form of du/dt = -v.grad(u) is M du/dt = K u + b, where M and K are the mass
    and advection matrices, and b describes the flow on the boundary. This can
    be written as a general ODE, du/dt = M^{-1} (K u + b), and this class is
    used to evaluate the right-hand side. */
class FE_Evolution : public TimeDependentOperator
{
private:
   BilinearForm &M, &K;
   const Vector &b;
//   Solver *M_prec;
//   CGSolver M_solver;
//   DG_Solver *dg_solver;

   mutable Vector z;
    
   // New code
   std::vector<DenseMatrix> invmass; // local inverse mass
   FiniteElementSpace &vfes; // vector finite element space

public:
   FE_Evolution(BilinearForm &M_, BilinearForm &K_, const Vector &b_, FiniteElementSpace &vfes_);

   void Mult(const Vector &x, Vector &y) const override;
   ~FE_Evolution() override;
};


int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   problem = 1;
   const char *mesh_file = "../data/cube2d_2.MFEM";
   //const char *mesh_file = "../data/segment-4.MFEM";
   int ref_levels = 0;
   int order = 1;
   real_t sigma = -1.0;
   //real_t kappa = -1.0;
   real_t kappa = (order+1)*(order+1);
   bool pa = false;
   bool ea = false;
   bool fa = false;
   const char *device_config = "cpu";
   int ode_solver_type = 4; //  4 RK4 method   21 implicit backward euler
   real_t t_final = 2.0;
   real_t dt = 0.000001;
   bool visualization = true;
   bool visit = false;
   bool paraview = false;
   bool binary = false;
   int vis_steps = 5;

   int precision = 8;
   cout.precision(precision);

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&problem, "-p", "--problem",
                  "Problem setup to use. See options in velocity_function().");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly.");
   args.AddOption(&order, "-o", "--order",
                  "Order (degree) of the finite elements.");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&ea, "-ea", "--element-assembly", "-no-ea",
                  "--no-element-assembly", "Enable Element Assembly.");
   args.AddOption(&fa, "-fa", "--full-assembly", "-no-fa",
                  "--no-full-assembly", "Enable Full Assembly.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  ODESolver::Types.c_str());
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Time step.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                  "--no-visit-datafiles",
                  "Save data files for VisIt (visit.llnl.gov) visualization.");
   args.AddOption(&paraview, "-paraview", "--paraview-datafiles", "-no-paraview",
                  "--no-paraview-datafiles",
                  "Save data files for ParaView (paraview.org) visualization.");
   args.AddOption(&binary, "-binary", "--binary-datafiles", "-ascii",
                  "--ascii-datafiles",
                  "Use binary (Sidre) or ascii format for VisIt data files.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   Device device(device_config);
   device.Print();

   // 2. Read the mesh from the given mesh file. We can handle geometrically
   //    periodic meshes in this code.
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   // 3. Define the ODE solver used for time integration. Several explicit
   //    Runge-Kutta methods are available.
   unique_ptr<ODESolver> ode_solver = ODESolver::Select(ode_solver_type);

   // 4. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement, where 'ref_levels' is a
   //    command-line parameter. If the mesh is of NURBS type, we convert it to
   //    a (piecewise-polynomial) high-order mesh.
   for (int lev = 0; lev < ref_levels; lev++)
   {
      mesh.UniformRefinement();
   }
   if (mesh.NURBSext)
   {
      mesh.SetCurvature(max(order, 1));
   }
   mesh.GetBoundingBox(bb_min, bb_max, max(order, 1));

   // 5. Define the discontinuous DG finite element space of the given
   //    polynomial order on the refined mesh.
   //DG_FECollection fec(order, dim, BasisType::GaussLobatto);
   DG_FECollection fec(order, dim);
   FiniteElementSpace fes(&mesh, &fec, dim, Ordering::byNODES);

   cout << "Number of unknowns: " << fes.GetVSize() << endl;

   // 6. Set up and assemble the bilinear and linear forms corresponding to the
   //    DG discretization. The DGTraceIntegrator involves integrals over mesh
   //    interior faces.
   VectorFunctionCoefficient velocity(dim, velocity_function);
   VectorFunctionCoefficient u_coeff(dim, u0_function); // Vector Diffusion
    
   FunctionCoefficient inflow(inflow_function);
   //VectorFunctionCoefficient u0(dim, u0_function);
    
   ConstantCoefficient one(1.0);
   ConstantCoefficient neg_one(-1.0);
   ConstantCoefficient diff_coeff(-1.0);
   ConstantCoefficient zero(0.0);
   
   // Forcing Function Coefficent 
   //FunctionCoefficient f(f_exact);
   VectorFunctionCoefficient f_coeff(dim, f_exact); // vector diffusion


   BilinearForm m(&fes);
   BilinearForm k(&fes);
   LinearForm b(&fes);

   if (pa)
   {
      m.SetAssemblyLevel(AssemblyLevel::PARTIAL);
      k.SetAssemblyLevel(AssemblyLevel::PARTIAL);
   }
   else if (ea)
   {
      m.SetAssemblyLevel(AssemblyLevel::ELEMENT);
      k.SetAssemblyLevel(AssemblyLevel::ELEMENT);
   }
   else if (fa)
   {
      m.SetAssemblyLevel(AssemblyLevel::FULL);
      k.SetAssemblyLevel(AssemblyLevel::FULL);
   }
    
   if (problem == 1) // Convection Only
   {
       std::cout << "Convection Only" << std::endl;
       m.AddDomainIntegrator(new MassIntegrator);
       constexpr real_t alpha = -1.0;
       k.AddDomainIntegrator(new ConvectionIntegrator(velocity, alpha));
       k.AddInteriorFaceIntegrator(
          new NonconservativeDGTraceIntegrator(velocity, alpha));
       k.AddBdrFaceIntegrator(
          new NonconservativeDGTraceIntegrator(velocity, alpha));

        
        // RHS Forcing Term
        b.AddDomainIntegrator(new VectorDomainLFIntegrator(f_coeff));
        //b.AddBdrFaceIntegrator(new BoundaryFlowIntegrator(inflow, velocity, -alpha)); // ?????

   }
   if (problem == 2) // Diffusion Only
   {
        std::cout << "Diffusion Only" << std::endl;
       
        // Time Dependant Term
        m.AddDomainIntegrator(new MassIntegrator);
                
        // Diffusion Integrators
        k.AddDomainIntegrator(new VectorDiffusionIntegrator(neg_one));
        k.AddInteriorFaceIntegrator(new VectorDGDiffusionIntegrator(neg_one, sigma, kappa, dim));
        k.AddBdrFaceIntegrator(new VectorDGDiffusionIntegrator(neg_one, sigma, kappa, dim));
       
        // RHS Forcing Term
        b.AddDomainIntegrator(new VectorDomainLFIntegrator(f_coeff));
        b.AddBdrFaceIntegrator(new VectorDGDirichletLFIntegrator(u_coeff, neg_one, sigma, kappa));
   }
    
   if (problem == 3) // Vector Convection Diffusion
   {
         std::cout << "DG Vector Diffusion + Convection" << std::endl;
       
       m.AddDomainIntegrator(new MassIntegrator);
       constexpr real_t alpha = -1.0;
       
       k.AddDomainIntegrator(new ConvectionIntegrator(velocity, alpha));
       k.AddInteriorFaceIntegrator(
          new NonconservativeDGTraceIntegrator(velocity, alpha));
       k.AddBdrFaceIntegrator(
          new NonconservativeDGTraceIntegrator(velocity, alpha));

        
        //Diffusion Part
        k.AddDomainIntegrator(new VectorDiffusionIntegrator(diff_coeff));
        k.AddInteriorFaceIntegrator(new VectorDGDiffusionIntegrator(diff_coeff, sigma, kappa, dim));
        k.AddBdrFaceIntegrator(new VectorDGDiffusionIntegrator(diff_coeff, sigma, kappa, dim));
       
       // RHS Forcing Term
       b.AddDomainIntegrator(new VectorDomainLFIntegrator(f_coeff));
       b.AddBdrFaceIntegrator(
           new VectorDGDirichletLFIntegrator(u_coeff, diff_coeff, sigma, kappa));

    }

    

   int skip_zeros = 0;
   m.Assemble(skip_zeros);
   m.Finalize(skip_zeros);

   k.Assemble(skip_zeros);
   k.Finalize(skip_zeros);
    
   std::ofstream linear_sys_m("DG-Vector-Diffusion-Sys_m.txt");
   m.PrintMatlab(linear_sys_m);
//   std::ofstream linear_sys_k("DG-Vector-Diffusion-Sys_k.txt");
//   k.PrintMatlab(linear_sys_k);
    
   b.Assemble();


   // 7. Define the initial conditions, save the corresponding grid function to
   //    a file and (optionally) save data in the VisIt format and initialize
   //    GLVis visualization.
   GridFunction u(&fes);
   u.ProjectCoefficient(u_coeff);

   {
      ofstream omesh("ex9.mesh");
      omesh.precision(precision);
      mesh.Print(omesh);
      ofstream osol("ex9-init.gf");
      osol.precision(precision);
      u.Save(osol);
   }

   // Create data collection for solution output: either VisItDataCollection for
   // ascii data files, or SidreDataCollection for binary data files.
   DataCollection *dc = NULL;
   if (visit)
   {
      if (binary)
      {
#ifdef MFEM_USE_SIDRE
         dc = new SidreDataCollection("Example9", &mesh);
#else
         MFEM_ABORT("Must build with MFEM_USE_SIDRE=YES for binary output.");
#endif
      }
      else
      {
         dc = new VisItDataCollection("Example9", &mesh);
         dc->SetPrecision(precision);
      }
      dc->RegisterField("solution", &u);
      dc->SetCycle(0);
      dc->SetTime(0.0);
      dc->Save();
   }

   ParaViewDataCollection *pd = NULL;
   if (paraview)
   {
      pd = new ParaViewDataCollection("Example9", &mesh);
      pd->SetPrefixPath("ParaView");
      pd->RegisterField("solution", &u);
      pd->SetLevelsOfDetail(order);
      pd->SetDataFormat(VTKFormat::BINARY);
      pd->SetHighOrderOutput(true);
      pd->SetCycle(0);
      pd->SetTime(0.0);
      pd->Save();
   }

   socketstream sout;
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      sout.open(vishost, visport);
      if (!sout)
      {
         cout << "Unable to connect to GLVis server at "
              << vishost << ':' << visport << endl;
         visualization = false;
         cout << "GLVis visualization disabled.\n";
      }
      else
      {
         sout.precision(precision);
         sout << "solution\n" << mesh << u;
         sout << "pause\n";
         sout << flush;
         cout << "GLVis visualization paused."
              << " Press space (in the GLVis window) to resume it.\n";
      }
   }

   // 8. Define the time-dependent evolution operator describing the ODE
   //    right-hand side, and perform time-integration (looping over the time
   //    iterations, ti, with a time-step dt).
   //FE_Evolution adv(m, k, b, fes);
    DGVectorEvolution adv(k, b, fes);
    
    Vector z(u.Size());
    adv.Mult(u, z);

   real_t t = 0.0;
   adv.SetTime(t);
   ode_solver->Init(adv);

   bool done = false;
   for (int ti = 0; !done; )
   {
      real_t dt_real = min(dt, t_final - t);
      ode_solver->Step(u, t, dt_real);
      ti++;

      done = (t >= t_final - 1e-8*dt);

      if (done || ti % vis_steps == 0)
      {
         cout << "time step: " << ti << ", time: " << t << endl;

         if (visualization)
         {
            sout << "solution\n" << mesh << u << flush;
         }

         if (visit)
         {
            dc->SetCycle(ti);
            dc->SetTime(t);
            dc->Save();
         }

         if (paraview)
         {
            pd->SetCycle(ti);
            pd->SetTime(t);
            pd->Save();
         }
      }
   }

   // 9. Save the final solution. This output can be viewed later using GLVis:
   //    "glvis -m ex9.mesh -g ex9-final.gf".
   {
      ofstream osol("ex9-final.gf");
      osol.precision(precision);
      u.Save(osol);
   }
    
    // Compute and print the L^2 norm of the error.
//    const IntegrationRule* irs[Geometry::NumGeom];
//    for (int i = 0; i < Geometry::NumGeom; i++)
//    {
//        if (i == 4)
//        {
//            // Tet Int Rule
//            irs[i] = &(IntRules.Get(i, 10));
//        }else
//        {
//            // Everything else
//            irs[i] = &(IntRules.Get(i, 16));
//        }
//    }
    cout << "Number of unknowns: " << fes.GetVSize() << endl;
    cout << "\n|| u_h - u ||_{L^2} = " << u.ComputeL2Error(u_coeff) << '\n' << endl;

   // 10. Free the used memory.
   delete pd;
   delete dc;

   return 0;
}


// Implementation of class FE_Evolution
FE_Evolution::FE_Evolution(BilinearForm &M_, BilinearForm &K_, const Vector &b_, FiniteElementSpace &vfes_)
   : TimeDependentOperator(M_.FESpace()->GetTrueVSize()),
     M(M_), K(K_), b(b_), z(height), vfes(vfes_)
{
    std::ofstream linear_sys_m("DG-Vector-Diffusion-Sys_m.txt");
    M.PrintMatlab(linear_sys_m);
    std::ofstream linear_sys_k("DG-Vector-Diffusion-Sys_k.txt");
    K.PrintMatlab(linear_sys_k);
    
    
    // Compute inverse mass matrix per element
    InverseIntegrator inv_mass(new MassIntegrator());
    invmass.resize(vfes.GetNE());
    for (int i=0; i<vfes.GetNE(); i++)
    {
        int dof = vfes.GetFE(i)->GetDof();
        //int dof = vfes.GetVDofs();
        //int dof = 6;
        invmass[i].SetSize(dof);
        inv_mass.AssembleElementMatrix(*vfes.GetFE(i),
                                       *vfes.GetElementTransformation(i),
                                       invmass[i]);
    }
    DenseMatrix Block_dense_mat;
    // Print each block of inverse mass
    for (int block = 0; block<vfes.GetNE(); block++)
    {
        Block_dense_mat = invmass[block];
        Block_dense_mat.Print();
        Block_dense_mat = 0.0;
    }
    
}

void FE_Evolution::Mult(const Vector &x, Vector &y) const
{
    /*
     // y = M^{-1} (K x + b)
     K.Mult(x, z);
     //z.Neg(); // z = -z //New change
     z += b;
     M_solver.Mult(z, y);*/
    
    // New implementation
    K.Mult(x, z); // Multiply x by stiffness matrix K
    z += b; // Add linear form RHS
    
    // Apply block inverser mass M^(-1) * (z)
    // Apply block inverse mass
    Vector zval; // z_loc, dof*num_eq
    int num_equations = 2;
    DenseMatrix current_zmat; // view of element auxiliary result, dof x num_eq
    DenseMatrix current_ymat; // view of element result, dof x num_eq
    Array<int> vdofs;
    for (int i=0; i<vfes.GetNE(); i++)
    {
       int dof = vfes.GetFE(i)->GetDof();
       vfes.GetElementVDofs(i, vdofs);
       z.GetSubVector(vdofs, zval);
       current_zmat.UseExternalData(zval.GetData(), dof, num_equations);
       current_ymat.SetSize(dof, num_equations);
       mfem::Mult(invmass[i], current_zmat, current_ymat);
       y.SetSubVector(vdofs, current_ymat.GetData());
        //y.Print();
        //std::cout << "hellp";
    }
    
}



FE_Evolution::~FE_Evolution()
{
 // not used
}


// Velocity coefficient
void velocity_function(const Vector &x, Vector &v)
{
   int dim = x.Size();

   switch (problem)
   {
    case 1: // Convection Only
    {
          // Translations in 1D, 2D, and 3D
          switch (dim)
          {
             case 1: v(0) = 1.0; break;
             case 2: v(0) = 1.0; v(1) = 1.0; break;
             case 3: break;
          }
          break;
     }
     case 2: // Diffusion Only
     {
         switch (dim)
         {
             case 1: v(0) = 1.0; break;
             case 2: v(0) = 1.0; v(1) = 1.0; break;
             case 3: break;
         }
           break;
     }
     case 3: // Vector Convection Diffusion
     {
         switch (dim)
         {
             case 1: break;
             case 2: v(0) = 1.0; v(1) = 1.0; break;
             case 3: break;
         }
         break;
     }
           
   }
}

// Initial condition
void u0_function(const Vector &x, Vector &u0)
{
   int dim = x.Size();

   switch (problem)
   {

      case 1: // Convection Only
       {
           if (dim == 1)
           {
               real_t x_ = M_PI*x(0);
               u0(0) =  sin(x_);
           }
           if (dim == 2)
           {
               real_t x_ = M_PI*x(0), y_ = M_PI*x(1);
               u0(0) =  sin(x_)*sin(y_);
               u0(1) =  cos(x_)*cos(y_);
           }
       }
       case 2: // Diffusion Only 
       {
           if (dim == 1)
           {
               real_t x_ = M_PI*x(0);
               u0(0) =  sin(x_);
           }
           if (dim == 2)
           {
               real_t x_ = M_PI*x(0), y_ = M_PI*x(1);
               u0(0) =  sin(x_)*sin(y_);
               u0(1) =  cos(x_)*cos(y_);
           }
       }
       case 3: // Vector Convection Diffusion
            
            // Initial vector field
            switch (dim)
            {
                case 1:
                {
                    real_t x_ = M_PI*x(0);
                    u0(0) =  sin(x_);
                    break;
                }
               case 2:
                {
                    real_t x_ = M_PI*x(0), y_ = M_PI*x(1);
                    u0(0) =  sin(x_)*sin(y_);
                    u0(1) =  cos(x_)*cos(y_);
                    break;
                }
               case 3: mfem_error("Dimension not valid"); break;
            }
            break;
   }
}


// Inflow boundary condition (zero for the problems considered in this example)
real_t inflow_function(const Vector &x)
{
   switch (problem)
   {
      case 0:
      case 1:
      case 2:
      case 3: return 0.0;
   }
   return 0.0;
}

// Forcing Funtion
void f_exact(const Vector &x, Vector &f_vec)
{
    
    int dim = x.Size();

    switch (problem)
    {
        case 1: //Convection only
        {
            if (dim == 1)
            {
                double v_x = 1.0;
                real_t x_ = M_PI*x(0);

                double conv_x = v_x * ( M_PI * cos(x_));
                f_vec(0) = conv_x;
                
                break;
            }
            if (dim == 2)
            {
                double v_x = 1.0;
                double v_y = 1.0;
                real_t x_ = M_PI*x(0), y_ = M_PI*x(1);
                
                
                double conv_x = v_x * ( M_PI * cos(x_)*sin(y_) ) + v_y * ( M_PI * sin(x_)*cos(y_) );
                double conv_y = v_x * ( -1.0 * M_PI * sin(x_)*cos(y_) ) + v_y * ( -1.0 * M_PI * cos(x_)*sin(y_) );
                
                f_vec(0) = conv_x;
                f_vec(1) = conv_y;
                
                break;
            }
        }
        case 2: // Diffusion Only
        {
            if (dim == 1)
            {
                real_t x_ = M_PI*x(0);
                f_vec(0) = M_PI*M_PI*sin(x_);
                break;
            }
            if (dim == 2)
            {
                real_t x_ = M_PI*x(0), y_ = M_PI*x(1);
                f_vec(0) = 2*M_PI*M_PI*sin(x_)*sin(y_);
                f_vec(1) = 2*M_PI*M_PI*cos(x_)*cos(y_);
                break;
            }
        }
        case 3:
            switch (dim)
            {
                case 1:
                case 2:
                {
                    double v_x = 1.0;
                    double v_y = 1.0;
                    double x_ = x(0);
                    double y_ = x(1);
                    const real_t f = M_PI;
                    const real_t s = M_PI_2;
                    
                    double conv_x = v_x * (-sin(2*f*x_ - s)*cos(2*f*y_ -s)*2*f)
                    + v_y * (-cos(2*f*x_ - s)*sin(2*f*y_ -s)*2*f);
                    double conv_y = v_x * (cos(f*x_)*sin(f*y_)*f) + v_y * (sin(f*x_)*cos(f*y_)*f);
                    
                    double laplace_x = -8*f*f*cos(2*f*x_ -s)*cos(2*f*y_ - s);
                    double laplace_y = -2*f*f*sin(f*x_)*sin(f*y_);
                    
                    f_vec(0) = conv_x + laplace_x;
                    f_vec(1) = conv_y + laplace_y;
                    
                }
                case 3:
                case 4:
            }
    }
}



