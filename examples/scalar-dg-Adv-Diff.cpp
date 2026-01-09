//                                DG Scalar Advection Diffusion
//
//
//
//

//
// Description:  This example code solves the time-dependent advection diffusion equation
//               du/dt + v.grad(u) + delta(u) = f, where v is a given fluid velocity, and
//               u0(x)=u(0,x) is a given initial condition.
//
//               The example demonstrates the use of Discontinuous Galerkin (DG)
//               bilinear forms in MFEM (face integrators), the use of implicit
//               and explicit ODE time integrators.

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
//#include "vector-dg-diffusion.hpp"


using namespace std;
using namespace mfem;

// Choice for the problem setup. The fluid velocity, initial condition and
// inflow boundary condition are chosen based on this parameter.
int problem;

// Velocity coefficient
void velocity_function(const Vector &x, Vector &v);

// Initial condition
real_t u0_function(const Vector &x);

// Inflow boundary condition
real_t inflow_function(const Vector &x);

// Forcing Function
real_t f_exact(const Vector &x);

// Mesh bounding box
Vector bb_min, bb_max;

class DG_Solver : public Solver
{
private:
   SparseMatrix &M, &K, A;
   GMRESSolver linear_solver;
   BlockILU prec;
   real_t dt;
public:
   DG_Solver(SparseMatrix &M_, SparseMatrix &K_, const FiniteElementSpace &fes)
      : M(M_),
        K(K_),
        prec(fes.GetTypicalFE()->GetDof(),
             BlockILU::Reordering::MINIMUM_DISCARDED_FILL),
        dt(-1.0)
   {
      linear_solver.iterative_mode = true;
      linear_solver.SetRelTol(1e-9);
      linear_solver.SetAbsTol(0.0);
      linear_solver.SetMaxIter(100);
      linear_solver.SetPrintLevel(0);
      linear_solver.SetPreconditioner(prec);
   }

   void SetTimeStep(real_t dt_)
   {
      if (dt_ != dt)
      {
         dt = dt_;
         // Form operator A = M - dt*K
         A = K;
         A *= -dt;
         A += M;

         // this will also call SetOperator on the preconditioner
         linear_solver.SetOperator(A);
      }
   }

   void SetOperator(const Operator &op) override
   {
      linear_solver.SetOperator(op);
   }

   void Mult(const Vector &x, Vector &y) const override
   {
      linear_solver.Mult(x, y);
   }
};

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
   Solver *M_prec;
   CGSolver M_solver;
   DG_Solver *dg_solver;

   mutable Vector z;

public:
   FE_Evolution(BilinearForm &M_, BilinearForm &K_, const Vector &b_);

   void Mult(const Vector &x, Vector &y) const override;
   void ImplicitSolve(const real_t dt, const Vector &x, Vector &k) override;

   ~FE_Evolution() override;
};


int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   problem = 2;
   const char *mesh_file = "../data/cube2d_2.MFEM";
   int ref_levels = 0;
   int order = 1;
   real_t sigma = -1.0;
   //real_t kappa = -1.0;
   real_t kappa = (order+1)*(order+1);
   bool pa = false;
   bool ea = false;
   bool fa = false;
   const char *device_config = "cpu";
   int ode_solver_type = 4;
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
   DG_FECollection fec(order, dim, BasisType::GaussLobatto);
   FiniteElementSpace fes(&mesh, &fec);

   cout << "Number of unknowns: " << fes.GetVSize() << endl;

   // 6. Set up and assemble the bilinear and linear forms corresponding to the
   //    DG discretization. The DGTraceIntegrator involves integrals over mesh
   //    interior faces.
   VectorFunctionCoefficient velocity(dim, velocity_function);
   FunctionCoefficient inflow(inflow_function);
   FunctionCoefficient u0(u0_function);
    
   ConstantCoefficient one(1.0);
   ConstantCoefficient neg_one(-1.0);
   ConstantCoefficient diff_coeff(-1.0);
   ConstantCoefficient zero(0.0);
   
   // Forcing Function Coefficent 
   FunctionCoefficient f(f_exact);

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
    
   if (problem == 0) // Convection and Diffusion
   {
       std::cout << "Convection and Diffusion" << std::endl;
       // Time-dependant term
       m.AddDomainIntegrator(new MassIntegrator);
       
       // Convection term
       constexpr real_t alpha = -1.0;
       k.AddDomainIntegrator(new ConvectionIntegrator(velocity, alpha));
       k.AddInteriorFaceIntegrator(
          new NonconservativeDGTraceIntegrator(velocity, alpha));
       k.AddBdrFaceIntegrator(
          new NonconservativeDGTraceIntegrator(velocity, alpha));
        
       //Diffusion term
       k.AddDomainIntegrator(new DiffusionIntegrator(diff_coeff));
       k.AddInteriorFaceIntegrator(new DGDiffusionIntegrator(diff_coeff, sigma, kappa));
       k.AddBdrFaceIntegrator(new DGDiffusionIntegrator(diff_coeff, sigma, kappa));
       
       // RHS Forcing Term
       b.AddDomainIntegrator(new DomainLFIntegrator(f));
       b.AddBdrFaceIntegrator(
           new DGDirichletLFIntegrator(inflow, diff_coeff, sigma, kappa));
       
   }
   if (problem == 1) // Convection Only
   {
       std::cout << "Convection Only" << std::endl;
       // Time-dependant term
       m.AddDomainIntegrator(new MassIntegrator);
       // Convection Term
       constexpr real_t alpha = -1.0;
       k.AddDomainIntegrator(new ConvectionIntegrator(velocity, alpha));
       k.AddInteriorFaceIntegrator(
          new NonconservativeDGTraceIntegrator(velocity, alpha));
       k.AddBdrFaceIntegrator(
          new NonconservativeDGTraceIntegrator(velocity, alpha));
        
       // RHS Forcing Term
       b.AddDomainIntegrator(new DomainLFIntegrator(f));
       //b.AddBdrFaceIntegrator(new BoundaryFlowIntegrator(inflow, velocity, -alpha)); // ?????

   }
   if (problem == 2) // Diffusion Only
   {
        std::cout << "Diffusion Only" << std::endl;
        // Time-dependant term
        m.AddDomainIntegrator(new MassIntegrator);
            
        // Diffusion term
        k.AddDomainIntegrator(new DiffusionIntegrator(diff_coeff));
        k.AddInteriorFaceIntegrator(new DGDiffusionIntegrator(diff_coeff, sigma, kappa));
        k.AddBdrFaceIntegrator(new DGDiffusionIntegrator(diff_coeff, sigma, kappa));
       
        // RHS Forcing Term
        b.AddDomainIntegrator(new DomainLFIntegrator(f));
        b.AddBdrFaceIntegrator(new DGDirichletLFIntegrator(inflow, diff_coeff, sigma, kappa));
   }


   int skip_zeros = 0;
   m.Assemble(skip_zeros);
   m.Finalize(skip_zeros);

   k.Assemble(skip_zeros);
   k.Finalize(skip_zeros);
    
   SparseMatrix M = m.SpMat();
   std::ofstream linear_sys_m("DG-Scalar-Mass-Matrix.txt");
   m.PrintMatlab(linear_sys_m);
    
   b.Assemble();


   // 7. Define the initial conditions, save the corresponding grid function to
   //    a file and (optionally) save data in the VisIt format and initialize
   //    GLVis visualization.
   GridFunction u(&fes);
   u.ProjectCoefficient(u0);

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
   FE_Evolution adv(m, k, b);

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
    cout << "Number of unknowns: " << fes.GetVSize() << endl;
    cout << "\n|| u_h - u ||_{L^2} = " << u.ComputeL2Error(u0) << '\n' << endl;

   // 10. Free the used memory.
   delete pd;
   delete dc;

   return 0;
}


// Implementation of class FE_Evolution
FE_Evolution::FE_Evolution(BilinearForm &M_, BilinearForm &K_, const Vector &b_)
   : TimeDependentOperator(M_.FESpace()->GetTrueVSize()),
     M(M_), K(K_), b(b_), z(height)
{
   Array<int> ess_tdof_list;
   if (M.GetAssemblyLevel() == AssemblyLevel::LEGACY)
   {
      M_prec = new DSmoother(M.SpMat());
      M_solver.SetOperator(M.SpMat());
      dg_solver = new DG_Solver(M.SpMat(), K.SpMat(), *M.FESpace());
   }
   else
   {
      M_prec = new OperatorJacobiSmoother(M, ess_tdof_list);
      M_solver.SetOperator(M);
      dg_solver = NULL;
   }
   M_solver.SetPreconditioner(*M_prec);
   M_solver.iterative_mode = true;
   M_solver.SetRelTol(1e-9);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(100);
   M_solver.SetPrintLevel(0);
}

void FE_Evolution::Mult(const Vector &x, Vector &y) const
{
   // y = M^{-1} (K x + b)
   K.Mult(x, z);
   z += b;
   M_solver.Mult(z, y);
}

void FE_Evolution::ImplicitSolve(const real_t dt, const Vector &x, Vector &k)
{
   MFEM_VERIFY(dg_solver != NULL,
               "Implicit time integration is not supported with partial assembly");
   K.Mult(x, z);
   //z.Neg(); // z = -z //New change
   z += b;
   dg_solver->SetTimeStep(dt);
   dg_solver->Mult(z, k);
}

FE_Evolution::~FE_Evolution()
{
   delete M_prec;
   delete dg_solver;
}


// Velocity coefficient
void velocity_function(const Vector &x, Vector &v)
{
   int dim = x.Size();

   switch (problem)
   {
    case 0: // Convecton - Diffusion
    {
         // Translations in 1D, 2D, and 3D
         switch (dim)
         {
            case 1: v(0) = 1.0; break;
            case 2: v(0) = 1.0; v(1) = 1.0; break;
            case 3: v(0) = sqrt(3./6.); v(1) = sqrt(2./6.); v(2) = sqrt(1./6.);
               break;
         }
         break;
    }
    case 1: // Convection Only
    {
          // Translations in 1D, 2D, and 3D
          switch (dim)
          {
             case 1: v(0) = 1.0; break;
             case 2: v(0) = 1.0; v(1) = 1.0; break;
             case 3: v(0) = sqrt(3./6.); v(1) = sqrt(2./6.); v(2) = sqrt(1./6.); break;
          }
          break;
     }
     case 2: // Diffusion Only
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
real_t u0_function(const Vector &x)
{
   int dim = x.Size();

   switch (problem)
   {
      case 0: // Convection - Diffusion
       {
           const real_t f = M_PI;
           return sin(f*x(0))*sin(f*x(1));
       }
      case 1: // Convection Only
       {
           if (dim == 2)
           {
               const real_t f = M_PI;
               return sin(f*x(0))*sin(f*x(1));
           }
       }
       case 2: // Diffusion Only 
       {
           if (dim == 2)
           {
               return x(0)*(1.0 - x(0)) * x(1)*(1.0 - x(1));
           }
       }
   }
   return 0.0;
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
real_t f_exact(const Vector &x)
{

    int dim = x.Size();

    switch (problem)
    {
       case 0:
       {
          // Conv - diff
          switch (dim)
          {
             case 2: 
              {
                  double v_x = 1.0;
                  double v_y = 1.0;
                  double x_ = x(0);
                  double y_ = x(1);
                  const real_t f = M_PI;
                  
                  return ( v_x*(f*cos(f*x_)*sin(f*y_)) + v_y*(f*sin(f*x_)*cos(f*y_)) + 2*f*f*sin(f*x_)*sin(f*y_) );
                  break;
              }
          }
       }
       case 1: //Convection only
       {
           double v_x = 1.0;
           double v_y = 1.0;
           double x_ = x(0);
           double y_ = x(1);
           const real_t f = M_PI;
           
           return ( v_x*(f*cos(f*x_)*sin(f*y_)) + v_y*(f*sin(f*x_)*cos(f*y_)) );
           break;
       }
       case 2: // Diffusion Only
       {
            
           return  (2.0*( x(0)*(1. - x(0)) + x(1)*(1. - x(1))) );

           break;
       }
    }
}



