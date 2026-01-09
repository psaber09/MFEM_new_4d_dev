//                                DG Vector Advection Diffusion with Explicit Solver



#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>

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



int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   problem = 2;
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
   FiniteElementSpace fes(&mesh, &fec, dim);

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
    
//    std::ofstream ofs("mass_matrix.txt");
//    m.Print(ofs);
//    ofs.close();
    SparseMatrix M = m.SpMat();
    DenseMatrix elmat;
    //m.GetElementMatrix(e, elmat); // e is the element index


    //std::ofstream MassMatrix("DG_MassMatrix.txt");
    //Mdense.PrintMatlab(MassMatrix);

   std::ofstream linear_sys_m("DG-Vector-Diffusion-Sys_m.mat");
   m.PrintMatlab(linear_sys_m);
   std::ofstream linear_sys_k("DG-Vector-Diffusion-Sys_k.txt");
   k.PrintMatlab(linear_sys_k);
    
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
    SparseMatrix &Kmat = k.SpMat();
    SparseMatrix &Mmat = m.SpMat();
    Vector Bvec = b;

        // --- Precompute element-local mass inverses ---
        // For DG, dofs are local to each element (no sharing), so we can invert per element.
        int ne = mesh.GetNE();
        Array<int> el_dofs;
        DenseMatrix localM;
        DenseMatrix localMinv;
        Vector local_rhs;
        Vector local_sol;

        // We'll store inverted local mass matrices for each element
        std::vector<DenseMatrix> Minv_per_elem;
        Minv_per_elem.resize(ne);

        for (int el = 0; el < ne; ++el)
        {
            // Get the global DOF indices for element 'el'
            fes.GetElementDofs(el, el_dofs); // el_dofs.Size() = dofs_per_element * vdim
            int loc_nd = el_dofs.Size();

            localM.SetSize(loc_nd, loc_nd);

            // Fill localM from global Mmat
            for (int i = 0; i < loc_nd; ++i)
            {
                int gi = el_dofs[i];
                for (int j = 0; j < loc_nd; ++j)
                {
                    int gj = el_dofs[j];
                    // SparseMatrix::GetIJ entry access is not convenient; use operator() which is implemented
                    localM(i, j) = Mmat(gi, gj);
                }
            }

            // Invert localM (make a copy to invert)
            localMinv = localM;
            localMinv.Invert();
            // For stability, you may want to check condition number or add tiny regularization
//            if (!localMinv.Invert())
//            {
//                cerr << "Warning: local mass matrix inversion failed on element " << el
//                     << " (size " << loc_nd << "). Trying pseudo-inverse fallback." << endl;
//                // Fallback: use diagonal inversion
//                DenseMatrix fallback(loc_nd);
//                fallback = 0.0;
//                for (int ii = 0; ii < loc_nd; ++ii)
//                {
//                    double v = localM(ii, ii);
//                    fallback(ii, ii) = (fabs(v) > 1e-12) ? 1.0 / v : 0.0;
//                }
//                localMinv = fallback;
//            }

            Minv_per_elem[el].SetSize(loc_nd);
            Minv_per_elem[el] = localMinv;
        }

        cout << "Precomputed per-element mass inverses for " << ne << " elements." << endl;

        // Time stepping params
        int steps = max(1, int(t_final / dt));
        double t = 0.0;
        Vector dudt;

        // Temporary global vectors
        Vector z(fes.GetVSize());
        Vector temp(fes.GetVSize());

        // Forcing function f(x,t) — here zero for simplicity; user can set as needed
        // If you have a known forcing, fill 'rhs' each step with f and use that in RHS.
        //rhs = 0.0;

        // Main time loop (Forward Euler)
        for (int it = 0; it < steps; ++it)
        {
            t = it * dt;

            // Compute Ku = K * u (Kmat from bilinear form)
            Kmat.Mult(u, z);
            z *= -1.0;
            z += Bvec;

            // rhs_total = -Ku + f
            // reuse 'temp' to store rhs_total
            //temp = Ku;
            //temp *= -1.0;
            // add forcing if present: temp += rhs (if rhs is time-dependent, compute it)
            // temp += rhs; // currently rhs is zero

            // Now apply M^{-1} block-locally: dudt = M^{-1} * temp
            dudt.SetSize(fes.GetVSize());
            dudt = 0.0; // ensure cleared

            for (int el = 0; el < ne; ++el)
            {
                fes.GetElementDofs(el, el_dofs);
                int loc_nd = el_dofs.Size();

                local_rhs.SetSize(loc_nd);
                local_sol.SetSize(loc_nd);

                // fill local_rhs from global temp
                for (int i = 0; i < loc_nd; ++i)
                {
                    local_rhs(i) = z(el_dofs[i]);
                }

                // local_sol = localMinv * local_rhs
                local_sol.SetSize(loc_nd);
                Minv_per_elem[el].Mult(local_rhs, local_sol);

                // scatter local_sol into dudt
                for (int i = 0; i < loc_nd; ++i)
                {
                    dudt(el_dofs[i]) = local_sol(i);
                }
            } // end elements loop

            // Forward Euler update: u^{n+1} = u^n + dt * dudt
            u.Add(dt, dudt);

            // (Optional) output or diagnostics
            if ((it % 10) == 0)
            {
                cout << "Step " << it << " t=" << t << endl;
            }
        } // end time loop

        // Save final solution
        {
            ofstream sol_ofs("u_final.gf");
            u.Save(sol_ofs);
            cout << "Saved final GridFunction to u_final.gf" << endl;
        }

        // Clean up
        //delete Mmat;
        //delete Kmat;
    

   //FE_Evolution adv(m, k, b);
   /*
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
   }*/

   // 9. Save the final solution. This output can be viewed later using GLVis:
   //    "glvis -m ex9.mesh -g ex9-final.gf".
   {
      ofstream osol("ex9-final.gf");
      osol.precision(precision);
      u.Save(osol);
   }
    
    
    cout << "Number of unknowns: " << fes.GetVSize() << endl;
    cout << "\n|| u_h - u ||_{L^2} = " << u.ComputeL2Error(u_coeff) << '\n' << endl;

   // 10. Free the used memory.
   delete pd;
   delete dc;

   return 0;
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



