//                       MFEM Example 3 - Parallel Version (SkwGrad, 4D)
//
// Compile with: make ex3p_SkwGrad
//
// Sample runs:  mpirun -np 4 ex3p_SkwGrad -m ../data/cube4d_96.MFEM -ref 1
//               mpirun -np 4 ex3p_SkwGrad -m ../data/cube4d_96.MFEM -o 2 -pr 1
//               mpirun -np 4 ex3p_SkwGrad -m ../data/cube4d_96.MFEM -sc -check
//
// Description:  Parallel version of ex3_SkwGrad.cpp. Solves the second order
//               definite problem
//                  SkwGrad SkwGrad E + E = f
//               with boundary condition E x n = <given tangential field>.
//               A known exact solution E (given in 4D, with 3D/2D variants
//               kept from the serial example) is used to compute f. We
//               discretize with the SkwGrad_FECollection, using the
//               SkwGradSkwGradIntegrator and the vector FE mass integrator.
//
//               The structure follows ex3p.cpp. The hypre AMS preconditioner
//               used there is specific to Nedelec spaces in 2D/3D and cannot
//               be used here; this example uses CG preconditioned with an
//               l1-Jacobi smoother by default (see -prec).

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <unistd.h>

using namespace std;
using namespace mfem;

// Exact solution, E, and r.h.s., f. See below for implementation.
void E_exact(const Vector &, Vector &);
void f_exact(const Vector &, Vector &);
real_t freq = 1.0, kappa;
int dim;

int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
        
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();

//   if (myid == 3) // Change 0 to whichever rank you want to debug
//   {
//       std::cerr << "Rank 3 PID: " << getpid() << " - waiting to attach..." << std::endl;
//       sleep(15); // Gives you 15 seconds to attach Xcode
//   }
//   MPI_Barrier(MPI_COMM_WORLD);
    
//    if (myid == 0)
//    {
//       for (int o = 0; o < 24; o++)
//       {
//          const DenseMatrix &M = ND_PentDofTransformation::GetFaceTransform(o);
//          cout << "tet orientation " << o << "  size " << M.Height() << "x" << M.Width()
//               << "  det " << M.Det() << endl;
//          M.Print(cout);
//       }
//    }

   Hypre::Init();


   // 2. Parse command-line options.
   const char *mesh_file = "../data/cube4d_96.MFEM";
   int order = 1;
   int ref_levels = 0;       // sequential (serial) refinement levels
   int par_ref_levels = 0;   // parallel refinement levels
   bool static_cond = false;
   bool pa = false;
   bool nc = false;
   const char *device_config = "cpu";
   bool visualization = false;
   double tol = 1e-12;
   int max_iter = 10000;
   int prec_type = 0;        // 0/1: l1-Jacobi, 2: BoomerAMG

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&freq, "-f", "--frequency", "Set the frequency for the exact"
                  " solution (used for the 2D/3D variants).");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&nc, "-nc", "--non-conforming", "-c",
                  "--conforming",
                  "Mark the mesh as nonconforming before partitioning.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization (2D/3D only).");
   args.AddOption(&ref_levels, "-ref", "--ref-levels",
                  "Number of sequential (serial) uniform refinements.");
   args.AddOption(&par_ref_levels, "-pr", "--parrefinement",
                  "Number of parallel uniform refinements.");
   args.AddOption(&tol, "-tol", "--tol", "Relative tolerance for PCG.");
   args.AddOption(&max_iter, "-mi", "--max-iter", "Maximum PCG iterations.");
   args.AddOption(&prec_type, "-prec", "--preconditioner",
                  "Preconditioner: 0 or 1 = l1-Jacobi, 2 = BoomerAMG.");


   args.Parse();
   if (!args.Good())
   {
      if (myid == 0)
      {
         args.PrintUsage(cout);
      }
      return 1;
   }
   if (myid == 0)
   {
      args.PrintOptions(cout);
   }
   kappa = 1.0;   // as in the serial ex3_SkwGrad.cpp

   // 3. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   if (myid == 0) { device.Print(); }

   // 4. Read the (serial) mesh from the given mesh file on all processors.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   dim = mesh->Dimension();
   int sdim = mesh->SpaceDimension();
   MFEM_VERIFY(dim >= 2 && dim <= 4 && sdim == dim,
               "E_exact/f_exact are only defined for 2D, 3D and 4D meshes.");
   if (nc)
   {
      // Can set to false to use conformal refinement for simplices.
      mesh->EnsureNCMesh(true);
   }

   // 5. Refine the serial mesh on all processors to increase the resolution.
   if (myid == 0) { cout << "ref_levels = " << ref_levels << endl; }
   for (int l = 0; l < ref_levels; l++)
   {
      mesh->UniformRefinement();
   }

   // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution. Once the
   //    parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   for (int l = 0; l < par_ref_levels; l++)
   {
      pmesh->UniformRefinement();
   }
   pmesh->PrintInfo(cout);

   // 7. Define a parallel finite element space on the parallel mesh. Here we
   //    use the SkwGrad finite elements of the specified order.
   FiniteElementCollection *fec = new SkwGrad_FECollection(order, dim);
   ParFiniteElementSpace *fespace = new ParFiniteElementSpace(pmesh, fec);
   HYPRE_BigInt size = fespace->GlobalTrueVSize();
   if (myid == 0)
   {
      cout << "Number of finite element unknowns: " << size << endl;
   }

   // 7a. Optional consistency test of the shared dofs: project the exact
   //     solution, restrict to true dofs and prolongate back. The result
   //     must reproduce the projected function on every processor.
   {
      VectorFunctionCoefficient Echk(sdim, E_exact);
      ParGridFunction g(fespace), g2(fespace);
      g.ProjectCoefficient(Echk);
      Vector tv(fespace->GetTrueVSize());
      g.GetTrueDofs(tv);
      g2.SetFromTrueDofs(tv);
      g2 -= g;

//      int nbad = 0;
//      for (int i = 0; i < g2.Size(); i++)
//      {
//         if (fabs(g2(i)) > 1e-10) { nbad++; }
//      }
//      int nbad_tot;
//      MPI_Allreduce(&nbad, &nbad_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//      if (myid == 0)
//      {
//         cout << "Shared dof check: inconsistent dofs = " << nbad_tot << endl;
//      }
       
       int cnt[4] = {0, 0, 0, 0};   // vertex, edge, triangle, tet
       const int e0 = fespace->GetNVDofs();
       const int p0 = e0 + fespace->GetNEDofs();
       const int f0 = p0 + fespace->GetNPDofs();
       for (int i = 0; i < g2.Size(); i++)
       {
          if (fabs(g2(i)) > 1e-10) { cnt[i < e0 ? 0 : i < p0 ? 1 : i < f0 ? 2 : 3]++; }
       }
       MPI_Allreduce(MPI_IN_PLACE, cnt, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
       if (myid == 0) { cout << "bad dofs  edge:" << cnt[1] << " tri:" << cnt[2]
                             << " tet:" << cnt[3] << endl; }
   }

   // 8. Determine the list of true (i.e. parallel conforming) essential
   //    boundary dofs. In this example, the boundary conditions are defined
   //    by marking all the boundary attributes from the mesh as essential
   //    (Dirichlet) and converting them to a list of true dofs.
   Array<int> ess_tdof_list;
   Array<int> ess_bdr;
   if (pmesh->bdr_attributes.Size())
   {
      ess_bdr.SetSize(pmesh->bdr_attributes.Max());
      ess_bdr = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 9. Set up the parallel linear form b(.) which corresponds to the
   //    right-hand side of the FEM linear system, which in this case is
   //    (f,phi_i) where f is given by the function f_exact and phi_i are the
   //    basis functions in the finite element fespace.
   VectorFunctionCoefficient f(sdim, f_exact);
   ParLinearForm *b = new ParLinearForm(fespace);
   b->AddDomainIntegrator(new VectorFEDomainLFIntegrator(f));
   b->Assemble();

   // 10. Define the solution vector x as a parallel finite element grid
   //     function corresponding to fespace. Initialize x by projecting the
   //     exact solution. Note that only values from the boundary will be used
   //     when eliminating the non-homogeneous boundary condition to modify the
   //     r.h.s. vector b.
   ParGridFunction x(fespace);
   VectorFunctionCoefficient E(sdim, E_exact);
   x.ProjectCoefficient(E);

   // Integration rules for the error computation, as in the serial example.
   // (In the serial code the triangle rule set for geometry 2 is overwritten
   // by the "everything else" branch, so the effective rules are: order 10
   // for tetrahedra, order 16 for everything else. That is kept here.)
   const IntegrationRule *irs[Geometry::NumGeom];
   for (int i = 0; i < Geometry::NumGeom; i++)
   {
      if (i == 4) { irs[i] = &(IntRules.Get(i, 10)); }
      else        { irs[i] = &(IntRules.Get(i, 16)); }
   }

   {
      real_t err0 = x.ComputeL2Error(E, irs);
      if (myid == 0)
      {
         cout << "\n Initial || E_h - E ||_{L^2} = " << err0 << '\n' << endl;
      }
   }

   // 11. Set up the parallel bilinear form corresponding to the operator
   //     SkwGrad muinv SkwGrad + sigma I, by adding the SkwGrad-SkwGrad and
   //     the mass domain integrators.
   Coefficient *muinv = new ConstantCoefficient(1.0);
   Coefficient *sigma = new ConstantCoefficient(1.0);
   ParBilinearForm *a = new ParBilinearForm(fespace);
   if (pa) { a->SetAssemblyLevel(AssemblyLevel::PARTIAL); }
   a->AddDomainIntegrator(new SkwGradSkwGradIntegrator(*muinv));
   a->AddDomainIntegrator(new VectorFEMassIntegrator(*sigma));

   // 12. Assemble the parallel bilinear form and the corresponding linear
   //     system, applying any necessary transformations such as: parallel
   //     assembly, eliminating boundary conditions, applying conforming
   //     constraints for non-conforming AMR, static condensation, etc.
   if (static_cond) { a->EnableStaticCondensation(); }
   a->Assemble();

   OperatorPtr A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);

   // 13. Solve the system AX=B using PCG. There is no AMS-type preconditioner
   //     for this space, so use Jacobi (default) or BoomerAMG. In the partial
   //     assembly case use Jacobi preconditioning.
   CGSolver cg(MPI_COMM_WORLD);
   cg.SetRelTol(tol);
   cg.SetMaxIter(max_iter);
   cg.SetPrintLevel(1);
   cg.SetOperator(*A);

   Solver *prec = NULL;
   if (pa)
   {
      prec = new OperatorJacobiSmoother(*a, ess_tdof_list);
   }
   else
   {
      HypreParMatrix &Ah = *A.As<HypreParMatrix>();
      if (myid == 0)
      {
         cout << "Size of linear system: " << Ah.GetGlobalNumRows() << endl;
      }
      if (prec_type == 2)
      {
         HypreBoomerAMG *amg = new HypreBoomerAMG(Ah);
         amg->SetPrintLevel(0);
         prec = amg;
      }
      else
      {
         prec = new HypreSmoother(Ah, HypreSmoother::l1Jacobi);
      }
   }
   cg.SetPreconditioner(*prec);
   cg.Mult(B, X);
   if (myid == 0) { cout << "PCG iterations: " << cg.GetNumIterations() << endl; }

   // 14. Recover the parallel grid function corresponding to X. This is the
   //     local finite element solution on each processor.
   a->RecoverFEMSolution(X, *b, x);

   // 15. Compute and print the L^2 norm of the error.
   {
      real_t error = x.ComputeL2Error(E, irs);
      if (myid == 0)
      {
         cout << "\n|| E_h - E ||_{L^2} = " << error << '\n' << endl;
      }
   }

   // 16. Save the refined mesh and the solution in parallel. This output can
   //     be viewed later using GLVis: "glvis -np <np> -m mesh -g sol".
   {
      ostringstream mesh_name, sol_name;
      mesh_name << "mesh." << setfill('0') << setw(6) << myid;
      sol_name << "sol." << setfill('0') << setw(6) << myid;

      ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(8);
      pmesh->Print(mesh_ofs);

      ofstream sol_ofs(sol_name.str().c_str());
      sol_ofs.precision(8);
      x.Save(sol_ofs);
   }

   // 17. Send the solution by socket to a GLVis server (2D/3D only).
   if (visualization && dim <= 3)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock.precision(8);
      sol_sock << "solution\n" << *pmesh << x << flush;
   }

   // 18. Free the used memory.
   delete prec;
   delete a;
   delete sigma;
   delete muinv;
   delete b;
   delete fespace;
   delete fec;
   delete pmesh;

   return 0;
}


void E_exact(const Vector &x, Vector &E)
{
   if (dim == 4)
   {
      E(0) =  sin(M_PI*x(0))*cos(M_PI*x(1))*cos(M_PI*x(2))*cos(M_PI*x(3));
      E(1) = -cos(M_PI*x(0))*sin(M_PI*x(1))*cos(M_PI*x(2))*cos(M_PI*x(3));
      E(2) =  cos(M_PI*x(0))*cos(M_PI*x(1))*sin(M_PI*x(2))*cos(M_PI*x(3));
      E(3) = -cos(M_PI*x(0))*cos(M_PI*x(1))*cos(M_PI*x(2))*sin(M_PI*x(3));
   }
   else if (dim == 3)
   {
      E(0) = sin(kappa * x(1));
      E(1) = sin(kappa * x(2));
      E(2) = sin(kappa * x(0));
   }
   else
   {
      E(0) = 2.0*(x(0) + x(1));
      E(1) = 0.5*(x(0) + x(1));
   }
}

void f_exact(const Vector &x, Vector &f)
{
   if (dim == 4)
   {
      const real_t s_factor = 0.5;
      const real_t c = 1.0 + 4.0*s_factor*M_PI*M_PI;

      f(0) =  c*sin(M_PI*x(0))*cos(M_PI*x(1))*cos(M_PI*x(2))*cos(M_PI*x(3));
      f(1) = -c*cos(M_PI*x(0))*sin(M_PI*x(1))*cos(M_PI*x(2))*cos(M_PI*x(3));
      f(2) =  c*cos(M_PI*x(0))*cos(M_PI*x(1))*sin(M_PI*x(2))*cos(M_PI*x(3));
      f(3) = -c*cos(M_PI*x(0))*cos(M_PI*x(1))*cos(M_PI*x(2))*sin(M_PI*x(3));
   }
   else if (dim == 3)
   {
      f(0) = (1. + kappa * kappa) * sin(kappa * x(1));
      f(1) = (1. + kappa * kappa) * sin(kappa * x(2));
      f(2) = (1. + kappa * kappa) * sin(kappa * x(0));
   }
   else
   {
      f(0) = (1. + kappa * kappa) * sin(kappa * x(1));
      f(1) = (1. + kappa * kappa) * sin(kappa * x(0));
      if (x.Size() == 3) { f(2) = 0.0; }
   }
}
