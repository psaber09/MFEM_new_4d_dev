//                                MFEM Example 14 - Parallel Version (L2 error)
//
// Compile with: make ex1p_L2
//
// Sample runs:  mpirun -np 4 ex1p_L2 -m ../data/inline-quad.mesh -o 0
//               mpirun -np 4 ex1p_L2 -m ../data/star.mesh -r 2 -rp 1 -o 2
//               mpirun -np 4 ex1p_L2 -m ../data/inline-hex.mesh -r 1 -rp 1 -o 2 -s 1
//               mpirun -np 4 ex1p_L2 -pa -r 1 -o 3
//
// Description:  DG discretization of the Poisson problem -Delta u = f with
//               homogeneous Dirichlet boundary conditions (imposed weakly),
//               where the exact solution is the product of x_i(1-x_i).
//
//               Problem setup and solution follow MFEM's parallel example
//               ex14p: the system is assembled with ParallelAssemble (no
//               essential dof elimination), preconditioned with BoomerAMG,
//               and solved with CG if sigma == -1 (symmetric) or GMRES
//               otherwise. The L^2 error of the DG solution is computed on
//               a custom high-order integration rule and reported globally.

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <unistd.h>


using namespace std;
using namespace mfem;

// Exact solution, E, and r.h.s., f. See below for implementation.
double E_exact(const Vector &);
double f_exact(const Vector &);

int dim;

int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
    
//    if (myid == 0) // Change 0 to whichever rank you want to debug
//    {
//        std::cerr << "Rank 0 PID: " << getpid() << " - waiting to attach..." << std::endl;
//        sleep(15); // Gives you 15 seconds to attach Xcode
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
    
   Hypre::Init();


   // 2. Parse command-line options.
   const char *mesh_file = "../data/star.mesh";
   int ref_levels = 0;     // serial refinements
   int par_ref_levels = 0; // parallel refinements
   int order = 1;
   real_t sigma = -1.0;  // SIPM sigma = -1; NSIPM sigma = 1;
   real_t kappa = -1.0;
   real_t eta = 0.0;
   real_t rel_tol = 1e-12;
   int max_iter = 500;
   int kdim = 10;
   bool pa = false;
   bool visualization = true;
   const char *device_config = "cpu";

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) >= 0.");
   args.AddOption(&sigma, "-s", "--sigma",
                  "One of the three DG penalty parameters, typically +1/-1."
                  " See the documentation of class DGDiffusionIntegrator.");
   args.AddOption(&kappa, "-k", "--kappa",
                  "One of the three DG penalty parameters, should be positive."
                  " Negative values are replaced with 3(order+1)(order+2).");
   args.AddOption(&eta, "-e", "--eta", "BR2 penalty parameter.");
   args.AddOption(&rel_tol, "-rtol", "--rel-tol",
                  "Relative tolerance for the Krylov solver.");
   args.AddOption(&max_iter, "-mi", "--max-iter",
                  "Maximum number of Krylov iterations.");
   args.AddOption(&kdim, "-kdim", "--krylov-dim",
                  "GMRES restart size (used when sigma != -1).");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.Parse();
   if (!args.Good())
   {
      if (myid == 0) { args.PrintUsage(cout); }
      return 1;
   }

   if (kappa < 0)
   {
      kappa = 3*(order+1)*(order+2);
   }
   if (myid == 0) { args.PrintOptions(cout); }

   // 3. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   if (myid == 0) { device.Print(); }

   // 4. Read the (serial) mesh from the given mesh file on all processors.
   Mesh mesh(mesh_file);
   dim = mesh.Dimension();

   // 5. Refine the serial mesh, then convert NURBS meshes to curved meshes.
   for (int l = 0; l < ref_levels; l++)
   {
      mesh.UniformRefinement();
   }
   if (mesh.NURBSext)
   {
      mesh.SetCurvature(max(order, 1));
   }

   // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel. Once the parallel mesh is defined, the
   //    serial mesh can be deleted.
   ParMesh pmesh(MPI_COMM_WORLD, mesh);
   mesh.Clear();
    
    
   for (int l = 0; l < par_ref_levels; l++)
   {
      pmesh.UniformRefinement();
   }
    
    // Number of elements per rank.
    {
       int ne = pmesh.GetNE();
       Array<int> all(num_procs);
       MPI_Gather(&ne, 1, MPI_INT, all.GetData(), 1, MPI_INT, 0, MPI_COMM_WORLD);
       if (myid == 0)
       {
          cout << "Elements per rank:";
          for (int r = 0; r < num_procs; r++) { cout << " " << all[r]; }
          cout << endl;
       }
    }
    
//    int nsf = pmesh.GetNSharedFaces(), tot = 0;
//    MPI_Allreduce(&nsf, &tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//    if (myid == 0) { cout << "Total shared faces: " << tot << endl; }

   // 7. Define a parallel finite element space on the parallel mesh. Here we
   //    use discontinuous finite elements of the specified order >= 0.
   const auto bt = pa ? BasisType::GaussLobatto : BasisType::GaussLegendre;
   DG_FECollection fec(order, dim, bt);
   ParFiniteElementSpace fespace(&pmesh, &fec);
   HYPRE_BigInt size = fespace.GlobalTrueVSize();
   if (myid == 0)
   {
      cout << "Number of unknowns: " << size << endl;
   }

   // 8. Set up the parallel linear form b(.) which corresponds to the
   //    right-hand side of the FEM linear system.
   ParLinearForm b(&fespace);
   ConstantCoefficient one(1.0);
   ConstantCoefficient beta(1.0);
   ConstantCoefficient zero(0.0);
   FunctionCoefficient f(f_exact);
   FunctionCoefficient E(E_exact);

   b.AddDomainIntegrator(new DomainLFIntegrator(f));
   b.AddBdrFaceIntegrator(
      new DGDirichletLFIntegrator(zero, one, sigma, kappa));
   b.Assemble();

   // 9. Define the solution vector x as a parallel finite element grid
   //    function corresponding to fespace. Initialize x with initial guess
   //    of zero.
   ParGridFunction x(&fespace);
   x.ProjectCoefficient(E);

   // Integration rules used for the L^2 error computation.
   const IntegrationRule *irs[Geometry::NumGeom];
   for (int i = 0; i < Geometry::NumGeom; i++)
   {
      if (i == Geometry::TETRAHEDRON)
      {
         irs[i] = &(IntRules.Get(i, 10)); // Tet rule
      }
      else
      {
         irs[i] = &(IntRules.Get(i, 16)); // Everything else
      }
   }

   // 10. Set up the bilinear form a(.,.) on the finite element space
   //     corresponding to the Laplacian operator -Delta, by adding the
   //     Diffusion domain integrator and the interior and boundary DG face
   //     integrators. Boundary conditions are imposed weakly in the form, so
   //     there is no need for dof elimination. After serial and parallel
   //     assembly we extract the corresponding parallel matrix A.
   ParBilinearForm a(&fespace);
   a.AddDomainIntegrator(new DiffusionIntegrator(beta));
   a.AddInteriorFaceIntegrator(new DGDiffusionIntegrator(beta, sigma, kappa));
   a.AddBdrFaceIntegrator(new DGDiffusionIntegrator(beta, sigma, kappa));
   if (eta > 0)
   {
      MFEM_VERIFY(!pa, "BR2 not yet compatible with partial assembly.");
      a.AddInteriorFaceIntegrator(new DGDiffusionBR2Integrator(fespace, eta));
      a.AddBdrFaceIntegrator(new DGDiffusionBR2Integrator(fespace, eta));
   }
   if (pa) { a.SetAssemblyLevel(AssemblyLevel::PARTIAL); }
   a.Assemble();
   a.Finalize();

   // 11. Define the parallel (hypre) matrix representing a(.,.) and, in the
   //     non-PA case, the BoomerAMG preconditioner.
   OperatorHandle A;

   std::unique_ptr<HypreBoomerAMG> amg;
   if (pa)
   {
      A.Reset(&a, false);
   }
   else
   {
      A.SetType(Operator::Hypre_ParCSR);
      a.ParallelAssemble(A);
      amg.reset(new HypreBoomerAMG(*A.As<HypreParMatrix>()));
      amg->SetPrintLevel(0);
   }

   // 12. Depending on the symmetry of A, define and apply a parallel PCG or
   //     GMRES solver for A x = b using the BoomerAMG preconditioner.
   bool converged = false;
   int num_iter = 0;
   if (sigma == -1.0)
   {
      CGSolver cg(MPI_COMM_WORLD);
      cg.SetRelTol(rel_tol);
      cg.SetAbsTol(0.0);
      cg.SetMaxIter(max_iter);
      cg.SetPrintLevel(1);
      cg.SetOperator(*A);
      if (amg) { cg.SetPreconditioner(*amg); }
      cg.Mult(b, x);
      converged = cg.GetConverged();
      num_iter = cg.GetNumIterations();
   }
   else
   {
      GMRESSolver gmres(MPI_COMM_WORLD);
      gmres.SetAbsTol(0.0);
      gmres.SetRelTol(rel_tol);
      gmres.SetMaxIter(max_iter);
      gmres.SetKDim(kdim);
      gmres.SetPrintLevel(1);
      gmres.SetOperator(*A);
      if (amg) { gmres.SetPreconditioner(*amg); }
      gmres.Mult(b, x);
      converged = gmres.GetConverged();
      num_iter = gmres.GetNumIterations();
   }

   // 13. Report solver convergence and the true relative residual.
   {
      Vector R(b.Size());
      A->Mult(x, R);
      R -= b;
      real_t rnorm = sqrt(InnerProduct(MPI_COMM_WORLD, R, R));
      real_t bnorm = sqrt(InnerProduct(MPI_COMM_WORLD, b, b));
      if (myid == 0)
      {
         cout << "Solver converged: " << (converged ? "yes" : "NO")
              << ", iterations: " << num_iter
              << ", true rel. residual: " << rnorm/bnorm << endl;
      }
   }

   // 14. Compute and print the global L^2 norm of the error.
   real_t err = x.ComputeL2Error(E, irs);
   if (myid == 0)
   {
      cout << "\n || E_h - E ||_{L^2} = " << err << '\n' << endl;
   }

   // 15. Save the refined mesh and the solution in parallel. This output can
   //     be viewed later using GLVis: "glvis -np <np> -m refined -g sol".
   {
      ostringstream mesh_name, sol_name;
      mesh_name << "refined." << setfill('0') << setw(6) << myid;
      sol_name << "sol." << setfill('0') << setw(6) << myid;

      ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(8);
      pmesh.Print(mesh_ofs);

      ofstream sol_ofs(sol_name.str().c_str());
      sol_ofs.precision(8);
      x.Save(sol_ofs);
   }

   // 16. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock.precision(8);
      sol_sock << "solution\n" << pmesh << x << flush;
   }

   return 0;
}



double E_exact(const Vector &x)
{
   double E_out = 0.0;

   if (dim == 2)
   {
      E_out = x(0)*(1. - x(0))*x(1)*(1. - x(1));
   }

   if (dim == 3)
   {
      E_out = x(0)*(1. - x(0)) * x(1)*(1. - x(1)) * x(2)*(1. - x(2));
   }

   if (dim == 4)
   {
      E_out = x(0)*(1. - x(0)) * x(1)*(1. - x(1)) * x(2)*(1. - x(2)) *
              x(3)*(1. - x(3));
   }

   return E_out;
}


double f_exact(const Vector &x)
{
   double f_out = 0.0;

   if (dim == 2)
   {
      f_out = 2.0*(x(0)*(1. - x(0)) + x(1)*(1. - x(1)));
   }

   if (dim == 3)
   {
      f_out = 2.0*((x(1)*(1. - x(1))*x(2)*(1. - x(2))) +
                   (x(0)*(1. - x(0))*x(2)*(1. - x(2))) +
                   (x(0)*(1. - x(0))*x(1)*(1. - x(1))));
   }

   if (dim == 4)
   {
      f_out = 2.0*((x(1)*(1. - x(1))*x(2)*(1. - x(2))*x(3)*(1. - x(3))) +
                   (x(0)*(1. - x(0))*x(2)*(1. - x(2))*x(3)*(1. - x(3))) +
                   (x(0)*(1. - x(0))*x(1)*(1. - x(1))*x(3)*(1. - x(3))) +
                   (x(0)*(1. - x(0))*x(1)*(1. - x(1))*x(2)*(1. - x(2))));
   }

   return f_out;
}
