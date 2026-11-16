//                       MFEM Example 3 (Curl) - Parallel Version
//
// Compile with: make ex3p_curl
//
// Sample runs:  mpirun -np 4 ex3p_curl -m <4D mesh> -o 1
//               mpirun -np 4 ex3p_curl -m <4D mesh> -o 2 -ref 1 -rp 1
//               mpirun -np 4 ex3p_curl -m <4D mesh> -o 1 -prec 2
//               mpirun -np 4 ex3p_curl -m <4D mesh> -o 2 -pa
//
// Description:  Parallel version of ex3_Curl.cpp. Solves the second order
//               definite Maxwell-type problem curl curl E + E = f with an
//               exact solution E (a 4D antisymmetric-matrix-valued field
//               stored as 12 components), using H(curl) finite elements.
//               The structure follows the parallel MFEM example ex3p.
//
//               Solver (non-PA mode): PCG with a parallel preconditioner,
//               selected by -prec:
//                 0 : HypreSmoother, l1 Gauss-Seidel   (default, analogue of
//                     the serial GSSmoother)
//                 1 : HypreSmoother, l1 Jacobi
//                 2 : HypreBoomerAMG
//               In PA mode, CG with OperatorJacobiSmoother is used.
//
//               Note: HypreAMS (used in ex3p) is only defined for 3D, so it
//               is not used here.

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;
using namespace mfem;

// Exact solution, E, and r.h.s., f. See below for implementation.
void E_exact_vec(const Vector &x, Vector &E);
void E_exact(const Vector &, DenseMatrix &);
void f_exact(const Vector &, DenseMatrix &);
real_t freq = 1.0, kappa;
int dim;

int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();

   // 2. Parse command-line options.
   const char *mesh_file = "../data/beam-tet.mesh";
   int order = 1;
   int ref_levels = 0;     // serial refinements
   int par_ref_levels = 0; // parallel refinements
   int prec_type = 0;
   real_t rel_tol = 1e-9;  // serial code: PCG(..., 1e-18, ...) => sqrt = 1e-9
   int max_iter = 5000;
   bool static_cond = false;
   bool pa = false;
   bool nc = false;
   const char *device_config = "cpu";
   bool visualization = true;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&freq, "-f", "--frequency", "Set the frequency for the exact"
                  " solution.");
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
                  "Enable or disable GLVis visualization.");
   args.AddOption(&ref_levels, "-ref", "--ref-levels",
                  "Number of uniform refinements of the serial mesh.");
   args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                  "Number of uniform refinements of the parallel mesh.");
   args.AddOption(&prec_type, "-prec", "--preconditioner",
                  "Preconditioner: 0 - l1 Gauss-Seidel, 1 - l1 Jacobi, "
                  "2 - BoomerAMG.");
   args.AddOption(&rel_tol, "-rtol", "--rel-tol",
                  "Relative tolerance for the Krylov solver.");
   args.AddOption(&max_iter, "-mi", "--max-iter",
                  "Maximum number of Krylov iterations.");

   args.Parse();
   if (!args.Good())
   {
      if (myid == 0) { args.PrintUsage(cout); }
      return 1;
   }
   if (myid == 0) { args.PrintOptions(cout); }
   kappa = 1.0;

   // 3. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   if (myid == 0) { device.Print(); }

   // 4. Read the (serial) mesh from the given mesh file on all processors.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   dim = mesh->Dimension();
   int sdim = mesh->SpaceDimension();
   if (nc)
   {
      if (myid == 0)
      {
         cout << "Inside of Non-Conforming Flag" << endl;
      }
      mesh->EnsureNCMesh(true);
   }

   // 5. Refine the serial mesh on all processors.
   for (int l = 0; l < ref_levels; l++)
   {
      mesh->UniformRefinement();
   }

   // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel. Once the parallel mesh is defined, the
   //    serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   for (int l = 0; l < par_ref_levels; l++)
   {
      pmesh->UniformRefinement();
   }

   // 7. Define a parallel finite element space on the parallel mesh. Here we
   //    use the H(curl) finite elements of the specified order.
   FiniteElementCollection *fec = new HCurl_FECollection(order, dim);
   ParFiniteElementSpace *fespace = new ParFiniteElementSpace(pmesh, fec);
   HYPRE_BigInt size = fespace->GlobalTrueVSize();
   if (myid == 0)
   {
      cout << "Number of finite element unknowns: " << size << endl;
   }

   // 8. Determine the list of true (i.e. parallel conforming) essential
   //    boundary dofs, marking all boundary attributes as essential.
   Array<int> ess_tdof_list;
   if (pmesh->bdr_attributes.Size())
   {
      Array<int> ess_bdr(pmesh->bdr_attributes.Max());
      ess_bdr = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 9. Set up the parallel linear form b(.) which corresponds to the
   //    right-hand side of the FEM linear system.
   MatrixFunctionCoefficient f(sdim, f_exact);
   VectorFunctionCoefficient solVec(12, E_exact_vec);

   ParLinearForm *b = new ParLinearForm(fespace);
   b->AddDomainIntegrator(new MatFEDomainLFIntegrator(f));
   b->Assemble();

   // 10. Define the solution vector x as a parallel finite element grid
   //     function. Initialize x by projecting the exact solution. Only the
   //     boundary values are used when eliminating the non-homogeneous
   //     boundary condition from the r.h.s. vector b.
   ParGridFunction x(fespace);
   x.ProjectCoefficient(solVec);

   // Integration rules used for the L^2 error computation. (In the serial
   // code the i == 2 branch was overwritten by the final else, so every
   // geometry other than i == 4 effectively used order 16; kept equivalent.)
   const IntegrationRule *irs[Geometry::NumGeom];
   for (int i = 0; i < Geometry::NumGeom; i++)
   {
      if (i == 4)
      {
         irs[i] = &(IntRules.Get(i, 10)); // Tet rule
      }
      else
      {
         irs[i] = &(IntRules.Get(i, 16)); // Everything else
      }
   }

   {
      real_t err0 = x.ComputeL2Error(solVec, irs);
      if (myid == 0)
      {
         cout << "\n Initial || E_h - E ||_{L^2} = " << err0 << '\n' << endl;
      }
   }

   // 11. Set up the parallel bilinear form corresponding to the EM diffusion
   //     operator curl muinv curl + sigma I.
   Coefficient *muinv = new ConstantCoefficient(1.0);
   Coefficient *sigma = new ConstantCoefficient(1.0);
   ParBilinearForm *a = new ParBilinearForm(fespace);
   if (pa) { a->SetAssemblyLevel(AssemblyLevel::PARTIAL); }
   a->AddDomainIntegrator(new CurlCurlIntegrator(*muinv));
   a->AddDomainIntegrator(new VectorFE_CurlMassIntegrator(*sigma));

   // 12. Assemble the parallel bilinear form and the corresponding linear
   //     system, applying parallel assembly, boundary condition elimination,
   //     conforming constraints for non-conforming meshes, static
   //     condensation, etc.
   if (static_cond) { a->EnableStaticCondensation(); }
   a->Assemble();

   OperatorPtr A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);

   // 13. Solve the system A X = B.
   if (pa)
   {
      // Jacobi preconditioning in partial assembly mode.
      OperatorJacobiSmoother M(*a, ess_tdof_list);
      CGSolver cg(MPI_COMM_WORLD);
      cg.SetOperator(*A);
      cg.SetPreconditioner(M);
      cg.SetRelTol(rel_tol);
      cg.SetAbsTol(0.0);
      cg.SetMaxIter(max_iter);
      cg.SetPrintLevel(1);
      cg.Mult(B, X);
   }
   else
   {
      HypreParMatrix &Ah = *A.As<HypreParMatrix>();
      if (myid == 0)
      {
         cout << "Size of linear system: " << Ah.GetGlobalNumRows() << endl;
      }

       Solver *prec = NULL;
       if (prec_type == 2)
       {
          HypreBoomerAMG *amg = new HypreBoomerAMG(Ah);
          amg->SetPrintLevel(0);
          prec = amg;
       }
       else
       {
          HypreSmoother *sm = new HypreSmoother;
          sm->SetType(prec_type == 1 ? HypreSmoother::l1Jacobi
                      : HypreSmoother::l1GS);
          sm->SetOperator(Ah);
          prec = sm;
       }

       CGSolver cg(MPI_COMM_WORLD);
       cg.SetOperator(Ah);
       cg.SetPreconditioner(*prec);
       cg.SetRelTol(rel_tol);
       cg.SetAbsTol(0.0);
       cg.SetMaxIter(max_iter);
       cg.SetPrintLevel(1);
       cg.Mult(B, X);

       delete prec;
   }

   // 14. Recover the parallel grid function corresponding to X. This is the
   //     local finite element solution on each processor.
   a->RecoverFEMSolution(X, *b, x);

   // 15. Compute and print the global L^2 norm of the error.
   {
      real_t err = x.ComputeL2Error(solVec, irs);
      if (myid == 0)
      {
         cout << "\n Final || E_h - E ||_{L^2} = " << err << '\n' << endl;
      }
   }

   // 16. Save the refined mesh and the solution in parallel. This output can
   //     be viewed later using GLVis: "glvis -np <np> -m refined -g sol_curl".
   {
      ostringstream mesh_name, sol_name;
      mesh_name << "refined." << setfill('0') << setw(6) << myid;
      sol_name << "sol_curl." << setfill('0') << setw(6) << myid;

      ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(8);
      pmesh->Print(mesh_ofs);

      ofstream sol_ofs(sol_name.str().c_str());
      sol_ofs.precision(8);
      x.Save(sol_ofs);
   }

   // 17. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock.precision(8);
      sol_sock << "solution\n" << *pmesh << x << flush;
   }

   // 18. Free the used memory.
   delete a;
   delete sigma;
   delete muinv;
   delete b;
   delete fespace;
   delete fec;
   delete pmesh;

   return 0;
}

void E_exact_vec(const Vector &x, Vector &E)
{
   int dim = x.Size();

   if (dim == 4)
   {
      E.SetSize(12);
      double sf = 0.5;

      double s0 = sin(M_PI*x(0)), s1 = sin(M_PI*x(1)), s2 = sin(M_PI*x(2)),
             s3 = sin(M_PI*x(3));
      double c0 = cos(M_PI*x(0)), c1 = cos(M_PI*x(1)), c2 = cos(M_PI*x(2)),
             c3 = cos(M_PI*x(3));

      // Current high order solution
      E(5) =  c0*c1*s2*s3;
      E(4) =  c0*s1*c2*s3;
      E(3) =  c0*s1*s2*c3;
      E(2) =  s0*c1*c2*s3;
      E(1) =  s0*c1*s2*c3;
      E(0) =  s0*s1*c2*c3;

      E(0) =  sf * E(0);
      E(1) =  sf * E(1);
      E(2) =  sf * E(2);
      E(3) =  sf * E(3);
      E(4) =  sf * E(4);
      E(5) =  sf * E(5);

      E(6)  = -1.0 * E(0);
      E(7)  = -1.0 * E(1);
      E(8)  = -1.0 * E(3);
      E(9)  = -1.0 * E(2);
      E(10) = -1.0 * E(4);
      E(11) = -1.0 * E(5);
   }
}

void E_exact(const Vector &x, DenseMatrix &E)
{
   int dim = x.Size();

   E.SetSize(dim*dim);

   if (dim == 4)
   {
      Vector vecE; E_exact_vec(x, vecE);

      E = 0.0;

      E(0,1) = vecE(0);
      E(0,2) = vecE(1);
      E(0,3) = vecE(2);
      E(1,2) = vecE(3);
      E(1,3) = vecE(4);
      E(2,3) = vecE(5);

      E(1,0) = -E(0,1);
      E(2,0) = -E(0,2);
      E(3,0) = -E(0,3);
      E(2,1) = -E(1,2);
      E(3,1) = -E(1,3);
      E(3,2) = -E(2,3);
   }
}

// f_exact = E + (Curl curl E), where P is the 4d permutation operator
void f_exact(const Vector &x, DenseMatrix &f)
{
   int dim = x.Size();

   f.SetSize(dim, dim);

   if (dim == 4)
   {
      f = 0.0;

      double s0 = sin(M_PI*x(0)), s1 = sin(M_PI*x(1)), s2 = sin(M_PI*x(2)),
             s3 = sin(M_PI*x(3));
      double c0 = cos(M_PI*x(0)), c1 = cos(M_PI*x(1)), c2 = cos(M_PI*x(2)),
             c3 = cos(M_PI*x(3));

      f(2,3) = (0.5 + 2.0 * M_PI*M_PI)*c0*c1*s2*s3;
      f(1,3) = (0.5 + 0.0 * M_PI*M_PI)*c0*s1*c2*s3;
      f(1,2) = (0.5 + 2.0 * M_PI*M_PI)*c0*s1*s2*c3;
      f(0,3) = (0.5 - 2.0 * M_PI*M_PI)*s0*c1*c2*s3;
      f(0,2) = (0.5 + 0.0 * M_PI*M_PI)*s0*c1*s2*c3;
      f(0,1) = (0.5 + 2.0 * M_PI*M_PI)*s0*s1*c2*c3;

      f(1,0) = -f(0,1);
      f(2,0) = -f(0,2);
      f(3,0) = -f(0,3);
      f(2,1) = -f(1,2);
      f(3,1) = -f(1,3);
      f(3,2) = -f(2,3);
   }
}
