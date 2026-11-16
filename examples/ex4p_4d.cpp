//                       MFEM Example 4 - Parallel Version (2D/3D/4D)
//
// Compile with: make ex4p_4d
//
// Sample runs:  mpirun -np 4 ex4p_4d -m ../data/cube4d_96.MFEM -ref 1
//               mpirun -np 4 ex4p_4d -m ../data/cube4d_96.MFEM -o 2 -pr 1
//               mpirun -np 4 ex4p_4d -m ../data/star.mesh
//               mpirun -np 4 ex4p_4d -m ../data/fichera.mesh -o 2
//
// Description:  Parallel H(div) diffusion problem
//                  -grad(alpha div F) + beta F = f,   F.n = <exact normal>
//               with a known exact solution F (2D, 3D or 4D, taken from the
//               serial ex4.cpp) and the corresponding right-hand side f.
//               Discretized with Raviart-Thomas elements, RT_FECollection
//               (order-1, dim). The structure follows ex1p_H1.cpp.
//
//               Preconditioning: hypre has an AMS (2D) and ADS (3D) solver for
//               H(div), but nothing for 4D. In 4D we default to an l1-Jacobi
//               smoothed CG. BoomerAMG can be selected with -prec 2 but is
//               not an effective preconditioner for H(div) problems.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// Exact solution, F, and r.h.s., f. See below for implementation.
void F_exact(const Vector &, Vector &);
void f_exact(const Vector &, Vector &);
real_t freq = 1.0, kappa;

int main(int argc, char *argv[])
{
   // 1. Initialize MPI.
   int num_procs, myid;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   bool verbose = (myid == 0);

   // 2. Parse command-line options.
   const char *mesh_file = "../data/cube4d_96.MFEM";
   int order = 1;
   bool set_bc = true;
   bool static_cond = false;
   bool visualization = false;
   int ref_levels = 0;       // sequential (serial) refinement levels
   int par_ref_levels = 0;   // parallel refinement levels
   double tol = 1e-12;
   int max_iter = 10000;
   int prec_type = 0;        // 0: auto, 1: l1-Jacobi, 2: BoomerAMG
   bool check_shared = false;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree); RT space of "
                  "degree order-1.");
   args.AddOption(&set_bc, "-bc", "--impose-bc", "-no-bc", "--dont-impose-bc",
                  "Impose or not essential boundary conditions.");
   args.AddOption(&freq, "-f", "--frequency", "Set the frequency for the exact"
                  " solution.");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
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
                  "Preconditioner: 0 = auto (AMS 2D, ADS 3D, Jacobi 4D), "
                  "1 = l1-Jacobi, 2 = BoomerAMG.");
   args.AddOption(&check_shared, "-check", "--check-shared", "-no-check",
                  "--no-check-shared",
                  "Check consistency of shared dofs (restriction followed "
                  "by prolongation).");
   args.Parse();
   if (!args.Good())
   {
      if (verbose) { args.PrintUsage(cout); }
      MPI_Finalize();
      return 1;
   }
   if (verbose) { args.PrintOptions(cout); }
   kappa = freq * M_PI;

   // 3. Read the serial mesh from the given mesh file.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();
   int sdim = mesh->SpaceDimension();
   MFEM_VERIFY(dim >= 2 && dim <= 4 && sdim == dim,
               "F_exact/f_exact are only defined for 2D, 3D and 4D meshes.");

   // 4. Refine the serial mesh, then partition it.
   if (verbose) { cout << "ref_levels = " << ref_levels << endl; }
   for (int l = 0; l < ref_levels; l++)
   {
      mesh->UniformRefinement();
   }

   if (verbose) { cout << "now we partition the mesh..." << endl << endl; }
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;

   for (int l = 0; l < par_ref_levels; l++)
   {
      pmesh->UniformRefinement();
   }
   pmesh->PrintInfo(cout);
   if (verbose) { cout << endl; }

   // 5. Define a parallel finite element space on the parallel mesh: the
   //    Raviart-Thomas elements of the specified order.
   FiniteElementCollection *fec = new RT_FECollection(order-1, dim);
   ParFiniteElementSpace *fespace = new ParFiniteElementSpace(pmesh, fec);
   HYPRE_BigInt size = fespace->GlobalTrueVSize();
   if (verbose)
   {
      cout << "Number of finite element unknowns: " << size << endl;
   }

   // 5a. Optional consistency test of the shared dofs: project the exact
   //     solution, restrict to true dofs and prolongate back. The result
   //     must reproduce the projected function on every processor.
   if (check_shared)
   {
      VectorFunctionCoefficient Fchk(sdim, F_exact);
      ParGridFunction g(fespace), g2(fespace);
      g.ProjectCoefficient(Fchk);
      Vector tv(fespace->GetTrueVSize());
      g.GetTrueDofs(tv);
      g2.SetFromTrueDofs(tv);
      g2 -= g;

      int nbad = 0;
      for (int i = 0; i < g2.Size(); i++)
      {
         if (fabs(g2(i)) > 1e-10) { nbad++; }
      }
      int nbad_tot;
      MPI_Allreduce(&nbad, &nbad_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if (verbose)
      {
         cout << "Shared dof check: inconsistent dofs = " << nbad_tot << endl;
      }
   }

   // 6. Determine the list of true (parallel conforming) essential boundary
   //    dofs: all boundary attributes are essential (Dirichlet).
   Array<int> ess_tdof_list;
   if (pmesh->bdr_attributes.Size())
   {
      Array<int> ess_bdr(pmesh->bdr_attributes.Max());
      ess_bdr = set_bc ? 1 : 0;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 7. Set up the parallel linear form b(.) with f_exact.
   VectorFunctionCoefficient f(sdim, f_exact);
   ParLinearForm *b = new ParLinearForm(fespace);
   b->AddDomainIntegrator(new VectorFEDomainLFIntegrator(f));
   b->Assemble();

   // 8. Define the solution vector x as a parallel grid function. Initialize
   //    it by projecting the exact solution (only boundary values are used
   //    when eliminating the essential BCs).
   ParGridFunction x(fespace);
   VectorFunctionCoefficient F(sdim, F_exact);
   x.ProjectCoefficient(F);

   // High-order integration rules for the error computation (as in the
   // serial version): geometry 4 is the tetrahedron.
   const IntegrationRule *irs[Geometry::NumGeom];
   for (int i = 0; i < Geometry::NumGeom; i++)
   {
      if (i == 4) { irs[i] = &(IntRules.Get(i, 10)); }
      else        { irs[i] = &(IntRules.Get(i, 16)); }
   }

   // 9. Set up the parallel bilinear form for the H(div) diffusion operator
   //    -grad(alpha div) + beta I: div-div plus the vector FE mass term.
   Coefficient *alpha = new ConstantCoefficient(1.0);
   Coefficient *beta  = new ConstantCoefficient(1.0);
   ParBilinearForm *a = new ParBilinearForm(fespace);
   a->AddDomainIntegrator(new DivDivIntegrator(*alpha));
   a->AddDomainIntegrator(new VectorFEMassIntegrator(*beta));

   // 10. Assemble the parallel bilinear form and the linear system.
   if (static_cond) { a->EnableStaticCondensation(); }
   a->Assemble();

   HypreParMatrix A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);

   if (verbose)
   {
      cout << "Size of linear system: " << A.GetGlobalNumRows() << endl;
   }

   // 11. Solve A X = B with PCG.
   Solver *prec = NULL;
   if (prec_type == 0)
   {
      ParFiniteElementSpace *prec_fespace =
         (a->StaticCondensationIsEnabled() ? a->SCParFESpace() : fespace);
      if (dim == 2)      { prec = new HypreAMS(A, prec_fespace); }
      else if (dim == 3) { prec = new HypreADS(A, prec_fespace); }
      else                { prec = new HypreSmoother(A, HypreSmoother::l1Jacobi); }
   }
   else if (prec_type == 1)
   {
      prec = new HypreSmoother(A, HypreSmoother::l1Jacobi);
   }
   else
   {
      HypreBoomerAMG *amg = new HypreBoomerAMG(A);
      amg->SetPrintLevel(0);
      prec = amg;
   }

   CGSolver *pcg = new CGSolver(MPI_COMM_WORLD);
   pcg->SetOperator(A);
   pcg->SetRelTol(tol);
   pcg->SetMaxIter(max_iter);
   pcg->SetPrintLevel(1);
   pcg->SetPreconditioner(*prec);
   pcg->Mult(B, X);
   if (verbose) { cout << "PCG iterations: " << pcg->GetNumIterations() << endl; }

   // 12. Recover the parallel grid function corresponding to X.
   a->RecoverFEMSolution(X, *b, x);

   // 13. Compute and print the L^2 norm of the error.
   {
      real_t err = x.ComputeL2Error(F, irs);
      if (verbose)
      {
         cout << "\n|| F_h - F ||_{L^2} = " << err << '\n' << endl;
      }
   }

   // 14. Save the refined mesh and the solution in parallel. View with
   //     "glvis -np <np> -m refined -g sol".
   {
      ostringstream mesh_name, sol_name;
      mesh_name << "refined." << setfill('0') << setw(6) << myid;
      sol_name << "sol." << setfill('0') << setw(6) << myid;

      ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(8);
      pmesh->Print(mesh_ofs);

      ofstream sol_ofs(sol_name.str().c_str());
      sol_ofs.precision(8);
      x.Save(sol_ofs);
   }

   // 15. Send the solution by socket to a GLVis server (2D/3D only).
   if (visualization && dim <= 3)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock.precision(8);
      sol_sock << "solution\n" << *pmesh << x << flush;
   }

   // 16. Free the used memory.
   delete pcg;
   delete prec;
   delete a;
   delete alpha;
   delete beta;
   delete b;
   delete fespace;
   delete fec;
   delete pmesh;

   MPI_Finalize();
   return 0;
}


// The exact solution
void F_exact(const Vector &p, Vector &F)
{
   int dim = p.Size();
   if (dim == 2)
   {
      real_t x = p(0);
      real_t y = p(1);

      F(0) = cos(kappa*x)*sin(kappa*y);
      F(1) = cos(kappa*y)*sin(kappa*x);
   }

   if (dim == 3)
   {
      real_t x = p(0);
      real_t y = p(1);
      real_t z = p(2);
      F(0) = cos(kappa*x)*sin(kappa*y)*sin(kappa*z);
      F(1) = cos(kappa*y)*sin(kappa*z)*sin(kappa*x);
      F(2) = cos(kappa*z)*sin(kappa*x)*sin(kappa*y);
   }

   if (dim == 4)
   {
      real_t x = p(0);
      real_t y = p(1);
      real_t z = p(2);
      real_t t = p(3);

      F(0) = cos(kappa*x)*sin(kappa*y)*sin(kappa*z)*sin(kappa*t);
      F(1) = cos(kappa*y)*sin(kappa*z)*sin(kappa*t)*sin(kappa*x);
      F(2) = cos(kappa*z)*sin(kappa*t)*sin(kappa*x)*sin(kappa*y);
      F(3) = cos(kappa*t)*sin(kappa*x)*sin(kappa*y)*sin(kappa*z);
   }
}

// The right hand side
void f_exact(const Vector &p, Vector &f)
{
   int dim = p.Size();

   if (dim == 2)
   {
      real_t x = p(0);
      real_t y = p(1);

      real_t temp = 1 + 2*kappa*kappa;

      f(0) = temp*cos(kappa*x)*sin(kappa*y);
      f(1) = temp*cos(kappa*y)*sin(kappa*x);
   }

   if (dim == 3)
   {
      real_t x = p(0);
      real_t y = p(1);
      real_t z = p(2);

      real_t temp = 1 + 3*kappa*kappa;
      f(0) = temp*cos(kappa*x)*sin(kappa*y)*sin(kappa*z);
      f(1) = temp*sin(kappa*x)*cos(kappa*y)*sin(kappa*z);
      f(2) = temp*sin(kappa*x)*sin(kappa*y)*cos(kappa*z);
   }

   if (dim == 4)
   {
      real_t x = p(0);
      real_t y = p(1);
      real_t z = p(2);
      real_t t = p(3);

      real_t temp = 1 + 4*kappa*kappa;
      f(0) = temp*cos(kappa*x)*sin(kappa*y)*sin(kappa*z)*sin(kappa*t);
      f(1) = temp*sin(kappa*x)*cos(kappa*y)*sin(kappa*z)*sin(kappa*t);
      f(2) = temp*sin(kappa*x)*sin(kappa*y)*cos(kappa*z)*sin(kappa*t);
      f(3) = temp*sin(kappa*x)*sin(kappa*y)*sin(kappa*z)*cos(kappa*t);
   }
}
