//                       MFEM Example 1 - Parallel Version (H1)
//
// Compile with: make ex1p_H1
//
// Sample runs:  mpirun -np 4 ex1p_H1 -m ../data/star.mesh -ref 2 -pr 1
//               mpirun -np 4 ex1p_H1 -m ../data/fichera.mesh -o 2
//               mpirun -np 4 ex1p_H1 -m ../data/cube4d_96.MFEM -sr 1
//
// Description:  Parallel version of ex1_H1.cpp. Solves -Delta u = f with
//               homogeneous Dirichlet boundary conditions, where the exact
//               solution is the product of x_i(1-x_i) over all coordinates
//               (2D, 3D, or 4D) and f is the corresponding forcing term.
//               The structure follows ex1p_4d.cpp: parallel mesh, parallel
//               FE space, HypreBoomerAMG-preconditioned CG.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// Exact solution, E, and r.h.s., f. See below for implementation.
double E_exact(const Vector &);
double f_exact(const Vector &);
int dim;

int main(int argc, char *argv[])
{
   // 1. Initialize MPI.
   int num_procs, myid;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   bool verbose = (myid == 0);

   // 2. Parse command-line options.
   const char *mesh_file = "../data/star.mesh";
   int order = 1;
   bool static_cond = false;
   bool visualization = false;
   int ref_levels = 0;       // sequential (serial) refinement levels
   int par_ref_levels = 0;   // parallel refinement levels
   double tol = 1e-12;
   int max_iter = 5000;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                  " isoparametric space.");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&ref_levels, "-ref", "--ref-levels",
                  "Number of sequential (serial) uniform refinements.");
   args.AddOption(&par_ref_levels, "-pr", "--parrefinement",
                  "Number of parallel uniform refinements.");
   args.AddOption(&tol, "-tol", "--tol", "Relative tolerance for PCG.");
   args.AddOption(&max_iter, "-mi", "--max-iter", "Maximum PCG iterations.");
   args.Parse();
   if (!args.Good())
   {
      if (verbose) { args.PrintUsage(cout); }
      MPI_Finalize();
      return 1;
   }
   if (verbose) { args.PrintOptions(cout); }

   // 3. Read the serial mesh from the given mesh file.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   dim = mesh->Dimension();
   int sdim = mesh->SpaceDimension();
   MFEM_VERIFY(dim >= 2 && dim <= 4 && sdim == dim,
               "E_exact/f_exact are only defined for 2D, 3D and 4D meshes.");

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

   // 5. Define a parallel finite element space on the parallel mesh.
   FiniteElementCollection *fec;
   bool delete_fec;
   if (order > 0)
   {
      fec = new H1_FECollection(order, dim);
      delete_fec = true;
   }
   else if (pmesh->GetNodes())
   {
      fec = pmesh->GetNodes()->OwnFEC();
      delete_fec = false;
      if (verbose)
      {
         cout << "Using isoparametric FEs: " << fec->Name() << endl;
      }
   }
   else
   {
      if (dim == 4) { fec = new LinearFECollection; }
      else { fec = new H1_FECollection(order = 1, dim); }
      delete_fec = true;
   }
   ParFiniteElementSpace *fespace = new ParFiniteElementSpace(pmesh, fec);
   HYPRE_BigInt size = fespace->GlobalTrueVSize();
   if (verbose)
   {
      cout << "Number of finite element unknowns: " << size << endl;
   }
    
    // Consistency test: g -> R g -> P R g should reproduce g
    {
       ParGridFunction g(fespace), g2(fespace);
       Vector tv(fespace->GetTrueVSize());
       for (int d = 0; d < dim; d++)
       {
          FunctionCoefficient cd([d](const Vector &p) { return p(d); });
          g.ProjectCoefficient(cd);      // nodal interpolation of x_d
          g.GetTrueDofs(tv);             // keep owner values (R)
          g2.SetFromTrueDofs(tv);        // communicate to shared dofs (P)
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
             cout << "coord " << d << ": inconsistent dofs = " << nbad_tot << endl;
          }
       }
    }

   // 6. Determine the list of true (parallel conforming) essential boundary
   //    dofs: all boundary attributes are Dirichlet.
   Array<int> ess_tdof_list;
   if (pmesh->bdr_attributes.Size())
   {
      Array<int> ess_bdr(pmesh->bdr_attributes.Max());
      ess_bdr = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 7. Set up the parallel linear form b(.) with f_exact.
   FunctionCoefficient f(f_exact);
   ParLinearForm *b = new ParLinearForm(fespace);
   b->AddDomainIntegrator(new DomainLFIntegrator(f));
   b->Assemble();

   // 8. Define the solution vector x as a parallel grid function. Initialize
   //    it by projecting the exact solution (only boundary values are used
   //    when eliminating essential BCs).
   ParGridFunction x(fespace);
   ConstantCoefficient one(1.0);
   FunctionCoefficient E(E_exact);
   x.ProjectCoefficient(E);

   // High-order integration rules for the error computation (as in the
   // serial version): geometry 4 is the tetrahedron.
   const IntegrationRule *irs[Geometry::NumGeom];
   for (int i = 0; i < Geometry::NumGeom; i++)
   {
      if (i == 4) { irs[i] = &(IntRules.Get(i, 10)); }
      else        { irs[i] = &(IntRules.Get(i, 16)); }
   }

   // 8a. Compute and print the initial L^2 error.
   {
      double err0 = x.ComputeL2Error(E, irs);
      if (verbose)
      {
         cout << "\n Initial || E_h - E ||_{L^2} = " << err0 << '\n' << endl;
      }
   }

   // 9. Set up the parallel bilinear form a(.,.) for the Laplacian.
   ParBilinearForm *a = new ParBilinearForm(fespace);
   a->AddDomainIntegrator(new DiffusionIntegrator(one)); // -Delta u

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

   // 11. Solve A X = B with PCG preconditioned by BoomerAMG.
   HypreBoomerAMG *amg = new HypreBoomerAMG(A);
   amg->SetPrintLevel(0);

   CGSolver *pcg = new CGSolver(MPI_COMM_WORLD);
   pcg->SetOperator(A);
   pcg->SetRelTol(tol);
   pcg->SetMaxIter(max_iter);
   pcg->SetPrintLevel(1);
   pcg->SetPreconditioner(*amg);
   pcg->Mult(B, X);
   int iter = pcg->GetNumIterations();
   if (verbose) { cout << "PCG iterations: " << iter << endl; }
   delete pcg;
   delete amg;

   // 12. Recover the parallel grid function corresponding to X.
   a->RecoverFEMSolution(X, *b, x);

   // 13. Compute and print the L^2 norm of the error.
   {
      double err = x.ComputeL2Error(E, irs);
      if (verbose)
      {
         cout << "\n|| E_h - E ||_{L^2} = " << err << '\n' << endl;
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
   delete a;
   delete b;
   delete fespace;
   if (delete_fec) { delete fec; }
   delete pmesh;

   MPI_Finalize();
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
      E_out = x(0)*(1. - x(0)) * x(1)*(1. - x(1)) * x(2)*(1. - x(2))
              * x(3)*(1. - x(3));
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
      f_out = 2.0*((x(1)*(1. - x(1))*x(2)*(1. - x(2)))
                   + (x(0)*(1. - x(0))*x(2)*(1. - x(2)))
                   + (x(0)*(1. - x(0))*x(1)*(1. - x(1))));
   }
   if (dim == 4)
   {
      f_out = 2.0*((x(1)*(1. - x(1))*x(2)*(1. - x(2))*x(3)*(1. - x(3)))
                   + (x(0)*(1. - x(0))*x(2)*(1. - x(2))*x(3)*(1. - x(3)))
                   + (x(0)*(1. - x(0))*x(1)*(1. - x(1))*x(3)*(1. - x(3)))
                   + (x(0)*(1. - x(0))*x(1)*(1. - x(1))*x(2)*(1. - x(2))));
   }
   return f_out;
}
