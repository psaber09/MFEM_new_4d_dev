#include "mfem.hpp"
//#include "vector-dg-diffusion.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace mfem;

class RepeatedCoefficient : public VectorCoefficient
{
   Coefficient &coeff;
public:
   RepeatedCoefficient(int dim, Coefficient &coeff_)
      : VectorCoefficient(dim), coeff(coeff_)
   { }
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip)
   {
      V.SetSize(vdim);
      V = coeff.Eval(T, ip);
   }
};



void u_fn(const Vector &xvec, Vector &u_sol);
void f_fn(const Vector &xvec, Vector &f_vec);

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "../data/cube2d_2.MFEM";
   //const char *mesh_file = "../data/segment-4.MFEM";
   int ref_levels = 0;
   int order = 1;
   real_t sigma = -1.0;
   real_t kappa = -1.0;
   const char *device_config = "cpu";

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly, -1 for auto.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) >= 0.");
   args.AddOption(&sigma, "-s", "--sigma",
                  "One of the three DG penalty parameters, typically +1/-1."
                  " See the documentation of class DGDiffusionIntegrator.");
   args.AddOption(&kappa, "-k", "--kappa",
                  "One of the three DG penalty parameters, should be positive."
                  " Negative values are replaced with (order+1)^2.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.ParseCheck();

   if (kappa < 0)
   {
      kappa = (order+1)*(order+1);
   }

   Device device(device_config);
   device.Print();

   Mesh mesh(mesh_file);
   const int dim = mesh.Dimension();

   {
//      if (ref_levels < 0)
//      {
//         ref_levels = (int)floor(log(50000./mesh.GetNE())/log(2.)/dim);
//      }
      for (int l = 0; l < ref_levels; l++)
      {
         mesh.UniformRefinement();
      }
   }

   DG_FECollection fec(order, dim);
   FiniteElementSpace fespace(&mesh, &fec, dim);
   cout << "Number of unknowns: " << fespace.GetVSize() << endl;


//   FunctionCoefficient scalar_f_coeff(f_fn);
//   FunctionCoefficient scalar_u_coeff(u_fn);
//   RepeatedCoefficient f_coeff(dim, scalar_f_coeff);
//   RepeatedCoefficient u_coeff(dim, scalar_u_coeff);
   VectorFunctionCoefficient f_coeff(dim, f_fn);
   VectorFunctionCoefficient u_coeff(dim, u_fn);

   ConstantCoefficient one(1.0);
   ConstantCoefficient zero(5.0);
   RepeatedCoefficient zero_vec(dim, zero);

   LinearForm b(&fespace);
   b.AddDomainIntegrator(new VectorDomainLFIntegrator(f_coeff));
   b.AddBdrFaceIntegrator(
      new VectorDGDirichletLFIntegrator(u_coeff, one, sigma, kappa));
   b.Assemble();

   GridFunction x(&fespace);
   x = 0.0;

   BilinearForm a(&fespace);
   a.AddDomainIntegrator(new VectorDiffusionIntegrator(one));
   a.AddInteriorFaceIntegrator(new VectorDGDiffusionIntegrator(
                                  one, sigma, kappa, dim));
   a.AddBdrFaceIntegrator(new VectorDGDiffusionIntegrator(
                             one, sigma, kappa, dim));
   a.Assemble();
   a.Finalize();
   std::ofstream linear_sys("DG-Vector-Diffusion-Sys.txt");
   a.PrintMatlab(linear_sys);

   const SparseMatrix &A = a.SpMat();
#ifndef MFEM_USE_SUITESPARSE
   GSSmoother M(A);
   if (sigma == -1.0)
   {
      PCG(A, M, b, x, 1, 500, 1e-12, 0.0);
   }
   else
   {
      GMRES(A, M, b, x, 1, 500, 10, 1e-12, 0.0);
   }
#else
   UMFPackSolver umf_solver;
   umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
   umf_solver.SetOperator(A);
   umf_solver.Mult(b, x);
#endif

   ParaViewDataCollection pv("DGDiffusion", &mesh);
   pv.SetPrefixPath("ParaView");
   pv.SetHighOrderOutput(true);
   pv.SetLevelsOfDetail(order);
   pv.RegisterField("u", &x);
   pv.SetCycle(0);
   pv.SetTime(0.0);
   pv.Save();

   cout << "L2 error: " << x.ComputeL2Error(u_coeff) << '\n';

   return 0;
}

constexpr real_t pi = M_PI;
constexpr real_t pi2 = pi*pi;

//real_t u_fn(const Vector &xvec)
//{
//   int dim = xvec.Size();
//   real_t x = pi*xvec[0], y = pi*xvec[1];
//   if (dim == 2) { return sin(x)*sin(y); }
//   else { real_t z = pi*xvec[2]; return sin(x)*sin(y)*sin(z); }
//}
//
//real_t f_fn(const Vector &xvec)
//{
//   int dim = xvec.Size();
//   real_t x = pi*xvec[0], y = pi*xvec[1];
//
//   if (dim == 2)
//   {
//      return 2*pi2*sin(x)*sin(y);
//   }
//   else // dim == 3
//   {
//      real_t z = pi*xvec[2];
//      return 3*pi2*sin(x)*sin(y)*sin(z);
//   }
//}

void u_fn(const Vector &xvec, Vector &u_sol)
{
   int dim = xvec.Size();
   real_t x = pi*xvec[0], y = pi*xvec[1];
   if (dim == 1)
   {
       u_sol(0) =  sin(x);

   }
   if (dim == 2)
   {
       u_sol(0) =  sin(x)*sin(y);
       u_sol(1) =  cos(x)*cos(y);
   }
}

void f_fn(const Vector &xvec, Vector &f_vec)
{
   int dim = xvec.Size();
   real_t x = pi*xvec[0], y = pi*xvec[1];
   if (dim == 1)
   {
      f_vec(0) = pi2*sin(x);

   }
   if (dim == 2)
   {
      f_vec(0) = 2*pi2*sin(x)*sin(y);
      f_vec(1) = 2*pi2*cos(x)*cos(y);

   }
}
