// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// Implementation of IntegrationRule(s) classes

// Acknowledgment: Some of the high-precision triangular and tetrahedral
// quadrature rules below were obtained from the Encyclopaedia of Cubature
// Formulas at http://nines.cs.kuleuven.be/research/ecf/ecf.html

#include "fem.hpp"
#include "../mesh/nurbs.hpp"
#include <cmath>

#ifdef MFEM_USE_MPFR
#include <mpfr.h>
#endif

using namespace std;

namespace mfem
{

IntegrationRule::IntegrationRule(IntegrationRule &irx, IntegrationRule &iry)
{
   int i, j, nx, ny;

   nx = irx.GetNPoints();
   ny = iry.GetNPoints();
   SetSize(nx * ny);
   SetPointIndices();
   Order = std::min(irx.GetOrder(), iry.GetOrder());

   for (j = 0; j < ny; j++)
   {
      IntegrationPoint &ipy = iry.IntPoint(j);
      for (i = 0; i < nx; i++)
      {
         IntegrationPoint &ipx = irx.IntPoint(i);
         IntegrationPoint &ip  = IntPoint(j*nx+i);

         ip.x = ipx.x;
         ip.y = ipy.x;
         ip.weight = ipx.weight * ipy.weight;
      }
   }
}

IntegrationRule::IntegrationRule(IntegrationRule &irx, IntegrationRule &iry,
                                 IntegrationRule &irz)
{
   const int nx = irx.GetNPoints();
   const int ny = iry.GetNPoints();
   const int nz = irz.GetNPoints();
   SetSize(nx*ny*nz);
   SetPointIndices();
   Order = std::min({irx.GetOrder(), iry.GetOrder(), irz.GetOrder()});

   for (int iz = 0; iz < nz; ++iz)
   {
      IntegrationPoint &ipz = irz.IntPoint(iz);
      for (int iy = 0; iy < ny; ++iy)
      {
         IntegrationPoint &ipy = iry.IntPoint(iy);
         for (int ix = 0; ix < nx; ++ix)
         {
            IntegrationPoint &ipx = irx.IntPoint(ix);
            IntegrationPoint &ip  = IntPoint(iz*nx*ny + iy*nx + ix);

            ip.x = ipx.x;
            ip.y = ipy.x;
            ip.z = ipz.x;
            ip.weight = ipx.weight*ipy.weight*ipz.weight;
         }
      }
   }
}

const Array<real_t> &IntegrationRule::GetWeights() const
{
   if (weights.Size() != GetNPoints())
   {
      weights.SetSize(GetNPoints());
      for (int i = 0; i < GetNPoints(); i++)
      {
         weights[i] = IntPoint(i).weight;
      }
   }
   return weights;
}

void IntegrationRule::SetPointIndices()
{
   for (int i = 0; i < Size(); i++)
   {
      IntPoint(i).index = i;
   }
}

void IntegrationRule::GrundmannMollerSimplexRule(int s, int n)
{
   // for pow on older compilers
   using std::pow;
   const int d = 2*s + 1;
   Vector fact(d + n + 1);
   Array<int> beta(n), sums(n);

   fact(0) = 1.;
   for (int i = 1; i < fact.Size(); i++)
   {
      fact(i) = fact(i - 1)*i;
   }

   // number of points is \binom{n + s + 1}{n + 1}
   int np = 1, f = 1;
   for (int i = 0; i <= n; i++)
   {
      np *= (s + i + 1), f *= (i + 1);
   }
   np /= f;
   SetSize(np);
   SetPointIndices();
   Order = 2*s + 1;

   int pt = 0;
   for (int i = 0; i <= s; i++)
   {
      real_t weight;

      weight = pow(2., -2*s)*pow(static_cast<real_t>(d + n - 2*i),
                                 d)/fact(i)/fact(d + n - i);
      if (i%2)
      {
         weight = -weight;
      }

      // loop over all beta : beta_0 + ... + beta_{n-1} <= s - i
      int k = s - i;
      beta = 0;
      sums = 0;
      while (true)
      {
         IntegrationPoint &ip = IntPoint(pt++);
         ip.weight = weight;
         ip.x = real_t(2*beta[0] + 1)/(d + n - 2*i);
         ip.y = real_t(2*beta[1] + 1)/(d + n - 2*i);
         if (n == 3)
         {
            ip.z = real_t(2*beta[2] + 1)/(d + n - 2*i);
         }

         int j = 0;
         while (sums[j] == k)
         {
            beta[j++] = 0;
            if (j == n)
            {
               goto done_beta;
            }
         }
         beta[j]++;
         sums[j]++;
         for (j--; j >= 0; j--)
         {
            sums[j] = sums[j+1];
         }
      }
   done_beta:
      ;
   }
}

IntegrationRule*
IntegrationRule::ApplyToKnotIntervals(KnotVector const& kv) const
{
   const int np = this->GetNPoints();
   const int ne = kv.GetNE();

   IntegrationRule *kvir = new IntegrationRule(ne * np);
   kvir->SetOrder(GetOrder());

   real_t x0 = kv[0];
   real_t x1 = x0;

   int id = 0;
   for (int e=0; e<ne; ++e)
   {
      x0 = x1;

      if (e == ne-1)
      {
         x1 = kv[kv.Size() - 1];
      }
      else
      {
         // Find the next unique knot
         while (id < kv.Size() - 1)
         {
            id++;
            if (kv[id] != x0)
            {
               x1 = kv[id];
               break;
            }
         }
      }

      const real_t s = x1 - x0;

      for (int j=0; j<this->GetNPoints(); ++j)
      {
         const real_t x = x0 + (s * (*this)[j].x);
         (*kvir)[(e * np) + j].Set1w(x, (*this)[j].weight);
      }
   }

   return kvir;
}

#ifdef MFEM_USE_MPFR

// Class for computing hi-precision (HP) quadrature in 1D
class HP_Quadrature1D
{
protected:
   mpfr_t pi, z, pp, p1, p2, p3, dz, w, rtol;

public:
   static const mpfr_rnd_t rnd = GMP_RNDN;
   static const int default_prec = 128;

   // prec = MPFR precision in bits
   HP_Quadrature1D(const int prec = default_prec)
   {
      mpfr_inits2(prec, pi, z, pp, p1, p2, p3, dz, w, rtol, (mpfr_ptr) 0);
      mpfr_const_pi(pi, rnd);
      mpfr_set_si_2exp(rtol, 1, -32, rnd); // 2^(-32) < 2.33e-10
   }

   // set rtol = 2^exponent
   // this is a tolerance for the last correction of x_i in Newton's algorithm;
   // this gives roughly rtol^2 accuracy for the final x_i.
   void SetRelTol(const int exponent = -32)
   {
      mpfr_set_si_2exp(rtol, 1, exponent, rnd);
   }

   // n - number of quadrature points
   // k - index of the point to compute, 0 <= k < n
   // see also: QuadratureFunctions1D::GaussLegendre
   void ComputeGaussLegendrePoint(const int n, const int k)
   {
      MFEM_ASSERT(n > 0 && 0 <= k && k < n, "invalid n = " << n
                  << " and/or k = " << k);

      int i = (k < (n+1)/2) ? k+1 : n-k;

      // Initial guess for the x-coordinate:
      // set z = cos(pi * (i - 0.25) / (n + 0.5)) =
      //       = sin(pi * ((n+1-2*i) / (2*n+1)))
      mpfr_set_si(z, n+1-2*i, rnd);
      mpfr_div_si(z, z, 2*n+1, rnd);
      mpfr_mul(z, z, pi, rnd);
      mpfr_sin(z, z, rnd);

      bool done = false;
      while (1)
      {
         mpfr_set_si(p2, 1, rnd);
         mpfr_set(p1, z, rnd);
         for (int j = 2; j <= n; j++)
         {
            mpfr_set(p3, p2, rnd);
            mpfr_set(p2, p1, rnd);
            // p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
            mpfr_mul_si(p1, z, 2*j-1, rnd);
            mpfr_mul_si(p3, p3, j-1, rnd);
            mpfr_fms(p1, p1, p2, p3, rnd);
            mpfr_div_si(p1, p1, j, rnd);
         }
         // p1 is Legendre polynomial

         // derivative of the Legendre polynomial:
         // pp = n * (z*p1-p2) / (z*z - 1);
         mpfr_fms(pp, z, p1, p2, rnd);
         mpfr_mul_si(pp, pp, n, rnd);
         mpfr_sqr(p2, z, rnd);
         mpfr_sub_si(p2, p2, 1, rnd);
         mpfr_div(pp, pp, p2, rnd);

         if (done) { break; }

         // set delta_z: dz = p1/pp;
         mpfr_div(dz, p1, pp, rnd);
         // compute absolute tolerance: atol = rtol*(1-z)
         mpfr_t &atol = w;
         mpfr_si_sub(atol, 1, z, rnd);
         mpfr_mul(atol, atol, rtol, rnd);
         if (mpfr_cmpabs(dz, atol) <= 0)
         {
            done = true;
            // continue the computation: get pp at the new point, then exit
         }
         // update z = z - dz
         mpfr_sub(z, z, dz, rnd);
      }

      // map z to (0,1): z = (1 - z)/2
      mpfr_si_sub(z, 1, z, rnd);
      mpfr_div_2si(z, z, 1, rnd);

      // weight: w = 1/(4*z*(1 - z)*pp*pp)
      mpfr_sqr(w, pp, rnd);
      mpfr_mul_2si(w, w, 2, rnd);
      mpfr_mul(w, w, z, rnd);
      mpfr_si_sub(p1, 1, z, rnd); // p1 = 1-z
      mpfr_mul(w, w, p1, rnd);
      mpfr_si_div(w, 1, w, rnd);

      if (k >= (n+1)/2) { mpfr_swap(z, p1); }
   }

   // n - number of quadrature points
   // k - index of the point to compute, 0 <= k < n
   // see also: QuadratureFunctions1D::GaussLobatto
   void ComputeGaussLobattoPoint(const int n, const int k)
   {
      MFEM_ASSERT(n > 1 && 0 <= k && k < n, "invalid n = " << n
                  << " and/or k = " << k);

      int i = (k < (n+1)/2) ? k : n-1-k;

      if (i == 0)
      {
         mpfr_set_si(z, 0, rnd);
         mpfr_set_si(p1, 1, rnd);
         mpfr_set_si(w, n*(n-1), rnd);
         mpfr_si_div(w, 1, w, rnd); // weight = 1/(n*(n-1))
         return;
      }
      // initial guess is the corresponding Chebyshev point, z:
      //    z = -cos(pi * i/(n-1)) = sin(pi * (2*i-n+1)/(2*n-2))
      mpfr_set_si(z, 2*i-n+1, rnd);
      mpfr_div_si(z, z, 2*(n-1), rnd);
      mpfr_mul(z, pi, z, rnd);
      mpfr_sin(z, z, rnd);
      bool done = false;
      for (int iter = 0 ; true ; ++iter)
      {
         // build Legendre polynomials, up to P_{n}(z)
         mpfr_set_si(p1, 1, rnd);
         mpfr_set(p2, z, rnd);

         for (int l = 1 ; l < (n-1) ; ++l)
         {
            // P_{l+1}(x) = [ (2*l+1)*x*P_l(x) - l*P_{l-1}(x) ]/(l+1)
            mpfr_mul_si(p1, p1, l, rnd);
            mpfr_mul_si(p3, z, 2*l+1, rnd);
            mpfr_fms(p3, p3, p2, p1, rnd);
            mpfr_div_si(p3, p3, l+1, rnd);

            mpfr_set(p1, p2, rnd);
            mpfr_set(p2, p3, rnd);
         }
         if (done) { break; }
         // compute dz = resid/deriv = (z*p2 - p1) / (n*p2);
         mpfr_fms(dz, z, p2, p1, rnd);
         mpfr_mul_si(p3, p2, n, rnd);
         mpfr_div(dz, dz, p3, rnd);
         // update: z = z - dz
         mpfr_sub(z, z, dz, rnd);
         // compute absolute tolerance: atol = rtol*(1 + z)
         mpfr_t &atol = w;
         mpfr_add_si(atol, z, 1, rnd);
         mpfr_mul(atol, atol, rtol, rnd);
         // check for convergence
         if (mpfr_cmpabs(dz, atol) <= 0)
         {
            done = true;
            // continue the computation: get p2 at the new point, then exit
         }
         // If the iteration does not converge fast, something is wrong.
         MFEM_VERIFY(iter < 8, "n = " << n << ", i = " << i
                     << ", dz = " << mpfr_get_d(dz, rnd));
      }
      // Map to the interval [0,1] and scale the weights
      mpfr_add_si(z, z, 1, rnd);
      mpfr_div_2si(z, z, 1, rnd);
      // set the symmetric point
      mpfr_si_sub(p1, 1, z, rnd);
      // w = 1/[ n*(n-1)*[P_{n-1}(z)]^2 ]
      mpfr_sqr(w, p2, rnd);
      mpfr_mul_si(w, w, n*(n-1), rnd);
      mpfr_si_div(w, 1, w, rnd);

      if (k >= (n+1)/2) { mpfr_swap(z, p1); }
   }

   real_t GetPoint() const { return mpfr_get_d(z, rnd); }
   real_t GetSymmPoint() const { return mpfr_get_d(p1, rnd); }
   real_t GetWeight() const { return mpfr_get_d(w, rnd); }

   const mpfr_t &GetHPPoint() const { return z; }
   const mpfr_t &GetHPSymmPoint() const { return p1; }
   const mpfr_t &GetHPWeight() const { return w; }

   ~HP_Quadrature1D()
   {
      mpfr_clears(pi, z, pp, p1, p2, p3, dz, w, rtol, (mpfr_ptr) 0);
      mpfr_free_cache();
   }
};

#endif // MFEM_USE_MPFR


void QuadratureFunctions1D::GaussLegendre(const int np, IntegrationRule* ir)
{
   ir->SetSize(np);
   ir->SetPointIndices();
   ir->SetOrder(2*np - 1);

   switch (np)
   {
      case 1:
         ir->IntPoint(0).Set1w(0.5, 1.0);
         return;
      case 2:
         ir->IntPoint(0).Set1w(0.21132486540518711775, 0.5);
         ir->IntPoint(1).Set1w(0.78867513459481288225, 0.5);
         return;
      case 3:
         ir->IntPoint(0).Set1w(0.11270166537925831148, 5./18.);
         ir->IntPoint(1).Set1w(0.5, 4./9.);
         ir->IntPoint(2).Set1w(0.88729833462074168852, 5./18.);
         return;
   }

   const int n = np;
   const int m = (n+1)/2;

#ifndef MFEM_USE_MPFR

   for (int i = 1; i <= m; i++)
   {
      real_t z = cos(M_PI * (i - 0.25) / (n + 0.5));
      real_t pp, p1, dz, xi = 0.;
      bool done = false;
      while (1)
      {
         real_t p2 = 1;
         p1 = z;
         for (int j = 2; j <= n; j++)
         {
            real_t p3 = p2;
            p2 = p1;
            p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
         }
         // p1 is Legendre polynomial

         pp = n * (z*p1-p2) / (z*z - 1);
         if (done) { break; }

         dz = p1/pp;
#ifdef MFEM_USE_SINGLE
         if (std::abs(dz) < 1e-7)
#elif defined MFEM_USE_DOUBLE
         if (std::abs(dz) < 1e-16)
#else
         MFEM_ABORT("Floating point type undefined");
         if (std::abs(dz) < 1e-16)
#endif
         {
            done = true;
            // map the new point (z-dz) to (0,1):
            xi = ((1 - z) + dz)/2; // (1 - (z - dz))/2 has bad round-off
            // continue the computation: get pp at the new point, then exit
         }
         // update: z = z - dz
         z -= dz;
      }

      ir->IntPoint(i-1).x = xi;
      ir->IntPoint(n-i).x = 1 - xi;
      ir->IntPoint(i-1).weight =
         ir->IntPoint(n-i).weight = 1./(4*xi*(1 - xi)*pp*pp);
   }

#else // MFEM_USE_MPFR is defined

   HP_Quadrature1D hp_quad;
   for (int i = 1; i <= m; i++)
   {
      hp_quad.ComputeGaussLegendrePoint(n, i-1);

      ir->IntPoint(i-1).x = hp_quad.GetPoint();
      ir->IntPoint(n-i).x = hp_quad.GetSymmPoint();
      ir->IntPoint(i-1).weight = ir->IntPoint(n-i).weight = hp_quad.GetWeight();
   }

#endif // MFEM_USE_MPFR

}

void QuadratureFunctions1D::GaussLobatto(const int np, IntegrationRule* ir)
{
   /* An np point Gauss-Lobatto quadrature has (np - 2) free abscissa the other
      (2) abscissa are the interval endpoints.

      The interior x_i are the zeros of P'_{np-1}(x). The weights of the
      interior points on the interval [-1,1] are:

      w_i = 2/(np*(np-1)*[P_{np-1}(x_i)]^2)

      The end point weights (on [-1,1]) are: w_{end} = 2/(np*(np-1)).

      The interior abscissa are found via a nonlinear solve, the initial guess
      for each point is the corresponding Chebyshev point.

      After we find all points on the interval [-1,1], we will map and scale the
      points and weights to the MFEM natural interval [0,1].

      References:
      [1] E. E. Lewis and W. F. Millier, "Computational Methods of Neutron
          Transport", Appendix A
      [2] the QUADRULE software by John Burkardt,
          https://people.sc.fsu.edu/~jburkardt/cpp_src/quadrule/quadrule.cpp
   */

   ir->SetSize(np);
   ir->SetPointIndices();
   if ( np == 1 )
   {
      ir->IntPoint(0).Set1w(0.5, 1.0);
      ir->SetOrder(1);
   }
   else
   {
      ir->SetOrder(2*np - 3);

#ifndef MFEM_USE_MPFR

      // endpoints and respective weights
      ir->IntPoint(0).x = 0.0;
      ir->IntPoint(np-1).x = 1.0;
      ir->IntPoint(0).weight = ir->IntPoint(np-1).weight = 1.0/(np*(np-1));

      // interior points and weights
      // use symmetry and compute just half of the points
      for (int i = 1 ; i <= (np-1)/2 ; ++i)
      {
         // initial guess is the corresponding Chebyshev point, x_i:
         //    x_i = -cos(\pi * (i / (np-1)))
         real_t x_i = std::sin(M_PI * ((real_t)(i)/(np-1) - 0.5));
         real_t z_i = 0., p_l;
         bool done = false;
         for (int iter = 0 ; true ; ++iter)
         {
            // build Legendre polynomials, up to P_{np}(x_i)
            real_t p_lm1 = 1.0;
            p_l = x_i;

            for (int l = 1 ; l < (np-1) ; ++l)
            {
               // The Legendre polynomials can be built by recursion:
               // x * P_l(x) = 1/(2*l+1)*[ (l+1)*P_{l+1}(x) + l*P_{l-1} ], i.e.
               // P_{l+1}(x) = [ (2*l+1)*x*P_l(x) - l*P_{l-1} ]/(l+1)
               real_t p_lp1 = ( (2*l + 1)*x_i*p_l - l*p_lm1)/(l + 1);

               p_lm1 = p_l;
               p_l = p_lp1;
            }
            if (done) { break; }
            // after this loop, p_l holds P_{np-1}(x_i)
            // resid = (x^2-1)*P'_{np-1}(x_i)
            // but use the recurrence relationship
            // (x^2 -1)P'_l(x) = l*[ x*P_l(x) - P_{l-1}(x) ]
            // thus, resid = (np-1) * (x_i*p_l - p_lm1)

            // The derivative of the residual is:
            // \frac{d}{d x} \left[ (x^2 -1)P'_l(x) ] \right] =
            // l * (l+1) * P_l(x), with l = np-1,
            // therefore, deriv = np * (np-1) * p_l;

            // compute dx = resid/deriv
            real_t dx = (x_i*p_l - p_lm1) / (np*p_l);
#ifdef MFEM_USE_SINGLE
            if (std::abs(dx) < 1e-7)
#elif defined MFEM_USE_DOUBLE
            if (std::abs(dx) < 1e-16)
#else
            MFEM_ABORT("Floating point type undefined");
            if (std::abs(dx) < 1e-16)
#endif
            {
               done = true;
               // Map the point to the interval [0,1]
               z_i = ((1.0 + x_i) - dx)/2;
               // continue the computation: get p_l at the new point, then exit
            }
            // If the iteration does not converge fast, something is wrong.
            MFEM_VERIFY(iter < 8, "np = " << np << ", i = " << i
                        << ", dx = " << dx);
            // update x_i:
            x_i -= dx;
         }
         // Map to the interval [0,1] and scale the weights
         IntegrationPoint &ip = ir->IntPoint(i);
         ip.x = z_i;
         // w_i = (2/[ n*(n-1)*[P_{n-1}(x_i)]^2 ]) / 2
         ip.weight = (real_t)(1.0 / (np*(np-1)*p_l*p_l));

         // set the symmetric point
         IntegrationPoint &symm_ip = ir->IntPoint(np-1-i);
         symm_ip.x = 1.0 - z_i;
         symm_ip.weight = ip.weight;
      }

#else // MFEM_USE_MPFR is defined

      HP_Quadrature1D hp_quad;
      // use symmetry and compute just half of the points
      for (int i = 0 ; i <= (np-1)/2 ; ++i)
      {
         hp_quad.ComputeGaussLobattoPoint(np, i);
         ir->IntPoint(i).x = hp_quad.GetPoint();
         ir->IntPoint(np-1-i).x = hp_quad.GetSymmPoint();
         ir->IntPoint(i).weight =
            ir->IntPoint(np-1-i).weight = hp_quad.GetWeight();
      }

#endif // MFEM_USE_MPFR

   }
}

void QuadratureFunctions1D::OpenUniform(const int np, IntegrationRule* ir)
{
   ir->SetSize(np);
   ir->SetPointIndices();
   ir->SetOrder(np - 1 + np%2);

   // The Newton-Cotes quadrature is based on weights that integrate exactly the
   // interpolatory polynomial through the equally spaced quadrature points.
   for (int i = 0; i < np ; ++i)
   {
      ir->IntPoint(i).x = real_t(i+1) / real_t(np + 1);
   }

   CalculateUniformWeights(ir, Quadrature1D::OpenUniform);
}

void QuadratureFunctions1D::ClosedUniform(const int np,
                                          IntegrationRule* ir)
{
   ir->SetSize(np);
   ir->SetPointIndices();
   ir->SetOrder(np - 1 + np%2);
   if ( np == 1 ) // allow this case as "closed"
   {
      ir->IntPoint(0).Set1w(0.5, 1.0);
      return;
   }

   for (int i = 0; i < np ; ++i)
   {
      ir->IntPoint(i).x = real_t(i) / (np-1);
   }

   CalculateUniformWeights(ir, Quadrature1D::ClosedUniform);
}

void QuadratureFunctions1D::OpenHalfUniform(const int np, IntegrationRule* ir)
{
   ir->SetSize(np);
   ir->SetPointIndices();
   ir->SetOrder(np - 1 + np%2);

   // Open half points: the centers of np uniform intervals
   for (int i = 0; i < np ; ++i)
   {
      ir->IntPoint(i).x = real_t(2*i+1) / (2*np);
   }

   CalculateUniformWeights(ir, Quadrature1D::OpenHalfUniform);
}

void QuadratureFunctions1D::ClosedGL(const int np, IntegrationRule* ir)
{
   ir->SetSize(np);
   ir->SetPointIndices();
   ir->IntPoint(0).x = 0.0;
   ir->IntPoint(np-1).x = 1.0;
   ir->SetOrder(np - 1 + np%2); // Is this the correct order?

   if ( np > 2 )
   {
      IntegrationRule gl_ir;
      GaussLegendre(np-1, &gl_ir);

      for (int i = 1; i < np-1; ++i)
      {
         ir->IntPoint(i).x = (gl_ir.IntPoint(i-1).x + gl_ir.IntPoint(i).x)/2;
      }
   }

   CalculateUniformWeights(ir, Quadrature1D::ClosedGL);
}

void QuadratureFunctions1D::GivePolyPoints(const int np, real_t *pts,
                                           const int type)
{
   IntegrationRule ir(np);

   switch (type)
   {
      case Quadrature1D::GaussLegendre:
      {
         GaussLegendre(np,&ir);
         break;
      }
      case Quadrature1D::GaussLobatto:
      {
         GaussLobatto(np, &ir);
         break;
      }
      case Quadrature1D::OpenUniform:
      {
         OpenUniform(np,&ir);
         break;
      }
      case Quadrature1D::ClosedUniform:
      {
         ClosedUniform(np,&ir);
         break;
      }
      case Quadrature1D::OpenHalfUniform:
      {
         OpenHalfUniform(np, &ir);
         break;
      }
      case Quadrature1D::ClosedGL:
      {
         ClosedGL(np, &ir);
         break;
      }
      case Quadrature1D::Invalid:
      {
         MFEM_ABORT("Asking for an unknown type of 1D Quadrature points, "
                    "type = " << type);
      }
   }

   for (int i = 0 ; i < np ; ++i)
   {
      pts[i] = ir.IntPoint(i).x;
   }
}

void QuadratureFunctions1D::CalculateUniformWeights(IntegrationRule *ir,
                                                    const int type)
{
   /* The Lagrange polynomials are:
           p_i = \prod_{j \neq i} {\frac{x - x_j }{x_i - x_j}}

      The weight associated with each abscissa is the integral of p_i over
      [0,1]. To calculate the integral of p_i, we use a Gauss-Legendre
      quadrature rule. This approach does not suffer from bad round-off/
      cancellation errors for large number of points.
   */
   const int n = ir->Size();
   switch (n)
   {
      case 1:
         ir->IntPoint(0).weight = 1.;
         return;
      case 2:
         ir->IntPoint(0).weight = .5;
         ir->IntPoint(1).weight = .5;
         return;
   }

#ifndef MFEM_USE_MPFR

   // This algorithm should work for any set of points, not just uniform
   const IntegrationRule &glob_ir = IntRules.Get(Geometry::SEGMENT, n-1);
   const int m = glob_ir.GetNPoints();
   Vector xv(n);
   for (int j = 0; j < n; j++)
   {
      xv(j) = ir->IntPoint(j).x;
   }
   Poly_1D::Basis basis(n-1, xv.GetData()); // nodal basis, with nodes at 'xv'
   Vector w(n);
   // Integrate all nodal basis functions using 'glob_ir':
   w = 0.0;
   for (int i = 0; i < m; i++)
   {
      const IntegrationPoint &ip = glob_ir.IntPoint(i);
      basis.Eval(ip.x, xv);
      w.Add(ip.weight, xv); // w += ip.weight * xv
   }
   for (int j = 0; j < n; j++)
   {
      ir->IntPoint(j).weight = w(j);
   }

#else // MFEM_USE_MPFR is defined

   static const mpfr_rnd_t rnd = HP_Quadrature1D::rnd;
   HP_Quadrature1D hp_quad;
   mpfr_t l, lk, w0, wi, tmp, *weights;
   mpfr_inits2(hp_quad.default_prec, l, lk, w0, wi, tmp, (mpfr_ptr) 0);
   weights = new mpfr_t[n];
   for (int i = 0; i < n; i++)
   {
      mpfr_init2(weights[i], hp_quad.default_prec);
      mpfr_set_si(weights[i], 0, rnd);
   }
   hp_quad.SetRelTol(-48); // rtol = 2^(-48) ~ 3.5e-15
   const int p = n-1;
   const int m = p/2+1; // number of points for Gauss-Legendre quadrature
   int hinv = 0, ihoffset = 0; // x_i = (i+ihoffset/2)/hinv
   switch (type)
   {
      case Quadrature1D::ClosedUniform:
         // x_i = i/p, i=0,...,p
         hinv = p;
         ihoffset = 0;
         break;
      case Quadrature1D::OpenUniform:
         // x_i = (i+1)/(p+2), i=0,...,p
         hinv = p+2;
         ihoffset = 2;
         break;
      case Quadrature1D::OpenHalfUniform:
         // x_i = (i+1/2)/(p+1), i=0,...,p
         hinv = p+1;
         ihoffset = 1;
         break;
      case Quadrature1D::GaussLegendre:
      case Quadrature1D::GaussLobatto:
      case Quadrature1D::ClosedGL:
      case Quadrature1D::Invalid:
         MFEM_ABORT("invalid Quadrature1D type: " << type);
   }
   // set w0 = (-1)^p*(p!)/(hinv^p)
   mpfr_fac_ui(w0, p, rnd);
   mpfr_ui_pow_ui(tmp, hinv, p, rnd);
   mpfr_div(w0, w0, tmp, rnd);
   if (p%2) { mpfr_neg(w0, w0, rnd); }

   for (int j = 0; j < m; j++)
   {
      hp_quad.ComputeGaussLegendrePoint(m, j);

      // Compute l = \prod_{i=0}^p (x-x_i) and lk = l/(x-x_k), where
      // x = hp_quad.GetHPPoint(), x_i = (i+ihoffset/2)/hinv, and x_k is the
      // node closest to x, i.e. k = min(max(round(x*hinv-ihoffset/2),0),p)
      mpfr_mul_si(tmp, hp_quad.GetHPPoint(), hinv, rnd);
      mpfr_sub_d(tmp, tmp, 0.5*ihoffset, rnd);
      mpfr_round(tmp, tmp);
      int k = min(max((int)mpfr_get_si(tmp, rnd), 0), p);
      mpfr_set_si(lk, 1, rnd);
      for (int i = 0; i <= p; i++)
      {
         mpfr_set_si(tmp, 2*i+ihoffset, rnd);
         mpfr_div_si(tmp, tmp, 2*hinv, rnd);
         mpfr_sub(tmp, hp_quad.GetHPPoint(), tmp, rnd);
         if (i != k)
         {
            mpfr_mul(lk, lk, tmp, rnd);
         }
         else
         {
            mpfr_set(l, tmp, rnd);
         }
      }
      mpfr_mul(l, l, lk, rnd);
      mpfr_set(wi, w0, rnd);
      for (int i = 0; true; i++)
      {
         if (i != k)
         {
            // tmp = l/(wi*(x - x_i))
            mpfr_set_si(tmp, 2*i+ihoffset, rnd);
            mpfr_div_si(tmp, tmp, 2*hinv, rnd);
            mpfr_sub(tmp, hp_quad.GetHPPoint(), tmp, rnd);
            mpfr_mul(tmp, tmp, wi, rnd);
            mpfr_div(tmp, l, tmp, rnd);
         }
         else
         {
            // tmp = lk/wi
            mpfr_div(tmp, lk, wi, rnd);
         }
         // weights[i] += hp_quad.weight*tmp
         mpfr_mul(tmp, tmp, hp_quad.GetHPWeight(), rnd);
         mpfr_add(weights[i], weights[i], tmp, rnd);

         if (i == p) { break; }

         // update wi *= (i+1)/(i-p)
         mpfr_mul_si(wi, wi, i+1, rnd);
         mpfr_div_si(wi, wi, i-p, rnd);
      }
   }
   for (int i = 0; i < n; i++)
   {
      ir->IntPoint(i).weight = mpfr_get_d(weights[i], rnd);
      mpfr_clear(weights[i]);
   }
   delete [] weights;
   mpfr_clears(l, lk, w0, wi, tmp, (mpfr_ptr) 0);

#endif // MFEM_USE_MPFR

}


int Quadrature1D::CheckClosed(int type)
{
   switch (type)
   {
      case GaussLobatto:
      case ClosedUniform:
      case ClosedGL:
         return type;
      default:
         return Invalid;
   }
}

int Quadrature1D::CheckOpen(int type)
{
   switch (type)
   {
      case GaussLegendre:
      case GaussLobatto:
      case OpenUniform:
      case ClosedUniform:
      case OpenHalfUniform:
         return type; // all types can work as open
      default:
         return Invalid;
   }
}


IntegrationRules IntRules(0, Quadrature1D::GaussLegendre);

IntegrationRules RefinedIntRules(1, Quadrature1D::GaussLegendre);

IntegrationRules::IntegrationRules(int ref, int type)
   : quad_type(type)
{
   refined = ref;

   if (refined < 0) { own_rules = 0; return; }

   own_rules = 1;

   const MemoryType h_mt = MemoryType::HOST;
   PointIntRules.SetSize(2, h_mt);
   PointIntRules = NULL;

   SegmentIntRules.SetSize(32, h_mt);
   SegmentIntRules = NULL;

   // TriangleIntegrationRule() assumes that this size is >= 26
   TriangleIntRules.SetSize(32, h_mt);
   TriangleIntRules = NULL;

   SquareIntRules.SetSize(32, h_mt);
   SquareIntRules = NULL;

   // TetrahedronIntegrationRule() assumes that this size is >= 10
   TetrahedronIntRules.SetSize(32, h_mt);
   TetrahedronIntRules = NULL;

   PyramidIntRules.SetSize(32, h_mt);
   PyramidIntRules = NULL;

   PrismIntRules.SetSize(32, h_mt);
   PrismIntRules = NULL;

   CubeIntRules.SetSize(32, h_mt);
   CubeIntRules = NULL;
    
   PentatopeIntRules.SetSize(32);
   PentatopeIntRules = NULL;


#if defined(MFEM_THREAD_SAFE) && defined(MFEM_USE_OPENMP)
   IntRuleLocks.SetSize(Geometry::NUM_GEOMETRIES, h_mt);
   for (int i = 0; i < Geometry::NUM_GEOMETRIES; i++)
   {
      omp_init_lock(&IntRuleLocks[i]);
   }
#endif
}

const IntegrationRule &IntegrationRules::Get(int GeomType, int Order)
{
   Array<IntegrationRule *> *ir_array = NULL;

   switch (GeomType)
   {
      case Geometry::POINT:       ir_array = &PointIntRules; Order = 0; break;
      case Geometry::SEGMENT:     ir_array = &SegmentIntRules; break;
      case Geometry::TRIANGLE:    ir_array = &TriangleIntRules; break;
      case Geometry::SQUARE:      ir_array = &SquareIntRules; break;
      case Geometry::TETRAHEDRON: ir_array = &TetrahedronIntRules; break;
      case Geometry::CUBE:        ir_array = &CubeIntRules; break;
      case Geometry::PRISM:       ir_array = &PrismIntRules; break;
      case Geometry::PYRAMID:     ir_array = &PyramidIntRules; break;
      case Geometry::PENTATOPE:   ir_array = &PentatopeIntRules; break;
      case Geometry::INVALID:
      case Geometry::NUM_GEOMETRIES:
         MFEM_ABORT("Unknown type of reference element!");
   }

   if (Order < 0)
   {
      Order = 0;
   }

#if defined(MFEM_THREAD_SAFE) && defined(MFEM_USE_OPENMP)
   omp_set_lock(&IntRuleLocks[GeomType]);
#endif

   if (!HaveIntRule(*ir_array, Order))
   {
      IntegrationRule *ir = GenerateIntegrationRule(GeomType, Order);
#ifdef MFEM_DEBUG
      int RealOrder = Order;
      while (RealOrder+1 < ir_array->Size() && (*ir_array)[RealOrder+1] == ir)
      {
         RealOrder++;
      }
      MFEM_VERIFY(RealOrder == ir->GetOrder(), "internal error");
#else
      MFEM_CONTRACT_VAR(ir);
#endif
   }

#if defined(MFEM_THREAD_SAFE) && defined(MFEM_USE_OPENMP)
   omp_unset_lock(&IntRuleLocks[GeomType]);
#endif

   return *(*ir_array)[Order];
}

void IntegrationRules::Set(int GeomType, int Order, IntegrationRule &IntRule)
{
   Array<IntegrationRule *> *ir_array = NULL;

   switch (GeomType)
   {
      case Geometry::POINT:       ir_array = &PointIntRules; break;
      case Geometry::SEGMENT:     ir_array = &SegmentIntRules; break;
      case Geometry::TRIANGLE:    ir_array = &TriangleIntRules; break;
      case Geometry::SQUARE:      ir_array = &SquareIntRules; break;
      case Geometry::TETRAHEDRON: ir_array = &TetrahedronIntRules; break;
      case Geometry::CUBE:        ir_array = &CubeIntRules; break;
      case Geometry::PRISM:       ir_array = &PrismIntRules; break;
      case Geometry::PYRAMID:     ir_array = &PyramidIntRules; break;
      case Geometry::PENTATOPE:   ir_array = &PentatopeIntRules; break;
      case Geometry::INVALID:
      case Geometry::NUM_GEOMETRIES:
         MFEM_ABORT("Unknown type of reference element!");
   }

#if defined(MFEM_THREAD_SAFE) && defined(MFEM_USE_OPENMP)
   omp_set_lock(&IntRuleLocks[GeomType]);
#endif

   if (HaveIntRule(*ir_array, Order))
   {
      MFEM_ABORT("Overwriting set rules is not supported!");
   }

   AllocIntRule(*ir_array, Order);

   (*ir_array)[Order] = &IntRule;

#if defined(MFEM_THREAD_SAFE) && defined(MFEM_USE_OPENMP)
   omp_unset_lock(&IntRuleLocks[GeomType]);
#endif
}

void IntegrationRules::DeleteIntRuleArray(
   Array<IntegrationRule *> &ir_array) const
{
   // Many of the intrules have multiple contiguous copies in the ir_array
   // so we have to be careful to not delete them twice.
   IntegrationRule *ir = NULL;
   for (int i = 0; i < ir_array.Size(); i++)
   {
      if (ir_array[i] != NULL && ir_array[i] != ir)
      {
         ir = ir_array[i];
         delete ir;
      }
   }
}

IntegrationRules::~IntegrationRules()
{
#if defined(MFEM_THREAD_SAFE) && defined(MFEM_USE_OPENMP)
   for (int i = 0; i < Geometry::NUM_GEOMETRIES; i++)
   {
      omp_destroy_lock(&IntRuleLocks[i]);
   }
#endif

   if (!own_rules) { return; }

   DeleteIntRuleArray(PointIntRules);
   DeleteIntRuleArray(SegmentIntRules);
   DeleteIntRuleArray(TriangleIntRules);
   DeleteIntRuleArray(SquareIntRules);
   DeleteIntRuleArray(TetrahedronIntRules);
   DeleteIntRuleArray(CubeIntRules);
   DeleteIntRuleArray(PrismIntRules);
   DeleteIntRuleArray(PyramidIntRules);
   DeleteIntRuleArray(PentatopeIntRules);
}


IntegrationRule *IntegrationRules::GenerateIntegrationRule(int GeomType,
                                                           int Order)
{
   switch (GeomType)
   {
      case Geometry::POINT:
         return PointIntegrationRule(Order);
      case Geometry::SEGMENT:
         return SegmentIntegrationRule(Order);
      case Geometry::TRIANGLE:
         return TriangleIntegrationRule(Order);
      case Geometry::SQUARE:
         return SquareIntegrationRule(Order);
      case Geometry::TETRAHEDRON:
         return TetrahedronIntegrationRule(Order);
      case Geometry::CUBE:
         return CubeIntegrationRule(Order);
      case Geometry::PRISM:
         return PrismIntegrationRule(Order);
      case Geometry::PYRAMID:
         return PyramidIntegrationRule(Order);
      case Geometry::PENTATOPE:
         return PentatopeIntegrationRule(Order);
      case Geometry::INVALID:
      case Geometry::NUM_GEOMETRIES:
         MFEM_ABORT("Unknown type of reference element!");
   }
   return NULL;
}


// Integration rules for a point
IntegrationRule *IntegrationRules::PointIntegrationRule(int Order)
{
   if (Order > 1)
   {
      MFEM_ABORT("Point Integration Rule of Order > 1 not defined");
      return NULL;
   }

   IntegrationRule *ir = new IntegrationRule(1);
   ir->IntPoint(0).x = .0;
   ir->IntPoint(0).weight = 1.;
   ir->SetOrder(1);

   PointIntRules[1] = PointIntRules[0] = ir;

   return ir;
}

// Integration rules for line segment [0,1]
IntegrationRule *IntegrationRules::SegmentIntegrationRule(int Order)
{
   int RealOrder = GetSegmentRealOrder(Order); // RealOrder >= Order
   // Order is one of {RealOrder-1,RealOrder}
   AllocIntRule(SegmentIntRules, RealOrder);

   IntegrationRule *ir = new IntegrationRule;

   int n = 0;
   // n is the number of points to achieve the exact integral of a
   // degree Order polynomial
   switch (quad_type)
   {
      case Quadrature1D::GaussLegendre:
      {
         // Gauss-Legendre is exact for 2*n-1
         n = Order/2 + 1;
         QuadratureFunctions1D::GaussLegendre(n, ir);
         break;
      }
      case Quadrature1D::GaussLobatto:
      {
         // Gauss-Lobatto is exact for 2*n-3
         n = Order/2 + 2;
         QuadratureFunctions1D::GaussLobatto(n, ir);
         break;
      }
      case Quadrature1D::OpenUniform:
      {
         // Open Newton Cotes is exact for n-(n+1)%2 = n-1+n%2
         n = Order | 1; // n is always odd
         QuadratureFunctions1D::OpenUniform(n, ir);
         break;
      }
      case Quadrature1D::ClosedUniform:
      {
         // Closed Newton Cotes is exact for n-(n+1)%2 = n-1+n%2
         n = Order | 1; // n is always odd
         QuadratureFunctions1D::ClosedUniform(n, ir);
         break;
      }
      case Quadrature1D::OpenHalfUniform:
      {
         // Open half Newton Cotes is exact for n-(n+1)%2 = n-1+n%2
         n = Order | 1; // n is always odd
         QuadratureFunctions1D::OpenHalfUniform(n, ir);
         break;
      }
      case Quadrature1D::Invalid:
      {
         MFEM_ABORT("unknown Quadrature1D type: " << quad_type);
      }
   }
   if (refined)
   {
      // Effectively passing memory management to SegmentIntegrationRules
      IntegrationRule *refined_ir = new IntegrationRule(2*n);
      refined_ir->SetOrder(ir->GetOrder());
      for (int j = 0; j < n; j++)
      {
         refined_ir->IntPoint(j).x = ir->IntPoint(j).x/2.0;
         refined_ir->IntPoint(j).weight = ir->IntPoint(j).weight/2.0;
         refined_ir->IntPoint(j+n).x = 0.5 + ir->IntPoint(j).x/2.0;
         refined_ir->IntPoint(j+n).weight = ir->IntPoint(j).weight/2.0;
      }
      delete ir;
      ir = refined_ir;
   }
   SegmentIntRules[RealOrder-1] = SegmentIntRules[RealOrder] = ir;
   return ir;
}

// Integration rules for reference triangle {[0,0],[1,0],[0,1]}
IntegrationRule *IntegrationRules::TriangleIntegrationRule(int Order)
{
   IntegrationRule *ir = NULL;
   // Note: Set TriangleIntRules[*] to ir only *after* ir is fully constructed.
   // This is needed in multithreaded environment.

   // assuming that orders <= 25 are pre-allocated
   switch (Order)
   {
      case 0:  // 1 point - degree 1
      case 1:
         ir = new IntegrationRule(1);
         ir->AddTriMidPoint(0, 0.5);
         ir->SetOrder(1);
         TriangleIntRules[0] = TriangleIntRules[1] = ir;
         return ir;

      case 2:  // 3 point - 2 degree
         ir = new IntegrationRule(3);
         ir->AddTriPoints3(0, 1./6., 1./6.);
         ir->SetOrder(2);
         TriangleIntRules[2] = ir;
         // interior points
         return ir;

      case 3:  // 4 point - 3 degree (has one negative weight)
         ir = new IntegrationRule(4);
         ir->AddTriMidPoint(0, -0.28125); // -9./32.
         ir->AddTriPoints3(1, 0.2, 25./96.);
         ir->SetOrder(3);
         TriangleIntRules[3] = ir;
         return ir;

      case 4:  // 6 point - 4 degree
         ir = new IntegrationRule(6);
         ir->AddTriPoints3(0, 0.091576213509770743460, 0.054975871827660933819);
         ir->AddTriPoints3(3, 0.44594849091596488632, 0.11169079483900573285);
         ir->SetOrder(4);
         TriangleIntRules[4] = ir;
         return ir;

      case 5:  // 7 point - 5 degree
         ir = new IntegrationRule(7);
         ir->AddTriMidPoint(0, 0.1125);
         ir->AddTriPoints3(1, 0.10128650732345633880, 0.062969590272413576298);
         ir->AddTriPoints3(4, 0.47014206410511508977, 0.066197076394253090369);
         ir->SetOrder(5);
         TriangleIntRules[5] = ir;
         return ir;

      case 6:  // 12 point - 6 degree
         ir = new IntegrationRule(12);
         ir->AddTriPoints3(0, 0.063089014491502228340, 0.025422453185103408460);
         ir->AddTriPoints3(3, 0.24928674517091042129, 0.058393137863189683013);
         ir->AddTriPoints6(6, 0.053145049844816947353, 0.31035245103378440542,
                           0.041425537809186787597);
         ir->SetOrder(6);
         TriangleIntRules[6] = ir;
         return ir;

      case 7:  // 12 point - degree 7
         ir = new IntegrationRule(12);
         ir->AddTriPoints3R(0, 0.062382265094402118174, 0.067517867073916085443,
                            0.026517028157436251429);
         ir->AddTriPoints3R(3, 0.055225456656926611737, 0.32150249385198182267,
                            0.043881408714446055037);
         // slightly better with explicit 3rd area coordinate
         ir->AddTriPoints3R(6, 0.034324302945097146470, 0.66094919618673565761,
                            0.30472650086816719592, 0.028775042784981585738);
         ir->AddTriPoints3R(9, 0.51584233435359177926, 0.27771616697639178257,
                            0.20644149867001643817, 0.067493187009802774463);
         ir->SetOrder(7);
         TriangleIntRules[7] = ir;
         return ir;

      case 8:  // 16 point - 8 degree
         ir = new IntegrationRule(16);
         ir->AddTriMidPoint(0, 0.0721578038388935841255455552445323);
         ir->AddTriPoints3(1, 0.170569307751760206622293501491464,
                           0.0516086852673591251408957751460645);
         ir->AddTriPoints3(4, 0.0505472283170309754584235505965989,
                           0.0162292488115990401554629641708902);
         ir->AddTriPoints3(7, 0.459292588292723156028815514494169,
                           0.0475458171336423123969480521942921);
         ir->AddTriPoints6(10, 0.008394777409957605337213834539296,
                           0.263112829634638113421785786284643,
                           0.0136151570872174971324223450369544);
         ir->SetOrder(8);
         TriangleIntRules[8] = ir;
         return ir;

      case 9:  // 19 point - 9 degree
         ir = new IntegrationRule(19);
         ir->AddTriMidPoint(0, 0.0485678981413994169096209912536443);
         ir->AddTriPoints3b(1, 0.020634961602524744433,
                            0.0156673501135695352684274156436046);
         ir->AddTriPoints3b(4, 0.12582081701412672546,
                            0.0389137705023871396583696781497019);
         ir->AddTriPoints3(7, 0.188203535619032730240961280467335,
                           0.0398238694636051265164458871320226);
         ir->AddTriPoints3(10, 0.0447295133944527098651065899662763,
                           0.0127888378293490156308393992794999);
         ir->AddTriPoints6(13, 0.0368384120547362836348175987833851,
                           0.2219629891607656956751025276931919,
                           0.0216417696886446886446886446886446);
         ir->SetOrder(9);
         TriangleIntRules[9] = ir;
         return ir;

      case 10:  // 25 point - 10 degree
         ir = new IntegrationRule(25);
         ir->AddTriMidPoint(0, 0.0454089951913767900476432975500142);
         ir->AddTriPoints3b(1, 0.028844733232685245264984935583748,
                            0.0183629788782333523585030359456832);
         ir->AddTriPoints3(4, 0.109481575485037054795458631340522,
                           0.0226605297177639673913028223692986);
         ir->AddTriPoints6(7, 0.141707219414879954756683250476361,
                           0.307939838764120950165155022930631,
                           0.0363789584227100543021575883096803);
         ir->AddTriPoints6(13, 0.025003534762686386073988481007746,
                           0.246672560639902693917276465411176,
                           0.0141636212655287424183685307910495);
         ir->AddTriPoints6(19, 0.0095408154002994575801528096228873,
                           0.0668032510122002657735402127620247,
                           4.71083348186641172996373548344341E-03);
         ir->SetOrder(10);
         TriangleIntRules[10] = ir;
         return ir;

      case 11: // 28 point -- 11 degree
         ir = new IntegrationRule(28);
         ir->AddTriPoints6(0, 0.0,
                           0.141129718717363295960826061941652,
                           3.68119189165027713212944752369032E-03);
         ir->AddTriMidPoint(6, 0.0439886505811161193990465846607278);
         ir->AddTriPoints3(7, 0.0259891409282873952600324854988407,
                           4.37215577686801152475821439991262E-03);
         ir->AddTriPoints3(10, 0.0942875026479224956305697762754049,
                           0.0190407859969674687575121697178070);
         ir->AddTriPoints3b(13, 0.010726449965572372516734795387128,
                            9.42772402806564602923839129555767E-03);
         ir->AddTriPoints3(16, 0.207343382614511333452934024112966,
                           0.0360798487723697630620149942932315);
         ir->AddTriPoints3b(19, 0.122184388599015809877869236727746,
                            0.0346645693527679499208828254519072);
         ir->AddTriPoints6(22, 0.0448416775891304433090523914688007,
                           0.2772206675282791551488214673424523,
                           0.0205281577146442833208261574536469);
         ir->SetOrder(11);
         TriangleIntRules[11] = ir;
         return ir;

      case 12: // 33 point - 12 degree
         ir = new IntegrationRule(33);
         ir->AddTriPoints3b(0, 2.35652204523900E-02, 1.28655332202275E-02);
         ir->AddTriPoints3b(3, 1.20551215411079E-01, 2.18462722690190E-02);
         ir->AddTriPoints3(6, 2.71210385012116E-01, 3.14291121089425E-02);
         ir->AddTriPoints3(9, 1.27576145541586E-01, 1.73980564653545E-02);
         ir->AddTriPoints3(12, 2.13173504532100E-02, 3.08313052577950E-03);
         ir->AddTriPoints6(15, 1.15343494534698E-01, 2.75713269685514E-01,
                           2.01857788831905E-02);
         ir->AddTriPoints6(21, 2.28383322222570E-02, 2.81325580989940E-01,
                           1.11783866011515E-02);
         ir->AddTriPoints6(27, 2.57340505483300E-02, 1.16251915907597E-01,
                           8.65811555432950E-03);
         ir->SetOrder(12);
         TriangleIntRules[12] = ir;
         return ir;

      case 13: // 37 point - 13 degree
         ir = new IntegrationRule(37);
         ir->AddTriPoints3b(0, 0.0,
                            2.67845189554543044455908674650066E-03);
         ir->AddTriMidPoint(3, 0.0293480398063595158995969648597808);
         ir->AddTriPoints3(4, 0.0246071886432302181878499494124643,
                           3.92538414805004016372590903990464E-03);
         ir->AddTriPoints3b(7, 0.159382493797610632566158925635800,
                            0.0253344765879434817105476355306468);
         ir->AddTriPoints3(10, 0.227900255506160619646298948153592,
                           0.0250401630452545330803738542916538);
         ir->AddTriPoints3(13, 0.116213058883517905247155321839271,
                           0.0158235572961491595176634480481793);
         ir->AddTriPoints3b(16, 0.046794039901841694097491569577008,
                            0.0157462815379843978450278590138683);
         ir->AddTriPoints6(19, 0.0227978945382486125477207592747430,
                           0.1254265183163409177176192369310890,
                           7.90126610763037567956187298486575E-03);
         ir->AddTriPoints6(25, 0.0162757709910885409437036075960413,
                           0.2909269114422506044621801030055257,
                           7.99081889046420266145965132482933E-03);
         ir->AddTriPoints6(31, 0.0897330604516053590796290561145196,
                           0.2723110556841851025078181617634414,
                           0.0182757511120486476280967518782978);
         ir->SetOrder(13);
         TriangleIntRules[13] = ir;
         return ir;

      case 14: // 42 point - 14 degree
         ir = new IntegrationRule(42);
         ir->AddTriPoints3b(0, 2.20721792756430E-02, 1.09417906847145E-02);
         ir->AddTriPoints3b(3, 1.64710561319092E-01, 1.63941767720625E-02);
         ir->AddTriPoints3(6, 2.73477528308839E-01, 2.58870522536460E-02);
         ir->AddTriPoints3(9, 1.77205532412543E-01, 2.10812943684965E-02);
         ir->AddTriPoints3(12, 6.17998830908730E-02, 7.21684983488850E-03);
         ir->AddTriPoints3(15, 1.93909612487010E-02, 2.46170180120000E-03);
         ir->AddTriPoints6(18, 5.71247574036480E-02, 1.72266687821356E-01,
                           1.23328766062820E-02);
         ir->AddTriPoints6(24, 9.29162493569720E-02, 3.36861459796345E-01,
                           1.92857553935305E-02);
         ir->AddTriPoints6(30, 1.46469500556540E-02, 2.98372882136258E-01,
                           7.21815405676700E-03);
         ir->AddTriPoints6(36, 1.26833093287200E-03, 1.18974497696957E-01,
                           2.50511441925050E-03);
         ir->SetOrder(14);
         TriangleIntRules[14] = ir;
         return ir;

      case 15: // 54 point - 15 degree
         ir = new IntegrationRule(54);
         ir->AddTriPoints3b(0, 0.0834384072617499333, 0.016330909424402645);
         ir->AddTriPoints3b(3, 0.192779070841738867, 0.01370640901568218);
         ir->AddTriPoints3(6, 0.293197167913025367, 0.01325501829935165);
         ir->AddTriPoints3(9, 0.146467786942772933, 0.014607981068243055);
         ir->AddTriPoints3(12, 0.0563628676656034333, 0.005292304033121995);
         ir->AddTriPoints3(15, 0.0165751268583703333, 0.0018073215320460175);
         ir->AddTriPoints6(18, 0.0099122033092248, 0.239534554154794445,
                           0.004263874050854718);
         ir->AddTriPoints6(24, 0.015803770630228, 0.404878807318339958,
                           0.006958088258345965);
         ir->AddTriPoints6(30, 0.00514360881697066667, 0.0950021131130448885,
                           0.0021459664703674175);
         ir->AddTriPoints6(36, 0.0489223257529888, 0.149753107322273969,
                           0.008117664640887445);
         ir->AddTriPoints6(42, 0.0687687486325192, 0.286919612441334979,
                           0.012803670460631195);
         ir->AddTriPoints6(48, 0.1684044181246992, 0.281835668099084562,
                           0.016544097765822835);
         ir->SetOrder(15);
         TriangleIntRules[15] = ir;
         return ir;

      case 16:  // 61 point - 17 degree (used for 16 as well)
      case 17:
         ir = new IntegrationRule(61);
         ir->AddTriMidPoint(0,  1.67185996454015E-02);
         ir->AddTriPoints3b(1,  5.65891888645200E-03, 2.54670772025350E-03);
         ir->AddTriPoints3b(4,  3.56473547507510E-02, 7.33543226381900E-03);
         ir->AddTriPoints3b(7,  9.95200619584370E-02, 1.21754391768360E-02);
         ir->AddTriPoints3b(10, 1.99467521245206E-01, 1.55537754344845E-02);
         ir->AddTriPoints3 (13, 2.52141267970953E-01, 1.56285556093100E-02);
         ir->AddTriPoints3 (16, 1.62047004658461E-01, 1.24078271698325E-02);
         ir->AddTriPoints3 (19, 7.58758822607460E-02, 7.02803653527850E-03);
         ir->AddTriPoints3 (22, 1.56547269678220E-02, 1.59733808688950E-03);
         ir->AddTriPoints6 (25, 1.01869288269190E-02, 3.34319867363658E-01,
                            4.05982765949650E-03);
         ir->AddTriPoints6 (31, 1.35440871671036E-01, 2.92221537796944E-01,
                            1.34028711415815E-02);
         ir->AddTriPoints6 (37, 5.44239242905830E-02, 3.19574885423190E-01,
                            9.22999660541100E-03);
         ir->AddTriPoints6 (43, 1.28685608336370E-02, 1.90704224192292E-01,
                            4.23843426716400E-03);
         ir->AddTriPoints6 (49, 6.71657824135240E-02, 1.80483211648746E-01,
                            9.14639838501250E-03);
         ir->AddTriPoints6 (55, 1.46631822248280E-02, 8.07113136795640E-02,
                            3.33281600208250E-03);
         ir->SetOrder(17);
         TriangleIntRules[16] = TriangleIntRules[17] = ir;
         return ir;

      case 18: // 73 point - 19 degree (used for 18 as well)
      case 19:
         ir = new IntegrationRule(73);
         ir->AddTriMidPoint(0,  0.0164531656944595);
         ir->AddTriPoints3b(1,  0.020780025853987, 0.005165365945636);
         ir->AddTriPoints3b(4,  0.090926214604215, 0.011193623631508);
         ir->AddTriPoints3b(7,  0.197166638701138, 0.015133062934734);
         ir->AddTriPoints3 (10, 0.255551654403098, 0.015245483901099);
         ir->AddTriPoints3 (13, 0.17707794215213,  0.0120796063708205);
         ir->AddTriPoints3 (16, 0.110061053227952, 0.0080254017934005);
         ir->AddTriPoints3 (19, 0.05552862425184,  0.004042290130892);
         ir->AddTriPoints3 (22, 0.012621863777229, 0.0010396810137425);
         ir->AddTriPoints6 (25, 0.003611417848412, 0.395754787356943,
                            0.0019424384524905);
         ir->AddTriPoints6 (31, 0.13446675453078, 0.307929983880436,
                            0.012787080306011);
         ir->AddTriPoints6 (37, 0.014446025776115, 0.26456694840652,
                            0.004440451786669);
         ir->AddTriPoints6 (43, 0.046933578838178, 0.358539352205951,
                            0.0080622733808655);
         ir->AddTriPoints6 (49, 0.002861120350567, 0.157807405968595,
                            0.0012459709087455);
         ir->AddTriPoints6 (55, 0.075050596975911, 0.223861424097916,
                            0.0091214200594755);
         ir->AddTriPoints6 (61, 0.03464707481676, 0.142421601113383,
                            0.0051292818680995);
         ir->AddTriPoints6 (67, 0.065494628082938, 0.010161119296278,
                            0.001899964427651);
         ir->SetOrder(19);
         TriangleIntRules[18] = TriangleIntRules[19] = ir;
         return ir;

      case 20: // 85 point - 20 degree
         ir = new IntegrationRule(85);
         ir->AddTriMidPoint(0, 0.01380521349884976);
         ir->AddTriPoints3b(1, 0.001500649324429,     0.00088951477366337);
         ir->AddTriPoints3b(4, 0.0941397519389508667, 0.010056199056980585);
         ir->AddTriPoints3b(7, 0.2044721240895264,    0.013408923629665785);
         ir->AddTriPoints3(10, 0.264500202532787333,  0.012261566900751005);
         ir->AddTriPoints3(13, 0.211018964092076767,  0.008197289205347695);
         ir->AddTriPoints3(16, 0.107735607171271333,  0.0073979536993248);
         ir->AddTriPoints3(19, 0.0390690878378026667, 0.0022896411388521255);
         ir->AddTriPoints3(22, 0.0111743797293296333, 0.0008259132577881085);
         ir->AddTriPoints6(25, 0.00534961818733726667, 0.0635496659083522206,
                           0.001174585454287792);
         ir->AddTriPoints6(31, 0.00795481706619893333, 0.157106918940706982,
                           0.0022329628770908965);
         ir->AddTriPoints6(37, 0.0104223982812638,     0.395642114364374018,
                           0.003049783403953986);
         ir->AddTriPoints6(43, 0.0109644147961233333,  0.273167570712910522,
                           0.0034455406635941015);
         ir->AddTriPoints6(49, 0.0385667120854623333,  0.101785382485017108,
                           0.0039987375362390815);
         ir->AddTriPoints6(55, 0.0355805078172182,     0.446658549176413815,
                           0.003693067142668012);
         ir->AddTriPoints6(61, 0.0496708163627641333,  0.199010794149503095,
                           0.00639966593932413);
         ir->AddTriPoints6(67, 0.0585197250843317333,  0.3242611836922827,
                           0.008629035587848275);
         ir->AddTriPoints6(73, 0.121497787004394267,   0.208531363210132855,
                           0.009336472951467735);
         ir->AddTriPoints6(79, 0.140710844943938733,   0.323170566536257485,
                           0.01140911202919763);
         ir->SetOrder(20);
         TriangleIntRules[20] = ir;
         return ir;

      case 21: // 126 point - 25 degree (used also for degrees from 21 to 24)
      case 22:
      case 23:
      case 24:
      case 25:
         ir = new IntegrationRule(126);
         ir->AddTriPoints3b(0, 0.0279464830731742,   0.0040027909400102085);
         ir->AddTriPoints3b(3, 0.131178601327651467, 0.00797353841619525);
         ir->AddTriPoints3b(6, 0.220221729512072267, 0.006554570615397765);
         ir->AddTriPoints3 (9, 0.298443234019804467,   0.00979150048281781);
         ir->AddTriPoints3(12, 0.2340441723373718,     0.008235442720768635);
         ir->AddTriPoints3(15, 0.151468334609017567,   0.00427363953704605);
         ir->AddTriPoints3(18, 0.112733893545993667,   0.004080942928613246);
         ir->AddTriPoints3(21, 0.0777156920915263,     0.0030605732699918895);
         ir->AddTriPoints3(24, 0.034893093614297,      0.0014542491324683325);
         ir->AddTriPoints3(27, 0.00725818462093236667, 0.00034613762283099815);
         ir->AddTriPoints6(30,  0.0012923527044422,     0.227214452153364077,
                           0.0006241445996386985);
         ir->AddTriPoints6(36,  0.0053997012721162,     0.435010554853571706,
                           0.001702376454401511);
         ir->AddTriPoints6(42,  0.006384003033975,      0.320309599272204437,
                           0.0016798271630320255);
         ir->AddTriPoints6(48,  0.00502821150199306667, 0.0917503222800051889,
                           0.000858078269748377);
         ir->AddTriPoints6(54,  0.00682675862178186667, 0.0380108358587243835,
                           0.000740428158357803);
         ir->AddTriPoints6(60,  0.0100161996399295333,  0.157425218485311668,
                           0.0017556563053643425);
         ir->AddTriPoints6(66,  0.02575781317339,       0.239889659778533193,
                           0.003696775074853242);
         ir->AddTriPoints6(72,  0.0302278981199158,     0.361943118126060531,
                           0.003991543738688279);
         ir->AddTriPoints6(78,  0.0305049901071620667,  0.0835519609548285602,
                           0.0021779813065790205);
         ir->AddTriPoints6(84,  0.0459565473625693333,  0.148443220732418205,
                           0.003682528350708916);
         ir->AddTriPoints6(90,  0.0674428005402775333,  0.283739708727534955,
                           0.005481786423209775);
         ir->AddTriPoints6(96,  0.0700450914159106,     0.406899375118787573,
                           0.00587498087177056);
         ir->AddTriPoints6(102, 0.0839115246401166,     0.194113987024892542,
                           0.005007800356899285);
         ir->AddTriPoints6(108, 0.120375535677152667,   0.32413434700070316,
                           0.00665482039381434);
         ir->AddTriPoints6(114, 0.148066899157366667,   0.229277483555980969,
                           0.00707722325261307);
         ir->AddTriPoints6(120, 0.191771865867325067,   0.325618122595983752,
                           0.007440689780584005);
         ir->SetOrder(25);
         TriangleIntRules[21] =
            TriangleIntRules[22] =
               TriangleIntRules[23] =
                  TriangleIntRules[24] =
                     TriangleIntRules[25] = ir;
         return ir;

      default:
         // Grundmann-Moller rules
         int i = (Order / 2) * 2 + 1;   // Get closest odd # >= Order
         AllocIntRule(TriangleIntRules, i);
         ir = new IntegrationRule;
         ir->GrundmannMollerSimplexRule(i/2,2);
         TriangleIntRules[i-1] = TriangleIntRules[i] = ir;
         return ir;
   }
}

// Integration rules for unit square
IntegrationRule *IntegrationRules::SquareIntegrationRule(int Order)
{
   int RealOrder = GetSegmentRealOrder(Order);
   // Order is one of {RealOrder-1,RealOrder}
   if (!HaveIntRule(SegmentIntRules, RealOrder))
   {
      SegmentIntegrationRule(RealOrder);
   }
   AllocIntRule(SquareIntRules, RealOrder); // RealOrder >= Order
   SquareIntRules[RealOrder-1] =
      SquareIntRules[RealOrder] =
         new IntegrationRule(*SegmentIntRules[RealOrder],
                             *SegmentIntRules[RealOrder]);
   return SquareIntRules[Order];
}

/** Integration rules for reference tetrahedron
 {[0,0,0],[1,0,0],[0,1,0],[0,0,1]}          */
IntegrationRule *IntegrationRules::TetrahedronIntegrationRule(int Order)
{
    IntegrationRule *ir;
    // Note: Set TetrahedronIntRules[*] to ir only *after* ir is fully
    // constructed. This is needed in multithreaded environment.
    
    // assuming that orders <= 9 are pre-allocated
    switch (Order)
    {
        case 0:  // 1 point - degree 1
        case 1:
            ir = new IntegrationRule(1);
            ir->AddTetMidPoint(0, 1./6.);
            TetrahedronIntRules[0] = TetrahedronIntRules[1] = ir;
            return ir;
            
        case 2:  // 4 points - degree 2
            ir = new IntegrationRule(4);
            // ir->AddTetPoints4(0, 0.13819660112501051518, 1./24.);
            ir->AddTetPoints4b(0, 0.58541019662496845446, 1./24.);
            TetrahedronIntRules[2] = ir;
            return ir;
            
        case 3:  // 5 points - degree 3 (negative weight)
            ir = new IntegrationRule(5);
            ir->AddTetMidPoint(0, -2./15.);
            ir->AddTetPoints4b(1, 0.5, 0.075);
            TetrahedronIntRules[3] = ir;
            return ir;
            
        case 4:  // 11 points - degree 4 (negative weight)
            ir = new IntegrationRule(11);
            ir->AddTetPoints4(0, 1./14., 343./45000.);
            ir->AddTetMidPoint(4, -74./5625.);
            ir->AddTetPoints6(5, 0.10059642383320079500, 28./1125.);
            TetrahedronIntRules[4] = ir;
            return ir;
            
        case 5:  // 14 points - degree 5
            ir = new IntegrationRule(14);
            ir->AddTetPoints6(0, 0.045503704125649649492,
                              7.0910034628469110730E-03);
            ir->AddTetPoints4(6, 0.092735250310891226402, 0.012248840519393658257);
            ir->AddTetPoints4b(10, 0.067342242210098170608,
                               0.018781320953002641800);
            TetrahedronIntRules[5] = ir;
            return ir;
            
        case 6:  // 24 points - degree 6
            ir = new IntegrationRule(24);
            ir->AddTetPoints4(0, 0.21460287125915202929,
                              6.6537917096945820166E-03);
            ir->AddTetPoints4(4, 0.040673958534611353116,
                              1.6795351758867738247E-03);
            ir->AddTetPoints4b(8, 0.032986329573173468968,
                               9.2261969239424536825E-03);
            ir->AddTetPoints12(12, 0.063661001875017525299, 0.26967233145831580803,
                               8.0357142857142857143E-03);
            TetrahedronIntRules[6] = ir;
            return ir;
            
        case 7:  // 31 points - degree 7 (negative weight)
            ir = new IntegrationRule(31);
            ir->AddTetPoints6(0, 0.0, 9.7001763668430335097E-04);
            ir->AddTetMidPoint(6, 0.018264223466108820291);
            ir->AddTetPoints4(7, 0.078213192330318064374, 0.010599941524413686916);
            ir->AddTetPoints4(11, 0.12184321666390517465,
                              -0.062517740114331851691);
            ir->AddTetPoints4b(15, 2.3825066607381275412E-03,
                               4.8914252630734993858E-03);
            ir->AddTetPoints12(19, 0.1, 0.2, 0.027557319223985890653);
            TetrahedronIntRules[7] = ir;
            return ir;
            
        case 8:  // 43 points - degree 8 (negative weight)
            ir = new IntegrationRule(43);
            ir->AddTetPoints4(0, 5.7819505051979972532E-03,
                              1.6983410909288737984E-04);
            ir->AddTetPoints4(4, 0.082103588310546723091,
                              1.9670333131339009876E-03);
            ir->AddTetPoints12(8, 0.036607749553197423679, 0.19048604193463345570,
                               2.1405191411620925965E-03);
            ir->AddTetPoints6(20, 0.050532740018894224426,
                              4.5796838244672818007E-03);
            ir->AddTetPoints12(26, 0.22906653611681113960, 0.035639582788534043717,
                               5.7044858086819185068E-03);
            ir->AddTetPoints4(38, 0.20682993161067320408, 0.014250305822866901248);
            ir->AddTetMidPoint(42, -0.020500188658639915841);
            TetrahedronIntRules[8] = ir;
            return ir;
            
        case 9: // orders 9 and higher -- Grundmann-Moller rules
            ir = new IntegrationRule;
            ir->GrundmannMollerSimplexRule(4,3);
            TetrahedronIntRules[9] = ir;
            return ir;
            
        case 10: // order 10 points 84
            ir = new IntegrationRule(84);
            std::cout << "MMCP Tet Int Rule" << std::endl;
            ir->AddTetPoint( 0 , 0.0268784744148170 , 0.0268784744148170 , 0.0268784744148170 , 0.0021449351443160 );
            ir->AddTetPoint( 1 , 0.9193645767555489 , 0.0268784744148170 , 0.0268784744148170 , 0.0021449351443160 );
            ir->AddTetPoint( 2 , 0.0268784744148170 , 0.9193645767555489 , 0.0268784744148170 , 0.0021449351443160 );
            ir->AddTetPoint( 3 , 0.0268784744148170 , 0.0268784744148170 , 0.9193645767555489 , 0.0021449351443160 );
            ir->AddTetPoint( 4 , 0.1871406758034700 , 0.1871406758034700 , 0.1871406758034700 , 0.0208266416907690 );
            ir->AddTetPoint( 5 , 0.4385779725895910 , 0.1871406758034700 , 0.1871406758034700 , 0.0208266416907690 );
            ir->AddTetPoint( 6 , 0.1871406758034700 , 0.4385779725895910 , 0.1871406758034700 , 0.0208266416907690 );
            ir->AddTetPoint( 7 , 0.1871406758034700 , 0.1871406758034700 , 0.4385779725895910 , 0.0208266416907690 );
            ir->AddTetPoint( 8 , 0.4735758351279370 , 0.0264241648720630 , 0.0264241648720630 , 0.0072101360644550 );
            ir->AddTetPoint( 9 , 0.0264241648720630 , 0.4735758351279370 , 0.0264241648720630 , 0.0072101360644550 );
            ir->AddTetPoint( 10 , 0.0264241648720630 , 0.0264241648720630 , 0.4735758351279370 , 0.0072101360644550 );
            ir->AddTetPoint( 11 , 0.4735758351279370 , 0.4735758351279370 , 0.0264241648720630 , 0.0072101360644550 );
            ir->AddTetPoint( 12 , 0.4735758351279370 , 0.0264241648720630 , 0.4735758351279370 , 0.0072101360644550 );
            ir->AddTetPoint( 13 , 0.0264241648720630 , 0.4735758351279370 , 0.4735758351279370 , 0.0072101360644550 );
            ir->AddTetPoint( 14 , 0.3520452620273560 , 0.1479547379726440 , 0.1479547379726440 , 0.0307989191597120 );
            ir->AddTetPoint( 15 , 0.1479547379726440 , 0.3520452620273560 , 0.1479547379726440 , 0.0307989191597120 );
            ir->AddTetPoint( 16 , 0.1479547379726440 , 0.1479547379726440 , 0.3520452620273560 , 0.0307989191597120 );
            ir->AddTetPoint( 17 , 0.3520452620273560 , 0.3520452620273560 , 0.1479547379726440 , 0.0307989191597120 );
            ir->AddTetPoint( 18 , 0.3520452620273560 , 0.1479547379726440 , 0.3520452620273560 , 0.0307989191597120 );
            ir->AddTetPoint( 19 , 0.1479547379726440 , 0.3520452620273560 , 0.3520452620273560 , 0.0307989191597120 );
            ir->AddTetPoint( 20 , 0.7323099096929470 , 0.0209534422200560 , 0.0209534422200560 , 0.0043578448138640 );
            ir->AddTetPoint( 21 , 0.0209534422200560 , 0.7323099096929470 , 0.0209534422200560 , 0.0043578448138640 );
            ir->AddTetPoint( 22 , 0.0209534422200560 , 0.0209534422200560 , 0.7323099096929470 , 0.0043578448138640 );
            ir->AddTetPoint( 23 , 0.2257832058669400 , 0.7323099096929470 , 0.0209534422200560 , 0.0043578448138640 );
            ir->AddTetPoint( 24 , 0.2257832058669400 , 0.0209534422200560 , 0.7323099096929470 , 0.0043578448138640 );
            ir->AddTetPoint( 25 , 0.0209534422200560 , 0.2257832058669400 , 0.7323099096929470 , 0.0043578448138640 );
            ir->AddTetPoint( 26 , 0.2257832058669400 , 0.0209534422200560 , 0.0209534422200560 , 0.0043578448138640 );
            ir->AddTetPoint( 27 , 0.0209534422200560 , 0.2257832058669400 , 0.0209534422200560 , 0.0043578448138640 );
            ir->AddTetPoint( 28 , 0.0209534422200560 , 0.0209534422200560 , 0.2257832058669400 , 0.0043578448138640 );
            ir->AddTetPoint( 29 , 0.7323099096929470 , 0.2257832058669400 , 0.0209534422200560 , 0.0043578448138640 );
            ir->AddTetPoint( 30 , 0.7323099096929470 , 0.0209534422200560 , 0.2257832058669400 , 0.0043578448138640 );
            ir->AddTetPoint( 31 , 0.0209534422200560 , 0.7323099096929470 , 0.2257832058669400 , 0.0043578448138640 );
            ir->AddTetPoint( 32 , 0.6475575940869750 , 0.0969897331234660 , 0.0969897331234660 , 0.0085935306778330 );
            ir->AddTetPoint( 33 , 0.0969897331234660 , 0.6475575940869750 , 0.0969897331234660 , 0.0085935306778330 );
            ir->AddTetPoint( 34 , 0.0969897331234660 , 0.0969897331234660 , 0.6475575940869750 , 0.0085935306778330 );
            ir->AddTetPoint( 35 , 0.1584629396660920 , 0.6475575940869750 , 0.0969897331234660 , 0.0085935306778330 );
            ir->AddTetPoint( 36 , 0.1584629396660920 , 0.0969897331234660 , 0.6475575940869750 , 0.0085935306778330 );
            ir->AddTetPoint( 37 , 0.0969897331234660 , 0.1584629396660920 , 0.6475575940869750 , 0.0085935306778330 );
            ir->AddTetPoint( 38 , 0.1584629396660920 , 0.0969897331234660 , 0.0969897331234660 , 0.0085935306778330 );
            ir->AddTetPoint( 39 , 0.0969897331234660 , 0.1584629396660920 , 0.0969897331234660 , 0.0085935306778330 );
            ir->AddTetPoint( 40 , 0.0969897331234660 , 0.0969897331234660 , 0.1584629396660920 , 0.0085935306778330 );
            ir->AddTetPoint( 41 , 0.6475575940869750 , 0.1584629396660920 , 0.0969897331234660 , 0.0085935306778330 );
            ir->AddTetPoint( 42 , 0.6475575940869750 , 0.0969897331234660 , 0.1584629396660920 , 0.0085935306778330 );
            ir->AddTetPoint( 43 , 0.0969897331234660 , 0.6475575940869750 , 0.1584629396660920 , 0.0085935306778330 );
            ir->AddTetPoint( 44 , 0.3221114318308570 , 0.3221114318308570 , 0.0336657045074290 , 0.0230006816692860 );
            ir->AddTetPoint( 45 , 0.3221114318308570 , 0.0336657045074290 , 0.3221114318308570 , 0.0230006816692860 );
            ir->AddTetPoint( 46 , 0.0336657045074290 , 0.3221114318308570 , 0.3221114318308570 , 0.0230006816692860 );
            ir->AddTetPoint( 47 , 0.3221114318308570 , 0.3221114318308570 , 0.3221114318308570 , 0.0230006816692860 );
            ir->AddTetPoint( 48 , 0.0976081628904420 , 0.0976081628904420 , 0.0118444177494980 , 0.0048630639049120 );
            ir->AddTetPoint( 49 , 0.7929392564696180 , 0.0976081628904420 , 0.0118444177494980 , 0.0048630639049120 );
            ir->AddTetPoint( 50 , 0.0976081628904420 , 0.7929392564696180 , 0.0118444177494980 , 0.0048630639049120 );
            ir->AddTetPoint( 51 , 0.0976081628904420 , 0.0118444177494980 , 0.0976081628904420 , 0.0048630639049120 );
            ir->AddTetPoint( 52 , 0.7929392564696180 , 0.0118444177494980 , 0.0976081628904420 , 0.0048630639049120 );
            ir->AddTetPoint( 53 , 0.0976081628904420 , 0.0118444177494980 , 0.7929392564696180 , 0.0048630639049120 );
            ir->AddTetPoint( 54 , 0.0118444177494980 , 0.0976081628904420 , 0.0976081628904420 , 0.0048630639049120 );
            ir->AddTetPoint( 55 , 0.0118444177494980 , 0.7929392564696180 , 0.0976081628904420 , 0.0048630639049120 );
            ir->AddTetPoint( 56 , 0.0118444177494980 , 0.0976081628904420 , 0.7929392564696180 , 0.0048630639049120 );
            ir->AddTetPoint( 57 , 0.7929392564696180 , 0.0976081628904420 , 0.0976081628904420 , 0.0048630639049120 );
            ir->AddTetPoint( 58 , 0.0976081628904420 , 0.7929392564696180 , 0.0976081628904420 , 0.0048630639049120 );
            ir->AddTetPoint( 59 , 0.0976081628904420 , 0.0976081628904420 , 0.7929392564696180 , 0.0048630639049120 );
            ir->AddTetPoint( 60 , 0.1335581607035680 , 0.2965010205431240 , 0.0287564059530710 , 0.0155951400782590 );
            ir->AddTetPoint( 61 , 0.5411844128002370 , 0.1335581607035680 , 0.0287564059530710 , 0.0155951400782590 );
            ir->AddTetPoint( 62 , 0.2965010205431240 , 0.5411844128002370 , 0.0287564059530710 , 0.0155951400782590 );
            ir->AddTetPoint( 63 , 0.2965010205431240 , 0.1335581607035680 , 0.0287564059530710 , 0.0155951400782590 );
            ir->AddTetPoint( 64 , 0.5411844128002370 , 0.2965010205431240 , 0.0287564059530710 , 0.0155951400782590 );
            ir->AddTetPoint( 65 , 0.1335581607035680 , 0.5411844128002370 , 0.0287564059530710 , 0.0155951400782590 );
            ir->AddTetPoint( 66 , 0.1335581607035680 , 0.0287564059530710 , 0.2965010205431240 , 0.0155951400782590 );
            ir->AddTetPoint( 67 , 0.5411844128002370 , 0.0287564059530710 , 0.1335581607035680 , 0.0155951400782590 );
            ir->AddTetPoint( 68 , 0.2965010205431240 , 0.0287564059530710 , 0.5411844128002370 , 0.0155951400782590 );
            ir->AddTetPoint( 69 , 0.2965010205431240 , 0.0287564059530710 , 0.1335581607035680 , 0.0155951400782590 );
            ir->AddTetPoint( 70 , 0.5411844128002370 , 0.0287564059530710 , 0.2965010205431240 , 0.0155951400782590 );
            ir->AddTetPoint( 71 , 0.1335581607035680 , 0.0287564059530710 , 0.5411844128002370 , 0.0155951400782590 );
            ir->AddTetPoint( 72 , 0.0287564059530710 , 0.1335581607035680 , 0.2965010205431240 , 0.0155951400782590 );
            ir->AddTetPoint( 73 , 0.0287564059530710 , 0.5411844128002370 , 0.1335581607035680 , 0.0155951400782590 );
            ir->AddTetPoint( 74 , 0.0287564059530710 , 0.2965010205431240 , 0.5411844128002370 , 0.0155951400782590 );
            ir->AddTetPoint( 75 , 0.0287564059530710 , 0.2965010205431240 , 0.1335581607035680 , 0.0155951400782590 );
            ir->AddTetPoint( 76 , 0.0287564059530710 , 0.5411844128002370 , 0.2965010205431240 , 0.0155951400782590 );
            ir->AddTetPoint( 77 , 0.0287564059530710 , 0.1335581607035680 , 0.5411844128002370 , 0.0155951400782590 );
            ir->AddTetPoint( 78 , 0.5411844128002370 , 0.1335581607035680 , 0.2965010205431240 , 0.0155951400782590 );
            ir->AddTetPoint( 79 , 0.2965010205431240 , 0.5411844128002370 , 0.1335581607035680 , 0.0155951400782590 );
            ir->AddTetPoint( 80 , 0.1335581607035680 , 0.2965010205431240 , 0.5411844128002370 , 0.0155951400782590 );
            ir->AddTetPoint( 81 , 0.5411844128002370 , 0.2965010205431240 , 0.1335581607035680 , 0.0155951400782590 );
            ir->AddTetPoint( 82 , 0.1335581607035680 , 0.5411844128002370 , 0.2965010205431240 , 0.0155951400782590 );
            ir->AddTetPoint( 83 , 0.2965010205431240 , 0.1335581607035680 , 0.5411844128002370 , 0.0155951400782590 );
            TetrahedronIntRules[10] = ir;
            return ir;
            
        default: // Grundmann-Moller rules
            int i = (Order / 2) * 2 + 1;   // Get closest odd # >= Order
            AllocIntRule(TetrahedronIntRules, i);
            ir = new IntegrationRule;
            ir->GrundmannMollerSimplexRule(i/2,3);
            TetrahedronIntRules[i-1] = TetrahedronIntRules[i] = ir;
            return ir;
            
    }
}
// Integration rules for reference pyramid
IntegrationRule *IntegrationRules::PyramidIntegrationRule(int Order)
{
   // This is a simple integration rule adapted from an integration
   // rule for a cube which seems to be adequate for now. We should continue
   // to search for a more appropriate integration rule designed specifically
   // for pyramid elements.
   const IntegrationRule &irc = Get(Geometry::CUBE, Order);
   int npts = irc.GetNPoints();
   AllocIntRule(PyramidIntRules, Order);
   PyramidIntRules[Order] = new IntegrationRule(npts);
   PyramidIntRules[Order]->SetOrder(Order);

   if (npts == 1)
   {
      // We handle this as a special case because with only one integration
      // point we cannot accurately integrate the quadratic factor
      // pow(1.0 - ipc.z, 2) and the resulting weight does not match the volume
      // of the reference element.
      IntegrationPoint &ipp = PyramidIntRules[Order]->IntPoint(0);
      ipp.x = 0.375;
      ipp.y = 0.375;
      ipp.z = 0.25;
      ipp.weight = 1.0 / 3.0;
   }
   else
   {
      for (int k=0; k<npts; k++)
      {
         const IntegrationPoint &ipc = irc.IntPoint(k);
         IntegrationPoint &ipp = PyramidIntRules[Order]->IntPoint(k);
         ipp.x = ipc.x * (1.0 - ipc.z);
         ipp.y = ipc.y * (1.0 - ipc.z);
         ipp.z = ipc.z;
         ipp.weight = ipc.weight * pow(1.0 - ipc.z, 2);
      }
   }
   return PyramidIntRules[Order];
}

// Integration rules for reference prism
IntegrationRule *IntegrationRules::PrismIntegrationRule(int Order)
{
   const IntegrationRule &irt = Get(Geometry::TRIANGLE, Order);
   const IntegrationRule &irs = Get(Geometry::SEGMENT, Order);
   int nt = irt.GetNPoints();
   int ns = irs.GetNPoints();
   AllocIntRule(PrismIntRules, Order);
   PrismIntRules[Order] = new IntegrationRule(nt * ns);
   PrismIntRules[Order]->SetOrder(std::min(irt.GetOrder(), irs.GetOrder()));
   while (Order < std::min(irt.GetOrder(), irs.GetOrder()))
   {
      AllocIntRule(PrismIntRules, ++Order);
      PrismIntRules[Order] = PrismIntRules[Order-1];
   }

   for (int ks=0; ks<ns; ks++)
   {
      const IntegrationPoint &ips = irs.IntPoint(ks);
      for (int kt=0; kt<nt; kt++)
      {
         int kp = ks * nt + kt;
         const IntegrationPoint &ipt = irt.IntPoint(kt);
         IntegrationPoint &ipp = PrismIntRules[Order]->IntPoint(kp);
         ipp.x = ipt.x;
         ipp.y = ipt.y;
         ipp.z = ips.x;
         ipp.weight = ipt.weight * ips.weight;
      }
   }
   return PrismIntRules[Order];
}

// Integration rules for reference cube
IntegrationRule *IntegrationRules::CubeIntegrationRule(int Order)
{
   int RealOrder = GetSegmentRealOrder(Order);
   if (!HaveIntRule(SegmentIntRules, RealOrder))
   {
      SegmentIntegrationRule(RealOrder);
   }
   AllocIntRule(CubeIntRules, RealOrder);
   CubeIntRules[RealOrder-1] =
      CubeIntRules[RealOrder] =
         new IntegrationRule(*SegmentIntRules[RealOrder],
                             *SegmentIntRules[RealOrder],
                             *SegmentIntRules[RealOrder]);
   return CubeIntRules[Order];
}

IntegrationRule *IntegrationRules::PentatopeIntegrationRule(int Order)
{
    IntegrationRule *ir;
    
#ifdef MFEM_DEBUG_INTRULES
    mfem::out << "requesting integration rules for pentatopes ( order = " << Order << " )!" << endl;
#endif
    
    switch (Order)
    {
        case 0:  // 1 point - degree 1
        case 1:
            PentatopeIntRules[0] = PentatopeIntRules[1] = ir = new IntegrationRule(1);
            ir->AddPentMidPoint(0, 1./24.);
            return ir;
            
        case 2:  // 5 points - degree 2 --
            PentatopeIntRules[2] = ir = new IntegrationRule(5);
            ir->AddPentPoint( 0 , 0.1183503419072274 , 0.1183503419072274 , 0.1183503419072274 , 0.5265986323710904 , 0.0083333333333333 );
            ir->AddPentPoint( 1 , 0.1183503419072274 , 0.1183503419072274 , 0.5265986323710904 , 0.1183503419072274 , 0.0083333333333333 );
            ir->AddPentPoint( 2 , 0.1183503419072274 , 0.5265986323710904 , 0.1183503419072274 , 0.1183503419072274 , 0.0083333333333333 );
            ir->AddPentPoint( 3 , 0.5265986323710904 , 0.1183503419072274 , 0.1183503419072274 , 0.1183503419072274 , 0.0083333333333333 );
            ir->AddPentPoint( 4 , 0.1183503419072274 , 0.1183503419072274 , 0.1183503419072274 , 0.1183503419072274 , 0.0083333333333333 );
            
            return ir;
            
        case 3:  // 15 points - degree 3 --
            PentatopeIntRules[3] = ir = new IntegrationRule(15);
            ir->AddPentPoint( 0 , 0.1100944859705240 , 0.1100944859705240 , 0.1100944859705240 , 0.5596220561179041 , 0.0036374776712251 );
            ir->AddPentPoint( 1 , 0.1100944859705240 , 0.1100944859705240 , 0.5596220561179041 , 0.1100944859705240 , 0.0036374776712251 );
            ir->AddPentPoint( 2 , 0.1100944859705240 , 0.5596220561179041 , 0.1100944859705240 , 0.1100944859705240 , 0.0036374776712251 );
            ir->AddPentPoint( 3 , 0.5596220561179041 , 0.1100944859705240 , 0.1100944859705240 , 0.1100944859705240 , 0.0036374776712251 );
            ir->AddPentPoint( 4 , 0.1100944859705240 , 0.1100944859705240 , 0.1100944859705240 , 0.1100944859705240 , 0.0036374776712251 );
            ir->AddPentPoint( 5 , 0.0781306353058967 , 0.0781306353058967 , 0.3828040470411550 , 0.3828040470411550 , 0.0023479278310541 );
            ir->AddPentPoint( 6 , 0.0781306353058967 , 0.3828040470411550 , 0.0781306353058967 , 0.3828040470411550 , 0.0023479278310541 );
            ir->AddPentPoint( 7 , 0.3828040470411550 , 0.0781306353058967 , 0.0781306353058967 , 0.3828040470411550 , 0.0023479278310541 );
            ir->AddPentPoint( 8 , 0.0781306353058967 , 0.0781306353058967 , 0.0781306353058967 , 0.3828040470411550 , 0.0023479278310541 );
            ir->AddPentPoint( 9 , 0.0781306353058967 , 0.3828040470411550 , 0.3828040470411550 , 0.0781306353058967 , 0.0023479278310541 );
            ir->AddPentPoint( 10 , 0.3828040470411550 , 0.0781306353058967 , 0.3828040470411550 , 0.0781306353058967 , 0.0023479278310541 );
            ir->AddPentPoint( 11 , 0.0781306353058967 , 0.0781306353058967 , 0.3828040470411550 , 0.0781306353058967 , 0.0023479278310541 );
            ir->AddPentPoint( 12 , 0.3828040470411550 , 0.3828040470411550 , 0.0781306353058967 , 0.0781306353058967 , 0.0023479278310541 );
            ir->AddPentPoint( 13 , 0.0781306353058967 , 0.3828040470411550 , 0.0781306353058967 , 0.0781306353058967 , 0.0023479278310541 );
            ir->AddPentPoint( 14 , 0.3828040470411550 , 0.0781306353058967 , 0.0781306353058967 , 0.0781306353058967 , 0.0023479278310541 );
            
            return ir;
            
        case 4:  // 20 points - degree 4 --
            PentatopeIntRules[4] = ir = new IntegrationRule(20);
            ir->AddPentPoint( 0 , 0.0675418538123058 , 0.0675418538123058 , 0.0675418538123058 , 0.7298325847507771 , 0.0011393080328237 );
            ir->AddPentPoint( 1 , 0.0675418538123058 , 0.0675418538123058 , 0.7298325847507771 , 0.0675418538123058 , 0.0011393080328237 );
            ir->AddPentPoint( 2 , 0.0675418538123058 , 0.7298325847507771 , 0.0675418538123058 , 0.0675418538123058 , 0.0011393080328237 );
            ir->AddPentPoint( 3 , 0.7298325847507771 , 0.0675418538123058 , 0.0675418538123058 , 0.0675418538123058 , 0.0011393080328237 );
            ir->AddPentPoint( 4 , 0.0675418538123058 , 0.0675418538123058 , 0.0675418538123058 , 0.0675418538123058 , 0.0011393080328237 );
            ir->AddPentPoint( 5 , 0.2485870141487720 , 0.2485870141487720 , 0.2485870141487720 , 0.0056519434049118 , 0.0011531957595205 );
            ir->AddPentPoint( 6 , 0.2485870141487720 , 0.2485870141487720 , 0.0056519434049118 , 0.2485870141487720 , 0.0011531957595205 );
            ir->AddPentPoint( 7 , 0.2485870141487720 , 0.0056519434049118 , 0.2485870141487720 , 0.2485870141487720 , 0.0011531957595205 );
            ir->AddPentPoint( 8 , 0.0056519434049118 , 0.2485870141487720 , 0.2485870141487720 , 0.2485870141487720 , 0.0011531957595205 );
            ir->AddPentPoint( 9 , 0.2485870141487720 , 0.2485870141487720 , 0.2485870141487720 , 0.2485870141487720 , 0.0011531957595205 );
            ir->AddPentPoint( 10 , 0.0795898176427519 , 0.0795898176427519 , 0.3806152735358721 , 0.3806152735358721 , 0.0030204147704946 );
            ir->AddPentPoint( 11 , 0.0795898176427519 , 0.3806152735358721 , 0.0795898176427519 , 0.3806152735358721 , 0.0030204147704946 );
            ir->AddPentPoint( 12 , 0.3806152735358721 , 0.0795898176427519 , 0.0795898176427519 , 0.3806152735358721 , 0.0030204147704946 );
            ir->AddPentPoint( 13 , 0.0795898176427519 , 0.0795898176427519 , 0.0795898176427519 , 0.3806152735358721 , 0.0030204147704946 );
            ir->AddPentPoint( 14 , 0.0795898176427519 , 0.3806152735358721 , 0.3806152735358721 , 0.0795898176427519 , 0.0030204147704946 );
            ir->AddPentPoint( 15 , 0.3806152735358721 , 0.0795898176427519 , 0.3806152735358721 , 0.0795898176427519 , 0.0030204147704946 );
            ir->AddPentPoint( 16 , 0.0795898176427519 , 0.0795898176427519 , 0.3806152735358721 , 0.0795898176427519 , 0.0030204147704946 );
            ir->AddPentPoint( 17 , 0.3806152735358721 , 0.3806152735358721 , 0.0795898176427519 , 0.0795898176427519 , 0.0030204147704946 );
            ir->AddPentPoint( 18 , 0.0795898176427519 , 0.3806152735358721 , 0.0795898176427519 , 0.0795898176427519 , 0.0030204147704946 );
            ir->AddPentPoint( 19 , 0.3806152735358721 , 0.0795898176427519 , 0.0795898176427519 , 0.0795898176427519 , 0.0030204147704946 );
            
            return ir;
            
        case 5:  // 30 points - degree 5 --
            PentatopeIntRules[5] = ir = new IntegrationRule(30);
            ir->AddPentPoint( 0 , 0.0710489106533826 , 0.0710489106533826 , 0.0710489106533826 , 0.7158043573864693 , 0.0001419101433670 );
            ir->AddPentPoint( 1 , 0.0710489106533826 , 0.0710489106533826 , 0.7158043573864693 , 0.0710489106533826 , 0.0001419101433670 );
            ir->AddPentPoint( 2 , 0.0710489106533826 , 0.7158043573864693 , 0.0710489106533826 , 0.0710489106533826 , 0.0001419101433670 );
            ir->AddPentPoint( 3 , 0.7158043573864693 , 0.0710489106533826 , 0.0710489106533826 , 0.0710489106533826 , 0.0001419101433670 );
            ir->AddPentPoint( 4 , 0.0710489106533826 , 0.0710489106533826 , 0.0710489106533826 , 0.0710489106533826 , 0.0001419101433670 );
            ir->AddPentPoint( 5 , 0.0883169605586097 , 0.0883169605586097 , 0.0883169605586097 , 0.6467321577655614 , 0.0020433720342941 );
            ir->AddPentPoint( 6 , 0.0883169605586097 , 0.0883169605586097 , 0.6467321577655614 , 0.0883169605586097 , 0.0020433720342941 );
            ir->AddPentPoint( 7 , 0.0883169605586097 , 0.6467321577655614 , 0.0883169605586097 , 0.0883169605586097 , 0.0020433720342941 );
            ir->AddPentPoint( 8 , 0.6467321577655614 , 0.0883169605586097 , 0.0883169605586097 , 0.0883169605586097 , 0.0020433720342941 );
            ir->AddPentPoint( 9 , 0.0883169605586097 , 0.0883169605586097 , 0.0883169605586097 , 0.0883169605586097 , 0.0020433720342941 );
            ir->AddPentPoint( 10 , 0.0220148177501016 , 0.0220148177501016 , 0.4669777733748475 , 0.4669777733748475 , 0.0004255237968879 );
            ir->AddPentPoint( 11 , 0.0220148177501016 , 0.4669777733748475 , 0.0220148177501016 , 0.4669777733748475 , 0.0004255237968879 );
            ir->AddPentPoint( 12 , 0.4669777733748475 , 0.0220148177501016 , 0.0220148177501016 , 0.4669777733748475 , 0.0004255237968879 );
            ir->AddPentPoint( 13 , 0.0220148177501016 , 0.0220148177501016 , 0.0220148177501016 , 0.4669777733748475 , 0.0004255237968879 );
            ir->AddPentPoint( 14 , 0.0220148177501016 , 0.4669777733748475 , 0.4669777733748475 , 0.0220148177501016 , 0.0004255237968879 );
            ir->AddPentPoint( 15 , 0.4669777733748475 , 0.0220148177501016 , 0.4669777733748475 , 0.0220148177501016 , 0.0004255237968879 );
            ir->AddPentPoint( 16 , 0.0220148177501016 , 0.0220148177501016 , 0.4669777733748475 , 0.0220148177501016 , 0.0004255237968879 );
            ir->AddPentPoint( 17 , 0.4669777733748475 , 0.4669777733748475 , 0.0220148177501016 , 0.0220148177501016 , 0.0004255237968879 );
            ir->AddPentPoint( 18 , 0.0220148177501016 , 0.4669777733748475 , 0.0220148177501016 , 0.0220148177501016 , 0.0004255237968879 );
            ir->AddPentPoint( 19 , 0.4669777733748475 , 0.0220148177501016 , 0.0220148177501016 , 0.0220148177501016 , 0.0004255237968879 );
            ir->AddPentPoint( 20 , 0.2941257011338569 , 0.2941257011338569 , 0.0588114482992146 , 0.0588114482992146 , 0.0026485017809482 );
            ir->AddPentPoint( 21 , 0.2941257011338569 , 0.0588114482992146 , 0.2941257011338569 , 0.0588114482992146 , 0.0026485017809482 );
            ir->AddPentPoint( 22 , 0.0588114482992146 , 0.2941257011338569 , 0.2941257011338569 , 0.0588114482992146 , 0.0026485017809482 );
            ir->AddPentPoint( 23 , 0.2941257011338569 , 0.2941257011338569 , 0.2941257011338569 , 0.0588114482992146 , 0.0026485017809482 );
            ir->AddPentPoint( 24 , 0.2941257011338569 , 0.0588114482992146 , 0.0588114482992146 , 0.2941257011338569 , 0.0026485017809482 );
            ir->AddPentPoint( 25 , 0.0588114482992146 , 0.2941257011338569 , 0.0588114482992146 , 0.2941257011338569 , 0.0026485017809482 );
            ir->AddPentPoint( 26 , 0.2941257011338569 , 0.2941257011338569 , 0.0588114482992146 , 0.2941257011338569 , 0.0026485017809482 );
            ir->AddPentPoint( 27 , 0.0588114482992146 , 0.0588114482992146 , 0.2941257011338569 , 0.2941257011338569 , 0.0026485017809482 );
            ir->AddPentPoint( 28 , 0.2941257011338569 , 0.0588114482992146 , 0.2941257011338569 , 0.2941257011338569 , 0.0026485017809482 );
            ir->AddPentPoint( 29 , 0.0588114482992146 , 0.2941257011338569 , 0.2941257011338569 , 0.2941257011338569 , 0.0026485017809482 );
            
            return ir;
            
        case 6: // 56 points - degree 6 --
            PentatopeIntRules[6] = ir = new IntegrationRule(56);
            ir->AddPentPoint( 0 , 0.2000000000000000 , 0.2000000000000000 , 0.2000000000000000 , 0.2000000000000000 , 0.0035911986465768 );
            ir->AddPentPoint( 1 , 0.0416622079805904 , 0.0416622079805904 , 0.0416622079805904 , 0.8333511680776382 , 0.0001922501317138 );
            ir->AddPentPoint( 2 , 0.0416622079805904 , 0.0416622079805904 , 0.8333511680776382 , 0.0416622079805904 , 0.0001922501317138 );
            ir->AddPentPoint( 3 , 0.0416622079805904 , 0.8333511680776382 , 0.0416622079805904 , 0.0416622079805904 , 0.0001922501317138 );
            ir->AddPentPoint( 4 , 0.8333511680776382 , 0.0416622079805904 , 0.0416622079805904 , 0.0416622079805904 , 0.0001922501317138 );
            ir->AddPentPoint( 5 , 0.0416622079805904 , 0.0416622079805904 , 0.0416622079805904 , 0.0416622079805904 , 0.0001922501317138 );
            ir->AddPentPoint( 6 , 0.2997649877085597 , 0.2997649877085597 , 0.0503525184371604 , 0.0503525184371604 , 0.0012597450269801 );
            ir->AddPentPoint( 7 , 0.2997649877085597 , 0.0503525184371604 , 0.2997649877085597 , 0.0503525184371604 , 0.0012597450269801 );
            ir->AddPentPoint( 8 , 0.0503525184371604 , 0.2997649877085597 , 0.2997649877085597 , 0.0503525184371604 , 0.0012597450269801 );
            ir->AddPentPoint( 9 , 0.2997649877085597 , 0.2997649877085597 , 0.2997649877085597 , 0.0503525184371604 , 0.0012597450269801 );
            ir->AddPentPoint( 10 , 0.2997649877085597 , 0.0503525184371604 , 0.0503525184371604 , 0.2997649877085597 , 0.0012597450269801 );
            ir->AddPentPoint( 11 , 0.0503525184371604 , 0.2997649877085597 , 0.0503525184371604 , 0.2997649877085597 , 0.0012597450269801 );
            ir->AddPentPoint( 12 , 0.2997649877085597 , 0.2997649877085597 , 0.0503525184371604 , 0.2997649877085597 , 0.0012597450269801 );
            ir->AddPentPoint( 13 , 0.0503525184371604 , 0.0503525184371604 , 0.2997649877085597 , 0.2997649877085597 , 0.0012597450269801 );
            ir->AddPentPoint( 14 , 0.2997649877085597 , 0.0503525184371604 , 0.2997649877085597 , 0.2997649877085597 , 0.0012597450269801 );
            ir->AddPentPoint( 15 , 0.0503525184371604 , 0.2997649877085597 , 0.2997649877085597 , 0.2997649877085597 , 0.0012597450269801 );
            ir->AddPentPoint( 16 , 0.0507519882160234 , 0.0507519882160234 , 0.2684660580011385 , 0.5792779773507911 , 0.0006608752665714 );
            ir->AddPentPoint( 17 , 0.0507519882160234 , 0.2684660580011385 , 0.0507519882160234 , 0.5792779773507911 , 0.0006608752665714 );
            ir->AddPentPoint( 18 , 0.2684660580011385 , 0.0507519882160234 , 0.0507519882160234 , 0.5792779773507911 , 0.0006608752665714 );
            ir->AddPentPoint( 19 , 0.0507519882160234 , 0.0507519882160234 , 0.0507519882160234 , 0.5792779773507911 , 0.0006608752665714 );
            ir->AddPentPoint( 20 , 0.0507519882160234 , 0.0507519882160234 , 0.5792779773507911 , 0.2684660580011385 , 0.0006608752665714 );
            ir->AddPentPoint( 21 , 0.0507519882160234 , 0.5792779773507911 , 0.0507519882160234 , 0.2684660580011385 , 0.0006608752665714 );
            ir->AddPentPoint( 22 , 0.5792779773507911 , 0.0507519882160234 , 0.0507519882160234 , 0.2684660580011385 , 0.0006608752665714 );
            ir->AddPentPoint( 23 , 0.0507519882160234 , 0.0507519882160234 , 0.0507519882160234 , 0.2684660580011385 , 0.0006608752665714 );
            ir->AddPentPoint( 24 , 0.0507519882160234 , 0.2684660580011385 , 0.5792779773507911 , 0.0507519882160234 , 0.0006608752665714 );
            ir->AddPentPoint( 25 , 0.2684660580011385 , 0.0507519882160234 , 0.5792779773507911 , 0.0507519882160234 , 0.0006608752665714 );
            ir->AddPentPoint( 26 , 0.0507519882160234 , 0.0507519882160234 , 0.5792779773507911 , 0.0507519882160234 , 0.0006608752665714 );
            ir->AddPentPoint( 27 , 0.0507519882160234 , 0.5792779773507911 , 0.2684660580011385 , 0.0507519882160234 , 0.0006608752665714 );
            ir->AddPentPoint( 28 , 0.5792779773507911 , 0.0507519882160234 , 0.2684660580011385 , 0.0507519882160234 , 0.0006608752665714 );
            ir->AddPentPoint( 29 , 0.0507519882160234 , 0.0507519882160234 , 0.2684660580011385 , 0.0507519882160234 , 0.0006608752665714 );
            ir->AddPentPoint( 30 , 0.2684660580011385 , 0.5792779773507911 , 0.0507519882160234 , 0.0507519882160234 , 0.0006608752665714 );
            ir->AddPentPoint( 31 , 0.0507519882160234 , 0.5792779773507911 , 0.0507519882160234 , 0.0507519882160234 , 0.0006608752665714 );
            ir->AddPentPoint( 32 , 0.5792779773507911 , 0.2684660580011385 , 0.0507519882160234 , 0.0507519882160234 , 0.0006608752665714 );
            ir->AddPentPoint( 33 , 0.0507519882160234 , 0.2684660580011385 , 0.0507519882160234 , 0.0507519882160234 , 0.0006608752665714 );
            ir->AddPentPoint( 34 , 0.5792779773507911 , 0.0507519882160234 , 0.0507519882160234 , 0.0507519882160234 , 0.0006608752665714 );
            ir->AddPentPoint( 35 , 0.2684660580011385 , 0.0507519882160234 , 0.0507519882160234 , 0.0507519882160234 , 0.0006608752665714 );
            ir->AddPentPoint( 36 , 0.1674479524682180 , 0.1674479524682180 , 0.0186976636164087 , 0.4789584789789373 , 0.0005649630880146 );
            ir->AddPentPoint( 37 , 0.1674479524682180 , 0.0186976636164087 , 0.1674479524682180 , 0.4789584789789373 , 0.0005649630880146 );
            ir->AddPentPoint( 38 , 0.0186976636164087 , 0.1674479524682180 , 0.1674479524682180 , 0.4789584789789373 , 0.0005649630880146 );
            ir->AddPentPoint( 39 , 0.1674479524682180 , 0.1674479524682180 , 0.1674479524682180 , 0.4789584789789373 , 0.0005649630880146 );
            ir->AddPentPoint( 40 , 0.1674479524682180 , 0.1674479524682180 , 0.4789584789789373 , 0.0186976636164087 , 0.0005649630880146 );
            ir->AddPentPoint( 41 , 0.1674479524682180 , 0.4789584789789373 , 0.1674479524682180 , 0.0186976636164087 , 0.0005649630880146 );
            ir->AddPentPoint( 42 , 0.4789584789789373 , 0.1674479524682180 , 0.1674479524682180 , 0.0186976636164087 , 0.0005649630880146 );
            ir->AddPentPoint( 43 , 0.1674479524682180 , 0.1674479524682180 , 0.1674479524682180 , 0.0186976636164087 , 0.0005649630880146 );
            ir->AddPentPoint( 44 , 0.1674479524682180 , 0.0186976636164087 , 0.4789584789789373 , 0.1674479524682180 , 0.0005649630880146 );
            ir->AddPentPoint( 45 , 0.0186976636164087 , 0.1674479524682180 , 0.4789584789789373 , 0.1674479524682180 , 0.0005649630880146 );
            ir->AddPentPoint( 46 , 0.1674479524682180 , 0.1674479524682180 , 0.4789584789789373 , 0.1674479524682180 , 0.0005649630880146 );
            ir->AddPentPoint( 47 , 0.1674479524682180 , 0.4789584789789373 , 0.0186976636164087 , 0.1674479524682180 , 0.0005649630880146 );
            ir->AddPentPoint( 48 , 0.4789584789789373 , 0.1674479524682180 , 0.0186976636164087 , 0.1674479524682180 , 0.0005649630880146 );
            ir->AddPentPoint( 49 , 0.1674479524682180 , 0.1674479524682180 , 0.0186976636164087 , 0.1674479524682180 , 0.0005649630880146 );
            ir->AddPentPoint( 50 , 0.0186976636164087 , 0.4789584789789373 , 0.1674479524682180 , 0.1674479524682180 , 0.0005649630880146 );
            ir->AddPentPoint( 51 , 0.1674479524682180 , 0.4789584789789373 , 0.1674479524682180 , 0.1674479524682180 , 0.0005649630880146 );
            ir->AddPentPoint( 52 , 0.4789584789789373 , 0.0186976636164087 , 0.1674479524682180 , 0.1674479524682180 , 0.0005649630880146 );
            ir->AddPentPoint( 53 , 0.1674479524682180 , 0.0186976636164087 , 0.1674479524682180 , 0.1674479524682180 , 0.0005649630880146 );
            ir->AddPentPoint( 54 , 0.4789584789789373 , 0.1674479524682180 , 0.1674479524682180 , 0.1674479524682180 , 0.0005649630880146 );
            ir->AddPentPoint( 55 , 0.0186976636164087 , 0.1674479524682180 , 0.1674479524682180 , 0.1674479524682180 , 0.0005649630880146 );
            return ir;
            
        case 7:  // 70 points - degree 7 --
            PentatopeIntRules[7] = ir = new IntegrationRule(70);
            ir->AddPentPoint( 0 , 0.0467168514448108 , 0.0467168514448108 , 0.0467168514448108 , 0.8131325942207568 , 0.0002232500865823 );
            ir->AddPentPoint( 1 , 0.0467168514448108 , 0.0467168514448108 , 0.8131325942207568 , 0.0467168514448108 , 0.0002232500865823 );
            ir->AddPentPoint( 2 , 0.0467168514448108 , 0.8131325942207568 , 0.0467168514448108 , 0.0467168514448108 , 0.0002232500865823 );
            ir->AddPentPoint( 3 , 0.8131325942207568 , 0.0467168514448108 , 0.0467168514448108 , 0.0467168514448108 , 0.0002232500865823 );
            ir->AddPentPoint( 4 , 0.0467168514448108 , 0.0467168514448108 , 0.0467168514448108 , 0.0467168514448108 , 0.0002232500865823 );
            ir->AddPentPoint( 5 , 0.2436479518537684 , 0.2436479518537684 , 0.2436479518537684 , 0.0254081925849265 , 0.0009020655120891 );
            ir->AddPentPoint( 6 , 0.2436479518537684 , 0.2436479518537684 , 0.0254081925849265 , 0.2436479518537684 , 0.0009020655120891 );
            ir->AddPentPoint( 7 , 0.2436479518537684 , 0.0254081925849265 , 0.2436479518537684 , 0.2436479518537684 , 0.0009020655120891 );
            ir->AddPentPoint( 8 , 0.0254081925849265 , 0.2436479518537684 , 0.2436479518537684 , 0.2436479518537684 , 0.0009020655120891 );
            ir->AddPentPoint( 9 , 0.2436479518537684 , 0.2436479518537684 , 0.2436479518537684 , 0.2436479518537684 , 0.0009020655120891 );
            ir->AddPentPoint( 10 , 0.1022854090355377 , 0.1022854090355377 , 0.3465718864466935 , 0.3465718864466935 , 0.0017432666832348 );
            ir->AddPentPoint( 11 , 0.1022854090355377 , 0.3465718864466935 , 0.1022854090355377 , 0.3465718864466935 , 0.0017432666832348 );
            ir->AddPentPoint( 12 , 0.3465718864466935 , 0.1022854090355377 , 0.1022854090355377 , 0.3465718864466935 , 0.0017432666832348 );
            ir->AddPentPoint( 13 , 0.1022854090355377 , 0.1022854090355377 , 0.1022854090355377 , 0.3465718864466935 , 0.0017432666832348 );
            ir->AddPentPoint( 14 , 0.1022854090355377 , 0.3465718864466935 , 0.3465718864466935 , 0.1022854090355377 , 0.0017432666832348 );
            ir->AddPentPoint( 15 , 0.3465718864466935 , 0.1022854090355377 , 0.3465718864466935 , 0.1022854090355377 , 0.0017432666832348 );
            ir->AddPentPoint( 16 , 0.1022854090355377 , 0.1022854090355377 , 0.3465718864466935 , 0.1022854090355377 , 0.0017432666832348 );
            ir->AddPentPoint( 17 , 0.3465718864466935 , 0.3465718864466935 , 0.1022854090355377 , 0.1022854090355377 , 0.0017432666832348 );
            ir->AddPentPoint( 18 , 0.1022854090355377 , 0.3465718864466935 , 0.1022854090355377 , 0.1022854090355377 , 0.0017432666832348 );
            ir->AddPentPoint( 19 , 0.3465718864466935 , 0.1022854090355377 , 0.1022854090355377 , 0.1022854090355377 , 0.0017432666832348 );
            ir->AddPentPoint( 20 , 0.3333333333333334 , 0.3333333333333334 , 0.0000000000000000 , 0.0000000000000000 , 0.0002900296721067 );
            ir->AddPentPoint( 21 , 0.3333333333333334 , 0.0000000000000000 , 0.3333333333333334 , 0.0000000000000000 , 0.0002900296721067 );
            ir->AddPentPoint( 22 , 0.0000000000000000 , 0.3333333333333334 , 0.3333333333333334 , 0.0000000000000000 , 0.0002900296721067 );
            ir->AddPentPoint( 23 , 0.3333333333333334 , 0.3333333333333334 , 0.3333333333333334 , 0.0000000000000000 , 0.0002900296721067 );
            ir->AddPentPoint( 24 , 0.3333333333333334 , 0.0000000000000000 , 0.0000000000000000 , 0.3333333333333334 , 0.0002900296721067 );
            ir->AddPentPoint( 25 , 0.0000000000000000 , 0.3333333333333334 , 0.0000000000000000 , 0.3333333333333334 , 0.0002900296721067 );
            ir->AddPentPoint( 26 , 0.3333333333333334 , 0.3333333333333334 , 0.0000000000000000 , 0.3333333333333334 , 0.0002900296721067 );
            ir->AddPentPoint( 27 , 0.0000000000000000 , 0.0000000000000000 , 0.3333333333333334 , 0.3333333333333334 , 0.0002900296721067 );
            ir->AddPentPoint( 28 , 0.3333333333333334 , 0.0000000000000000 , 0.3333333333333334 , 0.3333333333333334 , 0.0002900296721067 );
            ir->AddPentPoint( 29 , 0.0000000000000000 , 0.3333333333333334 , 0.3333333333333334 , 0.3333333333333334 , 0.0002900296721067 );
            ir->AddPentPoint( 30 , 0.0364597843329057 , 0.0364597843329057 , 0.2915694485578826 , 0.5990511984434004 , 0.0002783606971646 );
            ir->AddPentPoint( 31 , 0.0364597843329057 , 0.2915694485578826 , 0.0364597843329057 , 0.5990511984434004 , 0.0002783606971646 );
            ir->AddPentPoint( 32 , 0.2915694485578826 , 0.0364597843329057 , 0.0364597843329057 , 0.5990511984434004 , 0.0002783606971646 );
            ir->AddPentPoint( 33 , 0.0364597843329057 , 0.0364597843329057 , 0.0364597843329057 , 0.5990511984434004 , 0.0002783606971646 );
            ir->AddPentPoint( 34 , 0.0364597843329057 , 0.0364597843329057 , 0.5990511984434004 , 0.2915694485578826 , 0.0002783606971646 );
            ir->AddPentPoint( 35 , 0.0364597843329057 , 0.5990511984434004 , 0.0364597843329057 , 0.2915694485578826 , 0.0002783606971646 );
            ir->AddPentPoint( 36 , 0.5990511984434004 , 0.0364597843329057 , 0.0364597843329057 , 0.2915694485578826 , 0.0002783606971646 );
            ir->AddPentPoint( 37 , 0.0364597843329057 , 0.0364597843329057 , 0.0364597843329057 , 0.2915694485578826 , 0.0002783606971646 );
            ir->AddPentPoint( 38 , 0.0364597843329057 , 0.2915694485578826 , 0.5990511984434004 , 0.0364597843329057 , 0.0002783606971646 );
            ir->AddPentPoint( 39 , 0.2915694485578826 , 0.0364597843329057 , 0.5990511984434004 , 0.0364597843329057 , 0.0002783606971646 );
            ir->AddPentPoint( 40 , 0.0364597843329057 , 0.0364597843329057 , 0.5990511984434004 , 0.0364597843329057 , 0.0002783606971646 );
            ir->AddPentPoint( 41 , 0.0364597843329057 , 0.5990511984434004 , 0.2915694485578826 , 0.0364597843329057 , 0.0002783606971646 );
            ir->AddPentPoint( 42 , 0.5990511984434004 , 0.0364597843329057 , 0.2915694485578826 , 0.0364597843329057 , 0.0002783606971646 );
            ir->AddPentPoint( 43 , 0.0364597843329057 , 0.0364597843329057 , 0.2915694485578826 , 0.0364597843329057 , 0.0002783606971646 );
            ir->AddPentPoint( 44 , 0.2915694485578826 , 0.5990511984434004 , 0.0364597843329057 , 0.0364597843329057 , 0.0002783606971646 );
            ir->AddPentPoint( 45 , 0.0364597843329057 , 0.5990511984434004 , 0.0364597843329057 , 0.0364597843329057 , 0.0002783606971646 );
            ir->AddPentPoint( 46 , 0.5990511984434004 , 0.2915694485578826 , 0.0364597843329057 , 0.0364597843329057 , 0.0002783606971646 );
            ir->AddPentPoint( 47 , 0.0364597843329057 , 0.2915694485578826 , 0.0364597843329057 , 0.0364597843329057 , 0.0002783606971646 );
            ir->AddPentPoint( 48 , 0.5990511984434004 , 0.0364597843329057 , 0.0364597843329057 , 0.0364597843329057 , 0.0002783606971646 );
            ir->AddPentPoint( 49 , 0.2915694485578826 , 0.0364597843329057 , 0.0364597843329057 , 0.0364597843329057 , 0.0002783606971646 );
            ir->AddPentPoint( 50 , 0.1436015498551750 , 0.1436015498551750 , 0.0035133916533261 , 0.5656819587811487 , 0.0005069955588302 );
            ir->AddPentPoint( 51 , 0.1436015498551750 , 0.0035133916533261 , 0.1436015498551750 , 0.5656819587811487 , 0.0005069955588302 );
            ir->AddPentPoint( 52 , 0.0035133916533261 , 0.1436015498551750 , 0.1436015498551750 , 0.5656819587811487 , 0.0005069955588302 );
            ir->AddPentPoint( 53 , 0.1436015498551750 , 0.1436015498551750 , 0.1436015498551750 , 0.5656819587811487 , 0.0005069955588302 );
            ir->AddPentPoint( 54 , 0.1436015498551750 , 0.1436015498551750 , 0.5656819587811487 , 0.0035133916533261 , 0.0005069955588302 );
            ir->AddPentPoint( 55 , 0.1436015498551750 , 0.5656819587811487 , 0.1436015498551750 , 0.0035133916533261 , 0.0005069955588302 );
            ir->AddPentPoint( 56 , 0.5656819587811487 , 0.1436015498551750 , 0.1436015498551750 , 0.0035133916533261 , 0.0005069955588302 );
            ir->AddPentPoint( 57 , 0.1436015498551750 , 0.1436015498551750 , 0.1436015498551750 , 0.0035133916533261 , 0.0005069955588302 );
            ir->AddPentPoint( 58 , 0.1436015498551750 , 0.0035133916533261 , 0.5656819587811487 , 0.1436015498551750 , 0.0005069955588302 );
            ir->AddPentPoint( 59 , 0.0035133916533261 , 0.1436015498551750 , 0.5656819587811487 , 0.1436015498551750 , 0.0005069955588302 );
            ir->AddPentPoint( 60 , 0.1436015498551750 , 0.1436015498551750 , 0.5656819587811487 , 0.1436015498551750 , 0.0005069955588302 );
            ir->AddPentPoint( 61 , 0.1436015498551750 , 0.5656819587811487 , 0.0035133916533261 , 0.1436015498551750 , 0.0005069955588302 );
            ir->AddPentPoint( 62 , 0.5656819587811487 , 0.1436015498551750 , 0.0035133916533261 , 0.1436015498551750 , 0.0005069955588302 );
            ir->AddPentPoint( 63 , 0.1436015498551750 , 0.1436015498551750 , 0.0035133916533261 , 0.1436015498551750 , 0.0005069955588302 );
            ir->AddPentPoint( 64 , 0.0035133916533261 , 0.5656819587811487 , 0.1436015498551750 , 0.1436015498551750 , 0.0005069955588302 );
            ir->AddPentPoint( 65 , 0.1436015498551750 , 0.5656819587811487 , 0.1436015498551750 , 0.1436015498551750 , 0.0005069955588302 );
            ir->AddPentPoint( 66 , 0.5656819587811487 , 0.0035133916533261 , 0.1436015498551750 , 0.1436015498551750 , 0.0005069955588302 );
            ir->AddPentPoint( 67 , 0.1436015498551750 , 0.0035133916533261 , 0.1436015498551750 , 0.1436015498551750 , 0.0005069955588302 );
            ir->AddPentPoint( 68 , 0.5656819587811487 , 0.1436015498551750 , 0.1436015498551750 , 0.1436015498551750 , 0.0005069955588302 );
            ir->AddPentPoint( 69 , 0.0035133916533261 , 0.1436015498551750 , 0.1436015498551750 , 0.1436015498551750 , 0.0005069955588302 );
            
            return ir;
            
        case 8:  // 105 points - degree 8 --
            PentatopeIntRules[8] = ir = new IntegrationRule(105);
            ir->AddPentPoint( 0 , 0.0166503411338421 , 0.0166503411338421 , 0.0166503411338421 , 0.9333986354646318 , 0.0000190678723654 );
            ir->AddPentPoint( 1 , 0.0166503411338421 , 0.0166503411338421 , 0.9333986354646318 , 0.0166503411338421 , 0.0000190678723654 );
            ir->AddPentPoint( 2 , 0.0166503411338421 , 0.9333986354646318 , 0.0166503411338421 , 0.0166503411338421 , 0.0000190678723654 );
            ir->AddPentPoint( 3 , 0.9333986354646318 , 0.0166503411338421 , 0.0166503411338421 , 0.0166503411338421 , 0.0000190678723654 );
            ir->AddPentPoint( 4 , 0.0166503411338421 , 0.0166503411338421 , 0.0166503411338421 , 0.0166503411338421 , 0.0000190678723654 );
            ir->AddPentPoint( 5 , 0.0781784305909641 , 0.0781784305909641 , 0.0781784305909641 , 0.6872862776361435 , 0.0004710991020723 );
            ir->AddPentPoint( 6 , 0.0781784305909641 , 0.0781784305909641 , 0.6872862776361435 , 0.0781784305909641 , 0.0004710991020723 );
            ir->AddPentPoint( 7 , 0.0781784305909641 , 0.6872862776361435 , 0.0781784305909641 , 0.0781784305909641 , 0.0004710991020723 );
            ir->AddPentPoint( 8 , 0.6872862776361435 , 0.0781784305909641 , 0.0781784305909641 , 0.0781784305909641 , 0.0004710991020723 );
            ir->AddPentPoint( 9 , 0.0781784305909641 , 0.0781784305909641 , 0.0781784305909641 , 0.0781784305909641 , 0.0004710991020723 );
            ir->AddPentPoint( 10 , 0.1622138619205867 , 0.1622138619205867 , 0.1622138619205867 , 0.3511445523176532 , 0.0009532494110605 );
            ir->AddPentPoint( 11 , 0.1622138619205867 , 0.1622138619205867 , 0.3511445523176532 , 0.1622138619205867 , 0.0009532494110605 );
            ir->AddPentPoint( 12 , 0.1622138619205867 , 0.3511445523176532 , 0.1622138619205867 , 0.1622138619205867 , 0.0009532494110605 );
            ir->AddPentPoint( 13 , 0.3511445523176532 , 0.1622138619205867 , 0.1622138619205867 , 0.1622138619205867 , 0.0009532494110605 );
            ir->AddPentPoint( 14 , 0.1622138619205867 , 0.1622138619205867 , 0.1622138619205867 , 0.1622138619205867 , 0.0009532494110605 );
            ir->AddPentPoint( 15 , 0.0552897410438097 , 0.0552897410438097 , 0.4170653884342855 , 0.4170653884342855 , 0.0006636113686019 );
            ir->AddPentPoint( 16 , 0.0552897410438097 , 0.4170653884342855 , 0.0552897410438097 , 0.4170653884342855 , 0.0006636113686019 );
            ir->AddPentPoint( 17 , 0.4170653884342855 , 0.0552897410438097 , 0.0552897410438097 , 0.4170653884342855 , 0.0006636113686019 );
            ir->AddPentPoint( 18 , 0.0552897410438097 , 0.0552897410438097 , 0.0552897410438097 , 0.4170653884342855 , 0.0006636113686019 );
            ir->AddPentPoint( 19 , 0.0552897410438097 , 0.4170653884342855 , 0.4170653884342855 , 0.0552897410438097 , 0.0006636113686019 );
            ir->AddPentPoint( 20 , 0.4170653884342855 , 0.0552897410438097 , 0.4170653884342855 , 0.0552897410438097 , 0.0006636113686019 );
            ir->AddPentPoint( 21 , 0.0552897410438097 , 0.0552897410438097 , 0.4170653884342855 , 0.0552897410438097 , 0.0006636113686019 );
            ir->AddPentPoint( 22 , 0.4170653884342855 , 0.4170653884342855 , 0.0552897410438097 , 0.0552897410438097 , 0.0006636113686019 );
            ir->AddPentPoint( 23 , 0.0552897410438097 , 0.4170653884342855 , 0.0552897410438097 , 0.0552897410438097 , 0.0006636113686019 );
            ir->AddPentPoint( 24 , 0.4170653884342855 , 0.0552897410438097 , 0.0552897410438097 , 0.0552897410438097 , 0.0006636113686019 );
            ir->AddPentPoint( 25 , 0.2962849456389267 , 0.2962849456389267 , 0.0555725815416099 , 0.0555725815416099 , 0.0007344146470396 );
            ir->AddPentPoint( 26 , 0.2962849456389267 , 0.0555725815416099 , 0.2962849456389267 , 0.0555725815416099 , 0.0007344146470396 );
            ir->AddPentPoint( 27 , 0.0555725815416099 , 0.2962849456389267 , 0.2962849456389267 , 0.0555725815416099 , 0.0007344146470396 );
            ir->AddPentPoint( 28 , 0.2962849456389267 , 0.2962849456389267 , 0.2962849456389267 , 0.0555725815416099 , 0.0007344146470396 );
            ir->AddPentPoint( 29 , 0.2962849456389267 , 0.0555725815416099 , 0.0555725815416099 , 0.2962849456389267 , 0.0007344146470396 );
            ir->AddPentPoint( 30 , 0.0555725815416099 , 0.2962849456389267 , 0.0555725815416099 , 0.2962849456389267 , 0.0007344146470396 );
            ir->AddPentPoint( 31 , 0.2962849456389267 , 0.2962849456389267 , 0.0555725815416099 , 0.2962849456389267 , 0.0007344146470396 );
            ir->AddPentPoint( 32 , 0.0555725815416099 , 0.0555725815416099 , 0.2962849456389267 , 0.2962849456389267 , 0.0007344146470396 );
            ir->AddPentPoint( 33 , 0.2962849456389267 , 0.0555725815416099 , 0.2962849456389267 , 0.2962849456389267 , 0.0007344146470396 );
            ir->AddPentPoint( 34 , 0.0555725815416099 , 0.2962849456389267 , 0.2962849456389267 , 0.2962849456389267 , 0.0007344146470396 );
            ir->AddPentPoint( 35 , 0.0294607054555145 , 0.0294607054555145 , 0.1992981059127324 , 0.7123197777207242 , 0.0001252931679967 );
            ir->AddPentPoint( 36 , 0.0294607054555145 , 0.1992981059127324 , 0.0294607054555145 , 0.7123197777207242 , 0.0001252931679967 );
            ir->AddPentPoint( 37 , 0.1992981059127324 , 0.0294607054555145 , 0.0294607054555145 , 0.7123197777207242 , 0.0001252931679967 );
            ir->AddPentPoint( 38 , 0.0294607054555145 , 0.0294607054555145 , 0.0294607054555145 , 0.7123197777207242 , 0.0001252931679967 );
            ir->AddPentPoint( 39 , 0.0294607054555145 , 0.0294607054555145 , 0.7123197777207242 , 0.1992981059127324 , 0.0001252931679967 );
            ir->AddPentPoint( 40 , 0.0294607054555145 , 0.7123197777207242 , 0.0294607054555145 , 0.1992981059127324 , 0.0001252931679967 );
            ir->AddPentPoint( 41 , 0.7123197777207242 , 0.0294607054555145 , 0.0294607054555145 , 0.1992981059127324 , 0.0001252931679967 );
            ir->AddPentPoint( 42 , 0.0294607054555145 , 0.0294607054555145 , 0.0294607054555145 , 0.1992981059127324 , 0.0001252931679967 );
            ir->AddPentPoint( 43 , 0.0294607054555145 , 0.1992981059127324 , 0.7123197777207242 , 0.0294607054555145 , 0.0001252931679967 );
            ir->AddPentPoint( 44 , 0.1992981059127324 , 0.0294607054555145 , 0.7123197777207242 , 0.0294607054555145 , 0.0001252931679967 );
            ir->AddPentPoint( 45 , 0.0294607054555145 , 0.0294607054555145 , 0.7123197777207242 , 0.0294607054555145 , 0.0001252931679967 );
            ir->AddPentPoint( 46 , 0.0294607054555145 , 0.7123197777207242 , 0.1992981059127324 , 0.0294607054555145 , 0.0001252931679967 );
            ir->AddPentPoint( 47 , 0.7123197777207242 , 0.0294607054555145 , 0.1992981059127324 , 0.0294607054555145 , 0.0001252931679967 );
            ir->AddPentPoint( 48 , 0.0294607054555145 , 0.0294607054555145 , 0.1992981059127324 , 0.0294607054555145 , 0.0001252931679967 );
            ir->AddPentPoint( 49 , 0.1992981059127324 , 0.7123197777207242 , 0.0294607054555145 , 0.0294607054555145 , 0.0001252931679967 );
            ir->AddPentPoint( 50 , 0.0294607054555145 , 0.7123197777207242 , 0.0294607054555145 , 0.0294607054555145 , 0.0001252931679967 );
            ir->AddPentPoint( 51 , 0.7123197777207242 , 0.1992981059127324 , 0.0294607054555145 , 0.0294607054555145 , 0.0001252931679967 );
            ir->AddPentPoint( 52 , 0.0294607054555145 , 0.1992981059127324 , 0.0294607054555145 , 0.0294607054555145 , 0.0001252931679967 );
            ir->AddPentPoint( 53 , 0.7123197777207242 , 0.0294607054555145 , 0.0294607054555145 , 0.0294607054555145 , 0.0001252931679967 );
            ir->AddPentPoint( 54 , 0.1992981059127324 , 0.0294607054555145 , 0.0294607054555145 , 0.0294607054555145 , 0.0001252931679967 );
            ir->AddPentPoint( 55 , 0.1743045410275878 , 0.1743045410275878 , 0.0224608131707180 , 0.4546255637465186 , 0.0006724542822051 );
            ir->AddPentPoint( 56 , 0.1743045410275878 , 0.0224608131707180 , 0.1743045410275878 , 0.4546255637465186 , 0.0006724542822051 );
            ir->AddPentPoint( 57 , 0.0224608131707180 , 0.1743045410275878 , 0.1743045410275878 , 0.4546255637465186 , 0.0006724542822051 );
            ir->AddPentPoint( 58 , 0.1743045410275878 , 0.1743045410275878 , 0.1743045410275878 , 0.4546255637465186 , 0.0006724542822051 );
            ir->AddPentPoint( 59 , 0.1743045410275878 , 0.1743045410275878 , 0.4546255637465186 , 0.0224608131707180 , 0.0006724542822051 );
            ir->AddPentPoint( 60 , 0.1743045410275878 , 0.4546255637465186 , 0.1743045410275878 , 0.0224608131707180 , 0.0006724542822051 );
            ir->AddPentPoint( 61 , 0.4546255637465186 , 0.1743045410275878 , 0.1743045410275878 , 0.0224608131707180 , 0.0006724542822051 );
            ir->AddPentPoint( 62 , 0.1743045410275878 , 0.1743045410275878 , 0.1743045410275878 , 0.0224608131707180 , 0.0006724542822051 );
            ir->AddPentPoint( 63 , 0.1743045410275878 , 0.0224608131707180 , 0.4546255637465186 , 0.1743045410275878 , 0.0006724542822051 );
            ir->AddPentPoint( 64 , 0.0224608131707180 , 0.1743045410275878 , 0.4546255637465186 , 0.1743045410275878 , 0.0006724542822051 );
            ir->AddPentPoint( 65 , 0.1743045410275878 , 0.1743045410275878 , 0.4546255637465186 , 0.1743045410275878 , 0.0006724542822051 );
            ir->AddPentPoint( 66 , 0.1743045410275878 , 0.4546255637465186 , 0.0224608131707180 , 0.1743045410275878 , 0.0006724542822051 );
            ir->AddPentPoint( 67 , 0.4546255637465186 , 0.1743045410275878 , 0.0224608131707180 , 0.1743045410275878 , 0.0006724542822051 );
            ir->AddPentPoint( 68 , 0.1743045410275878 , 0.1743045410275878 , 0.0224608131707180 , 0.1743045410275878 , 0.0006724542822051 );
            ir->AddPentPoint( 69 , 0.0224608131707180 , 0.4546255637465186 , 0.1743045410275878 , 0.1743045410275878 , 0.0006724542822051 );
            ir->AddPentPoint( 70 , 0.1743045410275878 , 0.4546255637465186 , 0.1743045410275878 , 0.1743045410275878 , 0.0006724542822051 );
            ir->AddPentPoint( 71 , 0.4546255637465186 , 0.0224608131707180 , 0.1743045410275878 , 0.1743045410275878 , 0.0006724542822051 );
            ir->AddPentPoint( 72 , 0.1743045410275878 , 0.0224608131707180 , 0.1743045410275878 , 0.1743045410275878 , 0.0006724542822051 );
            ir->AddPentPoint( 73 , 0.4546255637465186 , 0.1743045410275878 , 0.1743045410275878 , 0.1743045410275878 , 0.0006724542822051 );
            ir->AddPentPoint( 74 , 0.0224608131707180 , 0.1743045410275878 , 0.1743045410275878 , 0.1743045410275878 , 0.0006724542822051 );
            ir->AddPentPoint( 75 , 0.0117056078942576 , 0.2085079973249413 , 0.2085079973249413 , 0.5595727895616022 , 0.0001504791859575 );
            ir->AddPentPoint( 76 , 0.2085079973249413 , 0.0117056078942576 , 0.2085079973249413 , 0.5595727895616022 , 0.0001504791859575 );
            ir->AddPentPoint( 77 , 0.0117056078942576 , 0.0117056078942576 , 0.2085079973249413 , 0.5595727895616022 , 0.0001504791859575 );
            ir->AddPentPoint( 78 , 0.2085079973249413 , 0.2085079973249413 , 0.0117056078942576 , 0.5595727895616022 , 0.0001504791859575 );
            ir->AddPentPoint( 79 , 0.0117056078942576 , 0.2085079973249413 , 0.0117056078942576 , 0.5595727895616022 , 0.0001504791859575 );
            ir->AddPentPoint( 80 , 0.2085079973249413 , 0.0117056078942576 , 0.0117056078942576 , 0.5595727895616022 , 0.0001504791859575 );
            ir->AddPentPoint( 81 , 0.0117056078942576 , 0.2085079973249413 , 0.5595727895616022 , 0.2085079973249413 , 0.0001504791859575 );
            ir->AddPentPoint( 82 , 0.2085079973249413 , 0.0117056078942576 , 0.5595727895616022 , 0.2085079973249413 , 0.0001504791859575 );
            ir->AddPentPoint( 83 , 0.0117056078942576 , 0.0117056078942576 , 0.5595727895616022 , 0.2085079973249413 , 0.0001504791859575 );
            ir->AddPentPoint( 84 , 0.0117056078942576 , 0.5595727895616022 , 0.2085079973249413 , 0.2085079973249413 , 0.0001504791859575 );
            ir->AddPentPoint( 85 , 0.5595727895616022 , 0.0117056078942576 , 0.2085079973249413 , 0.2085079973249413 , 0.0001504791859575 );
            ir->AddPentPoint( 86 , 0.0117056078942576 , 0.0117056078942576 , 0.2085079973249413 , 0.2085079973249413 , 0.0001504791859575 );
            ir->AddPentPoint( 87 , 0.2085079973249413 , 0.5595727895616022 , 0.0117056078942576 , 0.2085079973249413 , 0.0001504791859575 );
            ir->AddPentPoint( 88 , 0.0117056078942576 , 0.5595727895616022 , 0.0117056078942576 , 0.2085079973249413 , 0.0001504791859575 );
            ir->AddPentPoint( 89 , 0.5595727895616022 , 0.2085079973249413 , 0.0117056078942576 , 0.2085079973249413 , 0.0001504791859575 );
            ir->AddPentPoint( 90 , 0.0117056078942576 , 0.2085079973249413 , 0.0117056078942576 , 0.2085079973249413 , 0.0001504791859575 );
            ir->AddPentPoint( 91 , 0.5595727895616022 , 0.0117056078942576 , 0.0117056078942576 , 0.2085079973249413 , 0.0001504791859575 );
            ir->AddPentPoint( 92 , 0.2085079973249413 , 0.0117056078942576 , 0.0117056078942576 , 0.2085079973249413 , 0.0001504791859575 );
            ir->AddPentPoint( 93 , 0.2085079973249413 , 0.2085079973249413 , 0.5595727895616022 , 0.0117056078942576 , 0.0001504791859575 );
            ir->AddPentPoint( 94 , 0.0117056078942576 , 0.2085079973249413 , 0.5595727895616022 , 0.0117056078942576 , 0.0001504791859575 );
            ir->AddPentPoint( 95 , 0.2085079973249413 , 0.0117056078942576 , 0.5595727895616022 , 0.0117056078942576 , 0.0001504791859575 );
            ir->AddPentPoint( 96 , 0.2085079973249413 , 0.5595727895616022 , 0.2085079973249413 , 0.0117056078942576 , 0.0001504791859575 );
            ir->AddPentPoint( 97 , 0.0117056078942576 , 0.5595727895616022 , 0.2085079973249413 , 0.0117056078942576 , 0.0001504791859575 );
            ir->AddPentPoint( 98 , 0.5595727895616022 , 0.2085079973249413 , 0.2085079973249413 , 0.0117056078942576 , 0.0001504791859575 );
            ir->AddPentPoint( 99 , 0.0117056078942576 , 0.2085079973249413 , 0.2085079973249413 , 0.0117056078942576 , 0.0001504791859575 );
            ir->AddPentPoint( 100 , 0.5595727895616022 , 0.0117056078942576 , 0.2085079973249413 , 0.0117056078942576 , 0.0001504791859575 );
            ir->AddPentPoint( 101 , 0.2085079973249413 , 0.0117056078942576 , 0.2085079973249413 , 0.0117056078942576 , 0.0001504791859575 );
            ir->AddPentPoint( 102 , 0.2085079973249413 , 0.5595727895616022 , 0.0117056078942576 , 0.0117056078942576 , 0.0001504791859575 );
            ir->AddPentPoint( 103 , 0.5595727895616022 , 0.2085079973249413 , 0.0117056078942576 , 0.0117056078942576 , 0.0001504791859575 );
            ir->AddPentPoint( 104 , 0.2085079973249413 , 0.2085079973249413 , 0.0117056078942576 , 0.0117056078942576 , 0.0001504791859575 );
            
            return ir;
            
        case 9:  // 151 points - degree 9 --
            PentatopeIntRules[9] = ir = new IntegrationRule(151);
            ir->AddPentPoint( 0 , 0.2000000000000000 , 0.2000000000000000 , 0.2000000000000000 , 0.2000000000000000 , 0.0016427156665575 );
            ir->AddPentPoint( 1 , 0.0333960771949566 , 0.0333960771949566 , 0.0333960771949566 , 0.8664156912201735 , 0.0000594790909455 );
            ir->AddPentPoint( 2 , 0.0333960771949566 , 0.0333960771949566 , 0.8664156912201735 , 0.0333960771949566 , 0.0000594790909455 );
            ir->AddPentPoint( 3 , 0.0333960771949566 , 0.8664156912201735 , 0.0333960771949566 , 0.0333960771949566 , 0.0000594790909455 );
            ir->AddPentPoint( 4 , 0.8664156912201735 , 0.0333960771949566 , 0.0333960771949566 , 0.0333960771949566 , 0.0000594790909455 );
            ir->AddPentPoint( 5 , 0.0333960771949566 , 0.0333960771949566 , 0.0333960771949566 , 0.0333960771949566 , 0.0000594790909455 );
            ir->AddPentPoint( 6 , 0.2438966992434435 , 0.2438966992434435 , 0.2438966992434435 , 0.0244132030262262 , 0.0006660842122924 );
            ir->AddPentPoint( 7 , 0.2438966992434435 , 0.2438966992434435 , 0.0244132030262262 , 0.2438966992434435 , 0.0006660842122924 );
            ir->AddPentPoint( 8 , 0.2438966992434435 , 0.0244132030262262 , 0.2438966992434435 , 0.2438966992434435 , 0.0006660842122924 );
            ir->AddPentPoint( 9 , 0.0244132030262262 , 0.2438966992434435 , 0.2438966992434435 , 0.2438966992434435 , 0.0006660842122924 );
            ir->AddPentPoint( 10 , 0.2438966992434435 , 0.2438966992434435 , 0.2438966992434435 , 0.2438966992434435 , 0.0006660842122924 );
            ir->AddPentPoint( 11 , 0.0050791921061062 , 0.0050791921061062 , 0.4923812118408407 , 0.4923812118408407 , 0.0000196320910319 );
            ir->AddPentPoint( 12 , 0.0050791921061062 , 0.4923812118408407 , 0.0050791921061062 , 0.4923812118408407 , 0.0000196320910319 );
            ir->AddPentPoint( 13 , 0.4923812118408407 , 0.0050791921061062 , 0.0050791921061062 , 0.4923812118408407 , 0.0000196320910319 );
            ir->AddPentPoint( 14 , 0.0050791921061062 , 0.0050791921061062 , 0.0050791921061062 , 0.4923812118408407 , 0.0000196320910319 );
            ir->AddPentPoint( 15 , 0.0050791921061062 , 0.4923812118408407 , 0.4923812118408407 , 0.0050791921061062 , 0.0000196320910319 );
            ir->AddPentPoint( 16 , 0.4923812118408407 , 0.0050791921061062 , 0.4923812118408407 , 0.0050791921061062 , 0.0000196320910319 );
            ir->AddPentPoint( 17 , 0.0050791921061062 , 0.0050791921061062 , 0.4923812118408407 , 0.0050791921061062 , 0.0000196320910319 );
            ir->AddPentPoint( 18 , 0.4923812118408407 , 0.4923812118408407 , 0.0050791921061062 , 0.0050791921061062 , 0.0000196320910319 );
            ir->AddPentPoint( 19 , 0.0050791921061062 , 0.4923812118408407 , 0.0050791921061062 , 0.0050791921061062 , 0.0000196320910319 );
            ir->AddPentPoint( 20 , 0.4923812118408407 , 0.0050791921061062 , 0.0050791921061062 , 0.0050791921061062 , 0.0000196320910319 );
            ir->AddPentPoint( 21 , 0.2967495659603129 , 0.2967495659603129 , 0.0548756510595307 , 0.0548756510595307 , 0.0007083768588404 );
            ir->AddPentPoint( 22 , 0.2967495659603129 , 0.0548756510595307 , 0.2967495659603129 , 0.0548756510595307 , 0.0007083768588404 );
            ir->AddPentPoint( 23 , 0.0548756510595307 , 0.2967495659603129 , 0.2967495659603129 , 0.0548756510595307 , 0.0007083768588404 );
            ir->AddPentPoint( 24 , 0.2967495659603129 , 0.2967495659603129 , 0.2967495659603129 , 0.0548756510595307 , 0.0007083768588404 );
            ir->AddPentPoint( 25 , 0.2967495659603129 , 0.0548756510595307 , 0.0548756510595307 , 0.2967495659603129 , 0.0007083768588404 );
            ir->AddPentPoint( 26 , 0.0548756510595307 , 0.2967495659603129 , 0.0548756510595307 , 0.2967495659603129 , 0.0007083768588404 );
            ir->AddPentPoint( 27 , 0.2967495659603129 , 0.2967495659603129 , 0.0548756510595307 , 0.2967495659603129 , 0.0007083768588404 );
            ir->AddPentPoint( 28 , 0.0548756510595307 , 0.0548756510595307 , 0.2967495659603129 , 0.2967495659603129 , 0.0007083768588404 );
            ir->AddPentPoint( 29 , 0.2967495659603129 , 0.0548756510595307 , 0.2967495659603129 , 0.2967495659603129 , 0.0007083768588404 );
            ir->AddPentPoint( 30 , 0.0548756510595307 , 0.2967495659603129 , 0.2967495659603129 , 0.2967495659603129 , 0.0007083768588404 );
            ir->AddPentPoint( 31 , 0.0397035409474932 , 0.0397035409474932 , 0.1855665127074987 , 0.6953228644500216 , 0.0002160453626258 );
            ir->AddPentPoint( 32 , 0.0397035409474932 , 0.1855665127074987 , 0.0397035409474932 , 0.6953228644500216 , 0.0002160453626258 );
            ir->AddPentPoint( 33 , 0.1855665127074987 , 0.0397035409474932 , 0.0397035409474932 , 0.6953228644500216 , 0.0002160453626258 );
            ir->AddPentPoint( 34 , 0.0397035409474932 , 0.0397035409474932 , 0.0397035409474932 , 0.6953228644500216 , 0.0002160453626258 );
            ir->AddPentPoint( 35 , 0.0397035409474932 , 0.0397035409474932 , 0.6953228644500216 , 0.1855665127074987 , 0.0002160453626258 );
            ir->AddPentPoint( 36 , 0.0397035409474932 , 0.6953228644500216 , 0.0397035409474932 , 0.1855665127074987 , 0.0002160453626258 );
            ir->AddPentPoint( 37 , 0.6953228644500216 , 0.0397035409474932 , 0.0397035409474932 , 0.1855665127074987 , 0.0002160453626258 );
            ir->AddPentPoint( 38 , 0.0397035409474932 , 0.0397035409474932 , 0.0397035409474932 , 0.1855665127074987 , 0.0002160453626258 );
            ir->AddPentPoint( 39 , 0.0397035409474932 , 0.1855665127074987 , 0.6953228644500216 , 0.0397035409474932 , 0.0002160453626258 );
            ir->AddPentPoint( 40 , 0.1855665127074987 , 0.0397035409474932 , 0.6953228644500216 , 0.0397035409474932 , 0.0002160453626258 );
            ir->AddPentPoint( 41 , 0.0397035409474932 , 0.0397035409474932 , 0.6953228644500216 , 0.0397035409474932 , 0.0002160453626258 );
            ir->AddPentPoint( 42 , 0.0397035409474932 , 0.6953228644500216 , 0.1855665127074987 , 0.0397035409474932 , 0.0002160453626258 );
            ir->AddPentPoint( 43 , 0.6953228644500216 , 0.0397035409474932 , 0.1855665127074987 , 0.0397035409474932 , 0.0002160453626258 );
            ir->AddPentPoint( 44 , 0.0397035409474932 , 0.0397035409474932 , 0.1855665127074987 , 0.0397035409474932 , 0.0002160453626258 );
            ir->AddPentPoint( 45 , 0.1855665127074987 , 0.6953228644500216 , 0.0397035409474932 , 0.0397035409474932 , 0.0002160453626258 );
            ir->AddPentPoint( 46 , 0.0397035409474932 , 0.6953228644500216 , 0.0397035409474932 , 0.0397035409474932 , 0.0002160453626258 );
            ir->AddPentPoint( 47 , 0.6953228644500216 , 0.1855665127074987 , 0.0397035409474932 , 0.0397035409474932 , 0.0002160453626258 );
            ir->AddPentPoint( 48 , 0.0397035409474932 , 0.1855665127074987 , 0.0397035409474932 , 0.0397035409474932 , 0.0002160453626258 );
            ir->AddPentPoint( 49 , 0.6953228644500216 , 0.0397035409474932 , 0.0397035409474932 , 0.0397035409474932 , 0.0002160453626258 );
            ir->AddPentPoint( 50 , 0.1855665127074987 , 0.0397035409474932 , 0.0397035409474932 , 0.0397035409474932 , 0.0002160453626258 );
            ir->AddPentPoint( 51 , 0.0983990705650637 , 0.0983990705650637 , 0.2530927781459449 , 0.4517100101588639 , 0.0005593813535070 );
            ir->AddPentPoint( 52 , 0.0983990705650637 , 0.2530927781459449 , 0.0983990705650637 , 0.4517100101588639 , 0.0005593813535070 );
            ir->AddPentPoint( 53 , 0.2530927781459449 , 0.0983990705650637 , 0.0983990705650637 , 0.4517100101588639 , 0.0005593813535070 );
            ir->AddPentPoint( 54 , 0.0983990705650637 , 0.0983990705650637 , 0.0983990705650637 , 0.4517100101588639 , 0.0005593813535070 );
            ir->AddPentPoint( 55 , 0.0983990705650637 , 0.0983990705650637 , 0.4517100101588639 , 0.2530927781459449 , 0.0005593813535070 );
            ir->AddPentPoint( 56 , 0.0983990705650637 , 0.4517100101588639 , 0.0983990705650637 , 0.2530927781459449 , 0.0005593813535070 );
            ir->AddPentPoint( 57 , 0.4517100101588639 , 0.0983990705650637 , 0.0983990705650637 , 0.2530927781459449 , 0.0005593813535070 );
            ir->AddPentPoint( 58 , 0.0983990705650637 , 0.0983990705650637 , 0.0983990705650637 , 0.2530927781459449 , 0.0005593813535070 );
            ir->AddPentPoint( 59 , 0.0983990705650637 , 0.2530927781459449 , 0.4517100101588639 , 0.0983990705650637 , 0.0005593813535070 );
            ir->AddPentPoint( 60 , 0.2530927781459449 , 0.0983990705650637 , 0.4517100101588639 , 0.0983990705650637 , 0.0005593813535070 );
            ir->AddPentPoint( 61 , 0.0983990705650637 , 0.0983990705650637 , 0.4517100101588639 , 0.0983990705650637 , 0.0005593813535070 );
            ir->AddPentPoint( 62 , 0.0983990705650637 , 0.4517100101588639 , 0.2530927781459449 , 0.0983990705650637 , 0.0005593813535070 );
            ir->AddPentPoint( 63 , 0.4517100101588639 , 0.0983990705650637 , 0.2530927781459449 , 0.0983990705650637 , 0.0005593813535070 );
            ir->AddPentPoint( 64 , 0.0983990705650637 , 0.0983990705650637 , 0.2530927781459449 , 0.0983990705650637 , 0.0005593813535070 );
            ir->AddPentPoint( 65 , 0.2530927781459449 , 0.4517100101588639 , 0.0983990705650637 , 0.0983990705650637 , 0.0005593813535070 );
            ir->AddPentPoint( 66 , 0.0983990705650637 , 0.4517100101588639 , 0.0983990705650637 , 0.0983990705650637 , 0.0005593813535070 );
            ir->AddPentPoint( 67 , 0.4517100101588639 , 0.2530927781459449 , 0.0983990705650637 , 0.0983990705650637 , 0.0005593813535070 );
            ir->AddPentPoint( 68 , 0.0983990705650637 , 0.2530927781459449 , 0.0983990705650637 , 0.0983990705650637 , 0.0005593813535070 );
            ir->AddPentPoint( 69 , 0.4517100101588639 , 0.0983990705650637 , 0.0983990705650637 , 0.0983990705650637 , 0.0005593813535070 );
            ir->AddPentPoint( 70 , 0.2530927781459449 , 0.0983990705650637 , 0.0983990705650637 , 0.0983990705650637 , 0.0005593813535070 );
            ir->AddPentPoint( 71 , 0.1517781529659112 , 0.1517781529659112 , 0.0117376247998145 , 0.5329279163024520 , 0.0002650016066242 );
            ir->AddPentPoint( 72 , 0.1517781529659112 , 0.0117376247998145 , 0.1517781529659112 , 0.5329279163024520 , 0.0002650016066242 );
            ir->AddPentPoint( 73 , 0.0117376247998145 , 0.1517781529659112 , 0.1517781529659112 , 0.5329279163024520 , 0.0002650016066242 );
            ir->AddPentPoint( 74 , 0.1517781529659112 , 0.1517781529659112 , 0.1517781529659112 , 0.5329279163024520 , 0.0002650016066242 );
            ir->AddPentPoint( 75 , 0.1517781529659112 , 0.1517781529659112 , 0.5329279163024520 , 0.0117376247998145 , 0.0002650016066242 );
            ir->AddPentPoint( 76 , 0.1517781529659112 , 0.5329279163024520 , 0.1517781529659112 , 0.0117376247998145 , 0.0002650016066242 );
            ir->AddPentPoint( 77 , 0.5329279163024520 , 0.1517781529659112 , 0.1517781529659112 , 0.0117376247998145 , 0.0002650016066242 );
            ir->AddPentPoint( 78 , 0.1517781529659112 , 0.1517781529659112 , 0.1517781529659112 , 0.0117376247998145 , 0.0002650016066242 );
            ir->AddPentPoint( 79 , 0.1517781529659112 , 0.0117376247998145 , 0.5329279163024520 , 0.1517781529659112 , 0.0002650016066242 );
            ir->AddPentPoint( 80 , 0.0117376247998145 , 0.1517781529659112 , 0.5329279163024520 , 0.1517781529659112 , 0.0002650016066242 );
            ir->AddPentPoint( 81 , 0.1517781529659112 , 0.1517781529659112 , 0.5329279163024520 , 0.1517781529659112 , 0.0002650016066242 );
            ir->AddPentPoint( 82 , 0.1517781529659112 , 0.5329279163024520 , 0.0117376247998145 , 0.1517781529659112 , 0.0002650016066242 );
            ir->AddPentPoint( 83 , 0.5329279163024520 , 0.1517781529659112 , 0.0117376247998145 , 0.1517781529659112 , 0.0002650016066242 );
            ir->AddPentPoint( 84 , 0.1517781529659112 , 0.1517781529659112 , 0.0117376247998145 , 0.1517781529659112 , 0.0002650016066242 );
            ir->AddPentPoint( 85 , 0.0117376247998145 , 0.5329279163024520 , 0.1517781529659112 , 0.1517781529659112 , 0.0002650016066242 );
            ir->AddPentPoint( 86 , 0.1517781529659112 , 0.5329279163024520 , 0.1517781529659112 , 0.1517781529659112 , 0.0002650016066242 );
            ir->AddPentPoint( 87 , 0.5329279163024520 , 0.0117376247998145 , 0.1517781529659112 , 0.1517781529659112 , 0.0002650016066242 );
            ir->AddPentPoint( 88 , 0.1517781529659112 , 0.0117376247998145 , 0.1517781529659112 , 0.1517781529659112 , 0.0002650016066242 );
            ir->AddPentPoint( 89 , 0.5329279163024520 , 0.1517781529659112 , 0.1517781529659112 , 0.1517781529659112 , 0.0002650016066242 );
            ir->AddPentPoint( 90 , 0.0117376247998145 , 0.1517781529659112 , 0.1517781529659112 , 0.1517781529659112 , 0.0002650016066242 );
            ir->AddPentPoint( 91 , 0.0088561456625602 , 0.2163965033118229 , 0.2163965033118229 , 0.5494947020512340 , 0.0000970398164845 );
            ir->AddPentPoint( 92 , 0.2163965033118229 , 0.0088561456625602 , 0.2163965033118229 , 0.5494947020512340 , 0.0000970398164845 );
            ir->AddPentPoint( 93 , 0.0088561456625602 , 0.0088561456625602 , 0.2163965033118229 , 0.5494947020512340 , 0.0000970398164845 );
            ir->AddPentPoint( 94 , 0.2163965033118229 , 0.2163965033118229 , 0.0088561456625602 , 0.5494947020512340 , 0.0000970398164845 );
            ir->AddPentPoint( 95 , 0.0088561456625602 , 0.2163965033118229 , 0.0088561456625602 , 0.5494947020512340 , 0.0000970398164845 );
            ir->AddPentPoint( 96 , 0.2163965033118229 , 0.0088561456625602 , 0.0088561456625602 , 0.5494947020512340 , 0.0000970398164845 );
            ir->AddPentPoint( 97 , 0.0088561456625602 , 0.2163965033118229 , 0.5494947020512340 , 0.2163965033118229 , 0.0000970398164845 );
            ir->AddPentPoint( 98 , 0.2163965033118229 , 0.0088561456625602 , 0.5494947020512340 , 0.2163965033118229 , 0.0000970398164845 );
            ir->AddPentPoint( 99 , 0.0088561456625602 , 0.0088561456625602 , 0.5494947020512340 , 0.2163965033118229 , 0.0000970398164845 );
            ir->AddPentPoint( 100 , 0.0088561456625602 , 0.5494947020512340 , 0.2163965033118229 , 0.2163965033118229 , 0.0000970398164845 );
            ir->AddPentPoint( 101 , 0.5494947020512340 , 0.0088561456625602 , 0.2163965033118229 , 0.2163965033118229 , 0.0000970398164845 );
            ir->AddPentPoint( 102 , 0.0088561456625602 , 0.0088561456625602 , 0.2163965033118229 , 0.2163965033118229 , 0.0000970398164845 );
            ir->AddPentPoint( 103 , 0.2163965033118229 , 0.5494947020512340 , 0.0088561456625602 , 0.2163965033118229 , 0.0000970398164845 );
            ir->AddPentPoint( 104 , 0.0088561456625602 , 0.5494947020512340 , 0.0088561456625602 , 0.2163965033118229 , 0.0000970398164845 );
            ir->AddPentPoint( 105 , 0.5494947020512340 , 0.2163965033118229 , 0.0088561456625602 , 0.2163965033118229 , 0.0000970398164845 );
            ir->AddPentPoint( 106 , 0.0088561456625602 , 0.2163965033118229 , 0.0088561456625602 , 0.2163965033118229 , 0.0000970398164845 );
            ir->AddPentPoint( 107 , 0.5494947020512340 , 0.0088561456625602 , 0.0088561456625602 , 0.2163965033118229 , 0.0000970398164845 );
            ir->AddPentPoint( 108 , 0.2163965033118229 , 0.0088561456625602 , 0.0088561456625602 , 0.2163965033118229 , 0.0000970398164845 );
            ir->AddPentPoint( 109 , 0.2163965033118229 , 0.2163965033118229 , 0.5494947020512340 , 0.0088561456625602 , 0.0000970398164845 );
            ir->AddPentPoint( 110 , 0.0088561456625602 , 0.2163965033118229 , 0.5494947020512340 , 0.0088561456625602 , 0.0000970398164845 );
            ir->AddPentPoint( 111 , 0.2163965033118229 , 0.0088561456625602 , 0.5494947020512340 , 0.0088561456625602 , 0.0000970398164845 );
            ir->AddPentPoint( 112 , 0.2163965033118229 , 0.5494947020512340 , 0.2163965033118229 , 0.0088561456625602 , 0.0000970398164845 );
            ir->AddPentPoint( 113 , 0.0088561456625602 , 0.5494947020512340 , 0.2163965033118229 , 0.0088561456625602 , 0.0000970398164845 );
            ir->AddPentPoint( 114 , 0.5494947020512340 , 0.2163965033118229 , 0.2163965033118229 , 0.0088561456625602 , 0.0000970398164845 );
            ir->AddPentPoint( 115 , 0.0088561456625602 , 0.2163965033118229 , 0.2163965033118229 , 0.0088561456625602 , 0.0000970398164845 );
            ir->AddPentPoint( 116 , 0.5494947020512340 , 0.0088561456625602 , 0.2163965033118229 , 0.0088561456625602 , 0.0000970398164845 );
            ir->AddPentPoint( 117 , 0.2163965033118229 , 0.0088561456625602 , 0.2163965033118229 , 0.0088561456625602 , 0.0000970398164845 );
            ir->AddPentPoint( 118 , 0.2163965033118229 , 0.5494947020512340 , 0.0088561456625602 , 0.0088561456625602 , 0.0000970398164845 );
            ir->AddPentPoint( 119 , 0.5494947020512340 , 0.2163965033118229 , 0.0088561456625602 , 0.0088561456625602 , 0.0000970398164845 );
            ir->AddPentPoint( 120 , 0.2163965033118229 , 0.2163965033118229 , 0.0088561456625602 , 0.0088561456625602 , 0.0000970398164845 );
            ir->AddPentPoint( 121 , 0.0857812019328669 , 0.4135311897388892 , 0.4135311897388892 , 0.0013752166564878 , 0.0001798761345175 );
            ir->AddPentPoint( 122 , 0.4135311897388892 , 0.0857812019328669 , 0.4135311897388892 , 0.0013752166564878 , 0.0001798761345175 );
            ir->AddPentPoint( 123 , 0.0857812019328669 , 0.0857812019328669 , 0.4135311897388892 , 0.0013752166564878 , 0.0001798761345175 );
            ir->AddPentPoint( 124 , 0.4135311897388892 , 0.4135311897388892 , 0.0857812019328669 , 0.0013752166564878 , 0.0001798761345175 );
            ir->AddPentPoint( 125 , 0.0857812019328669 , 0.4135311897388892 , 0.0857812019328669 , 0.0013752166564878 , 0.0001798761345175 );
            ir->AddPentPoint( 126 , 0.4135311897388892 , 0.0857812019328669 , 0.0857812019328669 , 0.0013752166564878 , 0.0001798761345175 );
            ir->AddPentPoint( 127 , 0.0857812019328669 , 0.4135311897388892 , 0.0013752166564878 , 0.4135311897388892 , 0.0001798761345175 );
            ir->AddPentPoint( 128 , 0.4135311897388892 , 0.0857812019328669 , 0.0013752166564878 , 0.4135311897388892 , 0.0001798761345175 );
            ir->AddPentPoint( 129 , 0.0857812019328669 , 0.0857812019328669 , 0.0013752166564878 , 0.4135311897388892 , 0.0001798761345175 );
            ir->AddPentPoint( 130 , 0.0857812019328669 , 0.0013752166564878 , 0.4135311897388892 , 0.4135311897388892 , 0.0001798761345175 );
            ir->AddPentPoint( 131 , 0.0013752166564878 , 0.0857812019328669 , 0.4135311897388892 , 0.4135311897388892 , 0.0001798761345175 );
            ir->AddPentPoint( 132 , 0.0857812019328669 , 0.0857812019328669 , 0.4135311897388892 , 0.4135311897388892 , 0.0001798761345175 );
            ir->AddPentPoint( 133 , 0.4135311897388892 , 0.0013752166564878 , 0.0857812019328669 , 0.4135311897388892 , 0.0001798761345175 );
            ir->AddPentPoint( 134 , 0.0857812019328669 , 0.0013752166564878 , 0.0857812019328669 , 0.4135311897388892 , 0.0001798761345175 );
            ir->AddPentPoint( 135 , 0.0013752166564878 , 0.4135311897388892 , 0.0857812019328669 , 0.4135311897388892 , 0.0001798761345175 );
            ir->AddPentPoint( 136 , 0.0857812019328669 , 0.4135311897388892 , 0.0857812019328669 , 0.4135311897388892 , 0.0001798761345175 );
            ir->AddPentPoint( 137 , 0.0013752166564878 , 0.0857812019328669 , 0.0857812019328669 , 0.4135311897388892 , 0.0001798761345175 );
            ir->AddPentPoint( 138 , 0.4135311897388892 , 0.0857812019328669 , 0.0857812019328669 , 0.4135311897388892 , 0.0001798761345175 );
            ir->AddPentPoint( 139 , 0.4135311897388892 , 0.4135311897388892 , 0.0013752166564878 , 0.0857812019328669 , 0.0001798761345175 );
            ir->AddPentPoint( 140 , 0.0857812019328669 , 0.4135311897388892 , 0.0013752166564878 , 0.0857812019328669 , 0.0001798761345175 );
            ir->AddPentPoint( 141 , 0.4135311897388892 , 0.0857812019328669 , 0.0013752166564878 , 0.0857812019328669 , 0.0001798761345175 );
            ir->AddPentPoint( 142 , 0.4135311897388892 , 0.0013752166564878 , 0.4135311897388892 , 0.0857812019328669 , 0.0001798761345175 );
            ir->AddPentPoint( 143 , 0.0857812019328669 , 0.0013752166564878 , 0.4135311897388892 , 0.0857812019328669 , 0.0001798761345175 );
            ir->AddPentPoint( 144 , 0.0013752166564878 , 0.4135311897388892 , 0.4135311897388892 , 0.0857812019328669 , 0.0001798761345175 );
            ir->AddPentPoint( 145 , 0.0857812019328669 , 0.4135311897388892 , 0.4135311897388892 , 0.0857812019328669 , 0.0001798761345175 );
            ir->AddPentPoint( 146 , 0.0013752166564878 , 0.0857812019328669 , 0.4135311897388892 , 0.0857812019328669 , 0.0001798761345175 );
            ir->AddPentPoint( 147 , 0.4135311897388892 , 0.0857812019328669 , 0.4135311897388892 , 0.0857812019328669 , 0.0001798761345175 );
            ir->AddPentPoint( 148 , 0.4135311897388892 , 0.0013752166564878 , 0.0857812019328669 , 0.0857812019328669 , 0.0001798761345175 );
            ir->AddPentPoint( 149 , 0.0013752166564878 , 0.4135311897388892 , 0.0857812019328669 , 0.0857812019328669 , 0.0001798761345175 );
            ir->AddPentPoint( 150 , 0.4135311897388892 , 0.4135311897388892 , 0.0857812019328669 , 0.0857812019328669 , 0.0001798761345175 );
            
            return ir;
            
        case 10:  // 210 points - degree 10 --
            PentatopeIntRules[10] = ir = new IntegrationRule(210);
            ir->AddPentPoint( 0 , 0.0006800144927758 , 0.0006800144927758 , 0.0006800144927758 , 0.9972799420288967 , 0.0000020389157402 );
            ir->AddPentPoint( 1 , 0.0006800144927758 , 0.0006800144927758 , 0.9972799420288967 , 0.0006800144927758 , 0.0000020389157402 );
            ir->AddPentPoint( 2 , 0.0006800144927758 , 0.9972799420288967 , 0.0006800144927758 , 0.0006800144927758 , 0.0000020389157402 );
            ir->AddPentPoint( 3 , 0.9972799420288967 , 0.0006800144927758 , 0.0006800144927758 , 0.0006800144927758 , 0.0000020389157402 );
            ir->AddPentPoint( 4 , 0.0006800144927758 , 0.0006800144927758 , 0.0006800144927758 , 0.0006800144927758 , 0.0000020389157402 );
            ir->AddPentPoint( 5 , 0.0989918549880847 , 0.0989918549880847 , 0.0989918549880847 , 0.6040325800476612 , 0.0002561832835508 );
            ir->AddPentPoint( 6 , 0.0989918549880847 , 0.0989918549880847 , 0.6040325800476612 , 0.0989918549880847 , 0.0002561832835508 );
            ir->AddPentPoint( 7 , 0.0989918549880847 , 0.6040325800476612 , 0.0989918549880847 , 0.0989918549880847 , 0.0002561832835508 );
            ir->AddPentPoint( 8 , 0.6040325800476612 , 0.0989918549880847 , 0.0989918549880847 , 0.0989918549880847 , 0.0002561832835508 );
            ir->AddPentPoint( 9 , 0.0989918549880847 , 0.0989918549880847 , 0.0989918549880847 , 0.0989918549880847 , 0.0002561832835508 );
            ir->AddPentPoint( 10 , 0.1231204513985415 , 0.1231204513985415 , 0.1231204513985415 , 0.5075181944058341 , 0.0002541581830825 );
            ir->AddPentPoint( 11 , 0.1231204513985415 , 0.1231204513985415 , 0.5075181944058341 , 0.1231204513985415 , 0.0002541581830825 );
            ir->AddPentPoint( 12 , 0.1231204513985415 , 0.5075181944058341 , 0.1231204513985415 , 0.1231204513985415 , 0.0002541581830825 );
            ir->AddPentPoint( 13 , 0.5075181944058341 , 0.1231204513985415 , 0.1231204513985415 , 0.1231204513985415 , 0.0002541581830825 );
            ir->AddPentPoint( 14 , 0.1231204513985415 , 0.1231204513985415 , 0.1231204513985415 , 0.1231204513985415 , 0.0002541581830825 );
            ir->AddPentPoint( 15 , 0.2425841483279530 , 0.2425841483279530 , 0.2425841483279530 , 0.0296634066881880 , 0.0004061429331263 );
            ir->AddPentPoint( 16 , 0.2425841483279530 , 0.2425841483279530 , 0.0296634066881880 , 0.2425841483279530 , 0.0004061429331263 );
            ir->AddPentPoint( 17 , 0.2425841483279530 , 0.0296634066881880 , 0.2425841483279530 , 0.2425841483279530 , 0.0004061429331263 );
            ir->AddPentPoint( 18 , 0.0296634066881880 , 0.2425841483279530 , 0.2425841483279530 , 0.2425841483279530 , 0.0004061429331263 );
            ir->AddPentPoint( 19 , 0.2425841483279530 , 0.2425841483279530 , 0.2425841483279530 , 0.2425841483279530 , 0.0004061429331263 );
            ir->AddPentPoint( 20 , 0.1277741480154265 , 0.1277741480154265 , 0.3083387779768602 , 0.3083387779768602 , 0.0007365854844835 );
            ir->AddPentPoint( 21 , 0.1277741480154265 , 0.3083387779768602 , 0.1277741480154265 , 0.3083387779768602 , 0.0007365854844835 );
            ir->AddPentPoint( 22 , 0.3083387779768602 , 0.1277741480154265 , 0.1277741480154265 , 0.3083387779768602 , 0.0007365854844835 );
            ir->AddPentPoint( 23 , 0.1277741480154265 , 0.1277741480154265 , 0.1277741480154265 , 0.3083387779768602 , 0.0007365854844835 );
            ir->AddPentPoint( 24 , 0.1277741480154265 , 0.3083387779768602 , 0.3083387779768602 , 0.1277741480154265 , 0.0007365854844835 );
            ir->AddPentPoint( 25 , 0.3083387779768602 , 0.1277741480154265 , 0.3083387779768602 , 0.1277741480154265 , 0.0007365854844835 );
            ir->AddPentPoint( 26 , 0.1277741480154265 , 0.1277741480154265 , 0.3083387779768602 , 0.1277741480154265 , 0.0007365854844835 );
            ir->AddPentPoint( 27 , 0.3083387779768602 , 0.3083387779768602 , 0.1277741480154265 , 0.1277741480154265 , 0.0007365854844835 );
            ir->AddPentPoint( 28 , 0.1277741480154265 , 0.3083387779768602 , 0.1277741480154265 , 0.1277741480154265 , 0.0007365854844835 );
            ir->AddPentPoint( 29 , 0.3083387779768602 , 0.1277741480154265 , 0.1277741480154265 , 0.1277741480154265 , 0.0007365854844835 );
            ir->AddPentPoint( 30 , 0.3077308050572485 , 0.3077308050572485 , 0.0384037924141272 , 0.0384037924141272 , 0.0004223673935348 );
            ir->AddPentPoint( 31 , 0.3077308050572485 , 0.0384037924141272 , 0.3077308050572485 , 0.0384037924141272 , 0.0004223673935348 );
            ir->AddPentPoint( 32 , 0.0384037924141272 , 0.3077308050572485 , 0.3077308050572485 , 0.0384037924141272 , 0.0004223673935348 );
            ir->AddPentPoint( 33 , 0.3077308050572485 , 0.3077308050572485 , 0.3077308050572485 , 0.0384037924141272 , 0.0004223673935348 );
            ir->AddPentPoint( 34 , 0.3077308050572485 , 0.0384037924141272 , 0.0384037924141272 , 0.3077308050572485 , 0.0004223673935348 );
            ir->AddPentPoint( 35 , 0.0384037924141272 , 0.3077308050572485 , 0.0384037924141272 , 0.3077308050572485 , 0.0004223673935348 );
            ir->AddPentPoint( 36 , 0.3077308050572485 , 0.3077308050572485 , 0.0384037924141272 , 0.3077308050572485 , 0.0004223673935348 );
            ir->AddPentPoint( 37 , 0.0384037924141272 , 0.0384037924141272 , 0.3077308050572485 , 0.3077308050572485 , 0.0004223673935348 );
            ir->AddPentPoint( 38 , 0.3077308050572485 , 0.0384037924141272 , 0.3077308050572485 , 0.3077308050572485 , 0.0004223673935348 );
            ir->AddPentPoint( 39 , 0.0384037924141272 , 0.3077308050572485 , 0.3077308050572485 , 0.3077308050572485 , 0.0004223673935348 );
            ir->AddPentPoint( 40 , 0.0059989160658866 , 0.0059989160658866 , 0.2053177373179084 , 0.7766855144844318 , 0.0000102898818133 );
            ir->AddPentPoint( 41 , 0.0059989160658866 , 0.2053177373179084 , 0.0059989160658866 , 0.7766855144844318 , 0.0000102898818133 );
            ir->AddPentPoint( 42 , 0.2053177373179084 , 0.0059989160658866 , 0.0059989160658866 , 0.7766855144844318 , 0.0000102898818133 );
            ir->AddPentPoint( 43 , 0.0059989160658866 , 0.0059989160658866 , 0.0059989160658866 , 0.7766855144844318 , 0.0000102898818133 );
            ir->AddPentPoint( 44 , 0.0059989160658866 , 0.0059989160658866 , 0.7766855144844318 , 0.2053177373179084 , 0.0000102898818133 );
            ir->AddPentPoint( 45 , 0.0059989160658866 , 0.7766855144844318 , 0.0059989160658866 , 0.2053177373179084 , 0.0000102898818133 );
            ir->AddPentPoint( 46 , 0.7766855144844318 , 0.0059989160658866 , 0.0059989160658866 , 0.2053177373179084 , 0.0000102898818133 );
            ir->AddPentPoint( 47 , 0.0059989160658866 , 0.0059989160658866 , 0.0059989160658866 , 0.2053177373179084 , 0.0000102898818133 );
            ir->AddPentPoint( 48 , 0.0059989160658866 , 0.2053177373179084 , 0.7766855144844318 , 0.0059989160658866 , 0.0000102898818133 );
            ir->AddPentPoint( 49 , 0.2053177373179084 , 0.0059989160658866 , 0.7766855144844318 , 0.0059989160658866 , 0.0000102898818133 );
            ir->AddPentPoint( 50 , 0.0059989160658866 , 0.0059989160658866 , 0.7766855144844318 , 0.0059989160658866 , 0.0000102898818133 );
            ir->AddPentPoint( 51 , 0.0059989160658866 , 0.7766855144844318 , 0.2053177373179084 , 0.0059989160658866 , 0.0000102898818133 );
            ir->AddPentPoint( 52 , 0.7766855144844318 , 0.0059989160658866 , 0.2053177373179084 , 0.0059989160658866 , 0.0000102898818133 );
            ir->AddPentPoint( 53 , 0.0059989160658866 , 0.0059989160658866 , 0.2053177373179084 , 0.0059989160658866 , 0.0000102898818133 );
            ir->AddPentPoint( 54 , 0.2053177373179084 , 0.7766855144844318 , 0.0059989160658866 , 0.0059989160658866 , 0.0000102898818133 );
            ir->AddPentPoint( 55 , 0.0059989160658866 , 0.7766855144844318 , 0.0059989160658866 , 0.0059989160658866 , 0.0000102898818133 );
            ir->AddPentPoint( 56 , 0.7766855144844318 , 0.2053177373179084 , 0.0059989160658866 , 0.0059989160658866 , 0.0000102898818133 );
            ir->AddPentPoint( 57 , 0.0059989160658866 , 0.2053177373179084 , 0.0059989160658866 , 0.0059989160658866 , 0.0000102898818133 );
            ir->AddPentPoint( 58 , 0.7766855144844318 , 0.0059989160658866 , 0.0059989160658866 , 0.0059989160658866 , 0.0000102898818133 );
            ir->AddPentPoint( 59 , 0.2053177373179084 , 0.0059989160658866 , 0.0059989160658866 , 0.0059989160658866 , 0.0000102898818133 );
            ir->AddPentPoint( 60 , 0.0454971655794992 , 0.0454971655794992 , 0.2553045488676483 , 0.6082039543938539 , 0.0001971825981252 );
            ir->AddPentPoint( 61 , 0.0454971655794992 , 0.2553045488676483 , 0.0454971655794992 , 0.6082039543938539 , 0.0001971825981252 );
            ir->AddPentPoint( 62 , 0.2553045488676483 , 0.0454971655794992 , 0.0454971655794992 , 0.6082039543938539 , 0.0001971825981252 );
            ir->AddPentPoint( 63 , 0.0454971655794992 , 0.0454971655794992 , 0.0454971655794992 , 0.6082039543938539 , 0.0001971825981252 );
            ir->AddPentPoint( 64 , 0.0454971655794992 , 0.0454971655794992 , 0.6082039543938539 , 0.2553045488676483 , 0.0001971825981252 );
            ir->AddPentPoint( 65 , 0.0454971655794992 , 0.6082039543938539 , 0.0454971655794992 , 0.2553045488676483 , 0.0001971825981252 );
            ir->AddPentPoint( 66 , 0.6082039543938539 , 0.0454971655794992 , 0.0454971655794992 , 0.2553045488676483 , 0.0001971825981252 );
            ir->AddPentPoint( 67 , 0.0454971655794992 , 0.0454971655794992 , 0.0454971655794992 , 0.2553045488676483 , 0.0001971825981252 );
            ir->AddPentPoint( 68 , 0.0454971655794992 , 0.2553045488676483 , 0.6082039543938539 , 0.0454971655794992 , 0.0001971825981252 );
            ir->AddPentPoint( 69 , 0.2553045488676483 , 0.0454971655794992 , 0.6082039543938539 , 0.0454971655794992 , 0.0001971825981252 );
            ir->AddPentPoint( 70 , 0.0454971655794992 , 0.0454971655794992 , 0.6082039543938539 , 0.0454971655794992 , 0.0001971825981252 );
            ir->AddPentPoint( 71 , 0.0454971655794992 , 0.6082039543938539 , 0.2553045488676483 , 0.0454971655794992 , 0.0001971825981252 );
            ir->AddPentPoint( 72 , 0.6082039543938539 , 0.0454971655794992 , 0.2553045488676483 , 0.0454971655794992 , 0.0001971825981252 );
            ir->AddPentPoint( 73 , 0.0454971655794992 , 0.0454971655794992 , 0.2553045488676483 , 0.0454971655794992 , 0.0001971825981252 );
            ir->AddPentPoint( 74 , 0.2553045488676483 , 0.6082039543938539 , 0.0454971655794992 , 0.0454971655794992 , 0.0001971825981252 );
            ir->AddPentPoint( 75 , 0.0454971655794992 , 0.6082039543938539 , 0.0454971655794992 , 0.0454971655794992 , 0.0001971825981252 );
            ir->AddPentPoint( 76 , 0.6082039543938539 , 0.2553045488676483 , 0.0454971655794992 , 0.0454971655794992 , 0.0001971825981252 );
            ir->AddPentPoint( 77 , 0.0454971655794992 , 0.2553045488676483 , 0.0454971655794992 , 0.0454971655794992 , 0.0001971825981252 );
            ir->AddPentPoint( 78 , 0.6082039543938539 , 0.0454971655794992 , 0.0454971655794992 , 0.0454971655794992 , 0.0001971825981252 );
            ir->AddPentPoint( 79 , 0.2553045488676483 , 0.0454971655794992 , 0.0454971655794992 , 0.0454971655794992 , 0.0001971825981252 );
            ir->AddPentPoint( 80 , 0.0729383962178788 , 0.0729383962178788 , 0.0000000000000001 , 0.7811848113463635 , 0.0000723045529434 );
            ir->AddPentPoint( 81 , 0.0729383962178788 , 0.0000000000000001 , 0.0729383962178788 , 0.7811848113463635 , 0.0000723045529434 );
            ir->AddPentPoint( 82 , 0.0000000000000001 , 0.0729383962178788 , 0.0729383962178788 , 0.7811848113463635 , 0.0000723045529434 );
            ir->AddPentPoint( 83 , 0.0729383962178788 , 0.0729383962178788 , 0.0729383962178788 , 0.7811848113463635 , 0.0000723045529434 );
            ir->AddPentPoint( 84 , 0.0729383962178788 , 0.0729383962178788 , 0.7811848113463635 , 0.0000000000000001 , 0.0000723045529434 );
            ir->AddPentPoint( 85 , 0.0729383962178788 , 0.7811848113463635 , 0.0729383962178788 , 0.0000000000000001 , 0.0000723045529434 );
            ir->AddPentPoint( 86 , 0.7811848113463635 , 0.0729383962178788 , 0.0729383962178788 , 0.0000000000000001 , 0.0000723045529434 );
            ir->AddPentPoint( 87 , 0.0729383962178788 , 0.0729383962178788 , 0.0729383962178788 , 0.0000000000000001 , 0.0000723045529434 );
            ir->AddPentPoint( 88 , 0.0729383962178788 , 0.0000000000000001 , 0.7811848113463635 , 0.0729383962178788 , 0.0000723045529434 );
            ir->AddPentPoint( 89 , 0.0000000000000001 , 0.0729383962178788 , 0.7811848113463635 , 0.0729383962178788 , 0.0000723045529434 );
            ir->AddPentPoint( 90 , 0.0729383962178788 , 0.0729383962178788 , 0.7811848113463635 , 0.0729383962178788 , 0.0000723045529434 );
            ir->AddPentPoint( 91 , 0.0729383962178788 , 0.7811848113463635 , 0.0000000000000001 , 0.0729383962178788 , 0.0000723045529434 );
            ir->AddPentPoint( 92 , 0.7811848113463635 , 0.0729383962178788 , 0.0000000000000001 , 0.0729383962178788 , 0.0000723045529434 );
            ir->AddPentPoint( 93 , 0.0729383962178788 , 0.0729383962178788 , 0.0000000000000001 , 0.0729383962178788 , 0.0000723045529434 );
            ir->AddPentPoint( 94 , 0.0000000000000001 , 0.7811848113463635 , 0.0729383962178788 , 0.0729383962178788 , 0.0000723045529434 );
            ir->AddPentPoint( 95 , 0.0729383962178788 , 0.7811848113463635 , 0.0729383962178788 , 0.0729383962178788 , 0.0000723045529434 );
            ir->AddPentPoint( 96 , 0.7811848113463635 , 0.0000000000000001 , 0.0729383962178788 , 0.0729383962178788 , 0.0000723045529434 );
            ir->AddPentPoint( 97 , 0.0729383962178788 , 0.0000000000000001 , 0.0729383962178788 , 0.0729383962178788 , 0.0000723045529434 );
            ir->AddPentPoint( 98 , 0.7811848113463635 , 0.0729383962178788 , 0.0729383962178788 , 0.0729383962178788 , 0.0000723045529434 );
            ir->AddPentPoint( 99 , 0.0000000000000001 , 0.0729383962178788 , 0.0729383962178788 , 0.0729383962178788 , 0.0000723045529434 );
            ir->AddPentPoint( 100 , 0.1738181457692562 , 0.1738181457692562 , 0.0227203202344919 , 0.4558252424577396 , 0.0003272598432739 );
            ir->AddPentPoint( 101 , 0.1738181457692562 , 0.0227203202344919 , 0.1738181457692562 , 0.4558252424577396 , 0.0003272598432739 );
            ir->AddPentPoint( 102 , 0.0227203202344919 , 0.1738181457692562 , 0.1738181457692562 , 0.4558252424577396 , 0.0003272598432739 );
            ir->AddPentPoint( 103 , 0.1738181457692562 , 0.1738181457692562 , 0.1738181457692562 , 0.4558252424577396 , 0.0003272598432739 );
            ir->AddPentPoint( 104 , 0.1738181457692562 , 0.1738181457692562 , 0.4558252424577396 , 0.0227203202344919 , 0.0003272598432739 );
            ir->AddPentPoint( 105 , 0.1738181457692562 , 0.4558252424577396 , 0.1738181457692562 , 0.0227203202344919 , 0.0003272598432739 );
            ir->AddPentPoint( 106 , 0.4558252424577396 , 0.1738181457692562 , 0.1738181457692562 , 0.0227203202344919 , 0.0003272598432739 );
            ir->AddPentPoint( 107 , 0.1738181457692562 , 0.1738181457692562 , 0.1738181457692562 , 0.0227203202344919 , 0.0003272598432739 );
            ir->AddPentPoint( 108 , 0.1738181457692562 , 0.0227203202344919 , 0.4558252424577396 , 0.1738181457692562 , 0.0003272598432739 );
            ir->AddPentPoint( 109 , 0.0227203202344919 , 0.1738181457692562 , 0.4558252424577396 , 0.1738181457692562 , 0.0003272598432739 );
            ir->AddPentPoint( 110 , 0.1738181457692562 , 0.1738181457692562 , 0.4558252424577396 , 0.1738181457692562 , 0.0003272598432739 );
            ir->AddPentPoint( 111 , 0.1738181457692562 , 0.4558252424577396 , 0.0227203202344919 , 0.1738181457692562 , 0.0003272598432739 );
            ir->AddPentPoint( 112 , 0.4558252424577396 , 0.1738181457692562 , 0.0227203202344919 , 0.1738181457692562 , 0.0003272598432739 );
            ir->AddPentPoint( 113 , 0.1738181457692562 , 0.1738181457692562 , 0.0227203202344919 , 0.1738181457692562 , 0.0003272598432739 );
            ir->AddPentPoint( 114 , 0.0227203202344919 , 0.4558252424577396 , 0.1738181457692562 , 0.1738181457692562 , 0.0003272598432739 );
            ir->AddPentPoint( 115 , 0.1738181457692562 , 0.4558252424577396 , 0.1738181457692562 , 0.1738181457692562 , 0.0003272598432739 );
            ir->AddPentPoint( 116 , 0.4558252424577396 , 0.0227203202344919 , 0.1738181457692562 , 0.1738181457692562 , 0.0003272598432739 );
            ir->AddPentPoint( 117 , 0.1738181457692562 , 0.0227203202344919 , 0.1738181457692562 , 0.1738181457692562 , 0.0003272598432739 );
            ir->AddPentPoint( 118 , 0.4558252424577396 , 0.1738181457692562 , 0.1738181457692562 , 0.1738181457692562 , 0.0003272598432739 );
            ir->AddPentPoint( 119 , 0.0227203202344919 , 0.1738181457692562 , 0.1738181457692562 , 0.1738181457692562 , 0.0003272598432739 );
            ir->AddPentPoint( 120 , 0.0000004574705963 , 0.4568047468962194 , 0.4568047468962194 , 0.0863895912663686 , 0.0000298974252054 );
            ir->AddPentPoint( 121 , 0.4568047468962194 , 0.0000004574705963 , 0.4568047468962194 , 0.0863895912663686 , 0.0000298974252054 );
            ir->AddPentPoint( 122 , 0.0000004574705963 , 0.0000004574705963 , 0.4568047468962194 , 0.0863895912663686 , 0.0000298974252054 );
            ir->AddPentPoint( 123 , 0.4568047468962194 , 0.4568047468962194 , 0.0000004574705963 , 0.0863895912663686 , 0.0000298974252054 );
            ir->AddPentPoint( 124 , 0.0000004574705963 , 0.4568047468962194 , 0.0000004574705963 , 0.0863895912663686 , 0.0000298974252054 );
            ir->AddPentPoint( 125 , 0.4568047468962194 , 0.0000004574705963 , 0.0000004574705963 , 0.0863895912663686 , 0.0000298974252054 );
            ir->AddPentPoint( 126 , 0.0000004574705963 , 0.4568047468962194 , 0.0863895912663686 , 0.4568047468962194 , 0.0000298974252054 );
            ir->AddPentPoint( 127 , 0.4568047468962194 , 0.0000004574705963 , 0.0863895912663686 , 0.4568047468962194 , 0.0000298974252054 );
            ir->AddPentPoint( 128 , 0.0000004574705963 , 0.0000004574705963 , 0.0863895912663686 , 0.4568047468962194 , 0.0000298974252054 );
            ir->AddPentPoint( 129 , 0.0000004574705963 , 0.0863895912663686 , 0.4568047468962194 , 0.4568047468962194 , 0.0000298974252054 );
            ir->AddPentPoint( 130 , 0.0863895912663686 , 0.0000004574705963 , 0.4568047468962194 , 0.4568047468962194 , 0.0000298974252054 );
            ir->AddPentPoint( 131 , 0.0000004574705963 , 0.0000004574705963 , 0.4568047468962194 , 0.4568047468962194 , 0.0000298974252054 );
            ir->AddPentPoint( 132 , 0.4568047468962194 , 0.0863895912663686 , 0.0000004574705963 , 0.4568047468962194 , 0.0000298974252054 );
            ir->AddPentPoint( 133 , 0.0000004574705963 , 0.0863895912663686 , 0.0000004574705963 , 0.4568047468962194 , 0.0000298974252054 );
            ir->AddPentPoint( 134 , 0.0863895912663686 , 0.4568047468962194 , 0.0000004574705963 , 0.4568047468962194 , 0.0000298974252054 );
            ir->AddPentPoint( 135 , 0.0000004574705963 , 0.4568047468962194 , 0.0000004574705963 , 0.4568047468962194 , 0.0000298974252054 );
            ir->AddPentPoint( 136 , 0.0863895912663686 , 0.0000004574705963 , 0.0000004574705963 , 0.4568047468962194 , 0.0000298974252054 );
            ir->AddPentPoint( 137 , 0.4568047468962194 , 0.0000004574705963 , 0.0000004574705963 , 0.4568047468962194 , 0.0000298974252054 );
            ir->AddPentPoint( 138 , 0.4568047468962194 , 0.4568047468962194 , 0.0863895912663686 , 0.0000004574705963 , 0.0000298974252054 );
            ir->AddPentPoint( 139 , 0.0000004574705963 , 0.4568047468962194 , 0.0863895912663686 , 0.0000004574705963 , 0.0000298974252054 );
            ir->AddPentPoint( 140 , 0.4568047468962194 , 0.0000004574705963 , 0.0863895912663686 , 0.0000004574705963 , 0.0000298974252054 );
            ir->AddPentPoint( 141 , 0.4568047468962194 , 0.0863895912663686 , 0.4568047468962194 , 0.0000004574705963 , 0.0000298974252054 );
            ir->AddPentPoint( 142 , 0.0000004574705963 , 0.0863895912663686 , 0.4568047468962194 , 0.0000004574705963 , 0.0000298974252054 );
            ir->AddPentPoint( 143 , 0.0863895912663686 , 0.4568047468962194 , 0.4568047468962194 , 0.0000004574705963 , 0.0000298974252054 );
            ir->AddPentPoint( 144 , 0.0000004574705963 , 0.4568047468962194 , 0.4568047468962194 , 0.0000004574705963 , 0.0000298974252054 );
            ir->AddPentPoint( 145 , 0.0863895912663686 , 0.0000004574705963 , 0.4568047468962194 , 0.0000004574705963 , 0.0000298974252054 );
            ir->AddPentPoint( 146 , 0.4568047468962194 , 0.0000004574705963 , 0.4568047468962194 , 0.0000004574705963 , 0.0000298974252054 );
            ir->AddPentPoint( 147 , 0.4568047468962194 , 0.0863895912663686 , 0.0000004574705963 , 0.0000004574705963 , 0.0000298974252054 );
            ir->AddPentPoint( 148 , 0.0863895912663686 , 0.4568047468962194 , 0.0000004574705963 , 0.0000004574705963 , 0.0000298974252054 );
            ir->AddPentPoint( 149 , 0.4568047468962194 , 0.4568047468962194 , 0.0000004574705963 , 0.0000004574705963 , 0.0000298974252054 );
            ir->AddPentPoint( 150 , 0.0262918390784996 , 0.1915132492211721 , 0.1915132492211721 , 0.5643898234006565 , 0.0001813711012859 );
            ir->AddPentPoint( 151 , 0.1915132492211721 , 0.0262918390784996 , 0.1915132492211721 , 0.5643898234006565 , 0.0001813711012859 );
            ir->AddPentPoint( 152 , 0.0262918390784996 , 0.0262918390784996 , 0.1915132492211721 , 0.5643898234006565 , 0.0001813711012859 );
            ir->AddPentPoint( 153 , 0.1915132492211721 , 0.1915132492211721 , 0.0262918390784996 , 0.5643898234006565 , 0.0001813711012859 );
            ir->AddPentPoint( 154 , 0.0262918390784996 , 0.1915132492211721 , 0.0262918390784996 , 0.5643898234006565 , 0.0001813711012859 );
            ir->AddPentPoint( 155 , 0.1915132492211721 , 0.0262918390784996 , 0.0262918390784996 , 0.5643898234006565 , 0.0001813711012859 );
            ir->AddPentPoint( 156 , 0.0262918390784996 , 0.1915132492211721 , 0.5643898234006565 , 0.1915132492211721 , 0.0001813711012859 );
            ir->AddPentPoint( 157 , 0.1915132492211721 , 0.0262918390784996 , 0.5643898234006565 , 0.1915132492211721 , 0.0001813711012859 );
            ir->AddPentPoint( 158 , 0.0262918390784996 , 0.0262918390784996 , 0.5643898234006565 , 0.1915132492211721 , 0.0001813711012859 );
            ir->AddPentPoint( 159 , 0.0262918390784996 , 0.5643898234006565 , 0.1915132492211721 , 0.1915132492211721 , 0.0001813711012859 );
            ir->AddPentPoint( 160 , 0.5643898234006565 , 0.0262918390784996 , 0.1915132492211721 , 0.1915132492211721 , 0.0001813711012859 );
            ir->AddPentPoint( 161 , 0.0262918390784996 , 0.0262918390784996 , 0.1915132492211721 , 0.1915132492211721 , 0.0001813711012859 );
            ir->AddPentPoint( 162 , 0.1915132492211721 , 0.5643898234006565 , 0.0262918390784996 , 0.1915132492211721 , 0.0001813711012859 );
            ir->AddPentPoint( 163 , 0.0262918390784996 , 0.5643898234006565 , 0.0262918390784996 , 0.1915132492211721 , 0.0001813711012859 );
            ir->AddPentPoint( 164 , 0.5643898234006565 , 0.1915132492211721 , 0.0262918390784996 , 0.1915132492211721 , 0.0001813711012859 );
            ir->AddPentPoint( 165 , 0.0262918390784996 , 0.1915132492211721 , 0.0262918390784996 , 0.1915132492211721 , 0.0001813711012859 );
            ir->AddPentPoint( 166 , 0.5643898234006565 , 0.0262918390784996 , 0.0262918390784996 , 0.1915132492211721 , 0.0001813711012859 );
            ir->AddPentPoint( 167 , 0.1915132492211721 , 0.0262918390784996 , 0.0262918390784996 , 0.1915132492211721 , 0.0001813711012859 );
            ir->AddPentPoint( 168 , 0.1915132492211721 , 0.1915132492211721 , 0.5643898234006565 , 0.0262918390784996 , 0.0001813711012859 );
            ir->AddPentPoint( 169 , 0.0262918390784996 , 0.1915132492211721 , 0.5643898234006565 , 0.0262918390784996 , 0.0001813711012859 );
            ir->AddPentPoint( 170 , 0.1915132492211721 , 0.0262918390784996 , 0.5643898234006565 , 0.0262918390784996 , 0.0001813711012859 );
            ir->AddPentPoint( 171 , 0.1915132492211721 , 0.5643898234006565 , 0.1915132492211721 , 0.0262918390784996 , 0.0001813711012859 );
            ir->AddPentPoint( 172 , 0.0262918390784996 , 0.5643898234006565 , 0.1915132492211721 , 0.0262918390784996 , 0.0001813711012859 );
            ir->AddPentPoint( 173 , 0.5643898234006565 , 0.1915132492211721 , 0.1915132492211721 , 0.0262918390784996 , 0.0001813711012859 );
            ir->AddPentPoint( 174 , 0.0262918390784996 , 0.1915132492211721 , 0.1915132492211721 , 0.0262918390784996 , 0.0001813711012859 );
            ir->AddPentPoint( 175 , 0.5643898234006565 , 0.0262918390784996 , 0.1915132492211721 , 0.0262918390784996 , 0.0001813711012859 );
            ir->AddPentPoint( 176 , 0.1915132492211721 , 0.0262918390784996 , 0.1915132492211721 , 0.0262918390784996 , 0.0001813711012859 );
            ir->AddPentPoint( 177 , 0.1915132492211721 , 0.5643898234006565 , 0.0262918390784996 , 0.0262918390784996 , 0.0001813711012859 );
            ir->AddPentPoint( 178 , 0.5643898234006565 , 0.1915132492211721 , 0.0262918390784996 , 0.0262918390784996 , 0.0001813711012859 );
            ir->AddPentPoint( 179 , 0.1915132492211721 , 0.1915132492211721 , 0.0262918390784996 , 0.0262918390784996 , 0.0001813711012859 );
            ir->AddPentPoint( 180 , 0.1057514121324120 , 0.3860828233331977 , 0.3860828233331977 , 0.0163315290687808 , 0.0002335242663709 );
            ir->AddPentPoint( 181 , 0.3860828233331977 , 0.1057514121324120 , 0.3860828233331977 , 0.0163315290687808 , 0.0002335242663709 );
            ir->AddPentPoint( 182 , 0.1057514121324120 , 0.1057514121324120 , 0.3860828233331977 , 0.0163315290687808 , 0.0002335242663709 );
            ir->AddPentPoint( 183 , 0.3860828233331977 , 0.3860828233331977 , 0.1057514121324120 , 0.0163315290687808 , 0.0002335242663709 );
            ir->AddPentPoint( 184 , 0.1057514121324120 , 0.3860828233331977 , 0.1057514121324120 , 0.0163315290687808 , 0.0002335242663709 );
            ir->AddPentPoint( 185 , 0.3860828233331977 , 0.1057514121324120 , 0.1057514121324120 , 0.0163315290687808 , 0.0002335242663709 );
            ir->AddPentPoint( 186 , 0.1057514121324120 , 0.3860828233331977 , 0.0163315290687808 , 0.3860828233331977 , 0.0002335242663709 );
            ir->AddPentPoint( 187 , 0.3860828233331977 , 0.1057514121324120 , 0.0163315290687808 , 0.3860828233331977 , 0.0002335242663709 );
            ir->AddPentPoint( 188 , 0.1057514121324120 , 0.1057514121324120 , 0.0163315290687808 , 0.3860828233331977 , 0.0002335242663709 );
            ir->AddPentPoint( 189 , 0.1057514121324120 , 0.0163315290687808 , 0.3860828233331977 , 0.3860828233331977 , 0.0002335242663709 );
            ir->AddPentPoint( 190 , 0.0163315290687808 , 0.1057514121324120 , 0.3860828233331977 , 0.3860828233331977 , 0.0002335242663709 );
            ir->AddPentPoint( 191 , 0.1057514121324120 , 0.1057514121324120 , 0.3860828233331977 , 0.3860828233331977 , 0.0002335242663709 );
            ir->AddPentPoint( 192 , 0.3860828233331977 , 0.0163315290687808 , 0.1057514121324120 , 0.3860828233331977 , 0.0002335242663709 );
            ir->AddPentPoint( 193 , 0.1057514121324120 , 0.0163315290687808 , 0.1057514121324120 , 0.3860828233331977 , 0.0002335242663709 );
            ir->AddPentPoint( 194 , 0.0163315290687808 , 0.3860828233331977 , 0.1057514121324120 , 0.3860828233331977 , 0.0002335242663709 );
            ir->AddPentPoint( 195 , 0.1057514121324120 , 0.3860828233331977 , 0.1057514121324120 , 0.3860828233331977 , 0.0002335242663709 );
            ir->AddPentPoint( 196 , 0.0163315290687808 , 0.1057514121324120 , 0.1057514121324120 , 0.3860828233331977 , 0.0002335242663709 );
            ir->AddPentPoint( 197 , 0.3860828233331977 , 0.1057514121324120 , 0.1057514121324120 , 0.3860828233331977 , 0.0002335242663709 );
            ir->AddPentPoint( 198 , 0.3860828233331977 , 0.3860828233331977 , 0.0163315290687808 , 0.1057514121324120 , 0.0002335242663709 );
            ir->AddPentPoint( 199 , 0.1057514121324120 , 0.3860828233331977 , 0.0163315290687808 , 0.1057514121324120 , 0.0002335242663709 );
            ir->AddPentPoint( 200 , 0.3860828233331977 , 0.1057514121324120 , 0.0163315290687808 , 0.1057514121324120 , 0.0002335242663709 );
            ir->AddPentPoint( 201 , 0.3860828233331977 , 0.0163315290687808 , 0.3860828233331977 , 0.1057514121324120 , 0.0002335242663709 );
            ir->AddPentPoint( 202 , 0.1057514121324120 , 0.0163315290687808 , 0.3860828233331977 , 0.1057514121324120 , 0.0002335242663709 );
            ir->AddPentPoint( 203 , 0.0163315290687808 , 0.3860828233331977 , 0.3860828233331977 , 0.1057514121324120 , 0.0002335242663709 );
            ir->AddPentPoint( 204 , 0.1057514121324120 , 0.3860828233331977 , 0.3860828233331977 , 0.1057514121324120 , 0.0002335242663709 );
            ir->AddPentPoint( 205 , 0.0163315290687808 , 0.1057514121324120 , 0.3860828233331977 , 0.1057514121324120 , 0.0002335242663709 );
            ir->AddPentPoint( 206 , 0.3860828233331977 , 0.1057514121324120 , 0.3860828233331977 , 0.1057514121324120 , 0.0002335242663709 );
            ir->AddPentPoint( 207 , 0.3860828233331977 , 0.0163315290687808 , 0.1057514121324120 , 0.1057514121324120 , 0.0002335242663709 );
            ir->AddPentPoint( 208 , 0.0163315290687808 , 0.3860828233331977 , 0.1057514121324120 , 0.1057514121324120 , 0.0002335242663709 );
            ir->AddPentPoint( 209 , 0.3860828233331977 , 0.3860828233331977 , 0.1057514121324120 , 0.1057514121324120 , 0.0002335242663709 );
            
            return ir;
            
        case 11:  // 281 points - degree 11 --
            PentatopeIntRules[11] = ir = new IntegrationRule(281);
            ir->AddPentPoint( 0 , 0.2000000000000000 , 0.2000000000000000 , 0.2000000000000000 , 0.2000000000000000 , 0.0007490771117245 );
            ir->AddPentPoint( 1 , 0.0273293708859630 , 0.0273293708859630 , 0.0273293708859630 , 0.8906825164561480 , 0.0000249116974622 );
            ir->AddPentPoint( 2 , 0.0273293708859630 , 0.0273293708859630 , 0.8906825164561480 , 0.0273293708859630 , 0.0000249116974622 );
            ir->AddPentPoint( 3 , 0.0273293708859630 , 0.8906825164561480 , 0.0273293708859630 , 0.0273293708859630 , 0.0000249116974622 );
            ir->AddPentPoint( 4 , 0.8906825164561480 , 0.0273293708859630 , 0.0273293708859630 , 0.0273293708859630 , 0.0000249116974622 );
            ir->AddPentPoint( 5 , 0.0273293708859630 , 0.0273293708859630 , 0.0273293708859630 , 0.0273293708859630 , 0.0000249116974622 );
            ir->AddPentPoint( 6 , 0.0672416510291032 , 0.0672416510291032 , 0.0672416510291032 , 0.7310333958835870 , 0.0001982426851854 );
            ir->AddPentPoint( 7 , 0.0672416510291032 , 0.0672416510291032 , 0.7310333958835870 , 0.0672416510291032 , 0.0001982426851854 );
            ir->AddPentPoint( 8 , 0.0672416510291032 , 0.7310333958835870 , 0.0672416510291032 , 0.0672416510291032 , 0.0001982426851854 );
            ir->AddPentPoint( 9 , 0.7310333958835870 , 0.0672416510291032 , 0.0672416510291032 , 0.0672416510291032 , 0.0001982426851854 );
            ir->AddPentPoint( 10 , 0.0672416510291032 , 0.0672416510291032 , 0.0672416510291032 , 0.0672416510291032 , 0.0001982426851854 );
            ir->AddPentPoint( 11 , 0.1222724170020238 , 0.1222724170020238 , 0.1222724170020238 , 0.5109103319919047 , 0.0004260569890749 );
            ir->AddPentPoint( 12 , 0.1222724170020238 , 0.1222724170020238 , 0.5109103319919047 , 0.1222724170020238 , 0.0004260569890749 );
            ir->AddPentPoint( 13 , 0.1222724170020238 , 0.5109103319919047 , 0.1222724170020238 , 0.1222724170020238 , 0.0004260569890749 );
            ir->AddPentPoint( 14 , 0.5109103319919047 , 0.1222724170020238 , 0.1222724170020238 , 0.1222724170020238 , 0.0004260569890749 );
            ir->AddPentPoint( 15 , 0.1222724170020238 , 0.1222724170020238 , 0.1222724170020238 , 0.1222724170020238 , 0.0004260569890749 );
            ir->AddPentPoint( 16 , 0.2499981427637660 , 0.2499981427637660 , 0.2499981427637660 , 0.0000074289449360 , 0.0001099487523676 );
            ir->AddPentPoint( 17 , 0.2499981427637660 , 0.2499981427637660 , 0.0000074289449360 , 0.2499981427637660 , 0.0001099487523676 );
            ir->AddPentPoint( 18 , 0.2499981427637660 , 0.0000074289449360 , 0.2499981427637660 , 0.2499981427637660 , 0.0001099487523676 );
            ir->AddPentPoint( 19 , 0.0000074289449360 , 0.2499981427637660 , 0.2499981427637660 , 0.2499981427637660 , 0.0001099487523676 );
            ir->AddPentPoint( 20 , 0.2499981427637660 , 0.2499981427637660 , 0.2499981427637660 , 0.2499981427637660 , 0.0001099487523676 );
            ir->AddPentPoint( 21 , 0.0000000013178651 , 0.0000000013178651 , 0.4999999980232024 , 0.4999999980232024 , 0.0000056551263518 );
            ir->AddPentPoint( 22 , 0.0000000013178651 , 0.4999999980232024 , 0.0000000013178651 , 0.4999999980232024 , 0.0000056551263518 );
            ir->AddPentPoint( 23 , 0.4999999980232024 , 0.0000000013178651 , 0.0000000013178651 , 0.4999999980232024 , 0.0000056551263518 );
            ir->AddPentPoint( 24 , 0.0000000013178651 , 0.0000000013178651 , 0.0000000013178651 , 0.4999999980232024 , 0.0000056551263518 );
            ir->AddPentPoint( 25 , 0.0000000013178651 , 0.4999999980232024 , 0.4999999980232024 , 0.0000000013178651 , 0.0000056551263518 );
            ir->AddPentPoint( 26 , 0.4999999980232024 , 0.0000000013178651 , 0.4999999980232024 , 0.0000000013178651 , 0.0000056551263518 );
            ir->AddPentPoint( 27 , 0.0000000013178651 , 0.0000000013178651 , 0.4999999980232024 , 0.0000000013178651 , 0.0000056551263518 );
            ir->AddPentPoint( 28 , 0.4999999980232024 , 0.4999999980232024 , 0.0000000013178651 , 0.0000000013178651 , 0.0000056551263518 );
            ir->AddPentPoint( 29 , 0.0000000013178651 , 0.4999999980232024 , 0.0000000013178651 , 0.0000000013178651 , 0.0000056551263518 );
            ir->AddPentPoint( 30 , 0.4999999980232024 , 0.0000000013178651 , 0.0000000013178651 , 0.0000000013178651 , 0.0000056551263518 );
            ir->AddPentPoint( 31 , 0.0932486705850966 , 0.0932486705850966 , 0.3601269941223551 , 0.3601269941223551 , 0.0003732112383759 );
            ir->AddPentPoint( 32 , 0.0932486705850966 , 0.3601269941223551 , 0.0932486705850966 , 0.3601269941223551 , 0.0003732112383759 );
            ir->AddPentPoint( 33 , 0.3601269941223551 , 0.0932486705850966 , 0.0932486705850966 , 0.3601269941223551 , 0.0003732112383759 );
            ir->AddPentPoint( 34 , 0.0932486705850966 , 0.0932486705850966 , 0.0932486705850966 , 0.3601269941223551 , 0.0003732112383759 );
            ir->AddPentPoint( 35 , 0.0932486705850966 , 0.3601269941223551 , 0.3601269941223551 , 0.0932486705850966 , 0.0003732112383759 );
            ir->AddPentPoint( 36 , 0.3601269941223551 , 0.0932486705850966 , 0.3601269941223551 , 0.0932486705850966 , 0.0003732112383759 );
            ir->AddPentPoint( 37 , 0.0932486705850966 , 0.0932486705850966 , 0.3601269941223551 , 0.0932486705850966 , 0.0003732112383759 );
            ir->AddPentPoint( 38 , 0.3601269941223551 , 0.3601269941223551 , 0.0932486705850966 , 0.0932486705850966 , 0.0003732112383759 );
            ir->AddPentPoint( 39 , 0.0932486705850966 , 0.3601269941223551 , 0.0932486705850966 , 0.0932486705850966 , 0.0003732112383759 );
            ir->AddPentPoint( 40 , 0.3601269941223551 , 0.0932486705850966 , 0.0932486705850966 , 0.0932486705850966 , 0.0003732112383759 );
            ir->AddPentPoint( 41 , 0.2706108009024531 , 0.2706108009024531 , 0.0940837986463203 , 0.0940837986463203 , 0.0004865726490793 );
            ir->AddPentPoint( 42 , 0.2706108009024531 , 0.0940837986463203 , 0.2706108009024531 , 0.0940837986463203 , 0.0004865726490793 );
            ir->AddPentPoint( 43 , 0.0940837986463203 , 0.2706108009024531 , 0.2706108009024531 , 0.0940837986463203 , 0.0004865726490793 );
            ir->AddPentPoint( 44 , 0.2706108009024531 , 0.2706108009024531 , 0.2706108009024531 , 0.0940837986463203 , 0.0004865726490793 );
            ir->AddPentPoint( 45 , 0.2706108009024531 , 0.0940837986463203 , 0.0940837986463203 , 0.2706108009024531 , 0.0004865726490793 );
            ir->AddPentPoint( 46 , 0.0940837986463203 , 0.2706108009024531 , 0.0940837986463203 , 0.2706108009024531 , 0.0004865726490793 );
            ir->AddPentPoint( 47 , 0.2706108009024531 , 0.2706108009024531 , 0.0940837986463203 , 0.2706108009024531 , 0.0004865726490793 );
            ir->AddPentPoint( 48 , 0.0940837986463203 , 0.0940837986463203 , 0.2706108009024531 , 0.2706108009024531 , 0.0004865726490793 );
            ir->AddPentPoint( 49 , 0.2706108009024531 , 0.0940837986463203 , 0.2706108009024531 , 0.2706108009024531 , 0.0004865726490793 );
            ir->AddPentPoint( 50 , 0.0940837986463203 , 0.2706108009024531 , 0.2706108009024531 , 0.2706108009024531 , 0.0004865726490793 );
            ir->AddPentPoint( 51 , 0.0422640510205627 , 0.0422640510205627 , 0.2717998707940745 , 0.6014079761442372 , 0.0001808462102100 );
            ir->AddPentPoint( 52 , 0.0422640510205627 , 0.2717998707940745 , 0.0422640510205627 , 0.6014079761442372 , 0.0001808462102100 );
            ir->AddPentPoint( 53 , 0.2717998707940745 , 0.0422640510205627 , 0.0422640510205627 , 0.6014079761442372 , 0.0001808462102100 );
            ir->AddPentPoint( 54 , 0.0422640510205627 , 0.0422640510205627 , 0.0422640510205627 , 0.6014079761442372 , 0.0001808462102100 );
            ir->AddPentPoint( 55 , 0.0422640510205627 , 0.0422640510205627 , 0.6014079761442372 , 0.2717998707940745 , 0.0001808462102100 );
            ir->AddPentPoint( 56 , 0.0422640510205627 , 0.6014079761442372 , 0.0422640510205627 , 0.2717998707940745 , 0.0001808462102100 );
            ir->AddPentPoint( 57 , 0.6014079761442372 , 0.0422640510205627 , 0.0422640510205627 , 0.2717998707940745 , 0.0001808462102100 );
            ir->AddPentPoint( 58 , 0.0422640510205627 , 0.0422640510205627 , 0.0422640510205627 , 0.2717998707940745 , 0.0001808462102100 );
            ir->AddPentPoint( 59 , 0.0422640510205627 , 0.2717998707940745 , 0.6014079761442372 , 0.0422640510205627 , 0.0001808462102100 );
            ir->AddPentPoint( 60 , 0.2717998707940745 , 0.0422640510205627 , 0.6014079761442372 , 0.0422640510205627 , 0.0001808462102100 );
            ir->AddPentPoint( 61 , 0.0422640510205627 , 0.0422640510205627 , 0.6014079761442372 , 0.0422640510205627 , 0.0001808462102100 );
            ir->AddPentPoint( 62 , 0.0422640510205627 , 0.6014079761442372 , 0.2717998707940745 , 0.0422640510205627 , 0.0001808462102100 );
            ir->AddPentPoint( 63 , 0.6014079761442372 , 0.0422640510205627 , 0.2717998707940745 , 0.0422640510205627 , 0.0001808462102100 );
            ir->AddPentPoint( 64 , 0.0422640510205627 , 0.0422640510205627 , 0.2717998707940745 , 0.0422640510205627 , 0.0001808462102100 );
            ir->AddPentPoint( 65 , 0.2717998707940745 , 0.6014079761442372 , 0.0422640510205627 , 0.0422640510205627 , 0.0001808462102100 );
            ir->AddPentPoint( 66 , 0.0422640510205627 , 0.6014079761442372 , 0.0422640510205627 , 0.0422640510205627 , 0.0001808462102100 );
            ir->AddPentPoint( 67 , 0.6014079761442372 , 0.2717998707940745 , 0.0422640510205627 , 0.0422640510205627 , 0.0001808462102100 );
            ir->AddPentPoint( 68 , 0.0422640510205627 , 0.2717998707940745 , 0.0422640510205627 , 0.0422640510205627 , 0.0001808462102100 );
            ir->AddPentPoint( 69 , 0.6014079761442372 , 0.0422640510205627 , 0.0422640510205627 , 0.0422640510205627 , 0.0001808462102100 );
            ir->AddPentPoint( 70 , 0.2717998707940745 , 0.0422640510205627 , 0.0422640510205627 , 0.0422640510205627 , 0.0001808462102100 );
            ir->AddPentPoint( 71 , 0.1168661129986462 , 0.1168661129986462 , 0.0000000000000002 , 0.6494016610040614 , 0.0000825913954391 );
            ir->AddPentPoint( 72 , 0.1168661129986462 , 0.0000000000000002 , 0.1168661129986462 , 0.6494016610040614 , 0.0000825913954391 );
            ir->AddPentPoint( 73 , 0.0000000000000002 , 0.1168661129986462 , 0.1168661129986462 , 0.6494016610040614 , 0.0000825913954391 );
            ir->AddPentPoint( 74 , 0.1168661129986462 , 0.1168661129986462 , 0.1168661129986462 , 0.6494016610040614 , 0.0000825913954391 );
            ir->AddPentPoint( 75 , 0.1168661129986462 , 0.1168661129986462 , 0.6494016610040614 , 0.0000000000000002 , 0.0000825913954391 );
            ir->AddPentPoint( 76 , 0.1168661129986462 , 0.6494016610040614 , 0.1168661129986462 , 0.0000000000000002 , 0.0000825913954391 );
            ir->AddPentPoint( 77 , 0.6494016610040614 , 0.1168661129986462 , 0.1168661129986462 , 0.0000000000000002 , 0.0000825913954391 );
            ir->AddPentPoint( 78 , 0.1168661129986462 , 0.1168661129986462 , 0.1168661129986462 , 0.0000000000000002 , 0.0000825913954391 );
            ir->AddPentPoint( 79 , 0.1168661129986462 , 0.0000000000000002 , 0.6494016610040614 , 0.1168661129986462 , 0.0000825913954391 );
            ir->AddPentPoint( 80 , 0.0000000000000002 , 0.1168661129986462 , 0.6494016610040614 , 0.1168661129986462 , 0.0000825913954391 );
            ir->AddPentPoint( 81 , 0.1168661129986462 , 0.1168661129986462 , 0.6494016610040614 , 0.1168661129986462 , 0.0000825913954391 );
            ir->AddPentPoint( 82 , 0.1168661129986462 , 0.6494016610040614 , 0.0000000000000002 , 0.1168661129986462 , 0.0000825913954391 );
            ir->AddPentPoint( 83 , 0.6494016610040614 , 0.1168661129986462 , 0.0000000000000002 , 0.1168661129986462 , 0.0000825913954391 );
            ir->AddPentPoint( 84 , 0.1168661129986462 , 0.1168661129986462 , 0.0000000000000002 , 0.1168661129986462 , 0.0000825913954391 );
            ir->AddPentPoint( 85 , 0.0000000000000002 , 0.6494016610040614 , 0.1168661129986462 , 0.1168661129986462 , 0.0000825913954391 );
            ir->AddPentPoint( 86 , 0.1168661129986462 , 0.6494016610040614 , 0.1168661129986462 , 0.1168661129986462 , 0.0000825913954391 );
            ir->AddPentPoint( 87 , 0.6494016610040614 , 0.0000000000000002 , 0.1168661129986462 , 0.1168661129986462 , 0.0000825913954391 );
            ir->AddPentPoint( 88 , 0.1168661129986462 , 0.0000000000000002 , 0.1168661129986462 , 0.1168661129986462 , 0.0000825913954391 );
            ir->AddPentPoint( 89 , 0.6494016610040614 , 0.1168661129986462 , 0.1168661129986462 , 0.1168661129986462 , 0.0000825913954391 );
            ir->AddPentPoint( 90 , 0.0000000000000002 , 0.1168661129986462 , 0.1168661129986462 , 0.1168661129986462 , 0.0000825913954391 );
            ir->AddPentPoint( 91 , 0.1843113539422363 , 0.1843113539422363 , 0.0366827746698410 , 0.4103831635034501 , 0.0003655018208071 );
            ir->AddPentPoint( 92 , 0.1843113539422363 , 0.0366827746698410 , 0.1843113539422363 , 0.4103831635034501 , 0.0003655018208071 );
            ir->AddPentPoint( 93 , 0.0366827746698410 , 0.1843113539422363 , 0.1843113539422363 , 0.4103831635034501 , 0.0003655018208071 );
            ir->AddPentPoint( 94 , 0.1843113539422363 , 0.1843113539422363 , 0.1843113539422363 , 0.4103831635034501 , 0.0003655018208071 );
            ir->AddPentPoint( 95 , 0.1843113539422363 , 0.1843113539422363 , 0.4103831635034501 , 0.0366827746698410 , 0.0003655018208071 );
            ir->AddPentPoint( 96 , 0.1843113539422363 , 0.4103831635034501 , 0.1843113539422363 , 0.0366827746698410 , 0.0003655018208071 );
            ir->AddPentPoint( 97 , 0.4103831635034501 , 0.1843113539422363 , 0.1843113539422363 , 0.0366827746698410 , 0.0003655018208071 );
            ir->AddPentPoint( 98 , 0.1843113539422363 , 0.1843113539422363 , 0.1843113539422363 , 0.0366827746698410 , 0.0003655018208071 );
            ir->AddPentPoint( 99 , 0.1843113539422363 , 0.0366827746698410 , 0.4103831635034501 , 0.1843113539422363 , 0.0003655018208071 );
            ir->AddPentPoint( 100 , 0.0366827746698410 , 0.1843113539422363 , 0.4103831635034501 , 0.1843113539422363 , 0.0003655018208071 );
            ir->AddPentPoint( 101 , 0.1843113539422363 , 0.1843113539422363 , 0.4103831635034501 , 0.1843113539422363 , 0.0003655018208071 );
            ir->AddPentPoint( 102 , 0.1843113539422363 , 0.4103831635034501 , 0.0366827746698410 , 0.1843113539422363 , 0.0003655018208071 );
            ir->AddPentPoint( 103 , 0.4103831635034501 , 0.1843113539422363 , 0.0366827746698410 , 0.1843113539422363 , 0.0003655018208071 );
            ir->AddPentPoint( 104 , 0.1843113539422363 , 0.1843113539422363 , 0.0366827746698410 , 0.1843113539422363 , 0.0003655018208071 );
            ir->AddPentPoint( 105 , 0.0366827746698410 , 0.4103831635034501 , 0.1843113539422363 , 0.1843113539422363 , 0.0003655018208071 );
            ir->AddPentPoint( 106 , 0.1843113539422363 , 0.4103831635034501 , 0.1843113539422363 , 0.1843113539422363 , 0.0003655018208071 );
            ir->AddPentPoint( 107 , 0.4103831635034501 , 0.0366827746698410 , 0.1843113539422363 , 0.1843113539422363 , 0.0003655018208071 );
            ir->AddPentPoint( 108 , 0.1843113539422363 , 0.0366827746698410 , 0.1843113539422363 , 0.1843113539422363 , 0.0003655018208071 );
            ir->AddPentPoint( 109 , 0.4103831635034501 , 0.1843113539422363 , 0.1843113539422363 , 0.1843113539422363 , 0.0003655018208071 );
            ir->AddPentPoint( 110 , 0.0366827746698410 , 0.1843113539422363 , 0.1843113539422363 , 0.1843113539422363 , 0.0003655018208071 );
            ir->AddPentPoint( 111 , 0.3120944268112401 , 0.3120944268112401 , 0.0085885637602618 , 0.0551281558060177 , 0.0001212798478515 );
            ir->AddPentPoint( 112 , 0.3120944268112401 , 0.0085885637602618 , 0.3120944268112401 , 0.0551281558060177 , 0.0001212798478515 );
            ir->AddPentPoint( 113 , 0.0085885637602618 , 0.3120944268112401 , 0.3120944268112401 , 0.0551281558060177 , 0.0001212798478515 );
            ir->AddPentPoint( 114 , 0.3120944268112401 , 0.3120944268112401 , 0.3120944268112401 , 0.0551281558060177 , 0.0001212798478515 );
            ir->AddPentPoint( 115 , 0.3120944268112401 , 0.3120944268112401 , 0.0551281558060177 , 0.0085885637602618 , 0.0001212798478515 );
            ir->AddPentPoint( 116 , 0.3120944268112401 , 0.0551281558060177 , 0.3120944268112401 , 0.0085885637602618 , 0.0001212798478515 );
            ir->AddPentPoint( 117 , 0.0551281558060177 , 0.3120944268112401 , 0.3120944268112401 , 0.0085885637602618 , 0.0001212798478515 );
            ir->AddPentPoint( 118 , 0.3120944268112401 , 0.3120944268112401 , 0.3120944268112401 , 0.0085885637602618 , 0.0001212798478515 );
            ir->AddPentPoint( 119 , 0.3120944268112401 , 0.0085885637602618 , 0.0551281558060177 , 0.3120944268112401 , 0.0001212798478515 );
            ir->AddPentPoint( 120 , 0.0085885637602618 , 0.3120944268112401 , 0.0551281558060177 , 0.3120944268112401 , 0.0001212798478515 );
            ir->AddPentPoint( 121 , 0.3120944268112401 , 0.3120944268112401 , 0.0551281558060177 , 0.3120944268112401 , 0.0001212798478515 );
            ir->AddPentPoint( 122 , 0.3120944268112401 , 0.0551281558060177 , 0.0085885637602618 , 0.3120944268112401 , 0.0001212798478515 );
            ir->AddPentPoint( 123 , 0.0551281558060177 , 0.3120944268112401 , 0.0085885637602618 , 0.3120944268112401 , 0.0001212798478515 );
            ir->AddPentPoint( 124 , 0.3120944268112401 , 0.3120944268112401 , 0.0085885637602618 , 0.3120944268112401 , 0.0001212798478515 );
            ir->AddPentPoint( 125 , 0.0085885637602618 , 0.0551281558060177 , 0.3120944268112401 , 0.3120944268112401 , 0.0001212798478515 );
            ir->AddPentPoint( 126 , 0.3120944268112401 , 0.0551281558060177 , 0.3120944268112401 , 0.3120944268112401 , 0.0001212798478515 );
            ir->AddPentPoint( 127 , 0.0551281558060177 , 0.0085885637602618 , 0.3120944268112401 , 0.3120944268112401 , 0.0001212798478515 );
            ir->AddPentPoint( 128 , 0.3120944268112401 , 0.0085885637602618 , 0.3120944268112401 , 0.3120944268112401 , 0.0001212798478515 );
            ir->AddPentPoint( 129 , 0.0551281558060177 , 0.3120944268112401 , 0.3120944268112401 , 0.3120944268112401 , 0.0001212798478515 );
            ir->AddPentPoint( 130 , 0.0085885637602618 , 0.3120944268112401 , 0.3120944268112401 , 0.3120944268112401 , 0.0001212798478515 );
            ir->AddPentPoint( 131 , 0.0046052778761422 , 0.4407307016284600 , 0.4407307016284600 , 0.1093280409907956 , 0.0000363112853415 );
            ir->AddPentPoint( 132 , 0.4407307016284600 , 0.0046052778761422 , 0.4407307016284600 , 0.1093280409907956 , 0.0000363112853415 );
            ir->AddPentPoint( 133 , 0.0046052778761422 , 0.0046052778761422 , 0.4407307016284600 , 0.1093280409907956 , 0.0000363112853415 );
            ir->AddPentPoint( 134 , 0.4407307016284600 , 0.4407307016284600 , 0.0046052778761422 , 0.1093280409907956 , 0.0000363112853415 );
            ir->AddPentPoint( 135 , 0.0046052778761422 , 0.4407307016284600 , 0.0046052778761422 , 0.1093280409907956 , 0.0000363112853415 );
            ir->AddPentPoint( 136 , 0.4407307016284600 , 0.0046052778761422 , 0.0046052778761422 , 0.1093280409907956 , 0.0000363112853415 );
            ir->AddPentPoint( 137 , 0.0046052778761422 , 0.4407307016284600 , 0.1093280409907956 , 0.4407307016284600 , 0.0000363112853415 );
            ir->AddPentPoint( 138 , 0.4407307016284600 , 0.0046052778761422 , 0.1093280409907956 , 0.4407307016284600 , 0.0000363112853415 );
            ir->AddPentPoint( 139 , 0.0046052778761422 , 0.0046052778761422 , 0.1093280409907956 , 0.4407307016284600 , 0.0000363112853415 );
            ir->AddPentPoint( 140 , 0.0046052778761422 , 0.1093280409907956 , 0.4407307016284600 , 0.4407307016284600 , 0.0000363112853415 );
            ir->AddPentPoint( 141 , 0.1093280409907956 , 0.0046052778761422 , 0.4407307016284600 , 0.4407307016284600 , 0.0000363112853415 );
            ir->AddPentPoint( 142 , 0.0046052778761422 , 0.0046052778761422 , 0.4407307016284600 , 0.4407307016284600 , 0.0000363112853415 );
            ir->AddPentPoint( 143 , 0.4407307016284600 , 0.1093280409907956 , 0.0046052778761422 , 0.4407307016284600 , 0.0000363112853415 );
            ir->AddPentPoint( 144 , 0.0046052778761422 , 0.1093280409907956 , 0.0046052778761422 , 0.4407307016284600 , 0.0000363112853415 );
            ir->AddPentPoint( 145 , 0.1093280409907956 , 0.4407307016284600 , 0.0046052778761422 , 0.4407307016284600 , 0.0000363112853415 );
            ir->AddPentPoint( 146 , 0.0046052778761422 , 0.4407307016284600 , 0.0046052778761422 , 0.4407307016284600 , 0.0000363112853415 );
            ir->AddPentPoint( 147 , 0.1093280409907956 , 0.0046052778761422 , 0.0046052778761422 , 0.4407307016284600 , 0.0000363112853415 );
            ir->AddPentPoint( 148 , 0.4407307016284600 , 0.0046052778761422 , 0.0046052778761422 , 0.4407307016284600 , 0.0000363112853415 );
            ir->AddPentPoint( 149 , 0.4407307016284600 , 0.4407307016284600 , 0.1093280409907956 , 0.0046052778761422 , 0.0000363112853415 );
            ir->AddPentPoint( 150 , 0.0046052778761422 , 0.4407307016284600 , 0.1093280409907956 , 0.0046052778761422 , 0.0000363112853415 );
            ir->AddPentPoint( 151 , 0.4407307016284600 , 0.0046052778761422 , 0.1093280409907956 , 0.0046052778761422 , 0.0000363112853415 );
            ir->AddPentPoint( 152 , 0.4407307016284600 , 0.1093280409907956 , 0.4407307016284600 , 0.0046052778761422 , 0.0000363112853415 );
            ir->AddPentPoint( 153 , 0.0046052778761422 , 0.1093280409907956 , 0.4407307016284600 , 0.0046052778761422 , 0.0000363112853415 );
            ir->AddPentPoint( 154 , 0.1093280409907956 , 0.4407307016284600 , 0.4407307016284600 , 0.0046052778761422 , 0.0000363112853415 );
            ir->AddPentPoint( 155 , 0.0046052778761422 , 0.4407307016284600 , 0.4407307016284600 , 0.0046052778761422 , 0.0000363112853415 );
            ir->AddPentPoint( 156 , 0.1093280409907956 , 0.0046052778761422 , 0.4407307016284600 , 0.0046052778761422 , 0.0000363112853415 );
            ir->AddPentPoint( 157 , 0.4407307016284600 , 0.0046052778761422 , 0.4407307016284600 , 0.0046052778761422 , 0.0000363112853415 );
            ir->AddPentPoint( 158 , 0.4407307016284600 , 0.1093280409907956 , 0.0046052778761422 , 0.0046052778761422 , 0.0000363112853415 );
            ir->AddPentPoint( 159 , 0.1093280409907956 , 0.4407307016284600 , 0.0046052778761422 , 0.0046052778761422 , 0.0000363112853415 );
            ir->AddPentPoint( 160 , 0.4407307016284600 , 0.4407307016284600 , 0.0046052778761422 , 0.0046052778761422 , 0.0000363112853415 );
            ir->AddPentPoint( 161 , 0.0339903200688579 , 0.2021249233072104 , 0.2021249233072104 , 0.5277695132478635 , 0.0002317033311045 );
            ir->AddPentPoint( 162 , 0.2021249233072104 , 0.0339903200688579 , 0.2021249233072104 , 0.5277695132478635 , 0.0002317033311045 );
            ir->AddPentPoint( 163 , 0.0339903200688579 , 0.0339903200688579 , 0.2021249233072104 , 0.5277695132478635 , 0.0002317033311045 );
            ir->AddPentPoint( 164 , 0.2021249233072104 , 0.2021249233072104 , 0.0339903200688579 , 0.5277695132478635 , 0.0002317033311045 );
            ir->AddPentPoint( 165 , 0.0339903200688579 , 0.2021249233072104 , 0.0339903200688579 , 0.5277695132478635 , 0.0002317033311045 );
            ir->AddPentPoint( 166 , 0.2021249233072104 , 0.0339903200688579 , 0.0339903200688579 , 0.5277695132478635 , 0.0002317033311045 );
            ir->AddPentPoint( 167 , 0.0339903200688579 , 0.2021249233072104 , 0.5277695132478635 , 0.2021249233072104 , 0.0002317033311045 );
            ir->AddPentPoint( 168 , 0.2021249233072104 , 0.0339903200688579 , 0.5277695132478635 , 0.2021249233072104 , 0.0002317033311045 );
            ir->AddPentPoint( 169 , 0.0339903200688579 , 0.0339903200688579 , 0.5277695132478635 , 0.2021249233072104 , 0.0002317033311045 );
            ir->AddPentPoint( 170 , 0.0339903200688579 , 0.5277695132478635 , 0.2021249233072104 , 0.2021249233072104 , 0.0002317033311045 );
            ir->AddPentPoint( 171 , 0.5277695132478635 , 0.0339903200688579 , 0.2021249233072104 , 0.2021249233072104 , 0.0002317033311045 );
            ir->AddPentPoint( 172 , 0.0339903200688579 , 0.0339903200688579 , 0.2021249233072104 , 0.2021249233072104 , 0.0002317033311045 );
            ir->AddPentPoint( 173 , 0.2021249233072104 , 0.5277695132478635 , 0.0339903200688579 , 0.2021249233072104 , 0.0002317033311045 );
            ir->AddPentPoint( 174 , 0.0339903200688579 , 0.5277695132478635 , 0.0339903200688579 , 0.2021249233072104 , 0.0002317033311045 );
            ir->AddPentPoint( 175 , 0.5277695132478635 , 0.2021249233072104 , 0.0339903200688579 , 0.2021249233072104 , 0.0002317033311045 );
            ir->AddPentPoint( 176 , 0.0339903200688579 , 0.2021249233072104 , 0.0339903200688579 , 0.2021249233072104 , 0.0002317033311045 );
            ir->AddPentPoint( 177 , 0.5277695132478635 , 0.0339903200688579 , 0.0339903200688579 , 0.2021249233072104 , 0.0002317033311045 );
            ir->AddPentPoint( 178 , 0.2021249233072104 , 0.0339903200688579 , 0.0339903200688579 , 0.2021249233072104 , 0.0002317033311045 );
            ir->AddPentPoint( 179 , 0.2021249233072104 , 0.2021249233072104 , 0.5277695132478635 , 0.0339903200688579 , 0.0002317033311045 );
            ir->AddPentPoint( 180 , 0.0339903200688579 , 0.2021249233072104 , 0.5277695132478635 , 0.0339903200688579 , 0.0002317033311045 );
            ir->AddPentPoint( 181 , 0.2021249233072104 , 0.0339903200688579 , 0.5277695132478635 , 0.0339903200688579 , 0.0002317033311045 );
            ir->AddPentPoint( 182 , 0.2021249233072104 , 0.5277695132478635 , 0.2021249233072104 , 0.0339903200688579 , 0.0002317033311045 );
            ir->AddPentPoint( 183 , 0.0339903200688579 , 0.5277695132478635 , 0.2021249233072104 , 0.0339903200688579 , 0.0002317033311045 );
            ir->AddPentPoint( 184 , 0.5277695132478635 , 0.2021249233072104 , 0.2021249233072104 , 0.0339903200688579 , 0.0002317033311045 );
            ir->AddPentPoint( 185 , 0.0339903200688579 , 0.2021249233072104 , 0.2021249233072104 , 0.0339903200688579 , 0.0002317033311045 );
            ir->AddPentPoint( 186 , 0.5277695132478635 , 0.0339903200688579 , 0.2021249233072104 , 0.0339903200688579 , 0.0002317033311045 );
            ir->AddPentPoint( 187 , 0.2021249233072104 , 0.0339903200688579 , 0.2021249233072104 , 0.0339903200688579 , 0.0002317033311045 );
            ir->AddPentPoint( 188 , 0.2021249233072104 , 0.5277695132478635 , 0.0339903200688579 , 0.0339903200688579 , 0.0002317033311045 );
            ir->AddPentPoint( 189 , 0.5277695132478635 , 0.2021249233072104 , 0.0339903200688579 , 0.0339903200688579 , 0.0002317033311045 );
            ir->AddPentPoint( 190 , 0.2021249233072104 , 0.2021249233072104 , 0.0339903200688579 , 0.0339903200688579 , 0.0002317033311045 );
            ir->AddPentPoint( 191 , 0.1089486003638299 , 0.3867810636690229 , 0.3867810636690229 , 0.0085406719342944 , 0.0001529932409945 );
            ir->AddPentPoint( 192 , 0.3867810636690229 , 0.1089486003638299 , 0.3867810636690229 , 0.0085406719342944 , 0.0001529932409945 );
            ir->AddPentPoint( 193 , 0.1089486003638299 , 0.1089486003638299 , 0.3867810636690229 , 0.0085406719342944 , 0.0001529932409945 );
            ir->AddPentPoint( 194 , 0.3867810636690229 , 0.3867810636690229 , 0.1089486003638299 , 0.0085406719342944 , 0.0001529932409945 );
            ir->AddPentPoint( 195 , 0.1089486003638299 , 0.3867810636690229 , 0.1089486003638299 , 0.0085406719342944 , 0.0001529932409945 );
            ir->AddPentPoint( 196 , 0.3867810636690229 , 0.1089486003638299 , 0.1089486003638299 , 0.0085406719342944 , 0.0001529932409945 );
            ir->AddPentPoint( 197 , 0.1089486003638299 , 0.3867810636690229 , 0.0085406719342944 , 0.3867810636690229 , 0.0001529932409945 );
            ir->AddPentPoint( 198 , 0.3867810636690229 , 0.1089486003638299 , 0.0085406719342944 , 0.3867810636690229 , 0.0001529932409945 );
            ir->AddPentPoint( 199 , 0.1089486003638299 , 0.1089486003638299 , 0.0085406719342944 , 0.3867810636690229 , 0.0001529932409945 );
            ir->AddPentPoint( 200 , 0.1089486003638299 , 0.0085406719342944 , 0.3867810636690229 , 0.3867810636690229 , 0.0001529932409945 );
            ir->AddPentPoint( 201 , 0.0085406719342944 , 0.1089486003638299 , 0.3867810636690229 , 0.3867810636690229 , 0.0001529932409945 );
            ir->AddPentPoint( 202 , 0.1089486003638299 , 0.1089486003638299 , 0.3867810636690229 , 0.3867810636690229 , 0.0001529932409945 );
            ir->AddPentPoint( 203 , 0.3867810636690229 , 0.0085406719342944 , 0.1089486003638299 , 0.3867810636690229 , 0.0001529932409945 );
            ir->AddPentPoint( 204 , 0.1089486003638299 , 0.0085406719342944 , 0.1089486003638299 , 0.3867810636690229 , 0.0001529932409945 );
            ir->AddPentPoint( 205 , 0.0085406719342944 , 0.3867810636690229 , 0.1089486003638299 , 0.3867810636690229 , 0.0001529932409945 );
            ir->AddPentPoint( 206 , 0.1089486003638299 , 0.3867810636690229 , 0.1089486003638299 , 0.3867810636690229 , 0.0001529932409945 );
            ir->AddPentPoint( 207 , 0.0085406719342944 , 0.1089486003638299 , 0.1089486003638299 , 0.3867810636690229 , 0.0001529932409945 );
            ir->AddPentPoint( 208 , 0.3867810636690229 , 0.1089486003638299 , 0.1089486003638299 , 0.3867810636690229 , 0.0001529932409945 );
            ir->AddPentPoint( 209 , 0.3867810636690229 , 0.3867810636690229 , 0.0085406719342944 , 0.1089486003638299 , 0.0001529932409945 );
            ir->AddPentPoint( 210 , 0.1089486003638299 , 0.3867810636690229 , 0.0085406719342944 , 0.1089486003638299 , 0.0001529932409945 );
            ir->AddPentPoint( 211 , 0.3867810636690229 , 0.1089486003638299 , 0.0085406719342944 , 0.1089486003638299 , 0.0001529932409945 );
            ir->AddPentPoint( 212 , 0.3867810636690229 , 0.0085406719342944 , 0.3867810636690229 , 0.1089486003638299 , 0.0001529932409945 );
            ir->AddPentPoint( 213 , 0.1089486003638299 , 0.0085406719342944 , 0.3867810636690229 , 0.1089486003638299 , 0.0001529932409945 );
            ir->AddPentPoint( 214 , 0.0085406719342944 , 0.3867810636690229 , 0.3867810636690229 , 0.1089486003638299 , 0.0001529932409945 );
            ir->AddPentPoint( 215 , 0.1089486003638299 , 0.3867810636690229 , 0.3867810636690229 , 0.1089486003638299 , 0.0001529932409945 );
            ir->AddPentPoint( 216 , 0.0085406719342944 , 0.1089486003638299 , 0.3867810636690229 , 0.1089486003638299 , 0.0001529932409945 );
            ir->AddPentPoint( 217 , 0.3867810636690229 , 0.1089486003638299 , 0.3867810636690229 , 0.1089486003638299 , 0.0001529932409945 );
            ir->AddPentPoint( 218 , 0.3867810636690229 , 0.0085406719342944 , 0.1089486003638299 , 0.1089486003638299 , 0.0001529932409945 );
            ir->AddPentPoint( 219 , 0.0085406719342944 , 0.3867810636690229 , 0.1089486003638299 , 0.1089486003638299 , 0.0001529932409945 );
            ir->AddPentPoint( 220 , 0.3867810636690229 , 0.3867810636690229 , 0.1089486003638299 , 0.1089486003638299 , 0.0001529932409945 );
            ir->AddPentPoint( 221 , 0.0031786137103847 , 0.0666738318226067 , 0.1528608418201177 , 0.7741080989365061 , 0.0000138796264509 );
            ir->AddPentPoint( 222 , 0.0666738318226067 , 0.0031786137103847 , 0.1528608418201177 , 0.7741080989365061 , 0.0000138796264509 );
            ir->AddPentPoint( 223 , 0.0031786137103847 , 0.0031786137103847 , 0.1528608418201177 , 0.7741080989365061 , 0.0000138796264509 );
            ir->AddPentPoint( 224 , 0.0031786137103847 , 0.1528608418201177 , 0.0666738318226067 , 0.7741080989365061 , 0.0000138796264509 );
            ir->AddPentPoint( 225 , 0.1528608418201177 , 0.0031786137103847 , 0.0666738318226067 , 0.7741080989365061 , 0.0000138796264509 );
            ir->AddPentPoint( 226 , 0.0031786137103847 , 0.0031786137103847 , 0.0666738318226067 , 0.7741080989365061 , 0.0000138796264509 );
            ir->AddPentPoint( 227 , 0.0666738318226067 , 0.1528608418201177 , 0.0031786137103847 , 0.7741080989365061 , 0.0000138796264509 );
            ir->AddPentPoint( 228 , 0.0031786137103847 , 0.1528608418201177 , 0.0031786137103847 , 0.7741080989365061 , 0.0000138796264509 );
            ir->AddPentPoint( 229 , 0.1528608418201177 , 0.0666738318226067 , 0.0031786137103847 , 0.7741080989365061 , 0.0000138796264509 );
            ir->AddPentPoint( 230 , 0.0031786137103847 , 0.0666738318226067 , 0.0031786137103847 , 0.7741080989365061 , 0.0000138796264509 );
            ir->AddPentPoint( 231 , 0.1528608418201177 , 0.0031786137103847 , 0.0031786137103847 , 0.7741080989365061 , 0.0000138796264509 );
            ir->AddPentPoint( 232 , 0.0666738318226067 , 0.0031786137103847 , 0.0031786137103847 , 0.7741080989365061 , 0.0000138796264509 );
            ir->AddPentPoint( 233 , 0.0031786137103847 , 0.0666738318226067 , 0.7741080989365061 , 0.1528608418201177 , 0.0000138796264509 );
            ir->AddPentPoint( 234 , 0.0666738318226067 , 0.0031786137103847 , 0.7741080989365061 , 0.1528608418201177 , 0.0000138796264509 );
            ir->AddPentPoint( 235 , 0.0031786137103847 , 0.0031786137103847 , 0.7741080989365061 , 0.1528608418201177 , 0.0000138796264509 );
            ir->AddPentPoint( 236 , 0.0031786137103847 , 0.7741080989365061 , 0.0666738318226067 , 0.1528608418201177 , 0.0000138796264509 );
            ir->AddPentPoint( 237 , 0.7741080989365061 , 0.0031786137103847 , 0.0666738318226067 , 0.1528608418201177 , 0.0000138796264509 );
            ir->AddPentPoint( 238 , 0.0031786137103847 , 0.0031786137103847 , 0.0666738318226067 , 0.1528608418201177 , 0.0000138796264509 );
            ir->AddPentPoint( 239 , 0.0666738318226067 , 0.7741080989365061 , 0.0031786137103847 , 0.1528608418201177 , 0.0000138796264509 );
            ir->AddPentPoint( 240 , 0.0031786137103847 , 0.7741080989365061 , 0.0031786137103847 , 0.1528608418201177 , 0.0000138796264509 );
            ir->AddPentPoint( 241 , 0.7741080989365061 , 0.0666738318226067 , 0.0031786137103847 , 0.1528608418201177 , 0.0000138796264509 );
            ir->AddPentPoint( 242 , 0.0031786137103847 , 0.0666738318226067 , 0.0031786137103847 , 0.1528608418201177 , 0.0000138796264509 );
            ir->AddPentPoint( 243 , 0.7741080989365061 , 0.0031786137103847 , 0.0031786137103847 , 0.1528608418201177 , 0.0000138796264509 );
            ir->AddPentPoint( 244 , 0.0666738318226067 , 0.0031786137103847 , 0.0031786137103847 , 0.1528608418201177 , 0.0000138796264509 );
            ir->AddPentPoint( 245 , 0.0031786137103847 , 0.1528608418201177 , 0.7741080989365061 , 0.0666738318226067 , 0.0000138796264509 );
            ir->AddPentPoint( 246 , 0.1528608418201177 , 0.0031786137103847 , 0.7741080989365061 , 0.0666738318226067 , 0.0000138796264509 );
            ir->AddPentPoint( 247 , 0.0031786137103847 , 0.0031786137103847 , 0.7741080989365061 , 0.0666738318226067 , 0.0000138796264509 );
            ir->AddPentPoint( 248 , 0.0031786137103847 , 0.7741080989365061 , 0.1528608418201177 , 0.0666738318226067 , 0.0000138796264509 );
            ir->AddPentPoint( 249 , 0.7741080989365061 , 0.0031786137103847 , 0.1528608418201177 , 0.0666738318226067 , 0.0000138796264509 );
            ir->AddPentPoint( 250 , 0.0031786137103847 , 0.0031786137103847 , 0.1528608418201177 , 0.0666738318226067 , 0.0000138796264509 );
            ir->AddPentPoint( 251 , 0.1528608418201177 , 0.7741080989365061 , 0.0031786137103847 , 0.0666738318226067 , 0.0000138796264509 );
            ir->AddPentPoint( 252 , 0.0031786137103847 , 0.7741080989365061 , 0.0031786137103847 , 0.0666738318226067 , 0.0000138796264509 );
            ir->AddPentPoint( 253 , 0.7741080989365061 , 0.1528608418201177 , 0.0031786137103847 , 0.0666738318226067 , 0.0000138796264509 );
            ir->AddPentPoint( 254 , 0.0031786137103847 , 0.1528608418201177 , 0.0031786137103847 , 0.0666738318226067 , 0.0000138796264509 );
            ir->AddPentPoint( 255 , 0.7741080989365061 , 0.0031786137103847 , 0.0031786137103847 , 0.0666738318226067 , 0.0000138796264509 );
            ir->AddPentPoint( 256 , 0.1528608418201177 , 0.0031786137103847 , 0.0031786137103847 , 0.0666738318226067 , 0.0000138796264509 );
            ir->AddPentPoint( 257 , 0.0666738318226067 , 0.1528608418201177 , 0.7741080989365061 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 258 , 0.0031786137103847 , 0.1528608418201177 , 0.7741080989365061 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 259 , 0.1528608418201177 , 0.0666738318226067 , 0.7741080989365061 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 260 , 0.0031786137103847 , 0.0666738318226067 , 0.7741080989365061 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 261 , 0.1528608418201177 , 0.0031786137103847 , 0.7741080989365061 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 262 , 0.0666738318226067 , 0.0031786137103847 , 0.7741080989365061 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 263 , 0.0666738318226067 , 0.7741080989365061 , 0.1528608418201177 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 264 , 0.0031786137103847 , 0.7741080989365061 , 0.1528608418201177 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 265 , 0.7741080989365061 , 0.0666738318226067 , 0.1528608418201177 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 266 , 0.0031786137103847 , 0.0666738318226067 , 0.1528608418201177 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 267 , 0.7741080989365061 , 0.0031786137103847 , 0.1528608418201177 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 268 , 0.0666738318226067 , 0.0031786137103847 , 0.1528608418201177 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 269 , 0.1528608418201177 , 0.7741080989365061 , 0.0666738318226067 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 270 , 0.0031786137103847 , 0.7741080989365061 , 0.0666738318226067 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 271 , 0.7741080989365061 , 0.1528608418201177 , 0.0666738318226067 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 272 , 0.0031786137103847 , 0.1528608418201177 , 0.0666738318226067 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 273 , 0.7741080989365061 , 0.0031786137103847 , 0.0666738318226067 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 274 , 0.1528608418201177 , 0.0031786137103847 , 0.0666738318226067 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 275 , 0.1528608418201177 , 0.7741080989365061 , 0.0031786137103847 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 276 , 0.0666738318226067 , 0.7741080989365061 , 0.0031786137103847 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 277 , 0.7741080989365061 , 0.1528608418201177 , 0.0031786137103847 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 278 , 0.0666738318226067 , 0.1528608418201177 , 0.0031786137103847 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 279 , 0.7741080989365061 , 0.0666738318226067 , 0.0031786137103847 , 0.0031786137103847 , 0.0000138796264509 );
            ir->AddPentPoint( 280 , 0.1528608418201177 , 0.0666738318226067 , 0.0031786137103847 , 0.0031786137103847 , 0.0000138796264509 );
            
            return ir;
            
        case 12:  // 445 points - degree 12 --
            PentatopeIntRules[12] = ir = new IntegrationRule(445);
            ir->AddPentPoint( 0 , 0.0187405455280847 , 0.0187405455280847 , 0.0187405455280847 , 0.9250378178876613 , 0.0000068131746956 );
            ir->AddPentPoint( 1 , 0.0187405455280847 , 0.0187405455280847 , 0.9250378178876613 , 0.0187405455280847 , 0.0000068131746956 );
            ir->AddPentPoint( 2 , 0.0187405455280847 , 0.9250378178876613 , 0.0187405455280847 , 0.0187405455280847 , 0.0000068131746956 );
            ir->AddPentPoint( 3 , 0.9250378178876613 , 0.0187405455280847 , 0.0187405455280847 , 0.0187405455280847 , 0.0000068131746956 );
            ir->AddPentPoint( 4 , 0.0187405455280847 , 0.0187405455280847 , 0.0187405455280847 , 0.0187405455280847 , 0.0000068131746956 );
            ir->AddPentPoint( 5 , 0.1616906732413428 , 0.1616906732413428 , 0.1616906732413428 , 0.3532373070346286 , 0.0005714071565234 );
            ir->AddPentPoint( 6 , 0.1616906732413428 , 0.1616906732413428 , 0.3532373070346286 , 0.1616906732413428 , 0.0005714071565234 );
            ir->AddPentPoint( 7 , 0.1616906732413428 , 0.3532373070346286 , 0.1616906732413428 , 0.1616906732413428 , 0.0005714071565234 );
            ir->AddPentPoint( 8 , 0.3532373070346286 , 0.1616906732413428 , 0.1616906732413428 , 0.1616906732413428 , 0.0005714071565234 );
            ir->AddPentPoint( 9 , 0.1616906732413428 , 0.1616906732413428 , 0.1616906732413428 , 0.1616906732413428 , 0.0005714071565234 );
            ir->AddPentPoint( 10 , 0.2416166418750493 , 0.2416166418750493 , 0.2416166418750493 , 0.0335334324998028 , 0.0003386291960807 );
            ir->AddPentPoint( 11 , 0.2416166418750493 , 0.2416166418750493 , 0.0335334324998028 , 0.2416166418750493 , 0.0003386291960807 );
            ir->AddPentPoint( 12 , 0.2416166418750493 , 0.0335334324998028 , 0.2416166418750493 , 0.2416166418750493 , 0.0003386291960807 );
            ir->AddPentPoint( 13 , 0.0335334324998028 , 0.2416166418750493 , 0.2416166418750493 , 0.2416166418750493 , 0.0003386291960807 );
            ir->AddPentPoint( 14 , 0.2416166418750493 , 0.2416166418750493 , 0.2416166418750493 , 0.2416166418750493 , 0.0003386291960807 );
            ir->AddPentPoint( 15 , 0.0238708641272240 , 0.0238708641272240 , 0.4641937038091640 , 0.4641937038091640 , 0.0000446911422738 );
            ir->AddPentPoint( 16 , 0.0238708641272240 , 0.4641937038091640 , 0.0238708641272240 , 0.4641937038091640 , 0.0000446911422738 );
            ir->AddPentPoint( 17 , 0.4641937038091640 , 0.0238708641272240 , 0.0238708641272240 , 0.4641937038091640 , 0.0000446911422738 );
            ir->AddPentPoint( 18 , 0.0238708641272240 , 0.0238708641272240 , 0.0238708641272240 , 0.4641937038091640 , 0.0000446911422738 );
            ir->AddPentPoint( 19 , 0.0238708641272240 , 0.4641937038091640 , 0.4641937038091640 , 0.0238708641272240 , 0.0000446911422738 );
            ir->AddPentPoint( 20 , 0.4641937038091640 , 0.0238708641272240 , 0.4641937038091640 , 0.0238708641272240 , 0.0000446911422738 );
            ir->AddPentPoint( 21 , 0.0238708641272240 , 0.0238708641272240 , 0.4641937038091640 , 0.0238708641272240 , 0.0000446911422738 );
            ir->AddPentPoint( 22 , 0.4641937038091640 , 0.4641937038091640 , 0.0238708641272240 , 0.0238708641272240 , 0.0000446911422738 );
            ir->AddPentPoint( 23 , 0.0238708641272240 , 0.4641937038091640 , 0.0238708641272240 , 0.0238708641272240 , 0.0000446911422738 );
            ir->AddPentPoint( 24 , 0.4641937038091640 , 0.0238708641272240 , 0.0238708641272240 , 0.0238708641272240 , 0.0000446911422738 );
            ir->AddPentPoint( 25 , 0.0869616529443132 , 0.0869616529443132 , 0.3695575205835303 , 0.3695575205835303 , 0.0002872047109476 );
            ir->AddPentPoint( 26 , 0.0869616529443132 , 0.3695575205835303 , 0.0869616529443132 , 0.3695575205835303 , 0.0002872047109476 );
            ir->AddPentPoint( 27 , 0.3695575205835303 , 0.0869616529443132 , 0.0869616529443132 , 0.3695575205835303 , 0.0002872047109476 );
            ir->AddPentPoint( 28 , 0.0869616529443132 , 0.0869616529443132 , 0.0869616529443132 , 0.3695575205835303 , 0.0002872047109476 );
            ir->AddPentPoint( 29 , 0.0869616529443132 , 0.3695575205835303 , 0.3695575205835303 , 0.0869616529443132 , 0.0002872047109476 );
            ir->AddPentPoint( 30 , 0.3695575205835303 , 0.0869616529443132 , 0.3695575205835303 , 0.0869616529443132 , 0.0002872047109476 );
            ir->AddPentPoint( 31 , 0.0869616529443132 , 0.0869616529443132 , 0.3695575205835303 , 0.0869616529443132 , 0.0002872047109476 );
            ir->AddPentPoint( 32 , 0.3695575205835303 , 0.3695575205835303 , 0.0869616529443132 , 0.0869616529443132 , 0.0002872047109476 );
            ir->AddPentPoint( 33 , 0.0869616529443132 , 0.3695575205835303 , 0.0869616529443132 , 0.0869616529443132 , 0.0002872047109476 );
            ir->AddPentPoint( 34 , 0.3695575205835303 , 0.0869616529443132 , 0.0869616529443132 , 0.0869616529443132 , 0.0002872047109476 );
            ir->AddPentPoint( 35 , 0.2754833302419382 , 0.2754833302419382 , 0.0867750046370928 , 0.0867750046370928 , 0.0004347276084874 );
            ir->AddPentPoint( 36 , 0.2754833302419382 , 0.0867750046370928 , 0.2754833302419382 , 0.0867750046370928 , 0.0004347276084874 );
            ir->AddPentPoint( 37 , 0.0867750046370928 , 0.2754833302419382 , 0.2754833302419382 , 0.0867750046370928 , 0.0004347276084874 );
            ir->AddPentPoint( 38 , 0.2754833302419382 , 0.2754833302419382 , 0.2754833302419382 , 0.0867750046370928 , 0.0004347276084874 );
            ir->AddPentPoint( 39 , 0.2754833302419382 , 0.0867750046370928 , 0.0867750046370928 , 0.2754833302419382 , 0.0004347276084874 );
            ir->AddPentPoint( 40 , 0.0867750046370928 , 0.2754833302419382 , 0.0867750046370928 , 0.2754833302419382 , 0.0004347276084874 );
            ir->AddPentPoint( 41 , 0.2754833302419382 , 0.2754833302419382 , 0.0867750046370928 , 0.2754833302419382 , 0.0004347276084874 );
            ir->AddPentPoint( 42 , 0.0867750046370928 , 0.0867750046370928 , 0.2754833302419382 , 0.2754833302419382 , 0.0004347276084874 );
            ir->AddPentPoint( 43 , 0.2754833302419382 , 0.0867750046370928 , 0.2754833302419382 , 0.2754833302419382 , 0.0004347276084874 );
            ir->AddPentPoint( 44 , 0.0867750046370928 , 0.2754833302419382 , 0.2754833302419382 , 0.2754833302419382 , 0.0004347276084874 );
            ir->AddPentPoint( 45 , 0.0000000018557990 , 0.0000000018557990 , 0.2789044452126671 , 0.7210955492199358 , 0.0000027744431139 );
            ir->AddPentPoint( 46 , 0.0000000018557990 , 0.2789044452126671 , 0.0000000018557990 , 0.7210955492199358 , 0.0000027744431139 );
            ir->AddPentPoint( 47 , 0.2789044452126671 , 0.0000000018557990 , 0.0000000018557990 , 0.7210955492199358 , 0.0000027744431139 );
            ir->AddPentPoint( 48 , 0.0000000018557990 , 0.0000000018557990 , 0.0000000018557990 , 0.7210955492199358 , 0.0000027744431139 );
            ir->AddPentPoint( 49 , 0.0000000018557990 , 0.0000000018557990 , 0.7210955492199358 , 0.2789044452126671 , 0.0000027744431139 );
            ir->AddPentPoint( 50 , 0.0000000018557990 , 0.7210955492199358 , 0.0000000018557990 , 0.2789044452126671 , 0.0000027744431139 );
            ir->AddPentPoint( 51 , 0.7210955492199358 , 0.0000000018557990 , 0.0000000018557990 , 0.2789044452126671 , 0.0000027744431139 );
            ir->AddPentPoint( 52 , 0.0000000018557990 , 0.0000000018557990 , 0.0000000018557990 , 0.2789044452126671 , 0.0000027744431139 );
            ir->AddPentPoint( 53 , 0.0000000018557990 , 0.2789044452126671 , 0.7210955492199358 , 0.0000000018557990 , 0.0000027744431139 );
            ir->AddPentPoint( 54 , 0.2789044452126671 , 0.0000000018557990 , 0.7210955492199358 , 0.0000000018557990 , 0.0000027744431139 );
            ir->AddPentPoint( 55 , 0.0000000018557990 , 0.0000000018557990 , 0.7210955492199358 , 0.0000000018557990 , 0.0000027744431139 );
            ir->AddPentPoint( 56 , 0.0000000018557990 , 0.7210955492199358 , 0.2789044452126671 , 0.0000000018557990 , 0.0000027744431139 );
            ir->AddPentPoint( 57 , 0.7210955492199358 , 0.0000000018557990 , 0.2789044452126671 , 0.0000000018557990 , 0.0000027744431139 );
            ir->AddPentPoint( 58 , 0.0000000018557990 , 0.0000000018557990 , 0.2789044452126671 , 0.0000000018557990 , 0.0000027744431139 );
            ir->AddPentPoint( 59 , 0.2789044452126671 , 0.7210955492199358 , 0.0000000018557990 , 0.0000000018557990 , 0.0000027744431139 );
            ir->AddPentPoint( 60 , 0.0000000018557990 , 0.7210955492199358 , 0.0000000018557990 , 0.0000000018557990 , 0.0000027744431139 );
            ir->AddPentPoint( 61 , 0.7210955492199358 , 0.2789044452126671 , 0.0000000018557990 , 0.0000000018557990 , 0.0000027744431139 );
            ir->AddPentPoint( 62 , 0.0000000018557990 , 0.2789044452126671 , 0.0000000018557990 , 0.0000000018557990 , 0.0000027744431139 );
            ir->AddPentPoint( 63 , 0.7210955492199358 , 0.0000000018557990 , 0.0000000018557990 , 0.0000000018557990 , 0.0000027744431139 );
            ir->AddPentPoint( 64 , 0.2789044452126671 , 0.0000000018557990 , 0.0000000018557990 , 0.0000000018557990 , 0.0000027744431139 );
            ir->AddPentPoint( 65 , 0.0260584033343716 , 0.0260584033343716 , 0.1110382975784194 , 0.8107864924184656 , 0.0000350244076538 );
            ir->AddPentPoint( 66 , 0.0260584033343716 , 0.1110382975784194 , 0.0260584033343716 , 0.8107864924184656 , 0.0000350244076538 );
            ir->AddPentPoint( 67 , 0.1110382975784194 , 0.0260584033343716 , 0.0260584033343716 , 0.8107864924184656 , 0.0000350244076538 );
            ir->AddPentPoint( 68 , 0.0260584033343716 , 0.0260584033343716 , 0.0260584033343716 , 0.8107864924184656 , 0.0000350244076538 );
            ir->AddPentPoint( 69 , 0.0260584033343716 , 0.0260584033343716 , 0.8107864924184656 , 0.1110382975784194 , 0.0000350244076538 );
            ir->AddPentPoint( 70 , 0.0260584033343716 , 0.8107864924184656 , 0.0260584033343716 , 0.1110382975784194 , 0.0000350244076538 );
            ir->AddPentPoint( 71 , 0.8107864924184656 , 0.0260584033343716 , 0.0260584033343716 , 0.1110382975784194 , 0.0000350244076538 );
            ir->AddPentPoint( 72 , 0.0260584033343716 , 0.0260584033343716 , 0.0260584033343716 , 0.1110382975784194 , 0.0000350244076538 );
            ir->AddPentPoint( 73 , 0.0260584033343716 , 0.1110382975784194 , 0.8107864924184656 , 0.0260584033343716 , 0.0000350244076538 );
            ir->AddPentPoint( 74 , 0.1110382975784194 , 0.0260584033343716 , 0.8107864924184656 , 0.0260584033343716 , 0.0000350244076538 );
            ir->AddPentPoint( 75 , 0.0260584033343716 , 0.0260584033343716 , 0.8107864924184656 , 0.0260584033343716 , 0.0000350244076538 );
            ir->AddPentPoint( 76 , 0.0260584033343716 , 0.8107864924184656 , 0.1110382975784194 , 0.0260584033343716 , 0.0000350244076538 );
            ir->AddPentPoint( 77 , 0.8107864924184656 , 0.0260584033343716 , 0.1110382975784194 , 0.0260584033343716 , 0.0000350244076538 );
            ir->AddPentPoint( 78 , 0.0260584033343716 , 0.0260584033343716 , 0.1110382975784194 , 0.0260584033343716 , 0.0000350244076538 );
            ir->AddPentPoint( 79 , 0.1110382975784194 , 0.8107864924184656 , 0.0260584033343716 , 0.0260584033343716 , 0.0000350244076538 );
            ir->AddPentPoint( 80 , 0.0260584033343716 , 0.8107864924184656 , 0.0260584033343716 , 0.0260584033343716 , 0.0000350244076538 );
            ir->AddPentPoint( 81 , 0.8107864924184656 , 0.1110382975784194 , 0.0260584033343716 , 0.0260584033343716 , 0.0000350244076538 );
            ir->AddPentPoint( 82 , 0.0260584033343716 , 0.1110382975784194 , 0.0260584033343716 , 0.0260584033343716 , 0.0000350244076538 );
            ir->AddPentPoint( 83 , 0.8107864924184656 , 0.0260584033343716 , 0.0260584033343716 , 0.0260584033343716 , 0.0000350244076538 );
            ir->AddPentPoint( 84 , 0.1110382975784194 , 0.0260584033343716 , 0.0260584033343716 , 0.0260584033343716 , 0.0000350244076538 );
            ir->AddPentPoint( 85 , 0.0580662618708326 , 0.0580662618708326 , 0.1734657070267037 , 0.6523355073607985 , 0.0000984644006275 );
            ir->AddPentPoint( 86 , 0.0580662618708326 , 0.1734657070267037 , 0.0580662618708326 , 0.6523355073607985 , 0.0000984644006275 );
            ir->AddPentPoint( 87 , 0.1734657070267037 , 0.0580662618708326 , 0.0580662618708326 , 0.6523355073607985 , 0.0000984644006275 );
            ir->AddPentPoint( 88 , 0.0580662618708326 , 0.0580662618708326 , 0.0580662618708326 , 0.6523355073607985 , 0.0000984644006275 );
            ir->AddPentPoint( 89 , 0.0580662618708326 , 0.0580662618708326 , 0.6523355073607985 , 0.1734657070267037 , 0.0000984644006275 );
            ir->AddPentPoint( 90 , 0.0580662618708326 , 0.6523355073607985 , 0.0580662618708326 , 0.1734657070267037 , 0.0000984644006275 );
            ir->AddPentPoint( 91 , 0.6523355073607985 , 0.0580662618708326 , 0.0580662618708326 , 0.1734657070267037 , 0.0000984644006275 );
            ir->AddPentPoint( 92 , 0.0580662618708326 , 0.0580662618708326 , 0.0580662618708326 , 0.1734657070267037 , 0.0000984644006275 );
            ir->AddPentPoint( 93 , 0.0580662618708326 , 0.1734657070267037 , 0.6523355073607985 , 0.0580662618708326 , 0.0000984644006275 );
            ir->AddPentPoint( 94 , 0.1734657070267037 , 0.0580662618708326 , 0.6523355073607985 , 0.0580662618708326 , 0.0000984644006275 );
            ir->AddPentPoint( 95 , 0.0580662618708326 , 0.0580662618708326 , 0.6523355073607985 , 0.0580662618708326 , 0.0000984644006275 );
            ir->AddPentPoint( 96 , 0.0580662618708326 , 0.6523355073607985 , 0.1734657070267037 , 0.0580662618708326 , 0.0000984644006275 );
            ir->AddPentPoint( 97 , 0.6523355073607985 , 0.0580662618708326 , 0.1734657070267037 , 0.0580662618708326 , 0.0000984644006275 );
            ir->AddPentPoint( 98 , 0.0580662618708326 , 0.0580662618708326 , 0.1734657070267037 , 0.0580662618708326 , 0.0000984644006275 );
            ir->AddPentPoint( 99 , 0.1734657070267037 , 0.6523355073607985 , 0.0580662618708326 , 0.0580662618708326 , 0.0000984644006275 );
            ir->AddPentPoint( 100 , 0.0580662618708326 , 0.6523355073607985 , 0.0580662618708326 , 0.0580662618708326 , 0.0000984644006275 );
            ir->AddPentPoint( 101 , 0.6523355073607985 , 0.1734657070267037 , 0.0580662618708326 , 0.0580662618708326 , 0.0000984644006275 );
            ir->AddPentPoint( 102 , 0.0580662618708326 , 0.1734657070267037 , 0.0580662618708326 , 0.0580662618708326 , 0.0000984644006275 );
            ir->AddPentPoint( 103 , 0.6523355073607985 , 0.0580662618708326 , 0.0580662618708326 , 0.0580662618708326 , 0.0000984644006275 );
            ir->AddPentPoint( 104 , 0.1734657070267037 , 0.0580662618708326 , 0.0580662618708326 , 0.0580662618708326 , 0.0000984644006275 );
            ir->AddPentPoint( 105 , 0.1043482687149817 , 0.1043482687149817 , 0.0093049390716579 , 0.6776502547833971 , 0.0000686204897549 );
            ir->AddPentPoint( 106 , 0.1043482687149817 , 0.0093049390716579 , 0.1043482687149817 , 0.6776502547833971 , 0.0000686204897549 );
            ir->AddPentPoint( 107 , 0.0093049390716579 , 0.1043482687149817 , 0.1043482687149817 , 0.6776502547833971 , 0.0000686204897549 );
            ir->AddPentPoint( 108 , 0.1043482687149817 , 0.1043482687149817 , 0.1043482687149817 , 0.6776502547833971 , 0.0000686204897549 );
            ir->AddPentPoint( 109 , 0.1043482687149817 , 0.1043482687149817 , 0.6776502547833971 , 0.0093049390716579 , 0.0000686204897549 );
            ir->AddPentPoint( 110 , 0.1043482687149817 , 0.6776502547833971 , 0.1043482687149817 , 0.0093049390716579 , 0.0000686204897549 );
            ir->AddPentPoint( 111 , 0.6776502547833971 , 0.1043482687149817 , 0.1043482687149817 , 0.0093049390716579 , 0.0000686204897549 );
            ir->AddPentPoint( 112 , 0.1043482687149817 , 0.1043482687149817 , 0.1043482687149817 , 0.0093049390716579 , 0.0000686204897549 );
            ir->AddPentPoint( 113 , 0.1043482687149817 , 0.0093049390716579 , 0.6776502547833971 , 0.1043482687149817 , 0.0000686204897549 );
            ir->AddPentPoint( 114 , 0.0093049390716579 , 0.1043482687149817 , 0.6776502547833971 , 0.1043482687149817 , 0.0000686204897549 );
            ir->AddPentPoint( 115 , 0.1043482687149817 , 0.1043482687149817 , 0.6776502547833971 , 0.1043482687149817 , 0.0000686204897549 );
            ir->AddPentPoint( 116 , 0.1043482687149817 , 0.6776502547833971 , 0.0093049390716579 , 0.1043482687149817 , 0.0000686204897549 );
            ir->AddPentPoint( 117 , 0.6776502547833971 , 0.1043482687149817 , 0.0093049390716579 , 0.1043482687149817 , 0.0000686204897549 );
            ir->AddPentPoint( 118 , 0.1043482687149817 , 0.1043482687149817 , 0.0093049390716579 , 0.1043482687149817 , 0.0000686204897549 );
            ir->AddPentPoint( 119 , 0.0093049390716579 , 0.6776502547833971 , 0.1043482687149817 , 0.1043482687149817 , 0.0000686204897549 );
            ir->AddPentPoint( 120 , 0.1043482687149817 , 0.6776502547833971 , 0.1043482687149817 , 0.1043482687149817 , 0.0000686204897549 );
            ir->AddPentPoint( 121 , 0.6776502547833971 , 0.0093049390716579 , 0.1043482687149817 , 0.1043482687149817 , 0.0000686204897549 );
            ir->AddPentPoint( 122 , 0.1043482687149817 , 0.0093049390716579 , 0.1043482687149817 , 0.1043482687149817 , 0.0000686204897549 );
            ir->AddPentPoint( 123 , 0.6776502547833971 , 0.1043482687149817 , 0.1043482687149817 , 0.1043482687149817 , 0.0000686204897549 );
            ir->AddPentPoint( 124 , 0.0093049390716579 , 0.1043482687149817 , 0.1043482687149817 , 0.1043482687149817 , 0.0000686204897549 );
            ir->AddPentPoint( 125 , 0.3138484211983653 , 0.3138484211983653 , 0.0058960507492770 , 0.0525586856556272 , 0.0000924807728581 );
            ir->AddPentPoint( 126 , 0.3138484211983653 , 0.0058960507492770 , 0.3138484211983653 , 0.0525586856556272 , 0.0000924807728581 );
            ir->AddPentPoint( 127 , 0.0058960507492770 , 0.3138484211983653 , 0.3138484211983653 , 0.0525586856556272 , 0.0000924807728581 );
            ir->AddPentPoint( 128 , 0.3138484211983653 , 0.3138484211983653 , 0.3138484211983653 , 0.0525586856556272 , 0.0000924807728581 );
            ir->AddPentPoint( 129 , 0.3138484211983653 , 0.3138484211983653 , 0.0525586856556272 , 0.0058960507492770 , 0.0000924807728581 );
            ir->AddPentPoint( 130 , 0.3138484211983653 , 0.0525586856556272 , 0.3138484211983653 , 0.0058960507492770 , 0.0000924807728581 );
            ir->AddPentPoint( 131 , 0.0525586856556272 , 0.3138484211983653 , 0.3138484211983653 , 0.0058960507492770 , 0.0000924807728581 );
            ir->AddPentPoint( 132 , 0.3138484211983653 , 0.3138484211983653 , 0.3138484211983653 , 0.0058960507492770 , 0.0000924807728581 );
            ir->AddPentPoint( 133 , 0.3138484211983653 , 0.0058960507492770 , 0.0525586856556272 , 0.3138484211983653 , 0.0000924807728581 );
            ir->AddPentPoint( 134 , 0.0058960507492770 , 0.3138484211983653 , 0.0525586856556272 , 0.3138484211983653 , 0.0000924807728581 );
            ir->AddPentPoint( 135 , 0.3138484211983653 , 0.3138484211983653 , 0.0525586856556272 , 0.3138484211983653 , 0.0000924807728581 );
            ir->AddPentPoint( 136 , 0.3138484211983653 , 0.0525586856556272 , 0.0058960507492770 , 0.3138484211983653 , 0.0000924807728581 );
            ir->AddPentPoint( 137 , 0.0525586856556272 , 0.3138484211983653 , 0.0058960507492770 , 0.3138484211983653 , 0.0000924807728581 );
            ir->AddPentPoint( 138 , 0.3138484211983653 , 0.3138484211983653 , 0.0058960507492770 , 0.3138484211983653 , 0.0000924807728581 );
            ir->AddPentPoint( 139 , 0.0058960507492770 , 0.0525586856556272 , 0.3138484211983653 , 0.3138484211983653 , 0.0000924807728581 );
            ir->AddPentPoint( 140 , 0.3138484211983653 , 0.0525586856556272 , 0.3138484211983653 , 0.3138484211983653 , 0.0000924807728581 );
            ir->AddPentPoint( 141 , 0.0525586856556272 , 0.0058960507492770 , 0.3138484211983653 , 0.3138484211983653 , 0.0000924807728581 );
            ir->AddPentPoint( 142 , 0.3138484211983653 , 0.0058960507492770 , 0.3138484211983653 , 0.3138484211983653 , 0.0000924807728581 );
            ir->AddPentPoint( 143 , 0.0525586856556272 , 0.3138484211983653 , 0.3138484211983653 , 0.3138484211983653 , 0.0000924807728581 );
            ir->AddPentPoint( 144 , 0.0058960507492770 , 0.3138484211983653 , 0.3138484211983653 , 0.3138484211983653 , 0.0000924807728581 );
            ir->AddPentPoint( 145 , 0.0002635287125677 , 0.1328788841147149 , 0.1328788841147149 , 0.7337151743454349 , 0.0000120687092035 );
            ir->AddPentPoint( 146 , 0.1328788841147149 , 0.0002635287125677 , 0.1328788841147149 , 0.7337151743454349 , 0.0000120687092035 );
            ir->AddPentPoint( 147 , 0.0002635287125677 , 0.0002635287125677 , 0.1328788841147149 , 0.7337151743454349 , 0.0000120687092035 );
            ir->AddPentPoint( 148 , 0.1328788841147149 , 0.1328788841147149 , 0.0002635287125677 , 0.7337151743454349 , 0.0000120687092035 );
            ir->AddPentPoint( 149 , 0.0002635287125677 , 0.1328788841147149 , 0.0002635287125677 , 0.7337151743454349 , 0.0000120687092035 );
            ir->AddPentPoint( 150 , 0.1328788841147149 , 0.0002635287125677 , 0.0002635287125677 , 0.7337151743454349 , 0.0000120687092035 );
            ir->AddPentPoint( 151 , 0.0002635287125677 , 0.1328788841147149 , 0.7337151743454349 , 0.1328788841147149 , 0.0000120687092035 );
            ir->AddPentPoint( 152 , 0.1328788841147149 , 0.0002635287125677 , 0.7337151743454349 , 0.1328788841147149 , 0.0000120687092035 );
            ir->AddPentPoint( 153 , 0.0002635287125677 , 0.0002635287125677 , 0.7337151743454349 , 0.1328788841147149 , 0.0000120687092035 );
            ir->AddPentPoint( 154 , 0.0002635287125677 , 0.7337151743454349 , 0.1328788841147149 , 0.1328788841147149 , 0.0000120687092035 );
            ir->AddPentPoint( 155 , 0.7337151743454349 , 0.0002635287125677 , 0.1328788841147149 , 0.1328788841147149 , 0.0000120687092035 );
            ir->AddPentPoint( 156 , 0.0002635287125677 , 0.0002635287125677 , 0.1328788841147149 , 0.1328788841147149 , 0.0000120687092035 );
            ir->AddPentPoint( 157 , 0.1328788841147149 , 0.7337151743454349 , 0.0002635287125677 , 0.1328788841147149 , 0.0000120687092035 );
            ir->AddPentPoint( 158 , 0.0002635287125677 , 0.7337151743454349 , 0.0002635287125677 , 0.1328788841147149 , 0.0000120687092035 );
            ir->AddPentPoint( 159 , 0.7337151743454349 , 0.1328788841147149 , 0.0002635287125677 , 0.1328788841147149 , 0.0000120687092035 );
            ir->AddPentPoint( 160 , 0.0002635287125677 , 0.1328788841147149 , 0.0002635287125677 , 0.1328788841147149 , 0.0000120687092035 );
            ir->AddPentPoint( 161 , 0.7337151743454349 , 0.0002635287125677 , 0.0002635287125677 , 0.1328788841147149 , 0.0000120687092035 );
            ir->AddPentPoint( 162 , 0.1328788841147149 , 0.0002635287125677 , 0.0002635287125677 , 0.1328788841147149 , 0.0000120687092035 );
            ir->AddPentPoint( 163 , 0.1328788841147149 , 0.1328788841147149 , 0.7337151743454349 , 0.0002635287125677 , 0.0000120687092035 );
            ir->AddPentPoint( 164 , 0.0002635287125677 , 0.1328788841147149 , 0.7337151743454349 , 0.0002635287125677 , 0.0000120687092035 );
            ir->AddPentPoint( 165 , 0.1328788841147149 , 0.0002635287125677 , 0.7337151743454349 , 0.0002635287125677 , 0.0000120687092035 );
            ir->AddPentPoint( 166 , 0.1328788841147149 , 0.7337151743454349 , 0.1328788841147149 , 0.0002635287125677 , 0.0000120687092035 );
            ir->AddPentPoint( 167 , 0.0002635287125677 , 0.7337151743454349 , 0.1328788841147149 , 0.0002635287125677 , 0.0000120687092035 );
            ir->AddPentPoint( 168 , 0.7337151743454349 , 0.1328788841147149 , 0.1328788841147149 , 0.0002635287125677 , 0.0000120687092035 );
            ir->AddPentPoint( 169 , 0.0002635287125677 , 0.1328788841147149 , 0.1328788841147149 , 0.0002635287125677 , 0.0000120687092035 );
            ir->AddPentPoint( 170 , 0.7337151743454349 , 0.0002635287125677 , 0.1328788841147149 , 0.0002635287125677 , 0.0000120687092035 );
            ir->AddPentPoint( 171 , 0.1328788841147149 , 0.0002635287125677 , 0.1328788841147149 , 0.0002635287125677 , 0.0000120687092035 );
            ir->AddPentPoint( 172 , 0.1328788841147149 , 0.7337151743454349 , 0.0002635287125677 , 0.0002635287125677 , 0.0000120687092035 );
            ir->AddPentPoint( 173 , 0.7337151743454349 , 0.1328788841147149 , 0.0002635287125677 , 0.0002635287125677 , 0.0000120687092035 );
            ir->AddPentPoint( 174 , 0.1328788841147149 , 0.1328788841147149 , 0.0002635287125677 , 0.0002635287125677 , 0.0000120687092035 );
            ir->AddPentPoint( 175 , 0.0641103815880093 , 0.1835311602251138 , 0.1835311602251138 , 0.5047169163737538 , 0.0002192584644786 );
            ir->AddPentPoint( 176 , 0.1835311602251138 , 0.0641103815880093 , 0.1835311602251138 , 0.5047169163737538 , 0.0002192584644786 );
            ir->AddPentPoint( 177 , 0.0641103815880093 , 0.0641103815880093 , 0.1835311602251138 , 0.5047169163737538 , 0.0002192584644786 );
            ir->AddPentPoint( 178 , 0.1835311602251138 , 0.1835311602251138 , 0.0641103815880093 , 0.5047169163737538 , 0.0002192584644786 );
            ir->AddPentPoint( 179 , 0.0641103815880093 , 0.1835311602251138 , 0.0641103815880093 , 0.5047169163737538 , 0.0002192584644786 );
            ir->AddPentPoint( 180 , 0.1835311602251138 , 0.0641103815880093 , 0.0641103815880093 , 0.5047169163737538 , 0.0002192584644786 );
            ir->AddPentPoint( 181 , 0.0641103815880093 , 0.1835311602251138 , 0.5047169163737538 , 0.1835311602251138 , 0.0002192584644786 );
            ir->AddPentPoint( 182 , 0.1835311602251138 , 0.0641103815880093 , 0.5047169163737538 , 0.1835311602251138 , 0.0002192584644786 );
            ir->AddPentPoint( 183 , 0.0641103815880093 , 0.0641103815880093 , 0.5047169163737538 , 0.1835311602251138 , 0.0002192584644786 );
            ir->AddPentPoint( 184 , 0.0641103815880093 , 0.5047169163737538 , 0.1835311602251138 , 0.1835311602251138 , 0.0002192584644786 );
            ir->AddPentPoint( 185 , 0.5047169163737538 , 0.0641103815880093 , 0.1835311602251138 , 0.1835311602251138 , 0.0002192584644786 );
            ir->AddPentPoint( 186 , 0.0641103815880093 , 0.0641103815880093 , 0.1835311602251138 , 0.1835311602251138 , 0.0002192584644786 );
            ir->AddPentPoint( 187 , 0.1835311602251138 , 0.5047169163737538 , 0.0641103815880093 , 0.1835311602251138 , 0.0002192584644786 );
            ir->AddPentPoint( 188 , 0.0641103815880093 , 0.5047169163737538 , 0.0641103815880093 , 0.1835311602251138 , 0.0002192584644786 );
            ir->AddPentPoint( 189 , 0.5047169163737538 , 0.1835311602251138 , 0.0641103815880093 , 0.1835311602251138 , 0.0002192584644786 );
            ir->AddPentPoint( 190 , 0.0641103815880093 , 0.1835311602251138 , 0.0641103815880093 , 0.1835311602251138 , 0.0002192584644786 );
            ir->AddPentPoint( 191 , 0.5047169163737538 , 0.0641103815880093 , 0.0641103815880093 , 0.1835311602251138 , 0.0002192584644786 );
            ir->AddPentPoint( 192 , 0.1835311602251138 , 0.0641103815880093 , 0.0641103815880093 , 0.1835311602251138 , 0.0002192584644786 );
            ir->AddPentPoint( 193 , 0.1835311602251138 , 0.1835311602251138 , 0.5047169163737538 , 0.0641103815880093 , 0.0002192584644786 );
            ir->AddPentPoint( 194 , 0.0641103815880093 , 0.1835311602251138 , 0.5047169163737538 , 0.0641103815880093 , 0.0002192584644786 );
            ir->AddPentPoint( 195 , 0.1835311602251138 , 0.0641103815880093 , 0.5047169163737538 , 0.0641103815880093 , 0.0002192584644786 );
            ir->AddPentPoint( 196 , 0.1835311602251138 , 0.5047169163737538 , 0.1835311602251138 , 0.0641103815880093 , 0.0002192584644786 );
            ir->AddPentPoint( 197 , 0.0641103815880093 , 0.5047169163737538 , 0.1835311602251138 , 0.0641103815880093 , 0.0002192584644786 );
            ir->AddPentPoint( 198 , 0.5047169163737538 , 0.1835311602251138 , 0.1835311602251138 , 0.0641103815880093 , 0.0002192584644786 );
            ir->AddPentPoint( 199 , 0.0641103815880093 , 0.1835311602251138 , 0.1835311602251138 , 0.0641103815880093 , 0.0002192584644786 );
            ir->AddPentPoint( 200 , 0.5047169163737538 , 0.0641103815880093 , 0.1835311602251138 , 0.0641103815880093 , 0.0002192584644786 );
            ir->AddPentPoint( 201 , 0.1835311602251138 , 0.0641103815880093 , 0.1835311602251138 , 0.0641103815880093 , 0.0002192584644786 );
            ir->AddPentPoint( 202 , 0.1835311602251138 , 0.5047169163737538 , 0.0641103815880093 , 0.0641103815880093 , 0.0002192584644786 );
            ir->AddPentPoint( 203 , 0.5047169163737538 , 0.1835311602251138 , 0.0641103815880093 , 0.0641103815880093 , 0.0002192584644786 );
            ir->AddPentPoint( 204 , 0.1835311602251138 , 0.1835311602251138 , 0.0641103815880093 , 0.0641103815880093 , 0.0002192584644786 );
            ir->AddPentPoint( 205 , 0.0229833630026369 , 0.1342981602908115 , 0.3585906320589499 , 0.4611444816449648 , 0.0000595930802549 );
            ir->AddPentPoint( 206 , 0.1342981602908115 , 0.0229833630026369 , 0.3585906320589499 , 0.4611444816449648 , 0.0000595930802549 );
            ir->AddPentPoint( 207 , 0.0229833630026369 , 0.0229833630026369 , 0.3585906320589499 , 0.4611444816449648 , 0.0000595930802549 );
            ir->AddPentPoint( 208 , 0.0229833630026369 , 0.3585906320589499 , 0.1342981602908115 , 0.4611444816449648 , 0.0000595930802549 );
            ir->AddPentPoint( 209 , 0.3585906320589499 , 0.0229833630026369 , 0.1342981602908115 , 0.4611444816449648 , 0.0000595930802549 );
            ir->AddPentPoint( 210 , 0.0229833630026369 , 0.0229833630026369 , 0.1342981602908115 , 0.4611444816449648 , 0.0000595930802549 );
            ir->AddPentPoint( 211 , 0.1342981602908115 , 0.3585906320589499 , 0.0229833630026369 , 0.4611444816449648 , 0.0000595930802549 );
            ir->AddPentPoint( 212 , 0.0229833630026369 , 0.3585906320589499 , 0.0229833630026369 , 0.4611444816449648 , 0.0000595930802549 );
            ir->AddPentPoint( 213 , 0.3585906320589499 , 0.1342981602908115 , 0.0229833630026369 , 0.4611444816449648 , 0.0000595930802549 );
            ir->AddPentPoint( 214 , 0.0229833630026369 , 0.1342981602908115 , 0.0229833630026369 , 0.4611444816449648 , 0.0000595930802549 );
            ir->AddPentPoint( 215 , 0.3585906320589499 , 0.0229833630026369 , 0.0229833630026369 , 0.4611444816449648 , 0.0000595930802549 );
            ir->AddPentPoint( 216 , 0.1342981602908115 , 0.0229833630026369 , 0.0229833630026369 , 0.4611444816449648 , 0.0000595930802549 );
            ir->AddPentPoint( 217 , 0.0229833630026369 , 0.1342981602908115 , 0.4611444816449648 , 0.3585906320589499 , 0.0000595930802549 );
            ir->AddPentPoint( 218 , 0.1342981602908115 , 0.0229833630026369 , 0.4611444816449648 , 0.3585906320589499 , 0.0000595930802549 );
            ir->AddPentPoint( 219 , 0.0229833630026369 , 0.0229833630026369 , 0.4611444816449648 , 0.3585906320589499 , 0.0000595930802549 );
            ir->AddPentPoint( 220 , 0.0229833630026369 , 0.4611444816449648 , 0.1342981602908115 , 0.3585906320589499 , 0.0000595930802549 );
            ir->AddPentPoint( 221 , 0.4611444816449648 , 0.0229833630026369 , 0.1342981602908115 , 0.3585906320589499 , 0.0000595930802549 );
            ir->AddPentPoint( 222 , 0.0229833630026369 , 0.0229833630026369 , 0.1342981602908115 , 0.3585906320589499 , 0.0000595930802549 );
            ir->AddPentPoint( 223 , 0.1342981602908115 , 0.4611444816449648 , 0.0229833630026369 , 0.3585906320589499 , 0.0000595930802549 );
            ir->AddPentPoint( 224 , 0.0229833630026369 , 0.4611444816449648 , 0.0229833630026369 , 0.3585906320589499 , 0.0000595930802549 );
            ir->AddPentPoint( 225 , 0.4611444816449648 , 0.1342981602908115 , 0.0229833630026369 , 0.3585906320589499 , 0.0000595930802549 );
            ir->AddPentPoint( 226 , 0.0229833630026369 , 0.1342981602908115 , 0.0229833630026369 , 0.3585906320589499 , 0.0000595930802549 );
            ir->AddPentPoint( 227 , 0.4611444816449648 , 0.0229833630026369 , 0.0229833630026369 , 0.3585906320589499 , 0.0000595930802549 );
            ir->AddPentPoint( 228 , 0.1342981602908115 , 0.0229833630026369 , 0.0229833630026369 , 0.3585906320589499 , 0.0000595930802549 );
            ir->AddPentPoint( 229 , 0.0229833630026369 , 0.3585906320589499 , 0.4611444816449648 , 0.1342981602908115 , 0.0000595930802549 );
            ir->AddPentPoint( 230 , 0.3585906320589499 , 0.0229833630026369 , 0.4611444816449648 , 0.1342981602908115 , 0.0000595930802549 );
            ir->AddPentPoint( 231 , 0.0229833630026369 , 0.0229833630026369 , 0.4611444816449648 , 0.1342981602908115 , 0.0000595930802549 );
            ir->AddPentPoint( 232 , 0.0229833630026369 , 0.4611444816449648 , 0.3585906320589499 , 0.1342981602908115 , 0.0000595930802549 );
            ir->AddPentPoint( 233 , 0.4611444816449648 , 0.0229833630026369 , 0.3585906320589499 , 0.1342981602908115 , 0.0000595930802549 );
            ir->AddPentPoint( 234 , 0.0229833630026369 , 0.0229833630026369 , 0.3585906320589499 , 0.1342981602908115 , 0.0000595930802549 );
            ir->AddPentPoint( 235 , 0.3585906320589499 , 0.4611444816449648 , 0.0229833630026369 , 0.1342981602908115 , 0.0000595930802549 );
            ir->AddPentPoint( 236 , 0.0229833630026369 , 0.4611444816449648 , 0.0229833630026369 , 0.1342981602908115 , 0.0000595930802549 );
            ir->AddPentPoint( 237 , 0.4611444816449648 , 0.3585906320589499 , 0.0229833630026369 , 0.1342981602908115 , 0.0000595930802549 );
            ir->AddPentPoint( 238 , 0.0229833630026369 , 0.3585906320589499 , 0.0229833630026369 , 0.1342981602908115 , 0.0000595930802549 );
            ir->AddPentPoint( 239 , 0.4611444816449648 , 0.0229833630026369 , 0.0229833630026369 , 0.1342981602908115 , 0.0000595930802549 );
            ir->AddPentPoint( 240 , 0.3585906320589499 , 0.0229833630026369 , 0.0229833630026369 , 0.1342981602908115 , 0.0000595930802549 );
            ir->AddPentPoint( 241 , 0.1342981602908115 , 0.3585906320589499 , 0.4611444816449648 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 242 , 0.0229833630026369 , 0.3585906320589499 , 0.4611444816449648 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 243 , 0.3585906320589499 , 0.1342981602908115 , 0.4611444816449648 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 244 , 0.0229833630026369 , 0.1342981602908115 , 0.4611444816449648 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 245 , 0.3585906320589499 , 0.0229833630026369 , 0.4611444816449648 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 246 , 0.1342981602908115 , 0.0229833630026369 , 0.4611444816449648 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 247 , 0.1342981602908115 , 0.4611444816449648 , 0.3585906320589499 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 248 , 0.0229833630026369 , 0.4611444816449648 , 0.3585906320589499 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 249 , 0.4611444816449648 , 0.1342981602908115 , 0.3585906320589499 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 250 , 0.0229833630026369 , 0.1342981602908115 , 0.3585906320589499 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 251 , 0.4611444816449648 , 0.0229833630026369 , 0.3585906320589499 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 252 , 0.1342981602908115 , 0.0229833630026369 , 0.3585906320589499 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 253 , 0.3585906320589499 , 0.4611444816449648 , 0.1342981602908115 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 254 , 0.0229833630026369 , 0.4611444816449648 , 0.1342981602908115 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 255 , 0.4611444816449648 , 0.3585906320589499 , 0.1342981602908115 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 256 , 0.0229833630026369 , 0.3585906320589499 , 0.1342981602908115 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 257 , 0.4611444816449648 , 0.0229833630026369 , 0.1342981602908115 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 258 , 0.3585906320589499 , 0.0229833630026369 , 0.1342981602908115 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 259 , 0.3585906320589499 , 0.4611444816449648 , 0.0229833630026369 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 260 , 0.1342981602908115 , 0.4611444816449648 , 0.0229833630026369 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 261 , 0.4611444816449648 , 0.3585906320589499 , 0.0229833630026369 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 262 , 0.1342981602908115 , 0.3585906320589499 , 0.0229833630026369 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 263 , 0.4611444816449648 , 0.1342981602908115 , 0.0229833630026369 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 264 , 0.3585906320589499 , 0.1342981602908115 , 0.0229833630026369 , 0.0229833630026369 , 0.0000595930802549 );
            ir->AddPentPoint( 265 , 0.0577612562875865 , 0.0047476858698278 , 0.2741151471837904 , 0.6056146543712088 , 0.0000436847082236 );
            ir->AddPentPoint( 266 , 0.0047476858698278 , 0.0577612562875865 , 0.2741151471837904 , 0.6056146543712088 , 0.0000436847082236 );
            ir->AddPentPoint( 267 , 0.0577612562875865 , 0.0577612562875865 , 0.2741151471837904 , 0.6056146543712088 , 0.0000436847082236 );
            ir->AddPentPoint( 268 , 0.0577612562875865 , 0.2741151471837904 , 0.0047476858698278 , 0.6056146543712088 , 0.0000436847082236 );
            ir->AddPentPoint( 269 , 0.2741151471837904 , 0.0577612562875865 , 0.0047476858698278 , 0.6056146543712088 , 0.0000436847082236 );
            ir->AddPentPoint( 270 , 0.0577612562875865 , 0.0577612562875865 , 0.0047476858698278 , 0.6056146543712088 , 0.0000436847082236 );
            ir->AddPentPoint( 271 , 0.0047476858698278 , 0.2741151471837904 , 0.0577612562875865 , 0.6056146543712088 , 0.0000436847082236 );
            ir->AddPentPoint( 272 , 0.0577612562875865 , 0.2741151471837904 , 0.0577612562875865 , 0.6056146543712088 , 0.0000436847082236 );
            ir->AddPentPoint( 273 , 0.2741151471837904 , 0.0047476858698278 , 0.0577612562875865 , 0.6056146543712088 , 0.0000436847082236 );
            ir->AddPentPoint( 274 , 0.0577612562875865 , 0.0047476858698278 , 0.0577612562875865 , 0.6056146543712088 , 0.0000436847082236 );
            ir->AddPentPoint( 275 , 0.2741151471837904 , 0.0577612562875865 , 0.0577612562875865 , 0.6056146543712088 , 0.0000436847082236 );
            ir->AddPentPoint( 276 , 0.0047476858698278 , 0.0577612562875865 , 0.0577612562875865 , 0.6056146543712088 , 0.0000436847082236 );
            ir->AddPentPoint( 277 , 0.0577612562875865 , 0.0047476858698278 , 0.6056146543712088 , 0.2741151471837904 , 0.0000436847082236 );
            ir->AddPentPoint( 278 , 0.0047476858698278 , 0.0577612562875865 , 0.6056146543712088 , 0.2741151471837904 , 0.0000436847082236 );
            ir->AddPentPoint( 279 , 0.0577612562875865 , 0.0577612562875865 , 0.6056146543712088 , 0.2741151471837904 , 0.0000436847082236 );
            ir->AddPentPoint( 280 , 0.0577612562875865 , 0.6056146543712088 , 0.0047476858698278 , 0.2741151471837904 , 0.0000436847082236 );
            ir->AddPentPoint( 281 , 0.6056146543712088 , 0.0577612562875865 , 0.0047476858698278 , 0.2741151471837904 , 0.0000436847082236 );
            ir->AddPentPoint( 282 , 0.0577612562875865 , 0.0577612562875865 , 0.0047476858698278 , 0.2741151471837904 , 0.0000436847082236 );
            ir->AddPentPoint( 283 , 0.0047476858698278 , 0.6056146543712088 , 0.0577612562875865 , 0.2741151471837904 , 0.0000436847082236 );
            ir->AddPentPoint( 284 , 0.0577612562875865 , 0.6056146543712088 , 0.0577612562875865 , 0.2741151471837904 , 0.0000436847082236 );
            ir->AddPentPoint( 285 , 0.6056146543712088 , 0.0047476858698278 , 0.0577612562875865 , 0.2741151471837904 , 0.0000436847082236 );
            ir->AddPentPoint( 286 , 0.0577612562875865 , 0.0047476858698278 , 0.0577612562875865 , 0.2741151471837904 , 0.0000436847082236 );
            ir->AddPentPoint( 287 , 0.6056146543712088 , 0.0577612562875865 , 0.0577612562875865 , 0.2741151471837904 , 0.0000436847082236 );
            ir->AddPentPoint( 288 , 0.0047476858698278 , 0.0577612562875865 , 0.0577612562875865 , 0.2741151471837904 , 0.0000436847082236 );
            ir->AddPentPoint( 289 , 0.0577612562875865 , 0.2741151471837904 , 0.6056146543712088 , 0.0047476858698278 , 0.0000436847082236 );
            ir->AddPentPoint( 290 , 0.2741151471837904 , 0.0577612562875865 , 0.6056146543712088 , 0.0047476858698278 , 0.0000436847082236 );
            ir->AddPentPoint( 291 , 0.0577612562875865 , 0.0577612562875865 , 0.6056146543712088 , 0.0047476858698278 , 0.0000436847082236 );
            ir->AddPentPoint( 292 , 0.0577612562875865 , 0.6056146543712088 , 0.2741151471837904 , 0.0047476858698278 , 0.0000436847082236 );
            ir->AddPentPoint( 293 , 0.6056146543712088 , 0.0577612562875865 , 0.2741151471837904 , 0.0047476858698278 , 0.0000436847082236 );
            ir->AddPentPoint( 294 , 0.0577612562875865 , 0.0577612562875865 , 0.2741151471837904 , 0.0047476858698278 , 0.0000436847082236 );
            ir->AddPentPoint( 295 , 0.2741151471837904 , 0.6056146543712088 , 0.0577612562875865 , 0.0047476858698278 , 0.0000436847082236 );
            ir->AddPentPoint( 296 , 0.0577612562875865 , 0.6056146543712088 , 0.0577612562875865 , 0.0047476858698278 , 0.0000436847082236 );
            ir->AddPentPoint( 297 , 0.6056146543712088 , 0.2741151471837904 , 0.0577612562875865 , 0.0047476858698278 , 0.0000436847082236 );
            ir->AddPentPoint( 298 , 0.0577612562875865 , 0.2741151471837904 , 0.0577612562875865 , 0.0047476858698278 , 0.0000436847082236 );
            ir->AddPentPoint( 299 , 0.6056146543712088 , 0.0577612562875865 , 0.0577612562875865 , 0.0047476858698278 , 0.0000436847082236 );
            ir->AddPentPoint( 300 , 0.2741151471837904 , 0.0577612562875865 , 0.0577612562875865 , 0.0047476858698278 , 0.0000436847082236 );
            ir->AddPentPoint( 301 , 0.0047476858698278 , 0.2741151471837904 , 0.6056146543712088 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 302 , 0.0577612562875865 , 0.2741151471837904 , 0.6056146543712088 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 303 , 0.2741151471837904 , 0.0047476858698278 , 0.6056146543712088 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 304 , 0.0577612562875865 , 0.0047476858698278 , 0.6056146543712088 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 305 , 0.2741151471837904 , 0.0577612562875865 , 0.6056146543712088 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 306 , 0.0047476858698278 , 0.0577612562875865 , 0.6056146543712088 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 307 , 0.0047476858698278 , 0.6056146543712088 , 0.2741151471837904 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 308 , 0.0577612562875865 , 0.6056146543712088 , 0.2741151471837904 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 309 , 0.6056146543712088 , 0.0047476858698278 , 0.2741151471837904 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 310 , 0.0577612562875865 , 0.0047476858698278 , 0.2741151471837904 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 311 , 0.6056146543712088 , 0.0577612562875865 , 0.2741151471837904 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 312 , 0.0047476858698278 , 0.0577612562875865 , 0.2741151471837904 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 313 , 0.2741151471837904 , 0.6056146543712088 , 0.0047476858698278 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 314 , 0.0577612562875865 , 0.6056146543712088 , 0.0047476858698278 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 315 , 0.6056146543712088 , 0.2741151471837904 , 0.0047476858698278 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 316 , 0.0577612562875865 , 0.2741151471837904 , 0.0047476858698278 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 317 , 0.6056146543712088 , 0.0577612562875865 , 0.0047476858698278 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 318 , 0.2741151471837904 , 0.0577612562875865 , 0.0047476858698278 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 319 , 0.2741151471837904 , 0.6056146543712088 , 0.0577612562875865 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 320 , 0.0047476858698278 , 0.6056146543712088 , 0.0577612562875865 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 321 , 0.6056146543712088 , 0.2741151471837904 , 0.0577612562875865 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 322 , 0.0047476858698278 , 0.2741151471837904 , 0.0577612562875865 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 323 , 0.6056146543712088 , 0.0047476858698278 , 0.0577612562875865 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 324 , 0.2741151471837904 , 0.0047476858698278 , 0.0577612562875865 , 0.0577612562875865 , 0.0000436847082236 );
            ir->AddPentPoint( 325 , 0.1420256473819541 , 0.0171187111516862 , 0.2762455471960004 , 0.4225844468884052 , 0.0001395478845176 );
            ir->AddPentPoint( 326 , 0.0171187111516862 , 0.1420256473819541 , 0.2762455471960004 , 0.4225844468884052 , 0.0001395478845176 );
            ir->AddPentPoint( 327 , 0.1420256473819541 , 0.1420256473819541 , 0.2762455471960004 , 0.4225844468884052 , 0.0001395478845176 );
            ir->AddPentPoint( 328 , 0.1420256473819541 , 0.2762455471960004 , 0.0171187111516862 , 0.4225844468884052 , 0.0001395478845176 );
            ir->AddPentPoint( 329 , 0.2762455471960004 , 0.1420256473819541 , 0.0171187111516862 , 0.4225844468884052 , 0.0001395478845176 );
            ir->AddPentPoint( 330 , 0.1420256473819541 , 0.1420256473819541 , 0.0171187111516862 , 0.4225844468884052 , 0.0001395478845176 );
            ir->AddPentPoint( 331 , 0.0171187111516862 , 0.2762455471960004 , 0.1420256473819541 , 0.4225844468884052 , 0.0001395478845176 );
            ir->AddPentPoint( 332 , 0.1420256473819541 , 0.2762455471960004 , 0.1420256473819541 , 0.4225844468884052 , 0.0001395478845176 );
            ir->AddPentPoint( 333 , 0.2762455471960004 , 0.0171187111516862 , 0.1420256473819541 , 0.4225844468884052 , 0.0001395478845176 );
            ir->AddPentPoint( 334 , 0.1420256473819541 , 0.0171187111516862 , 0.1420256473819541 , 0.4225844468884052 , 0.0001395478845176 );
            ir->AddPentPoint( 335 , 0.2762455471960004 , 0.1420256473819541 , 0.1420256473819541 , 0.4225844468884052 , 0.0001395478845176 );
            ir->AddPentPoint( 336 , 0.0171187111516862 , 0.1420256473819541 , 0.1420256473819541 , 0.4225844468884052 , 0.0001395478845176 );
            ir->AddPentPoint( 337 , 0.1420256473819541 , 0.0171187111516862 , 0.4225844468884052 , 0.2762455471960004 , 0.0001395478845176 );
            ir->AddPentPoint( 338 , 0.0171187111516862 , 0.1420256473819541 , 0.4225844468884052 , 0.2762455471960004 , 0.0001395478845176 );
            ir->AddPentPoint( 339 , 0.1420256473819541 , 0.1420256473819541 , 0.4225844468884052 , 0.2762455471960004 , 0.0001395478845176 );
            ir->AddPentPoint( 340 , 0.1420256473819541 , 0.4225844468884052 , 0.0171187111516862 , 0.2762455471960004 , 0.0001395478845176 );
            ir->AddPentPoint( 341 , 0.4225844468884052 , 0.1420256473819541 , 0.0171187111516862 , 0.2762455471960004 , 0.0001395478845176 );
            ir->AddPentPoint( 342 , 0.1420256473819541 , 0.1420256473819541 , 0.0171187111516862 , 0.2762455471960004 , 0.0001395478845176 );
            ir->AddPentPoint( 343 , 0.0171187111516862 , 0.4225844468884052 , 0.1420256473819541 , 0.2762455471960004 , 0.0001395478845176 );
            ir->AddPentPoint( 344 , 0.1420256473819541 , 0.4225844468884052 , 0.1420256473819541 , 0.2762455471960004 , 0.0001395478845176 );
            ir->AddPentPoint( 345 , 0.4225844468884052 , 0.0171187111516862 , 0.1420256473819541 , 0.2762455471960004 , 0.0001395478845176 );
            ir->AddPentPoint( 346 , 0.1420256473819541 , 0.0171187111516862 , 0.1420256473819541 , 0.2762455471960004 , 0.0001395478845176 );
            ir->AddPentPoint( 347 , 0.4225844468884052 , 0.1420256473819541 , 0.1420256473819541 , 0.2762455471960004 , 0.0001395478845176 );
            ir->AddPentPoint( 348 , 0.0171187111516862 , 0.1420256473819541 , 0.1420256473819541 , 0.2762455471960004 , 0.0001395478845176 );
            ir->AddPentPoint( 349 , 0.1420256473819541 , 0.2762455471960004 , 0.4225844468884052 , 0.0171187111516862 , 0.0001395478845176 );
            ir->AddPentPoint( 350 , 0.2762455471960004 , 0.1420256473819541 , 0.4225844468884052 , 0.0171187111516862 , 0.0001395478845176 );
            ir->AddPentPoint( 351 , 0.1420256473819541 , 0.1420256473819541 , 0.4225844468884052 , 0.0171187111516862 , 0.0001395478845176 );
            ir->AddPentPoint( 352 , 0.1420256473819541 , 0.4225844468884052 , 0.2762455471960004 , 0.0171187111516862 , 0.0001395478845176 );
            ir->AddPentPoint( 353 , 0.4225844468884052 , 0.1420256473819541 , 0.2762455471960004 , 0.0171187111516862 , 0.0001395478845176 );
            ir->AddPentPoint( 354 , 0.1420256473819541 , 0.1420256473819541 , 0.2762455471960004 , 0.0171187111516862 , 0.0001395478845176 );
            ir->AddPentPoint( 355 , 0.2762455471960004 , 0.4225844468884052 , 0.1420256473819541 , 0.0171187111516862 , 0.0001395478845176 );
            ir->AddPentPoint( 356 , 0.1420256473819541 , 0.4225844468884052 , 0.1420256473819541 , 0.0171187111516862 , 0.0001395478845176 );
            ir->AddPentPoint( 357 , 0.4225844468884052 , 0.2762455471960004 , 0.1420256473819541 , 0.0171187111516862 , 0.0001395478845176 );
            ir->AddPentPoint( 358 , 0.1420256473819541 , 0.2762455471960004 , 0.1420256473819541 , 0.0171187111516862 , 0.0001395478845176 );
            ir->AddPentPoint( 359 , 0.4225844468884052 , 0.1420256473819541 , 0.1420256473819541 , 0.0171187111516862 , 0.0001395478845176 );
            ir->AddPentPoint( 360 , 0.2762455471960004 , 0.1420256473819541 , 0.1420256473819541 , 0.0171187111516862 , 0.0001395478845176 );
            ir->AddPentPoint( 361 , 0.0171187111516862 , 0.2762455471960004 , 0.4225844468884052 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 362 , 0.1420256473819541 , 0.2762455471960004 , 0.4225844468884052 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 363 , 0.2762455471960004 , 0.0171187111516862 , 0.4225844468884052 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 364 , 0.1420256473819541 , 0.0171187111516862 , 0.4225844468884052 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 365 , 0.2762455471960004 , 0.1420256473819541 , 0.4225844468884052 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 366 , 0.0171187111516862 , 0.1420256473819541 , 0.4225844468884052 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 367 , 0.0171187111516862 , 0.4225844468884052 , 0.2762455471960004 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 368 , 0.1420256473819541 , 0.4225844468884052 , 0.2762455471960004 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 369 , 0.4225844468884052 , 0.0171187111516862 , 0.2762455471960004 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 370 , 0.1420256473819541 , 0.0171187111516862 , 0.2762455471960004 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 371 , 0.4225844468884052 , 0.1420256473819541 , 0.2762455471960004 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 372 , 0.0171187111516862 , 0.1420256473819541 , 0.2762455471960004 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 373 , 0.2762455471960004 , 0.4225844468884052 , 0.0171187111516862 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 374 , 0.1420256473819541 , 0.4225844468884052 , 0.0171187111516862 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 375 , 0.4225844468884052 , 0.2762455471960004 , 0.0171187111516862 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 376 , 0.1420256473819541 , 0.2762455471960004 , 0.0171187111516862 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 377 , 0.4225844468884052 , 0.1420256473819541 , 0.0171187111516862 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 378 , 0.2762455471960004 , 0.1420256473819541 , 0.0171187111516862 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 379 , 0.2762455471960004 , 0.4225844468884052 , 0.1420256473819541 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 380 , 0.0171187111516862 , 0.4225844468884052 , 0.1420256473819541 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 381 , 0.4225844468884052 , 0.2762455471960004 , 0.1420256473819541 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 382 , 0.0171187111516862 , 0.2762455471960004 , 0.1420256473819541 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 383 , 0.4225844468884052 , 0.0171187111516862 , 0.1420256473819541 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 384 , 0.2762455471960004 , 0.0171187111516862 , 0.1420256473819541 , 0.1420256473819541 , 0.0001395478845176 );
            ir->AddPentPoint( 385 , 0.2096943241276536 , 0.0062247920880195 , 0.0364941472831919 , 0.5378924123734813 , 0.0000326589757115 );
            ir->AddPentPoint( 386 , 0.0062247920880195 , 0.2096943241276536 , 0.0364941472831919 , 0.5378924123734813 , 0.0000326589757115 );
            ir->AddPentPoint( 387 , 0.2096943241276536 , 0.2096943241276536 , 0.0364941472831919 , 0.5378924123734813 , 0.0000326589757115 );
            ir->AddPentPoint( 388 , 0.2096943241276536 , 0.0364941472831919 , 0.0062247920880195 , 0.5378924123734813 , 0.0000326589757115 );
            ir->AddPentPoint( 389 , 0.0364941472831919 , 0.2096943241276536 , 0.0062247920880195 , 0.5378924123734813 , 0.0000326589757115 );
            ir->AddPentPoint( 390 , 0.2096943241276536 , 0.2096943241276536 , 0.0062247920880195 , 0.5378924123734813 , 0.0000326589757115 );
            ir->AddPentPoint( 391 , 0.0062247920880195 , 0.0364941472831919 , 0.2096943241276536 , 0.5378924123734813 , 0.0000326589757115 );
            ir->AddPentPoint( 392 , 0.2096943241276536 , 0.0364941472831919 , 0.2096943241276536 , 0.5378924123734813 , 0.0000326589757115 );
            ir->AddPentPoint( 393 , 0.0364941472831919 , 0.0062247920880195 , 0.2096943241276536 , 0.5378924123734813 , 0.0000326589757115 );
            ir->AddPentPoint( 394 , 0.2096943241276536 , 0.0062247920880195 , 0.2096943241276536 , 0.5378924123734813 , 0.0000326589757115 );
            ir->AddPentPoint( 395 , 0.0364941472831919 , 0.2096943241276536 , 0.2096943241276536 , 0.5378924123734813 , 0.0000326589757115 );
            ir->AddPentPoint( 396 , 0.0062247920880195 , 0.2096943241276536 , 0.2096943241276536 , 0.5378924123734813 , 0.0000326589757115 );
            ir->AddPentPoint( 397 , 0.2096943241276536 , 0.0062247920880195 , 0.5378924123734813 , 0.0364941472831919 , 0.0000326589757115 );
            ir->AddPentPoint( 398 , 0.0062247920880195 , 0.2096943241276536 , 0.5378924123734813 , 0.0364941472831919 , 0.0000326589757115 );
            ir->AddPentPoint( 399 , 0.2096943241276536 , 0.2096943241276536 , 0.5378924123734813 , 0.0364941472831919 , 0.0000326589757115 );
            ir->AddPentPoint( 400 , 0.2096943241276536 , 0.5378924123734813 , 0.0062247920880195 , 0.0364941472831919 , 0.0000326589757115 );
            ir->AddPentPoint( 401 , 0.5378924123734813 , 0.2096943241276536 , 0.0062247920880195 , 0.0364941472831919 , 0.0000326589757115 );
            ir->AddPentPoint( 402 , 0.2096943241276536 , 0.2096943241276536 , 0.0062247920880195 , 0.0364941472831919 , 0.0000326589757115 );
            ir->AddPentPoint( 403 , 0.0062247920880195 , 0.5378924123734813 , 0.2096943241276536 , 0.0364941472831919 , 0.0000326589757115 );
            ir->AddPentPoint( 404 , 0.2096943241276536 , 0.5378924123734813 , 0.2096943241276536 , 0.0364941472831919 , 0.0000326589757115 );
            ir->AddPentPoint( 405 , 0.5378924123734813 , 0.0062247920880195 , 0.2096943241276536 , 0.0364941472831919 , 0.0000326589757115 );
            ir->AddPentPoint( 406 , 0.2096943241276536 , 0.0062247920880195 , 0.2096943241276536 , 0.0364941472831919 , 0.0000326589757115 );
            ir->AddPentPoint( 407 , 0.5378924123734813 , 0.2096943241276536 , 0.2096943241276536 , 0.0364941472831919 , 0.0000326589757115 );
            ir->AddPentPoint( 408 , 0.0062247920880195 , 0.2096943241276536 , 0.2096943241276536 , 0.0364941472831919 , 0.0000326589757115 );
            ir->AddPentPoint( 409 , 0.2096943241276536 , 0.0364941472831919 , 0.5378924123734813 , 0.0062247920880195 , 0.0000326589757115 );
            ir->AddPentPoint( 410 , 0.0364941472831919 , 0.2096943241276536 , 0.5378924123734813 , 0.0062247920880195 , 0.0000326589757115 );
            ir->AddPentPoint( 411 , 0.2096943241276536 , 0.2096943241276536 , 0.5378924123734813 , 0.0062247920880195 , 0.0000326589757115 );
            ir->AddPentPoint( 412 , 0.2096943241276536 , 0.5378924123734813 , 0.0364941472831919 , 0.0062247920880195 , 0.0000326589757115 );
            ir->AddPentPoint( 413 , 0.5378924123734813 , 0.2096943241276536 , 0.0364941472831919 , 0.0062247920880195 , 0.0000326589757115 );
            ir->AddPentPoint( 414 , 0.2096943241276536 , 0.2096943241276536 , 0.0364941472831919 , 0.0062247920880195 , 0.0000326589757115 );
            ir->AddPentPoint( 415 , 0.0364941472831919 , 0.5378924123734813 , 0.2096943241276536 , 0.0062247920880195 , 0.0000326589757115 );
            ir->AddPentPoint( 416 , 0.2096943241276536 , 0.5378924123734813 , 0.2096943241276536 , 0.0062247920880195 , 0.0000326589757115 );
            ir->AddPentPoint( 417 , 0.5378924123734813 , 0.0364941472831919 , 0.2096943241276536 , 0.0062247920880195 , 0.0000326589757115 );
            ir->AddPentPoint( 418 , 0.2096943241276536 , 0.0364941472831919 , 0.2096943241276536 , 0.0062247920880195 , 0.0000326589757115 );
            ir->AddPentPoint( 419 , 0.5378924123734813 , 0.2096943241276536 , 0.2096943241276536 , 0.0062247920880195 , 0.0000326589757115 );
            ir->AddPentPoint( 420 , 0.0364941472831919 , 0.2096943241276536 , 0.2096943241276536 , 0.0062247920880195 , 0.0000326589757115 );
            ir->AddPentPoint( 421 , 0.0062247920880195 , 0.0364941472831919 , 0.5378924123734813 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 422 , 0.2096943241276536 , 0.0364941472831919 , 0.5378924123734813 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 423 , 0.0364941472831919 , 0.0062247920880195 , 0.5378924123734813 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 424 , 0.2096943241276536 , 0.0062247920880195 , 0.5378924123734813 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 425 , 0.0364941472831919 , 0.2096943241276536 , 0.5378924123734813 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 426 , 0.0062247920880195 , 0.2096943241276536 , 0.5378924123734813 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 427 , 0.0062247920880195 , 0.5378924123734813 , 0.0364941472831919 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 428 , 0.2096943241276536 , 0.5378924123734813 , 0.0364941472831919 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 429 , 0.5378924123734813 , 0.0062247920880195 , 0.0364941472831919 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 430 , 0.2096943241276536 , 0.0062247920880195 , 0.0364941472831919 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 431 , 0.5378924123734813 , 0.2096943241276536 , 0.0364941472831919 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 432 , 0.0062247920880195 , 0.2096943241276536 , 0.0364941472831919 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 433 , 0.0364941472831919 , 0.5378924123734813 , 0.0062247920880195 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 434 , 0.2096943241276536 , 0.5378924123734813 , 0.0062247920880195 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 435 , 0.5378924123734813 , 0.0364941472831919 , 0.0062247920880195 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 436 , 0.2096943241276536 , 0.0364941472831919 , 0.0062247920880195 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 437 , 0.5378924123734813 , 0.2096943241276536 , 0.0062247920880195 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 438 , 0.0364941472831919 , 0.2096943241276536 , 0.0062247920880195 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 439 , 0.0364941472831919 , 0.5378924123734813 , 0.2096943241276536 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 440 , 0.0062247920880195 , 0.5378924123734813 , 0.2096943241276536 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 441 , 0.5378924123734813 , 0.0364941472831919 , 0.2096943241276536 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 442 , 0.0062247920880195 , 0.0364941472831919 , 0.2096943241276536 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 443 , 0.5378924123734813 , 0.0062247920880195 , 0.2096943241276536 , 0.2096943241276536 , 0.0000326589757115 );
            ir->AddPentPoint( 444 , 0.0364941472831919 , 0.0062247920880195 , 0.2096943241276536 , 0.2096943241276536 , 0.0000326589757115 );
            
            return ir;
            
            
        case 13:  // 555 points - degree 13 --
            PentatopeIntRules[13] = ir = new IntegrationRule(555);
            ir->AddPentPoint( 0 , 0.0107574287908706 , 0.0107574287908706 , 0.0107574287908706 , 0.9569702848365175 , 0.0000010772832712 );
            ir->AddPentPoint( 1 , 0.0107574287908706 , 0.0107574287908706 , 0.9569702848365175 , 0.0107574287908706 , 0.0000010772832712 );
            ir->AddPentPoint( 2 , 0.0107574287908706 , 0.9569702848365175 , 0.0107574287908706 , 0.0107574287908706 , 0.0000010772832712 );
            ir->AddPentPoint( 3 , 0.9569702848365175 , 0.0107574287908706 , 0.0107574287908706 , 0.0107574287908706 , 0.0000010772832712 );
            ir->AddPentPoint( 4 , 0.0107574287908706 , 0.0107574287908706 , 0.0107574287908706 , 0.0107574287908706 , 0.0000010772832712 );
            ir->AddPentPoint( 5 , 0.0788267222478506 , 0.0788267222478506 , 0.0788267222478506 , 0.6846931110085976 , 0.0001005898411824 );
            ir->AddPentPoint( 6 , 0.0788267222478506 , 0.0788267222478506 , 0.6846931110085976 , 0.0788267222478506 , 0.0001005898411824 );
            ir->AddPentPoint( 7 , 0.0788267222478506 , 0.6846931110085976 , 0.0788267222478506 , 0.0788267222478506 , 0.0001005898411824 );
            ir->AddPentPoint( 8 , 0.6846931110085976 , 0.0788267222478506 , 0.0788267222478506 , 0.0788267222478506 , 0.0001005898411824 );
            ir->AddPentPoint( 9 , 0.0788267222478506 , 0.0788267222478506 , 0.0788267222478506 , 0.0788267222478506 , 0.0001005898411824 );
            ir->AddPentPoint( 10 , 0.1453928459440898 , 0.1453928459440898 , 0.1453928459440898 , 0.4184286162236409 , 0.0003936217726832 );
            ir->AddPentPoint( 11 , 0.1453928459440898 , 0.1453928459440898 , 0.4184286162236409 , 0.1453928459440898 , 0.0003936217726832 );
            ir->AddPentPoint( 12 , 0.1453928459440898 , 0.4184286162236409 , 0.1453928459440898 , 0.1453928459440898 , 0.0003936217726832 );
            ir->AddPentPoint( 13 , 0.4184286162236409 , 0.1453928459440898 , 0.1453928459440898 , 0.1453928459440898 , 0.0003936217726832 );
            ir->AddPentPoint( 14 , 0.1453928459440898 , 0.1453928459440898 , 0.1453928459440898 , 0.1453928459440898 , 0.0003936217726832 );
            ir->AddPentPoint( 15 , 0.2266042893834491 , 0.2266042893834491 , 0.2266042893834491 , 0.0935828424662035 , 0.0004333886530915 );
            ir->AddPentPoint( 16 , 0.2266042893834491 , 0.2266042893834491 , 0.0935828424662035 , 0.2266042893834491 , 0.0004333886530915 );
            ir->AddPentPoint( 17 , 0.2266042893834491 , 0.0935828424662035 , 0.2266042893834491 , 0.2266042893834491 , 0.0004333886530915 );
            ir->AddPentPoint( 18 , 0.0935828424662035 , 0.2266042893834491 , 0.2266042893834491 , 0.2266042893834491 , 0.0004333886530915 );
            ir->AddPentPoint( 19 , 0.2266042893834491 , 0.2266042893834491 , 0.2266042893834491 , 0.2266042893834491 , 0.0004333886530915 );
            ir->AddPentPoint( 20 , 0.2477780735524920 , 0.2477780735524920 , 0.2477780735524920 , 0.0088877057900320 , 0.0001454345900369 );
            ir->AddPentPoint( 21 , 0.2477780735524920 , 0.2477780735524920 , 0.0088877057900320 , 0.2477780735524920 , 0.0001454345900369 );
            ir->AddPentPoint( 22 , 0.2477780735524920 , 0.0088877057900320 , 0.2477780735524920 , 0.2477780735524920 , 0.0001454345900369 );
            ir->AddPentPoint( 23 , 0.0088877057900320 , 0.2477780735524920 , 0.2477780735524920 , 0.2477780735524920 , 0.0001454345900369 );
            ir->AddPentPoint( 24 , 0.2477780735524920 , 0.2477780735524920 , 0.2477780735524920 , 0.2477780735524920 , 0.0001454345900369 );
            ir->AddPentPoint( 25 , 0.0430542476412466 , 0.0430542476412466 , 0.4354186285381301 , 0.4354186285381301 , 0.0000864578511064 );
            ir->AddPentPoint( 26 , 0.0430542476412466 , 0.4354186285381301 , 0.0430542476412466 , 0.4354186285381301 , 0.0000864578511064 );
            ir->AddPentPoint( 27 , 0.4354186285381301 , 0.0430542476412466 , 0.0430542476412466 , 0.4354186285381301 , 0.0000864578511064 );
            ir->AddPentPoint( 28 , 0.0430542476412466 , 0.0430542476412466 , 0.0430542476412466 , 0.4354186285381301 , 0.0000864578511064 );
            ir->AddPentPoint( 29 , 0.0430542476412466 , 0.4354186285381301 , 0.4354186285381301 , 0.0430542476412466 , 0.0000864578511064 );
            ir->AddPentPoint( 30 , 0.4354186285381301 , 0.0430542476412466 , 0.4354186285381301 , 0.0430542476412466 , 0.0000864578511064 );
            ir->AddPentPoint( 31 , 0.0430542476412466 , 0.0430542476412466 , 0.4354186285381301 , 0.0430542476412466 , 0.0000864578511064 );
            ir->AddPentPoint( 32 , 0.4354186285381301 , 0.4354186285381301 , 0.0430542476412466 , 0.0430542476412466 , 0.0000864578511064 );
            ir->AddPentPoint( 33 , 0.0430542476412466 , 0.4354186285381301 , 0.0430542476412466 , 0.0430542476412466 , 0.0000864578511064 );
            ir->AddPentPoint( 34 , 0.4354186285381301 , 0.0430542476412466 , 0.0430542476412466 , 0.0430542476412466 , 0.0000864578511064 );
            ir->AddPentPoint( 35 , 0.1088288352150819 , 0.1088288352150819 , 0.3367567471773771 , 0.3367567471773771 , 0.0002598226299977 );
            ir->AddPentPoint( 36 , 0.1088288352150819 , 0.3367567471773771 , 0.1088288352150819 , 0.3367567471773771 , 0.0002598226299977 );
            ir->AddPentPoint( 37 , 0.3367567471773771 , 0.1088288352150819 , 0.1088288352150819 , 0.3367567471773771 , 0.0002598226299977 );
            ir->AddPentPoint( 38 , 0.1088288352150819 , 0.1088288352150819 , 0.1088288352150819 , 0.3367567471773771 , 0.0002598226299977 );
            ir->AddPentPoint( 39 , 0.1088288352150819 , 0.3367567471773771 , 0.3367567471773771 , 0.1088288352150819 , 0.0002598226299977 );
            ir->AddPentPoint( 40 , 0.3367567471773771 , 0.1088288352150819 , 0.3367567471773771 , 0.1088288352150819 , 0.0002598226299977 );
            ir->AddPentPoint( 41 , 0.1088288352150819 , 0.1088288352150819 , 0.3367567471773771 , 0.1088288352150819 , 0.0002598226299977 );
            ir->AddPentPoint( 42 , 0.3367567471773771 , 0.3367567471773771 , 0.1088288352150819 , 0.1088288352150819 , 0.0002598226299977 );
            ir->AddPentPoint( 43 , 0.1088288352150819 , 0.3367567471773771 , 0.1088288352150819 , 0.1088288352150819 , 0.0002598226299977 );
            ir->AddPentPoint( 44 , 0.3367567471773771 , 0.1088288352150819 , 0.1088288352150819 , 0.1088288352150819 , 0.0002598226299977 );
            ir->AddPentPoint( 45 , 0.2836128758790833 , 0.2836128758790833 , 0.0745806861813750 , 0.0745806861813750 , 0.0003071106754142 );
            ir->AddPentPoint( 46 , 0.2836128758790833 , 0.0745806861813750 , 0.2836128758790833 , 0.0745806861813750 , 0.0003071106754142 );
            ir->AddPentPoint( 47 , 0.0745806861813750 , 0.2836128758790833 , 0.2836128758790833 , 0.0745806861813750 , 0.0003071106754142 );
            ir->AddPentPoint( 48 , 0.2836128758790833 , 0.2836128758790833 , 0.2836128758790833 , 0.0745806861813750 , 0.0003071106754142 );
            ir->AddPentPoint( 49 , 0.2836128758790833 , 0.0745806861813750 , 0.0745806861813750 , 0.2836128758790833 , 0.0003071106754142 );
            ir->AddPentPoint( 50 , 0.0745806861813750 , 0.2836128758790833 , 0.0745806861813750 , 0.2836128758790833 , 0.0003071106754142 );
            ir->AddPentPoint( 51 , 0.2836128758790833 , 0.2836128758790833 , 0.0745806861813750 , 0.2836128758790833 , 0.0003071106754142 );
            ir->AddPentPoint( 52 , 0.0745806861813750 , 0.0745806861813750 , 0.2836128758790833 , 0.2836128758790833 , 0.0003071106754142 );
            ir->AddPentPoint( 53 , 0.2836128758790833 , 0.0745806861813750 , 0.2836128758790833 , 0.2836128758790833 , 0.0003071106754142 );
            ir->AddPentPoint( 54 , 0.0745806861813750 , 0.2836128758790833 , 0.2836128758790833 , 0.2836128758790833 , 0.0003071106754142 );
            ir->AddPentPoint( 55 , 0.0205766324183731 , 0.0205766324183731 , 0.1097851969171234 , 0.8284849058277572 , 0.0000177303737100 );
            ir->AddPentPoint( 56 , 0.0205766324183731 , 0.1097851969171234 , 0.0205766324183731 , 0.8284849058277572 , 0.0000177303737100 );
            ir->AddPentPoint( 57 , 0.1097851969171234 , 0.0205766324183731 , 0.0205766324183731 , 0.8284849058277572 , 0.0000177303737100 );
            ir->AddPentPoint( 58 , 0.0205766324183731 , 0.0205766324183731 , 0.0205766324183731 , 0.8284849058277572 , 0.0000177303737100 );
            ir->AddPentPoint( 59 , 0.0205766324183731 , 0.0205766324183731 , 0.8284849058277572 , 0.1097851969171234 , 0.0000177303737100 );
            ir->AddPentPoint( 60 , 0.0205766324183731 , 0.8284849058277572 , 0.0205766324183731 , 0.1097851969171234 , 0.0000177303737100 );
            ir->AddPentPoint( 61 , 0.8284849058277572 , 0.0205766324183731 , 0.0205766324183731 , 0.1097851969171234 , 0.0000177303737100 );
            ir->AddPentPoint( 62 , 0.0205766324183731 , 0.0205766324183731 , 0.0205766324183731 , 0.1097851969171234 , 0.0000177303737100 );
            ir->AddPentPoint( 63 , 0.0205766324183731 , 0.1097851969171234 , 0.8284849058277572 , 0.0205766324183731 , 0.0000177303737100 );
            ir->AddPentPoint( 64 , 0.1097851969171234 , 0.0205766324183731 , 0.8284849058277572 , 0.0205766324183731 , 0.0000177303737100 );
            ir->AddPentPoint( 65 , 0.0205766324183731 , 0.0205766324183731 , 0.8284849058277572 , 0.0205766324183731 , 0.0000177303737100 );
            ir->AddPentPoint( 66 , 0.0205766324183731 , 0.8284849058277572 , 0.1097851969171234 , 0.0205766324183731 , 0.0000177303737100 );
            ir->AddPentPoint( 67 , 0.8284849058277572 , 0.0205766324183731 , 0.1097851969171234 , 0.0205766324183731 , 0.0000177303737100 );
            ir->AddPentPoint( 68 , 0.0205766324183731 , 0.0205766324183731 , 0.1097851969171234 , 0.0205766324183731 , 0.0000177303737100 );
            ir->AddPentPoint( 69 , 0.1097851969171234 , 0.8284849058277572 , 0.0205766324183731 , 0.0205766324183731 , 0.0000177303737100 );
            ir->AddPentPoint( 70 , 0.0205766324183731 , 0.8284849058277572 , 0.0205766324183731 , 0.0205766324183731 , 0.0000177303737100 );
            ir->AddPentPoint( 71 , 0.8284849058277572 , 0.1097851969171234 , 0.0205766324183731 , 0.0205766324183731 , 0.0000177303737100 );
            ir->AddPentPoint( 72 , 0.0205766324183731 , 0.1097851969171234 , 0.0205766324183731 , 0.0205766324183731 , 0.0000177303737100 );
            ir->AddPentPoint( 73 , 0.8284849058277572 , 0.0205766324183731 , 0.0205766324183731 , 0.0205766324183731 , 0.0000177303737100 );
            ir->AddPentPoint( 74 , 0.1097851969171234 , 0.0205766324183731 , 0.0205766324183731 , 0.0205766324183731 , 0.0000177303737100 );
            ir->AddPentPoint( 75 , 0.0222178991799967 , 0.0222178991799967 , 0.0492959792306157 , 0.8840503232293943 , 0.0000034588412869 );
            ir->AddPentPoint( 76 , 0.0222178991799967 , 0.0492959792306157 , 0.0222178991799967 , 0.8840503232293943 , 0.0000034588412869 );
            ir->AddPentPoint( 77 , 0.0492959792306157 , 0.0222178991799967 , 0.0222178991799967 , 0.8840503232293943 , 0.0000034588412869 );
            ir->AddPentPoint( 78 , 0.0222178991799967 , 0.0222178991799967 , 0.0222178991799967 , 0.8840503232293943 , 0.0000034588412869 );
            ir->AddPentPoint( 79 , 0.0222178991799967 , 0.0222178991799967 , 0.8840503232293943 , 0.0492959792306157 , 0.0000034588412869 );
            ir->AddPentPoint( 80 , 0.0222178991799967 , 0.8840503232293943 , 0.0222178991799967 , 0.0492959792306157 , 0.0000034588412869 );
            ir->AddPentPoint( 81 , 0.8840503232293943 , 0.0222178991799967 , 0.0222178991799967 , 0.0492959792306157 , 0.0000034588412869 );
            ir->AddPentPoint( 82 , 0.0222178991799967 , 0.0222178991799967 , 0.0222178991799967 , 0.0492959792306157 , 0.0000034588412869 );
            ir->AddPentPoint( 83 , 0.0222178991799967 , 0.0492959792306157 , 0.8840503232293943 , 0.0222178991799967 , 0.0000034588412869 );
            ir->AddPentPoint( 84 , 0.0492959792306157 , 0.0222178991799967 , 0.8840503232293943 , 0.0222178991799967 , 0.0000034588412869 );
            ir->AddPentPoint( 85 , 0.0222178991799967 , 0.0222178991799967 , 0.8840503232293943 , 0.0222178991799967 , 0.0000034588412869 );
            ir->AddPentPoint( 86 , 0.0222178991799967 , 0.8840503232293943 , 0.0492959792306157 , 0.0222178991799967 , 0.0000034588412869 );
            ir->AddPentPoint( 87 , 0.8840503232293943 , 0.0222178991799967 , 0.0492959792306157 , 0.0222178991799967 , 0.0000034588412869 );
            ir->AddPentPoint( 88 , 0.0222178991799967 , 0.0222178991799967 , 0.0492959792306157 , 0.0222178991799967 , 0.0000034588412869 );
            ir->AddPentPoint( 89 , 0.0492959792306157 , 0.8840503232293943 , 0.0222178991799967 , 0.0222178991799967 , 0.0000034588412869 );
            ir->AddPentPoint( 90 , 0.0222178991799967 , 0.8840503232293943 , 0.0222178991799967 , 0.0222178991799967 , 0.0000034588412869 );
            ir->AddPentPoint( 91 , 0.8840503232293943 , 0.0492959792306157 , 0.0222178991799967 , 0.0222178991799967 , 0.0000034588412869 );
            ir->AddPentPoint( 92 , 0.0222178991799967 , 0.0492959792306157 , 0.0222178991799967 , 0.0222178991799967 , 0.0000034588412869 );
            ir->AddPentPoint( 93 , 0.8840503232293943 , 0.0222178991799967 , 0.0222178991799967 , 0.0222178991799967 , 0.0000034588412869 );
            ir->AddPentPoint( 94 , 0.0492959792306157 , 0.0222178991799967 , 0.0222178991799967 , 0.0222178991799967 , 0.0000034588412869 );
            ir->AddPentPoint( 95 , 0.0230693073596855 , 0.0230693073596855 , 0.2532573583869335 , 0.6775347195340101 , 0.0000336483272322 );
            ir->AddPentPoint( 96 , 0.0230693073596855 , 0.2532573583869335 , 0.0230693073596855 , 0.6775347195340101 , 0.0000336483272322 );
            ir->AddPentPoint( 97 , 0.2532573583869335 , 0.0230693073596855 , 0.0230693073596855 , 0.6775347195340101 , 0.0000336483272322 );
            ir->AddPentPoint( 98 , 0.0230693073596855 , 0.0230693073596855 , 0.0230693073596855 , 0.6775347195340101 , 0.0000336483272322 );
            ir->AddPentPoint( 99 , 0.0230693073596855 , 0.0230693073596855 , 0.6775347195340101 , 0.2532573583869335 , 0.0000336483272322 );
            ir->AddPentPoint( 100 , 0.0230693073596855 , 0.6775347195340101 , 0.0230693073596855 , 0.2532573583869335 , 0.0000336483272322 );
            ir->AddPentPoint( 101 , 0.6775347195340101 , 0.0230693073596855 , 0.0230693073596855 , 0.2532573583869335 , 0.0000336483272322 );
            ir->AddPentPoint( 102 , 0.0230693073596855 , 0.0230693073596855 , 0.0230693073596855 , 0.2532573583869335 , 0.0000336483272322 );
            ir->AddPentPoint( 103 , 0.0230693073596855 , 0.2532573583869335 , 0.6775347195340101 , 0.0230693073596855 , 0.0000336483272322 );
            ir->AddPentPoint( 104 , 0.2532573583869335 , 0.0230693073596855 , 0.6775347195340101 , 0.0230693073596855 , 0.0000336483272322 );
            ir->AddPentPoint( 105 , 0.0230693073596855 , 0.0230693073596855 , 0.6775347195340101 , 0.0230693073596855 , 0.0000336483272322 );
            ir->AddPentPoint( 106 , 0.0230693073596855 , 0.6775347195340101 , 0.2532573583869335 , 0.0230693073596855 , 0.0000336483272322 );
            ir->AddPentPoint( 107 , 0.6775347195340101 , 0.0230693073596855 , 0.2532573583869335 , 0.0230693073596855 , 0.0000336483272322 );
            ir->AddPentPoint( 108 , 0.0230693073596855 , 0.0230693073596855 , 0.2532573583869335 , 0.0230693073596855 , 0.0000336483272322 );
            ir->AddPentPoint( 109 , 0.2532573583869335 , 0.6775347195340101 , 0.0230693073596855 , 0.0230693073596855 , 0.0000336483272322 );
            ir->AddPentPoint( 110 , 0.0230693073596855 , 0.6775347195340101 , 0.0230693073596855 , 0.0230693073596855 , 0.0000336483272322 );
            ir->AddPentPoint( 111 , 0.6775347195340101 , 0.2532573583869335 , 0.0230693073596855 , 0.0230693073596855 , 0.0000336483272322 );
            ir->AddPentPoint( 112 , 0.0230693073596855 , 0.2532573583869335 , 0.0230693073596855 , 0.0230693073596855 , 0.0000336483272322 );
            ir->AddPentPoint( 113 , 0.6775347195340101 , 0.0230693073596855 , 0.0230693073596855 , 0.0230693073596855 , 0.0000336483272322 );
            ir->AddPentPoint( 114 , 0.2532573583869335 , 0.0230693073596855 , 0.0230693073596855 , 0.0230693073596855 , 0.0000336483272322 );
            ir->AddPentPoint( 115 , 0.0722111547180188 , 0.0722111547180188 , 0.2162388270579805 , 0.5671277087879630 , 0.0001806112584012 );
            ir->AddPentPoint( 116 , 0.0722111547180188 , 0.2162388270579805 , 0.0722111547180188 , 0.5671277087879630 , 0.0001806112584012 );
            ir->AddPentPoint( 117 , 0.2162388270579805 , 0.0722111547180188 , 0.0722111547180188 , 0.5671277087879630 , 0.0001806112584012 );
            ir->AddPentPoint( 118 , 0.0722111547180188 , 0.0722111547180188 , 0.0722111547180188 , 0.5671277087879630 , 0.0001806112584012 );
            ir->AddPentPoint( 119 , 0.0722111547180188 , 0.0722111547180188 , 0.5671277087879630 , 0.2162388270579805 , 0.0001806112584012 );
            ir->AddPentPoint( 120 , 0.0722111547180188 , 0.5671277087879630 , 0.0722111547180188 , 0.2162388270579805 , 0.0001806112584012 );
            ir->AddPentPoint( 121 , 0.5671277087879630 , 0.0722111547180188 , 0.0722111547180188 , 0.2162388270579805 , 0.0001806112584012 );
            ir->AddPentPoint( 122 , 0.0722111547180188 , 0.0722111547180188 , 0.0722111547180188 , 0.2162388270579805 , 0.0001806112584012 );
            ir->AddPentPoint( 123 , 0.0722111547180188 , 0.2162388270579805 , 0.5671277087879630 , 0.0722111547180188 , 0.0001806112584012 );
            ir->AddPentPoint( 124 , 0.2162388270579805 , 0.0722111547180188 , 0.5671277087879630 , 0.0722111547180188 , 0.0001806112584012 );
            ir->AddPentPoint( 125 , 0.0722111547180188 , 0.0722111547180188 , 0.5671277087879630 , 0.0722111547180188 , 0.0001806112584012 );
            ir->AddPentPoint( 126 , 0.0722111547180188 , 0.5671277087879630 , 0.2162388270579805 , 0.0722111547180188 , 0.0001806112584012 );
            ir->AddPentPoint( 127 , 0.5671277087879630 , 0.0722111547180188 , 0.2162388270579805 , 0.0722111547180188 , 0.0001806112584012 );
            ir->AddPentPoint( 128 , 0.0722111547180188 , 0.0722111547180188 , 0.2162388270579805 , 0.0722111547180188 , 0.0001806112584012 );
            ir->AddPentPoint( 129 , 0.2162388270579805 , 0.5671277087879630 , 0.0722111547180188 , 0.0722111547180188 , 0.0001806112584012 );
            ir->AddPentPoint( 130 , 0.0722111547180188 , 0.5671277087879630 , 0.0722111547180188 , 0.0722111547180188 , 0.0001806112584012 );
            ir->AddPentPoint( 131 , 0.5671277087879630 , 0.2162388270579805 , 0.0722111547180188 , 0.0722111547180188 , 0.0001806112584012 );
            ir->AddPentPoint( 132 , 0.0722111547180188 , 0.2162388270579805 , 0.0722111547180188 , 0.0722111547180188 , 0.0001806112584012 );
            ir->AddPentPoint( 133 , 0.5671277087879630 , 0.0722111547180188 , 0.0722111547180188 , 0.0722111547180188 , 0.0001806112584012 );
            ir->AddPentPoint( 134 , 0.2162388270579805 , 0.0722111547180188 , 0.0722111547180188 , 0.0722111547180188 , 0.0001806112584012 );
            ir->AddPentPoint( 135 , 0.0888837602751494 , 0.0888837602751494 , 0.0009872476689666 , 0.7323614715055851 , 0.0000136103311491 );
            ir->AddPentPoint( 136 , 0.0888837602751494 , 0.0009872476689666 , 0.0888837602751494 , 0.7323614715055851 , 0.0000136103311491 );
            ir->AddPentPoint( 137 , 0.0009872476689666 , 0.0888837602751494 , 0.0888837602751494 , 0.7323614715055851 , 0.0000136103311491 );
            ir->AddPentPoint( 138 , 0.0888837602751494 , 0.0888837602751494 , 0.0888837602751494 , 0.7323614715055851 , 0.0000136103311491 );
            ir->AddPentPoint( 139 , 0.0888837602751494 , 0.0888837602751494 , 0.7323614715055851 , 0.0009872476689666 , 0.0000136103311491 );
            ir->AddPentPoint( 140 , 0.0888837602751494 , 0.7323614715055851 , 0.0888837602751494 , 0.0009872476689666 , 0.0000136103311491 );
            ir->AddPentPoint( 141 , 0.7323614715055851 , 0.0888837602751494 , 0.0888837602751494 , 0.0009872476689666 , 0.0000136103311491 );
            ir->AddPentPoint( 142 , 0.0888837602751494 , 0.0888837602751494 , 0.0888837602751494 , 0.0009872476689666 , 0.0000136103311491 );
            ir->AddPentPoint( 143 , 0.0888837602751494 , 0.0009872476689666 , 0.7323614715055851 , 0.0888837602751494 , 0.0000136103311491 );
            ir->AddPentPoint( 144 , 0.0009872476689666 , 0.0888837602751494 , 0.7323614715055851 , 0.0888837602751494 , 0.0000136103311491 );
            ir->AddPentPoint( 145 , 0.0888837602751494 , 0.0888837602751494 , 0.7323614715055851 , 0.0888837602751494 , 0.0000136103311491 );
            ir->AddPentPoint( 146 , 0.0888837602751494 , 0.7323614715055851 , 0.0009872476689666 , 0.0888837602751494 , 0.0000136103311491 );
            ir->AddPentPoint( 147 , 0.7323614715055851 , 0.0888837602751494 , 0.0009872476689666 , 0.0888837602751494 , 0.0000136103311491 );
            ir->AddPentPoint( 148 , 0.0888837602751494 , 0.0888837602751494 , 0.0009872476689666 , 0.0888837602751494 , 0.0000136103311491 );
            ir->AddPentPoint( 149 , 0.0009872476689666 , 0.7323614715055851 , 0.0888837602751494 , 0.0888837602751494 , 0.0000136103311491 );
            ir->AddPentPoint( 150 , 0.0888837602751494 , 0.7323614715055851 , 0.0888837602751494 , 0.0888837602751494 , 0.0000136103311491 );
            ir->AddPentPoint( 151 , 0.7323614715055851 , 0.0009872476689666 , 0.0888837602751494 , 0.0888837602751494 , 0.0000136103311491 );
            ir->AddPentPoint( 152 , 0.0888837602751494 , 0.0009872476689666 , 0.0888837602751494 , 0.0888837602751494 , 0.0000136103311491 );
            ir->AddPentPoint( 153 , 0.7323614715055851 , 0.0888837602751494 , 0.0888837602751494 , 0.0888837602751494 , 0.0000136103311491 );
            ir->AddPentPoint( 154 , 0.0009872476689666 , 0.0888837602751494 , 0.0888837602751494 , 0.0888837602751494 , 0.0000136103311491 );
            ir->AddPentPoint( 155 , 0.1350153362854520 , 0.1350153362854520 , 0.0104784532131386 , 0.5844755379305052 , 0.0000650682536048 );
            ir->AddPentPoint( 156 , 0.1350153362854520 , 0.0104784532131386 , 0.1350153362854520 , 0.5844755379305052 , 0.0000650682536048 );
            ir->AddPentPoint( 157 , 0.0104784532131386 , 0.1350153362854520 , 0.1350153362854520 , 0.5844755379305052 , 0.0000650682536048 );
            ir->AddPentPoint( 158 , 0.1350153362854520 , 0.1350153362854520 , 0.1350153362854520 , 0.5844755379305052 , 0.0000650682536048 );
            ir->AddPentPoint( 159 , 0.1350153362854520 , 0.1350153362854520 , 0.5844755379305052 , 0.0104784532131386 , 0.0000650682536048 );
            ir->AddPentPoint( 160 , 0.1350153362854520 , 0.5844755379305052 , 0.1350153362854520 , 0.0104784532131386 , 0.0000650682536048 );
            ir->AddPentPoint( 161 , 0.5844755379305052 , 0.1350153362854520 , 0.1350153362854520 , 0.0104784532131386 , 0.0000650682536048 );
            ir->AddPentPoint( 162 , 0.1350153362854520 , 0.1350153362854520 , 0.1350153362854520 , 0.0104784532131386 , 0.0000650682536048 );
            ir->AddPentPoint( 163 , 0.1350153362854520 , 0.0104784532131386 , 0.5844755379305052 , 0.1350153362854520 , 0.0000650682536048 );
            ir->AddPentPoint( 164 , 0.0104784532131386 , 0.1350153362854520 , 0.5844755379305052 , 0.1350153362854520 , 0.0000650682536048 );
            ir->AddPentPoint( 165 , 0.1350153362854520 , 0.1350153362854520 , 0.5844755379305052 , 0.1350153362854520 , 0.0000650682536048 );
            ir->AddPentPoint( 166 , 0.1350153362854520 , 0.5844755379305052 , 0.0104784532131386 , 0.1350153362854520 , 0.0000650682536048 );
            ir->AddPentPoint( 167 , 0.5844755379305052 , 0.1350153362854520 , 0.0104784532131386 , 0.1350153362854520 , 0.0000650682536048 );
            ir->AddPentPoint( 168 , 0.1350153362854520 , 0.1350153362854520 , 0.0104784532131386 , 0.1350153362854520 , 0.0000650682536048 );
            ir->AddPentPoint( 169 , 0.0104784532131386 , 0.5844755379305052 , 0.1350153362854520 , 0.1350153362854520 , 0.0000650682536048 );
            ir->AddPentPoint( 170 , 0.1350153362854520 , 0.5844755379305052 , 0.1350153362854520 , 0.1350153362854520 , 0.0000650682536048 );
            ir->AddPentPoint( 171 , 0.5844755379305052 , 0.0104784532131386 , 0.1350153362854520 , 0.1350153362854520 , 0.0000650682536048 );
            ir->AddPentPoint( 172 , 0.1350153362854520 , 0.0104784532131386 , 0.1350153362854520 , 0.1350153362854520 , 0.0000650682536048 );
            ir->AddPentPoint( 173 , 0.5844755379305052 , 0.1350153362854520 , 0.1350153362854520 , 0.1350153362854520 , 0.0000650682536048 );
            ir->AddPentPoint( 174 , 0.0104784532131386 , 0.1350153362854520 , 0.1350153362854520 , 0.1350153362854520 , 0.0000650682536048 );
            ir->AddPentPoint( 175 , 0.3284242055668443 , 0.3284242055668443 , 0.0000409931017634 , 0.0146863901977038 , 0.0000165186146256 );
            ir->AddPentPoint( 176 , 0.3284242055668443 , 0.0000409931017634 , 0.3284242055668443 , 0.0146863901977038 , 0.0000165186146256 );
            ir->AddPentPoint( 177 , 0.0000409931017634 , 0.3284242055668443 , 0.3284242055668443 , 0.0146863901977038 , 0.0000165186146256 );
            ir->AddPentPoint( 178 , 0.3284242055668443 , 0.3284242055668443 , 0.3284242055668443 , 0.0146863901977038 , 0.0000165186146256 );
            ir->AddPentPoint( 179 , 0.3284242055668443 , 0.3284242055668443 , 0.0146863901977038 , 0.0000409931017634 , 0.0000165186146256 );
            ir->AddPentPoint( 180 , 0.3284242055668443 , 0.0146863901977038 , 0.3284242055668443 , 0.0000409931017634 , 0.0000165186146256 );
            ir->AddPentPoint( 181 , 0.0146863901977038 , 0.3284242055668443 , 0.3284242055668443 , 0.0000409931017634 , 0.0000165186146256 );
            ir->AddPentPoint( 182 , 0.3284242055668443 , 0.3284242055668443 , 0.3284242055668443 , 0.0000409931017634 , 0.0000165186146256 );
            ir->AddPentPoint( 183 , 0.3284242055668443 , 0.0000409931017634 , 0.0146863901977038 , 0.3284242055668443 , 0.0000165186146256 );
            ir->AddPentPoint( 184 , 0.0000409931017634 , 0.3284242055668443 , 0.0146863901977038 , 0.3284242055668443 , 0.0000165186146256 );
            ir->AddPentPoint( 185 , 0.3284242055668443 , 0.3284242055668443 , 0.0146863901977038 , 0.3284242055668443 , 0.0000165186146256 );
            ir->AddPentPoint( 186 , 0.3284242055668443 , 0.0146863901977038 , 0.0000409931017634 , 0.3284242055668443 , 0.0000165186146256 );
            ir->AddPentPoint( 187 , 0.0146863901977038 , 0.3284242055668443 , 0.0000409931017634 , 0.3284242055668443 , 0.0000165186146256 );
            ir->AddPentPoint( 188 , 0.3284242055668443 , 0.3284242055668443 , 0.0000409931017634 , 0.3284242055668443 , 0.0000165186146256 );
            ir->AddPentPoint( 189 , 0.0000409931017634 , 0.0146863901977038 , 0.3284242055668443 , 0.3284242055668443 , 0.0000165186146256 );
            ir->AddPentPoint( 190 , 0.3284242055668443 , 0.0146863901977038 , 0.3284242055668443 , 0.3284242055668443 , 0.0000165186146256 );
            ir->AddPentPoint( 191 , 0.0146863901977038 , 0.0000409931017634 , 0.3284242055668443 , 0.3284242055668443 , 0.0000165186146256 );
            ir->AddPentPoint( 192 , 0.3284242055668443 , 0.0000409931017634 , 0.3284242055668443 , 0.3284242055668443 , 0.0000165186146256 );
            ir->AddPentPoint( 193 , 0.0146863901977038 , 0.3284242055668443 , 0.3284242055668443 , 0.3284242055668443 , 0.0000165186146256 );
            ir->AddPentPoint( 194 , 0.0000409931017634 , 0.3284242055668443 , 0.3284242055668443 , 0.3284242055668443 , 0.0000165186146256 );
            ir->AddPentPoint( 195 , 0.0042040139986841 , 0.4779107324155366 , 0.4779107324155366 , 0.0357705071715587 , 0.0000056154038485 );
            ir->AddPentPoint( 196 , 0.4779107324155366 , 0.0042040139986841 , 0.4779107324155366 , 0.0357705071715587 , 0.0000056154038485 );
            ir->AddPentPoint( 197 , 0.0042040139986841 , 0.0042040139986841 , 0.4779107324155366 , 0.0357705071715587 , 0.0000056154038485 );
            ir->AddPentPoint( 198 , 0.4779107324155366 , 0.4779107324155366 , 0.0042040139986841 , 0.0357705071715587 , 0.0000056154038485 );
            ir->AddPentPoint( 199 , 0.0042040139986841 , 0.4779107324155366 , 0.0042040139986841 , 0.0357705071715587 , 0.0000056154038485 );
            ir->AddPentPoint( 200 , 0.4779107324155366 , 0.0042040139986841 , 0.0042040139986841 , 0.0357705071715587 , 0.0000056154038485 );
            ir->AddPentPoint( 201 , 0.0042040139986841 , 0.4779107324155366 , 0.0357705071715587 , 0.4779107324155366 , 0.0000056154038485 );
            ir->AddPentPoint( 202 , 0.4779107324155366 , 0.0042040139986841 , 0.0357705071715587 , 0.4779107324155366 , 0.0000056154038485 );
            ir->AddPentPoint( 203 , 0.0042040139986841 , 0.0042040139986841 , 0.0357705071715587 , 0.4779107324155366 , 0.0000056154038485 );
            ir->AddPentPoint( 204 , 0.0042040139986841 , 0.0357705071715587 , 0.4779107324155366 , 0.4779107324155366 , 0.0000056154038485 );
            ir->AddPentPoint( 205 , 0.0357705071715587 , 0.0042040139986841 , 0.4779107324155366 , 0.4779107324155366 , 0.0000056154038485 );
            ir->AddPentPoint( 206 , 0.0042040139986841 , 0.0042040139986841 , 0.4779107324155366 , 0.4779107324155366 , 0.0000056154038485 );
            ir->AddPentPoint( 207 , 0.4779107324155366 , 0.0357705071715587 , 0.0042040139986841 , 0.4779107324155366 , 0.0000056154038485 );
            ir->AddPentPoint( 208 , 0.0042040139986841 , 0.0357705071715587 , 0.0042040139986841 , 0.4779107324155366 , 0.0000056154038485 );
            ir->AddPentPoint( 209 , 0.0357705071715587 , 0.4779107324155366 , 0.0042040139986841 , 0.4779107324155366 , 0.0000056154038485 );
            ir->AddPentPoint( 210 , 0.0042040139986841 , 0.4779107324155366 , 0.0042040139986841 , 0.4779107324155366 , 0.0000056154038485 );
            ir->AddPentPoint( 211 , 0.0357705071715587 , 0.0042040139986841 , 0.0042040139986841 , 0.4779107324155366 , 0.0000056154038485 );
            ir->AddPentPoint( 212 , 0.4779107324155366 , 0.0042040139986841 , 0.0042040139986841 , 0.4779107324155366 , 0.0000056154038485 );
            ir->AddPentPoint( 213 , 0.4779107324155366 , 0.4779107324155366 , 0.0357705071715587 , 0.0042040139986841 , 0.0000056154038485 );
            ir->AddPentPoint( 214 , 0.0042040139986841 , 0.4779107324155366 , 0.0357705071715587 , 0.0042040139986841 , 0.0000056154038485 );
            ir->AddPentPoint( 215 , 0.4779107324155366 , 0.0042040139986841 , 0.0357705071715587 , 0.0042040139986841 , 0.0000056154038485 );
            ir->AddPentPoint( 216 , 0.4779107324155366 , 0.0357705071715587 , 0.4779107324155366 , 0.0042040139986841 , 0.0000056154038485 );
            ir->AddPentPoint( 217 , 0.0042040139986841 , 0.0357705071715587 , 0.4779107324155366 , 0.0042040139986841 , 0.0000056154038485 );
            ir->AddPentPoint( 218 , 0.0357705071715587 , 0.4779107324155366 , 0.4779107324155366 , 0.0042040139986841 , 0.0000056154038485 );
            ir->AddPentPoint( 219 , 0.0042040139986841 , 0.4779107324155366 , 0.4779107324155366 , 0.0042040139986841 , 0.0000056154038485 );
            ir->AddPentPoint( 220 , 0.0357705071715587 , 0.0042040139986841 , 0.4779107324155366 , 0.0042040139986841 , 0.0000056154038485 );
            ir->AddPentPoint( 221 , 0.4779107324155366 , 0.0042040139986841 , 0.4779107324155366 , 0.0042040139986841 , 0.0000056154038485 );
            ir->AddPentPoint( 222 , 0.4779107324155366 , 0.0357705071715587 , 0.0042040139986841 , 0.0042040139986841 , 0.0000056154038485 );
            ir->AddPentPoint( 223 , 0.0357705071715587 , 0.4779107324155366 , 0.0042040139986841 , 0.0042040139986841 , 0.0000056154038485 );
            ir->AddPentPoint( 224 , 0.4779107324155366 , 0.4779107324155366 , 0.0042040139986841 , 0.0042040139986841 , 0.0000056154038485 );
            ir->AddPentPoint( 225 , 0.0209198170627071 , 0.1196871103129186 , 0.1196871103129186 , 0.7187861452487485 , 0.0000508812549516 );
            ir->AddPentPoint( 226 , 0.1196871103129186 , 0.0209198170627071 , 0.1196871103129186 , 0.7187861452487485 , 0.0000508812549516 );
            ir->AddPentPoint( 227 , 0.0209198170627071 , 0.0209198170627071 , 0.1196871103129186 , 0.7187861452487485 , 0.0000508812549516 );
            ir->AddPentPoint( 228 , 0.1196871103129186 , 0.1196871103129186 , 0.0209198170627071 , 0.7187861452487485 , 0.0000508812549516 );
            ir->AddPentPoint( 229 , 0.0209198170627071 , 0.1196871103129186 , 0.0209198170627071 , 0.7187861452487485 , 0.0000508812549516 );
            ir->AddPentPoint( 230 , 0.1196871103129186 , 0.0209198170627071 , 0.0209198170627071 , 0.7187861452487485 , 0.0000508812549516 );
            ir->AddPentPoint( 231 , 0.0209198170627071 , 0.1196871103129186 , 0.7187861452487485 , 0.1196871103129186 , 0.0000508812549516 );
            ir->AddPentPoint( 232 , 0.1196871103129186 , 0.0209198170627071 , 0.7187861452487485 , 0.1196871103129186 , 0.0000508812549516 );
            ir->AddPentPoint( 233 , 0.0209198170627071 , 0.0209198170627071 , 0.7187861452487485 , 0.1196871103129186 , 0.0000508812549516 );
            ir->AddPentPoint( 234 , 0.0209198170627071 , 0.7187861452487485 , 0.1196871103129186 , 0.1196871103129186 , 0.0000508812549516 );
            ir->AddPentPoint( 235 , 0.7187861452487485 , 0.0209198170627071 , 0.1196871103129186 , 0.1196871103129186 , 0.0000508812549516 );
            ir->AddPentPoint( 236 , 0.0209198170627071 , 0.0209198170627071 , 0.1196871103129186 , 0.1196871103129186 , 0.0000508812549516 );
            ir->AddPentPoint( 237 , 0.1196871103129186 , 0.7187861452487485 , 0.0209198170627071 , 0.1196871103129186 , 0.0000508812549516 );
            ir->AddPentPoint( 238 , 0.0209198170627071 , 0.7187861452487485 , 0.0209198170627071 , 0.1196871103129186 , 0.0000508812549516 );
            ir->AddPentPoint( 239 , 0.7187861452487485 , 0.1196871103129186 , 0.0209198170627071 , 0.1196871103129186 , 0.0000508812549516 );
            ir->AddPentPoint( 240 , 0.0209198170627071 , 0.1196871103129186 , 0.0209198170627071 , 0.1196871103129186 , 0.0000508812549516 );
            ir->AddPentPoint( 241 , 0.7187861452487485 , 0.0209198170627071 , 0.0209198170627071 , 0.1196871103129186 , 0.0000508812549516 );
            ir->AddPentPoint( 242 , 0.1196871103129186 , 0.0209198170627071 , 0.0209198170627071 , 0.1196871103129186 , 0.0000508812549516 );
            ir->AddPentPoint( 243 , 0.1196871103129186 , 0.1196871103129186 , 0.7187861452487485 , 0.0209198170627071 , 0.0000508812549516 );
            ir->AddPentPoint( 244 , 0.0209198170627071 , 0.1196871103129186 , 0.7187861452487485 , 0.0209198170627071 , 0.0000508812549516 );
            ir->AddPentPoint( 245 , 0.1196871103129186 , 0.0209198170627071 , 0.7187861452487485 , 0.0209198170627071 , 0.0000508812549516 );
            ir->AddPentPoint( 246 , 0.1196871103129186 , 0.7187861452487485 , 0.1196871103129186 , 0.0209198170627071 , 0.0000508812549516 );
            ir->AddPentPoint( 247 , 0.0209198170627071 , 0.7187861452487485 , 0.1196871103129186 , 0.0209198170627071 , 0.0000508812549516 );
            ir->AddPentPoint( 248 , 0.7187861452487485 , 0.1196871103129186 , 0.1196871103129186 , 0.0209198170627071 , 0.0000508812549516 );
            ir->AddPentPoint( 249 , 0.0209198170627071 , 0.1196871103129186 , 0.1196871103129186 , 0.0209198170627071 , 0.0000508812549516 );
            ir->AddPentPoint( 250 , 0.7187861452487485 , 0.0209198170627071 , 0.1196871103129186 , 0.0209198170627071 , 0.0000508812549516 );
            ir->AddPentPoint( 251 , 0.1196871103129186 , 0.0209198170627071 , 0.1196871103129186 , 0.0209198170627071 , 0.0000508812549516 );
            ir->AddPentPoint( 252 , 0.1196871103129186 , 0.7187861452487485 , 0.0209198170627071 , 0.0209198170627071 , 0.0000508812549516 );
            ir->AddPentPoint( 253 , 0.7187861452487485 , 0.1196871103129186 , 0.0209198170627071 , 0.0209198170627071 , 0.0000508812549516 );
            ir->AddPentPoint( 254 , 0.1196871103129186 , 0.1196871103129186 , 0.0209198170627071 , 0.0209198170627071 , 0.0000508812549516 );
            ir->AddPentPoint( 255 , 0.0297565060157150 , 0.4027982497493572 , 0.4027982497493572 , 0.1348904884698556 , 0.0000890696152729 );
            ir->AddPentPoint( 256 , 0.4027982497493572 , 0.0297565060157150 , 0.4027982497493572 , 0.1348904884698556 , 0.0000890696152729 );
            ir->AddPentPoint( 257 , 0.0297565060157150 , 0.0297565060157150 , 0.4027982497493572 , 0.1348904884698556 , 0.0000890696152729 );
            ir->AddPentPoint( 258 , 0.4027982497493572 , 0.4027982497493572 , 0.0297565060157150 , 0.1348904884698556 , 0.0000890696152729 );
            ir->AddPentPoint( 259 , 0.0297565060157150 , 0.4027982497493572 , 0.0297565060157150 , 0.1348904884698556 , 0.0000890696152729 );
            ir->AddPentPoint( 260 , 0.4027982497493572 , 0.0297565060157150 , 0.0297565060157150 , 0.1348904884698556 , 0.0000890696152729 );
            ir->AddPentPoint( 261 , 0.0297565060157150 , 0.4027982497493572 , 0.1348904884698556 , 0.4027982497493572 , 0.0000890696152729 );
            ir->AddPentPoint( 262 , 0.4027982497493572 , 0.0297565060157150 , 0.1348904884698556 , 0.4027982497493572 , 0.0000890696152729 );
            ir->AddPentPoint( 263 , 0.0297565060157150 , 0.0297565060157150 , 0.1348904884698556 , 0.4027982497493572 , 0.0000890696152729 );
            ir->AddPentPoint( 264 , 0.0297565060157150 , 0.1348904884698556 , 0.4027982497493572 , 0.4027982497493572 , 0.0000890696152729 );
            ir->AddPentPoint( 265 , 0.1348904884698556 , 0.0297565060157150 , 0.4027982497493572 , 0.4027982497493572 , 0.0000890696152729 );
            ir->AddPentPoint( 266 , 0.0297565060157150 , 0.0297565060157150 , 0.4027982497493572 , 0.4027982497493572 , 0.0000890696152729 );
            ir->AddPentPoint( 267 , 0.4027982497493572 , 0.1348904884698556 , 0.0297565060157150 , 0.4027982497493572 , 0.0000890696152729 );
            ir->AddPentPoint( 268 , 0.0297565060157150 , 0.1348904884698556 , 0.0297565060157150 , 0.4027982497493572 , 0.0000890696152729 );
            ir->AddPentPoint( 269 , 0.1348904884698556 , 0.4027982497493572 , 0.0297565060157150 , 0.4027982497493572 , 0.0000890696152729 );
            ir->AddPentPoint( 270 , 0.0297565060157150 , 0.4027982497493572 , 0.0297565060157150 , 0.4027982497493572 , 0.0000890696152729 );
            ir->AddPentPoint( 271 , 0.1348904884698556 , 0.0297565060157150 , 0.0297565060157150 , 0.4027982497493572 , 0.0000890696152729 );
            ir->AddPentPoint( 272 , 0.4027982497493572 , 0.0297565060157150 , 0.0297565060157150 , 0.4027982497493572 , 0.0000890696152729 );
            ir->AddPentPoint( 273 , 0.4027982497493572 , 0.4027982497493572 , 0.1348904884698556 , 0.0297565060157150 , 0.0000890696152729 );
            ir->AddPentPoint( 274 , 0.0297565060157150 , 0.4027982497493572 , 0.1348904884698556 , 0.0297565060157150 , 0.0000890696152729 );
            ir->AddPentPoint( 275 , 0.4027982497493572 , 0.0297565060157150 , 0.1348904884698556 , 0.0297565060157150 , 0.0000890696152729 );
            ir->AddPentPoint( 276 , 0.4027982497493572 , 0.1348904884698556 , 0.4027982497493572 , 0.0297565060157150 , 0.0000890696152729 );
            ir->AddPentPoint( 277 , 0.0297565060157150 , 0.1348904884698556 , 0.4027982497493572 , 0.0297565060157150 , 0.0000890696152729 );
            ir->AddPentPoint( 278 , 0.1348904884698556 , 0.4027982497493572 , 0.4027982497493572 , 0.0297565060157150 , 0.0000890696152729 );
            ir->AddPentPoint( 279 , 0.0297565060157150 , 0.4027982497493572 , 0.4027982497493572 , 0.0297565060157150 , 0.0000890696152729 );
            ir->AddPentPoint( 280 , 0.1348904884698556 , 0.0297565060157150 , 0.4027982497493572 , 0.0297565060157150 , 0.0000890696152729 );
            ir->AddPentPoint( 281 , 0.4027982497493572 , 0.0297565060157150 , 0.4027982497493572 , 0.0297565060157150 , 0.0000890696152729 );
            ir->AddPentPoint( 282 , 0.4027982497493572 , 0.1348904884698556 , 0.0297565060157150 , 0.0297565060157150 , 0.0000890696152729 );
            ir->AddPentPoint( 283 , 0.1348904884698556 , 0.4027982497493572 , 0.0297565060157150 , 0.0297565060157150 , 0.0000890696152729 );
            ir->AddPentPoint( 284 , 0.4027982497493572 , 0.4027982497493572 , 0.0297565060157150 , 0.0297565060157150 , 0.0000890696152729 );
            ir->AddPentPoint( 285 , 0.0325068593353589 , 0.2112903026922174 , 0.2112903026922174 , 0.5124056759448475 , 0.0001329751563098 );
            ir->AddPentPoint( 286 , 0.2112903026922174 , 0.0325068593353589 , 0.2112903026922174 , 0.5124056759448475 , 0.0001329751563098 );
            ir->AddPentPoint( 287 , 0.0325068593353589 , 0.0325068593353589 , 0.2112903026922174 , 0.5124056759448475 , 0.0001329751563098 );
            ir->AddPentPoint( 288 , 0.2112903026922174 , 0.2112903026922174 , 0.0325068593353589 , 0.5124056759448475 , 0.0001329751563098 );
            ir->AddPentPoint( 289 , 0.0325068593353589 , 0.2112903026922174 , 0.0325068593353589 , 0.5124056759448475 , 0.0001329751563098 );
            ir->AddPentPoint( 290 , 0.2112903026922174 , 0.0325068593353589 , 0.0325068593353589 , 0.5124056759448475 , 0.0001329751563098 );
            ir->AddPentPoint( 291 , 0.0325068593353589 , 0.2112903026922174 , 0.5124056759448475 , 0.2112903026922174 , 0.0001329751563098 );
            ir->AddPentPoint( 292 , 0.2112903026922174 , 0.0325068593353589 , 0.5124056759448475 , 0.2112903026922174 , 0.0001329751563098 );
            ir->AddPentPoint( 293 , 0.0325068593353589 , 0.0325068593353589 , 0.5124056759448475 , 0.2112903026922174 , 0.0001329751563098 );
            ir->AddPentPoint( 294 , 0.0325068593353589 , 0.5124056759448475 , 0.2112903026922174 , 0.2112903026922174 , 0.0001329751563098 );
            ir->AddPentPoint( 295 , 0.5124056759448475 , 0.0325068593353589 , 0.2112903026922174 , 0.2112903026922174 , 0.0001329751563098 );
            ir->AddPentPoint( 296 , 0.0325068593353589 , 0.0325068593353589 , 0.2112903026922174 , 0.2112903026922174 , 0.0001329751563098 );
            ir->AddPentPoint( 297 , 0.2112903026922174 , 0.5124056759448475 , 0.0325068593353589 , 0.2112903026922174 , 0.0001329751563098 );
            ir->AddPentPoint( 298 , 0.0325068593353589 , 0.5124056759448475 , 0.0325068593353589 , 0.2112903026922174 , 0.0001329751563098 );
            ir->AddPentPoint( 299 , 0.5124056759448475 , 0.2112903026922174 , 0.0325068593353589 , 0.2112903026922174 , 0.0001329751563098 );
            ir->AddPentPoint( 300 , 0.0325068593353589 , 0.2112903026922174 , 0.0325068593353589 , 0.2112903026922174 , 0.0001329751563098 );
            ir->AddPentPoint( 301 , 0.5124056759448475 , 0.0325068593353589 , 0.0325068593353589 , 0.2112903026922174 , 0.0001329751563098 );
            ir->AddPentPoint( 302 , 0.2112903026922174 , 0.0325068593353589 , 0.0325068593353589 , 0.2112903026922174 , 0.0001329751563098 );
            ir->AddPentPoint( 303 , 0.2112903026922174 , 0.2112903026922174 , 0.5124056759448475 , 0.0325068593353589 , 0.0001329751563098 );
            ir->AddPentPoint( 304 , 0.0325068593353589 , 0.2112903026922174 , 0.5124056759448475 , 0.0325068593353589 , 0.0001329751563098 );
            ir->AddPentPoint( 305 , 0.2112903026922174 , 0.0325068593353589 , 0.5124056759448475 , 0.0325068593353589 , 0.0001329751563098 );
            ir->AddPentPoint( 306 , 0.2112903026922174 , 0.5124056759448475 , 0.2112903026922174 , 0.0325068593353589 , 0.0001329751563098 );
            ir->AddPentPoint( 307 , 0.0325068593353589 , 0.5124056759448475 , 0.2112903026922174 , 0.0325068593353589 , 0.0001329751563098 );
            ir->AddPentPoint( 308 , 0.5124056759448475 , 0.2112903026922174 , 0.2112903026922174 , 0.0325068593353589 , 0.0001329751563098 );
            ir->AddPentPoint( 309 , 0.0325068593353589 , 0.2112903026922174 , 0.2112903026922174 , 0.0325068593353589 , 0.0001329751563098 );
            ir->AddPentPoint( 310 , 0.5124056759448475 , 0.0325068593353589 , 0.2112903026922174 , 0.0325068593353589 , 0.0001329751563098 );
            ir->AddPentPoint( 311 , 0.2112903026922174 , 0.0325068593353589 , 0.2112903026922174 , 0.0325068593353589 , 0.0001329751563098 );
            ir->AddPentPoint( 312 , 0.2112903026922174 , 0.5124056759448475 , 0.0325068593353589 , 0.0325068593353589 , 0.0001329751563098 );
            ir->AddPentPoint( 313 , 0.5124056759448475 , 0.2112903026922174 , 0.0325068593353589 , 0.0325068593353589 , 0.0001329751563098 );
            ir->AddPentPoint( 314 , 0.2112903026922174 , 0.2112903026922174 , 0.0325068593353589 , 0.0325068593353589 , 0.0001329751563098 );
            ir->AddPentPoint( 315 , 0.0002012214607435 , 0.1210433476083445 , 0.3020268264169776 , 0.5765273830531908 , 0.0000083626065639 );
            ir->AddPentPoint( 316 , 0.1210433476083445 , 0.0002012214607435 , 0.3020268264169776 , 0.5765273830531908 , 0.0000083626065639 );
            ir->AddPentPoint( 317 , 0.0002012214607435 , 0.0002012214607435 , 0.3020268264169776 , 0.5765273830531908 , 0.0000083626065639 );
            ir->AddPentPoint( 318 , 0.0002012214607435 , 0.3020268264169776 , 0.1210433476083445 , 0.5765273830531908 , 0.0000083626065639 );
            ir->AddPentPoint( 319 , 0.3020268264169776 , 0.0002012214607435 , 0.1210433476083445 , 0.5765273830531908 , 0.0000083626065639 );
            ir->AddPentPoint( 320 , 0.0002012214607435 , 0.0002012214607435 , 0.1210433476083445 , 0.5765273830531908 , 0.0000083626065639 );
            ir->AddPentPoint( 321 , 0.1210433476083445 , 0.3020268264169776 , 0.0002012214607435 , 0.5765273830531908 , 0.0000083626065639 );
            ir->AddPentPoint( 322 , 0.0002012214607435 , 0.3020268264169776 , 0.0002012214607435 , 0.5765273830531908 , 0.0000083626065639 );
            ir->AddPentPoint( 323 , 0.3020268264169776 , 0.1210433476083445 , 0.0002012214607435 , 0.5765273830531908 , 0.0000083626065639 );
            ir->AddPentPoint( 324 , 0.0002012214607435 , 0.1210433476083445 , 0.0002012214607435 , 0.5765273830531908 , 0.0000083626065639 );
            ir->AddPentPoint( 325 , 0.3020268264169776 , 0.0002012214607435 , 0.0002012214607435 , 0.5765273830531908 , 0.0000083626065639 );
            ir->AddPentPoint( 326 , 0.1210433476083445 , 0.0002012214607435 , 0.0002012214607435 , 0.5765273830531908 , 0.0000083626065639 );
            ir->AddPentPoint( 327 , 0.0002012214607435 , 0.1210433476083445 , 0.5765273830531908 , 0.3020268264169776 , 0.0000083626065639 );
            ir->AddPentPoint( 328 , 0.1210433476083445 , 0.0002012214607435 , 0.5765273830531908 , 0.3020268264169776 , 0.0000083626065639 );
            ir->AddPentPoint( 329 , 0.0002012214607435 , 0.0002012214607435 , 0.5765273830531908 , 0.3020268264169776 , 0.0000083626065639 );
            ir->AddPentPoint( 330 , 0.0002012214607435 , 0.5765273830531908 , 0.1210433476083445 , 0.3020268264169776 , 0.0000083626065639 );
            ir->AddPentPoint( 331 , 0.5765273830531908 , 0.0002012214607435 , 0.1210433476083445 , 0.3020268264169776 , 0.0000083626065639 );
            ir->AddPentPoint( 332 , 0.0002012214607435 , 0.0002012214607435 , 0.1210433476083445 , 0.3020268264169776 , 0.0000083626065639 );
            ir->AddPentPoint( 333 , 0.1210433476083445 , 0.5765273830531908 , 0.0002012214607435 , 0.3020268264169776 , 0.0000083626065639 );
            ir->AddPentPoint( 334 , 0.0002012214607435 , 0.5765273830531908 , 0.0002012214607435 , 0.3020268264169776 , 0.0000083626065639 );
            ir->AddPentPoint( 335 , 0.5765273830531908 , 0.1210433476083445 , 0.0002012214607435 , 0.3020268264169776 , 0.0000083626065639 );
            ir->AddPentPoint( 336 , 0.0002012214607435 , 0.1210433476083445 , 0.0002012214607435 , 0.3020268264169776 , 0.0000083626065639 );
            ir->AddPentPoint( 337 , 0.5765273830531908 , 0.0002012214607435 , 0.0002012214607435 , 0.3020268264169776 , 0.0000083626065639 );
            ir->AddPentPoint( 338 , 0.1210433476083445 , 0.0002012214607435 , 0.0002012214607435 , 0.3020268264169776 , 0.0000083626065639 );
            ir->AddPentPoint( 339 , 0.0002012214607435 , 0.3020268264169776 , 0.5765273830531908 , 0.1210433476083445 , 0.0000083626065639 );
            ir->AddPentPoint( 340 , 0.3020268264169776 , 0.0002012214607435 , 0.5765273830531908 , 0.1210433476083445 , 0.0000083626065639 );
            ir->AddPentPoint( 341 , 0.0002012214607435 , 0.0002012214607435 , 0.5765273830531908 , 0.1210433476083445 , 0.0000083626065639 );
            ir->AddPentPoint( 342 , 0.0002012214607435 , 0.5765273830531908 , 0.3020268264169776 , 0.1210433476083445 , 0.0000083626065639 );
            ir->AddPentPoint( 343 , 0.5765273830531908 , 0.0002012214607435 , 0.3020268264169776 , 0.1210433476083445 , 0.0000083626065639 );
            ir->AddPentPoint( 344 , 0.0002012214607435 , 0.0002012214607435 , 0.3020268264169776 , 0.1210433476083445 , 0.0000083626065639 );
            ir->AddPentPoint( 345 , 0.3020268264169776 , 0.5765273830531908 , 0.0002012214607435 , 0.1210433476083445 , 0.0000083626065639 );
            ir->AddPentPoint( 346 , 0.0002012214607435 , 0.5765273830531908 , 0.0002012214607435 , 0.1210433476083445 , 0.0000083626065639 );
            ir->AddPentPoint( 347 , 0.5765273830531908 , 0.3020268264169776 , 0.0002012214607435 , 0.1210433476083445 , 0.0000083626065639 );
            ir->AddPentPoint( 348 , 0.0002012214607435 , 0.3020268264169776 , 0.0002012214607435 , 0.1210433476083445 , 0.0000083626065639 );
            ir->AddPentPoint( 349 , 0.5765273830531908 , 0.0002012214607435 , 0.0002012214607435 , 0.1210433476083445 , 0.0000083626065639 );
            ir->AddPentPoint( 350 , 0.3020268264169776 , 0.0002012214607435 , 0.0002012214607435 , 0.1210433476083445 , 0.0000083626065639 );
            ir->AddPentPoint( 351 , 0.1210433476083445 , 0.3020268264169776 , 0.5765273830531908 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 352 , 0.0002012214607435 , 0.3020268264169776 , 0.5765273830531908 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 353 , 0.3020268264169776 , 0.1210433476083445 , 0.5765273830531908 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 354 , 0.0002012214607435 , 0.1210433476083445 , 0.5765273830531908 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 355 , 0.3020268264169776 , 0.0002012214607435 , 0.5765273830531908 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 356 , 0.1210433476083445 , 0.0002012214607435 , 0.5765273830531908 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 357 , 0.1210433476083445 , 0.5765273830531908 , 0.3020268264169776 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 358 , 0.0002012214607435 , 0.5765273830531908 , 0.3020268264169776 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 359 , 0.5765273830531908 , 0.1210433476083445 , 0.3020268264169776 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 360 , 0.0002012214607435 , 0.1210433476083445 , 0.3020268264169776 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 361 , 0.5765273830531908 , 0.0002012214607435 , 0.3020268264169776 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 362 , 0.1210433476083445 , 0.0002012214607435 , 0.3020268264169776 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 363 , 0.3020268264169776 , 0.5765273830531908 , 0.1210433476083445 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 364 , 0.0002012214607435 , 0.5765273830531908 , 0.1210433476083445 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 365 , 0.5765273830531908 , 0.3020268264169776 , 0.1210433476083445 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 366 , 0.0002012214607435 , 0.3020268264169776 , 0.1210433476083445 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 367 , 0.5765273830531908 , 0.0002012214607435 , 0.1210433476083445 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 368 , 0.3020268264169776 , 0.0002012214607435 , 0.1210433476083445 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 369 , 0.3020268264169776 , 0.5765273830531908 , 0.0002012214607435 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 370 , 0.1210433476083445 , 0.5765273830531908 , 0.0002012214607435 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 371 , 0.5765273830531908 , 0.3020268264169776 , 0.0002012214607435 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 372 , 0.1210433476083445 , 0.3020268264169776 , 0.0002012214607435 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 373 , 0.5765273830531908 , 0.1210433476083445 , 0.0002012214607435 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 374 , 0.3020268264169776 , 0.1210433476083445 , 0.0002012214607435 , 0.0002012214607435 , 0.0000083626065639 );
            ir->AddPentPoint( 375 , 0.0739363213733936 , 0.0000000000000000 , 0.2809706979708689 , 0.5711566592823440 , 0.0000366901472403 );
            ir->AddPentPoint( 376 , 0.0000000000000000 , 0.0739363213733936 , 0.2809706979708689 , 0.5711566592823440 , 0.0000366901472403 );
            ir->AddPentPoint( 377 , 0.0739363213733936 , 0.0739363213733936 , 0.2809706979708689 , 0.5711566592823440 , 0.0000366901472403 );
            ir->AddPentPoint( 378 , 0.0739363213733936 , 0.2809706979708689 , 0.0000000000000000 , 0.5711566592823440 , 0.0000366901472403 );
            ir->AddPentPoint( 379 , 0.2809706979708689 , 0.0739363213733936 , 0.0000000000000000 , 0.5711566592823440 , 0.0000366901472403 );
            ir->AddPentPoint( 380 , 0.0739363213733936 , 0.0739363213733936 , 0.0000000000000000 , 0.5711566592823440 , 0.0000366901472403 );
            ir->AddPentPoint( 381 , 0.0000000000000000 , 0.2809706979708689 , 0.0739363213733936 , 0.5711566592823440 , 0.0000366901472403 );
            ir->AddPentPoint( 382 , 0.0739363213733936 , 0.2809706979708689 , 0.0739363213733936 , 0.5711566592823440 , 0.0000366901472403 );
            ir->AddPentPoint( 383 , 0.2809706979708689 , 0.0000000000000000 , 0.0739363213733936 , 0.5711566592823440 , 0.0000366901472403 );
            ir->AddPentPoint( 384 , 0.0739363213733936 , 0.0000000000000000 , 0.0739363213733936 , 0.5711566592823440 , 0.0000366901472403 );
            ir->AddPentPoint( 385 , 0.2809706979708689 , 0.0739363213733936 , 0.0739363213733936 , 0.5711566592823440 , 0.0000366901472403 );
            ir->AddPentPoint( 386 , 0.0000000000000000 , 0.0739363213733936 , 0.0739363213733936 , 0.5711566592823440 , 0.0000366901472403 );
            ir->AddPentPoint( 387 , 0.0739363213733936 , 0.0000000000000000 , 0.5711566592823440 , 0.2809706979708689 , 0.0000366901472403 );
            ir->AddPentPoint( 388 , 0.0000000000000000 , 0.0739363213733936 , 0.5711566592823440 , 0.2809706979708689 , 0.0000366901472403 );
            ir->AddPentPoint( 389 , 0.0739363213733936 , 0.0739363213733936 , 0.5711566592823440 , 0.2809706979708689 , 0.0000366901472403 );
            ir->AddPentPoint( 390 , 0.0739363213733936 , 0.5711566592823440 , 0.0000000000000000 , 0.2809706979708689 , 0.0000366901472403 );
            ir->AddPentPoint( 391 , 0.5711566592823440 , 0.0739363213733936 , 0.0000000000000000 , 0.2809706979708689 , 0.0000366901472403 );
            ir->AddPentPoint( 392 , 0.0739363213733936 , 0.0739363213733936 , 0.0000000000000000 , 0.2809706979708689 , 0.0000366901472403 );
            ir->AddPentPoint( 393 , 0.0000000000000000 , 0.5711566592823440 , 0.0739363213733936 , 0.2809706979708689 , 0.0000366901472403 );
            ir->AddPentPoint( 394 , 0.0739363213733936 , 0.5711566592823440 , 0.0739363213733936 , 0.2809706979708689 , 0.0000366901472403 );
            ir->AddPentPoint( 395 , 0.5711566592823440 , 0.0000000000000000 , 0.0739363213733936 , 0.2809706979708689 , 0.0000366901472403 );
            ir->AddPentPoint( 396 , 0.0739363213733936 , 0.0000000000000000 , 0.0739363213733936 , 0.2809706979708689 , 0.0000366901472403 );
            ir->AddPentPoint( 397 , 0.5711566592823440 , 0.0739363213733936 , 0.0739363213733936 , 0.2809706979708689 , 0.0000366901472403 );
            ir->AddPentPoint( 398 , 0.0000000000000000 , 0.0739363213733936 , 0.0739363213733936 , 0.2809706979708689 , 0.0000366901472403 );
            ir->AddPentPoint( 399 , 0.0739363213733936 , 0.2809706979708689 , 0.5711566592823440 , 0.0000000000000000 , 0.0000366901472403 );
            ir->AddPentPoint( 400 , 0.2809706979708689 , 0.0739363213733936 , 0.5711566592823440 , 0.0000000000000000 , 0.0000366901472403 );
            ir->AddPentPoint( 401 , 0.0739363213733936 , 0.0739363213733936 , 0.5711566592823440 , 0.0000000000000000 , 0.0000366901472403 );
            ir->AddPentPoint( 402 , 0.0739363213733936 , 0.5711566592823440 , 0.2809706979708689 , 0.0000000000000000 , 0.0000366901472403 );
            ir->AddPentPoint( 403 , 0.5711566592823440 , 0.0739363213733936 , 0.2809706979708689 , 0.0000000000000000 , 0.0000366901472403 );
            ir->AddPentPoint( 404 , 0.0739363213733936 , 0.0739363213733936 , 0.2809706979708689 , 0.0000000000000000 , 0.0000366901472403 );
            ir->AddPentPoint( 405 , 0.2809706979708689 , 0.5711566592823440 , 0.0739363213733936 , 0.0000000000000000 , 0.0000366901472403 );
            ir->AddPentPoint( 406 , 0.0739363213733936 , 0.5711566592823440 , 0.0739363213733936 , 0.0000000000000000 , 0.0000366901472403 );
            ir->AddPentPoint( 407 , 0.5711566592823440 , 0.2809706979708689 , 0.0739363213733936 , 0.0000000000000000 , 0.0000366901472403 );
            ir->AddPentPoint( 408 , 0.0739363213733936 , 0.2809706979708689 , 0.0739363213733936 , 0.0000000000000000 , 0.0000366901472403 );
            ir->AddPentPoint( 409 , 0.5711566592823440 , 0.0739363213733936 , 0.0739363213733936 , 0.0000000000000000 , 0.0000366901472403 );
            ir->AddPentPoint( 410 , 0.2809706979708689 , 0.0739363213733936 , 0.0739363213733936 , 0.0000000000000000 , 0.0000366901472403 );
            ir->AddPentPoint( 411 , 0.0000000000000000 , 0.2809706979708689 , 0.5711566592823440 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 412 , 0.0739363213733936 , 0.2809706979708689 , 0.5711566592823440 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 413 , 0.2809706979708689 , 0.0000000000000000 , 0.5711566592823440 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 414 , 0.0739363213733936 , 0.0000000000000000 , 0.5711566592823440 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 415 , 0.2809706979708689 , 0.0739363213733936 , 0.5711566592823440 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 416 , 0.0000000000000000 , 0.0739363213733936 , 0.5711566592823440 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 417 , 0.0000000000000000 , 0.5711566592823440 , 0.2809706979708689 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 418 , 0.0739363213733936 , 0.5711566592823440 , 0.2809706979708689 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 419 , 0.5711566592823440 , 0.0000000000000000 , 0.2809706979708689 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 420 , 0.0739363213733936 , 0.0000000000000000 , 0.2809706979708689 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 421 , 0.5711566592823440 , 0.0739363213733936 , 0.2809706979708689 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 422 , 0.0000000000000000 , 0.0739363213733936 , 0.2809706979708689 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 423 , 0.2809706979708689 , 0.5711566592823440 , 0.0000000000000000 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 424 , 0.0739363213733936 , 0.5711566592823440 , 0.0000000000000000 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 425 , 0.5711566592823440 , 0.2809706979708689 , 0.0000000000000000 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 426 , 0.0739363213733936 , 0.2809706979708689 , 0.0000000000000000 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 427 , 0.5711566592823440 , 0.0739363213733936 , 0.0000000000000000 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 428 , 0.2809706979708689 , 0.0739363213733936 , 0.0000000000000000 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 429 , 0.2809706979708689 , 0.5711566592823440 , 0.0739363213733936 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 430 , 0.0000000000000000 , 0.5711566592823440 , 0.0739363213733936 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 431 , 0.5711566592823440 , 0.2809706979708689 , 0.0739363213733936 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 432 , 0.0000000000000000 , 0.2809706979708689 , 0.0739363213733936 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 433 , 0.5711566592823440 , 0.0000000000000000 , 0.0739363213733936 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 434 , 0.2809706979708689 , 0.0000000000000000 , 0.0739363213733936 , 0.0739363213733936 , 0.0000366901472403 );
            ir->AddPentPoint( 435 , 0.1463142972627170 , 0.0323129857983894 , 0.2546936679720009 , 0.4203647517041757 , 0.0001463144167375 );
            ir->AddPentPoint( 436 , 0.0323129857983894 , 0.1463142972627170 , 0.2546936679720009 , 0.4203647517041757 , 0.0001463144167375 );
            ir->AddPentPoint( 437 , 0.1463142972627170 , 0.1463142972627170 , 0.2546936679720009 , 0.4203647517041757 , 0.0001463144167375 );
            ir->AddPentPoint( 438 , 0.1463142972627170 , 0.2546936679720009 , 0.0323129857983894 , 0.4203647517041757 , 0.0001463144167375 );
            ir->AddPentPoint( 439 , 0.2546936679720009 , 0.1463142972627170 , 0.0323129857983894 , 0.4203647517041757 , 0.0001463144167375 );
            ir->AddPentPoint( 440 , 0.1463142972627170 , 0.1463142972627170 , 0.0323129857983894 , 0.4203647517041757 , 0.0001463144167375 );
            ir->AddPentPoint( 441 , 0.0323129857983894 , 0.2546936679720009 , 0.1463142972627170 , 0.4203647517041757 , 0.0001463144167375 );
            ir->AddPentPoint( 442 , 0.1463142972627170 , 0.2546936679720009 , 0.1463142972627170 , 0.4203647517041757 , 0.0001463144167375 );
            ir->AddPentPoint( 443 , 0.2546936679720009 , 0.0323129857983894 , 0.1463142972627170 , 0.4203647517041757 , 0.0001463144167375 );
            ir->AddPentPoint( 444 , 0.1463142972627170 , 0.0323129857983894 , 0.1463142972627170 , 0.4203647517041757 , 0.0001463144167375 );
            ir->AddPentPoint( 445 , 0.2546936679720009 , 0.1463142972627170 , 0.1463142972627170 , 0.4203647517041757 , 0.0001463144167375 );
            ir->AddPentPoint( 446 , 0.0323129857983894 , 0.1463142972627170 , 0.1463142972627170 , 0.4203647517041757 , 0.0001463144167375 );
            ir->AddPentPoint( 447 , 0.1463142972627170 , 0.0323129857983894 , 0.4203647517041757 , 0.2546936679720009 , 0.0001463144167375 );
            ir->AddPentPoint( 448 , 0.0323129857983894 , 0.1463142972627170 , 0.4203647517041757 , 0.2546936679720009 , 0.0001463144167375 );
            ir->AddPentPoint( 449 , 0.1463142972627170 , 0.1463142972627170 , 0.4203647517041757 , 0.2546936679720009 , 0.0001463144167375 );
            ir->AddPentPoint( 450 , 0.1463142972627170 , 0.4203647517041757 , 0.0323129857983894 , 0.2546936679720009 , 0.0001463144167375 );
            ir->AddPentPoint( 451 , 0.4203647517041757 , 0.1463142972627170 , 0.0323129857983894 , 0.2546936679720009 , 0.0001463144167375 );
            ir->AddPentPoint( 452 , 0.1463142972627170 , 0.1463142972627170 , 0.0323129857983894 , 0.2546936679720009 , 0.0001463144167375 );
            ir->AddPentPoint( 453 , 0.0323129857983894 , 0.4203647517041757 , 0.1463142972627170 , 0.2546936679720009 , 0.0001463144167375 );
            ir->AddPentPoint( 454 , 0.1463142972627170 , 0.4203647517041757 , 0.1463142972627170 , 0.2546936679720009 , 0.0001463144167375 );
            ir->AddPentPoint( 455 , 0.4203647517041757 , 0.0323129857983894 , 0.1463142972627170 , 0.2546936679720009 , 0.0001463144167375 );
            ir->AddPentPoint( 456 , 0.1463142972627170 , 0.0323129857983894 , 0.1463142972627170 , 0.2546936679720009 , 0.0001463144167375 );
            ir->AddPentPoint( 457 , 0.4203647517041757 , 0.1463142972627170 , 0.1463142972627170 , 0.2546936679720009 , 0.0001463144167375 );
            ir->AddPentPoint( 458 , 0.0323129857983894 , 0.1463142972627170 , 0.1463142972627170 , 0.2546936679720009 , 0.0001463144167375 );
            ir->AddPentPoint( 459 , 0.1463142972627170 , 0.2546936679720009 , 0.4203647517041757 , 0.0323129857983894 , 0.0001463144167375 );
            ir->AddPentPoint( 460 , 0.2546936679720009 , 0.1463142972627170 , 0.4203647517041757 , 0.0323129857983894 , 0.0001463144167375 );
            ir->AddPentPoint( 461 , 0.1463142972627170 , 0.1463142972627170 , 0.4203647517041757 , 0.0323129857983894 , 0.0001463144167375 );
            ir->AddPentPoint( 462 , 0.1463142972627170 , 0.4203647517041757 , 0.2546936679720009 , 0.0323129857983894 , 0.0001463144167375 );
            ir->AddPentPoint( 463 , 0.4203647517041757 , 0.1463142972627170 , 0.2546936679720009 , 0.0323129857983894 , 0.0001463144167375 );
            ir->AddPentPoint( 464 , 0.1463142972627170 , 0.1463142972627170 , 0.2546936679720009 , 0.0323129857983894 , 0.0001463144167375 );
            ir->AddPentPoint( 465 , 0.2546936679720009 , 0.4203647517041757 , 0.1463142972627170 , 0.0323129857983894 , 0.0001463144167375 );
            ir->AddPentPoint( 466 , 0.1463142972627170 , 0.4203647517041757 , 0.1463142972627170 , 0.0323129857983894 , 0.0001463144167375 );
            ir->AddPentPoint( 467 , 0.4203647517041757 , 0.2546936679720009 , 0.1463142972627170 , 0.0323129857983894 , 0.0001463144167375 );
            ir->AddPentPoint( 468 , 0.1463142972627170 , 0.2546936679720009 , 0.1463142972627170 , 0.0323129857983894 , 0.0001463144167375 );
            ir->AddPentPoint( 469 , 0.4203647517041757 , 0.1463142972627170 , 0.1463142972627170 , 0.0323129857983894 , 0.0001463144167375 );
            ir->AddPentPoint( 470 , 0.2546936679720009 , 0.1463142972627170 , 0.1463142972627170 , 0.0323129857983894 , 0.0001463144167375 );
            ir->AddPentPoint( 471 , 0.0323129857983894 , 0.2546936679720009 , 0.4203647517041757 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 472 , 0.1463142972627170 , 0.2546936679720009 , 0.4203647517041757 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 473 , 0.2546936679720009 , 0.0323129857983894 , 0.4203647517041757 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 474 , 0.1463142972627170 , 0.0323129857983894 , 0.4203647517041757 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 475 , 0.2546936679720009 , 0.1463142972627170 , 0.4203647517041757 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 476 , 0.0323129857983894 , 0.1463142972627170 , 0.4203647517041757 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 477 , 0.0323129857983894 , 0.4203647517041757 , 0.2546936679720009 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 478 , 0.1463142972627170 , 0.4203647517041757 , 0.2546936679720009 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 479 , 0.4203647517041757 , 0.0323129857983894 , 0.2546936679720009 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 480 , 0.1463142972627170 , 0.0323129857983894 , 0.2546936679720009 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 481 , 0.4203647517041757 , 0.1463142972627170 , 0.2546936679720009 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 482 , 0.0323129857983894 , 0.1463142972627170 , 0.2546936679720009 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 483 , 0.2546936679720009 , 0.4203647517041757 , 0.0323129857983894 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 484 , 0.1463142972627170 , 0.4203647517041757 , 0.0323129857983894 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 485 , 0.4203647517041757 , 0.2546936679720009 , 0.0323129857983894 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 486 , 0.1463142972627170 , 0.2546936679720009 , 0.0323129857983894 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 487 , 0.4203647517041757 , 0.1463142972627170 , 0.0323129857983894 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 488 , 0.2546936679720009 , 0.1463142972627170 , 0.0323129857983894 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 489 , 0.2546936679720009 , 0.4203647517041757 , 0.1463142972627170 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 490 , 0.0323129857983894 , 0.4203647517041757 , 0.1463142972627170 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 491 , 0.4203647517041757 , 0.2546936679720009 , 0.1463142972627170 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 492 , 0.0323129857983894 , 0.2546936679720009 , 0.1463142972627170 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 493 , 0.4203647517041757 , 0.0323129857983894 , 0.1463142972627170 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 494 , 0.2546936679720009 , 0.0323129857983894 , 0.1463142972627170 , 0.1463142972627170 , 0.0001463144167375 );
            ir->AddPentPoint( 495 , 0.3534152321648770 , 0.0056219112112306 , 0.0821056872809434 , 0.2054419371780721 , 0.0000551833542664 );
            ir->AddPentPoint( 496 , 0.0056219112112306 , 0.3534152321648770 , 0.0821056872809434 , 0.2054419371780721 , 0.0000551833542664 );
            ir->AddPentPoint( 497 , 0.3534152321648770 , 0.3534152321648770 , 0.0821056872809434 , 0.2054419371780721 , 0.0000551833542664 );
            ir->AddPentPoint( 498 , 0.3534152321648770 , 0.0821056872809434 , 0.0056219112112306 , 0.2054419371780721 , 0.0000551833542664 );
            ir->AddPentPoint( 499 , 0.0821056872809434 , 0.3534152321648770 , 0.0056219112112306 , 0.2054419371780721 , 0.0000551833542664 );
            ir->AddPentPoint( 500 , 0.3534152321648770 , 0.3534152321648770 , 0.0056219112112306 , 0.2054419371780721 , 0.0000551833542664 );
            ir->AddPentPoint( 501 , 0.0056219112112306 , 0.0821056872809434 , 0.3534152321648770 , 0.2054419371780721 , 0.0000551833542664 );
            ir->AddPentPoint( 502 , 0.3534152321648770 , 0.0821056872809434 , 0.3534152321648770 , 0.2054419371780721 , 0.0000551833542664 );
            ir->AddPentPoint( 503 , 0.0821056872809434 , 0.0056219112112306 , 0.3534152321648770 , 0.2054419371780721 , 0.0000551833542664 );
            ir->AddPentPoint( 504 , 0.3534152321648770 , 0.0056219112112306 , 0.3534152321648770 , 0.2054419371780721 , 0.0000551833542664 );
            ir->AddPentPoint( 505 , 0.0821056872809434 , 0.3534152321648770 , 0.3534152321648770 , 0.2054419371780721 , 0.0000551833542664 );
            ir->AddPentPoint( 506 , 0.0056219112112306 , 0.3534152321648770 , 0.3534152321648770 , 0.2054419371780721 , 0.0000551833542664 );
            ir->AddPentPoint( 507 , 0.3534152321648770 , 0.0056219112112306 , 0.2054419371780721 , 0.0821056872809434 , 0.0000551833542664 );
            ir->AddPentPoint( 508 , 0.0056219112112306 , 0.3534152321648770 , 0.2054419371780721 , 0.0821056872809434 , 0.0000551833542664 );
            ir->AddPentPoint( 509 , 0.3534152321648770 , 0.3534152321648770 , 0.2054419371780721 , 0.0821056872809434 , 0.0000551833542664 );
            ir->AddPentPoint( 510 , 0.3534152321648770 , 0.2054419371780721 , 0.0056219112112306 , 0.0821056872809434 , 0.0000551833542664 );
            ir->AddPentPoint( 511 , 0.2054419371780721 , 0.3534152321648770 , 0.0056219112112306 , 0.0821056872809434 , 0.0000551833542664 );
            ir->AddPentPoint( 512 , 0.3534152321648770 , 0.3534152321648770 , 0.0056219112112306 , 0.0821056872809434 , 0.0000551833542664 );
            ir->AddPentPoint( 513 , 0.0056219112112306 , 0.2054419371780721 , 0.3534152321648770 , 0.0821056872809434 , 0.0000551833542664 );
            ir->AddPentPoint( 514 , 0.3534152321648770 , 0.2054419371780721 , 0.3534152321648770 , 0.0821056872809434 , 0.0000551833542664 );
            ir->AddPentPoint( 515 , 0.2054419371780721 , 0.0056219112112306 , 0.3534152321648770 , 0.0821056872809434 , 0.0000551833542664 );
            ir->AddPentPoint( 516 , 0.3534152321648770 , 0.0056219112112306 , 0.3534152321648770 , 0.0821056872809434 , 0.0000551833542664 );
            ir->AddPentPoint( 517 , 0.2054419371780721 , 0.3534152321648770 , 0.3534152321648770 , 0.0821056872809434 , 0.0000551833542664 );
            ir->AddPentPoint( 518 , 0.0056219112112306 , 0.3534152321648770 , 0.3534152321648770 , 0.0821056872809434 , 0.0000551833542664 );
            ir->AddPentPoint( 519 , 0.3534152321648770 , 0.0821056872809434 , 0.2054419371780721 , 0.0056219112112306 , 0.0000551833542664 );
            ir->AddPentPoint( 520 , 0.0821056872809434 , 0.3534152321648770 , 0.2054419371780721 , 0.0056219112112306 , 0.0000551833542664 );
            ir->AddPentPoint( 521 , 0.3534152321648770 , 0.3534152321648770 , 0.2054419371780721 , 0.0056219112112306 , 0.0000551833542664 );
            ir->AddPentPoint( 522 , 0.3534152321648770 , 0.2054419371780721 , 0.0821056872809434 , 0.0056219112112306 , 0.0000551833542664 );
            ir->AddPentPoint( 523 , 0.2054419371780721 , 0.3534152321648770 , 0.0821056872809434 , 0.0056219112112306 , 0.0000551833542664 );
            ir->AddPentPoint( 524 , 0.3534152321648770 , 0.3534152321648770 , 0.0821056872809434 , 0.0056219112112306 , 0.0000551833542664 );
            ir->AddPentPoint( 525 , 0.0821056872809434 , 0.2054419371780721 , 0.3534152321648770 , 0.0056219112112306 , 0.0000551833542664 );
            ir->AddPentPoint( 526 , 0.3534152321648770 , 0.2054419371780721 , 0.3534152321648770 , 0.0056219112112306 , 0.0000551833542664 );
            ir->AddPentPoint( 527 , 0.2054419371780721 , 0.0821056872809434 , 0.3534152321648770 , 0.0056219112112306 , 0.0000551833542664 );
            ir->AddPentPoint( 528 , 0.3534152321648770 , 0.0821056872809434 , 0.3534152321648770 , 0.0056219112112306 , 0.0000551833542664 );
            ir->AddPentPoint( 529 , 0.2054419371780721 , 0.3534152321648770 , 0.3534152321648770 , 0.0056219112112306 , 0.0000551833542664 );
            ir->AddPentPoint( 530 , 0.0821056872809434 , 0.3534152321648770 , 0.3534152321648770 , 0.0056219112112306 , 0.0000551833542664 );
            ir->AddPentPoint( 531 , 0.0056219112112306 , 0.0821056872809434 , 0.2054419371780721 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 532 , 0.3534152321648770 , 0.0821056872809434 , 0.2054419371780721 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 533 , 0.0821056872809434 , 0.0056219112112306 , 0.2054419371780721 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 534 , 0.3534152321648770 , 0.0056219112112306 , 0.2054419371780721 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 535 , 0.0821056872809434 , 0.3534152321648770 , 0.2054419371780721 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 536 , 0.0056219112112306 , 0.3534152321648770 , 0.2054419371780721 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 537 , 0.0056219112112306 , 0.2054419371780721 , 0.0821056872809434 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 538 , 0.3534152321648770 , 0.2054419371780721 , 0.0821056872809434 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 539 , 0.2054419371780721 , 0.0056219112112306 , 0.0821056872809434 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 540 , 0.3534152321648770 , 0.0056219112112306 , 0.0821056872809434 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 541 , 0.2054419371780721 , 0.3534152321648770 , 0.0821056872809434 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 542 , 0.0056219112112306 , 0.3534152321648770 , 0.0821056872809434 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 543 , 0.0821056872809434 , 0.2054419371780721 , 0.0056219112112306 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 544 , 0.3534152321648770 , 0.2054419371780721 , 0.0056219112112306 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 545 , 0.2054419371780721 , 0.0821056872809434 , 0.0056219112112306 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 546 , 0.3534152321648770 , 0.0821056872809434 , 0.0056219112112306 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 547 , 0.2054419371780721 , 0.3534152321648770 , 0.0056219112112306 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 548 , 0.0821056872809434 , 0.3534152321648770 , 0.0056219112112306 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 549 , 0.0821056872809434 , 0.2054419371780721 , 0.3534152321648770 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 550 , 0.0056219112112306 , 0.2054419371780721 , 0.3534152321648770 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 551 , 0.2054419371780721 , 0.0821056872809434 , 0.3534152321648770 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 552 , 0.0056219112112306 , 0.0821056872809434 , 0.3534152321648770 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 553 , 0.2054419371780721 , 0.0056219112112306 , 0.3534152321648770 , 0.3534152321648770 , 0.0000551833542664 );
            ir->AddPentPoint( 554 , 0.0821056872809434 , 0.0056219112112306 , 0.3534152321648770 , 0.3534152321648770 , 0.0000551833542664 );
            
            return ir;
            
        case 16:  // 1055 points - degree 16 --
            PentatopeIntRules[16] = ir = new IntegrationRule(1055);
            ir->AddPentPoint( 0 , 0.0707613942979452 , 0.0707613942979452 , 0.0707613942979452 , 0.7169544228082192 , 0.0000643856425658 );
            ir->AddPentPoint( 1 , 0.0707613942979452 , 0.0707613942979452 , 0.7169544228082192 , 0.0707613942979452 , 0.0000643856425658 );
            ir->AddPentPoint( 2 , 0.0707613942979452 , 0.7169544228082192 , 0.0707613942979452 , 0.0707613942979452 , 0.0000643856425658 );
            ir->AddPentPoint( 3 , 0.7169544228082192 , 0.0707613942979452 , 0.0707613942979452 , 0.0707613942979452 , 0.0000643856425658 );
            ir->AddPentPoint( 4 , 0.0707613942979452 , 0.0707613942979452 , 0.0707613942979452 , 0.0707613942979452 , 0.0000643856425658 );
            ir->AddPentPoint( 5 , 0.0164924596382321 , 0.0164924596382321 , 0.4752613105426517 , 0.4752613105426517 , 0.0000107669979120 );
            ir->AddPentPoint( 6 , 0.0164924596382321 , 0.4752613105426517 , 0.0164924596382321 , 0.4752613105426517 , 0.0000107669979120 );
            ir->AddPentPoint( 7 , 0.4752613105426517 , 0.0164924596382321 , 0.0164924596382321 , 0.4752613105426517 , 0.0000107669979120 );
            ir->AddPentPoint( 8 , 0.0164924596382321 , 0.0164924596382321 , 0.0164924596382321 , 0.4752613105426517 , 0.0000107669979120 );
            ir->AddPentPoint( 9 , 0.0164924596382321 , 0.4752613105426517 , 0.4752613105426517 , 0.0164924596382321 , 0.0000107669979120 );
            ir->AddPentPoint( 10 , 0.4752613105426517 , 0.0164924596382321 , 0.4752613105426517 , 0.0164924596382321 , 0.0000107669979120 );
            ir->AddPentPoint( 11 , 0.0164924596382321 , 0.0164924596382321 , 0.4752613105426517 , 0.0164924596382321 , 0.0000107669979120 );
            ir->AddPentPoint( 12 , 0.4752613105426517 , 0.4752613105426517 , 0.0164924596382321 , 0.0164924596382321 , 0.0000107669979120 );
            ir->AddPentPoint( 13 , 0.0164924596382321 , 0.4752613105426517 , 0.0164924596382321 , 0.0164924596382321 , 0.0000107669979120 );
            ir->AddPentPoint( 14 , 0.4752613105426517 , 0.0164924596382321 , 0.0164924596382321 , 0.0164924596382321 , 0.0000107669979120 );
            ir->AddPentPoint( 15 , 0.0536522941130353 , 0.0536522941130353 , 0.4195215588304471 , 0.4195215588304471 , 0.0000733749748583 );
            ir->AddPentPoint( 16 , 0.0536522941130353 , 0.4195215588304471 , 0.0536522941130353 , 0.4195215588304471 , 0.0000733749748583 );
            ir->AddPentPoint( 17 , 0.4195215588304471 , 0.0536522941130353 , 0.0536522941130353 , 0.4195215588304471 , 0.0000733749748583 );
            ir->AddPentPoint( 18 , 0.0536522941130353 , 0.0536522941130353 , 0.0536522941130353 , 0.4195215588304471 , 0.0000733749748583 );
            ir->AddPentPoint( 19 , 0.0536522941130353 , 0.4195215588304471 , 0.4195215588304471 , 0.0536522941130353 , 0.0000733749748583 );
            ir->AddPentPoint( 20 , 0.4195215588304471 , 0.0536522941130353 , 0.4195215588304471 , 0.0536522941130353 , 0.0000733749748583 );
            ir->AddPentPoint( 21 , 0.0536522941130353 , 0.0536522941130353 , 0.4195215588304471 , 0.0536522941130353 , 0.0000733749748583 );
            ir->AddPentPoint( 22 , 0.4195215588304471 , 0.4195215588304471 , 0.0536522941130353 , 0.0536522941130353 , 0.0000733749748583 );
            ir->AddPentPoint( 23 , 0.0536522941130353 , 0.4195215588304471 , 0.0536522941130353 , 0.0536522941130353 , 0.0000733749748583 );
            ir->AddPentPoint( 24 , 0.4195215588304471 , 0.0536522941130353 , 0.0536522941130353 , 0.0536522941130353 , 0.0000733749748583 );
            ir->AddPentPoint( 25 , 0.2423101553883049 , 0.2423101553883049 , 0.1365347669175427 , 0.1365347669175427 , 0.0001509931908024 );
            ir->AddPentPoint( 26 , 0.2423101553883049 , 0.1365347669175427 , 0.2423101553883049 , 0.1365347669175427 , 0.0001509931908024 );
            ir->AddPentPoint( 27 , 0.1365347669175427 , 0.2423101553883049 , 0.2423101553883049 , 0.1365347669175427 , 0.0001509931908024 );
            ir->AddPentPoint( 28 , 0.2423101553883049 , 0.2423101553883049 , 0.2423101553883049 , 0.1365347669175427 , 0.0001509931908024 );
            ir->AddPentPoint( 29 , 0.2423101553883049 , 0.1365347669175427 , 0.1365347669175427 , 0.2423101553883049 , 0.0001509931908024 );
            ir->AddPentPoint( 30 , 0.1365347669175427 , 0.2423101553883049 , 0.1365347669175427 , 0.2423101553883049 , 0.0001509931908024 );
            ir->AddPentPoint( 31 , 0.2423101553883049 , 0.2423101553883049 , 0.1365347669175427 , 0.2423101553883049 , 0.0001509931908024 );
            ir->AddPentPoint( 32 , 0.1365347669175427 , 0.1365347669175427 , 0.2423101553883049 , 0.2423101553883049 , 0.0001509931908024 );
            ir->AddPentPoint( 33 , 0.2423101553883049 , 0.1365347669175427 , 0.2423101553883049 , 0.2423101553883049 , 0.0001509931908024 );
            ir->AddPentPoint( 34 , 0.1365347669175427 , 0.2423101553883049 , 0.2423101553883049 , 0.2423101553883049 , 0.0001509931908024 );
            ir->AddPentPoint( 35 , 0.0057103322101029 , 0.0057103322101029 , 0.0526564020408650 , 0.9302126013288265 , 0.0000009301306146 );
            ir->AddPentPoint( 36 , 0.0057103322101029 , 0.0526564020408650 , 0.0057103322101029 , 0.9302126013288265 , 0.0000009301306146 );
            ir->AddPentPoint( 37 , 0.0526564020408650 , 0.0057103322101029 , 0.0057103322101029 , 0.9302126013288265 , 0.0000009301306146 );
            ir->AddPentPoint( 38 , 0.0057103322101029 , 0.0057103322101029 , 0.0057103322101029 , 0.9302126013288265 , 0.0000009301306146 );
            ir->AddPentPoint( 39 , 0.0057103322101029 , 0.0057103322101029 , 0.9302126013288265 , 0.0526564020408650 , 0.0000009301306146 );
            ir->AddPentPoint( 40 , 0.0057103322101029 , 0.9302126013288265 , 0.0057103322101029 , 0.0526564020408650 , 0.0000009301306146 );
            ir->AddPentPoint( 41 , 0.9302126013288265 , 0.0057103322101029 , 0.0057103322101029 , 0.0526564020408650 , 0.0000009301306146 );
            ir->AddPentPoint( 42 , 0.0057103322101029 , 0.0057103322101029 , 0.0057103322101029 , 0.0526564020408650 , 0.0000009301306146 );
            ir->AddPentPoint( 43 , 0.0057103322101029 , 0.0526564020408650 , 0.9302126013288265 , 0.0057103322101029 , 0.0000009301306146 );
            ir->AddPentPoint( 44 , 0.0526564020408650 , 0.0057103322101029 , 0.9302126013288265 , 0.0057103322101029 , 0.0000009301306146 );
            ir->AddPentPoint( 45 , 0.0057103322101029 , 0.0057103322101029 , 0.9302126013288265 , 0.0057103322101029 , 0.0000009301306146 );
            ir->AddPentPoint( 46 , 0.0057103322101029 , 0.9302126013288265 , 0.0526564020408650 , 0.0057103322101029 , 0.0000009301306146 );
            ir->AddPentPoint( 47 , 0.9302126013288265 , 0.0057103322101029 , 0.0526564020408650 , 0.0057103322101029 , 0.0000009301306146 );
            ir->AddPentPoint( 48 , 0.0057103322101029 , 0.0057103322101029 , 0.0526564020408650 , 0.0057103322101029 , 0.0000009301306146 );
            ir->AddPentPoint( 49 , 0.0526564020408650 , 0.9302126013288265 , 0.0057103322101029 , 0.0057103322101029 , 0.0000009301306146 );
            ir->AddPentPoint( 50 , 0.0057103322101029 , 0.9302126013288265 , 0.0057103322101029 , 0.0057103322101029 , 0.0000009301306146 );
            ir->AddPentPoint( 51 , 0.9302126013288265 , 0.0526564020408650 , 0.0057103322101029 , 0.0057103322101029 , 0.0000009301306146 );
            ir->AddPentPoint( 52 , 0.0057103322101029 , 0.0526564020408650 , 0.0057103322101029 , 0.0057103322101029 , 0.0000009301306146 );
            ir->AddPentPoint( 53 , 0.9302126013288265 , 0.0057103322101029 , 0.0057103322101029 , 0.0057103322101029 , 0.0000009301306146 );
            ir->AddPentPoint( 54 , 0.0526564020408650 , 0.0057103322101029 , 0.0057103322101029 , 0.0057103322101029 , 0.0000009301306146 );
            ir->AddPentPoint( 55 , 0.0146009267283278 , 0.0146009267283278 , 0.3210162946000336 , 0.6351809252149829 , 0.0000079791371467 );
            ir->AddPentPoint( 56 , 0.0146009267283278 , 0.3210162946000336 , 0.0146009267283278 , 0.6351809252149829 , 0.0000079791371467 );
            ir->AddPentPoint( 57 , 0.3210162946000336 , 0.0146009267283278 , 0.0146009267283278 , 0.6351809252149829 , 0.0000079791371467 );
            ir->AddPentPoint( 58 , 0.0146009267283278 , 0.0146009267283278 , 0.0146009267283278 , 0.6351809252149829 , 0.0000079791371467 );
            ir->AddPentPoint( 59 , 0.0146009267283278 , 0.0146009267283278 , 0.6351809252149829 , 0.3210162946000336 , 0.0000079791371467 );
            ir->AddPentPoint( 60 , 0.0146009267283278 , 0.6351809252149829 , 0.0146009267283278 , 0.3210162946000336 , 0.0000079791371467 );
            ir->AddPentPoint( 61 , 0.6351809252149829 , 0.0146009267283278 , 0.0146009267283278 , 0.3210162946000336 , 0.0000079791371467 );
            ir->AddPentPoint( 62 , 0.0146009267283278 , 0.0146009267283278 , 0.0146009267283278 , 0.3210162946000336 , 0.0000079791371467 );
            ir->AddPentPoint( 63 , 0.0146009267283278 , 0.3210162946000336 , 0.6351809252149829 , 0.0146009267283278 , 0.0000079791371467 );
            ir->AddPentPoint( 64 , 0.3210162946000336 , 0.0146009267283278 , 0.6351809252149829 , 0.0146009267283278 , 0.0000079791371467 );
            ir->AddPentPoint( 65 , 0.0146009267283278 , 0.0146009267283278 , 0.6351809252149829 , 0.0146009267283278 , 0.0000079791371467 );
            ir->AddPentPoint( 66 , 0.0146009267283278 , 0.6351809252149829 , 0.3210162946000336 , 0.0146009267283278 , 0.0000079791371467 );
            ir->AddPentPoint( 67 , 0.6351809252149829 , 0.0146009267283278 , 0.3210162946000336 , 0.0146009267283278 , 0.0000079791371467 );
            ir->AddPentPoint( 68 , 0.0146009267283278 , 0.0146009267283278 , 0.3210162946000336 , 0.0146009267283278 , 0.0000079791371467 );
            ir->AddPentPoint( 69 , 0.3210162946000336 , 0.6351809252149829 , 0.0146009267283278 , 0.0146009267283278 , 0.0000079791371467 );
            ir->AddPentPoint( 70 , 0.0146009267283278 , 0.6351809252149829 , 0.0146009267283278 , 0.0146009267283278 , 0.0000079791371467 );
            ir->AddPentPoint( 71 , 0.6351809252149829 , 0.3210162946000336 , 0.0146009267283278 , 0.0146009267283278 , 0.0000079791371467 );
            ir->AddPentPoint( 72 , 0.0146009267283278 , 0.3210162946000336 , 0.0146009267283278 , 0.0146009267283278 , 0.0000079791371467 );
            ir->AddPentPoint( 73 , 0.6351809252149829 , 0.0146009267283278 , 0.0146009267283278 , 0.0146009267283278 , 0.0000079791371467 );
            ir->AddPentPoint( 74 , 0.3210162946000336 , 0.0146009267283278 , 0.0146009267283278 , 0.0146009267283278 , 0.0000079791371467 );
            ir->AddPentPoint( 75 , 0.0167720160970359 , 0.0167720160970359 , 0.1553333404869812 , 0.7943506112219112 , 0.0000102707205107 );
            ir->AddPentPoint( 76 , 0.0167720160970359 , 0.1553333404869812 , 0.0167720160970359 , 0.7943506112219112 , 0.0000102707205107 );
            ir->AddPentPoint( 77 , 0.1553333404869812 , 0.0167720160970359 , 0.0167720160970359 , 0.7943506112219112 , 0.0000102707205107 );
            ir->AddPentPoint( 78 , 0.0167720160970359 , 0.0167720160970359 , 0.0167720160970359 , 0.7943506112219112 , 0.0000102707205107 );
            ir->AddPentPoint( 79 , 0.0167720160970359 , 0.0167720160970359 , 0.7943506112219112 , 0.1553333404869812 , 0.0000102707205107 );
            ir->AddPentPoint( 80 , 0.0167720160970359 , 0.7943506112219112 , 0.0167720160970359 , 0.1553333404869812 , 0.0000102707205107 );
            ir->AddPentPoint( 81 , 0.7943506112219112 , 0.0167720160970359 , 0.0167720160970359 , 0.1553333404869812 , 0.0000102707205107 );
            ir->AddPentPoint( 82 , 0.0167720160970359 , 0.0167720160970359 , 0.0167720160970359 , 0.1553333404869812 , 0.0000102707205107 );
            ir->AddPentPoint( 83 , 0.0167720160970359 , 0.1553333404869812 , 0.7943506112219112 , 0.0167720160970359 , 0.0000102707205107 );
            ir->AddPentPoint( 84 , 0.1553333404869812 , 0.0167720160970359 , 0.7943506112219112 , 0.0167720160970359 , 0.0000102707205107 );
            ir->AddPentPoint( 85 , 0.0167720160970359 , 0.0167720160970359 , 0.7943506112219112 , 0.0167720160970359 , 0.0000102707205107 );
            ir->AddPentPoint( 86 , 0.0167720160970359 , 0.7943506112219112 , 0.1553333404869812 , 0.0167720160970359 , 0.0000102707205107 );
            ir->AddPentPoint( 87 , 0.7943506112219112 , 0.0167720160970359 , 0.1553333404869812 , 0.0167720160970359 , 0.0000102707205107 );
            ir->AddPentPoint( 88 , 0.0167720160970359 , 0.0167720160970359 , 0.1553333404869812 , 0.0167720160970359 , 0.0000102707205107 );
            ir->AddPentPoint( 89 , 0.1553333404869812 , 0.7943506112219112 , 0.0167720160970359 , 0.0167720160970359 , 0.0000102707205107 );
            ir->AddPentPoint( 90 , 0.0167720160970359 , 0.7943506112219112 , 0.0167720160970359 , 0.0167720160970359 , 0.0000102707205107 );
            ir->AddPentPoint( 91 , 0.7943506112219112 , 0.1553333404869812 , 0.0167720160970359 , 0.0167720160970359 , 0.0000102707205107 );
            ir->AddPentPoint( 92 , 0.0167720160970359 , 0.1553333404869812 , 0.0167720160970359 , 0.0167720160970359 , 0.0000102707205107 );
            ir->AddPentPoint( 93 , 0.7943506112219112 , 0.0167720160970359 , 0.0167720160970359 , 0.0167720160970359 , 0.0000102707205107 );
            ir->AddPentPoint( 94 , 0.1553333404869812 , 0.0167720160970359 , 0.0167720160970359 , 0.0167720160970359 , 0.0000102707205107 );
            ir->AddPentPoint( 95 , 0.0333716873923710 , 0.0333716873923710 , 0.2927459459394552 , 0.6071389918834320 , 0.0000222604721338 );
            ir->AddPentPoint( 96 , 0.0333716873923710 , 0.2927459459394552 , 0.0333716873923710 , 0.6071389918834320 , 0.0000222604721338 );
            ir->AddPentPoint( 97 , 0.2927459459394552 , 0.0333716873923710 , 0.0333716873923710 , 0.6071389918834320 , 0.0000222604721338 );
            ir->AddPentPoint( 98 , 0.0333716873923710 , 0.0333716873923710 , 0.0333716873923710 , 0.6071389918834320 , 0.0000222604721338 );
            ir->AddPentPoint( 99 , 0.0333716873923710 , 0.0333716873923710 , 0.6071389918834320 , 0.2927459459394552 , 0.0000222604721338 );
            ir->AddPentPoint( 100 , 0.0333716873923710 , 0.6071389918834320 , 0.0333716873923710 , 0.2927459459394552 , 0.0000222604721338 );
            ir->AddPentPoint( 101 , 0.6071389918834320 , 0.0333716873923710 , 0.0333716873923710 , 0.2927459459394552 , 0.0000222604721338 );
            ir->AddPentPoint( 102 , 0.0333716873923710 , 0.0333716873923710 , 0.0333716873923710 , 0.2927459459394552 , 0.0000222604721338 );
            ir->AddPentPoint( 103 , 0.0333716873923710 , 0.2927459459394552 , 0.6071389918834320 , 0.0333716873923710 , 0.0000222604721338 );
            ir->AddPentPoint( 104 , 0.2927459459394552 , 0.0333716873923710 , 0.6071389918834320 , 0.0333716873923710 , 0.0000222604721338 );
            ir->AddPentPoint( 105 , 0.0333716873923710 , 0.0333716873923710 , 0.6071389918834320 , 0.0333716873923710 , 0.0000222604721338 );
            ir->AddPentPoint( 106 , 0.0333716873923710 , 0.6071389918834320 , 0.2927459459394552 , 0.0333716873923710 , 0.0000222604721338 );
            ir->AddPentPoint( 107 , 0.6071389918834320 , 0.0333716873923710 , 0.2927459459394552 , 0.0333716873923710 , 0.0000222604721338 );
            ir->AddPentPoint( 108 , 0.0333716873923710 , 0.0333716873923710 , 0.2927459459394552 , 0.0333716873923710 , 0.0000222604721338 );
            ir->AddPentPoint( 109 , 0.2927459459394552 , 0.6071389918834320 , 0.0333716873923710 , 0.0333716873923710 , 0.0000222604721338 );
            ir->AddPentPoint( 110 , 0.0333716873923710 , 0.6071389918834320 , 0.0333716873923710 , 0.0333716873923710 , 0.0000222604721338 );
            ir->AddPentPoint( 111 , 0.6071389918834320 , 0.2927459459394552 , 0.0333716873923710 , 0.0333716873923710 , 0.0000222604721338 );
            ir->AddPentPoint( 112 , 0.0333716873923710 , 0.2927459459394552 , 0.0333716873923710 , 0.0333716873923710 , 0.0000222604721338 );
            ir->AddPentPoint( 113 , 0.6071389918834320 , 0.0333716873923710 , 0.0333716873923710 , 0.0333716873923710 , 0.0000222604721338 );
            ir->AddPentPoint( 114 , 0.2927459459394552 , 0.0333716873923710 , 0.0333716873923710 , 0.0333716873923710 , 0.0000222604721338 );
            ir->AddPentPoint( 115 , 0.0396663043431213 , 0.0396663043431213 , 0.0072990563678513 , 0.8737020306027848 , 0.0000042150803680 );
            ir->AddPentPoint( 116 , 0.0396663043431213 , 0.0072990563678513 , 0.0396663043431213 , 0.8737020306027848 , 0.0000042150803680 );
            ir->AddPentPoint( 117 , 0.0072990563678513 , 0.0396663043431213 , 0.0396663043431213 , 0.8737020306027848 , 0.0000042150803680 );
            ir->AddPentPoint( 118 , 0.0396663043431213 , 0.0396663043431213 , 0.0396663043431213 , 0.8737020306027848 , 0.0000042150803680 );
            ir->AddPentPoint( 119 , 0.0396663043431213 , 0.0396663043431213 , 0.8737020306027848 , 0.0072990563678513 , 0.0000042150803680 );
            ir->AddPentPoint( 120 , 0.0396663043431213 , 0.8737020306027848 , 0.0396663043431213 , 0.0072990563678513 , 0.0000042150803680 );
            ir->AddPentPoint( 121 , 0.8737020306027848 , 0.0396663043431213 , 0.0396663043431213 , 0.0072990563678513 , 0.0000042150803680 );
            ir->AddPentPoint( 122 , 0.0396663043431213 , 0.0396663043431213 , 0.0396663043431213 , 0.0072990563678513 , 0.0000042150803680 );
            ir->AddPentPoint( 123 , 0.0396663043431213 , 0.0072990563678513 , 0.8737020306027848 , 0.0396663043431213 , 0.0000042150803680 );
            ir->AddPentPoint( 124 , 0.0072990563678513 , 0.0396663043431213 , 0.8737020306027848 , 0.0396663043431213 , 0.0000042150803680 );
            ir->AddPentPoint( 125 , 0.0396663043431213 , 0.0396663043431213 , 0.8737020306027848 , 0.0396663043431213 , 0.0000042150803680 );
            ir->AddPentPoint( 126 , 0.0396663043431213 , 0.8737020306027848 , 0.0072990563678513 , 0.0396663043431213 , 0.0000042150803680 );
            ir->AddPentPoint( 127 , 0.8737020306027848 , 0.0396663043431213 , 0.0072990563678513 , 0.0396663043431213 , 0.0000042150803680 );
            ir->AddPentPoint( 128 , 0.0396663043431213 , 0.0396663043431213 , 0.0072990563678513 , 0.0396663043431213 , 0.0000042150803680 );
            ir->AddPentPoint( 129 , 0.0072990563678513 , 0.8737020306027848 , 0.0396663043431213 , 0.0396663043431213 , 0.0000042150803680 );
            ir->AddPentPoint( 130 , 0.0396663043431213 , 0.8737020306027848 , 0.0396663043431213 , 0.0396663043431213 , 0.0000042150803680 );
            ir->AddPentPoint( 131 , 0.8737020306027848 , 0.0072990563678513 , 0.0396663043431213 , 0.0396663043431213 , 0.0000042150803680 );
            ir->AddPentPoint( 132 , 0.0396663043431213 , 0.0072990563678513 , 0.0396663043431213 , 0.0396663043431213 , 0.0000042150803680 );
            ir->AddPentPoint( 133 , 0.8737020306027848 , 0.0396663043431213 , 0.0396663043431213 , 0.0396663043431213 , 0.0000042150803680 );
            ir->AddPentPoint( 134 , 0.0072990563678513 , 0.0396663043431213 , 0.0396663043431213 , 0.0396663043431213 , 0.0000042150803680 );
            ir->AddPentPoint( 135 , 0.0624034595127863 , 0.0624034595127863 , 0.2344030462701640 , 0.5783865751914772 , 0.0000299311351699 );
            ir->AddPentPoint( 136 , 0.0624034595127863 , 0.2344030462701640 , 0.0624034595127863 , 0.5783865751914772 , 0.0000299311351699 );
            ir->AddPentPoint( 137 , 0.2344030462701640 , 0.0624034595127863 , 0.0624034595127863 , 0.5783865751914772 , 0.0000299311351699 );
            ir->AddPentPoint( 138 , 0.0624034595127863 , 0.0624034595127863 , 0.0624034595127863 , 0.5783865751914772 , 0.0000299311351699 );
            ir->AddPentPoint( 139 , 0.0624034595127863 , 0.0624034595127863 , 0.5783865751914772 , 0.2344030462701640 , 0.0000299311351699 );
            ir->AddPentPoint( 140 , 0.0624034595127863 , 0.5783865751914772 , 0.0624034595127863 , 0.2344030462701640 , 0.0000299311351699 );
            ir->AddPentPoint( 141 , 0.5783865751914772 , 0.0624034595127863 , 0.0624034595127863 , 0.2344030462701640 , 0.0000299311351699 );
            ir->AddPentPoint( 142 , 0.0624034595127863 , 0.0624034595127863 , 0.0624034595127863 , 0.2344030462701640 , 0.0000299311351699 );
            ir->AddPentPoint( 143 , 0.0624034595127863 , 0.2344030462701640 , 0.5783865751914772 , 0.0624034595127863 , 0.0000299311351699 );
            ir->AddPentPoint( 144 , 0.2344030462701640 , 0.0624034595127863 , 0.5783865751914772 , 0.0624034595127863 , 0.0000299311351699 );
            ir->AddPentPoint( 145 , 0.0624034595127863 , 0.0624034595127863 , 0.5783865751914772 , 0.0624034595127863 , 0.0000299311351699 );
            ir->AddPentPoint( 146 , 0.0624034595127863 , 0.5783865751914772 , 0.2344030462701640 , 0.0624034595127863 , 0.0000299311351699 );
            ir->AddPentPoint( 147 , 0.5783865751914772 , 0.0624034595127863 , 0.2344030462701640 , 0.0624034595127863 , 0.0000299311351699 );
            ir->AddPentPoint( 148 , 0.0624034595127863 , 0.0624034595127863 , 0.2344030462701640 , 0.0624034595127863 , 0.0000299311351699 );
            ir->AddPentPoint( 149 , 0.2344030462701640 , 0.5783865751914772 , 0.0624034595127863 , 0.0624034595127863 , 0.0000299311351699 );
            ir->AddPentPoint( 150 , 0.0624034595127863 , 0.5783865751914772 , 0.0624034595127863 , 0.0624034595127863 , 0.0000299311351699 );
            ir->AddPentPoint( 151 , 0.5783865751914772 , 0.2344030462701640 , 0.0624034595127863 , 0.0624034595127863 , 0.0000299311351699 );
            ir->AddPentPoint( 152 , 0.0624034595127863 , 0.2344030462701640 , 0.0624034595127863 , 0.0624034595127863 , 0.0000299311351699 );
            ir->AddPentPoint( 153 , 0.5783865751914772 , 0.0624034595127863 , 0.0624034595127863 , 0.0624034595127863 , 0.0000299311351699 );
            ir->AddPentPoint( 154 , 0.2344030462701640 , 0.0624034595127863 , 0.0624034595127863 , 0.0624034595127863 , 0.0000299311351699 );
            ir->AddPentPoint( 155 , 0.0832903691214466 , 0.0832903691214466 , 0.0144990370333865 , 0.7356298556022737 , 0.0000125836469210 );
            ir->AddPentPoint( 156 , 0.0832903691214466 , 0.0144990370333865 , 0.0832903691214466 , 0.7356298556022737 , 0.0000125836469210 );
            ir->AddPentPoint( 157 , 0.0144990370333865 , 0.0832903691214466 , 0.0832903691214466 , 0.7356298556022737 , 0.0000125836469210 );
            ir->AddPentPoint( 158 , 0.0832903691214466 , 0.0832903691214466 , 0.0832903691214466 , 0.7356298556022737 , 0.0000125836469210 );
            ir->AddPentPoint( 159 , 0.0832903691214466 , 0.0832903691214466 , 0.7356298556022737 , 0.0144990370333865 , 0.0000125836469210 );
            ir->AddPentPoint( 160 , 0.0832903691214466 , 0.7356298556022737 , 0.0832903691214466 , 0.0144990370333865 , 0.0000125836469210 );
            ir->AddPentPoint( 161 , 0.7356298556022737 , 0.0832903691214466 , 0.0832903691214466 , 0.0144990370333865 , 0.0000125836469210 );
            ir->AddPentPoint( 162 , 0.0832903691214466 , 0.0832903691214466 , 0.0832903691214466 , 0.0144990370333865 , 0.0000125836469210 );
            ir->AddPentPoint( 163 , 0.0832903691214466 , 0.0144990370333865 , 0.7356298556022737 , 0.0832903691214466 , 0.0000125836469210 );
            ir->AddPentPoint( 164 , 0.0144990370333865 , 0.0832903691214466 , 0.7356298556022737 , 0.0832903691214466 , 0.0000125836469210 );
            ir->AddPentPoint( 165 , 0.0832903691214466 , 0.0832903691214466 , 0.7356298556022737 , 0.0832903691214466 , 0.0000125836469210 );
            ir->AddPentPoint( 166 , 0.0832903691214466 , 0.7356298556022737 , 0.0144990370333865 , 0.0832903691214466 , 0.0000125836469210 );
            ir->AddPentPoint( 167 , 0.7356298556022737 , 0.0832903691214466 , 0.0144990370333865 , 0.0832903691214466 , 0.0000125836469210 );
            ir->AddPentPoint( 168 , 0.0832903691214466 , 0.0832903691214466 , 0.0144990370333865 , 0.0832903691214466 , 0.0000125836469210 );
            ir->AddPentPoint( 169 , 0.0144990370333865 , 0.7356298556022737 , 0.0832903691214466 , 0.0832903691214466 , 0.0000125836469210 );
            ir->AddPentPoint( 170 , 0.0832903691214466 , 0.7356298556022737 , 0.0832903691214466 , 0.0832903691214466 , 0.0000125836469210 );
            ir->AddPentPoint( 171 , 0.7356298556022737 , 0.0144990370333865 , 0.0832903691214466 , 0.0832903691214466 , 0.0000125836469210 );
            ir->AddPentPoint( 172 , 0.0832903691214466 , 0.0144990370333865 , 0.0832903691214466 , 0.0832903691214466 , 0.0000125836469210 );
            ir->AddPentPoint( 173 , 0.7356298556022737 , 0.0832903691214466 , 0.0832903691214466 , 0.0832903691214466 , 0.0000125836469210 );
            ir->AddPentPoint( 174 , 0.0144990370333865 , 0.0832903691214466 , 0.0832903691214466 , 0.0832903691214466 , 0.0000125836469210 );
            ir->AddPentPoint( 175 , 0.0934050772864460 , 0.0934050772864460 , 0.1779572463206686 , 0.5418275218199933 , 0.0000665877792863 );
            ir->AddPentPoint( 176 , 0.0934050772864460 , 0.1779572463206686 , 0.0934050772864460 , 0.5418275218199933 , 0.0000665877792863 );
            ir->AddPentPoint( 177 , 0.1779572463206686 , 0.0934050772864460 , 0.0934050772864460 , 0.5418275218199933 , 0.0000665877792863 );
            ir->AddPentPoint( 178 , 0.0934050772864460 , 0.0934050772864460 , 0.0934050772864460 , 0.5418275218199933 , 0.0000665877792863 );
            ir->AddPentPoint( 179 , 0.0934050772864460 , 0.0934050772864460 , 0.5418275218199933 , 0.1779572463206686 , 0.0000665877792863 );
            ir->AddPentPoint( 180 , 0.0934050772864460 , 0.5418275218199933 , 0.0934050772864460 , 0.1779572463206686 , 0.0000665877792863 );
            ir->AddPentPoint( 181 , 0.5418275218199933 , 0.0934050772864460 , 0.0934050772864460 , 0.1779572463206686 , 0.0000665877792863 );
            ir->AddPentPoint( 182 , 0.0934050772864460 , 0.0934050772864460 , 0.0934050772864460 , 0.1779572463206686 , 0.0000665877792863 );
            ir->AddPentPoint( 183 , 0.0934050772864460 , 0.1779572463206686 , 0.5418275218199933 , 0.0934050772864460 , 0.0000665877792863 );
            ir->AddPentPoint( 184 , 0.1779572463206686 , 0.0934050772864460 , 0.5418275218199933 , 0.0934050772864460 , 0.0000665877792863 );
            ir->AddPentPoint( 185 , 0.0934050772864460 , 0.0934050772864460 , 0.5418275218199933 , 0.0934050772864460 , 0.0000665877792863 );
            ir->AddPentPoint( 186 , 0.0934050772864460 , 0.5418275218199933 , 0.1779572463206686 , 0.0934050772864460 , 0.0000665877792863 );
            ir->AddPentPoint( 187 , 0.5418275218199933 , 0.0934050772864460 , 0.1779572463206686 , 0.0934050772864460 , 0.0000665877792863 );
            ir->AddPentPoint( 188 , 0.0934050772864460 , 0.0934050772864460 , 0.1779572463206686 , 0.0934050772864460 , 0.0000665877792863 );
            ir->AddPentPoint( 189 , 0.1779572463206686 , 0.5418275218199933 , 0.0934050772864460 , 0.0934050772864460 , 0.0000665877792863 );
            ir->AddPentPoint( 190 , 0.0934050772864460 , 0.5418275218199933 , 0.0934050772864460 , 0.0934050772864460 , 0.0000665877792863 );
            ir->AddPentPoint( 191 , 0.5418275218199933 , 0.1779572463206686 , 0.0934050772864460 , 0.0934050772864460 , 0.0000665877792863 );
            ir->AddPentPoint( 192 , 0.0934050772864460 , 0.1779572463206686 , 0.0934050772864460 , 0.0934050772864460 , 0.0000665877792863 );
            ir->AddPentPoint( 193 , 0.5418275218199933 , 0.0934050772864460 , 0.0934050772864460 , 0.0934050772864460 , 0.0000665877792863 );
            ir->AddPentPoint( 194 , 0.1779572463206686 , 0.0934050772864460 , 0.0934050772864460 , 0.0934050772864460 , 0.0000665877792863 );
            ir->AddPentPoint( 195 , 0.1362299703081878 , 0.1362299703081878 , 0.0220798971966119 , 0.5692301918788247 , 0.0000530688483067 );
            ir->AddPentPoint( 196 , 0.1362299703081878 , 0.0220798971966119 , 0.1362299703081878 , 0.5692301918788247 , 0.0000530688483067 );
            ir->AddPentPoint( 197 , 0.0220798971966119 , 0.1362299703081878 , 0.1362299703081878 , 0.5692301918788247 , 0.0000530688483067 );
            ir->AddPentPoint( 198 , 0.1362299703081878 , 0.1362299703081878 , 0.1362299703081878 , 0.5692301918788247 , 0.0000530688483067 );
            ir->AddPentPoint( 199 , 0.1362299703081878 , 0.1362299703081878 , 0.5692301918788247 , 0.0220798971966119 , 0.0000530688483067 );
            ir->AddPentPoint( 200 , 0.1362299703081878 , 0.5692301918788247 , 0.1362299703081878 , 0.0220798971966119 , 0.0000530688483067 );
            ir->AddPentPoint( 201 , 0.5692301918788247 , 0.1362299703081878 , 0.1362299703081878 , 0.0220798971966119 , 0.0000530688483067 );
            ir->AddPentPoint( 202 , 0.1362299703081878 , 0.1362299703081878 , 0.1362299703081878 , 0.0220798971966119 , 0.0000530688483067 );
            ir->AddPentPoint( 203 , 0.1362299703081878 , 0.0220798971966119 , 0.5692301918788247 , 0.1362299703081878 , 0.0000530688483067 );
            ir->AddPentPoint( 204 , 0.0220798971966119 , 0.1362299703081878 , 0.5692301918788247 , 0.1362299703081878 , 0.0000530688483067 );
            ir->AddPentPoint( 205 , 0.1362299703081878 , 0.1362299703081878 , 0.5692301918788247 , 0.1362299703081878 , 0.0000530688483067 );
            ir->AddPentPoint( 206 , 0.1362299703081878 , 0.5692301918788247 , 0.0220798971966119 , 0.1362299703081878 , 0.0000530688483067 );
            ir->AddPentPoint( 207 , 0.5692301918788247 , 0.1362299703081878 , 0.0220798971966119 , 0.1362299703081878 , 0.0000530688483067 );
            ir->AddPentPoint( 208 , 0.1362299703081878 , 0.1362299703081878 , 0.0220798971966119 , 0.1362299703081878 , 0.0000530688483067 );
            ir->AddPentPoint( 209 , 0.0220798971966119 , 0.5692301918788247 , 0.1362299703081878 , 0.1362299703081878 , 0.0000530688483067 );
            ir->AddPentPoint( 210 , 0.1362299703081878 , 0.5692301918788247 , 0.1362299703081878 , 0.1362299703081878 , 0.0000530688483067 );
            ir->AddPentPoint( 211 , 0.5692301918788247 , 0.0220798971966119 , 0.1362299703081878 , 0.1362299703081878 , 0.0000530688483067 );
            ir->AddPentPoint( 212 , 0.1362299703081878 , 0.0220798971966119 , 0.1362299703081878 , 0.1362299703081878 , 0.0000530688483067 );
            ir->AddPentPoint( 213 , 0.5692301918788247 , 0.1362299703081878 , 0.1362299703081878 , 0.1362299703081878 , 0.0000530688483067 );
            ir->AddPentPoint( 214 , 0.0220798971966119 , 0.1362299703081878 , 0.1362299703081878 , 0.1362299703081878 , 0.0000530688483067 );
            ir->AddPentPoint( 215 , 0.1732979951188451 , 0.1732979951188451 , 0.0841021374124987 , 0.3960038772309660 , 0.0001527784571035 );
            ir->AddPentPoint( 216 , 0.1732979951188451 , 0.0841021374124987 , 0.1732979951188451 , 0.3960038772309660 , 0.0001527784571035 );
            ir->AddPentPoint( 217 , 0.0841021374124987 , 0.1732979951188451 , 0.1732979951188451 , 0.3960038772309660 , 0.0001527784571035 );
            ir->AddPentPoint( 218 , 0.1732979951188451 , 0.1732979951188451 , 0.1732979951188451 , 0.3960038772309660 , 0.0001527784571035 );
            ir->AddPentPoint( 219 , 0.1732979951188451 , 0.1732979951188451 , 0.3960038772309660 , 0.0841021374124987 , 0.0001527784571035 );
            ir->AddPentPoint( 220 , 0.1732979951188451 , 0.3960038772309660 , 0.1732979951188451 , 0.0841021374124987 , 0.0001527784571035 );
            ir->AddPentPoint( 221 , 0.3960038772309660 , 0.1732979951188451 , 0.1732979951188451 , 0.0841021374124987 , 0.0001527784571035 );
            ir->AddPentPoint( 222 , 0.1732979951188451 , 0.1732979951188451 , 0.1732979951188451 , 0.0841021374124987 , 0.0001527784571035 );
            ir->AddPentPoint( 223 , 0.1732979951188451 , 0.0841021374124987 , 0.3960038772309660 , 0.1732979951188451 , 0.0001527784571035 );
            ir->AddPentPoint( 224 , 0.0841021374124987 , 0.1732979951188451 , 0.3960038772309660 , 0.1732979951188451 , 0.0001527784571035 );
            ir->AddPentPoint( 225 , 0.1732979951188451 , 0.1732979951188451 , 0.3960038772309660 , 0.1732979951188451 , 0.0001527784571035 );
            ir->AddPentPoint( 226 , 0.1732979951188451 , 0.3960038772309660 , 0.0841021374124987 , 0.1732979951188451 , 0.0001527784571035 );
            ir->AddPentPoint( 227 , 0.3960038772309660 , 0.1732979951188451 , 0.0841021374124987 , 0.1732979951188451 , 0.0001527784571035 );
            ir->AddPentPoint( 228 , 0.1732979951188451 , 0.1732979951188451 , 0.0841021374124987 , 0.1732979951188451 , 0.0001527784571035 );
            ir->AddPentPoint( 229 , 0.0841021374124987 , 0.3960038772309660 , 0.1732979951188451 , 0.1732979951188451 , 0.0001527784571035 );
            ir->AddPentPoint( 230 , 0.1732979951188451 , 0.3960038772309660 , 0.1732979951188451 , 0.1732979951188451 , 0.0001527784571035 );
            ir->AddPentPoint( 231 , 0.3960038772309660 , 0.0841021374124987 , 0.1732979951188451 , 0.1732979951188451 , 0.0001527784571035 );
            ir->AddPentPoint( 232 , 0.1732979951188451 , 0.0841021374124987 , 0.1732979951188451 , 0.1732979951188451 , 0.0001527784571035 );
            ir->AddPentPoint( 233 , 0.3960038772309660 , 0.1732979951188451 , 0.1732979951188451 , 0.1732979951188451 , 0.0001527784571035 );
            ir->AddPentPoint( 234 , 0.0841021374124987 , 0.1732979951188451 , 0.1732979951188451 , 0.1732979951188451 , 0.0001527784571035 );
            ir->AddPentPoint( 235 , 0.2683851232214430 , 0.2683851232214430 , 0.0510056350893379 , 0.1438389952463330 , 0.0001224244766359 );
            ir->AddPentPoint( 236 , 0.2683851232214430 , 0.0510056350893379 , 0.2683851232214430 , 0.1438389952463330 , 0.0001224244766359 );
            ir->AddPentPoint( 237 , 0.0510056350893379 , 0.2683851232214430 , 0.2683851232214430 , 0.1438389952463330 , 0.0001224244766359 );
            ir->AddPentPoint( 238 , 0.2683851232214430 , 0.2683851232214430 , 0.2683851232214430 , 0.1438389952463330 , 0.0001224244766359 );
            ir->AddPentPoint( 239 , 0.2683851232214430 , 0.2683851232214430 , 0.1438389952463330 , 0.0510056350893379 , 0.0001224244766359 );
            ir->AddPentPoint( 240 , 0.2683851232214430 , 0.1438389952463330 , 0.2683851232214430 , 0.0510056350893379 , 0.0001224244766359 );
            ir->AddPentPoint( 241 , 0.1438389952463330 , 0.2683851232214430 , 0.2683851232214430 , 0.0510056350893379 , 0.0001224244766359 );
            ir->AddPentPoint( 242 , 0.2683851232214430 , 0.2683851232214430 , 0.2683851232214430 , 0.0510056350893379 , 0.0001224244766359 );
            ir->AddPentPoint( 243 , 0.2683851232214430 , 0.0510056350893379 , 0.1438389952463330 , 0.2683851232214430 , 0.0001224244766359 );
            ir->AddPentPoint( 244 , 0.0510056350893379 , 0.2683851232214430 , 0.1438389952463330 , 0.2683851232214430 , 0.0001224244766359 );
            ir->AddPentPoint( 245 , 0.2683851232214430 , 0.2683851232214430 , 0.1438389952463330 , 0.2683851232214430 , 0.0001224244766359 );
            ir->AddPentPoint( 246 , 0.2683851232214430 , 0.1438389952463330 , 0.0510056350893379 , 0.2683851232214430 , 0.0001224244766359 );
            ir->AddPentPoint( 247 , 0.1438389952463330 , 0.2683851232214430 , 0.0510056350893379 , 0.2683851232214430 , 0.0001224244766359 );
            ir->AddPentPoint( 248 , 0.2683851232214430 , 0.2683851232214430 , 0.0510056350893379 , 0.2683851232214430 , 0.0001224244766359 );
            ir->AddPentPoint( 249 , 0.0510056350893379 , 0.1438389952463330 , 0.2683851232214430 , 0.2683851232214430 , 0.0001224244766359 );
            ir->AddPentPoint( 250 , 0.2683851232214430 , 0.1438389952463330 , 0.2683851232214430 , 0.2683851232214430 , 0.0001224244766359 );
            ir->AddPentPoint( 251 , 0.1438389952463330 , 0.0510056350893379 , 0.2683851232214430 , 0.2683851232214430 , 0.0001224244766359 );
            ir->AddPentPoint( 252 , 0.2683851232214430 , 0.0510056350893379 , 0.2683851232214430 , 0.2683851232214430 , 0.0001224244766359 );
            ir->AddPentPoint( 253 , 0.1438389952463330 , 0.2683851232214430 , 0.2683851232214430 , 0.2683851232214430 , 0.0001224244766359 );
            ir->AddPentPoint( 254 , 0.0510056350893379 , 0.2683851232214430 , 0.2683851232214430 , 0.2683851232214430 , 0.0001224244766359 );
            ir->AddPentPoint( 255 , 0.3142975959592881 , 0.3142975959592881 , 0.0098158564884535 , 0.0472913556336823 , 0.0000374546635434 );
            ir->AddPentPoint( 256 , 0.3142975959592881 , 0.0098158564884535 , 0.3142975959592881 , 0.0472913556336823 , 0.0000374546635434 );
            ir->AddPentPoint( 257 , 0.0098158564884535 , 0.3142975959592881 , 0.3142975959592881 , 0.0472913556336823 , 0.0000374546635434 );
            ir->AddPentPoint( 258 , 0.3142975959592881 , 0.3142975959592881 , 0.3142975959592881 , 0.0472913556336823 , 0.0000374546635434 );
            ir->AddPentPoint( 259 , 0.3142975959592881 , 0.3142975959592881 , 0.0472913556336823 , 0.0098158564884535 , 0.0000374546635434 );
            ir->AddPentPoint( 260 , 0.3142975959592881 , 0.0472913556336823 , 0.3142975959592881 , 0.0098158564884535 , 0.0000374546635434 );
            ir->AddPentPoint( 261 , 0.0472913556336823 , 0.3142975959592881 , 0.3142975959592881 , 0.0098158564884535 , 0.0000374546635434 );
            ir->AddPentPoint( 262 , 0.3142975959592881 , 0.3142975959592881 , 0.3142975959592881 , 0.0098158564884535 , 0.0000374546635434 );
            ir->AddPentPoint( 263 , 0.3142975959592881 , 0.0098158564884535 , 0.0472913556336823 , 0.3142975959592881 , 0.0000374546635434 );
            ir->AddPentPoint( 264 , 0.0098158564884535 , 0.3142975959592881 , 0.0472913556336823 , 0.3142975959592881 , 0.0000374546635434 );
            ir->AddPentPoint( 265 , 0.3142975959592881 , 0.3142975959592881 , 0.0472913556336823 , 0.3142975959592881 , 0.0000374546635434 );
            ir->AddPentPoint( 266 , 0.3142975959592881 , 0.0472913556336823 , 0.0098158564884535 , 0.3142975959592881 , 0.0000374546635434 );
            ir->AddPentPoint( 267 , 0.0472913556336823 , 0.3142975959592881 , 0.0098158564884535 , 0.3142975959592881 , 0.0000374546635434 );
            ir->AddPentPoint( 268 , 0.3142975959592881 , 0.3142975959592881 , 0.0098158564884535 , 0.3142975959592881 , 0.0000374546635434 );
            ir->AddPentPoint( 269 , 0.0098158564884535 , 0.0472913556336823 , 0.3142975959592881 , 0.3142975959592881 , 0.0000374546635434 );
            ir->AddPentPoint( 270 , 0.3142975959592881 , 0.0472913556336823 , 0.3142975959592881 , 0.3142975959592881 , 0.0000374546635434 );
            ir->AddPentPoint( 271 , 0.0472913556336823 , 0.0098158564884535 , 0.3142975959592881 , 0.3142975959592881 , 0.0000374546635434 );
            ir->AddPentPoint( 272 , 0.3142975959592881 , 0.0098158564884535 , 0.3142975959592881 , 0.3142975959592881 , 0.0000374546635434 );
            ir->AddPentPoint( 273 , 0.0472913556336823 , 0.3142975959592881 , 0.3142975959592881 , 0.3142975959592881 , 0.0000374546635434 );
            ir->AddPentPoint( 274 , 0.0098158564884535 , 0.3142975959592881 , 0.3142975959592881 , 0.3142975959592881 , 0.0000374546635434 );
            ir->AddPentPoint( 275 , 0.0089460554813759 , 0.4509983952392004 , 0.4509983952392004 , 0.0801110985588474 , 0.0000110932015523 );
            ir->AddPentPoint( 276 , 0.4509983952392004 , 0.0089460554813759 , 0.4509983952392004 , 0.0801110985588474 , 0.0000110932015523 );
            ir->AddPentPoint( 277 , 0.0089460554813759 , 0.0089460554813759 , 0.4509983952392004 , 0.0801110985588474 , 0.0000110932015523 );
            ir->AddPentPoint( 278 , 0.4509983952392004 , 0.4509983952392004 , 0.0089460554813759 , 0.0801110985588474 , 0.0000110932015523 );
            ir->AddPentPoint( 279 , 0.0089460554813759 , 0.4509983952392004 , 0.0089460554813759 , 0.0801110985588474 , 0.0000110932015523 );
            ir->AddPentPoint( 280 , 0.4509983952392004 , 0.0089460554813759 , 0.0089460554813759 , 0.0801110985588474 , 0.0000110932015523 );
            ir->AddPentPoint( 281 , 0.0089460554813759 , 0.4509983952392004 , 0.0801110985588474 , 0.4509983952392004 , 0.0000110932015523 );
            ir->AddPentPoint( 282 , 0.4509983952392004 , 0.0089460554813759 , 0.0801110985588474 , 0.4509983952392004 , 0.0000110932015523 );
            ir->AddPentPoint( 283 , 0.0089460554813759 , 0.0089460554813759 , 0.0801110985588474 , 0.4509983952392004 , 0.0000110932015523 );
            ir->AddPentPoint( 284 , 0.0089460554813759 , 0.0801110985588474 , 0.4509983952392004 , 0.4509983952392004 , 0.0000110932015523 );
            ir->AddPentPoint( 285 , 0.0801110985588474 , 0.0089460554813759 , 0.4509983952392004 , 0.4509983952392004 , 0.0000110932015523 );
            ir->AddPentPoint( 286 , 0.0089460554813759 , 0.0089460554813759 , 0.4509983952392004 , 0.4509983952392004 , 0.0000110932015523 );
            ir->AddPentPoint( 287 , 0.4509983952392004 , 0.0801110985588474 , 0.0089460554813759 , 0.4509983952392004 , 0.0000110932015523 );
            ir->AddPentPoint( 288 , 0.0089460554813759 , 0.0801110985588474 , 0.0089460554813759 , 0.4509983952392004 , 0.0000110932015523 );
            ir->AddPentPoint( 289 , 0.0801110985588474 , 0.4509983952392004 , 0.0089460554813759 , 0.4509983952392004 , 0.0000110932015523 );
            ir->AddPentPoint( 290 , 0.0089460554813759 , 0.4509983952392004 , 0.0089460554813759 , 0.4509983952392004 , 0.0000110932015523 );
            ir->AddPentPoint( 291 , 0.0801110985588474 , 0.0089460554813759 , 0.0089460554813759 , 0.4509983952392004 , 0.0000110932015523 );
            ir->AddPentPoint( 292 , 0.4509983952392004 , 0.0089460554813759 , 0.0089460554813759 , 0.4509983952392004 , 0.0000110932015523 );
            ir->AddPentPoint( 293 , 0.4509983952392004 , 0.4509983952392004 , 0.0801110985588474 , 0.0089460554813759 , 0.0000110932015523 );
            ir->AddPentPoint( 294 , 0.0089460554813759 , 0.4509983952392004 , 0.0801110985588474 , 0.0089460554813759 , 0.0000110932015523 );
            ir->AddPentPoint( 295 , 0.4509983952392004 , 0.0089460554813759 , 0.0801110985588474 , 0.0089460554813759 , 0.0000110932015523 );
            ir->AddPentPoint( 296 , 0.4509983952392004 , 0.0801110985588474 , 0.4509983952392004 , 0.0089460554813759 , 0.0000110932015523 );
            ir->AddPentPoint( 297 , 0.0089460554813759 , 0.0801110985588474 , 0.4509983952392004 , 0.0089460554813759 , 0.0000110932015523 );
            ir->AddPentPoint( 298 , 0.0801110985588474 , 0.4509983952392004 , 0.4509983952392004 , 0.0089460554813759 , 0.0000110932015523 );
            ir->AddPentPoint( 299 , 0.0089460554813759 , 0.4509983952392004 , 0.4509983952392004 , 0.0089460554813759 , 0.0000110932015523 );
            ir->AddPentPoint( 300 , 0.0801110985588474 , 0.0089460554813759 , 0.4509983952392004 , 0.0089460554813759 , 0.0000110932015523 );
            ir->AddPentPoint( 301 , 0.4509983952392004 , 0.0089460554813759 , 0.4509983952392004 , 0.0089460554813759 , 0.0000110932015523 );
            ir->AddPentPoint( 302 , 0.4509983952392004 , 0.0801110985588474 , 0.0089460554813759 , 0.0089460554813759 , 0.0000110932015523 );
            ir->AddPentPoint( 303 , 0.0801110985588474 , 0.4509983952392004 , 0.0089460554813759 , 0.0089460554813759 , 0.0000110932015523 );
            ir->AddPentPoint( 304 , 0.4509983952392004 , 0.4509983952392004 , 0.0089460554813759 , 0.0089460554813759 , 0.0000110932015523 );
            ir->AddPentPoint( 305 , 0.0099797990160544 , 0.3886935638562216 , 0.3886935638562216 , 0.2026532742554480 , 0.0000140359239615 );
            ir->AddPentPoint( 306 , 0.3886935638562216 , 0.0099797990160544 , 0.3886935638562216 , 0.2026532742554480 , 0.0000140359239615 );
            ir->AddPentPoint( 307 , 0.0099797990160544 , 0.0099797990160544 , 0.3886935638562216 , 0.2026532742554480 , 0.0000140359239615 );
            ir->AddPentPoint( 308 , 0.3886935638562216 , 0.3886935638562216 , 0.0099797990160544 , 0.2026532742554480 , 0.0000140359239615 );
            ir->AddPentPoint( 309 , 0.0099797990160544 , 0.3886935638562216 , 0.0099797990160544 , 0.2026532742554480 , 0.0000140359239615 );
            ir->AddPentPoint( 310 , 0.3886935638562216 , 0.0099797990160544 , 0.0099797990160544 , 0.2026532742554480 , 0.0000140359239615 );
            ir->AddPentPoint( 311 , 0.0099797990160544 , 0.3886935638562216 , 0.2026532742554480 , 0.3886935638562216 , 0.0000140359239615 );
            ir->AddPentPoint( 312 , 0.3886935638562216 , 0.0099797990160544 , 0.2026532742554480 , 0.3886935638562216 , 0.0000140359239615 );
            ir->AddPentPoint( 313 , 0.0099797990160544 , 0.0099797990160544 , 0.2026532742554480 , 0.3886935638562216 , 0.0000140359239615 );
            ir->AddPentPoint( 314 , 0.0099797990160544 , 0.2026532742554480 , 0.3886935638562216 , 0.3886935638562216 , 0.0000140359239615 );
            ir->AddPentPoint( 315 , 0.2026532742554480 , 0.0099797990160544 , 0.3886935638562216 , 0.3886935638562216 , 0.0000140359239615 );
            ir->AddPentPoint( 316 , 0.0099797990160544 , 0.0099797990160544 , 0.3886935638562216 , 0.3886935638562216 , 0.0000140359239615 );
            ir->AddPentPoint( 317 , 0.3886935638562216 , 0.2026532742554480 , 0.0099797990160544 , 0.3886935638562216 , 0.0000140359239615 );
            ir->AddPentPoint( 318 , 0.0099797990160544 , 0.2026532742554480 , 0.0099797990160544 , 0.3886935638562216 , 0.0000140359239615 );
            ir->AddPentPoint( 319 , 0.2026532742554480 , 0.3886935638562216 , 0.0099797990160544 , 0.3886935638562216 , 0.0000140359239615 );
            ir->AddPentPoint( 320 , 0.0099797990160544 , 0.3886935638562216 , 0.0099797990160544 , 0.3886935638562216 , 0.0000140359239615 );
            ir->AddPentPoint( 321 , 0.2026532742554480 , 0.0099797990160544 , 0.0099797990160544 , 0.3886935638562216 , 0.0000140359239615 );
            ir->AddPentPoint( 322 , 0.3886935638562216 , 0.0099797990160544 , 0.0099797990160544 , 0.3886935638562216 , 0.0000140359239615 );
            ir->AddPentPoint( 323 , 0.3886935638562216 , 0.3886935638562216 , 0.2026532742554480 , 0.0099797990160544 , 0.0000140359239615 );
            ir->AddPentPoint( 324 , 0.0099797990160544 , 0.3886935638562216 , 0.2026532742554480 , 0.0099797990160544 , 0.0000140359239615 );
            ir->AddPentPoint( 325 , 0.3886935638562216 , 0.0099797990160544 , 0.2026532742554480 , 0.0099797990160544 , 0.0000140359239615 );
            ir->AddPentPoint( 326 , 0.3886935638562216 , 0.2026532742554480 , 0.3886935638562216 , 0.0099797990160544 , 0.0000140359239615 );
            ir->AddPentPoint( 327 , 0.0099797990160544 , 0.2026532742554480 , 0.3886935638562216 , 0.0099797990160544 , 0.0000140359239615 );
            ir->AddPentPoint( 328 , 0.2026532742554480 , 0.3886935638562216 , 0.3886935638562216 , 0.0099797990160544 , 0.0000140359239615 );
            ir->AddPentPoint( 329 , 0.0099797990160544 , 0.3886935638562216 , 0.3886935638562216 , 0.0099797990160544 , 0.0000140359239615 );
            ir->AddPentPoint( 330 , 0.2026532742554480 , 0.0099797990160544 , 0.3886935638562216 , 0.0099797990160544 , 0.0000140359239615 );
            ir->AddPentPoint( 331 , 0.3886935638562216 , 0.0099797990160544 , 0.3886935638562216 , 0.0099797990160544 , 0.0000140359239615 );
            ir->AddPentPoint( 332 , 0.3886935638562216 , 0.2026532742554480 , 0.0099797990160544 , 0.0099797990160544 , 0.0000140359239615 );
            ir->AddPentPoint( 333 , 0.2026532742554480 , 0.3886935638562216 , 0.0099797990160544 , 0.0099797990160544 , 0.0000140359239615 );
            ir->AddPentPoint( 334 , 0.3886935638562216 , 0.3886935638562216 , 0.0099797990160544 , 0.0099797990160544 , 0.0000140359239615 );
            ir->AddPentPoint( 335 , 0.0124135340163404 , 0.2456525404904195 , 0.2456525404904195 , 0.4838678509864802 , 0.0000193340058483 );
            ir->AddPentPoint( 336 , 0.2456525404904195 , 0.0124135340163404 , 0.2456525404904195 , 0.4838678509864802 , 0.0000193340058483 );
            ir->AddPentPoint( 337 , 0.0124135340163404 , 0.0124135340163404 , 0.2456525404904195 , 0.4838678509864802 , 0.0000193340058483 );
            ir->AddPentPoint( 338 , 0.2456525404904195 , 0.2456525404904195 , 0.0124135340163404 , 0.4838678509864802 , 0.0000193340058483 );
            ir->AddPentPoint( 339 , 0.0124135340163404 , 0.2456525404904195 , 0.0124135340163404 , 0.4838678509864802 , 0.0000193340058483 );
            ir->AddPentPoint( 340 , 0.2456525404904195 , 0.0124135340163404 , 0.0124135340163404 , 0.4838678509864802 , 0.0000193340058483 );
            ir->AddPentPoint( 341 , 0.0124135340163404 , 0.2456525404904195 , 0.4838678509864802 , 0.2456525404904195 , 0.0000193340058483 );
            ir->AddPentPoint( 342 , 0.2456525404904195 , 0.0124135340163404 , 0.4838678509864802 , 0.2456525404904195 , 0.0000193340058483 );
            ir->AddPentPoint( 343 , 0.0124135340163404 , 0.0124135340163404 , 0.4838678509864802 , 0.2456525404904195 , 0.0000193340058483 );
            ir->AddPentPoint( 344 , 0.0124135340163404 , 0.4838678509864802 , 0.2456525404904195 , 0.2456525404904195 , 0.0000193340058483 );
            ir->AddPentPoint( 345 , 0.4838678509864802 , 0.0124135340163404 , 0.2456525404904195 , 0.2456525404904195 , 0.0000193340058483 );
            ir->AddPentPoint( 346 , 0.0124135340163404 , 0.0124135340163404 , 0.2456525404904195 , 0.2456525404904195 , 0.0000193340058483 );
            ir->AddPentPoint( 347 , 0.2456525404904195 , 0.4838678509864802 , 0.0124135340163404 , 0.2456525404904195 , 0.0000193340058483 );
            ir->AddPentPoint( 348 , 0.0124135340163404 , 0.4838678509864802 , 0.0124135340163404 , 0.2456525404904195 , 0.0000193340058483 );
            ir->AddPentPoint( 349 , 0.4838678509864802 , 0.2456525404904195 , 0.0124135340163404 , 0.2456525404904195 , 0.0000193340058483 );
            ir->AddPentPoint( 350 , 0.0124135340163404 , 0.2456525404904195 , 0.0124135340163404 , 0.2456525404904195 , 0.0000193340058483 );
            ir->AddPentPoint( 351 , 0.4838678509864802 , 0.0124135340163404 , 0.0124135340163404 , 0.2456525404904195 , 0.0000193340058483 );
            ir->AddPentPoint( 352 , 0.2456525404904195 , 0.0124135340163404 , 0.0124135340163404 , 0.2456525404904195 , 0.0000193340058483 );
            ir->AddPentPoint( 353 , 0.2456525404904195 , 0.2456525404904195 , 0.4838678509864802 , 0.0124135340163404 , 0.0000193340058483 );
            ir->AddPentPoint( 354 , 0.0124135340163404 , 0.2456525404904195 , 0.4838678509864802 , 0.0124135340163404 , 0.0000193340058483 );
            ir->AddPentPoint( 355 , 0.2456525404904195 , 0.0124135340163404 , 0.4838678509864802 , 0.0124135340163404 , 0.0000193340058483 );
            ir->AddPentPoint( 356 , 0.2456525404904195 , 0.4838678509864802 , 0.2456525404904195 , 0.0124135340163404 , 0.0000193340058483 );
            ir->AddPentPoint( 357 , 0.0124135340163404 , 0.4838678509864802 , 0.2456525404904195 , 0.0124135340163404 , 0.0000193340058483 );
            ir->AddPentPoint( 358 , 0.4838678509864802 , 0.2456525404904195 , 0.2456525404904195 , 0.0124135340163404 , 0.0000193340058483 );
            ir->AddPentPoint( 359 , 0.0124135340163404 , 0.2456525404904195 , 0.2456525404904195 , 0.0124135340163404 , 0.0000193340058483 );
            ir->AddPentPoint( 360 , 0.4838678509864802 , 0.0124135340163404 , 0.2456525404904195 , 0.0124135340163404 , 0.0000193340058483 );
            ir->AddPentPoint( 361 , 0.2456525404904195 , 0.0124135340163404 , 0.2456525404904195 , 0.0124135340163404 , 0.0000193340058483 );
            ir->AddPentPoint( 362 , 0.2456525404904195 , 0.4838678509864802 , 0.0124135340163404 , 0.0124135340163404 , 0.0000193340058483 );
            ir->AddPentPoint( 363 , 0.4838678509864802 , 0.2456525404904195 , 0.0124135340163404 , 0.0124135340163404 , 0.0000193340058483 );
            ir->AddPentPoint( 364 , 0.2456525404904195 , 0.2456525404904195 , 0.0124135340163404 , 0.0124135340163404 , 0.0000193340058483 );
            ir->AddPentPoint( 365 , 0.0145879385947127 , 0.0815578081256928 , 0.0815578081256928 , 0.8077085065591889 , 0.0000126033861441 );
            ir->AddPentPoint( 366 , 0.0815578081256928 , 0.0145879385947127 , 0.0815578081256928 , 0.8077085065591889 , 0.0000126033861441 );
            ir->AddPentPoint( 367 , 0.0145879385947127 , 0.0145879385947127 , 0.0815578081256928 , 0.8077085065591889 , 0.0000126033861441 );
            ir->AddPentPoint( 368 , 0.0815578081256928 , 0.0815578081256928 , 0.0145879385947127 , 0.8077085065591889 , 0.0000126033861441 );
            ir->AddPentPoint( 369 , 0.0145879385947127 , 0.0815578081256928 , 0.0145879385947127 , 0.8077085065591889 , 0.0000126033861441 );
            ir->AddPentPoint( 370 , 0.0815578081256928 , 0.0145879385947127 , 0.0145879385947127 , 0.8077085065591889 , 0.0000126033861441 );
            ir->AddPentPoint( 371 , 0.0145879385947127 , 0.0815578081256928 , 0.8077085065591889 , 0.0815578081256928 , 0.0000126033861441 );
            ir->AddPentPoint( 372 , 0.0815578081256928 , 0.0145879385947127 , 0.8077085065591889 , 0.0815578081256928 , 0.0000126033861441 );
            ir->AddPentPoint( 373 , 0.0145879385947127 , 0.0145879385947127 , 0.8077085065591889 , 0.0815578081256928 , 0.0000126033861441 );
            ir->AddPentPoint( 374 , 0.0145879385947127 , 0.8077085065591889 , 0.0815578081256928 , 0.0815578081256928 , 0.0000126033861441 );
            ir->AddPentPoint( 375 , 0.8077085065591889 , 0.0145879385947127 , 0.0815578081256928 , 0.0815578081256928 , 0.0000126033861441 );
            ir->AddPentPoint( 376 , 0.0145879385947127 , 0.0145879385947127 , 0.0815578081256928 , 0.0815578081256928 , 0.0000126033861441 );
            ir->AddPentPoint( 377 , 0.0815578081256928 , 0.8077085065591889 , 0.0145879385947127 , 0.0815578081256928 , 0.0000126033861441 );
            ir->AddPentPoint( 378 , 0.0145879385947127 , 0.8077085065591889 , 0.0145879385947127 , 0.0815578081256928 , 0.0000126033861441 );
            ir->AddPentPoint( 379 , 0.8077085065591889 , 0.0815578081256928 , 0.0145879385947127 , 0.0815578081256928 , 0.0000126033861441 );
            ir->AddPentPoint( 380 , 0.0145879385947127 , 0.0815578081256928 , 0.0145879385947127 , 0.0815578081256928 , 0.0000126033861441 );
            ir->AddPentPoint( 381 , 0.8077085065591889 , 0.0145879385947127 , 0.0145879385947127 , 0.0815578081256928 , 0.0000126033861441 );
            ir->AddPentPoint( 382 , 0.0815578081256928 , 0.0145879385947127 , 0.0145879385947127 , 0.0815578081256928 , 0.0000126033861441 );
            ir->AddPentPoint( 383 , 0.0815578081256928 , 0.0815578081256928 , 0.8077085065591889 , 0.0145879385947127 , 0.0000126033861441 );
            ir->AddPentPoint( 384 , 0.0145879385947127 , 0.0815578081256928 , 0.8077085065591889 , 0.0145879385947127 , 0.0000126033861441 );
            ir->AddPentPoint( 385 , 0.0815578081256928 , 0.0145879385947127 , 0.8077085065591889 , 0.0145879385947127 , 0.0000126033861441 );
            ir->AddPentPoint( 386 , 0.0815578081256928 , 0.8077085065591889 , 0.0815578081256928 , 0.0145879385947127 , 0.0000126033861441 );
            ir->AddPentPoint( 387 , 0.0145879385947127 , 0.8077085065591889 , 0.0815578081256928 , 0.0145879385947127 , 0.0000126033861441 );
            ir->AddPentPoint( 388 , 0.8077085065591889 , 0.0815578081256928 , 0.0815578081256928 , 0.0145879385947127 , 0.0000126033861441 );
            ir->AddPentPoint( 389 , 0.0145879385947127 , 0.0815578081256928 , 0.0815578081256928 , 0.0145879385947127 , 0.0000126033861441 );
            ir->AddPentPoint( 390 , 0.8077085065591889 , 0.0145879385947127 , 0.0815578081256928 , 0.0145879385947127 , 0.0000126033861441 );
            ir->AddPentPoint( 391 , 0.0815578081256928 , 0.0145879385947127 , 0.0815578081256928 , 0.0145879385947127 , 0.0000126033861441 );
            ir->AddPentPoint( 392 , 0.0815578081256928 , 0.8077085065591889 , 0.0145879385947127 , 0.0145879385947127 , 0.0000126033861441 );
            ir->AddPentPoint( 393 , 0.8077085065591889 , 0.0815578081256928 , 0.0145879385947127 , 0.0145879385947127 , 0.0000126033861441 );
            ir->AddPentPoint( 394 , 0.0815578081256928 , 0.0815578081256928 , 0.0145879385947127 , 0.0145879385947127 , 0.0000126033861441 );
            ir->AddPentPoint( 395 , 0.0146849968005364 , 0.1527973568303998 , 0.1527973568303998 , 0.6650352927381276 , 0.0000164938495139 );
            ir->AddPentPoint( 396 , 0.1527973568303998 , 0.0146849968005364 , 0.1527973568303998 , 0.6650352927381276 , 0.0000164938495139 );
            ir->AddPentPoint( 397 , 0.0146849968005364 , 0.0146849968005364 , 0.1527973568303998 , 0.6650352927381276 , 0.0000164938495139 );
            ir->AddPentPoint( 398 , 0.1527973568303998 , 0.1527973568303998 , 0.0146849968005364 , 0.6650352927381276 , 0.0000164938495139 );
            ir->AddPentPoint( 399 , 0.0146849968005364 , 0.1527973568303998 , 0.0146849968005364 , 0.6650352927381276 , 0.0000164938495139 );
            ir->AddPentPoint( 400 , 0.1527973568303998 , 0.0146849968005364 , 0.0146849968005364 , 0.6650352927381276 , 0.0000164938495139 );
            ir->AddPentPoint( 401 , 0.0146849968005364 , 0.1527973568303998 , 0.6650352927381276 , 0.1527973568303998 , 0.0000164938495139 );
            ir->AddPentPoint( 402 , 0.1527973568303998 , 0.0146849968005364 , 0.6650352927381276 , 0.1527973568303998 , 0.0000164938495139 );
            ir->AddPentPoint( 403 , 0.0146849968005364 , 0.0146849968005364 , 0.6650352927381276 , 0.1527973568303998 , 0.0000164938495139 );
            ir->AddPentPoint( 404 , 0.0146849968005364 , 0.6650352927381276 , 0.1527973568303998 , 0.1527973568303998 , 0.0000164938495139 );
            ir->AddPentPoint( 405 , 0.6650352927381276 , 0.0146849968005364 , 0.1527973568303998 , 0.1527973568303998 , 0.0000164938495139 );
            ir->AddPentPoint( 406 , 0.0146849968005364 , 0.0146849968005364 , 0.1527973568303998 , 0.1527973568303998 , 0.0000164938495139 );
            ir->AddPentPoint( 407 , 0.1527973568303998 , 0.6650352927381276 , 0.0146849968005364 , 0.1527973568303998 , 0.0000164938495139 );
            ir->AddPentPoint( 408 , 0.0146849968005364 , 0.6650352927381276 , 0.0146849968005364 , 0.1527973568303998 , 0.0000164938495139 );
            ir->AddPentPoint( 409 , 0.6650352927381276 , 0.1527973568303998 , 0.0146849968005364 , 0.1527973568303998 , 0.0000164938495139 );
            ir->AddPentPoint( 410 , 0.0146849968005364 , 0.1527973568303998 , 0.0146849968005364 , 0.1527973568303998 , 0.0000164938495139 );
            ir->AddPentPoint( 411 , 0.6650352927381276 , 0.0146849968005364 , 0.0146849968005364 , 0.1527973568303998 , 0.0000164938495139 );
            ir->AddPentPoint( 412 , 0.1527973568303998 , 0.0146849968005364 , 0.0146849968005364 , 0.1527973568303998 , 0.0000164938495139 );
            ir->AddPentPoint( 413 , 0.1527973568303998 , 0.1527973568303998 , 0.6650352927381276 , 0.0146849968005364 , 0.0000164938495139 );
            ir->AddPentPoint( 414 , 0.0146849968005364 , 0.1527973568303998 , 0.6650352927381276 , 0.0146849968005364 , 0.0000164938495139 );
            ir->AddPentPoint( 415 , 0.1527973568303998 , 0.0146849968005364 , 0.6650352927381276 , 0.0146849968005364 , 0.0000164938495139 );
            ir->AddPentPoint( 416 , 0.1527973568303998 , 0.6650352927381276 , 0.1527973568303998 , 0.0146849968005364 , 0.0000164938495139 );
            ir->AddPentPoint( 417 , 0.0146849968005364 , 0.6650352927381276 , 0.1527973568303998 , 0.0146849968005364 , 0.0000164938495139 );
            ir->AddPentPoint( 418 , 0.6650352927381276 , 0.1527973568303998 , 0.1527973568303998 , 0.0146849968005364 , 0.0000164938495139 );
            ir->AddPentPoint( 419 , 0.0146849968005364 , 0.1527973568303998 , 0.1527973568303998 , 0.0146849968005364 , 0.0000164938495139 );
            ir->AddPentPoint( 420 , 0.6650352927381276 , 0.0146849968005364 , 0.1527973568303998 , 0.0146849968005364 , 0.0000164938495139 );
            ir->AddPentPoint( 421 , 0.1527973568303998 , 0.0146849968005364 , 0.1527973568303998 , 0.0146849968005364 , 0.0000164938495139 );
            ir->AddPentPoint( 422 , 0.1527973568303998 , 0.6650352927381276 , 0.0146849968005364 , 0.0146849968005364 , 0.0000164938495139 );
            ir->AddPentPoint( 423 , 0.6650352927381276 , 0.1527973568303998 , 0.0146849968005364 , 0.0146849968005364 , 0.0000164938495139 );
            ir->AddPentPoint( 424 , 0.1527973568303998 , 0.1527973568303998 , 0.0146849968005364 , 0.0146849968005364 , 0.0000164938495139 );
            ir->AddPentPoint( 425 , 0.0403754775891479 , 0.3921976592941051 , 0.3921976592941051 , 0.1348537262334940 , 0.0000325275043721 );
            ir->AddPentPoint( 426 , 0.3921976592941051 , 0.0403754775891479 , 0.3921976592941051 , 0.1348537262334940 , 0.0000325275043721 );
            ir->AddPentPoint( 427 , 0.0403754775891479 , 0.0403754775891479 , 0.3921976592941051 , 0.1348537262334940 , 0.0000325275043721 );
            ir->AddPentPoint( 428 , 0.3921976592941051 , 0.3921976592941051 , 0.0403754775891479 , 0.1348537262334940 , 0.0000325275043721 );
            ir->AddPentPoint( 429 , 0.0403754775891479 , 0.3921976592941051 , 0.0403754775891479 , 0.1348537262334940 , 0.0000325275043721 );
            ir->AddPentPoint( 430 , 0.3921976592941051 , 0.0403754775891479 , 0.0403754775891479 , 0.1348537262334940 , 0.0000325275043721 );
            ir->AddPentPoint( 431 , 0.0403754775891479 , 0.3921976592941051 , 0.1348537262334940 , 0.3921976592941051 , 0.0000325275043721 );
            ir->AddPentPoint( 432 , 0.3921976592941051 , 0.0403754775891479 , 0.1348537262334940 , 0.3921976592941051 , 0.0000325275043721 );
            ir->AddPentPoint( 433 , 0.0403754775891479 , 0.0403754775891479 , 0.1348537262334940 , 0.3921976592941051 , 0.0000325275043721 );
            ir->AddPentPoint( 434 , 0.0403754775891479 , 0.1348537262334940 , 0.3921976592941051 , 0.3921976592941051 , 0.0000325275043721 );
            ir->AddPentPoint( 435 , 0.1348537262334940 , 0.0403754775891479 , 0.3921976592941051 , 0.3921976592941051 , 0.0000325275043721 );
            ir->AddPentPoint( 436 , 0.0403754775891479 , 0.0403754775891479 , 0.3921976592941051 , 0.3921976592941051 , 0.0000325275043721 );
            ir->AddPentPoint( 437 , 0.3921976592941051 , 0.1348537262334940 , 0.0403754775891479 , 0.3921976592941051 , 0.0000325275043721 );
            ir->AddPentPoint( 438 , 0.0403754775891479 , 0.1348537262334940 , 0.0403754775891479 , 0.3921976592941051 , 0.0000325275043721 );
            ir->AddPentPoint( 439 , 0.1348537262334940 , 0.3921976592941051 , 0.0403754775891479 , 0.3921976592941051 , 0.0000325275043721 );
            ir->AddPentPoint( 440 , 0.0403754775891479 , 0.3921976592941051 , 0.0403754775891479 , 0.3921976592941051 , 0.0000325275043721 );
            ir->AddPentPoint( 441 , 0.1348537262334940 , 0.0403754775891479 , 0.0403754775891479 , 0.3921976592941051 , 0.0000325275043721 );
            ir->AddPentPoint( 442 , 0.3921976592941051 , 0.0403754775891479 , 0.0403754775891479 , 0.3921976592941051 , 0.0000325275043721 );
            ir->AddPentPoint( 443 , 0.3921976592941051 , 0.3921976592941051 , 0.1348537262334940 , 0.0403754775891479 , 0.0000325275043721 );
            ir->AddPentPoint( 444 , 0.0403754775891479 , 0.3921976592941051 , 0.1348537262334940 , 0.0403754775891479 , 0.0000325275043721 );
            ir->AddPentPoint( 445 , 0.3921976592941051 , 0.0403754775891479 , 0.1348537262334940 , 0.0403754775891479 , 0.0000325275043721 );
            ir->AddPentPoint( 446 , 0.3921976592941051 , 0.1348537262334940 , 0.3921976592941051 , 0.0403754775891479 , 0.0000325275043721 );
            ir->AddPentPoint( 447 , 0.0403754775891479 , 0.1348537262334940 , 0.3921976592941051 , 0.0403754775891479 , 0.0000325275043721 );
            ir->AddPentPoint( 448 , 0.1348537262334940 , 0.3921976592941051 , 0.3921976592941051 , 0.0403754775891479 , 0.0000325275043721 );
            ir->AddPentPoint( 449 , 0.0403754775891479 , 0.3921976592941051 , 0.3921976592941051 , 0.0403754775891479 , 0.0000325275043721 );
            ir->AddPentPoint( 450 , 0.1348537262334940 , 0.0403754775891479 , 0.3921976592941051 , 0.0403754775891479 , 0.0000325275043721 );
            ir->AddPentPoint( 451 , 0.3921976592941051 , 0.0403754775891479 , 0.3921976592941051 , 0.0403754775891479 , 0.0000325275043721 );
            ir->AddPentPoint( 452 , 0.3921976592941051 , 0.1348537262334940 , 0.0403754775891479 , 0.0403754775891479 , 0.0000325275043721 );
            ir->AddPentPoint( 453 , 0.1348537262334940 , 0.3921976592941051 , 0.0403754775891479 , 0.0403754775891479 , 0.0000325275043721 );
            ir->AddPentPoint( 454 , 0.3921976592941051 , 0.3921976592941051 , 0.0403754775891479 , 0.0403754775891479 , 0.0000325275043721 );
            ir->AddPentPoint( 455 , 0.0449581283040572 , 0.1543547700531318 , 0.1543547700531318 , 0.6013742032856220 , 0.0000491662886531 );
            ir->AddPentPoint( 456 , 0.1543547700531318 , 0.0449581283040572 , 0.1543547700531318 , 0.6013742032856220 , 0.0000491662886531 );
            ir->AddPentPoint( 457 , 0.0449581283040572 , 0.0449581283040572 , 0.1543547700531318 , 0.6013742032856220 , 0.0000491662886531 );
            ir->AddPentPoint( 458 , 0.1543547700531318 , 0.1543547700531318 , 0.0449581283040572 , 0.6013742032856220 , 0.0000491662886531 );
            ir->AddPentPoint( 459 , 0.0449581283040572 , 0.1543547700531318 , 0.0449581283040572 , 0.6013742032856220 , 0.0000491662886531 );
            ir->AddPentPoint( 460 , 0.1543547700531318 , 0.0449581283040572 , 0.0449581283040572 , 0.6013742032856220 , 0.0000491662886531 );
            ir->AddPentPoint( 461 , 0.0449581283040572 , 0.1543547700531318 , 0.6013742032856220 , 0.1543547700531318 , 0.0000491662886531 );
            ir->AddPentPoint( 462 , 0.1543547700531318 , 0.0449581283040572 , 0.6013742032856220 , 0.1543547700531318 , 0.0000491662886531 );
            ir->AddPentPoint( 463 , 0.0449581283040572 , 0.0449581283040572 , 0.6013742032856220 , 0.1543547700531318 , 0.0000491662886531 );
            ir->AddPentPoint( 464 , 0.0449581283040572 , 0.6013742032856220 , 0.1543547700531318 , 0.1543547700531318 , 0.0000491662886531 );
            ir->AddPentPoint( 465 , 0.6013742032856220 , 0.0449581283040572 , 0.1543547700531318 , 0.1543547700531318 , 0.0000491662886531 );
            ir->AddPentPoint( 466 , 0.0449581283040572 , 0.0449581283040572 , 0.1543547700531318 , 0.1543547700531318 , 0.0000491662886531 );
            ir->AddPentPoint( 467 , 0.1543547700531318 , 0.6013742032856220 , 0.0449581283040572 , 0.1543547700531318 , 0.0000491662886531 );
            ir->AddPentPoint( 468 , 0.0449581283040572 , 0.6013742032856220 , 0.0449581283040572 , 0.1543547700531318 , 0.0000491662886531 );
            ir->AddPentPoint( 469 , 0.6013742032856220 , 0.1543547700531318 , 0.0449581283040572 , 0.1543547700531318 , 0.0000491662886531 );
            ir->AddPentPoint( 470 , 0.0449581283040572 , 0.1543547700531318 , 0.0449581283040572 , 0.1543547700531318 , 0.0000491662886531 );
            ir->AddPentPoint( 471 , 0.6013742032856220 , 0.0449581283040572 , 0.0449581283040572 , 0.1543547700531318 , 0.0000491662886531 );
            ir->AddPentPoint( 472 , 0.1543547700531318 , 0.0449581283040572 , 0.0449581283040572 , 0.1543547700531318 , 0.0000491662886531 );
            ir->AddPentPoint( 473 , 0.1543547700531318 , 0.1543547700531318 , 0.6013742032856220 , 0.0449581283040572 , 0.0000491662886531 );
            ir->AddPentPoint( 474 , 0.0449581283040572 , 0.1543547700531318 , 0.6013742032856220 , 0.0449581283040572 , 0.0000491662886531 );
            ir->AddPentPoint( 475 , 0.1543547700531318 , 0.0449581283040572 , 0.6013742032856220 , 0.0449581283040572 , 0.0000491662886531 );
            ir->AddPentPoint( 476 , 0.1543547700531318 , 0.6013742032856220 , 0.1543547700531318 , 0.0449581283040572 , 0.0000491662886531 );
            ir->AddPentPoint( 477 , 0.0449581283040572 , 0.6013742032856220 , 0.1543547700531318 , 0.0449581283040572 , 0.0000491662886531 );
            ir->AddPentPoint( 478 , 0.6013742032856220 , 0.1543547700531318 , 0.1543547700531318 , 0.0449581283040572 , 0.0000491662886531 );
            ir->AddPentPoint( 479 , 0.0449581283040572 , 0.1543547700531318 , 0.1543547700531318 , 0.0449581283040572 , 0.0000491662886531 );
            ir->AddPentPoint( 480 , 0.6013742032856220 , 0.0449581283040572 , 0.1543547700531318 , 0.0449581283040572 , 0.0000491662886531 );
            ir->AddPentPoint( 481 , 0.1543547700531318 , 0.0449581283040572 , 0.1543547700531318 , 0.0449581283040572 , 0.0000491662886531 );
            ir->AddPentPoint( 482 , 0.1543547700531318 , 0.6013742032856220 , 0.0449581283040572 , 0.0449581283040572 , 0.0000491662886531 );
            ir->AddPentPoint( 483 , 0.6013742032856220 , 0.1543547700531318 , 0.0449581283040572 , 0.0449581283040572 , 0.0000491662886531 );
            ir->AddPentPoint( 484 , 0.1543547700531318 , 0.1543547700531318 , 0.0449581283040572 , 0.0449581283040572 , 0.0000491662886531 );
            ir->AddPentPoint( 485 , 0.0576496002944590 , 0.2361915441598070 , 0.2361915441598070 , 0.4123177110914679 , 0.0001221121983753 );
            ir->AddPentPoint( 486 , 0.2361915441598070 , 0.0576496002944590 , 0.2361915441598070 , 0.4123177110914679 , 0.0001221121983753 );
            ir->AddPentPoint( 487 , 0.0576496002944590 , 0.0576496002944590 , 0.2361915441598070 , 0.4123177110914679 , 0.0001221121983753 );
            ir->AddPentPoint( 488 , 0.2361915441598070 , 0.2361915441598070 , 0.0576496002944590 , 0.4123177110914679 , 0.0001221121983753 );
            ir->AddPentPoint( 489 , 0.0576496002944590 , 0.2361915441598070 , 0.0576496002944590 , 0.4123177110914679 , 0.0001221121983753 );
            ir->AddPentPoint( 490 , 0.2361915441598070 , 0.0576496002944590 , 0.0576496002944590 , 0.4123177110914679 , 0.0001221121983753 );
            ir->AddPentPoint( 491 , 0.0576496002944590 , 0.2361915441598070 , 0.4123177110914679 , 0.2361915441598070 , 0.0001221121983753 );
            ir->AddPentPoint( 492 , 0.2361915441598070 , 0.0576496002944590 , 0.4123177110914679 , 0.2361915441598070 , 0.0001221121983753 );
            ir->AddPentPoint( 493 , 0.0576496002944590 , 0.0576496002944590 , 0.4123177110914679 , 0.2361915441598070 , 0.0001221121983753 );
            ir->AddPentPoint( 494 , 0.0576496002944590 , 0.4123177110914679 , 0.2361915441598070 , 0.2361915441598070 , 0.0001221121983753 );
            ir->AddPentPoint( 495 , 0.4123177110914679 , 0.0576496002944590 , 0.2361915441598070 , 0.2361915441598070 , 0.0001221121983753 );
            ir->AddPentPoint( 496 , 0.0576496002944590 , 0.0576496002944590 , 0.2361915441598070 , 0.2361915441598070 , 0.0001221121983753 );
            ir->AddPentPoint( 497 , 0.2361915441598070 , 0.4123177110914679 , 0.0576496002944590 , 0.2361915441598070 , 0.0001221121983753 );
            ir->AddPentPoint( 498 , 0.0576496002944590 , 0.4123177110914679 , 0.0576496002944590 , 0.2361915441598070 , 0.0001221121983753 );
            ir->AddPentPoint( 499 , 0.4123177110914679 , 0.2361915441598070 , 0.0576496002944590 , 0.2361915441598070 , 0.0001221121983753 );
            ir->AddPentPoint( 500 , 0.0576496002944590 , 0.2361915441598070 , 0.0576496002944590 , 0.2361915441598070 , 0.0001221121983753 );
            ir->AddPentPoint( 501 , 0.4123177110914679 , 0.0576496002944590 , 0.0576496002944590 , 0.2361915441598070 , 0.0001221121983753 );
            ir->AddPentPoint( 502 , 0.2361915441598070 , 0.0576496002944590 , 0.0576496002944590 , 0.2361915441598070 , 0.0001221121983753 );
            ir->AddPentPoint( 503 , 0.2361915441598070 , 0.2361915441598070 , 0.4123177110914679 , 0.0576496002944590 , 0.0001221121983753 );
            ir->AddPentPoint( 504 , 0.0576496002944590 , 0.2361915441598070 , 0.4123177110914679 , 0.0576496002944590 , 0.0001221121983753 );
            ir->AddPentPoint( 505 , 0.2361915441598070 , 0.0576496002944590 , 0.4123177110914679 , 0.0576496002944590 , 0.0001221121983753 );
            ir->AddPentPoint( 506 , 0.2361915441598070 , 0.4123177110914679 , 0.2361915441598070 , 0.0576496002944590 , 0.0001221121983753 );
            ir->AddPentPoint( 507 , 0.0576496002944590 , 0.4123177110914679 , 0.2361915441598070 , 0.0576496002944590 , 0.0001221121983753 );
            ir->AddPentPoint( 508 , 0.4123177110914679 , 0.2361915441598070 , 0.2361915441598070 , 0.0576496002944590 , 0.0001221121983753 );
            ir->AddPentPoint( 509 , 0.0576496002944590 , 0.2361915441598070 , 0.2361915441598070 , 0.0576496002944590 , 0.0001221121983753 );
            ir->AddPentPoint( 510 , 0.4123177110914679 , 0.0576496002944590 , 0.2361915441598070 , 0.0576496002944590 , 0.0001221121983753 );
            ir->AddPentPoint( 511 , 0.2361915441598070 , 0.0576496002944590 , 0.2361915441598070 , 0.0576496002944590 , 0.0001221121983753 );
            ir->AddPentPoint( 512 , 0.2361915441598070 , 0.4123177110914679 , 0.0576496002944590 , 0.0576496002944590 , 0.0001221121983753 );
            ir->AddPentPoint( 513 , 0.4123177110914679 , 0.2361915441598070 , 0.0576496002944590 , 0.0576496002944590 , 0.0001221121983753 );
            ir->AddPentPoint( 514 , 0.2361915441598070 , 0.2361915441598070 , 0.0576496002944590 , 0.0576496002944590 , 0.0001221121983753 );
            ir->AddPentPoint( 515 , 0.1334996251965174 , 0.3284191033922808 , 0.3284191033922808 , 0.0761625428224036 , 0.0000995317498740 );
            ir->AddPentPoint( 516 , 0.3284191033922808 , 0.1334996251965174 , 0.3284191033922808 , 0.0761625428224036 , 0.0000995317498740 );
            ir->AddPentPoint( 517 , 0.1334996251965174 , 0.1334996251965174 , 0.3284191033922808 , 0.0761625428224036 , 0.0000995317498740 );
            ir->AddPentPoint( 518 , 0.3284191033922808 , 0.3284191033922808 , 0.1334996251965174 , 0.0761625428224036 , 0.0000995317498740 );
            ir->AddPentPoint( 519 , 0.1334996251965174 , 0.3284191033922808 , 0.1334996251965174 , 0.0761625428224036 , 0.0000995317498740 );
            ir->AddPentPoint( 520 , 0.3284191033922808 , 0.1334996251965174 , 0.1334996251965174 , 0.0761625428224036 , 0.0000995317498740 );
            ir->AddPentPoint( 521 , 0.1334996251965174 , 0.3284191033922808 , 0.0761625428224036 , 0.3284191033922808 , 0.0000995317498740 );
            ir->AddPentPoint( 522 , 0.3284191033922808 , 0.1334996251965174 , 0.0761625428224036 , 0.3284191033922808 , 0.0000995317498740 );
            ir->AddPentPoint( 523 , 0.1334996251965174 , 0.1334996251965174 , 0.0761625428224036 , 0.3284191033922808 , 0.0000995317498740 );
            ir->AddPentPoint( 524 , 0.1334996251965174 , 0.0761625428224036 , 0.3284191033922808 , 0.3284191033922808 , 0.0000995317498740 );
            ir->AddPentPoint( 525 , 0.0761625428224036 , 0.1334996251965174 , 0.3284191033922808 , 0.3284191033922808 , 0.0000995317498740 );
            ir->AddPentPoint( 526 , 0.1334996251965174 , 0.1334996251965174 , 0.3284191033922808 , 0.3284191033922808 , 0.0000995317498740 );
            ir->AddPentPoint( 527 , 0.3284191033922808 , 0.0761625428224036 , 0.1334996251965174 , 0.3284191033922808 , 0.0000995317498740 );
            ir->AddPentPoint( 528 , 0.1334996251965174 , 0.0761625428224036 , 0.1334996251965174 , 0.3284191033922808 , 0.0000995317498740 );
            ir->AddPentPoint( 529 , 0.0761625428224036 , 0.3284191033922808 , 0.1334996251965174 , 0.3284191033922808 , 0.0000995317498740 );
            ir->AddPentPoint( 530 , 0.1334996251965174 , 0.3284191033922808 , 0.1334996251965174 , 0.3284191033922808 , 0.0000995317498740 );
            ir->AddPentPoint( 531 , 0.0761625428224036 , 0.1334996251965174 , 0.1334996251965174 , 0.3284191033922808 , 0.0000995317498740 );
            ir->AddPentPoint( 532 , 0.3284191033922808 , 0.1334996251965174 , 0.1334996251965174 , 0.3284191033922808 , 0.0000995317498740 );
            ir->AddPentPoint( 533 , 0.3284191033922808 , 0.3284191033922808 , 0.0761625428224036 , 0.1334996251965174 , 0.0000995317498740 );
            ir->AddPentPoint( 534 , 0.1334996251965174 , 0.3284191033922808 , 0.0761625428224036 , 0.1334996251965174 , 0.0000995317498740 );
            ir->AddPentPoint( 535 , 0.3284191033922808 , 0.1334996251965174 , 0.0761625428224036 , 0.1334996251965174 , 0.0000995317498740 );
            ir->AddPentPoint( 536 , 0.3284191033922808 , 0.0761625428224036 , 0.3284191033922808 , 0.1334996251965174 , 0.0000995317498740 );
            ir->AddPentPoint( 537 , 0.1334996251965174 , 0.0761625428224036 , 0.3284191033922808 , 0.1334996251965174 , 0.0000995317498740 );
            ir->AddPentPoint( 538 , 0.0761625428224036 , 0.3284191033922808 , 0.3284191033922808 , 0.1334996251965174 , 0.0000995317498740 );
            ir->AddPentPoint( 539 , 0.1334996251965174 , 0.3284191033922808 , 0.3284191033922808 , 0.1334996251965174 , 0.0000995317498740 );
            ir->AddPentPoint( 540 , 0.0761625428224036 , 0.1334996251965174 , 0.3284191033922808 , 0.1334996251965174 , 0.0000995317498740 );
            ir->AddPentPoint( 541 , 0.3284191033922808 , 0.1334996251965174 , 0.3284191033922808 , 0.1334996251965174 , 0.0000995317498740 );
            ir->AddPentPoint( 542 , 0.3284191033922808 , 0.0761625428224036 , 0.1334996251965174 , 0.1334996251965174 , 0.0000995317498740 );
            ir->AddPentPoint( 543 , 0.0761625428224036 , 0.3284191033922808 , 0.1334996251965174 , 0.1334996251965174 , 0.0000995317498740 );
            ir->AddPentPoint( 544 , 0.3284191033922808 , 0.3284191033922808 , 0.1334996251965174 , 0.1334996251965174 , 0.0000995317498740 );
            ir->AddPentPoint( 545 , 0.1964198071512829 , 0.3008745000470239 , 0.3008745000470239 , 0.0054113856033864 , 0.0000078137298661 );
            ir->AddPentPoint( 546 , 0.3008745000470239 , 0.1964198071512829 , 0.3008745000470239 , 0.0054113856033864 , 0.0000078137298661 );
            ir->AddPentPoint( 547 , 0.1964198071512829 , 0.1964198071512829 , 0.3008745000470239 , 0.0054113856033864 , 0.0000078137298661 );
            ir->AddPentPoint( 548 , 0.3008745000470239 , 0.3008745000470239 , 0.1964198071512829 , 0.0054113856033864 , 0.0000078137298661 );
            ir->AddPentPoint( 549 , 0.1964198071512829 , 0.3008745000470239 , 0.1964198071512829 , 0.0054113856033864 , 0.0000078137298661 );
            ir->AddPentPoint( 550 , 0.3008745000470239 , 0.1964198071512829 , 0.1964198071512829 , 0.0054113856033864 , 0.0000078137298661 );
            ir->AddPentPoint( 551 , 0.1964198071512829 , 0.3008745000470239 , 0.0054113856033864 , 0.3008745000470239 , 0.0000078137298661 );
            ir->AddPentPoint( 552 , 0.3008745000470239 , 0.1964198071512829 , 0.0054113856033864 , 0.3008745000470239 , 0.0000078137298661 );
            ir->AddPentPoint( 553 , 0.1964198071512829 , 0.1964198071512829 , 0.0054113856033864 , 0.3008745000470239 , 0.0000078137298661 );
            ir->AddPentPoint( 554 , 0.1964198071512829 , 0.0054113856033864 , 0.3008745000470239 , 0.3008745000470239 , 0.0000078137298661 );
            ir->AddPentPoint( 555 , 0.0054113856033864 , 0.1964198071512829 , 0.3008745000470239 , 0.3008745000470239 , 0.0000078137298661 );
            ir->AddPentPoint( 556 , 0.1964198071512829 , 0.1964198071512829 , 0.3008745000470239 , 0.3008745000470239 , 0.0000078137298661 );
            ir->AddPentPoint( 557 , 0.3008745000470239 , 0.0054113856033864 , 0.1964198071512829 , 0.3008745000470239 , 0.0000078137298661 );
            ir->AddPentPoint( 558 , 0.1964198071512829 , 0.0054113856033864 , 0.1964198071512829 , 0.3008745000470239 , 0.0000078137298661 );
            ir->AddPentPoint( 559 , 0.0054113856033864 , 0.3008745000470239 , 0.1964198071512829 , 0.3008745000470239 , 0.0000078137298661 );
            ir->AddPentPoint( 560 , 0.1964198071512829 , 0.3008745000470239 , 0.1964198071512829 , 0.3008745000470239 , 0.0000078137298661 );
            ir->AddPentPoint( 561 , 0.0054113856033864 , 0.1964198071512829 , 0.1964198071512829 , 0.3008745000470239 , 0.0000078137298661 );
            ir->AddPentPoint( 562 , 0.3008745000470239 , 0.1964198071512829 , 0.1964198071512829 , 0.3008745000470239 , 0.0000078137298661 );
            ir->AddPentPoint( 563 , 0.3008745000470239 , 0.3008745000470239 , 0.0054113856033864 , 0.1964198071512829 , 0.0000078137298661 );
            ir->AddPentPoint( 564 , 0.1964198071512829 , 0.3008745000470239 , 0.0054113856033864 , 0.1964198071512829 , 0.0000078137298661 );
            ir->AddPentPoint( 565 , 0.3008745000470239 , 0.1964198071512829 , 0.0054113856033864 , 0.1964198071512829 , 0.0000078137298661 );
            ir->AddPentPoint( 566 , 0.3008745000470239 , 0.0054113856033864 , 0.3008745000470239 , 0.1964198071512829 , 0.0000078137298661 );
            ir->AddPentPoint( 567 , 0.1964198071512829 , 0.0054113856033864 , 0.3008745000470239 , 0.1964198071512829 , 0.0000078137298661 );
            ir->AddPentPoint( 568 , 0.0054113856033864 , 0.3008745000470239 , 0.3008745000470239 , 0.1964198071512829 , 0.0000078137298661 );
            ir->AddPentPoint( 569 , 0.1964198071512829 , 0.3008745000470239 , 0.3008745000470239 , 0.1964198071512829 , 0.0000078137298661 );
            ir->AddPentPoint( 570 , 0.0054113856033864 , 0.1964198071512829 , 0.3008745000470239 , 0.1964198071512829 , 0.0000078137298661 );
            ir->AddPentPoint( 571 , 0.3008745000470239 , 0.1964198071512829 , 0.3008745000470239 , 0.1964198071512829 , 0.0000078137298661 );
            ir->AddPentPoint( 572 , 0.3008745000470239 , 0.0054113856033864 , 0.1964198071512829 , 0.1964198071512829 , 0.0000078137298661 );
            ir->AddPentPoint( 573 , 0.0054113856033864 , 0.3008745000470239 , 0.1964198071512829 , 0.1964198071512829 , 0.0000078137298661 );
            ir->AddPentPoint( 574 , 0.3008745000470239 , 0.3008745000470239 , 0.1964198071512829 , 0.1964198071512829 , 0.0000078137298661 );
            ir->AddPentPoint( 575 , 0.0000000001507476 , 0.0753397736760257 , 0.2432370253979445 , 0.6814232006245345 , 0.0000032345030169 );
            ir->AddPentPoint( 576 , 0.0753397736760257 , 0.0000000001507476 , 0.2432370253979445 , 0.6814232006245345 , 0.0000032345030169 );
            ir->AddPentPoint( 577 , 0.0000000001507476 , 0.0000000001507476 , 0.2432370253979445 , 0.6814232006245345 , 0.0000032345030169 );
            ir->AddPentPoint( 578 , 0.0000000001507476 , 0.2432370253979445 , 0.0753397736760257 , 0.6814232006245345 , 0.0000032345030169 );
            ir->AddPentPoint( 579 , 0.2432370253979445 , 0.0000000001507476 , 0.0753397736760257 , 0.6814232006245345 , 0.0000032345030169 );
            ir->AddPentPoint( 580 , 0.0000000001507476 , 0.0000000001507476 , 0.0753397736760257 , 0.6814232006245345 , 0.0000032345030169 );
            ir->AddPentPoint( 581 , 0.0753397736760257 , 0.2432370253979445 , 0.0000000001507476 , 0.6814232006245345 , 0.0000032345030169 );
            ir->AddPentPoint( 582 , 0.0000000001507476 , 0.2432370253979445 , 0.0000000001507476 , 0.6814232006245345 , 0.0000032345030169 );
            ir->AddPentPoint( 583 , 0.2432370253979445 , 0.0753397736760257 , 0.0000000001507476 , 0.6814232006245345 , 0.0000032345030169 );
            ir->AddPentPoint( 584 , 0.0000000001507476 , 0.0753397736760257 , 0.0000000001507476 , 0.6814232006245345 , 0.0000032345030169 );
            ir->AddPentPoint( 585 , 0.2432370253979445 , 0.0000000001507476 , 0.0000000001507476 , 0.6814232006245345 , 0.0000032345030169 );
            ir->AddPentPoint( 586 , 0.0753397736760257 , 0.0000000001507476 , 0.0000000001507476 , 0.6814232006245345 , 0.0000032345030169 );
            ir->AddPentPoint( 587 , 0.0000000001507476 , 0.0753397736760257 , 0.6814232006245345 , 0.2432370253979445 , 0.0000032345030169 );
            ir->AddPentPoint( 588 , 0.0753397736760257 , 0.0000000001507476 , 0.6814232006245345 , 0.2432370253979445 , 0.0000032345030169 );
            ir->AddPentPoint( 589 , 0.0000000001507476 , 0.0000000001507476 , 0.6814232006245345 , 0.2432370253979445 , 0.0000032345030169 );
            ir->AddPentPoint( 590 , 0.0000000001507476 , 0.6814232006245345 , 0.0753397736760257 , 0.2432370253979445 , 0.0000032345030169 );
            ir->AddPentPoint( 591 , 0.6814232006245345 , 0.0000000001507476 , 0.0753397736760257 , 0.2432370253979445 , 0.0000032345030169 );
            ir->AddPentPoint( 592 , 0.0000000001507476 , 0.0000000001507476 , 0.0753397736760257 , 0.2432370253979445 , 0.0000032345030169 );
            ir->AddPentPoint( 593 , 0.0753397736760257 , 0.6814232006245345 , 0.0000000001507476 , 0.2432370253979445 , 0.0000032345030169 );
            ir->AddPentPoint( 594 , 0.0000000001507476 , 0.6814232006245345 , 0.0000000001507476 , 0.2432370253979445 , 0.0000032345030169 );
            ir->AddPentPoint( 595 , 0.6814232006245345 , 0.0753397736760257 , 0.0000000001507476 , 0.2432370253979445 , 0.0000032345030169 );
            ir->AddPentPoint( 596 , 0.0000000001507476 , 0.0753397736760257 , 0.0000000001507476 , 0.2432370253979445 , 0.0000032345030169 );
            ir->AddPentPoint( 597 , 0.6814232006245345 , 0.0000000001507476 , 0.0000000001507476 , 0.2432370253979445 , 0.0000032345030169 );
            ir->AddPentPoint( 598 , 0.0753397736760257 , 0.0000000001507476 , 0.0000000001507476 , 0.2432370253979445 , 0.0000032345030169 );
            ir->AddPentPoint( 599 , 0.0000000001507476 , 0.2432370253979445 , 0.6814232006245345 , 0.0753397736760257 , 0.0000032345030169 );
            ir->AddPentPoint( 600 , 0.2432370253979445 , 0.0000000001507476 , 0.6814232006245345 , 0.0753397736760257 , 0.0000032345030169 );
            ir->AddPentPoint( 601 , 0.0000000001507476 , 0.0000000001507476 , 0.6814232006245345 , 0.0753397736760257 , 0.0000032345030169 );
            ir->AddPentPoint( 602 , 0.0000000001507476 , 0.6814232006245345 , 0.2432370253979445 , 0.0753397736760257 , 0.0000032345030169 );
            ir->AddPentPoint( 603 , 0.6814232006245345 , 0.0000000001507476 , 0.2432370253979445 , 0.0753397736760257 , 0.0000032345030169 );
            ir->AddPentPoint( 604 , 0.0000000001507476 , 0.0000000001507476 , 0.2432370253979445 , 0.0753397736760257 , 0.0000032345030169 );
            ir->AddPentPoint( 605 , 0.2432370253979445 , 0.6814232006245345 , 0.0000000001507476 , 0.0753397736760257 , 0.0000032345030169 );
            ir->AddPentPoint( 606 , 0.0000000001507476 , 0.6814232006245345 , 0.0000000001507476 , 0.0753397736760257 , 0.0000032345030169 );
            ir->AddPentPoint( 607 , 0.6814232006245345 , 0.2432370253979445 , 0.0000000001507476 , 0.0753397736760257 , 0.0000032345030169 );
            ir->AddPentPoint( 608 , 0.0000000001507476 , 0.2432370253979445 , 0.0000000001507476 , 0.0753397736760257 , 0.0000032345030169 );
            ir->AddPentPoint( 609 , 0.6814232006245345 , 0.0000000001507476 , 0.0000000001507476 , 0.0753397736760257 , 0.0000032345030169 );
            ir->AddPentPoint( 610 , 0.2432370253979445 , 0.0000000001507476 , 0.0000000001507476 , 0.0753397736760257 , 0.0000032345030169 );
            ir->AddPentPoint( 611 , 0.0753397736760257 , 0.2432370253979445 , 0.6814232006245345 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 612 , 0.0000000001507476 , 0.2432370253979445 , 0.6814232006245345 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 613 , 0.2432370253979445 , 0.0753397736760257 , 0.6814232006245345 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 614 , 0.0000000001507476 , 0.0753397736760257 , 0.6814232006245345 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 615 , 0.2432370253979445 , 0.0000000001507476 , 0.6814232006245345 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 616 , 0.0753397736760257 , 0.0000000001507476 , 0.6814232006245345 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 617 , 0.0753397736760257 , 0.6814232006245345 , 0.2432370253979445 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 618 , 0.0000000001507476 , 0.6814232006245345 , 0.2432370253979445 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 619 , 0.6814232006245345 , 0.0753397736760257 , 0.2432370253979445 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 620 , 0.0000000001507476 , 0.0753397736760257 , 0.2432370253979445 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 621 , 0.6814232006245345 , 0.0000000001507476 , 0.2432370253979445 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 622 , 0.0753397736760257 , 0.0000000001507476 , 0.2432370253979445 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 623 , 0.2432370253979445 , 0.6814232006245345 , 0.0753397736760257 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 624 , 0.0000000001507476 , 0.6814232006245345 , 0.0753397736760257 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 625 , 0.6814232006245345 , 0.2432370253979445 , 0.0753397736760257 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 626 , 0.0000000001507476 , 0.2432370253979445 , 0.0753397736760257 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 627 , 0.6814232006245345 , 0.0000000001507476 , 0.0753397736760257 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 628 , 0.2432370253979445 , 0.0000000001507476 , 0.0753397736760257 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 629 , 0.2432370253979445 , 0.6814232006245345 , 0.0000000001507476 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 630 , 0.0753397736760257 , 0.6814232006245345 , 0.0000000001507476 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 631 , 0.6814232006245345 , 0.2432370253979445 , 0.0000000001507476 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 632 , 0.0753397736760257 , 0.2432370253979445 , 0.0000000001507476 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 633 , 0.6814232006245345 , 0.0753397736760257 , 0.0000000001507476 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 634 , 0.2432370253979445 , 0.0753397736760257 , 0.0000000001507476 , 0.0000000001507476 , 0.0000032345030169 );
            ir->AddPentPoint( 635 , 0.0229518395598165 , 0.1240602798626259 , 0.2923308644281944 , 0.5377051765895465 , 0.0000375587155977 );
            ir->AddPentPoint( 636 , 0.1240602798626259 , 0.0229518395598165 , 0.2923308644281944 , 0.5377051765895465 , 0.0000375587155977 );
            ir->AddPentPoint( 637 , 0.0229518395598165 , 0.0229518395598165 , 0.2923308644281944 , 0.5377051765895465 , 0.0000375587155977 );
            ir->AddPentPoint( 638 , 0.0229518395598165 , 0.2923308644281944 , 0.1240602798626259 , 0.5377051765895465 , 0.0000375587155977 );
            ir->AddPentPoint( 639 , 0.2923308644281944 , 0.0229518395598165 , 0.1240602798626259 , 0.5377051765895465 , 0.0000375587155977 );
            ir->AddPentPoint( 640 , 0.0229518395598165 , 0.0229518395598165 , 0.1240602798626259 , 0.5377051765895465 , 0.0000375587155977 );
            ir->AddPentPoint( 641 , 0.1240602798626259 , 0.2923308644281944 , 0.0229518395598165 , 0.5377051765895465 , 0.0000375587155977 );
            ir->AddPentPoint( 642 , 0.0229518395598165 , 0.2923308644281944 , 0.0229518395598165 , 0.5377051765895465 , 0.0000375587155977 );
            ir->AddPentPoint( 643 , 0.2923308644281944 , 0.1240602798626259 , 0.0229518395598165 , 0.5377051765895465 , 0.0000375587155977 );
            ir->AddPentPoint( 644 , 0.0229518395598165 , 0.1240602798626259 , 0.0229518395598165 , 0.5377051765895465 , 0.0000375587155977 );
            ir->AddPentPoint( 645 , 0.2923308644281944 , 0.0229518395598165 , 0.0229518395598165 , 0.5377051765895465 , 0.0000375587155977 );
            ir->AddPentPoint( 646 , 0.1240602798626259 , 0.0229518395598165 , 0.0229518395598165 , 0.5377051765895465 , 0.0000375587155977 );
            ir->AddPentPoint( 647 , 0.0229518395598165 , 0.1240602798626259 , 0.5377051765895465 , 0.2923308644281944 , 0.0000375587155977 );
            ir->AddPentPoint( 648 , 0.1240602798626259 , 0.0229518395598165 , 0.5377051765895465 , 0.2923308644281944 , 0.0000375587155977 );
            ir->AddPentPoint( 649 , 0.0229518395598165 , 0.0229518395598165 , 0.5377051765895465 , 0.2923308644281944 , 0.0000375587155977 );
            ir->AddPentPoint( 650 , 0.0229518395598165 , 0.5377051765895465 , 0.1240602798626259 , 0.2923308644281944 , 0.0000375587155977 );
            ir->AddPentPoint( 651 , 0.5377051765895465 , 0.0229518395598165 , 0.1240602798626259 , 0.2923308644281944 , 0.0000375587155977 );
            ir->AddPentPoint( 652 , 0.0229518395598165 , 0.0229518395598165 , 0.1240602798626259 , 0.2923308644281944 , 0.0000375587155977 );
            ir->AddPentPoint( 653 , 0.1240602798626259 , 0.5377051765895465 , 0.0229518395598165 , 0.2923308644281944 , 0.0000375587155977 );
            ir->AddPentPoint( 654 , 0.0229518395598165 , 0.5377051765895465 , 0.0229518395598165 , 0.2923308644281944 , 0.0000375587155977 );
            ir->AddPentPoint( 655 , 0.5377051765895465 , 0.1240602798626259 , 0.0229518395598165 , 0.2923308644281944 , 0.0000375587155977 );
            ir->AddPentPoint( 656 , 0.0229518395598165 , 0.1240602798626259 , 0.0229518395598165 , 0.2923308644281944 , 0.0000375587155977 );
            ir->AddPentPoint( 657 , 0.5377051765895465 , 0.0229518395598165 , 0.0229518395598165 , 0.2923308644281944 , 0.0000375587155977 );
            ir->AddPentPoint( 658 , 0.1240602798626259 , 0.0229518395598165 , 0.0229518395598165 , 0.2923308644281944 , 0.0000375587155977 );
            ir->AddPentPoint( 659 , 0.0229518395598165 , 0.2923308644281944 , 0.5377051765895465 , 0.1240602798626259 , 0.0000375587155977 );
            ir->AddPentPoint( 660 , 0.2923308644281944 , 0.0229518395598165 , 0.5377051765895465 , 0.1240602798626259 , 0.0000375587155977 );
            ir->AddPentPoint( 661 , 0.0229518395598165 , 0.0229518395598165 , 0.5377051765895465 , 0.1240602798626259 , 0.0000375587155977 );
            ir->AddPentPoint( 662 , 0.0229518395598165 , 0.5377051765895465 , 0.2923308644281944 , 0.1240602798626259 , 0.0000375587155977 );
            ir->AddPentPoint( 663 , 0.5377051765895465 , 0.0229518395598165 , 0.2923308644281944 , 0.1240602798626259 , 0.0000375587155977 );
            ir->AddPentPoint( 664 , 0.0229518395598165 , 0.0229518395598165 , 0.2923308644281944 , 0.1240602798626259 , 0.0000375587155977 );
            ir->AddPentPoint( 665 , 0.2923308644281944 , 0.5377051765895465 , 0.0229518395598165 , 0.1240602798626259 , 0.0000375587155977 );
            ir->AddPentPoint( 666 , 0.0229518395598165 , 0.5377051765895465 , 0.0229518395598165 , 0.1240602798626259 , 0.0000375587155977 );
            ir->AddPentPoint( 667 , 0.5377051765895465 , 0.2923308644281944 , 0.0229518395598165 , 0.1240602798626259 , 0.0000375587155977 );
            ir->AddPentPoint( 668 , 0.0229518395598165 , 0.2923308644281944 , 0.0229518395598165 , 0.1240602798626259 , 0.0000375587155977 );
            ir->AddPentPoint( 669 , 0.5377051765895465 , 0.0229518395598165 , 0.0229518395598165 , 0.1240602798626259 , 0.0000375587155977 );
            ir->AddPentPoint( 670 , 0.2923308644281944 , 0.0229518395598165 , 0.0229518395598165 , 0.1240602798626259 , 0.0000375587155977 );
            ir->AddPentPoint( 671 , 0.1240602798626259 , 0.2923308644281944 , 0.5377051765895465 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 672 , 0.0229518395598165 , 0.2923308644281944 , 0.5377051765895465 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 673 , 0.2923308644281944 , 0.1240602798626259 , 0.5377051765895465 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 674 , 0.0229518395598165 , 0.1240602798626259 , 0.5377051765895465 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 675 , 0.2923308644281944 , 0.0229518395598165 , 0.5377051765895465 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 676 , 0.1240602798626259 , 0.0229518395598165 , 0.5377051765895465 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 677 , 0.1240602798626259 , 0.5377051765895465 , 0.2923308644281944 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 678 , 0.0229518395598165 , 0.5377051765895465 , 0.2923308644281944 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 679 , 0.5377051765895465 , 0.1240602798626259 , 0.2923308644281944 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 680 , 0.0229518395598165 , 0.1240602798626259 , 0.2923308644281944 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 681 , 0.5377051765895465 , 0.0229518395598165 , 0.2923308644281944 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 682 , 0.1240602798626259 , 0.0229518395598165 , 0.2923308644281944 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 683 , 0.2923308644281944 , 0.5377051765895465 , 0.1240602798626259 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 684 , 0.0229518395598165 , 0.5377051765895465 , 0.1240602798626259 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 685 , 0.5377051765895465 , 0.2923308644281944 , 0.1240602798626259 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 686 , 0.0229518395598165 , 0.2923308644281944 , 0.1240602798626259 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 687 , 0.5377051765895465 , 0.0229518395598165 , 0.1240602798626259 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 688 , 0.2923308644281944 , 0.0229518395598165 , 0.1240602798626259 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 689 , 0.2923308644281944 , 0.5377051765895465 , 0.0229518395598165 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 690 , 0.1240602798626259 , 0.5377051765895465 , 0.0229518395598165 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 691 , 0.5377051765895465 , 0.2923308644281944 , 0.0229518395598165 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 692 , 0.1240602798626259 , 0.2923308644281944 , 0.0229518395598165 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 693 , 0.5377051765895465 , 0.1240602798626259 , 0.0229518395598165 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 694 , 0.2923308644281944 , 0.1240602798626259 , 0.0229518395598165 , 0.0229518395598165 , 0.0000375587155977 );
            ir->AddPentPoint( 695 , 0.0637915197105361 , 0.0107373897791662 , 0.1765900596203117 , 0.6850895111794498 , 0.0000288081136079 );
            ir->AddPentPoint( 696 , 0.0107373897791662 , 0.0637915197105361 , 0.1765900596203117 , 0.6850895111794498 , 0.0000288081136079 );
            ir->AddPentPoint( 697 , 0.0637915197105361 , 0.0637915197105361 , 0.1765900596203117 , 0.6850895111794498 , 0.0000288081136079 );
            ir->AddPentPoint( 698 , 0.0637915197105361 , 0.1765900596203117 , 0.0107373897791662 , 0.6850895111794498 , 0.0000288081136079 );
            ir->AddPentPoint( 699 , 0.1765900596203117 , 0.0637915197105361 , 0.0107373897791662 , 0.6850895111794498 , 0.0000288081136079 );
            ir->AddPentPoint( 700 , 0.0637915197105361 , 0.0637915197105361 , 0.0107373897791662 , 0.6850895111794498 , 0.0000288081136079 );
            ir->AddPentPoint( 701 , 0.0107373897791662 , 0.1765900596203117 , 0.0637915197105361 , 0.6850895111794498 , 0.0000288081136079 );
            ir->AddPentPoint( 702 , 0.0637915197105361 , 0.1765900596203117 , 0.0637915197105361 , 0.6850895111794498 , 0.0000288081136079 );
            ir->AddPentPoint( 703 , 0.1765900596203117 , 0.0107373897791662 , 0.0637915197105361 , 0.6850895111794498 , 0.0000288081136079 );
            ir->AddPentPoint( 704 , 0.0637915197105361 , 0.0107373897791662 , 0.0637915197105361 , 0.6850895111794498 , 0.0000288081136079 );
            ir->AddPentPoint( 705 , 0.1765900596203117 , 0.0637915197105361 , 0.0637915197105361 , 0.6850895111794498 , 0.0000288081136079 );
            ir->AddPentPoint( 706 , 0.0107373897791662 , 0.0637915197105361 , 0.0637915197105361 , 0.6850895111794498 , 0.0000288081136079 );
            ir->AddPentPoint( 707 , 0.0637915197105361 , 0.0107373897791662 , 0.6850895111794498 , 0.1765900596203117 , 0.0000288081136079 );
            ir->AddPentPoint( 708 , 0.0107373897791662 , 0.0637915197105361 , 0.6850895111794498 , 0.1765900596203117 , 0.0000288081136079 );
            ir->AddPentPoint( 709 , 0.0637915197105361 , 0.0637915197105361 , 0.6850895111794498 , 0.1765900596203117 , 0.0000288081136079 );
            ir->AddPentPoint( 710 , 0.0637915197105361 , 0.6850895111794498 , 0.0107373897791662 , 0.1765900596203117 , 0.0000288081136079 );
            ir->AddPentPoint( 711 , 0.6850895111794498 , 0.0637915197105361 , 0.0107373897791662 , 0.1765900596203117 , 0.0000288081136079 );
            ir->AddPentPoint( 712 , 0.0637915197105361 , 0.0637915197105361 , 0.0107373897791662 , 0.1765900596203117 , 0.0000288081136079 );
            ir->AddPentPoint( 713 , 0.0107373897791662 , 0.6850895111794498 , 0.0637915197105361 , 0.1765900596203117 , 0.0000288081136079 );
            ir->AddPentPoint( 714 , 0.0637915197105361 , 0.6850895111794498 , 0.0637915197105361 , 0.1765900596203117 , 0.0000288081136079 );
            ir->AddPentPoint( 715 , 0.6850895111794498 , 0.0107373897791662 , 0.0637915197105361 , 0.1765900596203117 , 0.0000288081136079 );
            ir->AddPentPoint( 716 , 0.0637915197105361 , 0.0107373897791662 , 0.0637915197105361 , 0.1765900596203117 , 0.0000288081136079 );
            ir->AddPentPoint( 717 , 0.6850895111794498 , 0.0637915197105361 , 0.0637915197105361 , 0.1765900596203117 , 0.0000288081136079 );
            ir->AddPentPoint( 718 , 0.0107373897791662 , 0.0637915197105361 , 0.0637915197105361 , 0.1765900596203117 , 0.0000288081136079 );
            ir->AddPentPoint( 719 , 0.0637915197105361 , 0.1765900596203117 , 0.6850895111794498 , 0.0107373897791662 , 0.0000288081136079 );
            ir->AddPentPoint( 720 , 0.1765900596203117 , 0.0637915197105361 , 0.6850895111794498 , 0.0107373897791662 , 0.0000288081136079 );
            ir->AddPentPoint( 721 , 0.0637915197105361 , 0.0637915197105361 , 0.6850895111794498 , 0.0107373897791662 , 0.0000288081136079 );
            ir->AddPentPoint( 722 , 0.0637915197105361 , 0.6850895111794498 , 0.1765900596203117 , 0.0107373897791662 , 0.0000288081136079 );
            ir->AddPentPoint( 723 , 0.6850895111794498 , 0.0637915197105361 , 0.1765900596203117 , 0.0107373897791662 , 0.0000288081136079 );
            ir->AddPentPoint( 724 , 0.0637915197105361 , 0.0637915197105361 , 0.1765900596203117 , 0.0107373897791662 , 0.0000288081136079 );
            ir->AddPentPoint( 725 , 0.1765900596203117 , 0.6850895111794498 , 0.0637915197105361 , 0.0107373897791662 , 0.0000288081136079 );
            ir->AddPentPoint( 726 , 0.0637915197105361 , 0.6850895111794498 , 0.0637915197105361 , 0.0107373897791662 , 0.0000288081136079 );
            ir->AddPentPoint( 727 , 0.6850895111794498 , 0.1765900596203117 , 0.0637915197105361 , 0.0107373897791662 , 0.0000288081136079 );
            ir->AddPentPoint( 728 , 0.0637915197105361 , 0.1765900596203117 , 0.0637915197105361 , 0.0107373897791662 , 0.0000288081136079 );
            ir->AddPentPoint( 729 , 0.6850895111794498 , 0.0637915197105361 , 0.0637915197105361 , 0.0107373897791662 , 0.0000288081136079 );
            ir->AddPentPoint( 730 , 0.1765900596203117 , 0.0637915197105361 , 0.0637915197105361 , 0.0107373897791662 , 0.0000288081136079 );
            ir->AddPentPoint( 731 , 0.0107373897791662 , 0.1765900596203117 , 0.6850895111794498 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 732 , 0.0637915197105361 , 0.1765900596203117 , 0.6850895111794498 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 733 , 0.1765900596203117 , 0.0107373897791662 , 0.6850895111794498 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 734 , 0.0637915197105361 , 0.0107373897791662 , 0.6850895111794498 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 735 , 0.1765900596203117 , 0.0637915197105361 , 0.6850895111794498 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 736 , 0.0107373897791662 , 0.0637915197105361 , 0.6850895111794498 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 737 , 0.0107373897791662 , 0.6850895111794498 , 0.1765900596203117 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 738 , 0.0637915197105361 , 0.6850895111794498 , 0.1765900596203117 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 739 , 0.6850895111794498 , 0.0107373897791662 , 0.1765900596203117 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 740 , 0.0637915197105361 , 0.0107373897791662 , 0.1765900596203117 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 741 , 0.6850895111794498 , 0.0637915197105361 , 0.1765900596203117 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 742 , 0.0107373897791662 , 0.0637915197105361 , 0.1765900596203117 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 743 , 0.1765900596203117 , 0.6850895111794498 , 0.0107373897791662 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 744 , 0.0637915197105361 , 0.6850895111794498 , 0.0107373897791662 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 745 , 0.6850895111794498 , 0.1765900596203117 , 0.0107373897791662 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 746 , 0.0637915197105361 , 0.1765900596203117 , 0.0107373897791662 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 747 , 0.6850895111794498 , 0.0637915197105361 , 0.0107373897791662 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 748 , 0.1765900596203117 , 0.0637915197105361 , 0.0107373897791662 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 749 , 0.1765900596203117 , 0.6850895111794498 , 0.0637915197105361 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 750 , 0.0107373897791662 , 0.6850895111794498 , 0.0637915197105361 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 751 , 0.6850895111794498 , 0.1765900596203117 , 0.0637915197105361 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 752 , 0.0107373897791662 , 0.1765900596203117 , 0.0637915197105361 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 753 , 0.6850895111794498 , 0.0107373897791662 , 0.0637915197105361 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 754 , 0.1765900596203117 , 0.0107373897791662 , 0.0637915197105361 , 0.0637915197105361 , 0.0000288081136079 );
            ir->AddPentPoint( 755 , 0.0668415703475457 , 0.0000000001038131 , 0.3396879598598828 , 0.5266288993412127 , 0.0000151165915450 );
            ir->AddPentPoint( 756 , 0.0000000001038131 , 0.0668415703475457 , 0.3396879598598828 , 0.5266288993412127 , 0.0000151165915450 );
            ir->AddPentPoint( 757 , 0.0668415703475457 , 0.0668415703475457 , 0.3396879598598828 , 0.5266288993412127 , 0.0000151165915450 );
            ir->AddPentPoint( 758 , 0.0668415703475457 , 0.3396879598598828 , 0.0000000001038131 , 0.5266288993412127 , 0.0000151165915450 );
            ir->AddPentPoint( 759 , 0.3396879598598828 , 0.0668415703475457 , 0.0000000001038131 , 0.5266288993412127 , 0.0000151165915450 );
            ir->AddPentPoint( 760 , 0.0668415703475457 , 0.0668415703475457 , 0.0000000001038131 , 0.5266288993412127 , 0.0000151165915450 );
            ir->AddPentPoint( 761 , 0.0000000001038131 , 0.3396879598598828 , 0.0668415703475457 , 0.5266288993412127 , 0.0000151165915450 );
            ir->AddPentPoint( 762 , 0.0668415703475457 , 0.3396879598598828 , 0.0668415703475457 , 0.5266288993412127 , 0.0000151165915450 );
            ir->AddPentPoint( 763 , 0.3396879598598828 , 0.0000000001038131 , 0.0668415703475457 , 0.5266288993412127 , 0.0000151165915450 );
            ir->AddPentPoint( 764 , 0.0668415703475457 , 0.0000000001038131 , 0.0668415703475457 , 0.5266288993412127 , 0.0000151165915450 );
            ir->AddPentPoint( 765 , 0.3396879598598828 , 0.0668415703475457 , 0.0668415703475457 , 0.5266288993412127 , 0.0000151165915450 );
            ir->AddPentPoint( 766 , 0.0000000001038131 , 0.0668415703475457 , 0.0668415703475457 , 0.5266288993412127 , 0.0000151165915450 );
            ir->AddPentPoint( 767 , 0.0668415703475457 , 0.0000000001038131 , 0.5266288993412127 , 0.3396879598598828 , 0.0000151165915450 );
            ir->AddPentPoint( 768 , 0.0000000001038131 , 0.0668415703475457 , 0.5266288993412127 , 0.3396879598598828 , 0.0000151165915450 );
            ir->AddPentPoint( 769 , 0.0668415703475457 , 0.0668415703475457 , 0.5266288993412127 , 0.3396879598598828 , 0.0000151165915450 );
            ir->AddPentPoint( 770 , 0.0668415703475457 , 0.5266288993412127 , 0.0000000001038131 , 0.3396879598598828 , 0.0000151165915450 );
            ir->AddPentPoint( 771 , 0.5266288993412127 , 0.0668415703475457 , 0.0000000001038131 , 0.3396879598598828 , 0.0000151165915450 );
            ir->AddPentPoint( 772 , 0.0668415703475457 , 0.0668415703475457 , 0.0000000001038131 , 0.3396879598598828 , 0.0000151165915450 );
            ir->AddPentPoint( 773 , 0.0000000001038131 , 0.5266288993412127 , 0.0668415703475457 , 0.3396879598598828 , 0.0000151165915450 );
            ir->AddPentPoint( 774 , 0.0668415703475457 , 0.5266288993412127 , 0.0668415703475457 , 0.3396879598598828 , 0.0000151165915450 );
            ir->AddPentPoint( 775 , 0.5266288993412127 , 0.0000000001038131 , 0.0668415703475457 , 0.3396879598598828 , 0.0000151165915450 );
            ir->AddPentPoint( 776 , 0.0668415703475457 , 0.0000000001038131 , 0.0668415703475457 , 0.3396879598598828 , 0.0000151165915450 );
            ir->AddPentPoint( 777 , 0.5266288993412127 , 0.0668415703475457 , 0.0668415703475457 , 0.3396879598598828 , 0.0000151165915450 );
            ir->AddPentPoint( 778 , 0.0000000001038131 , 0.0668415703475457 , 0.0668415703475457 , 0.3396879598598828 , 0.0000151165915450 );
            ir->AddPentPoint( 779 , 0.0668415703475457 , 0.3396879598598828 , 0.5266288993412127 , 0.0000000001038131 , 0.0000151165915450 );
            ir->AddPentPoint( 780 , 0.3396879598598828 , 0.0668415703475457 , 0.5266288993412127 , 0.0000000001038131 , 0.0000151165915450 );
            ir->AddPentPoint( 781 , 0.0668415703475457 , 0.0668415703475457 , 0.5266288993412127 , 0.0000000001038131 , 0.0000151165915450 );
            ir->AddPentPoint( 782 , 0.0668415703475457 , 0.5266288993412127 , 0.3396879598598828 , 0.0000000001038131 , 0.0000151165915450 );
            ir->AddPentPoint( 783 , 0.5266288993412127 , 0.0668415703475457 , 0.3396879598598828 , 0.0000000001038131 , 0.0000151165915450 );
            ir->AddPentPoint( 784 , 0.0668415703475457 , 0.0668415703475457 , 0.3396879598598828 , 0.0000000001038131 , 0.0000151165915450 );
            ir->AddPentPoint( 785 , 0.3396879598598828 , 0.5266288993412127 , 0.0668415703475457 , 0.0000000001038131 , 0.0000151165915450 );
            ir->AddPentPoint( 786 , 0.0668415703475457 , 0.5266288993412127 , 0.0668415703475457 , 0.0000000001038131 , 0.0000151165915450 );
            ir->AddPentPoint( 787 , 0.5266288993412127 , 0.3396879598598828 , 0.0668415703475457 , 0.0000000001038131 , 0.0000151165915450 );
            ir->AddPentPoint( 788 , 0.0668415703475457 , 0.3396879598598828 , 0.0668415703475457 , 0.0000000001038131 , 0.0000151165915450 );
            ir->AddPentPoint( 789 , 0.5266288993412127 , 0.0668415703475457 , 0.0668415703475457 , 0.0000000001038131 , 0.0000151165915450 );
            ir->AddPentPoint( 790 , 0.3396879598598828 , 0.0668415703475457 , 0.0668415703475457 , 0.0000000001038131 , 0.0000151165915450 );
            ir->AddPentPoint( 791 , 0.0000000001038131 , 0.3396879598598828 , 0.5266288993412127 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 792 , 0.0668415703475457 , 0.3396879598598828 , 0.5266288993412127 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 793 , 0.3396879598598828 , 0.0000000001038131 , 0.5266288993412127 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 794 , 0.0668415703475457 , 0.0000000001038131 , 0.5266288993412127 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 795 , 0.3396879598598828 , 0.0668415703475457 , 0.5266288993412127 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 796 , 0.0000000001038131 , 0.0668415703475457 , 0.5266288993412127 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 797 , 0.0000000001038131 , 0.5266288993412127 , 0.3396879598598828 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 798 , 0.0668415703475457 , 0.5266288993412127 , 0.3396879598598828 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 799 , 0.5266288993412127 , 0.0000000001038131 , 0.3396879598598828 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 800 , 0.0668415703475457 , 0.0000000001038131 , 0.3396879598598828 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 801 , 0.5266288993412127 , 0.0668415703475457 , 0.3396879598598828 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 802 , 0.0000000001038131 , 0.0668415703475457 , 0.3396879598598828 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 803 , 0.3396879598598828 , 0.5266288993412127 , 0.0000000001038131 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 804 , 0.0668415703475457 , 0.5266288993412127 , 0.0000000001038131 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 805 , 0.5266288993412127 , 0.3396879598598828 , 0.0000000001038131 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 806 , 0.0668415703475457 , 0.3396879598598828 , 0.0000000001038131 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 807 , 0.5266288993412127 , 0.0668415703475457 , 0.0000000001038131 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 808 , 0.3396879598598828 , 0.0668415703475457 , 0.0000000001038131 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 809 , 0.3396879598598828 , 0.5266288993412127 , 0.0668415703475457 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 810 , 0.0000000001038131 , 0.5266288993412127 , 0.0668415703475457 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 811 , 0.5266288993412127 , 0.3396879598598828 , 0.0668415703475457 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 812 , 0.0000000001038131 , 0.3396879598598828 , 0.0668415703475457 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 813 , 0.5266288993412127 , 0.0000000001038131 , 0.0668415703475457 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 814 , 0.3396879598598828 , 0.0000000001038131 , 0.0668415703475457 , 0.0668415703475457 , 0.0000151165915450 );
            ir->AddPentPoint( 815 , 0.1153852675306786 , 0.0334572708282863 , 0.2633720195447322 , 0.4724001745656243 , 0.0000855167401672 );
            ir->AddPentPoint( 816 , 0.0334572708282863 , 0.1153852675306786 , 0.2633720195447322 , 0.4724001745656243 , 0.0000855167401672 );
            ir->AddPentPoint( 817 , 0.1153852675306786 , 0.1153852675306786 , 0.2633720195447322 , 0.4724001745656243 , 0.0000855167401672 );
            ir->AddPentPoint( 818 , 0.1153852675306786 , 0.2633720195447322 , 0.0334572708282863 , 0.4724001745656243 , 0.0000855167401672 );
            ir->AddPentPoint( 819 , 0.2633720195447322 , 0.1153852675306786 , 0.0334572708282863 , 0.4724001745656243 , 0.0000855167401672 );
            ir->AddPentPoint( 820 , 0.1153852675306786 , 0.1153852675306786 , 0.0334572708282863 , 0.4724001745656243 , 0.0000855167401672 );
            ir->AddPentPoint( 821 , 0.0334572708282863 , 0.2633720195447322 , 0.1153852675306786 , 0.4724001745656243 , 0.0000855167401672 );
            ir->AddPentPoint( 822 , 0.1153852675306786 , 0.2633720195447322 , 0.1153852675306786 , 0.4724001745656243 , 0.0000855167401672 );
            ir->AddPentPoint( 823 , 0.2633720195447322 , 0.0334572708282863 , 0.1153852675306786 , 0.4724001745656243 , 0.0000855167401672 );
            ir->AddPentPoint( 824 , 0.1153852675306786 , 0.0334572708282863 , 0.1153852675306786 , 0.4724001745656243 , 0.0000855167401672 );
            ir->AddPentPoint( 825 , 0.2633720195447322 , 0.1153852675306786 , 0.1153852675306786 , 0.4724001745656243 , 0.0000855167401672 );
            ir->AddPentPoint( 826 , 0.0334572708282863 , 0.1153852675306786 , 0.1153852675306786 , 0.4724001745656243 , 0.0000855167401672 );
            ir->AddPentPoint( 827 , 0.1153852675306786 , 0.0334572708282863 , 0.4724001745656243 , 0.2633720195447322 , 0.0000855167401672 );
            ir->AddPentPoint( 828 , 0.0334572708282863 , 0.1153852675306786 , 0.4724001745656243 , 0.2633720195447322 , 0.0000855167401672 );
            ir->AddPentPoint( 829 , 0.1153852675306786 , 0.1153852675306786 , 0.4724001745656243 , 0.2633720195447322 , 0.0000855167401672 );
            ir->AddPentPoint( 830 , 0.1153852675306786 , 0.4724001745656243 , 0.0334572708282863 , 0.2633720195447322 , 0.0000855167401672 );
            ir->AddPentPoint( 831 , 0.4724001745656243 , 0.1153852675306786 , 0.0334572708282863 , 0.2633720195447322 , 0.0000855167401672 );
            ir->AddPentPoint( 832 , 0.1153852675306786 , 0.1153852675306786 , 0.0334572708282863 , 0.2633720195447322 , 0.0000855167401672 );
            ir->AddPentPoint( 833 , 0.0334572708282863 , 0.4724001745656243 , 0.1153852675306786 , 0.2633720195447322 , 0.0000855167401672 );
            ir->AddPentPoint( 834 , 0.1153852675306786 , 0.4724001745656243 , 0.1153852675306786 , 0.2633720195447322 , 0.0000855167401672 );
            ir->AddPentPoint( 835 , 0.4724001745656243 , 0.0334572708282863 , 0.1153852675306786 , 0.2633720195447322 , 0.0000855167401672 );
            ir->AddPentPoint( 836 , 0.1153852675306786 , 0.0334572708282863 , 0.1153852675306786 , 0.2633720195447322 , 0.0000855167401672 );
            ir->AddPentPoint( 837 , 0.4724001745656243 , 0.1153852675306786 , 0.1153852675306786 , 0.2633720195447322 , 0.0000855167401672 );
            ir->AddPentPoint( 838 , 0.0334572708282863 , 0.1153852675306786 , 0.1153852675306786 , 0.2633720195447322 , 0.0000855167401672 );
            ir->AddPentPoint( 839 , 0.1153852675306786 , 0.2633720195447322 , 0.4724001745656243 , 0.0334572708282863 , 0.0000855167401672 );
            ir->AddPentPoint( 840 , 0.2633720195447322 , 0.1153852675306786 , 0.4724001745656243 , 0.0334572708282863 , 0.0000855167401672 );
            ir->AddPentPoint( 841 , 0.1153852675306786 , 0.1153852675306786 , 0.4724001745656243 , 0.0334572708282863 , 0.0000855167401672 );
            ir->AddPentPoint( 842 , 0.1153852675306786 , 0.4724001745656243 , 0.2633720195447322 , 0.0334572708282863 , 0.0000855167401672 );
            ir->AddPentPoint( 843 , 0.4724001745656243 , 0.1153852675306786 , 0.2633720195447322 , 0.0334572708282863 , 0.0000855167401672 );
            ir->AddPentPoint( 844 , 0.1153852675306786 , 0.1153852675306786 , 0.2633720195447322 , 0.0334572708282863 , 0.0000855167401672 );
            ir->AddPentPoint( 845 , 0.2633720195447322 , 0.4724001745656243 , 0.1153852675306786 , 0.0334572708282863 , 0.0000855167401672 );
            ir->AddPentPoint( 846 , 0.1153852675306786 , 0.4724001745656243 , 0.1153852675306786 , 0.0334572708282863 , 0.0000855167401672 );
            ir->AddPentPoint( 847 , 0.4724001745656243 , 0.2633720195447322 , 0.1153852675306786 , 0.0334572708282863 , 0.0000855167401672 );
            ir->AddPentPoint( 848 , 0.1153852675306786 , 0.2633720195447322 , 0.1153852675306786 , 0.0334572708282863 , 0.0000855167401672 );
            ir->AddPentPoint( 849 , 0.4724001745656243 , 0.1153852675306786 , 0.1153852675306786 , 0.0334572708282863 , 0.0000855167401672 );
            ir->AddPentPoint( 850 , 0.2633720195447322 , 0.1153852675306786 , 0.1153852675306786 , 0.0334572708282863 , 0.0000855167401672 );
            ir->AddPentPoint( 851 , 0.0334572708282863 , 0.2633720195447322 , 0.4724001745656243 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 852 , 0.1153852675306786 , 0.2633720195447322 , 0.4724001745656243 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 853 , 0.2633720195447322 , 0.0334572708282863 , 0.4724001745656243 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 854 , 0.1153852675306786 , 0.0334572708282863 , 0.4724001745656243 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 855 , 0.2633720195447322 , 0.1153852675306786 , 0.4724001745656243 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 856 , 0.0334572708282863 , 0.1153852675306786 , 0.4724001745656243 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 857 , 0.0334572708282863 , 0.4724001745656243 , 0.2633720195447322 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 858 , 0.1153852675306786 , 0.4724001745656243 , 0.2633720195447322 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 859 , 0.4724001745656243 , 0.0334572708282863 , 0.2633720195447322 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 860 , 0.1153852675306786 , 0.0334572708282863 , 0.2633720195447322 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 861 , 0.4724001745656243 , 0.1153852675306786 , 0.2633720195447322 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 862 , 0.0334572708282863 , 0.1153852675306786 , 0.2633720195447322 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 863 , 0.2633720195447322 , 0.4724001745656243 , 0.0334572708282863 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 864 , 0.1153852675306786 , 0.4724001745656243 , 0.0334572708282863 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 865 , 0.4724001745656243 , 0.2633720195447322 , 0.0334572708282863 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 866 , 0.1153852675306786 , 0.2633720195447322 , 0.0334572708282863 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 867 , 0.4724001745656243 , 0.1153852675306786 , 0.0334572708282863 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 868 , 0.2633720195447322 , 0.1153852675306786 , 0.0334572708282863 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 869 , 0.2633720195447322 , 0.4724001745656243 , 0.1153852675306786 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 870 , 0.0334572708282863 , 0.4724001745656243 , 0.1153852675306786 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 871 , 0.4724001745656243 , 0.2633720195447322 , 0.1153852675306786 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 872 , 0.0334572708282863 , 0.2633720195447322 , 0.1153852675306786 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 873 , 0.4724001745656243 , 0.0334572708282863 , 0.1153852675306786 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 874 , 0.2633720195447322 , 0.0334572708282863 , 0.1153852675306786 , 0.1153852675306786 , 0.0000855167401672 );
            ir->AddPentPoint( 875 , 0.1987579365227028 , 0.0000000000000000 , 0.0747748152891036 , 0.5277093116654908 , 0.0000227107898376 );
            ir->AddPentPoint( 876 , 0.0000000000000000 , 0.1987579365227028 , 0.0747748152891036 , 0.5277093116654908 , 0.0000227107898376 );
            ir->AddPentPoint( 877 , 0.1987579365227028 , 0.1987579365227028 , 0.0747748152891036 , 0.5277093116654908 , 0.0000227107898376 );
            ir->AddPentPoint( 878 , 0.1987579365227028 , 0.0747748152891036 , 0.0000000000000000 , 0.5277093116654908 , 0.0000227107898376 );
            ir->AddPentPoint( 879 , 0.0747748152891036 , 0.1987579365227028 , 0.0000000000000000 , 0.5277093116654908 , 0.0000227107898376 );
            ir->AddPentPoint( 880 , 0.1987579365227028 , 0.1987579365227028 , 0.0000000000000000 , 0.5277093116654908 , 0.0000227107898376 );
            ir->AddPentPoint( 881 , 0.0000000000000000 , 0.0747748152891036 , 0.1987579365227028 , 0.5277093116654908 , 0.0000227107898376 );
            ir->AddPentPoint( 882 , 0.1987579365227028 , 0.0747748152891036 , 0.1987579365227028 , 0.5277093116654908 , 0.0000227107898376 );
            ir->AddPentPoint( 883 , 0.0747748152891036 , 0.0000000000000000 , 0.1987579365227028 , 0.5277093116654908 , 0.0000227107898376 );
            ir->AddPentPoint( 884 , 0.1987579365227028 , 0.0000000000000000 , 0.1987579365227028 , 0.5277093116654908 , 0.0000227107898376 );
            ir->AddPentPoint( 885 , 0.0747748152891036 , 0.1987579365227028 , 0.1987579365227028 , 0.5277093116654908 , 0.0000227107898376 );
            ir->AddPentPoint( 886 , 0.0000000000000000 , 0.1987579365227028 , 0.1987579365227028 , 0.5277093116654908 , 0.0000227107898376 );
            ir->AddPentPoint( 887 , 0.1987579365227028 , 0.0000000000000000 , 0.5277093116654908 , 0.0747748152891036 , 0.0000227107898376 );
            ir->AddPentPoint( 888 , 0.0000000000000000 , 0.1987579365227028 , 0.5277093116654908 , 0.0747748152891036 , 0.0000227107898376 );
            ir->AddPentPoint( 889 , 0.1987579365227028 , 0.1987579365227028 , 0.5277093116654908 , 0.0747748152891036 , 0.0000227107898376 );
            ir->AddPentPoint( 890 , 0.1987579365227028 , 0.5277093116654908 , 0.0000000000000000 , 0.0747748152891036 , 0.0000227107898376 );
            ir->AddPentPoint( 891 , 0.5277093116654908 , 0.1987579365227028 , 0.0000000000000000 , 0.0747748152891036 , 0.0000227107898376 );
            ir->AddPentPoint( 892 , 0.1987579365227028 , 0.1987579365227028 , 0.0000000000000000 , 0.0747748152891036 , 0.0000227107898376 );
            ir->AddPentPoint( 893 , 0.0000000000000000 , 0.5277093116654908 , 0.1987579365227028 , 0.0747748152891036 , 0.0000227107898376 );
            ir->AddPentPoint( 894 , 0.1987579365227028 , 0.5277093116654908 , 0.1987579365227028 , 0.0747748152891036 , 0.0000227107898376 );
            ir->AddPentPoint( 895 , 0.5277093116654908 , 0.0000000000000000 , 0.1987579365227028 , 0.0747748152891036 , 0.0000227107898376 );
            ir->AddPentPoint( 896 , 0.1987579365227028 , 0.0000000000000000 , 0.1987579365227028 , 0.0747748152891036 , 0.0000227107898376 );
            ir->AddPentPoint( 897 , 0.5277093116654908 , 0.1987579365227028 , 0.1987579365227028 , 0.0747748152891036 , 0.0000227107898376 );
            ir->AddPentPoint( 898 , 0.0000000000000000 , 0.1987579365227028 , 0.1987579365227028 , 0.0747748152891036 , 0.0000227107898376 );
            ir->AddPentPoint( 899 , 0.1987579365227028 , 0.0747748152891036 , 0.5277093116654908 , 0.0000000000000000 , 0.0000227107898376 );
            ir->AddPentPoint( 900 , 0.0747748152891036 , 0.1987579365227028 , 0.5277093116654908 , 0.0000000000000000 , 0.0000227107898376 );
            ir->AddPentPoint( 901 , 0.1987579365227028 , 0.1987579365227028 , 0.5277093116654908 , 0.0000000000000000 , 0.0000227107898376 );
            ir->AddPentPoint( 902 , 0.1987579365227028 , 0.5277093116654908 , 0.0747748152891036 , 0.0000000000000000 , 0.0000227107898376 );
            ir->AddPentPoint( 903 , 0.5277093116654908 , 0.1987579365227028 , 0.0747748152891036 , 0.0000000000000000 , 0.0000227107898376 );
            ir->AddPentPoint( 904 , 0.1987579365227028 , 0.1987579365227028 , 0.0747748152891036 , 0.0000000000000000 , 0.0000227107898376 );
            ir->AddPentPoint( 905 , 0.0747748152891036 , 0.5277093116654908 , 0.1987579365227028 , 0.0000000000000000 , 0.0000227107898376 );
            ir->AddPentPoint( 906 , 0.1987579365227028 , 0.5277093116654908 , 0.1987579365227028 , 0.0000000000000000 , 0.0000227107898376 );
            ir->AddPentPoint( 907 , 0.5277093116654908 , 0.0747748152891036 , 0.1987579365227028 , 0.0000000000000000 , 0.0000227107898376 );
            ir->AddPentPoint( 908 , 0.1987579365227028 , 0.0747748152891036 , 0.1987579365227028 , 0.0000000000000000 , 0.0000227107898376 );
            ir->AddPentPoint( 909 , 0.5277093116654908 , 0.1987579365227028 , 0.1987579365227028 , 0.0000000000000000 , 0.0000227107898376 );
            ir->AddPentPoint( 910 , 0.0747748152891036 , 0.1987579365227028 , 0.1987579365227028 , 0.0000000000000000 , 0.0000227107898376 );
            ir->AddPentPoint( 911 , 0.0000000000000000 , 0.0747748152891036 , 0.5277093116654908 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 912 , 0.1987579365227028 , 0.0747748152891036 , 0.5277093116654908 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 913 , 0.0747748152891036 , 0.0000000000000000 , 0.5277093116654908 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 914 , 0.1987579365227028 , 0.0000000000000000 , 0.5277093116654908 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 915 , 0.0747748152891036 , 0.1987579365227028 , 0.5277093116654908 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 916 , 0.0000000000000000 , 0.1987579365227028 , 0.5277093116654908 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 917 , 0.0000000000000000 , 0.5277093116654908 , 0.0747748152891036 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 918 , 0.1987579365227028 , 0.5277093116654908 , 0.0747748152891036 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 919 , 0.5277093116654908 , 0.0000000000000000 , 0.0747748152891036 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 920 , 0.1987579365227028 , 0.0000000000000000 , 0.0747748152891036 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 921 , 0.5277093116654908 , 0.1987579365227028 , 0.0747748152891036 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 922 , 0.0000000000000000 , 0.1987579365227028 , 0.0747748152891036 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 923 , 0.0747748152891036 , 0.5277093116654908 , 0.0000000000000000 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 924 , 0.1987579365227028 , 0.5277093116654908 , 0.0000000000000000 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 925 , 0.5277093116654908 , 0.0747748152891036 , 0.0000000000000000 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 926 , 0.1987579365227028 , 0.0747748152891036 , 0.0000000000000000 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 927 , 0.5277093116654908 , 0.1987579365227028 , 0.0000000000000000 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 928 , 0.0747748152891036 , 0.1987579365227028 , 0.0000000000000000 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 929 , 0.0747748152891036 , 0.5277093116654908 , 0.1987579365227028 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 930 , 0.0000000000000000 , 0.5277093116654908 , 0.1987579365227028 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 931 , 0.5277093116654908 , 0.0747748152891036 , 0.1987579365227028 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 932 , 0.0000000000000000 , 0.0747748152891036 , 0.1987579365227028 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 933 , 0.5277093116654908 , 0.0000000000000000 , 0.1987579365227028 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 934 , 0.0747748152891036 , 0.0000000000000000 , 0.1987579365227028 , 0.1987579365227028 , 0.0000227107898376 );
            ir->AddPentPoint( 935 , 0.2331750657935490 , 0.0129306209538988 , 0.1473703181830719 , 0.3733489292759313 , 0.0000454609731552 );
            ir->AddPentPoint( 936 , 0.0129306209538988 , 0.2331750657935490 , 0.1473703181830719 , 0.3733489292759313 , 0.0000454609731552 );
            ir->AddPentPoint( 937 , 0.2331750657935490 , 0.2331750657935490 , 0.1473703181830719 , 0.3733489292759313 , 0.0000454609731552 );
            ir->AddPentPoint( 938 , 0.2331750657935490 , 0.1473703181830719 , 0.0129306209538988 , 0.3733489292759313 , 0.0000454609731552 );
            ir->AddPentPoint( 939 , 0.1473703181830719 , 0.2331750657935490 , 0.0129306209538988 , 0.3733489292759313 , 0.0000454609731552 );
            ir->AddPentPoint( 940 , 0.2331750657935490 , 0.2331750657935490 , 0.0129306209538988 , 0.3733489292759313 , 0.0000454609731552 );
            ir->AddPentPoint( 941 , 0.0129306209538988 , 0.1473703181830719 , 0.2331750657935490 , 0.3733489292759313 , 0.0000454609731552 );
            ir->AddPentPoint( 942 , 0.2331750657935490 , 0.1473703181830719 , 0.2331750657935490 , 0.3733489292759313 , 0.0000454609731552 );
            ir->AddPentPoint( 943 , 0.1473703181830719 , 0.0129306209538988 , 0.2331750657935490 , 0.3733489292759313 , 0.0000454609731552 );
            ir->AddPentPoint( 944 , 0.2331750657935490 , 0.0129306209538988 , 0.2331750657935490 , 0.3733489292759313 , 0.0000454609731552 );
            ir->AddPentPoint( 945 , 0.1473703181830719 , 0.2331750657935490 , 0.2331750657935490 , 0.3733489292759313 , 0.0000454609731552 );
            ir->AddPentPoint( 946 , 0.0129306209538988 , 0.2331750657935490 , 0.2331750657935490 , 0.3733489292759313 , 0.0000454609731552 );
            ir->AddPentPoint( 947 , 0.2331750657935490 , 0.0129306209538988 , 0.3733489292759313 , 0.1473703181830719 , 0.0000454609731552 );
            ir->AddPentPoint( 948 , 0.0129306209538988 , 0.2331750657935490 , 0.3733489292759313 , 0.1473703181830719 , 0.0000454609731552 );
            ir->AddPentPoint( 949 , 0.2331750657935490 , 0.2331750657935490 , 0.3733489292759313 , 0.1473703181830719 , 0.0000454609731552 );
            ir->AddPentPoint( 950 , 0.2331750657935490 , 0.3733489292759313 , 0.0129306209538988 , 0.1473703181830719 , 0.0000454609731552 );
            ir->AddPentPoint( 951 , 0.3733489292759313 , 0.2331750657935490 , 0.0129306209538988 , 0.1473703181830719 , 0.0000454609731552 );
            ir->AddPentPoint( 952 , 0.2331750657935490 , 0.2331750657935490 , 0.0129306209538988 , 0.1473703181830719 , 0.0000454609731552 );
            ir->AddPentPoint( 953 , 0.0129306209538988 , 0.3733489292759313 , 0.2331750657935490 , 0.1473703181830719 , 0.0000454609731552 );
            ir->AddPentPoint( 954 , 0.2331750657935490 , 0.3733489292759313 , 0.2331750657935490 , 0.1473703181830719 , 0.0000454609731552 );
            ir->AddPentPoint( 955 , 0.3733489292759313 , 0.0129306209538988 , 0.2331750657935490 , 0.1473703181830719 , 0.0000454609731552 );
            ir->AddPentPoint( 956 , 0.2331750657935490 , 0.0129306209538988 , 0.2331750657935490 , 0.1473703181830719 , 0.0000454609731552 );
            ir->AddPentPoint( 957 , 0.3733489292759313 , 0.2331750657935490 , 0.2331750657935490 , 0.1473703181830719 , 0.0000454609731552 );
            ir->AddPentPoint( 958 , 0.0129306209538988 , 0.2331750657935490 , 0.2331750657935490 , 0.1473703181830719 , 0.0000454609731552 );
            ir->AddPentPoint( 959 , 0.2331750657935490 , 0.1473703181830719 , 0.3733489292759313 , 0.0129306209538988 , 0.0000454609731552 );
            ir->AddPentPoint( 960 , 0.1473703181830719 , 0.2331750657935490 , 0.3733489292759313 , 0.0129306209538988 , 0.0000454609731552 );
            ir->AddPentPoint( 961 , 0.2331750657935490 , 0.2331750657935490 , 0.3733489292759313 , 0.0129306209538988 , 0.0000454609731552 );
            ir->AddPentPoint( 962 , 0.2331750657935490 , 0.3733489292759313 , 0.1473703181830719 , 0.0129306209538988 , 0.0000454609731552 );
            ir->AddPentPoint( 963 , 0.3733489292759313 , 0.2331750657935490 , 0.1473703181830719 , 0.0129306209538988 , 0.0000454609731552 );
            ir->AddPentPoint( 964 , 0.2331750657935490 , 0.2331750657935490 , 0.1473703181830719 , 0.0129306209538988 , 0.0000454609731552 );
            ir->AddPentPoint( 965 , 0.1473703181830719 , 0.3733489292759313 , 0.2331750657935490 , 0.0129306209538988 , 0.0000454609731552 );
            ir->AddPentPoint( 966 , 0.2331750657935490 , 0.3733489292759313 , 0.2331750657935490 , 0.0129306209538988 , 0.0000454609731552 );
            ir->AddPentPoint( 967 , 0.3733489292759313 , 0.1473703181830719 , 0.2331750657935490 , 0.0129306209538988 , 0.0000454609731552 );
            ir->AddPentPoint( 968 , 0.2331750657935490 , 0.1473703181830719 , 0.2331750657935490 , 0.0129306209538988 , 0.0000454609731552 );
            ir->AddPentPoint( 969 , 0.3733489292759313 , 0.2331750657935490 , 0.2331750657935490 , 0.0129306209538988 , 0.0000454609731552 );
            ir->AddPentPoint( 970 , 0.1473703181830719 , 0.2331750657935490 , 0.2331750657935490 , 0.0129306209538988 , 0.0000454609731552 );
            ir->AddPentPoint( 971 , 0.0129306209538988 , 0.1473703181830719 , 0.3733489292759313 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 972 , 0.2331750657935490 , 0.1473703181830719 , 0.3733489292759313 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 973 , 0.1473703181830719 , 0.0129306209538988 , 0.3733489292759313 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 974 , 0.2331750657935490 , 0.0129306209538988 , 0.3733489292759313 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 975 , 0.1473703181830719 , 0.2331750657935490 , 0.3733489292759313 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 976 , 0.0129306209538988 , 0.2331750657935490 , 0.3733489292759313 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 977 , 0.0129306209538988 , 0.3733489292759313 , 0.1473703181830719 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 978 , 0.2331750657935490 , 0.3733489292759313 , 0.1473703181830719 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 979 , 0.3733489292759313 , 0.0129306209538988 , 0.1473703181830719 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 980 , 0.2331750657935490 , 0.0129306209538988 , 0.1473703181830719 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 981 , 0.3733489292759313 , 0.2331750657935490 , 0.1473703181830719 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 982 , 0.0129306209538988 , 0.2331750657935490 , 0.1473703181830719 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 983 , 0.1473703181830719 , 0.3733489292759313 , 0.0129306209538988 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 984 , 0.2331750657935490 , 0.3733489292759313 , 0.0129306209538988 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 985 , 0.3733489292759313 , 0.1473703181830719 , 0.0129306209538988 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 986 , 0.2331750657935490 , 0.1473703181830719 , 0.0129306209538988 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 987 , 0.3733489292759313 , 0.2331750657935490 , 0.0129306209538988 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 988 , 0.1473703181830719 , 0.2331750657935490 , 0.0129306209538988 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 989 , 0.1473703181830719 , 0.3733489292759313 , 0.2331750657935490 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 990 , 0.0129306209538988 , 0.3733489292759313 , 0.2331750657935490 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 991 , 0.3733489292759313 , 0.1473703181830719 , 0.2331750657935490 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 992 , 0.0129306209538988 , 0.1473703181830719 , 0.2331750657935490 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 993 , 0.3733489292759313 , 0.0129306209538988 , 0.2331750657935490 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 994 , 0.1473703181830719 , 0.0129306209538988 , 0.2331750657935490 , 0.2331750657935490 , 0.0000454609731552 );
            ir->AddPentPoint( 995 , 0.3679171564486097 , 0.0103703935710663 , 0.0786829086131733 , 0.1751123849185410 , 0.0000456325850473 );
            ir->AddPentPoint( 996 , 0.0103703935710663 , 0.3679171564486097 , 0.0786829086131733 , 0.1751123849185410 , 0.0000456325850473 );
            ir->AddPentPoint( 997 , 0.3679171564486097 , 0.3679171564486097 , 0.0786829086131733 , 0.1751123849185410 , 0.0000456325850473 );
            ir->AddPentPoint( 998 , 0.3679171564486097 , 0.0786829086131733 , 0.0103703935710663 , 0.1751123849185410 , 0.0000456325850473 );
            ir->AddPentPoint( 999 , 0.0786829086131733 , 0.3679171564486097 , 0.0103703935710663 , 0.1751123849185410 , 0.0000456325850473 );
            ir->AddPentPoint( 1000 , 0.3679171564486097 , 0.3679171564486097 , 0.0103703935710663 , 0.1751123849185410 , 0.0000456325850473 );
            ir->AddPentPoint( 1001 , 0.0103703935710663 , 0.0786829086131733 , 0.3679171564486097 , 0.1751123849185410 , 0.0000456325850473 );
            ir->AddPentPoint( 1002 , 0.3679171564486097 , 0.0786829086131733 , 0.3679171564486097 , 0.1751123849185410 , 0.0000456325850473 );
            ir->AddPentPoint( 1003 , 0.0786829086131733 , 0.0103703935710663 , 0.3679171564486097 , 0.1751123849185410 , 0.0000456325850473 );
            ir->AddPentPoint( 1004 , 0.3679171564486097 , 0.0103703935710663 , 0.3679171564486097 , 0.1751123849185410 , 0.0000456325850473 );
            ir->AddPentPoint( 1005 , 0.0786829086131733 , 0.3679171564486097 , 0.3679171564486097 , 0.1751123849185410 , 0.0000456325850473 );
            ir->AddPentPoint( 1006 , 0.0103703935710663 , 0.3679171564486097 , 0.3679171564486097 , 0.1751123849185410 , 0.0000456325850473 );
            ir->AddPentPoint( 1007 , 0.3679171564486097 , 0.0103703935710663 , 0.1751123849185410 , 0.0786829086131733 , 0.0000456325850473 );
            ir->AddPentPoint( 1008 , 0.0103703935710663 , 0.3679171564486097 , 0.1751123849185410 , 0.0786829086131733 , 0.0000456325850473 );
            ir->AddPentPoint( 1009 , 0.3679171564486097 , 0.3679171564486097 , 0.1751123849185410 , 0.0786829086131733 , 0.0000456325850473 );
            ir->AddPentPoint( 1010 , 0.3679171564486097 , 0.1751123849185410 , 0.0103703935710663 , 0.0786829086131733 , 0.0000456325850473 );
            ir->AddPentPoint( 1011 , 0.1751123849185410 , 0.3679171564486097 , 0.0103703935710663 , 0.0786829086131733 , 0.0000456325850473 );
            ir->AddPentPoint( 1012 , 0.3679171564486097 , 0.3679171564486097 , 0.0103703935710663 , 0.0786829086131733 , 0.0000456325850473 );
            ir->AddPentPoint( 1013 , 0.0103703935710663 , 0.1751123849185410 , 0.3679171564486097 , 0.0786829086131733 , 0.0000456325850473 );
            ir->AddPentPoint( 1014 , 0.3679171564486097 , 0.1751123849185410 , 0.3679171564486097 , 0.0786829086131733 , 0.0000456325850473 );
            ir->AddPentPoint( 1015 , 0.1751123849185410 , 0.0103703935710663 , 0.3679171564486097 , 0.0786829086131733 , 0.0000456325850473 );
            ir->AddPentPoint( 1016 , 0.3679171564486097 , 0.0103703935710663 , 0.3679171564486097 , 0.0786829086131733 , 0.0000456325850473 );
            ir->AddPentPoint( 1017 , 0.1751123849185410 , 0.3679171564486097 , 0.3679171564486097 , 0.0786829086131733 , 0.0000456325850473 );
            ir->AddPentPoint( 1018 , 0.0103703935710663 , 0.3679171564486097 , 0.3679171564486097 , 0.0786829086131733 , 0.0000456325850473 );
            ir->AddPentPoint( 1019 , 0.3679171564486097 , 0.0786829086131733 , 0.1751123849185410 , 0.0103703935710663 , 0.0000456325850473 );
            ir->AddPentPoint( 1020 , 0.0786829086131733 , 0.3679171564486097 , 0.1751123849185410 , 0.0103703935710663 , 0.0000456325850473 );
            ir->AddPentPoint( 1021 , 0.3679171564486097 , 0.3679171564486097 , 0.1751123849185410 , 0.0103703935710663 , 0.0000456325850473 );
            ir->AddPentPoint( 1022 , 0.3679171564486097 , 0.1751123849185410 , 0.0786829086131733 , 0.0103703935710663 , 0.0000456325850473 );
            ir->AddPentPoint( 1023 , 0.1751123849185410 , 0.3679171564486097 , 0.0786829086131733 , 0.0103703935710663 , 0.0000456325850473 );
            ir->AddPentPoint( 1024 , 0.3679171564486097 , 0.3679171564486097 , 0.0786829086131733 , 0.0103703935710663 , 0.0000456325850473 );
            ir->AddPentPoint( 1025 , 0.0786829086131733 , 0.1751123849185410 , 0.3679171564486097 , 0.0103703935710663 , 0.0000456325850473 );
            ir->AddPentPoint( 1026 , 0.3679171564486097 , 0.1751123849185410 , 0.3679171564486097 , 0.0103703935710663 , 0.0000456325850473 );
            ir->AddPentPoint( 1027 , 0.1751123849185410 , 0.0786829086131733 , 0.3679171564486097 , 0.0103703935710663 , 0.0000456325850473 );
            ir->AddPentPoint( 1028 , 0.3679171564486097 , 0.0786829086131733 , 0.3679171564486097 , 0.0103703935710663 , 0.0000456325850473 );
            ir->AddPentPoint( 1029 , 0.1751123849185410 , 0.3679171564486097 , 0.3679171564486097 , 0.0103703935710663 , 0.0000456325850473 );
            ir->AddPentPoint( 1030 , 0.0786829086131733 , 0.3679171564486097 , 0.3679171564486097 , 0.0103703935710663 , 0.0000456325850473 );
            ir->AddPentPoint( 1031 , 0.0103703935710663 , 0.0786829086131733 , 0.1751123849185410 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1032 , 0.3679171564486097 , 0.0786829086131733 , 0.1751123849185410 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1033 , 0.0786829086131733 , 0.0103703935710663 , 0.1751123849185410 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1034 , 0.3679171564486097 , 0.0103703935710663 , 0.1751123849185410 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1035 , 0.0786829086131733 , 0.3679171564486097 , 0.1751123849185410 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1036 , 0.0103703935710663 , 0.3679171564486097 , 0.1751123849185410 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1037 , 0.0103703935710663 , 0.1751123849185410 , 0.0786829086131733 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1038 , 0.3679171564486097 , 0.1751123849185410 , 0.0786829086131733 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1039 , 0.1751123849185410 , 0.0103703935710663 , 0.0786829086131733 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1040 , 0.3679171564486097 , 0.0103703935710663 , 0.0786829086131733 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1041 , 0.1751123849185410 , 0.3679171564486097 , 0.0786829086131733 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1042 , 0.0103703935710663 , 0.3679171564486097 , 0.0786829086131733 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1043 , 0.0786829086131733 , 0.1751123849185410 , 0.0103703935710663 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1044 , 0.3679171564486097 , 0.1751123849185410 , 0.0103703935710663 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1045 , 0.1751123849185410 , 0.0786829086131733 , 0.0103703935710663 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1046 , 0.3679171564486097 , 0.0786829086131733 , 0.0103703935710663 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1047 , 0.1751123849185410 , 0.3679171564486097 , 0.0103703935710663 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1048 , 0.0786829086131733 , 0.3679171564486097 , 0.0103703935710663 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1049 , 0.0786829086131733 , 0.1751123849185410 , 0.3679171564486097 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1050 , 0.0103703935710663 , 0.1751123849185410 , 0.3679171564486097 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1051 , 0.1751123849185410 , 0.0786829086131733 , 0.3679171564486097 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1052 , 0.0103703935710663 , 0.0786829086131733 , 0.3679171564486097 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1053 , 0.1751123849185410 , 0.0103703935710663 , 0.3679171564486097 , 0.3679171564486097 , 0.0000456325850473 );
            ir->AddPentPoint( 1054 , 0.0786829086131733 , 0.0103703935710663 , 0.3679171564486097 , 0.3679171564486097 , 0.0000456325850473 );

            return ir;

        default:
        {
           int i = (Order / 2) * 2 + 1;   // Get closest odd # >= Order
           AllocIntRule(PentatopeIntRules, i);
           ir = new IntegrationRule;
           ir->GrundmannMollerSimplexRule(i/2,4);
           PentatopeIntRules[i-1] = PentatopeIntRules[i] = ir;
           return ir;
        }
            
    }
    
    return PentatopeIntRules[Order];
}


IntegrationRule& NURBSMeshRules::GetElementRule(const int elem,
                                                const int patch, const int *ijk,
                                                Array<const KnotVector*> const& kv) const
{
   // First check whether a rule has been assigned to element index elem.
   auto search = elementToRule.find(elem);
   if (search != elementToRule.end())
   {
      return *elementRule[search->second];
   }

#ifndef MFEM_THREAD_SAFE
   // If no prescribed rule is given for the current element, a temporary one is
   // formed by restricting a tensor-product of 1D rules to the element. The
   // ownership model for this temporary rule is not thread-safe.

   MFEM_VERIFY(patchRules1D.NumRows(),
               "Undefined rule in NURBSMeshRules::GetElementRule");

   // Use a tensor product of rules on the patch.
   MFEM_VERIFY(kv.Size() == dim, "");

   int np = 1;
   std::vector<std::vector<real_t>> el(dim);

   std::vector<int> npd;
   npd.assign(3, 0);

   for (int d=0; d<dim; ++d)
   {
      const int order = kv[d]->GetOrder();

      const real_t kv0 = (*kv[d])[order + ijk[d]];
      const real_t kv1 = (*kv[d])[order + ijk[d] + 1];

      const bool rightEnd = (order + ijk[d] + 1) == (kv[d]->Size() - 1);

      for (int i=0; i<patchRules1D(patch,d)->Size(); ++i)
      {
         const IntegrationPoint& ip = (*patchRules1D(patch,d))[i];
         if (kv0 <= ip.x && (ip.x < kv1 || rightEnd))
         {
            const real_t x = (ip.x - kv0) / (kv1 - kv0);
            el[d].push_back(x);
            el[d].push_back(ip.weight);
         }
      }

      npd[d] = static_cast<int>(el[d].size() / 2);
      np *= npd[d];
   }

   temporaryElementRule.SetSize(np);

   // Set temporaryElementRule[i + j*npd[0] + k*npd[0]*npd[1]] =
   //     (el[0][2*i], el[1][2*j], el[2][2*k])

   MFEM_VERIFY(npd[0] > 0 && npd[1] > 0, "Assuming 2D or 3D");

   for (int i = 0; i < npd[0]; ++i)
   {
      for (int j = 0; j < npd[1]; ++j)
      {
         for (int k = 0; k < std::max(npd[2], 1); ++k)
         {
            const int id = i + j*npd[0] + k*npd[0]*npd[1];
            temporaryElementRule[id].x = el[0][2*i];
            temporaryElementRule[id].y = el[1][2*j];

            temporaryElementRule[id].weight = el[0][(2*i)+1];
            temporaryElementRule[id].weight *= el[1][(2*j)+1];

            if (npd[2] > 0)
            {
               temporaryElementRule[id].z = el[2][2*k];
               temporaryElementRule[id].weight *= el[2][(2*k)+1];
            }
         }
      }
   }

   return temporaryElementRule;
#else
   MFEM_ABORT("Temporary integration rules on NURBS elements "
              "are not thread-safe.");
#endif
}

void NURBSMeshRules::GetIntegrationPointFrom1D(const int patch, int i, int j,
                                               int k, IntegrationPoint & ip)
{
   MFEM_VERIFY(patchRules1D.NumRows() > 0,
               "Assuming patchRules1D is set.");

   ip.weight = (*patchRules1D(patch,0))[i].weight;
   ip.x = (*patchRules1D(patch,0))[i].x;

   if (dim > 1)
   {
      ip.weight *= (*patchRules1D(patch,1))[j].weight;
      ip.y = (*patchRules1D(patch,1))[j].x;  // 1D rule only has x
   }

   if (dim > 2)
   {
      ip.weight *= (*patchRules1D(patch,2))[k].weight;
      ip.z = (*patchRules1D(patch,2))[k].x;  // 1D rule only has x
   }
}

void NURBSMeshRules::Finalize(Mesh const& mesh)
{
   if ((int) pointToElem.size() == npatches) { return; }  // Already set

   MFEM_VERIFY(elementToRule.empty() && patchRules1D.NumRows() > 0
               && npatches > 0, "Assuming patchRules1D is set.");
   MFEM_VERIFY(mesh.NURBSext, "");
   MFEM_VERIFY(mesh.Dimension() == dim, "");

   pointToElem.resize(npatches);
   patchRules1D_KnotSpan.resize(npatches);

   // First, find all the elements in each patch.
   std::vector<std::vector<int>> patchElements(npatches);

   for (int e=0; e<mesh.GetNE(); ++e)
   {
      patchElements[mesh.NURBSext->GetElementPatch(e)].push_back(e);
   }

   Array<int> ijk(3);
   Array<int> maxijk(3);
   Array<int> np(3);  // Number of points in each dimension
   ijk = 0;

   Array<const KnotVector*> pkv;

   for (int p=0; p<npatches; ++p)
   {
      patchRules1D_KnotSpan[p].resize(dim);

      // For each patch, get the range of ijk.
      mesh.NURBSext->GetPatchKnotVectors(p, pkv);
      MFEM_VERIFY((int) pkv.Size() == dim, "");

      maxijk = 1;
      np = 1;
      for (int d=0; d<dim; ++d)
      {
         maxijk[d] = pkv[d]->GetNKS();
         np[d] = patchRules1D(p,d)->Size();
      }

      // For each patch, set a map from ijk to element index.
      Array3D<int> ijk2elem(maxijk[0], maxijk[1], maxijk[2]);
      ijk2elem = -1;

      for (auto elem : patchElements[p])
      {
         mesh.NURBSext->GetElementIJK(elem, ijk);
         MFEM_VERIFY(ijk2elem(ijk[0], ijk[1], ijk[2]) == -1, "");
         ijk2elem(ijk[0], ijk[1], ijk[2]) = elem;
      }

      // For each point, find its ijk and from that its element index.
      // It is assumed here that the NURBSFiniteElement kv the same as the
      // patch kv.

      for (int d=0; d<dim; ++d)
      {
         patchRules1D_KnotSpan[p][d].SetSize(patchRules1D(p,d)->Size());

         for (int r=0; r<patchRules1D(p,d)->Size(); ++r)
         {
            const IntegrationPoint& ip = (*patchRules1D(p,d))[r];

            const int order = pkv[d]->GetOrder();

            // Find ijk_d such that ip.x is in the corresponding knot-span.
            int ijk_d = 0;
            bool found = false;
            while (!found)
            {
               const real_t kv0 = (*pkv[d])[order + ijk_d];
               const real_t kv1 = (*pkv[d])[order + ijk_d + 1];

               const bool rightEnd = (order + ijk_d + 1) == (pkv[d]->Size() - 1);

               if (kv0 <= ip.x && (ip.x < kv1 || rightEnd))
               {
                  found = true;
               }
               else
               {
                  ijk_d++;
               }
            }

            patchRules1D_KnotSpan[p][d][r] = ijk_d;
         }
      }

      pointToElem[p].SetSize(np[0], np[1], np[2]);
      for (int i=0; i<np[0]; ++i)
         for (int j=0; j<np[1]; ++j)
            for (int k=0; k<np[2]; ++k)
            {
               const int elem = ijk2elem(patchRules1D_KnotSpan[p][0][i],
                                         patchRules1D_KnotSpan[p][1][j],
                                         patchRules1D_KnotSpan[p][2][k]);
               MFEM_VERIFY(elem >= 0, "");
               pointToElem[p](i,j,k) = elem;
            }
   } // Loop (p) over patches
}

void NURBSMeshRules::SetPatchRules1D(const int patch,
                                     std::vector<const IntegrationRule*> & ir1D)
{
   MFEM_VERIFY((int) ir1D.size() == dim, "Wrong dimension");

   for (int i=0; i<dim; ++i)
   {
      patchRules1D(patch,i) = ir1D[i];
   }
}

NURBSMeshRules::~NURBSMeshRules()
{
   for (int i=0; i<patchRules1D.NumRows(); ++i)
      for (int j=0; j<patchRules1D.NumCols(); ++j)
      {
         delete patchRules1D(i, j);
      }
}

}
