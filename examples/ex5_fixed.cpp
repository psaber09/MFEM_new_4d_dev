//                                MFEM Example 5
//
// Compile with: make ex5
//
// Sample runs:  ex5 -m ../data/square-disc.mesh
//               ex5 -m ../data/star.mesh
//               ex5 -m ../data/star.mesh -pa
//               ex5 -m ../data/beam-tet.mesh
//               ex5 -m ../data/beam-hex.mesh
//               ex5 -m ../data/beam-hex.mesh -pa
//               ex5 -m ../data/escher.mesh
//               ex5 -m ../data/fichera.mesh
//
// Device sample runs:
//               ex5 -m ../data/star.mesh -pa -d cuda
//               ex5 -m ../data/star.mesh -pa -d raja-cuda
//               ex5 -m ../data/star.mesh -pa -d raja-omp
//               ex5 -m ../data/beam-hex.mesh -pa -d cuda
//
// Description:  This example code solves a simple 2D/3D mixed Darcy problem
//               corresponding to the saddle point system
//
//                                 k*u + grad p = f
//                                 - div u      = g
//
//               with natural boundary condition -p = <given pressure>.
//               Here, we use a given exact solution (u,p) and compute the
//               corresponding r.h.s. (f,g).  We discretize with Raviart-Thomas
//               finite elements (velocity u) and piecewise discontinuous
//               polynomials (pressure p).
//
//               The example demonstrates the use of the BlockOperator class, as
//               well as the collective saving of several grid functions in
//               VisIt (visit.llnl.gov) and ParaView (paraview.org) formats.
//
//               We recommend viewing examples 1-4 before viewing this example.

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace mfem;

// Define the analytical solution and forcing terms / boundary conditions
void uFun_ex(const Vector & x, Vector & u);
real_t pFun_ex(const Vector & x);
void fFun(const Vector & x, Vector & f);
real_t gFun(const Vector & x);
real_t f_natural(const Vector & x);
int dim;
//int ref_levels;

int main(int argc, char *argv[])
{
    // Parse command-line options
    const char *mesh_file = "../data/star.mesh";
    int order = 1;

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&order, "-o", "--order", "Finite element order.");
    //args.AddOption(&ref_levels, "-ref", "--ref-levels", "");

    args.Parse();
    if (!args.Good()) { args.PrintUsage(std::cout); return 1; }
    args.PrintOptions(std::cout);

    std::vector<int> ref_levels = {0, 1, 2, 3, 4};
    std::vector<double> u_errors, p_errors;

    for (int r : ref_levels)
    {
        Mesh mesh(mesh_file, 1, 1);
        dim = mesh.Dimension();
        for (int l = 0; l < r; l++)
        {
            mesh.UniformRefinement();
        }
        
        RT_FECollection hdiv_coll(order, dim);
        L2_FECollection l2_coll(order, dim);
        FiniteElementSpace R_space(&mesh, &hdiv_coll);
        FiniteElementSpace W_space(&mesh, &l2_coll);
        
        Array<int> block_offsets(3);
        block_offsets[0] = 0;
        block_offsets[1] = R_space.GetVSize();
        block_offsets[2] = W_space.GetVSize();
        block_offsets.PartialSum();
                
        BlockVector x(block_offsets), rhs(block_offsets);
        
        // Coefficients
        ConstantCoefficient k(1.0);
        VectorFunctionCoefficient fcoeff(dim, fFun);
        FunctionCoefficient fnatcoeff(f_natural);
        FunctionCoefficient gcoeff(gFun);
        VectorFunctionCoefficient ucoeff(dim, uFun_ex);
        FunctionCoefficient pcoeff(pFun_ex);
        
        // Linear forms with proper memory binding
        LinearForm fform(&R_space);
        fform.Update(&R_space, rhs.GetBlock(0), 0);
        fform.AddDomainIntegrator(new VectorFEDomainLFIntegrator(fcoeff));
        fform.AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator(fnatcoeff));
        fform.Assemble();
        fform.SyncAliasMemory(rhs.GetBlock(0));
        
        LinearForm gform(&W_space);
        gform.Update(&W_space, rhs.GetBlock(1), 0);
        gform.AddDomainIntegrator(new DomainLFIntegrator(gcoeff));
        gform.Assemble();
        gform.SyncAliasMemory(rhs.GetBlock(1));
        
        // Bilinear forms
        BilinearForm mVarf(&R_space);
        mVarf.AddDomainIntegrator(new VectorFEMassIntegrator(k));
        mVarf.Assemble();
        mVarf.Finalize();
        
        MixedBilinearForm bVarf(&R_space, &W_space);
        bVarf.AddDomainIntegrator(new VectorFEDivergenceIntegrator);
        bVarf.Assemble();
        bVarf.Finalize();
        
        // Block operator
        SparseMatrix &M = mVarf.SpMat();
        SparseMatrix &B = bVarf.SpMat();
        B *= -1.0;
        TransposeOperator Bt(&B);
        
        BlockOperator darcyOp(block_offsets);
        darcyOp.SetBlock(0, 0, &M);
        darcyOp.SetBlock(0, 1, &Bt);
        darcyOp.SetBlock(1, 0, &B);
        
        // Preconditioner
        Vector Md;
        M.GetDiag(Md);
        SparseMatrix *MinvBt = Transpose(B);
        for (int i = 0; i < Md.Size(); i++) {
            (*MinvBt).ScaleRow(i, 1.0 / Md(i));
        }
        SparseMatrix *S = Mult(B, *MinvBt);
        
        DSmoother invM(M);
#ifndef MFEM_USE_SUITESPARSE
        GSSmoother invS(*S);
#else
        UMFPackSolver invS(*S);
#endif
        invM.iterative_mode = false;
        invS.iterative_mode = false;
        
        BlockDiagonalPreconditioner darcyPrec(block_offsets);
        darcyPrec.SetDiagonalBlock(0, &invM);
        darcyPrec.SetDiagonalBlock(1, &invS);
        
        // Solve with MINRES
        MINRESSolver solver;
        solver.SetAbsTol(1e-14);
        solver.SetRelTol(1e-12);
        solver.SetMaxIter(2000);
        solver.SetPrintLevel(0);
        solver.SetOperator(darcyOp);
        solver.SetPreconditioner(darcyPrec);
        x = 0.0;
        solver.Mult(rhs, x);
        
        // Error computation
        GridFunction u, p;
        u.MakeRef(&R_space, x.GetBlock(0), 0);
        p.MakeRef(&W_space, x.GetBlock(1), 0);
        
//        int order_quad = std::max(2, 2*order + 1);
//        const IntegrationRule *irs[Geometry::NumGeom];
//        for (int i = 0; i < Geometry::NumGeom; i++) {
//            irs[i] = &IntRules.Get(i, order_quad);
//        }
        
        const IntegrationRule* irs[Geometry::NumGeom];
        for (int i = 0; i < Geometry::NumGeom; i++)
        {
            if (i == 4)
            {
                // Tet Int Rule
                irs[i] = &(IntRules.Get(i, 10));
            }else
            {
                // Everything else
                irs[i] = &(IntRules.Get(i, 16));
            }
        }
        
        double err_u = u.ComputeL2Error(ucoeff, irs);
        double norm_u = ComputeLpNorm(2.0, ucoeff, mesh, irs);
        double err_p = p.ComputeL2Error(pcoeff, irs);
        double norm_p = ComputeLpNorm(2.0, pcoeff, mesh, irs);
        
        u_errors.push_back(err_u / norm_u);
        p_errors.push_back(err_p / norm_p);
        
        std::cout << "RT dofs = " << block_offsets[1] - block_offsets[0] << "\n";
        std::cout << "L2 dofs = " << block_offsets[2] - block_offsets[1] << "\n";
        std::cout << "Total dofs = " << block_offsets.Last() << std::endl;
        
        //std::cout << "||u - u_ex|| = " << err_u << std::endl;
        //std::cout << "||p - p_ex|| =" << err_p << std::endl;
            cout << "ref = " << r
                 << ", DOFs = " << block_offsets.Last()
                 << ", ||u - u_ex||/||u_ex|| = " << err_u / norm_u
                 << ", ||p - p_ex||/||p_ex|| = " << err_p / norm_p << endl;
        
        delete MinvBt;
        delete S;
    }


    // Convergence rates
    cout << "\nConvergence Rates:\n";
    for (size_t i = 1; i < u_errors.size(); i++) {
        double rate_u = log(u_errors[i-1] / u_errors[i]) / log(2.0);
        double rate_p = log(p_errors[i-1] / p_errors[i]) / log(2.0);
        cout << "ref " << i << ": rate_u = " << rate_u << ", rate_p = " << rate_p << endl;
    }

    return 0;
}


// Exact pressure: p(x, y) = exp(x) * sin(y)
real_t pFun_ex(const Vector & x)
{
    if (dim == 2)
    {
        return exp(x(0)) * sin(x(1));
    }
    if (dim == 3)
    {
        return exp(x(0)) * sin(x(1)) * cos(x(2));
    }
}

// Exact velocity: u = -grad(p)
void uFun_ex(const Vector & x, Vector & u)
{
    if (dim == 2)
    {
        u.SetSize(2);
        double xi = x(0);
        double yi = x(1);
        double ex = exp(xi);
        u(0) = -ex * sin(yi);  // -dp/dx
        u(1) = -ex * cos(yi);  // -dp/dy
    }
    if (dim == 3)
    {
        u.SetSize(3);
        u(0) = - exp(x(0)) * sin(x(1)) * cos(x(2)); //-dp/dx
        u(1) = - exp(x(0)) * cos(x(1)) * cos(x(2)); // -dp/dy
        u(2) = exp(x(0)) * sin(x(1)) * sin(x(2)); // -dp/dz
    }
}


// Right-hand side: f = k*u + grad(p) = 0 for this exact solution
void fFun(const Vector & x, Vector & f)
{
    if (dim == 2)
    {
        f.SetSize(2);
        double xi = x(0);
        double yi = x(1);
        double ex = exp(xi);
        double siny = sin(yi);
        double cosy = cos(yi);
        
        double dpdx = ex * siny;
        double dpdy = ex * cosy;
        
        // f = k*u + grad(p) = -grad(p) + grad(p) = 0
        f(0) = -dpdx + dpdx;
        f(1) = -dpdy + dpdy;
    }
    if (dim == 3)
    {
        f.SetSize(3);
        double xi = x(0);
        double yi = x(1);
        double zi = x(2);
        double ex = exp(xi);
        double siny = sin(yi);
        double cosy = cos(yi);
        double cosz = cos(zi);
        double nsinz = -1.0 * sin(zi);
        
        double dpdx = ex * siny * cosz;
        double dpdy = ex * cosy * cosz;
        double dpdz = ex * siny * nsinz;
        
        // f = k*u + grad(p) = -grad(p) + grad(p) = 0
        f(0) = -dpdx + dpdx;
        f(1) = -dpdy + dpdy;
        f(2) = -dpdz + dpdz;

    }
}

// Divergence condition: g = -div(u) = -div(-grad(p)) = Δp = 0
real_t gFun(const Vector & x)
{
    if (dim == 2)
    {
        // Divergence condition: g = -div(u) = -div(-grad(p)) = Δp = 0
        return 0.0;
    }
    if (dim == 3)
    {
        // Divergence condition: g = -div(u) = -div(-grad(p)) = Δp = -p(x,y,z)
        double g = -1.0*(exp(x(0)) * sin(x(1)) * cos(x(2)));
        return g;
    }
}

// Natural boundary term (Neumann-like condition): -p = -p(x)
real_t f_natural(const Vector & x)
{
    return -pFun_ex(x);
}






