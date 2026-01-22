//
//  fe_skgrad.cpp
//  mfem
//
//  Created by psaber on 6/20/25.
//

#include "fe_curl.hpp"
#include "face_map_utils.hpp"
#include "../coefficient.hpp"
#include <fstream>


namespace mfem
{

using namespace std;

//const double HCurl_PentatopeElement::tk[40] =
//{ 0,0,-0.5,-0.5,  0,-0.5,0,0.5,  0,-0.5,-0.5,0,  -0.5,0,0,-0.5,  -0.5,0,-0.5,0,  -0.5,-0.5,0,0,  0.5,0.5,0.5,0,  0.5,0.5,0,0.5,  0.5,0,0.5,0.5,  0,0.5,0.5,0.5 };
//
//const double HCurl_PentatopeElement::tk1[40] =
//{1,0,0,0, 1,0,0,0, 1,0,0,0, 0,1,0,0, 0,1,0,0, 0,0,1,0, -1,1,0,0, -1,1,0,0, -1,0,1,0, 0,-1,1,0};
//
//const double HCurl_PentatopeElement::tk2[40] =
//{-1,1,0,0, -1,0,1,0, -1,0,0,1, 0,-1,1,0, 0,-1,0,1, 0,0,-1,1, 0,-1,1,0, 0,-1,0,1, 0,0,-1,1, 0,0,-1,1};

//const double HSkwGrad_PentatopeElement::tk[40] =
//{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1, -1,1,0,0, -1,0,1,0, -1,0,0,1, 0,-1,1,0, 0,-1,0,1, 0,0,-1,1};

//const double HCurl_PentatopeElement::c = 1./5.;

HCurl_PentatopeElement::HCurl_PentatopeElement(const int p)
: VectorFiniteElement(4, Geometry::PENTATOPE, p*(p + 2)*(p + 3)*(p + 4)/6,
                      p, H_CURL, FunctionSpace::Pk),
dof2tk(dof), doftrans(p)
{
    
//    const real_t *eop = poly1d.OpenPoints(p - 1);
//    const real_t *fop = (p > 1) ? poly1d.OpenPoints(p - 2) : NULL;
//    const real_t *ftop = (p > 2) ? poly1d.OpenPoints(p - 3) : NULL;
//    const real_t *iop = (p > 3) ? poly1d.OpenPoints(p - 4) : NULL;
    
    const real_t *fop = poly1d.OpenPoints(p - 1);
    const real_t *ftop = (p > 1) ? poly1d.OpenPoints(p - 2) : NULL;
    const real_t *iop = (p > 2) ? poly1d.OpenPoints(p - 3) : NULL;
    
    const int pm1 = p - 1, pm2 = p - 2, pm3 = p - 3, pm4 = p -4;
    
    // Lamda function for computing the scaled skew-symmetric-outer-product of two four vectors
    auto skw_sym_outerprod_fnc = [&](const std::vector<double>& vec1, const std::vector<double>& vec2)
    {

        DenseMatrix result(4,4);
        double s_factor = 1.0;
        //std::cout << "Skw-sym Prod" << std::endl;
        // "Forward" Product
        for (int i = 0; i < vec1.size(); ++i)
        {
            for (int j = 0; j < vec2.size(); ++j)
            {
                // Forward Product
                double outer_prod_1 = vec1[i] * vec2[j];
                
                // "Backward" Product
                double outer_prod_2 = vec2[i] * vec1[j];
                
                // Compute Result
                result(i,j) = s_factor*(outer_prod_1 - outer_prod_2);
                //std::cout << 0.5*(outer_prod_1 - outer_prod_2) << std::endl;
            }
            //std::cout << "New COl --------" << std::endl;
        }
        
        return result;
    }; // end of lamda function

    
#ifndef MFEM_THREAD_SAFE
    shape_x.SetSize(p);
    shape_y.SetSize(p);
    shape_z.SetSize(p);
    shape_t.SetSize(p);
    shape_l.SetSize(p);
    dshape_x.SetSize(p);
    dshape_y.SetSize(p);
    dshape_z.SetSize(p);
    dshape_t.SetSize(p);
    dshape_l.SetSize(p);
    u.SetSize(dof, dim);
#else
    Vector shape_x(p), shape_y(p), shape_z(p), shape_t(p),
    shape_l(p);
#endif
    // edges (see Tetrahedron::edges in mesh/tetrahedron.cpp)
    int o = 0;
    
    // faces (see Mesh::GeneratePlanars in mesh/mesh.cpp)
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (0,1,2)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
           double xcoord = fop[i]/w;
           double ycoord = fop[j]/w;
          Nodes.IntPoint(o).Set4(fop[i]/w, fop[j]/w, 0.0, 0.0);
          dof2tk[o++] = 0;
          //std::cout << "inside dof" << std::endl;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (0,1,3)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(fop[i]/w, 0.0, fop[j]/w, 0.0);
          dof2tk[o++] = 1;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (0,1,4)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(fop[i]/w, 0.0, 0.0, fop[j]/w);
          dof2tk[o++] = 2;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (0,2,3)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(0.0, fop[i]/w, fop[j]/w, 0.0);
          dof2tk[o++] = 3;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (0,2,4)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(0.0, fop[i]/w, 0.0, fop[j]/w);
          dof2tk[o++] = 4;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (0,3,4)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(0.0, 0.0, fop[i]/w, fop[j]/w);
          dof2tk[o++] = 5;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (1,2,3)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(fop[pm1-i-j]/w, fop[i]/w, fop[j]/w, 0.0);
          dof2tk[o++] = 6;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (1,2,4)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(fop[pm1-i-j]/w, fop[i]/w, 0.0, fop[j]/w);
          dof2tk[o++] = 7;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (1,3,4)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(fop[pm1-i-j]/w, 0.0, fop[i]/w, fop[j]/w);
          dof2tk[o++] = 8;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (2,3,4)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(0.0, fop[pm1-i-j]/w, fop[i]/w, fop[j]/w);
          dof2tk[o++] = 9;
       }

    // facets (see Mesh::GenerateFaces in mesh/mesh.cpp)
    for (int k = 0; k <= pm2; k++)
       for (int j = 0; j + k <= pm2; j++)
          for (int i = 0; i + j + k <= pm2; i++)  // (0,1,2,3)
          {
             double w = ftop[i] + ftop[j] + ftop[k] + ftop[pm2-i-j-k];
             Nodes.IntPoint(o).Set4(ftop[i]/w, ftop[j]/w, ftop[k]/w, 0.0);
             dof2tk[o++] = 0;
             Nodes.IntPoint(o).Set4(ftop[i]/w, ftop[j]/w, ftop[k]/w, 0.0);
             dof2tk[o++] = 1;

          }
    for (int k = 0; k <= pm2; k++)
       for (int j = 0; j + k <= pm2; j++)
          for (int i = 0; i + j + k <= pm2; i++)  // (0,2,1,4)
          {
             double w = ftop[i] + ftop[j] + ftop[k] + ftop[pm2-i-j-k];
             Nodes.IntPoint(o).Set4(ftop[j]/w, ftop[i]/w, 0.0, ftop[k]/w);
             dof2tk[o++] = 1;
             Nodes.IntPoint(o).Set4(ftop[j]/w, ftop[i]/w, 0.0, ftop[k]/w);
             dof2tk[o++] = 0;

          }
    for (int k = 0; k <= pm2; k++)
       for (int j = 0; j + k <= pm2; j++)
          for (int i = 0; i + j + k <= pm2; i++)  // (0,1,3,4)
          {
             double w = ftop[i] + ftop[j] + ftop[k] + ftop[pm2-i-j-k];
             Nodes.IntPoint(o).Set4(ftop[i]/w, 0.0, ftop[j]/w, ftop[k]/w);
             dof2tk[o++] = 0;
             Nodes.IntPoint(o).Set4(ftop[i]/w, 0.0, ftop[j]/w, ftop[k]/w);
             dof2tk[o++] = 2;

          }
    for (int k = 0; k <= pm2; k++)
       for (int j = 0; j + k <= pm2; j++)
          for (int i = 0; i + j + k <= pm2; i++)  // (0,3,2,4)
          {
             double w = ftop[i] + ftop[j] + ftop[k] + ftop[pm2-i-j-k];
             Nodes.IntPoint(o).Set4(0.0, ftop[j]/w, ftop[i]/w, ftop[k]/w);
             dof2tk[o++] = 2;
             Nodes.IntPoint(o).Set4(0.0, ftop[j]/w, ftop[i]/w, ftop[k]/w);
             dof2tk[o++] = 1;

          }
    for (int k = 0; k <= pm2; k++)
       for (int j = 0; j + k <= pm2; j++)
          for (int i = 0; i + j + k <= pm2; i++)  // (1,2,3,4)
          {
             double w = ftop[i] + ftop[j] + ftop[k] + ftop[pm2-i-j-k];
             Nodes.IntPoint(o).Set4(ftop[pm2-i-j-k]/w, ftop[i]/w, ftop[j]/w, ftop[k]/w);
             dof2tk[o++] = 4;
             Nodes.IntPoint(o).Set4(ftop[pm2-i-j-k]/w, ftop[i]/w, ftop[j]/w, ftop[k]/w);
             dof2tk[o++] = 5;

          }

    // interior bubbles
    for (int l = 0; l <= pm3; l++)
       for (int k = 0; k + l <= pm3; k++)
          for (int j = 0; j + k + l <= pm3; j++)
             for (int i = 0; i + j + k + l <= pm3; i++)
             {
                double w = iop[i] + iop[j] + iop[k] + iop[l] + iop[pm3-i-j-k-l];
                Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
                dof2tk[o++] = 0;
                Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
                dof2tk[o++] = 1;
                Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
                dof2tk[o++] = 2;
                Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
                dof2tk[o++] = 3;
                Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
                dof2tk[o++] = 4;
                Nodes.IntPoint(o).Set4(iop[i]/w, iop[j]/w, iop[k]/w, iop[l]/w);
                dof2tk[o++] = 5;

             }
    /*
    DenseMatrix T(dof);
    DenseMatrix B(4, 4);
    //std::vector<DenseMatrix> Mat_B(dof);
    int num_dof = dof;

    double Tensor_B[4][4][num_dof];
    double Mat_B[4][num_dof];
    
    for (int q = 0; q < dof; q++)
    {
        const IntegrationPoint &ip = Nodes.IntPoint(q);
        //std::cout << "Dof " << q << " = " << ip.x << ", " << ip.y << ", " << ip.z << ", " << ip.t << std::endl;
        
        //const double *nm = nk + 4*dof2nk[m];
        
        double tm1[4] = {tk1[4*dof2tk[q]], tk1[4*dof2tk[q]+1], tk1[4*dof2tk[q]+2], tk1[4*dof2tk[q]+3]};
        
        double tm2[4] = {tk2[4*dof2tk[q]], tk2[4*dof2tk[q]+1], tk2[4*dof2tk[q]+2], tk2[4*dof2tk[q]+3]};

        //const Vector tm({tk[4*dof2tk[q]], tk[4*dof2tk[q]+1], tk[4*dof2tk[q]+2], tk[4*dof2tk[q]+3]});
        //std::cout << "tm < " << tm(0) << "," << tm(1) << "," << tm(2) << "," << tm(3) << ">" <<std::endl;
        int o = 0;
        //int num_dof = dof;
        
        //compute barycentric coordinates as function of ip
        std::vector<double> bary_vector{ip.x, ip.y, ip.z, ip.t, (1.0 - ip.x - ip.y - ip.z - ip.t)};
        
        // compute the gradient of the barycentric coords
        std::vector<double> gradL1{1,0,0,0};
        std::vector<double> gradL2{0,1,0,0};
        std::vector<double> gradL3{0,0,1,0};
        std::vector<double> gradL4{0,0,0,1};
        std::vector<double> gradL5{-1,-1,-1,-1};
        
        std::vector<std::vector<double>> gradbary_vector{gradL1, gradL2, gradL3, gradL4, gradL5};
        
        double La, Lb, Lc, Ld, Le;
        
        std::vector<double> grad_La, grad_Lb, grad_Lc, grad_Ld, grad_Le;
        
        
        //Faces
        for(int i=0; i<p;i++)
        {
            for(int j=0; j<p;j++)
            {
                for(int a=0; a<5;a++)
                {
                    for(int b=0; b<5;b++)
                    {
                        for(int c=0; c<5;c++)
                        {
                            if((a<b)&&(b<c)&&((i+j)<p))
                            {
                                
                                La = bary_vector[a];
                                Lb = bary_vector[b];
                                Lc = bary_vector[c];
                                
                                grad_La = gradbary_vector[a];
                                grad_Lb = gradbary_vector[b];
                                grad_Lc = gradbary_vector[c];
                                
                                // compute polynomials
                                std::vector<double> Legendre_i;
                                double x = Lb;
                                double y = La + Lb;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                std::vector<double> Jacobi_j;
                                x = Lc;
                                y = La + Lb + Lc;
                                double alpha = 2*i + 1;
                                poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                                
                                // Scaled Skew-sym outer product
                                DenseMatrix skw_sym_outerprod_1 = skw_sym_outerprod_fnc(grad_Lb, grad_Lc);
                                DenseMatrix skw_sym_outerprod_2 = skw_sym_outerprod_fnc(grad_Lc, grad_La);
                                DenseMatrix skw_sym_outerprod_3 = skw_sym_outerprod_fnc(grad_La, grad_Lb);

                                
                                // Add Basis Funcitons
                                Tensor_B[0][0][o] = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(0,0) + Lb*skw_sym_outerprod_2(0,0) + Lc*skw_sym_outerprod_3(0,0));
                                Tensor_B[0][1][o] = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(0,1) + Lb*skw_sym_outerprod_2(0,1) + Lc*skw_sym_outerprod_3(0,1));
                                Tensor_B[0][2][o] = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(0,2) + Lb*skw_sym_outerprod_2(0,2) + Lc*skw_sym_outerprod_3(0,2));
                                Tensor_B[0][3][o] = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(0,3) + Lb*skw_sym_outerprod_2(0,3) + Lc*skw_sym_outerprod_3(0,3));
                                Tensor_B[1][0][o] = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(1,0) + Lb*skw_sym_outerprod_2(1,0) + Lc*skw_sym_outerprod_3(1,0));
                                Tensor_B[1][1][o] = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(1,1) + Lb*skw_sym_outerprod_2(1,1) + Lc*skw_sym_outerprod_3(1,1));
                                Tensor_B[1][2][o] = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(1,2) + Lb*skw_sym_outerprod_2(1,2) + Lc*skw_sym_outerprod_3(1,2));
                                Tensor_B[1][3][o] = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(1,3) + Lb*skw_sym_outerprod_2(1,3) + Lc*skw_sym_outerprod_3(1,3));
                                Tensor_B[2][0][o] = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(2,0) + Lb*skw_sym_outerprod_2(2,0) + Lc*skw_sym_outerprod_3(2,0));
                                Tensor_B[2][1][o] = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(2,1) + Lb*skw_sym_outerprod_2(2,1) + Lc*skw_sym_outerprod_3(2,1));
                                Tensor_B[2][2][o] = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(2,2) + Lb*skw_sym_outerprod_2(2,2) + Lc*skw_sym_outerprod_3(2,2));
                                Tensor_B[2][3][o] = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(2,3) + Lb*skw_sym_outerprod_2(2,3) + Lc*skw_sym_outerprod_3(2,3));
                                Tensor_B[3][0][o] = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(3,0) + Lb*skw_sym_outerprod_2(3,0) + Lc*skw_sym_outerprod_3(3,0));
                                Tensor_B[3][1][o] = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(3,1) + Lb*skw_sym_outerprod_2(3,1) + Lc*skw_sym_outerprod_3(3,1));
                                Tensor_B[3][2][o] = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(3,2) + Lb*skw_sym_outerprod_2(3,2) + Lc*skw_sym_outerprod_3(3,2));
                                Tensor_B[3][3][o] = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(3,3) + Lb*skw_sym_outerprod_2(3,3) + Lc*skw_sym_outerprod_3(3,3));
                                
                                //Mat_B.push_back(B);
                                //Mat_B[o] = B;
                                o++;
                                
                                
                            }
                        }
                    }
                }
            }
        }
        // end of faces
        

        // Initialize result to zero
        for (int j = 0; j < 4; j++) {
            for (int n = 0; n < num_dof; n++) {
                Mat_B[j][n] = 0.0;
            }
        }

        // Perform contraction over i
        for (int n = 0; n < num_dof; n++) {
            for (int j = 0; j < 4; j++) {
                double sum = 0.0;
                for (int i = 0; i < 4; i++) {
                    sum += Tensor_B[i][j][n] * tm1[i];
                }
                Mat_B[j][n] = sum;
            }
        }
        // Perfrom contraction over j
        std::vector<double> vec_B(num_dof);
        for (int n = 0; n < num_dof; n++)
        {
            std::cout << "New Dof" << std::endl;
            double sum = 0.0;
            for (int j = 0; j < 4; j++)
            {
                sum += tm2[j] * Mat_B[j][n];
            }
            std::cout << "sum value" << std::endl;
            std::cout << sum << std::endl;
            vec_B[n] = sum;  // This vector is of size (1 x N)
        }
        
        // Add to Vandermonde matrix
        for (int n = 0; n < num_dof; n++) 
        {
            T(n,q) = vec_B[n];
        }
        
        //vec_B = T.GetColumn(q);
    
        //B.Mult(tm, T.GetColumn(q));

    }*/
        
//    std::ofstream T_file("T_4DCurl.txt");
//    T.PrintMatlab(T_file);
    
    //Ti.Factor(T);
    
    
//    std::ofstream Ti_file("Ti_4DCurl.txt");
//    Ti.PrintMatlab(Ti_file);

    //mfem::out << "HCurl_PentatopeElement(" << p << ") : "; Ti.TestInversion();
    
}

const double HCurl_PentatopeElement::tk1[10][4] =
{{1,0,0,0},{1,0,0,0},{1,0,0,0},{0,1,0,0},{0,1,0,0},{0,0,1,0},{-1,1,0,0},{-1,1,0,0},{-1,0,1,0},{0,-1,1,0}};

const double HCurl_PentatopeElement::tk2[10][4] =
{{-1,1,0,0},{-1,0,1,0},{-1,0,0,1},{0,-1,1,0},{0,-1,0,1},{0,0,-1,1},{0,-1,1,0},{0,-1,0,1},{0,0,-1,1},{0,0,-1,1}};

void HCurl_PentatopeElement::CalcVShape(const IntegrationPoint &ip,
                                           DenseMatrix &shape) const
{
    
    /*
    const int p = order;
    //std::cout << "p_CalcShape = " << p << std::endl;
    int reset_o;
    
#ifdef MFEM_THREAD_SAFE
    Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p + 1),
    shape_l(p + 1);
    DenseMatrix u(Dof, Dim);
#endif
    int size_ip = Nodes.Size();
    //std::cout << "Inside" << std::endl;
    //std::cout << "Int Point = " << ip.x << ", " << ip.y << ", " << ip.z << ", " << ip.t << std::endl;
    
    int o = 0;
    int num_dof = dof;
    u.SetSize(num_dof,4*dim);
    
    //compute barycentric coordinates as function of ip
    std::vector<double> bary_vector{ip.x, ip.y, ip.z, ip.t, (1.0 - ip.x - ip.y - ip.z - ip.t)};
    
    // compute the gradient of the barycentric coords
    std::vector<double> gradL1{1,0,0,0};
    std::vector<double> gradL2{0,1,0,0};
    std::vector<double> gradL3{0,0,1,0};
    std::vector<double> gradL4{0,0,0,1};
    std::vector<double> gradL5{-1,-1,-1,-1};
    
    std::vector<std::vector<double>> gradbary_vector{gradL1, gradL2, gradL3, gradL4, gradL5};
    
    double La, Lb, Lc, Ld, Le;
    
    std::vector<double> grad_La, grad_Lb, grad_Lc, grad_Ld, grad_Le;
    
    
    // Lamda function for computing the scaled skew-symmetric-outer-product of two four vectors
    auto skw_sym_outerprod_fnc = [&](const std::vector<double>& vec1, const std::vector<double>& vec2)
    {

        DenseMatrix result(4,4);
        double s_factor = 1.0;
        //std::cout << "Skw-sym Prod" << std::endl;
        // "Forward" Product
        for (int i = 0; i < vec1.size(); ++i)
        {
            for (int j = 0; j < vec2.size(); ++j)
            {
                // Forward Product
                double outer_prod_1 = vec1[i] * vec2[j];
                
                // "Backward" Product
                double outer_prod_2 = vec2[i] * vec1[j];
                
                // Compute Result
                result(i,j) = s_factor*(outer_prod_1 - outer_prod_2);
                //std::cout << 0.5*(outer_prod_1 - outer_prod_2) << std::endl;
            }
            //std::cout << "New COl --------" << std::endl;
        }
        
        return result;
    }; // end of lamda function

    
    
    //Faces
    for(int i=0; i<p;i++)
    {
        for(int j=0; j<p;j++)
        {
            for(int a=0; a<5;a++)
            {
                for(int b=0; b<5;b++)
                {
                    for(int c=0; c<5;c++)
                    {
                        if((a<b)&&(b<c)&&((i+j)<p))
                        {
                            
                            La = bary_vector[a];
                            Lb = bary_vector[b];
                            Lc = bary_vector[c];
                            
                            grad_La = gradbary_vector[a];
                            grad_Lb = gradbary_vector[b];
                            grad_Lc = gradbary_vector[c];
                            
                            // compute polynomials
                            std::vector<double> Legendre_i;
                            double x = Lb;
                            double y = La + Lb;
                            poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                            
                            std::vector<double> Jacobi_j;
                            x = Lc;
                            y = La + Lb + Lc;
                            double alpha = 2*i + 1;
                            poly1d.CalcJacobi(j, x, y, alpha, Jacobi_j);
                            
                            // Scaled Skew-sym outer product
                            DenseMatrix skw_sym_outerprod_1 = skw_sym_outerprod_fnc(grad_Lb, grad_Lc);
                            DenseMatrix skw_sym_outerprod_2 = skw_sym_outerprod_fnc(grad_Lc, grad_La);
                            DenseMatrix skw_sym_outerprod_3 = skw_sym_outerprod_fnc(grad_La, grad_Lb);

                            
                            // Add Basis Funcitons
//                            u(o, 0) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(0,0) + Lb*skw_sym_outerprod_2(0,0) + Lc*skw_sym_outerprod_3(0,0));
//                            u(o, 1) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(0,1) + Lb*skw_sym_outerprod_2(0,1) + Lc*skw_sym_outerprod_3(0,1));
//                            u(o, 2) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(0,2) + Lb*skw_sym_outerprod_2(0,2) + Lc*skw_sym_outerprod_3(0,2));
//                            u(o, 3) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(0,3) + Lb*skw_sym_outerprod_2(0,3) + Lc*skw_sym_outerprod_3(0,3));
//                            u(o, 4) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(1,0) + Lb*skw_sym_outerprod_2(1,0) + Lc*skw_sym_outerprod_3(1,0));
//                            u(o, 5) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(1,1) + Lb*skw_sym_outerprod_2(1,1) + Lc*skw_sym_outerprod_3(1,1));
//                            u(o, 6) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(1,2) + Lb*skw_sym_outerprod_2(1,2) + Lc*skw_sym_outerprod_3(1,2));
//                            u(o, 7) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(1,3) + Lb*skw_sym_outerprod_2(1,3) + Lc*skw_sym_outerprod_3(1,3));
//                            u(o, 8) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(2,0) + Lb*skw_sym_outerprod_2(2,0) + Lc*skw_sym_outerprod_3(2,0));
//                            u(o, 9) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(2,1) + Lb*skw_sym_outerprod_2(2,1) + Lc*skw_sym_outerprod_3(2,1));
//                            u(o, 10) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(2,2) + Lb*skw_sym_outerprod_2(2,2) + Lc*skw_sym_outerprod_3(2,2));
//                            u(o, 11) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(2,3) + Lb*skw_sym_outerprod_2(2,3) + Lc*skw_sym_outerprod_3(2,3));
//                            u(o, 12) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(3,0) + Lb*skw_sym_outerprod_2(3,0) + Lc*skw_sym_outerprod_3(3,0));
//                            u(o, 13) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(3,1) + Lb*skw_sym_outerprod_2(3,1) + Lc*skw_sym_outerprod_3(3,1));
//                            u(o, 14) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(3,2) + Lb*skw_sym_outerprod_2(3,2) + Lc*skw_sym_outerprod_3(3,2));
//                            u(o, 15) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(3,3) + Lb*skw_sym_outerprod_2(3,3) + Lc*skw_sym_outerprod_3(3,3));

                            shape(o, 0) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(0,0) + Lb*skw_sym_outerprod_2(0,0) + Lc*skw_sym_outerprod_3(0,0));
                            shape(o, 1) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(0,1) + Lb*skw_sym_outerprod_2(0,1) + Lc*skw_sym_outerprod_3(0,1));
                            shape(o, 2) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(0,2) + Lb*skw_sym_outerprod_2(0,2) + Lc*skw_sym_outerprod_3(0,2));
                            shape(o, 3) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(0,3) + Lb*skw_sym_outerprod_2(0,3) + Lc*skw_sym_outerprod_3(0,3));
                            shape(o, 4) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(1,0) + Lb*skw_sym_outerprod_2(1,0) + Lc*skw_sym_outerprod_3(1,0));
                            shape(o, 5) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(1,1) + Lb*skw_sym_outerprod_2(1,1) + Lc*skw_sym_outerprod_3(1,1));
                            shape(o, 6) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(1,2) + Lb*skw_sym_outerprod_2(1,2) + Lc*skw_sym_outerprod_3(1,2));
                            shape(o, 7) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(1,3) + Lb*skw_sym_outerprod_2(1,3) + Lc*skw_sym_outerprod_3(1,3));
                            shape(o, 8) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(2,0) + Lb*skw_sym_outerprod_2(2,0) + Lc*skw_sym_outerprod_3(2,0));
                            shape(o, 9) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(2,1) + Lb*skw_sym_outerprod_2(2,1) + Lc*skw_sym_outerprod_3(2,1));
                            shape(o, 10) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(2,2) + Lb*skw_sym_outerprod_2(2,2) + Lc*skw_sym_outerprod_3(2,2));
                            shape(o, 11) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(2,3) + Lb*skw_sym_outerprod_2(2,3) + Lc*skw_sym_outerprod_3(2,3));
                            shape(o, 12) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(3,0) + Lb*skw_sym_outerprod_2(3,0) + Lc*skw_sym_outerprod_3(3,0));
                            shape(o, 13) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(3,1) + Lb*skw_sym_outerprod_2(3,1) + Lc*skw_sym_outerprod_3(3,1));
                            shape(o, 14) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(3,2) + Lb*skw_sym_outerprod_2(3,2) + Lc*skw_sym_outerprod_3(3,2));
                            shape(o, 15) = Legendre_i[Legendre_i.size()-1] * Jacobi_j[Jacobi_j.size()-1] * (La*skw_sym_outerprod_1(3,3) + Lb*skw_sym_outerprod_2(3,3) + Lc*skw_sym_outerprod_3(3,3));
                            o++;
                            
                            
                        }
                    }
                }
            }
        }
    }
    // end of faces
    
    for (int mod; mod < 10; mod++)
    {
        std::cout << "Start -----" << std::endl;
        std::cout << shape(mod,0) << std::endl;
        std::cout << shape(mod,1) << std::endl;
        std::cout << shape(mod,2) << std::endl;
        std::cout << shape(mod,3) << std::endl;
        std::cout << shape(mod,4) << std::endl;
        std::cout << shape(mod,5) << std::endl;
        std::cout << shape(mod,6) << std::endl;
        std::cout << shape(mod,7) << std::endl;
        std::cout << shape(mod,8) << std::endl;
        std::cout << shape(mod,9) << std::endl;
        std::cout << shape(mod,10) << std::endl;
        std::cout << shape(mod,11) << std::endl;
        std::cout << shape(mod,12) << std::endl;
        std::cout << shape(mod,13) << std::endl;
        std::cout << shape(mod,14) << std::endl;
        std::cout << shape(mod,15) << std::endl;


    }
    
    std::cout << "The End " << std::endl;*/
    
    //Ti.Mult(u, shape);

    double x1 = ip.x, x2 = ip.y, x3 = ip.z, x4 = ip.t;

    shape(0,0) = 0.;
    shape(0,1) = 0.;
    shape(0,2) = x1;
    shape(0,3) = -1.*x1;
    shape(0,4) = 0.;
    shape(0,5) = 0.;
    shape(0,6) = x2;
    shape(0,7) = -1.*x2;
    shape(0,8) = -1.*x1;
    shape(0,9) = -1.*x2;
    shape(0,10) = 0.;
    shape(0,11) = 1. - 1.*x3 - 1.*x4;
    shape(0,12) = x1;
    shape(0,13) = x2;
    shape(0,14) = -1. + x3 + x4;
    shape(0,15) = 0.;

    shape(1,0) = 0.;
    shape(1,1) = -1.*x1;
    shape(1,2) = 0.;
    shape(1,3) = x1;
    shape(1,4) = x1;
    shape(1,5) = 0.;
    shape(1,6) = x3;
    shape(1,7) = -1. + x2 + x4;
    shape(1,8) = 0.;
    shape(1,9) = -1.*x3;
    shape(1,10) = 0.;
    shape(1,11) = x3;
    shape(1,12) = -1.*x1;
    shape(1,13) = 1. - 1.*x2 - 1.*x4;
    shape(1,14) = -1.*x3;
    shape(1,15) = 0.;

    shape(2,0) = 0.;
    shape(2,1) = x1;
    shape(2,2) = -1.*x1;
    shape(2,3) = 0.;
    shape(2,4) = -1.*x1;
    shape(2,5) = 0.;
    shape(2,6) = 1. - 1.*x2 - 1.*x3;
    shape(2,7) = -1.*x4;
    shape(2,8) = x1;
    shape(2,9) = -1. + x2 + x3;
    shape(2,10) = 0.;
    shape(2,11) = x4;
    shape(2,12) = 0.;
    shape(2,13) = x4;
    shape(2,14) = -1.*x4;
    shape(2,15) = 0.;

    shape(3,0) = 0.;
    shape(3,1) = -1.*x2;
    shape(3,2) = -1.*x3;
    shape(3,3) = 1. - 1.*x1 - 1.*x4;
    shape(3,4) = x2;
    shape(3,5) = 0.;
    shape(3,6) = 0.;
    shape(3,7) = -1.*x2;
    shape(3,8) = x3;
    shape(3,9) = 0.;
    shape(3,10) = 0.;
    shape(3,11) = -1.*x3;
    shape(3,12) = -1. + x1 + x4;
    shape(3,13) = x2;
    shape(3,14) = x3;
    shape(3,15) = 0.;

    shape(4,0) = 0.;
    shape(4,1) = x2;
    shape(4,2) = -1. + x1 + x3;
    shape(4,3) = x4;
    shape(4,4) = -1.*x2;
    shape(4,5) = 0.;
    shape(4,6) = x2;
    shape(4,7) = 0.;
    shape(4,8) = 1. - 1.*x1 - 1.*x3;
    shape(4,9) = -1.*x2;
    shape(4,10) = 0.;
    shape(4,11) = -1.*x4;
    shape(4,12) = -1.*x4;
    shape(4,13) = 0.;
    shape(4,14) = x4;
    shape(4,15) = 0.;

    shape(5,0) = 0.;
    shape(5,1) = 1. - 1.*x1 - 1.*x2;
    shape(5,2) = -1.*x3;
    shape(5,3) = -1.*x4;
    shape(5,4) = -1. + x1 + x2;
    shape(5,5) = 0.;
    shape(5,6) = x3;
    shape(5,7) = x4;
    shape(5,8) = x3;
    shape(5,9) = -1.*x3;
    shape(5,10) = 0.;
    shape(5,11) = 0.;
    shape(5,12) = x4;
    shape(5,13) = -1.*x4;
    shape(5,14) = 0.;
    shape(5,15) = 0.;

    shape(6,0) = 0.;
    shape(6,1) = 0.;
    shape(6,2) = 0.;
    shape(6,3) = x1;
    shape(6,4) = 0.;
    shape(6,5) = 0.;
    shape(6,6) = 0.;
    shape(6,7) = x2;
    shape(6,8) = 0.;
    shape(6,9) = 0.;
    shape(6,10) = 0.;
    shape(6,11) = x3;
    shape(6,12) = -1.*x1;
    shape(6,13) = -1.*x2;
    shape(6,14) = -1.*x3;
    shape(6,15) = 0.;

    shape(7,0) = 0.;
    shape(7,1) = 0.;
    shape(7,2) = -1.*x1;
    shape(7,3) = 0.;
    shape(7,4) = 0.;
    shape(7,5) = 0.;
    shape(7,6) = -1.*x2;
    shape(7,7) = 0.;
    shape(7,8) = x1;
    shape(7,9) = x2;
    shape(7,10) = 0.;
    shape(7,11) = x4;
    shape(7,12) = 0.;
    shape(7,13) = 0.;
    shape(7,14) = -1.*x4;
    shape(7,15) = 0.;

    shape(8,0) = 0.;
    shape(8,1) = x1;
    shape(8,2) = 0.;
    shape(8,3) = 0.;
    shape(8,4) = -1.*x1;
    shape(8,5) = 0.;
    shape(8,6) = -1.*x3;
    shape(8,7) = -1.*x4;
    shape(8,8) = 0.;
    shape(8,9) = x3;
    shape(8,10) = 0.;
    shape(8,11) = 0.;
    shape(8,12) = 0.;
    shape(8,13) = x4;
    shape(8,14) = 0.;
    shape(8,15) = 0.;

    shape(9,0) = 0.;
    shape(9,1) = x2;
    shape(9,2) = x3;
    shape(9,3) = x4;
    shape(9,4) = -1.*x2;
    shape(9,5) = 0.;
    shape(9,6) = 0.;
    shape(9,7) = 0.;
    shape(9,8) = -1.*x3;
    shape(9,9) = 0.;
    shape(9,10) = 0.;
    shape(9,11) = 0.;
    shape(9,12) = -1.*x4;
    shape(9,13) = 0.;
    shape(9,14) = 0.;
    shape(9,15) = 0.;

}

void HCurl_PentatopeElement::CalcCurlShape(const IntegrationPoint &ip,
                                                 DenseMatrix &Curlshape) const
{
    const int p = order;
    
#ifdef MFEM_THREAD_SAFE
    Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_t(p + 1),
    shape_l(p + 1);
    DenseMatrix u(Dof, Dim);
#endif
    int size_ip = Nodes.Size();
    //int num_func =
    Curlu.SetSize(dof, int(6));

    //std::cout << "Int Point = " << ip.x << ", " << ip.y << ", " << ip.z << ", " << ip.t << std::endl;
    
    int o = 0;
    int num_dof = dof;
    
    //compute barycentric coordinates as function of ip
    std::vector<double> bary_vector{ip.x, ip.y, ip.z, ip.t, (1.0 - ip.x - ip.y - ip.z - ip.t)};
    
    // compute the gradient of the barycentric coords
    std::vector<double> gradL1{1,0,0,0};
    std::vector<double> gradL2{0,1,0,0};
    std::vector<double> gradL3{0,0,1,0};
    std::vector<double> gradL4{0,0,0,1};
    std::vector<double> gradL5{-1,-1,-1,-1};
    
    std::vector<std::vector<double>> gradbary_vector{gradL1, gradL2, gradL3, gradL4, gradL5};
    
    double La, Lb, Lc, Ld, Le;
    
    std::vector<double> grad_La, grad_Lb, grad_Lc, grad_Ld, grad_Le;
    

    double s_factor = 0.5;

    // Lamda function for computing the skew-symmetric-outer-product of two four vectors
    auto skw_sym_outerprod_fnc = [&](const std::vector<double>& vec1, const std::vector<double>& vec2)
    {

        DenseMatrix result(4,4);
        //std::cout << "Skw-sym Prod" << std::endl;
        // "Forward" Product
        for (int i = 0; i < vec1.size(); ++i)
        {
            for (int j = 0; j < vec2.size(); ++j)
            {
                // Forward Product
                double outer_prod_1 = vec1[i] * vec2[j];
                
                // "Backward" Product
                double outer_prod_2 = vec2[i] * vec1[j];
                
                // Compute Result
                result(i,j) = s_factor*(outer_prod_1 - outer_prod_2);
                //std::cout << 0.5*(outer_prod_1 - outer_prod_2) << std::endl;
            }
            //std::cout << "New COl --------" << std::endl;
        }
        
        return result;
    }; // end of lamda function


                    
   
    
    //Ti.Mult(SkwGradu, SkwGradshape);


}

void HCurl_PentatopeElement::Project (
   VectorCoefficient &vc, ElementTransformation &Trans,
   Vector &dofs) const
{
   Vector v(6);
   double t1[4]; Vector t1i(t1, 4);
   double t2[4]; Vector t2i(t2, 4);
   Vector Mt(4);
   DenseMatrix mat(4,4); mat = 0.0;

   dofs.SetSize(10); dofs = 0.0;
   for (int k = 0; k < 10; k++)
   {
      Trans.SetIntPoint (&Nodes.IntPoint (k));
      const DenseMatrix &J = Trans.Jacobian();

      vc.Eval(v, Trans, Nodes.IntPoint (k));

      mat(0,1) =  v(5); mat(0,2) = -v(4); mat(0,3) =  v(3);
      mat(1,0) = -v(5);                   mat(1,2) =  v(2); mat(1,3) = -v(1);
      mat(2,0) =  v(4); mat(2,1) = -v(2);                   mat(2,3) =  v(0);
      mat(3,0) = -v(3); mat(3,1) =  v(1); mat(3,2) = -v(0);

      J.Mult(tk1[k],t1);
      J.Mult(tk2[k],t2);

      mat.Mult(t2i, Mt);

      dofs(k) = t1i * Mt;
   }
}

//void HCurl_PentatopeElement::Project(const FiniteElement &fe,
//                                        ElementTransformation &Trans,
//                                        DenseMatrix &Id) const
//{
//   int dim = fe.GetDim();
//   int dof = fe.GetDof();
//
//   Id.SetSize(10,6*dof); Id = 0.0;
//
//   double t1[4]; Vector t1i(t1, 4);
//   double t2[4]; Vector t2i(t2, 4);
//   Vector Mt(4);
//   DenseMatrix mat(4,4); mat = 0.0;
//
//   int indI[6] = {0,0,0,1,1,2};
//   int indJ[6] = {1,2,3,2,3,3};
//
//   Vector shape(dof);
//   for (int k = 0; k < 10; k++)
//   {
//      Trans.SetIntPoint(&Nodes.IntPoint(k));
//      const DenseMatrix &J = Trans.Jacobian();
//      fe.CalcShape(Nodes.IntPoint(k), shape);
//
//      J.Mult(tk1[k],t1);
//      J.Mult(tk2[k],t2);
//
//      for (int di=0; di<6; di++)
//      {
//         for (int j=0; j<dof; j++)
//         {
//            mat = 0.0;
//            mat(indI[di], indJ[di]) =  shape(j);
//            mat(indJ[di], indI[di]) = -shape(j);
//
//            mat.Mult(t2i, Mt);
//
//            Id(k, di*dof+j) = t1i * Mt;
//         }
//      }
//   }
//}

void HCurl_PentatopeElement::ProjectDivSkew(const FiniteElement& fe,
                                               ElementTransformation& Trans, DenseMatrix& DivSkew)
{
    mfem_error("SkwGrad_PentatopeElement_Fuentes::ProjectDivSkew not implemented");
    
}
} // end of namespace
