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

const double HCurl_PentatopeElement::tk[40] =
{ 1,0,0,0,  0,1,0,0,  0,0,1,0,  0,0,0,1,  -1,1,0,0,  -1,0,1,0,  0,-1,1,0,  -1,0,0,1,  0,-1,0,1,  0,0,-1,1 };

//const double HSkwGrad_PentatopeElement::tk[40] =
//{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1, -1,1,0,0, -1,0,1,0, -1,0,0,1, 0,-1,1,0, 0,-1,0,1, 0,0,-1,1};

const double HCurl_PentatopeElement::c = 1./5.;

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
          Nodes.IntPoint(o).Set4(fop[i]/w, fop[j]/w, 0.0, 0.0);
          dof2tk[o++] = 0;
          //std::cout << "inside dof" << std::endl;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (0,1,3)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(fop[i]/w, 0.0, fop[j]/w, 0.0);
          dof2tk[o++] = 0;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (0,1,4)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(fop[i]/w, 0.0, 0.0, fop[j]/w);
          dof2tk[o++] = 0;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (0,2,3)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(0.0, fop[i]/w, fop[j]/w, 0.0);
          dof2tk[o++] = 1;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (0,2,4)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(0.0, fop[i]/w, 0.0, fop[j]/w);
          dof2tk[o++] = 1;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (0,3,4)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(0.0, 0.0, fop[i]/w, fop[j]/w);
          dof2tk[o++] = 2;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (1,2,3)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(fop[pm1-i-j]/w, fop[i]/w, fop[j]/w, 0.0);
          dof2tk[o++] = 4;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (1,2,4)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(fop[pm1-i-j]/w, fop[i]/w, 0.0, fop[j]/w);
          dof2tk[o++] = 4;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (1,3,4)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(fop[pm1-i-j]/w, 0.0, fop[i]/w, fop[j]/w);
          dof2tk[o++] = 5;
       }
    for (int j = 0; j < p; j++)
       for (int i=0; i + j < p; i++) // (2,3,4)
       {
          double w = fop[i] + fop[j] + fop[pm1-i-j];
          Nodes.IntPoint(o).Set4(0.0, fop[pm1-i-j]/w, fop[i]/w, fop[j]/w);
          dof2tk[o++] = 6;
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
    
    DenseMatrix T(dof);
    DenseMatrix B(dof, dim);
    
    for (int q = 0; q < dof; q++)
    {
        const IntegrationPoint &ip = Nodes.IntPoint(q);
        //std::cout << "Dof " << q << " = " << ip.x << ", " << ip.y << ", " << ip.z << ", " << ip.t << std::endl;
        
        //const double *nm = nk + 4*dof2nk[m];
        
        const Vector tm({tk[4*dof2tk[q]], tk[4*dof2tk[q]+1], tk[4*dof2tk[q]+2], tk[4*dof2tk[q]+3]});
        //std::cout << "tm < " << tm(0) << "," << tm(1) << "," << tm(2) << "," << tm(3) << ">" <<std::endl;
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
        
        
        int edge_counter = 0;
        //Edges
        for(int i=0; i<p; i++)
        {
            for(int a=0; a<5; a++)
            {
                for(int b=0; b<5; b++)
                {
                    if(a<b)
                    {
                        La = bary_vector[a];
                        Lb = bary_vector[b];
                        
                        edge_counter = edge_counter + 1;
                        
                        // compute polynomials
                        std::vector<double> Legendre_i;
                        double x = Lb;
                        double y = La + Lb;
                        poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                        
                        // Whitney Function
                        //std::vector<double> Whit_Vec = WH_func(a,b);
//                        
//                        // Add Basis Funcitons
//                        double x_comp = Legendre_i[Legendre_i.size()-1]*Whit_Vec[0];
//                        double Whitx = Whit_Vec[0];
//                        B(o, 0) = x_comp;
//                        
//                        double y_comp = Legendre_i[Legendre_i.size()-1]*Whit_Vec[1];
//                        double Whity = Whit_Vec[1];
//                        
//                        B(o, 1) = y_comp;
//                        
//                        double z_comp = Legendre_i[Legendre_i.size()-1]*Whit_Vec[2];
//                        double Whitz = Whit_Vec[2];
//                        
//                        B(o, 2) = z_comp;
//                        
//                        double t_comp = Legendre_i[Legendre_i.size()-1]*Whit_Vec[3];
//                        double Whitt = Whit_Vec[3];
//                        
//                        B(o, 3) = t_comp;
                        
                        o++;
                        
                    }
                }
            }
        }  // end of edges
    }
        

    
    Ti.Factor(T);
    
    
//    std::ofstream Ti_file("Ti_matrix.txt");
//    Ti.PrintMatlab(Ti_file);

    mfem::out << "HCurl_PentatopeElement(" << p << ") : "; Ti.TestInversion();
}

void HCurl_PentatopeElement::CalcVShape(const IntegrationPoint &ip,
                                           DenseMatrix &shape) const
{
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
        for(int j=1; j<p;j++)
        {
            for(int a=0; a<5;a++)
            {
                for(int b=0; b<5;b++)
                {
                    for(int c=0; c<5;c++)
                    {
                        if((a<b)&&(b<c)&&((i+j)<p))
                        {
                            // Define Family
                            int Family = 1;
                            
                            if (Family == 1)
                            {
                                
                                // Family I:
                                La = bary_vector[a];
                                Lb = bary_vector[b];
                                Lc = bary_vector[c];
                                
                                // compute polynomials
                                std::vector<double> Legendre_i;
                                double x = Lb;
                                double y = La + Lb;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                std::vector<double> Int_Jacobi_j;
                                x = Lc;
                                y = La + Lb + Lc;
                                double alpha = 2*i +1;
                                poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                                
                                // Whitney Function
                                //std::vector<double> Whit_Vec = WH_func(a,b);
//                                
//                                // Add Basis Funcitons
//                                u(o, 0) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Whit_Vec[0];
//                                
//                                u(o, 1) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Whit_Vec[1];
//                                
//                                u(o, 2) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Whit_Vec[2];
//                                
//                                u(o, 3) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Whit_Vec[3];
                                
                                o++;
                                
                            }
                            
                            Family++;
                            // Clean up storage for Polynomials
                            //Legendre_i.clear();
                            //Int_Jacobi_j.clear();
                            
                            if (Family == 2)
                            {
                                
                                // Family I:
                                La = bary_vector[b];
                                Lb = bary_vector[c];
                                Lc = bary_vector[a];
                                
                                // compute polynomials
                                std::vector<double> Legendre_i;
                                double x = Lb;
                                double y = La + Lb;
                                poly1d.CalcLegendreShifted(i, x, y, Legendre_i);
                                
                                std::vector<double> Int_Jacobi_j;
                                x = Lc;
                                y = La + Lb + Lc;
                                double alpha = 2*i +1;
                                poly1d.CalcIntJacobi(j, x, y, alpha, Int_Jacobi_j);
                                
//                                // Whitney Function
//                                std::vector<double> Whit_Vec = WH_func(b,c);
//                                
//                                // Add Basis Funcitons
//                                u(o, 0) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Whit_Vec[0];
//                                
//                                u(o, 1) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Whit_Vec[1];
//                                
//                                u(o, 2) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Whit_Vec[2];
//                                
//                                u(o, 3) = Legendre_i[Legendre_i.size()-1] * Int_Jacobi_j[Int_Jacobi_j.size()-1] * Whit_Vec[3];
//                                
                                o++;
                                
                            }
                            
                            
                        }
                    }
                }
            }
        }
    }
    // end of faces
    
    Ti.Mult(u, shape);

    
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

void HCurl_PentatopeElement::ProjectDivSkew(const FiniteElement& fe,
                                               ElementTransformation& Trans, DenseMatrix& DivSkew)
{
    mfem_error("SkwGrad_PentatopeElement_Fuentes::ProjectDivSkew not implemented");
    
}
} // end of namespace
