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

#include "fem.hpp"

namespace mfem
{

void DofTransformation::TransformPrimal(real_t *v) const
{
   MFEM_ASSERT(dof_trans_,
               "DofTransformation has no local transformation, call "
               "SetDofTransformation first!");
   int size = dof_trans_->Size();
   int ndim = dof_trans_->Dim();

   if (vdim_ == 1 || (Ordering::Type)ordering_ == Ordering::byNODES)
   {
      for (int i=0; i<vdim_; i++)
      {
         if (ndim == 4)
         {
             dof_trans_->TransformPrimal(Po_, Fo_, &v[i*size]);
         }
         else
         {
             dof_trans_->TransformPrimal(Fo_, &v[i*size]);
         }
      }
   }
   else
   {
      Vector vec(size);
      for (int i=0; i<vdim_; i++)
      {
         for (int j=0; j<size; j++)
         {
            vec(j) = v[j*vdim_+i];
         }
         dof_trans_->TransformPrimal(Fo_, vec);
         for (int j=0; j<size; j++)
         {
            v[j*vdim_+i] = vec(j);
         }
      }
   }
}

void DofTransformation::InvTransformPrimal(real_t *v) const
{
   MFEM_ASSERT(dof_trans_,
               "DofTransformation has no local transformation, call "
               "SetDofTransformation first!");
   int size = dof_trans_->Height();
   int dim = dof_trans_->Dim();

   if (vdim_ == 1 || (Ordering::Type)ordering_ == Ordering::byNODES)
   {
      for (int i=0; i<vdim_; i++)
      {
         if (dim == 4)
         {
             dof_trans_->InvTransformPrimal(Po_, Fo_, &v[i*size]);
         }
         else
         {
             dof_trans_->InvTransformPrimal(Fo_, &v[i*size]);
         }
      }
   }
   else
   {
      Vector vec(size);
      for (int i=0; i<vdim_; i++)
      {
         for (int j=0; j<size; j++)
         {
            vec(j) = v[j*vdim_+i];
         }
         dof_trans_->InvTransformPrimal(Fo_, vec);
         for (int j=0; j<size; j++)
         {
            v[j*vdim_+i] = vec(j);
         }
      }
   }
}

void DofTransformation::TransformDual(real_t *v) const
{
   MFEM_ASSERT(dof_trans_,
               "DofTransformation has no local transformation, call "
               "SetDofTransformation first!");
   int size = dof_trans_->Size();
   int ndim = dof_trans_->Dim();

   if (vdim_ == 1 || (Ordering::Type)ordering_ == Ordering::byNODES)
   {
      for (int i=0; i<vdim_; i++)
      {
          if (ndim == 4)
          {
              dof_trans_->TransformDual(Po_, Fo_, &v[i*size]);
          }
          else
          {
              dof_trans_->TransformDual(Fo_, &v[i*size]);
          }
      }
   }
   else
   {
      Vector vec(size);
      for (int i=0; i<vdim_; i++)
      {
         for (int j=0; j<size; j++)
         {
            vec(j) = v[j*vdim_+i];
         }
         dof_trans_->TransformDual(Fo_, vec);
         for (int j=0; j<size; j++)
         {
            v[j*vdim_+i] = vec(j);
         }
      }
   }
}

void DofTransformation::InvTransformDual(real_t *v) const
{
   MFEM_ASSERT(dof_trans_,
               "DofTransformation has no local transformation, call "
               "SetDofTransformation first!");
   int size = dof_trans_->Size();
   int ndim = dof_trans_->Dim();

   if (vdim_ == 1 || (Ordering::Type)ordering_ == Ordering::byNODES)
   {
      for (int i=0; i<vdim_; i++)
      {
          if (ndim == 4)
          {
              dof_trans_->InvTransformDual(Po_, Fo_, &v[i*size]);
          }
          else
          {
              dof_trans_->InvTransformDual(Fo_, &v[i*size]);
          }      }
   }
   else
   {
      Vector vec(size);
      for (int i=0; i<vdim_; i++)
      {
         for (int j=0; j<size; j++)
         {
            vec(j) = v[j*vdim_+i];
         }
         dof_trans_->InvTransformDual(Fo_, vec);
         for (int j=0; j<size; j++)
         {
            v[j*vdim_+i] = vec(j);
         }
      }
   }
}

void TransformPrimal(const DofTransformation *ran_dof_trans,
                     const DofTransformation *dom_dof_trans,
                     DenseMatrix &elmat)
{
   // No action if both transformations are NULL
   if (ran_dof_trans)
   {
      ran_dof_trans->TransformPrimalCols(elmat);
   }
   if (dom_dof_trans)
   {
      dom_dof_trans->TransformDualRows(elmat);
   }
}

void TransformDual(const DofTransformation *ran_dof_trans,
                   const DofTransformation *dom_dof_trans,
                   DenseMatrix &elmat)
{
   // No action if both transformations are NULL
   if (ran_dof_trans)
   {
      ran_dof_trans->TransformDualCols(elmat);
   }
   if (dom_dof_trans)
   {
      dom_dof_trans->TransformDualRows(elmat);
   }
}

// ordering (i0j0, i1j0, i0j1, i1j1), each row is a column major matrix
const real_t ND_DofTransformation::T_data[24] =
{
   1.0,  0.0,  0.0,  1.0,
   -1.0, -1.0,  0.0,  1.0,
   0.0,  1.0, -1.0, -1.0,
   1.0,  0.0, -1.0, -1.0,
   -1.0, -1.0,  1.0,  0.0,
   0.0,  1.0,  1.0,  0.0
};

const DenseTensor ND_DofTransformation
::T(const_cast<real_t *>(ND_DofTransformation::T_data), 2, 2, 6);

// ordering (i0j0, i1j0, i0j1, i1j1), each row is a column major matrix
const real_t ND_DofTransformation::TInv_data[24] =
{
   1.0,  0.0,  0.0,  1.0,
   -1.0, -1.0,  0.0,  1.0,
   -1.0, -1.0,  1.0,  0.0,
   1.0,  0.0, -1.0, -1.0,
   0.0,  1.0, -1.0, -1.0,
   0.0,  1.0,  1.0,  0.0
};

const DenseTensor ND_DofTransformation
::TInv(const_cast<real_t *>(TInv_data), 2, 2, 6);

// ordering (i0j0, i1j0, i0j1, i1j1), each row is a column major matrix
const real_t ND_DofTransformation::T_data4D[216] =
{
   1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0,  0.0,  1.0,
   1.0, 0.0, 0.0,  0.0, 0.0, 1.0,  0.0,  1.0,  0.0,
   0.0, 0.0, 1.0,  1.0, 0.0, 0.0,  0.0,  1.0,  0.0,
   0.0, 1.0, 0.0,  1.0, 0.0, 0.0,  0.0,  0.0,  1.0,
   0.0, 1.0, 0.0,  0.0, 0.0, 1.0,  1.0,  0.0,  0.0,
   0.0, 0.0, 1.0,  0.0, 1.0, 0.0,  1.0,  0.0,  0.0,
   -1.0, -1.0, -1.0,  1.0, 0.0, 0.0,  0.0,  0.0,  1.0,
   -1.0, -1.0, -1.0,  1.0, 0.0, 0.0,  0.0,  1.0,  0.0,
   -1.0, -1.0, -1.0,  0.0, 1.0, 0.0,  1.0,  0.0,  0.0,
   -1.0, -1.0, -1.0,  0.0, 0.0, 1.0,  1.0,  0.0,  0.0,
   -1.0, -1.0, -1.0,  0.0, 0.0, 1.0,  0.0,  1.0,  0.0,
   -1.0, -1.0, -1.0,  0.0, 1.0, 0.0,  0.0,  0.0,  1.0,
   0.0, 0.0, 1.0,  -1.0, -1.0, -1.0,  1.0,  0.0,  0.0,
   0.0, 1.0, 0.0,  -1.0, -1.0, -1.0,  1.0,  0.0,  0.0,
   0.0, 1.0, 0.0,  -1.0, -1.0, -1.0,  0.0,  0.0,  1.0,
   0.0, 0.0, 1.0,  -1.0, -1.0, -1.0,  0.0,  1.0,  0.0,
   1.0, 0.0, 0.0,  -1.0, -1.0, -1.0,  0.0,  1.0,  0.0,
   1.0, 0.0, 0.0,  -1.0, -1.0, -1.0,  0.0,  0.0,  1.0,
   0.0, 0.0, 1.0,  0.0, 1.0, 0.0,  -1.0,  -1.0,  -1.0,
   0.0, 1.0, 0.0,  0.0, 0.0, 1.0,  -1.0,  -1.0,  -1.0,
   1.0, 0.0, 0.0,  0.0, 0.0, 1.0,  -1.0,  -1.0,  -1.0,
   1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  -1.0,  -1.0,  -1.0,
   0.0, 1.0, 0.0,  1.0, 0.0, 0.0,  -1.0,  -1.0,  -1.0,
   0.0, 0.0, 1.0,  1.0, 0.0, 0.0,  -1.0,  -1.0,  -1.0,

};

const DenseTensor ND_DofTransformation
::T4D(const_cast<real_t *>(ND_DofTransformation::T_data4D), 3, 3, 24);

// ordering (i0j0, i1j0, i0j1, i1j1), each row is a column major matrix
const real_t ND_DofTransformation::TInv_data4D[216] =
{
    1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0,  0.0,  1.0,
    1.0, 0.0, 0.0,  0.0, 0.0, 1.0,  0.0,  1.0,  0.0,
    0.0, 1.0, 0.0,  0.0, 0.0, 1.0,  1.0,  0.0,  0.0,
    0.0, 1.0, 0.0,  1.0, 0.0, 0.0,  0.0,  0.0,  1.0,
    0.0, 0.0, 1.0,  1.0, 0.0, 0.0,  0.0,  1.0,  0.0,
    0.0, 0.0, 1.0,  0.0, 1.0, 0.0,  1.0,  0.0,  0.0,
    0.0, 1.0, 0.0,  -1.0, -1.0, -1.0,  0.0,  0.0,  1.0,
    0.0, 1.0, 0.0,  0.0, 0.0, 1.0,  -1.0,  -1.0,  -1.0,
    0.0, 0.0, 1.0,  0.0, 1.0, 0.0,  -1.0,  -1.0,  -1.0,
    0.0, 0.0, 1.0,  -1.0, -1.0, -1.0,  0.0,  1.0,  0.0,
    -1.0, -1.0, -1.0,  0.0, 0.0, 1.0,  0.0,  1.0,  0.0,
    -1.0, -1.0, -1.0,  0.0, 1.0, 0.0,  0.0,  0.0,  1.0,
    0.0, 0.0, 1.0,  -1.0, -1.0, -1.0,  1.0,  0.0,  0.0,
    0.0, 0.0, 1.0,  1.0, 0.0, 0.0,  -1.0,  -1.0,  -1.0,
    -1.0, -1.0, -1.0,  1.0, 0.0, 0.0,  0.0,  0.0,  1.0,
    -1.0, -1.0, -1.0,  0.0, 0.0, 1.0,  1.0,  0.0,  0.0,
    1.0, 0.0, 0.0,  0.0, 0.0, 1.0,  -1.0,  -1.0,  -1.0,
    1.0, 0.0, 0.0,  -1.0, -1.0, -1.0,  0.0,  0.0,  1.0,
    -1.0, -1.0, -1.0,  0.0, 1.0, 0.0,  1.0,  0.0,  0.0,
    -1.0, -1.0, -1.0,  1.0, 0.0, 0.0,  0.0,  1.0,  0.0,
    1.0, 0.0, 0.0,  -1.0, -1.0, -1.0,  0.0,  1.0,  0.0,
    1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  -1.0,  -1.0,  -1.0,
    0.0, 1.0, 0.0,  1.0, 0.0, 0.0,  -1.0,  -1.0,  -1.0,
    0.0, 1.0, 0.0,  -1.0, -1.0, -1.0,  1.0,  0.0,  0.0,
    
};

const DenseTensor ND_DofTransformation
::TInv4D(const_cast<real_t *>(TInv_data4D), 3, 3, 24);


ND_DofTransformation::ND_DofTransformation(int size, int ndim, int p, int num_edges,
                                           int num_faces,
                                           int face_types[])
   : StatelessDofTransformation(size, ndim)
   , order(p)
   , nedofs(p)
   , ntdofs(p*(p-1))
   //, ntetdofs(0) // not sure
   , ntetdofs(p*(p-1)*(p-2)*0.5)
   , nqdofs(2*p*(p-1))
   , nedges(num_edges)
   , nplanars(0) // not sure
   , nfaces(num_faces)
   , ftypes(face_types)
{
}

ND_DofTransformation::ND_DofTransformation(int size, int ndim, int p, int num_edges, int num_planars,
                                           int num_faces,
                                           int planar_types[], int face_types[])
   : StatelessDofTransformation(size, ndim)
   , order(p)
   , nedofs(p)
   , ntdofs(p*(p-1))
   , nqdofs(2*p*(p-1))
   , ntetdofs(p*(p-1)*(p-2)*0.5)
   , nedges(num_edges)
   , nplanars(num_planars)
   , nfaces(num_faces)
   , ptypes(planar_types)
   , ftypes(face_types)
{
}

void ND_DofTransformation::TransformPrimal(const Array<int> & Fo,
                                           real_t *v) const
{
   // Return immediately when no face DoFs are present
   if (IsIdentity()) { return; }

//   MFEM_VERIFY(Fo.Size() >= nfaces,
//               "Face orientation array is shorter than the number of faces in "
//               "ND_DofTransformation");

   int of = 0;
   real_t data[2];
   Vector v2(data, 2);
   DenseMatrix T2;
    
        std::cout << "Face 3D Oreint 1 " << Fo[0] << std::endl;
        std::cout << "Face 3D Oreint 2 " << Fo[1] << std::endl;
        std::cout << "Face 3D Oreint 3 " << Fo[2] << std::endl;
        std::cout << "Face 3D Oreint 4 " << Fo[3] << std::endl;

   // Transform face DoFs
   for (int f=0; f<nfaces; f++)
   {
      if (ftypes[f] == Geometry::TRIANGLE)
      {
         for (int i=0; i<ntdofs/2; i++)
         {
            v2 = &v[nedges*nedofs + of + 2*i];
            T2.UseExternalData(const_cast<real_t *>(T.GetData(Fo[f])), 2, 2);
            T2.Mult(v2, &v[nedges*nedofs + of + 2*i]);
         }
         of += ntdofs;
      }
      else
      {
         of += nqdofs;
      }
   }
}

void ND_DofTransformation::TransformPrimal(const Array<int> & Po, const Array<int> & Fo,
                                           real_t *v) const
{
   // Return immediately when no face DoFs are present
   //if (IsIdentity()) { return; }

//   MFEM_VERIFY(Fo.Size() >= nfaces,
//               "Face orientation array is shorter than the number of faces in "
//               "ND_DofTransformation");

   int of_pl = 0;
   int of_ft = 0;
   real_t data[2];
   real_t data4D[3];
   Vector v2(data, 2);
   Vector v3(data4D,3);
   DenseMatrix T2_planar;
   DenseMatrix T2_facet;
    
//    std::cout << "Planar Oreint 1 " << Po[0] << std::endl;
//    std::cout << "Planar Oreint 2 " << Po[1] << std::endl;
//    std::cout << "Planar Oreint 3 " << Po[2] << std::endl;
//    std::cout << "Planar Oreint 4 " << Po[3] << std::endl;
//    std::cout << "Planar Oreint 5 " << Po[4] << std::endl;
//    std::cout << "Planar Oreint 6 " << Po[5] << std::endl;
//    std::cout << "Planar Oreint 7 " << Po[6] << std::endl;
//    std::cout << "Planar Oreint 8 " << Po[7] << std::endl;
//    std::cout << "Planar Oreint 9 " << Po[8] << std::endl;
//    std::cout << "Planar Oreint 10 " << Po[9] << std::endl;

   // Transform Planar DoFs
   for (int pl=0; pl<nplanars; pl++)
   {
      if (ptypes[pl] == Geometry::TRIANGLE)
      {
         for (int i=0; i<ntdofs/2; i++)
         {
            v2 = &v[nedges*nedofs + of_pl + 2*i];
            T2_planar.UseExternalData(const_cast<real_t *>(T.GetData(Po[pl])), 2, 2);
            T2_planar.Mult(v2, &v[nedges*nedofs + of_pl + 2*i]);
         }
         of_pl += ntdofs;
      }
      else
      {
         of_pl += nqdofs;
      }
   }
    
    // Transform Facet DoFs
//    std::cout << "Primal Transform -----" << std::endl;
//    std::cout << "Facet Oreint 1 " << Fo[0] << std::endl;
//    std::cout << "Facet Oreint 2 " << Fo[1] << std::endl;
//    std::cout << "Facet Oreint 3 " << Fo[2] << std::endl;
//    std::cout << "Facet Oreint 4 " << Fo[3] << std::endl;
//    std::cout << "Facet Oreint 5 " << Fo[4] << std::endl;
//    if (Fo[0] > 23) {
//        mfem_error("Invalid Orientation index");
//    }
//    if (Fo[1] > 23) {
//        mfem_error("Invalid Orientation index");
//    }
//    if (Fo[2] > 23) {
//        mfem_error("Invalid Orientation index");
//    }
//    if (Fo[3] > 23) {
//        mfem_error("Invalid Orientation index");
//    }
//    if (Fo[4] > 23) {
//        mfem_error("Invalid Orientation index");
//    }

    for (int ft=0; ft<nfaces; ft++)
    {
       if (ftypes[ft] == Geometry::TETRAHEDRON)
       {
          for (int i=0; i<ntetdofs/3; i++)
          {
             v3 = &v[nedges*nedofs + nplanars*ntdofs + of_ft + 3*i];
             T2_facet.UseExternalData(const_cast<real_t *>(T4D.GetData(Fo[ft])), 3, 3);
             //T2_facet.Print(std::cout);
             T2_facet.Mult(v3, &v[nedges*nedofs + nplanars*ntdofs + of_ft + 3*i]);
          }
          of_ft += ntetdofs;
       }
       else
       {
           mfem_error("Unsupported Facet Type");
       }
    }
    
}

void ND_DofTransformation::InvTransformPrimal(const Array<int> & Fo,
                                              real_t *v) const
{
   // Return immediately when no face DoFs are present
   if (IsIdentity()) { return; }

   MFEM_VERIFY(Fo.Size() >= nfaces,
               "Face orientation array is shorter than the number of faces in "
               "ND_DofTransformation");

   int of = 0;
   real_t data[2];
   Vector v2(data, 2);
   DenseMatrix T2Inv;

   // Transform face DoFs
   for (int f=0; f<nfaces; f++)
   {
      if (ftypes[f] == Geometry::TRIANGLE)
      {
         for (int i=0; i<ntdofs/2; i++)
         {
            v2 = &v[nedges*nedofs + of + 2*i];
            T2Inv.UseExternalData(const_cast<real_t *>(TInv.GetData(Fo[f])), 2, 2);
            T2Inv.Mult(v2, &v[nedges*nedofs + of + 2*i]);
         }
         of += ntdofs;
      }
      else
      {
         of += nqdofs;
      }
   }
}

void ND_DofTransformation::InvTransformPrimal(const Array<int> & Po, const Array<int> & Fo,
                                              real_t *v) const
{
   // Return immediately when no face DoFs are present
   //if (IsIdentity()) { return; }
//
//   MFEM_VERIFY(Fo.Size() >= nfaces,
//               "Face orientation array is shorter than the number of faces in "
//               "ND_DofTransformation");

   int of_pl = 0;
   int of_ft = 0;
   real_t data[2];
   real_t data4D[3];
   Vector v2(data, 2);
   Vector v3(data4D, 3);
   DenseMatrix T2Inv_planar;
   DenseMatrix T2Inv_facet;

   // Transform planar DoFs
   for (int pl=0; pl<nplanars; pl++)
   {
      if (ptypes[pl] == Geometry::TRIANGLE)
      {
         for (int i=0; i<ntdofs/2; i++)
         {
            v2 = &v[nedges*nedofs + of_pl + 2*i];
            T2Inv_planar.UseExternalData(const_cast<real_t *>(TInv.GetData(Po[pl])), 2, 2);
            T2Inv_planar.Mult(v2, &v[nedges*nedofs + of_pl + 2*i]);
         }
         of_pl += ntdofs;
      }
      else
      {
         of_pl += nqdofs;
      }
   }
    
    // Transform facet DoFs
    for (int ft=0; ft<nfaces; ft++)
    {
       if (ftypes[ft] == Geometry::TETRAHEDRON)
       {
          for (int i=0; i<ntetdofs/3; i++)
          {
             v3 = &v[nedges*nedofs + nplanars*ntdofs + of_ft + 3*i];
             T2Inv_facet.UseExternalData(const_cast<real_t *>(TInv4D.GetData(Fo[ft])), 3, 3);
             T2Inv_facet.Mult(v3, &v[nedges*nedofs + nplanars*ntdofs + of_ft + 3*i]);
          }
          of_ft += ntetdofs;
       }
       else
       {
           mfem_error("Unsupported Facet Type");
       }
    }

}

void ND_DofTransformation::TransformDual(const Array<int> & Fo, real_t *v) const
{
   // Return immediately when no face DoFs are present
   if (IsIdentity()) { return; }

   MFEM_VERIFY(Fo.Size() >= nfaces,
               "Face orientation array is shorter than the number of faces in "
               "ND_DofTransformation");

   int of = 0;
   real_t data[2];
   Vector v2(data, 2);
   DenseMatrix T2Inv;

   // Transform face DoFs
   for (int f=0; f<nfaces; f++)
   {
      if (ftypes[f] == Geometry::TRIANGLE)
      {
         for (int i=0; i<ntdofs/2; i++)
         {
            v2 = &v[nedges*nedofs + of + 2*i];
            T2Inv.UseExternalData(const_cast<real_t *>(TInv.GetData(Fo[f])), 2, 2);
            T2Inv.MultTranspose(v2, &v[nedges*nedofs + of + 2*i]);
         }
         of += ntdofs;
      }
      else
      {
         of += nqdofs;
      }

   }
}

void ND_DofTransformation::TransformDual(const Array<int> & Po, const Array<int> & Fo, real_t *v) const
{
   // Return immediately when no face DoFs are present
   if (IsIdentity()) { return; }

   MFEM_VERIFY(Fo.Size() >= nfaces,
               "Face orientation array is shorter than the number of faces in "
               "ND_DofTransformation");

   int of_pl = 0;
   int of_ft = 0;
   real_t data[2];
   real_t data4D[3];
   Vector v2(data, 2);
   Vector v3(data4D, 3);
   DenseMatrix T2Inv_planars;
   DenseMatrix T2Inv_factets;

   // Transform planar DoFs
   for (int pl=0; pl<nplanars; pl++)
   {
      if (ptypes[pl] == Geometry::TRIANGLE)
      {
         for (int i=0; i<ntdofs/2; i++)
         {
            v2 = &v[nedges*nedofs + of_pl + 2*i];
            T2Inv_planars.UseExternalData(const_cast<real_t *>(TInv.GetData(Po[pl])), 2, 2);
            T2Inv_planars.MultTranspose(v2, &v[nedges*nedofs + of_pl + 2*i]);
         }
         of_pl += ntdofs;
      }
      else
      {
         of_pl += nqdofs;
      }

   }
    
    // Transform face DoFs
    for (int ft=0; ft<nfaces; ft++)
    {
       if (ftypes[ft] == Geometry::TETRAHEDRON)
       {
          for (int i=0; i<ntetdofs/3; i++)
          {
             v3 = &v[nedges*nedofs +nplanars*ntdofs + of_ft + 3*i];
             T2Inv_factets.UseExternalData(const_cast<real_t *>(TInv4D.GetData(Fo[ft])), 3, 3);
             T2Inv_factets.MultTranspose(v3, &v[nedges*nedofs + nplanars*ntdofs + of_ft + 3*i]);
          }
          of_ft += ntetdofs;
       }
       else
       {
           mfem_error("Unsupported Facet Type");
       }

    }
}

void ND_DofTransformation::InvTransformDual(const Array<int> & Fo,
                                            real_t *v) const
{
   // Return immediately when no face DoFs are present
   if (IsIdentity()) { return; }

   MFEM_VERIFY(Fo.Size() >= nfaces,
               "Face orientation array is shorter than the number of faces in "
               "ND_DofTransformation");

   int of = 0;
   real_t data[2];
   Vector v2(data, 2);
   DenseMatrix T2;

   // Transform face DoFs
   for (int f=0; f<nfaces; f++)
   {
      if (ftypes[f] == Geometry::TRIANGLE)
      {
         for (int i=0; i<ntdofs/2; i++)
         {
            v2 = &v[nedges*nedofs + of + 2*i];
            T2.UseExternalData(const_cast<real_t *>(T.GetData(Fo[f])), 2, 2);
            T2.MultTranspose(v2, &v[nedges*nedofs + of + 2*i]);
         }
         of += ntdofs;
      }
      else
      {
         of += nqdofs;
      }
   }
}

void ND_DofTransformation::InvTransformDual(const Array<int> & Po, const Array<int> & Fo,
                                            real_t *v) const
{
   // Return immediately when no face DoFs are present
   if (IsIdentity()) { return; }

   MFEM_VERIFY(Fo.Size() >= nfaces,
               "Face orientation array is shorter than the number of faces in "
               "ND_DofTransformation");

   int of_pl = 0;
   int of_ft = 0;
   real_t data[2];
   real_t data4D[3];
   Vector v2(data, 2);
   Vector v3(data4D, 3);
   DenseMatrix T2_planars;
   DenseMatrix T2_facets;

   // Transform planar DoFs
   for (int pl=0; pl<nplanars; pl++)
   {
      if (ptypes[pl] == Geometry::TRIANGLE)
      {
         for (int i=0; i<ntdofs/2; i++)
         {
            v2 = &v[nedges*nedofs + of_pl + 2*i];
            T2_planars.UseExternalData(const_cast<real_t *>(T.GetData(Po[pl])), 2, 2);
            T2_planars.MultTranspose(v2, &v[nedges*nedofs + of_pl + 2*i]);
         }
         of_pl += ntdofs;
      }
      else
      {
         of_pl += nqdofs;
      }
   }
    
    // Transform face DoFs
    for (int ft=0; ft<nfaces; ft++)
    {
       if (ftypes[ft] == Geometry::TETRAHEDRON)
       {
          for (int i=0; i<ntetdofs/3; i++)
          {
             v3 = &v[nedges*nedofs + nplanars*ntdofs + of_ft + 3*i];
             T2_facets.UseExternalData(const_cast<real_t *>(T4D.GetData(Fo[ft])), 3, 3);
             T2_facets.MultTranspose(v3, &v[nedges*nedofs +nplanars*ntdofs + of_ft + 3*i]);
          }
          of_ft += ntetdofs;
       }
       else
       {
           mfem_error("Unsupported Facet Type");
       }
    }
}


} // namespace mfem
