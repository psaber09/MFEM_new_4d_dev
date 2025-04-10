// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability see http://mfem.googlecode.com.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#ifndef MFEM_TESSERACT
#define MFEM_TESSERACT

#include "../config/config.hpp"
#include "element.hpp"

namespace mfem
{

/// Data type tesseract element
class Tesseract : public Element
{
protected:
   int indices[16];

public:
   static const int edges[32][2];
   static const int faces[8][8];  // same as Mesh::tes_faces

   Tesseract() : Element(Geometry::TESSERACT) { }

   /// Constructs hexahedron by specifying the indices and the attribute.
   Tesseract(const int *ind, int attr = 1);

   /// Constructs hexahedron by specifying the indices and the attribute.
   Tesseract(int ind1, int ind2, int ind3, int ind4,
             int ind5, int ind6, int ind7, int ind8,
             int ind9, int ind10, int ind11, int ind12,
             int ind13, int ind14, int ind15, int ind16, int attr = 1);

   /// Return element's type
   Type GetType() const override { return Element::TESSERACT; }

   /// Returns the indices of the element's vertices.
   void GetVertices(Array<int> &v) const override;
    
   /// Set the indices defining the vertices.
   void SetVertices(const int *ind) override;

   void SetVertices(const Array<int> &v) override; // Added

   int *GetVertices() override { return indices; }

   int GetNVertices() const override { return 16; }

   int GetNEdges() const override { return 32; }

   int GetNFaces() const override { return 8; }

   int GetNFaceVertices(int) const override { return 8; }

   const int *GetEdgeVertices(int ei) const override
   { return edges[ei]; }
    
//   const int *GetEdgeVertices(int ei) const override
//   { return geom_t::Edges[ei]; }

   int GetNFaces(int &nFaceVertices) const override
   { nFaceVertices = 8; return 8; }

   const int *GetFaceVertices(int fi) const override
   { return faces[fi]; }

   Element *Duplicate(Mesh *m) const override
   { return new Tesseract(indices, attribute); }

   virtual ~Tesseract() = default;
};

extern MFEM_EXPORT class QuadLinear4DFiniteElement TesseractFE;

}

#endif
