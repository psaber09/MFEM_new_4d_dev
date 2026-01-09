//
//  fe_skgrad.hpp
//  mfem
//
//  Created by psaber on 6/20/25.
//

#ifndef fe_curl_hpp
#define fe_curl_hpp

#include <stdio.h>
#include "fe_base.hpp"
#include "fe_h1.hpp"
#include "fe_nd.hpp"
#include "fe_l2.hpp"
#include "fe_pyramid.hpp"

namespace mfem
{

class HCurl_PentatopeElement : public VectorFiniteElement
{
    static const double tk[40], c;
    
#ifndef MFEM_THREAD_SAFE
    mutable Vector shape_x, shape_y, shape_z, shape_t, shape_l;
    mutable Vector dshape_x, dshape_y, dshape_z, dshape_t, dshape_l;
    mutable DenseMatrix u;
    mutable DenseMatrix Curlu;
#endif
    Array<int> dof2tk;
    DenseMatrixInverse Ti;
    
    ND_PentDofTransformation doftrans;

    
public:
    HCurl_PentatopeElement(const int p);
    virtual void CalcVShape(const IntegrationPoint &ip,
                            DenseMatrix &shape) const override;
    virtual void CalcVShape(ElementTransformation &Trans,
                            DenseMatrix &shape) const override
    { CalcVShape_DivSkew(Trans, shape); }
    virtual void CalcCurlShape(const IntegrationPoint &ip,
                                  DenseMatrix &SkwGradshape) const;
    virtual void GetLocalInterpolation(ElementTransformation &Trans,
                                       DenseMatrix &I) const override
    {mfem_error("GetLocalInterpolaton error");}
    //LocalInterpolation_RT(*this, nk, dof2nk, Trans, I); }
    virtual void GetLocalRestriction(ElementTransformation &Trans,
                                     DenseMatrix &R) const override
    {          mfem_error("GetLocalRestriction error");}
    //LocalRestriction_RT(nk, dof2nk, Trans, R); }
    virtual void GetTransferMatrix(const FiniteElement &fe,
                                   ElementTransformation &Trans,
                                   DenseMatrix &I) const override
    {          mfem_error("GetTransferMatrix error");}
    //LocalInterpolation_RT(CheckVectorFE(fe), nk, dof2nk, Trans, I); }
    const StatelessDofTransformation *GetDofTransformation() const override
    { return &doftrans; }
    using FiniteElement::Project;
    virtual void Project(VectorCoefficient &vc,
                         ElementTransformation &Trans, Vector &dofs) const override
    { mfem_error("Error in Project ND referance"); }//Project_ND(tk, dof2tk, vc, Trans, dofs); }
    virtual void ProjectMatrixCoefficient(
                                          MatrixCoefficient &mc, ElementTransformation &T, Vector &dofs) const override
    {          mfem_error("ProjectMatrixCoefficient error");}
    //ProjectMatrixCoefficient_RT(nk, dof2nk, mc, T, dofs); }
    virtual void Project(const FiniteElement &fe, ElementTransformation &Trans,
                         DenseMatrix &I) const override 
    {          mfem_error("Project error");}
    //Project_RT(nk, dof2nk, fe, Trans, I); }
    virtual void ProjectDivSkew(const FiniteElement &fe,
                                ElementTransformation &Trans,
                                DenseMatrix &DivSkew);
};

} // end of namespace
#endif /* fe_curl_hpp */
