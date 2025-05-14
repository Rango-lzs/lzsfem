#pragma once
#include "FEDomain.h"
#include "elements/FEElement2d.h"

class FEElement2D;

//-----------------------------------------------------------------------------
//! Abstract base class for shell elements
class FEM_EXPORT FEDomain2D : public FEDomain
{
    DECLARE_META_CLASS(FEDomain2D, FEDomain);

public:
    //! constructor
    FEDomain2D(FEModel* fem) : FEDomain(FE_DOMAIN_2D, fem) {}
    
    //! create storage for elements
	bool Create(int nsize, FE_Element_Spec espec) override;

    //! return nr of elements
    int Elements() const override;
    
    //! element access
    FEElement2D& Element(int n);
    FEElement& ElementRef(int n) override;
    const FEElement& ElementRef(int n) const override;

    int GetElementType();
    
    //! Initialize elements
    void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
    //! Reset element data
    void Reset() override;
    
    // inverse jacobian with respect to reference frame
    double invjac0(FEElement2D& el, double J[2][2], int n);
    
    // inverse jacobian with respect to current frame
    double invjact(FEElement2D& el, double J[2][2], int n);
    
    //! calculate in-plane gradient of function at integration points
    Vector2d gradient(FEElement2D& el, double* fn, int n);
    
    //! calculate in-plane gradient of function at integration points
    Vector2d gradient(FEElement2D& el, std::vector<double>& fn, int n);

    //! calculate in-plane gradient of vector function at integration points
    Matrix2d gradient(FEElement2D& el, Vector2d* fn, int n);
    
    //! calculate in-plane gradient of vector function at integration points
    Matrix3d gradient(FEElement2D& el, Vector3d* fn, int n);
    
    // jacobian with respect to reference frame
    double detJ0(FEElement2D& el, int n);
    
    //! calculate jacobian in current frame
    double detJt(FEElement2D& el, int n);
    
    //! calculates covariant basis vectors at an integration point
    void CoBaseVectors(FEElement2D& el, int j, Vector2d g[2]);
    
    //! calculates contravariant basis vectors at an integration point
    void ContraBaseVectors(FEElement2D& el, int j, Vector2d g[2]);
    
    //! calculates parametric derivatives of covariant basis vectors at an integration point
    void CoBaseVectorDerivatives(FEElement2D& el, int j, Vector2d dg[2][2]);
    
    //! calculates parametric derivatives of contravariant basis vectors at an integration point
    void ContraBaseVectorDerivatives(FEElement2D& el, int j, Vector2d dg[2][2]);
    
    //! calculate the laplacian of a vector function at an integration point
    Vector2d lapvec(FEElement2D& el, Vector2d* fn, int n);
    
    //! calculate the gradient of the divergence of a vector function at an integration point
    Vector2d gradivec(FEElement2D& el, Vector2d* fn, int n);
    
protected:
    std::vector<FEElement2D>	m_Elem;	//!< array of elements
};
