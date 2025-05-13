#pragma once
#include "FEDomain.h"
#include "femcore/FEDofList.h"
#include "femcore/FELinearSystem.h"
#include "elements/FESolidElement.h"
#include <functional>
#include "datastructure/tens3d.h"

//-----------------------------------------------------------------------------
// This typedef defines a surface integrand. 
// It evaluates the function at surface material point mp, and returns the value
// it the val std::vector. The size of the std::vector is determined by the field variable
// that is being integrated and is already set when the integrand is called.
// This is used in the FESurface::LoadVector function.
typedef std::function<void(FEMaterialPoint& mp, int node_a, std::vector<double>& val)> FEVolumeVectorIntegrand;

typedef std::function<void(FEMaterialPoint& mp, int node_a, int node_b, Matrix& val)> FEVolumeMatrixIntegrand;

//-----------------------------------------------------------------------------
//! abstract base class for 3D volumetric elements
class FEM_EXPORT FESolidDomain : public FEDomain
{
    DECLARE_META_CLASS(FESolidDomain, FEDomain);

public:
    //! constructor
    FESolidDomain(FEModel* pfem);
    
    //! create storage for elements
	bool Create(int nsize, FE_Element_Spec espec) override;

    //! return nr of elements
	int Elements() const override;

	//! initialize element data
	bool Init() override;
    
    //! reset data (overridden from FEDomain)
    void Reset() override;
    
    //! copy data from another domain (overridden from FEDomain)
    void CopyFrom(FEMeshPartition* pd) override;

    //! element access
	FESolidElement& Element(int n);
    FEElement& ElementRef(int n) override { return m_Elem[n]; }
	const FEElement& ElementRef(int n) const override { return m_Elem[n]; }

    ElementType GetElementType() const { return m_Elem[0].elementType(); }
    
    int GetElementShape() const { return m_Elem[0].ShapeFunctions(-1); }

	FE_Element_Spec GetElementSpec() const;
    
    //! find the element in which point y lies
    FESolidElement* FindElement(const Vector3d& y, double r[3]);

	//! find the element in which point y lies (reference configuration)
	FESolidElement* FindReferenceElement(const Vector3d& y, double r[3]);

	//! Project a point to an element and return natural coordinates
	void ProjectToElement(FESolidElement& el, const Vector3d& p, double r[3]);

	//! Project a point to an element in the reference frame and return natural coordinates
	//! returns true if the point lies in the element
	bool ProjectToReferenceElement(FESolidElement& el, const Vector3d& p, double r[3]);

    //! Calculate deformation gradient at integration point n
    double defgrad(FESolidElement& el, Matrix3d& F, int n);
	double defgrad(FESolidElement& el, Matrix3d& F, int n, Vector3d* r);

    //! Calculate deformation gradient at integration point n
    double defgrad(FESolidElement& el, Matrix3d& F, double r, double s, double t);
    
    //! Calculate deformation gradient at integration point n at previous time
    double defgradp(FESolidElement& el, Matrix3d& F, int n);
    
    //! Calculate GradJ at integration point n at current time
    Vector3d GradJ(FESolidElement& el, int n);
    
    //! Calculate GradJ at integration point n at prev time
    Vector3d GradJp(FESolidElement& el, int n);
    
    //! Calculate GradJ at integration point n at current time
    tens3dls Gradb(FESolidElement& el, int n);
    
    //! Calculate GradJ at integration point n at prev time
    tens3dls Gradbp(FESolidElement& el, int n);
    
    //! calculate inverse jacobian matrix w.r.t. reference frame
    double invjac0(const FESolidElement& el, double J[3][3], int n);
    
    //! calculate inverse jacobian matrix w.r.t. reference frame
    double invjac0(const FESolidElement& el, double J[3][3], double r, double s, double t);
    
    //! calculate inverse jacobian matrix w.r.t. reference frame
    double invjac0(const FESolidElement& el, double r, double s, double t, Matrix3d& J);
    
    //! calculate inverse jacobian matrix w.r.t. current frame
    double invjact(FESolidElement& el, double J[3][3], int n);
	double invjact(FESolidElement& el, double J[3][3], int n, const Vector3d* r);

    //! calculate inverse jacobian matrix w.r.t. reference frame
    double invjact(FESolidElement& el, double J[3][3], double r, double s, double t);
    
    //! calculate inverse jacobian matrix w.r.t. intermediate frame
    double invjact(FESolidElement& el, double J[3][3], int n, const double alpha);
    
    //! calculate inverse jacobian matrix w.r.t. previous time
    double invjactp(FESolidElement& el, double J[3][3], int n);
    
    //! calculate gradient of function at integration points
    Vector3d gradient(FESolidElement& el, double* fn, int n);
	Vector3d gradient(FESolidElement& el, int order, double* fn, int n);

    //! calculate gradient of function at integration points
    Vector3d gradient(FESolidElement& el, std::vector<double>& fn, int n);
    
    //! calculate spatial gradient of std::vector function at integration points
    Matrix3d gradient(FESolidElement& el, Vector3d* fn, int n);
    
    //! calculate gradient of function at integration points at intermediate time
    Vector3d gradient(FESolidElement& el, std::vector<double>& fn, int n, const double alpha);
    
    //! calculate spatial gradient of std::vector function at integration points at intermediate time
    Matrix3d gradient(FESolidElement& el, Vector3d* fn, int n, const double alpha);
    
    //! calculate spatial gradient of std::vector function at integration points
    //! at previous time
    Matrix3d gradientp(FESolidElement& el, Vector3d* fn, int n);
    
    //! calculate spatial gradient of std::vector function at integration points
    tens3dls gradient(FESolidElement& el, Matrix3ds* fn, int n);
    
    //! calculate spatial gradient of std::vector function at integration points
    //! at previous time
    tens3dls gradientp(FESolidElement& el, Matrix3ds* fn, int n);
    
    //! calculate material gradient of scalar function at integration points
    Vector3d Gradient(FESolidElement& el, double* fn, int n);
    
    //! calculate material gradient of scalar function at integration points
    Vector3d Gradient(FESolidElement& el, std::vector<double>& fn, int n);
    
    //! calculate material gradient of std::vector function at integration points
    Matrix3d Gradient(FESolidElement& el, Vector3d* fn, int n);
    
    //! calculate material gradient of std::vector function at integration points
    tens3dls Gradient(FESolidElement& el, Matrix3ds* fn, int n);
    
    //! calculate jacobian in reference frame
    double detJ0(FESolidElement& el, int n);
    
    //! calculate jacobian in current frame
    double detJt(FESolidElement& el, int n);
    
    //! calculate jacobian in current frame
    double detJt(FESolidElement& el, int n, const double alpha);
    
    //! calculates covariant basis vectors in reference configuration at an integration point
    void CoBaseVectors0(FESolidElement& el, int j, Vector3d g[3]);
    
    //! calculates covariant basis vectors at an integration point
    void CoBaseVectors(FESolidElement& el, int j, Vector3d g[3]);
    
    //! calculates covariant basis vectors at an integration point
    void CoBaseVectors(FESolidElement& el, int j, Vector3d g[3], const double alpha);
    
    //! calculates contravariant basis vectors in reference configuration at an integration point
    void ContraBaseVectors0(FESolidElement& el, int j, Vector3d g[3]);
    
    //! calculates contravariant basis vectors at an integration point
    void ContraBaseVectors(FESolidElement& el, int j, Vector3d g[3]);
    
    //! calculates contravariant basis vectors at an integration point
    void ContraBaseVectors(FESolidElement& el, int j, Vector3d g[3], const double alpha);
    
    //! calculates parametric derivatives of covariant basis vectors at an integration point
    void CoBaseVectorDerivatives(FESolidElement& el, int j, Vector3d dg[3][3]);
    
    //! calculates parametric derivatives of covariant basis vectors at an integration point
    void CoBaseVectorDerivatives0(FESolidElement& el, int j, Vector3d dg[3][3]);
    
    //! calculates parametric derivatives of covariant basis vectors at an integration point
    void CoBaseVectorDerivatives(FESolidElement& el, int j, Vector3d dg[3][3], const double alpha);
    
    //! calculates parametric derivatives of contravariant basis vectors at an integration point
    void ContraBaseVectorDerivatives(FESolidElement& el, int j, Vector3d dg[3][3]);
    
    //! calculates parametric derivatives of contravariant basis vectors at an integration point
    void ContraBaseVectorDerivatives0(FESolidElement& el, int j, Vector3d dg[3][3]);
    
    //! calculates parametric derivatives of contravariant basis vectors at an integration point
    void ContraBaseVectorDerivatives(FESolidElement& el, int j, Vector3d dg[3][3], const double alpha);
    
    //! calculate the laplacian of a std::vector function at an integration point
    Vector3d lapvec(FESolidElement& el, Vector3d* fn, int n);
    
    //! calculate the gradient of the divergence of a std::vector function at an integration point
    Vector3d gradivec(FESolidElement& el, Vector3d* fn, int n);
    
    //! calculate the transpose of the gradient of the shape function gradients at an integration point
    void gradTgradShape(FESolidElement& el, int j, std::vector<Matrix3d>& mn);
    
    //! calculate spatial gradient of shapefunctions at integration point (returns Jacobian determinant)
    double ShapeGradient(FESolidElement& el, int n, Vector3d* GradH);
    
    //! calculate spatial gradient of shapefunctions at integration point (returns Jacobian determinant)
    double ShapeGradient(FESolidElement& el, int n, Vector3d* GradH, const double alpha);
    
    //! calculate spatial gradient of shapefunctions at integration point in reference frame (returns Jacobian determinant)
    double ShapeGradient0(FESolidElement& el, int n, Vector3d* GradH);
    
    //! calculate spatial gradient of shapefunctions at integration point (returns Jacobian determinant)
    double ShapeGradient(FESolidElement& el, double r, double s, double t, Vector3d* GradH);
    
    //! calculate spatial gradient of shapefunctions at integration point in reference frame (returns Jacobian determinant)
    double ShapeGradient0(FESolidElement& el, double r, double s, double t, Vector3d* GradH);

	//! calculate the volume of an element in reference frame
	double Volume(FESolidElement& el);

	//! calculate the volume of an element in current frame
	double CurrentVolume(FESolidElement& el);

public:
	//! get the current nodal coordinates
	void GetCurrentNodalCoordinates(const FESolidElement& el, Vector3d* rt);
	void GetCurrentNodalCoordinates(const FESolidElement& el, Vector3d* rt, double alpha);

	//! get the reference nodal coordinates
	void GetReferenceNodalCoordinates(const FESolidElement& el, Vector3d* r0);

	//! get the nodal coordinates at previous state
	void GetPreviousNodalCoordinates(const FESolidElement& el, Vector3d* rp);

public:
	//! loop over elements
	void ForEachSolidElement(std::function<void(FESolidElement& el)> f);

	//! return the degrees of freedom of an element for this domain
	virtual int GetElementDofs(FESolidElement& el);

public:
	// Evaluate an integral over the domain and assemble into global load std::vector
	virtual void LoadVector(
		FEGlobalVector& R,			// the global std::vector to assembe the load std::vector in
		const FEDofList& dofList,	// the degree of freedom list
		FEVolumeVectorIntegrand f	// the actual integrand function
	);

	//! Evaluate the stiffness matrix of a load
	virtual void LoadStiffness(
		FELinearSystem& LS,			// The solver does the assembling
		const FEDofList& dofList_a,	// The degree of freedom list of node a
		const FEDofList& dofList_b,	// The degree of freedom list of node b
		FEVolumeMatrixIntegrand f	// the matrix function to evaluate
	);

protected:
    std::vector<FESolidElement>	m_Elem;		//!< array of elements
	FE_Element_Spec			m_elemSpec;	//!< the element spec

	FEDofList	m_dofU;
	FEDofList	m_dofSU;
};
