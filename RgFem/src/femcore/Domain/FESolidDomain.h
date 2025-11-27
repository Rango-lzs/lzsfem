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
    FEElement& ElementRef(int n) override { return *m_Elem[n]; }
	const FEElement& ElementRef(int n) const override { return *m_Elem[n]; }

    ElementType GetElementType() const { return m_Elem[0]->elementType(); }
    
    int GetElementShape() const { return m_Elem[0]->ShapeFunctions(-1); }

	FE_Element_Spec GetElementSpec() const;
    
    //! find the element in which point y lies
    FESolidElement* FindElement(const Vector3d& y, double r[3]);

	//! find the element in which point y lies (reference configuration)
	FESolidElement* FindReferenceElement(const Vector3d& y, double r[3]);


public:
	// Evaluate an integral over the domain and assemble into global load std::vector
	virtual void LoadVector(
		FEGlobalVector& R,			// the global std::vector to assembe the load std::vector in
		const FEDofList& dofList,	// the degree of freedom list
		FEVolumeVectorIntegrand f	// the actual integrand function
	);

	//! Evaluate the stiffness Matrix of a load
	virtual void LoadStiffness(
		FELinearSystem& LS,			// The solver does the assembling
		const FEDofList& dofList_a,	// The degree of freedom list of node a
		const FEDofList& dofList_b,	// The degree of freedom list of node b
		FEVolumeMatrixIntegrand f	// the Matrix function to evaluate
	);

protected:
    std::vector<FESolidElement*>	m_Elem;		//!< array of elements
	FE_Element_Spec			m_elemSpec;	//!< the element spec

	FEDofList	m_dofU;
	FEDofList	m_dofSU;
};
