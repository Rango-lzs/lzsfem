#pragma once
#include "materials/RgMaterialPoint.h"
#include "femcore/FEBodyLoad.h"

//-----------------------------------------------------------------------------
//! This class is the base class for body forces
//! Derived classes need to implement the force and stiffness functions.
//
class FEM_EXPORT FEBodyForce : public FEBodyLoad
{
public:
	//! constructor
	FEBodyForce();

public:
	//! calculate the body force at a material point
	virtual Vector3d force(FEMaterialPoint& pt) = 0;

    //! calculate the divergence of the body force at a material point
    virtual double divforce(FEMaterialPoint& pt) { return (stiffness(pt)).tr(); }
    
	//! calculate constribution to stiffness Matrix
	virtual Matrix3ds stiffness(FEMaterialPoint& pt) = 0;

public:
	void LoadVector(FEGlobalVector& R) override;
	void StiffnessMatrix(FELinearSystem& LS) override;
};
