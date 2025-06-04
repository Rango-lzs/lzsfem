#pragma once
#include "FEValuator.h"
#include "FEDataMap.h"
//#include "MathObject.h"

//---------------------------------------------------------------------------------------
// Base class for evaluating Vector3d parameters
class FEM_EXPORT FEVec3dValuator : public FEValuator
{	
	DECLARE_META_CLASS(FEVec3dValuator, FEValuator);

public:
	FEVec3dValuator();

public:
	// evaluate value at a material point
	virtual Vector3d operator()(const FEMaterialPoint& pt) = 0;

	// create a copy of the valuator
	virtual FEVec3dValuator* copy() = 0;

	// is the valuator constant
	virtual bool isConst() { return false; }

	// return the const value
	virtual Vector3d* constValue() { return nullptr; }

	// return a unit vector
	Vector3d unitVector(const FEMaterialPoint& pt)
	{
		Vector3d v = operator () (pt);
		return v.Normalized();
	}
};
