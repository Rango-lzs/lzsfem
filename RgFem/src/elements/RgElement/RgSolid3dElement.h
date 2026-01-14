#pragma once
// RgSolid3dElement.h
// Derived 3D solid element for RgFem
//
// Minimal interface derived from FESolidElement.
// Add/adjust methods to match your project's FE core.

#include "RgSolidElement.h"
#include <vector>

class FEElementMatrix;
class FEElementVector;
class DumpStream;
struct vec3d;


//定义了三维单元的行为，坐标维数为3
class FEM_EXPORT RgSolid3dElement : public RgSolidElement
{
public:
	// default constructor / destructor
	RgSolid3dElement() = default;
	virtual ~RgSolid3dElement() = default;

	// copy
	RgSolid3dElement(const RgSolid3dElement& el) = default;
	RgSolid3dElement& operator=(const RgSolid3dElement& el) = default;

	// return spatial dimension (3 for 3D solids)
	virtual int dim() override;                   
};