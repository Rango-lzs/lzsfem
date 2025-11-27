#pragma once
#include "elements/RgElement/RgElement.h"
#include <vector>

//! This class defines a solid element
class RgSolidElementTraits;

//! Corresponds to Abaqus Continuum(Solid) elements, which are multidimensional (3D)
class FEM_EXPORT RgSolidElement : public RgElement
{
public:
	//! default constructor
	RgSolidElement() = default;
	~RgSolidElement() = default;

	RgSolidElement(const RgSolidElement& el);
	RgSolidElement& operator = (const RgSolidElement& el);

	virtual int dim();
};
