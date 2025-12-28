#pragma once
#include "elements/RgElement/RgElement.h"

//-----------------------------------------------------------------------------
//!  This class defines a surface element

class FEM_EXPORT RgSurfaceElement : public RgElement
{
public:
	RgSurfaceElement();

	RgSurfaceElement(const RgSurfaceElement& el);

	RgSurfaceElement& operator = (const RgSurfaceElement& el);

public:
    //! local ID of surface element
	int		m_lid;

	// indices of solid or shell element this surface is a face of
	// For solids, a surface element can be connected to two elements 
	// if the surface is an inside surface. For boundary surfaces
	// the second element index is -1. 
	RgElement*		m_elem[2];
};

