
#pragma once
#include <vector>
#include "femcore/fecore_enum.h"
#include "femcore/fem_export.h"

class FEElement;
class FEElementTraits;
class FEElementShape;

//-----------------------------------------------------------------------------
//! This class stores the different element traits classes
//! The purpose of this class is to store all the different element traits classes.
//! It are these traits classes that define the different element types. All different
//! traits must be registered before they can be assigned to elements.

class FEM_EXPORT FEElementLibrary
{
public:
	//! destructor
	~FEElementLibrary();

	//! return the element library
	static FEElementLibrary* GetInstance();

	//! Assign a traits class to an element
	static void SetElementTraits(FEElement& el, int id);

	//! return element traits data
	static FEElementTraits* GetElementTraits(int ntype);

	//! return element shape class
	static FEElementShape* GetElementShapeClass(FE_Element_Shape eshape);

	//! return the element shape of a given element type
	static FE_Element_Shape GetElementShape(int ntype);

	//! return the element class of a given element type
	static FE_Element_Category GetElementClass(int ntype);

	//! initialize library
	static void Initialize();

private:
	//! constructor
	FEElementLibrary(){}
	FEElementLibrary(const FEElementLibrary&){}

	//! Function to register an element shape class
	int RegisterShape(FEElementShape* pshape);

	//! Function to register an element traits class
	int RegisterTraits(FEElementTraits* ptrait);

private:
	std::vector<FEElementTraits*>	m_Traits;	//!< pointer to registered element traits
	std::vector<FEElementShape*>	m_Shape;	//!< pointer to registered element shapes
	static FEElementLibrary*	m_pThis;
};
