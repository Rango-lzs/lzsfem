#pragma once

#include "elements/RgElement/RgSolid3dElement.h"

//!  This class defines a 2D element
class FEM_EXPORT RgSolid2dElement : public RgSolidElement
{
public:
    //! default constructor
    RgSolid2dElement() = default;

    //! copy constructor
    RgSolid2dElement(const RgSolid2dElement& el);

    //! assignment operator
    RgSolid2dElement& operator=(const RgSolid2dElement& el);

    virtual int dim() override;

};