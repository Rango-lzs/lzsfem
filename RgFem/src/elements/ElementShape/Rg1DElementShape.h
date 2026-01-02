#pragma once
#include "elements/ElementShape/RgElementShape.h"
#include "elements/NaturalCoord.h"
#include "../RgElemTypeDefine.h"
#include <vector>

namespace RgFem
{
    class NaturalCoord;
}

//=============================================================================
// Base class for defining element shape classes for (3D) solid elements
namespace RgFem {

class FEM_EXPORT Rg1DElementShape : public RgElementShape
{
public:
	Rg1DElementShape(ElementShape shape, int nodes) : RgElementShape(shape, nodes) {}

	//values of shape functions with size N
	virtual std::vector<double> evalH(const RgFem::NaturalCoord& coord) override;

	//values of shape function derivatives with size 2,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) override;

	//values of shape function second derivatives with size 3,N (for 2D surface elements)
    virtual std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) override;
};

} // namespace RgFem