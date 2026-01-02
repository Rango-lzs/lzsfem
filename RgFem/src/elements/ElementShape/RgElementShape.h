#pragma once
#include "elements/RgElemTypeDefine.h"
#include <vector>

namespace RgFem{
    class NaturalCoord;
}

// Base class for defining element shape classes for (3D) solid elements
class RgElementShape 
{
public:
    RgElementShape(ElementShape shape, int nodes);

	//values of shape functions with size N
	virtual std::vector<double> evalH(const RgFem::NaturalCoord& coord) = 0;

	//values of shape function derivatives with size 3,N (2,N for 2d)
    virtual  std::vector<std::vector<double>> evalDeriv(const RgFem::NaturalCoord& coord) = 0;

	//values of shape function second derivatives with size 6,N (3,N for 2d)
    virtual std::vector<std::vector<double>> evalDeriv2(const RgFem::NaturalCoord& coord) = 0;

    ElementShape shapeType() const
    {
        return mShpType;
    }

    int nodes() const
    {
        return mNodes;
    }

private:
    ElementShape mShpType;
    int mNodes;
};