#pragma once
#include "elements/ElementShape/RgElementShape.h"
#include "../RgElemTypeDefine.h"
#include <vector>

struct NaturalCoord;

//=============================================================================
// Base class for defining element shape classes for (3D) solid elements
// 只计算形状函数、以及形状函数对自然坐标的导数
class RgSolidElementShape : public RgElementShape
{
public:
	RgSolidElementShape(ElementShape shape, int nodes) : RgElementShape(shape, nodes) {}

	// values of shape functions with size N
    virtual std::vector<double> evalH(const NaturalCoord& coord) {};

    // values of shape function derivatives with size 3,N (2,N for 2d)
    virtual std::vector<std::vector<double>> evalDeriv(const NaturalCoord& coord)  {};

    // values of shape function second derivatives with size 6,N (3,N for 2d)
    virtual std::vector<std::vector<double>> evalDeriv2(const NaturalCoord& coord) {};
};
