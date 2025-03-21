#pragma once
#include "elements/RgElemTypeDefine.h"
#include <vector>

// 自然坐标点类型
struct NaturalCoord
{
    double r;
    double s;
    double t;
};

// Base class for defining element shape classes for (3D) solid elements
// 只计算形状函数、以及形状函数对自然坐标的导数
class RgElementShape 
{
public:
    RgElementShape(ElementShape shape, int nodes);

	//values of shape functions with size N
	virtual std::vector<double> evalH(NaturalCoord coord) = 0;

	//values of shape function derivatives with size 3,N (2,N for 2d)
    virtual  std::vector<std::vector<double>> evalDeriv(NaturalCoord coord) = 0;

	//values of shape function second derivatives with size 6,N (3,N for 2d)
    virtual std::vector<std::vector<double>> evalDeriv2(NaturalCoord coord) = 0;

private:
    ElementShape mShpType;
    int mNodes;
};
