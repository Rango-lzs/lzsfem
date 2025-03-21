#pragma once
#include "elements/RgElemTypeDefine.h"
#include <vector>

// ��Ȼ���������
struct NaturalCoord
{
    double r;
    double s;
    double t;
};

// Base class for defining element shape classes for (3D) solid elements
// ֻ������״�������Լ���״��������Ȼ����ĵ���
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
