#pragma once
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
    RgElementShape(ElementShapeType shape, int nodes);

	//! values of shape functions
	virtual std::vector<double> evalH(NaturalCoord coord) = 0;

	//! values of shape function derivatives
    virtual  std::vector<std::vector<double>> evalDeriv(NaturalCoord coord) = 0;

	//! values of shape function second derivatives
    virtual std::vector<std::vector<double>> evalDeriv2(NaturalCoord coord) = 0;

private:
    ElementShapeType  mShpType;
    int mNodes;
};
