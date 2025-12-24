#pragma once
#include "elements/RgElemTypeDefine.h"
#include <vector>

//class NaturalCoord
//{
//public:
//    virtual ~NaturalCoord() = default;
//    NaturalCoord() = default;
//};
//
//
//class NaturalCoord3d : public NaturalCoord
//{
//public:
//    NaturalCoord3d() : r(0.0), s(0.0), t(0.0) {}
//    NaturalCoord3d(double r, double s, double t) : r(r), s(s), t(t) {}
//    
//    double getR() const { return r; }
//    double getS() const { return s; }
//    double getT() const { return t; }
//    
//    void setR(double r_val) { r = r_val; }
//    void setS(double s_val) { s = s_val; }
//    void setT(double t_val) { t = t_val; }
//
//private:
//    double r;
//    double s;
//    double t;
//};
//
//class NaturalCoord2d : public NaturalCoord
//{
//public:
//    NaturalCoord2d() : r(0.0), s(0.0) {}
//    NaturalCoord2d(double r, double s) : r(r), s(s) {}
//    
//    double getR() const { return r; }
//    double getS() const { return s; }
//    
//    void setR(double r_val) { r = r_val; }
//    void setS(double s_val) { s = s_val; }
//
//private:
//    double r;
//    double s;
//};

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