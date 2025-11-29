#pragma once
// RgSolid3dElement.h
// Derived 3D solid element for RgFem
//
// Minimal interface derived from FESolidElement.
// Add/adjust methods to match your project's FE core.

#include "RgSolidElement.h"
#include <vector>

class FEElementMatrix;
class FEElementVector;
class FEMaterialPoint;
class DumpStream;
struct vec3d;

//定义了三维单元的行为，坐标维数为3
class RgSolid3dElement : public RgSolidElement
{
public:
	// default constructor / destructor
	RgSolid3dElement() = default;
	virtual ~RgSolid3dElement() = default;

	// copy
	RgSolid3dElement(const RgSolid3dElement& el) = default;
	RgSolid3dElement& operator=(const RgSolid3dElement& el) = default;

	// return spatial dimension (3 for 3D solids)
	int dim() override;

    // 高斯点的形函数相关计算
    // n : the n-th gauss point
    std::vector<double> gaussPoint(int n);
    
    // returns weights of integration points as a vector
    std::vector<double> GaussWeights() const;

    // dH/dr[n] -- return shape function derivative arrays as vectors
    // n : the n-th gauss point
    // return : N shape function derivative
    std::vector<double> Gr(int n) const;

    // shape function derivative to r
    std::vector<double> Gs(int n) const;
    
    // shape function derivative to s
    std::vector<double> Gt(int n) const;

    // dH2/dr2[n] -- second derivatives as vectors
    std::vector<double> Grr(int n) const;

    std::vector<double> Gsr(int n) const;

    std::vector<double> Gtr(int n) const;
  
    std::vector<double> Grs(int n) const;

    std::vector<double> Gss(int n) const;

    std::vector<double> Gts(int n) const;

    std::vector<double> Grt(int n) const;

    std::vector<double> Gst(int n) const;

    std::vector<double> Gtt(int n) const;

    //! values of shape functions (unchanged API)  这些接口计算任意点的形函数
    void shape_fnc(double* H, double r, double s, double t) const;

    //! values of shape function derivatives (unchanged API)
    void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t) const;

    //! values of shape function second derivatives (unchanged API)
    void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s,
                      double t) const;
};