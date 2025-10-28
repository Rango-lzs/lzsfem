#pragma once
#include "elements/FEElementShape.h"
#include "../RgElemTypeDefine.h"
#include <vector>

struct NaturalCoord;

//=============================================================================
// Base class for defining element shape classes for (3D) solid elements
// ֻ״Լ״Ȼĵ
class FESolidElementShape : public FEElementShape
{
public:
	FESolidElementShape(ElementShape shape, int nodes) : FEElementShape(shape, nodes) {}

	//! values of shape functions, size N(number of shape function)
	virtual void shape_fnc(double* H, double r, double s, double t) = 0;
    virtual std::vector<double> shape_fnc(const NaturalCoord& coor) = 0;

	//! values of shape function derivatives,size 3 x N
	virtual void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t) = 0;
    virtual std::vector<std::vector<double>> shape_deriv(const NaturalCoord& coor) = 0;
	 
	//! values of shape function second derivatives,size 6 x N
	virtual void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t) = 0;
    virtual std::vector<std::vector<double>> shape_deriv2(const NaturalCoord& coor) = 0;
};

//=============================================================================
class FETet4 : public FESolidElementShape
{
public:
	FETet4() : FESolidElementShape(ET_TET4, 4) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};

//=============================================================================
class FETet10 : public FESolidElementShape
{
public:
    FETet10()
        : FESolidElementShape(ET_TET10, 10)
    {
    }

    //! values of shape functions
    void shape_fnc(double* H, double r, double s, double t);

    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

    //! values of shape function second derivatives
    void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s,
                      double t);
};

//=============================================================================
class FEHex8 : public FESolidElementShape
{
public:
	FEHex8() : FESolidElementShape(ET_HEX8, 8) {}

	//! values of shape functions
	//! H<8>
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	//! H<3,8>
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	//! H<6,8>
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};

//=============================================================================
class FEHex20 : public FESolidElementShape
{
public:
	FEHex20() : FESolidElementShape(ET_HEX20, 20) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};


//=============================================================================
//! Base class for 27-node quadratic hexahedral element
class FEHex27 : public FESolidElementShape
{
public:
	FEHex27() : FESolidElementShape(ET_HEX27, 27) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};

//=============================================================================
class FEPenta6 : public FESolidElementShape
{
public:
    FEPenta6()
        : FESolidElementShape(ET_PENTA6, 6)
    {
    }

    //! values of shape functions
    void shape_fnc(double* H, double r, double s, double t);

    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

    //! values of shape function second derivatives
    void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s,
                      double t);
};

//=============================================================================
class FEPenta15 : public FESolidElementShape
{
public:
    FEPenta15()
        : FESolidElementShape(ET_PENTA15, 15)
    {
    }

    //! values of shape functions
    void shape_fnc(double* H, double r, double s, double t);

    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

    //! values of shape function second derivatives
    void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s,
                      double t);
};

//=============================================================================
class FEPyra5 : public FESolidElementShape
{
public:
	FEPyra5() : FESolidElementShape(ET_PYRA5, 5) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t);

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};


//=============================================================================
class FEPyra13 : public FESolidElementShape
{
public:
    FEPyra13() : FESolidElementShape(ET_PYRA13, 13) {}
    
    //! values of shape functions
    void shape_fnc(double* H, double r, double s, double t);
    
    //! values of shape function derivatives
    void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t);
    
    //! values of shape function second derivatives
    void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t);
};

