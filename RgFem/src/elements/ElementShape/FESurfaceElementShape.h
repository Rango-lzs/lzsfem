#pragma once
#include "FEElementShape.h"

//=============================================================================
// Base class for defining element shape classes for (3D) solid elements
class FESurfaceElementShape : public FEElementShape
{
public:
	FESurfaceElementShape(ElementShape shape, int nodes) : FEElementShape(shape, nodes) {}

	//! values of shape functions
	virtual void shape_fnc(double* H, double r, double s) = 0;

	//! values of shape function derivatives
	virtual void shape_deriv(double* Hr, double* Hs, double r, double s) = 0;

	//! values of shape function second derivatives
	virtual void shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s) = 0;
};

//=============================================================================
// Class for QUAD4 elements
class FEQuad4 : public FESurfaceElementShape
{
public:
	FEQuad4() : FESurfaceElementShape(ET_QUAD4, 4) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s) override;

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) override;

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s) override;
};

//=============================================================================
// Class for QUAD8 elements
class FEQuad8 : public FESurfaceElementShape
{
public:
	FEQuad8() : FESurfaceElementShape(ET_QUAD8, 8) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s) override;

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) override;

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s) override;
};

//=============================================================================
// Class for QUAD9 elements
class FEQuad9 : public FESurfaceElementShape
{
public:
	FEQuad9() : FESurfaceElementShape(ET_QUAD9, 9) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s) override;

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) override;

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s) override;
};

//=============================================================================
// Class for TRI3 elements
class FETri3 : public FESurfaceElementShape
{
public:
	FETri3() : FESurfaceElementShape(ET_TRI3, 4) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s) override;

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) override;

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s) override;
};

//=============================================================================
// Class for TRI6 elements
class FETri6 : public FESurfaceElementShape
{
public:
	FETri6() : FESurfaceElementShape(ET_TRI6, 6) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s) override;

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) override;

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s) override;
};

//=============================================================================
// Class for TRI7 elements
class FETri7 : public FESurfaceElementShape
{
public:
	FETri7() : FESurfaceElementShape(ET_TRI7, 7) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s) override;

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) override;

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s) override;
};

//=============================================================================
// Class for TRI10 elements
class FETri10 : public FESurfaceElementShape
{
public:
	FETri10() : FESurfaceElementShape(ET_TRI10, 10) {}

	//! values of shape functions
	void shape_fnc(double* H, double r, double s) override;

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) override;

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s) override;
};
