/*****************************************************************//**
 * \file   fei2d_line_linear.cpp
 * \brief
 *
 * \author Leizs
 * \date   November 2023
 *********************************************************************/

#include "fei2d_line_linear.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "gaussintegrationrule.h"
#include "fei_element_geometry_wrapper.h"

namespace fem
{

	FloatArrayF<2> FEI2dLineLin::evalN(double xi)
	{
		return { (1. - xi) * 0.5, (1. + xi) * 0.5 };
	}

	void FEI2dLineLin::evalN(FloatArray& answer, const FloatArray& lcoords) const
	{
		double xi = lcoords(0);
		answer.resize(2);
		answer.at(1) = (1. - xi) * 0.5;
		answer.at(2) = (1. + xi) * 0.5;
	}

	double FEI2dLineLin::evaldNdx(FloatMatrix& answer, const FloatArray& lcoords) const
	{
		// Not meaningful to return anything.
		answer.clear();
		return 0.;
	}

	void FEI2dLineLin::evaldNdxi(FloatMatrix& answer, const FloatArray& lcoords) const
	{
		answer.resize(2, 1);
		answer(0, 0) = -0.5;
		answer(1, 0) = 0.5;
	}

	void FEI2dLineLin::local2global(FloatArray& answer, const FloatArray& lcoords) const
	{
		FloatArray n;
		this->evalN(n, lcoords);
		answer.resize(max(xind, yind));
		answer.zero();
		const std::unique_ptr<FEIElementGeometry>& cellgeo = giveElemGeomety();
		answer.at(xind) = n(0) * cellgeo->giveVertexCoordinates(1).at(xind) +
			n(1) * cellgeo->giveVertexCoordinates(2).at(xind);
		answer.at(yind) = n(0) * cellgeo->giveVertexCoordinates(1).at(yind) +
			n(1) * cellgeo->giveVertexCoordinates(2).at(yind);
	}

	int FEI2dLineLin::global2local(FloatArray& answer, const FloatArray& gcoords) const
	{
		const std::unique_ptr<FEIElementGeometry>& cellgeo = giveElemGeomety();
		double x2_x1 = cellgeo->giveVertexCoordinates(2).at(xind) - cellgeo->giveVertexCoordinates(1).at(xind);
		double y2_y1 = cellgeo->giveVertexCoordinates(2).at(yind) - cellgeo->giveVertexCoordinates(1).at(yind);

		// Projection of the global coordinate gives the value interpolated in [0,1].
		double dx = gcoords(0) - cellgeo->giveVertexCoordinates(1).at(xind);
		double dy = gcoords(1) - cellgeo->giveVertexCoordinates(1).at(yind);
		double xi = (x2_x1) ? (sqrt(dx * dx * (1 + (y2_y1 / x2_x1) * (y2_y1 / x2_x1)))) / (sqrt(x2_x1 * x2_x1 + y2_y1 * y2_y1)) : sqrt(dy * dy) / (sqrt(x2_x1 * x2_x1 + y2_y1 * y2_y1));
		// Map to [-1,1] domain.
		xi = xi * 2 - 1;

		answer.resize(1);
		answer(0) = clamp(xi, -1., 1.);
		return false;
	}

	double FEI2dLineLin::giveTransformationJacobian(const FloatArray& lcoords) const
	{
		const std::unique_ptr<FEIElementGeometry>& cellgeo = giveElemGeomety();
		double x2_x1 = cellgeo->giveVertexCoordinates(2).at(xind) - cellgeo->giveVertexCoordinates(1).at(xind);
		double y2_y1 = cellgeo->giveVertexCoordinates(2).at(yind) - cellgeo->giveVertexCoordinates(1).at(yind);
		return sqrt(x2_x1 * x2_x1 + y2_y1 * y2_y1) / 2.0;
	}
} // end namespace oofem
