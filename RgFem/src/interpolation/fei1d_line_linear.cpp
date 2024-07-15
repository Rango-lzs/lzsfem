/*****************************************************************//**
 * \file   fei1d_line_linear.cpp
 * \brief
 *
 * \author Leizs
 * \date   November 2023
 *********************************************************************/
#include "fei1d_line_linear.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "floatmatrixf.h"
#include "floatarrayf.h"
#include "fei_element_geometry_wrapper.h"

namespace fem
{
	double FEI1dLinear::giveLength() const
	{
		const std::unique_ptr<FEIElementGeometry>& cellgeo = giveElemGeomety();
		return fabs(cellgeo->giveVertexCoordinates(2).at(cindx) - cellgeo->giveVertexCoordinates(1).at(cindx));
	}

	FloatArrayF<2> FEI1dLinear::evalN(double ksi)
	{
		return { (1. - ksi) * 0.5, (1. + ksi) * 0.5 };
	}

	void FEI1dLinear::evalN(FloatArray& answer, const FloatArray& lcoords) const
	{
		double ksi = lcoords.at(1);
		answer.resize(2);

		answer.at(1) = (1. - ksi) * 0.5;
		answer.at(2) = (1. + ksi) * 0.5;
	}

	std::pair<double, FloatMatrixF<1, 2>> FEI1dLinear::evaldNdx() const
	{
		const std::unique_ptr<FEIElementGeometry>& cellgeo = giveElemGeomety();
		double l = cellgeo->giveVertexCoordinates(2).at(cindx) - cellgeo->giveVertexCoordinates(1).at(cindx);
		return { 0.5 * l, {-1.0 / l, 1.0 / l} };
	}

	double FEI1dLinear::evaldNdx(FloatMatrix& answer, const FloatArray& lcoords) const
	{
		double l = giveElemGeomety()->giveVertexCoordinates(2).at(cindx) - giveElemGeomety()->giveVertexCoordinates(1).at(cindx);
		answer.resize(2, 1);

		answer.at(1, 1) = -1.0 / l;
		answer.at(2, 1) = 1.0 / l;
		return 0.5 * l;
	}

	//1 d has only one coordinate 'x'
	void FEI1dLinear::local2global(FloatArray& answer, const FloatArray& lcoords) const
	{
		FloatArray n;
		answer.resize(1);

		this->evalN(n, lcoords);
		const std::unique_ptr<FEIElementGeometry>& cellgeo = giveElemGeomety();
		answer.at(1) = n.at(1) * cellgeo->giveVertexCoordinates(1).at(cindx) +
			n.at(2) * cellgeo->giveVertexCoordinates(2).at(cindx);
	}

	int FEI1dLinear::global2local(FloatArray& answer, const FloatArray& coords) const
	{
		const std::unique_ptr<FEIElementGeometry>& cellgeo = giveElemGeomety();
		double x1 = cellgeo->giveVertexCoordinates(1).at(cindx);
		double x2 = cellgeo->giveVertexCoordinates(2).at(cindx);
		double ksi = (2.0 * coords.at(1) - (x1 + x2)) / (x2 - x1);
		answer.resize(1);
		//answer.at(1) = clamp(ksi, -1., 1.);
		return fabs(ksi) <= 1.0;
	}

	double FEI1dLinear::giveTransformationJacobian(const FloatArray& lcoords) const
	{
		const std::unique_ptr<FEIElementGeometry>& cellgeo = giveElemGeomety();
		//dn/dksi = 1/2(x2-x1)
		return 0.5 * (cellgeo->giveVertexCoordinates(2).at(cindx) - cellgeo->giveVertexCoordinates(1).at(cindx));
	}
} // end namespace fem
