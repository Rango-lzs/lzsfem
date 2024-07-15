/*****************************************************************//**
 * \file   fei1d_line_quad.cpp
 * \brief
 *
 * \author Leizs
 * \date   November 2023
 *********************************************************************/

#include "fei1d_line_quad.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "floatmatrixf.h"
#include "floatarrayf.h"
#include "fei_element_geometry_wrapper.h"

namespace fem
{
	FloatArrayF<3> FEI1dQuad::evalN(double ksi)
	{
		// 1.....3.....2
		return { ksi * (ksi - 1.) * 0.5, ksi * (1. + ksi) * 0.5, (1. - ksi * ksi) };
	}

	void FEI1dQuad::evalN(FloatArray& answer, const FloatArray& lcoords) const
	{
		double ksi = lcoords.at(1);
		answer.resize(3);
		answer.zero();

		answer.at(1) = ksi * (ksi - 1.) * 0.5;
		answer.at(2) = ksi * (1. + ksi) * 0.5;
		answer.at(3) = (1. - ksi * ksi);
	}

	std::pair<double, FloatMatrixF<1, 3>>
		FEI1dQuad::evaldNdx(double ksi) const
	{
		const std::unique_ptr<FEIElementGeometry>& cellgeo = giveElemGeomety();
		double x1 = cellgeo->giveVertexCoordinates(1).at(cindx);
		double x2 = cellgeo->giveVertexCoordinates(2).at(cindx);
		double x3 = cellgeo->giveVertexCoordinates(3).at(cindx);

		double J = 1. / 2. * (2 * ksi - 1) * x1 + 1. / 2. * (2 * ksi + 1) * x2 - 2. * ksi * x3;

		FloatMatrixF<1, 3> ans = {
			(-1. / 2. + ksi) / J,
			(1. / 2. + ksi) / J,
			-2. * ksi / J,
		};
		return { J, ans };
	}


	double
		FEI1dQuad::evaldNdx(FloatMatrix& answer, const FloatArray& lcoords) const
	{
		const std::unique_ptr<FEIElementGeometry>& cellgeo = giveElemGeomety();
		double J = this->giveTransformationJacobian(lcoords);
		double ksi = lcoords.at(1);
		answer.resize(1, 3);
		answer.zero();

		answer.at(1, 1) = (-1. / 2. + ksi) / J;
		answer.at(1, 2) = (1. / 2. + ksi) / J;
		answer.at(1, 3) = -2. * ksi / J;
		return J;
	}

	void
		FEI1dQuad::local2global(FloatArray& answer, const FloatArray& lcoords) const
	{
		FloatArray n;
		answer.resize(1);

		const std::unique_ptr<FEIElementGeometry>& cellgeo = giveElemGeomety();
		this->evalN(n, lcoords);
		answer.at(1) = n.at(1) * cellgeo->giveVertexCoordinates(1).at(cindx) +
			n.at(2) * cellgeo->giveVertexCoordinates(2).at(cindx) +
			n.at(3) * cellgeo->giveVertexCoordinates(3).at(cindx);
	}

	int
		FEI1dQuad::global2local(FloatArray& answer, const FloatArray& coords) const
	{
		const std::unique_ptr<FEIElementGeometry>& cellgeo = giveElemGeomety();
		double x1 = cellgeo->giveVertexCoordinates(1).at(cindx);
		double x2 = cellgeo->giveVertexCoordinates(2).at(cindx);
		double x3 = cellgeo->giveVertexCoordinates(3).at(cindx);

		double a = 0.5 * (x1 + x2) - x3;
		double b = 0.5 * (x2 - x1);
		double c = x3 - coords.at(1);

		answer.resize(1);
		if (fabs(a) < 1.e-6) {
			double ksi = (2.0 * coords.at(1) - (x1 + x2)) / (x2 - x1);
			answer.at(1) = clamp(ksi, -1., 1.);
			return fabs(ksi) <= 1.0;
		}
		else {
			double ksi1 = (-b + sqrt(b * b - 4. * a * c)) / (2. * a);
			double ksi2 = (-b - sqrt(b * b - 4. * a * c)) / (2. * a);

			if ((fabs(ksi1) <= 1.) && (fabs(ksi2) <= 1.)) { // Two roots, element must be bad
				answer.at(1) = 0.;
				return 0;
			}
			else if (fabs(ksi1) <= 1.) {
				answer.at(1) = ksi1;
				return 1;
			}
			else if (fabs(ksi2) <= 1.) {
				answer.at(1) = ksi2;
				return 1;
			}
			else {
				answer.at(1) = 0.;
				return 0;
			}
		}
	}

	double
		FEI1dQuad::giveTransformationJacobian(const FloatArray& lcoords) const
	{
		const std::unique_ptr<FEIElementGeometry>& cellgeo = giveElemGeomety();
		double x1 = cellgeo->giveVertexCoordinates(1).at(cindx);
		double x2 = cellgeo->giveVertexCoordinates(2).at(cindx);
		double x3 = cellgeo->giveVertexCoordinates(3).at(cindx);
		double ksi = lcoords.at(1);

		double J = 1. / 2. * (2 * ksi - 1) * x1 + 1. / 2. * (2 * ksi + 1) * x2 - 2. * ksi * x3;
		return J;
	}

	double FEI1dQuad::giveLength() const
	{
		const std::unique_ptr<FEIElementGeometry>& cellgeo = giveElemGeomety();
		return fabs(cellgeo->giveVertexCoordinates(2).at(cindx) - cellgeo->giveVertexCoordinates(1).at(cindx));
	}

} // end namespace fem
