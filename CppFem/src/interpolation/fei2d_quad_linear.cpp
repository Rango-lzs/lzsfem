/*****************************************************************//**
 * \file   fei2d_quad_linear.cpp
 * \brief
 *
 * \author Leizs
 * \date   October 2023
 *********************************************************************/

#include "fei_element_geometry_wrapper.h"
#include "fei2d_quad_linear.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "floatmatrixf.h"
#include "floatarrayf.h"
#include "gaussintegrationrule.h"

namespace fem
{
	double FEI2dQuadLin::giveArea() const
	{
		FEIElementGeometry* cellgeo = giveElemGeomety();
		const auto& node1 = cellgeo->giveVertexCoordinates(1);
		const auto& node2 = cellgeo->giveVertexCoordinates(2);
		const auto& node3 = cellgeo->giveVertexCoordinates(3);
		const auto& node4 = cellgeo->giveVertexCoordinates(4);

		double x13 = node1.at(xind) - node3.at(xind);
		double y13 = node1.at(yind) - node3.at(yind);
		double x24 = node2.at(xind) - node4.at(xind);
		double y24 = node2.at(yind) - node4.at(yind);

		return fabs(0.5 * (x13 * y24 - x24 * y13));
	}

	FloatArrayF<4> FEI2dQuadLin::evalN(const FloatArrayF<2>& lcoords)
	{
		double ksi = lcoords[0];
		double eta = lcoords[1];

		return {
			(1. + ksi) * (1. + eta) * 0.25,
			(1. - ksi) * (1. + eta) * 0.25,
			(1. - ksi) * (1. - eta) * 0.25,
			(1. + ksi) * (1. - eta) * 0.25
		};
	}

	void FEI2dQuadLin::evalN(FloatArray& answer, const FloatArray& lcoords) const
	{
		answer = evalN({ lcoords[0], lcoords[1] });
	}

	std::pair<double, FloatMatrixF<2, 4>> FEI2dQuadLin::evaldNdx(const FloatArrayF<2>& lcoords) const
	{
		auto dndu = this->evaldNdxi(lcoords);
		FEIElementGeometry* cellgeo = giveElemGeomety();

		FloatMatrixF<2, 2> jacT;
		for (std::size_t i = 0; i < dndu.cols(); i++) {
			double x = cellgeo->giveVertexCoordinates(i + 1).at(xind);
			double y = cellgeo->giveVertexCoordinates(i + 1).at(yind);

			jacT(0, 0) += dndu(0, i) * x;
			jacT(0, 1) += dndu(0, i) * y;
			jacT(1, 0) += dndu(1, i) * x;
			jacT(1, 1) += dndu(1, i) * y;
		}
		return { det(jacT), dot(inv(jacT), dndu) };
	}

	double FEI2dQuadLin::evaldNdx(FloatMatrix& answer, const FloatArray& lcoords) const
	{
		auto tmp = evaldNdx(lcoords);
		answer = transpose(tmp.second);
		return tmp.first;
	}

	void FEI2dQuadLin::local2global(FloatArray& answer, const FloatArray& lcoords) const
	{
		double ksi = lcoords.at(1);
		double eta = lcoords.at(2);

		double n1 = (1. + ksi) * (1. + eta) * 0.25;
		double n2 = (1. - ksi) * (1. + eta) * 0.25;
		double n3 = (1. - ksi) * (1. - eta) * 0.25;
		double n4 = (1. + ksi) * (1. - eta) * 0.25;

		FEIElementGeometry* cellgeo = giveElemGeomety();
		const auto& p1 = cellgeo->giveVertexCoordinates(1);
		const auto& p2 = cellgeo->giveVertexCoordinates(2);
		const auto& p3 = cellgeo->giveVertexCoordinates(3);
		const auto& p4 = cellgeo->giveVertexCoordinates(4);

		answer = { n1 * p1.at(xind) + n2 * p2.at(xind) +
				  n3 * p3.at(xind) + n4 * p4.at(xind),
				  n1 * p1.at(yind) + n2 * p2.at(yind) +
				  n3 * p3.at(yind) + n4 * p4.at(yind) };
	}

#define POINT_TOL 1.e-6

	int FEI2dQuadLin::global2local(FloatArray& answer, const FloatArray& coords) const
	{
		double x1, x2, x3, x4, y1, y2, y3, y4, a1, a2, a3, a4, b1, b2, b3, b4;
		double a, b, c, ksi1, ksi2, ksi3, eta1 = 0.0, eta2 = 0.0, denom;
		int nroot;

		answer.resize(2);

		FEIElementGeometry* cellgeo = giveElemGeomety();
		x1 = cellgeo->giveVertexCoordinates(1).at(xind);
		x2 = cellgeo->giveVertexCoordinates(2).at(xind);
		x3 = cellgeo->giveVertexCoordinates(3).at(xind);
		x4 = cellgeo->giveVertexCoordinates(4).at(xind);
					
		y1 = cellgeo->giveVertexCoordinates(1).at(yind);
		y2 = cellgeo->giveVertexCoordinates(2).at(yind);
		y3 = cellgeo->giveVertexCoordinates(3).at(yind);
		y4 = cellgeo->giveVertexCoordinates(4).at(yind);

		a1 = x1 + x2 + x3 + x4;
		a2 = x1 - x2 - x3 + x4;
		a3 = x1 + x2 - x3 - x4;
		a4 = x1 - x2 + x3 - x4;

		b1 = y1 + y2 + y3 + y4;
		b2 = y1 - y2 - y3 + y4;
		b3 = y1 + y2 - y3 - y4;
		b4 = y1 - y2 + y3 - y4;

		a = a2 * b4 - b2 * a4;
		b = a1 * b4 + a2 * b3 - a3 * b2 - b1 * a4 - b4 * 4.0 * coords.at(xind) + a4 * 4.0 * coords.at(yind);
		c = a1 * b3 - a3 * b1 - 4.0 * coords.at(xind) * b3 + 4.0 * coords.at(yind) * a3;

		// solve quadratic equation for ksi
		//cubic(0.0, a, b, c, &ksi1, &ksi2, &ksi3, &nroot);

		if (nroot == 0) {
			answer.zero();
			return false;
		}

		if (nroot) {
			denom = (b3 + ksi1 * b4);
			if (fabs(denom) <= 1.0e-10) {
				eta1 = (4.0 * coords.at(xind) - a1 - ksi1 * a2) / (a3 + ksi1 * a4);
			}
			else {
				eta1 = (4.0 * coords.at(yind) - b1 - ksi1 * b2) / denom;
			}
		}

		if (nroot > 1) {
			double diff_ksi1, diff_eta1, diff_ksi2, diff_eta2, diff1, diff2;

			denom = b3 + ksi2 * b4;
			if (fabs(denom) <= 1.0e-10) {
				eta2 = (4.0 * coords.at(xind) - a1 - ksi2 * a2) / (a3 + ksi2 * a4);
			}
			else {
				eta2 = (4.0 * coords.at(yind) - b1 - ksi2 * b2) / denom;
			}

			// choose the one which seems to be closer to the parametric space (square <-1;1>x<-1;1>)
			diff_ksi1 = 0.0;
			if (ksi1 > 1.0) {
				diff_ksi1 = ksi1 - 1.0;
			}

			if (ksi1 < -1.0) {
				diff_ksi1 = ksi1 + 1.0;
			}

			diff_eta1 = 0.0;
			if (eta1 > 1.0) {
				diff_eta1 = eta1 - 1.0;
			}

			if (eta1 < -1.0) {
				diff_eta1 = eta1 + 1.0;
			}

			diff_ksi2 = 0.0;
			if (ksi2 > 1.0) {
				diff_ksi2 = ksi2 - 1.0;
			}

			if (ksi2 < -1.0) {
				diff_ksi2 = ksi2 + 1.0;
			}

			diff_eta2 = 0.0;
			if (eta2 > 1.0) {
				diff_eta2 = eta2 - 1.0;
			}

			if (eta2 < -1.0) {
				diff_eta2 = eta2 + 1.0;
			}

			diff1 = diff_ksi1 * diff_ksi1 + diff_eta1 * diff_eta1;
			diff2 = diff_ksi2 * diff_ksi2 + diff_eta2 * diff_eta2;

			// ksi2, eta2 seems to be closer
			if (diff1 > diff2) {
				ksi1 = ksi2;
				eta1 = eta2;
			}
		}

		answer.at(1) = ksi1;
		answer.at(2) = eta1;

		// test if inside
		bool inside = true;
		for (int i = 1; i <= 2; i++) {
			if (answer.at(i) < (-1. - POINT_TOL)) {
				answer.at(i) = -1.;
				inside = false;
			}
			else if (answer.at(i) > (1. + POINT_TOL)) {
				answer.at(i) = 1.;
				inside = false;
			}
		}

		return inside;
	}

	bool FEI2dQuadLin::inside(const FloatArray& lcoords) const
	{
		const double point_tol = 1.0e-3;
		bool inside = true;
		for (int i = 1; i <= 2; i++) {
			if (lcoords.at(i) < (-1. - point_tol)) {
				inside = false;
			}
			else if (lcoords.at(i) > (1. + point_tol)) {
				inside = false;
			}
		}

		return inside;
	}

	FloatMatrixF<2, 4> FEI2dQuadLin::evaldNdxi(const FloatArrayF<2>& lcoords)
	{
		const double& ksi = lcoords[0];
		const double& eta = lcoords[1];

		FloatMatrixF<2, 4> answer;
		// dn/dxi
		answer.at(1, 1) = 0.25 * (1. + eta);
		answer.at(1, 2) = -0.25 * (1. + eta);
		answer.at(1, 3) = -0.25 * (1. - eta);
		answer.at(1, 4) = 0.25 * (1. - eta);

		// dn/deta
		answer.at(2, 1) = 0.25 * (1. + ksi);
		answer.at(2, 2) = 0.25 * (1. - ksi);
		answer.at(2, 3) = -0.25 * (1. - ksi);
		answer.at(2, 4) = -0.25 * (1. + ksi);
		return answer;
	}

	
	std::unique_ptr<IntegrationRule>
		FEI2dQuadLin::giveIntegrationRule(int order) const
	{
		auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
		int points = iRule->getRequiredNumberOfIntegrationPoints(_Square, order + 2);
		iRule->SetUpPointsOnSquare(points, _Unknown);
		return std::move(iRule);
	}


	/*
	 * FEI2dQuadlinAxi element
	 */
	double FEI2dQuadLinAxi::giveTransformationJacobian(const FloatArray& lcoords) const
	{
		FloatArray N;
		this->evalN(N, lcoords);

		double r = 0.0;
		FEIElementGeometry* cellgeo = giveElemGeomety();
		for (int i = 1; i <= 4; i++) {
			double x = cellgeo->giveVertexCoordinates(i).at(1);
			r += x * N.at(i);
		}

		return r * FEI2dQuadLin::giveTransformationJacobian(lcoords);
	}

} // end namespace oofem
