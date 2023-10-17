/*****************************************************************//**
 * \file   feinterpol2d.cpp
 * \brief
 *
 * \author Leizs
 * \date   October 2023
 *********************************************************************/

#include "feinterpol2d.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"

namespace fem
{
	double FEInterpolation2d::giveArea() const
	{
		FEM_ERROR("Not implemented in subclass.");
		return 0;
	}

#define POINT_TOL 1.e-3

	int FEInterpolation2d::global2local(FloatArray& answer, const FloatArray& gcoords) const
	{
		FloatArray res, delta, guess, lcoords_guess;
		FloatMatrix jac;
		double convergence_limit, error = 0.0;

		// find a suitable convergence limit
		convergence_limit = 1e-6 * this->giveCharacteristicLength();

		// setup initial guess
		lcoords_guess.resize(2);
		lcoords_guess.zero();

		// apply Newton-Raphson to solve the problem
		for (int nite = 0; nite < 10; nite++) {
			// compute the residual
			this->local2global(guess, lcoords_guess);
			res = { gcoords(0) - guess(0), gcoords(1) - guess(1) };

			// check for convergence
			error = res.computeNorm();
			if (error < convergence_limit) {
				break;
			}

			// compute the corrections
			this->giveJacobianMatrixAt(jac, lcoords_guess);
			jac.solveForRhs(res, delta);

			// update guess
			lcoords_guess.add(delta);
		}
		if (error > convergence_limit) { // Imperfect, could give false negatives.
			FEM_WARNING("Failed convergence");
			answer = { 1. / 3., 1. / 3. };
			return false;
		}

		answer = { lcoords_guess(0), lcoords_guess(1) };

		return inside(answer);
	}

	void FEInterpolation2d::giveJacobianMatrixAt(FloatMatrix& jacobianMatrix, const FloatArray& lcoords) const
		// Returns the jacobian matrix  J (x,y)/(ksi,eta)  of the receiver.
	{
		FloatMatrix dn;

		jacobianMatrix.resize(2, 2);
		jacobianMatrix.zero();

		this->evaldNdxi(dn, lcoords);

		for (int i = 1; i <= dn.giveNumberOfRows(); i++) {
			double x = 0;//cellgeo.giveVertexCoordinates(i).at(xind);
			double y = 0; //cellgeo.giveVertexCoordinates(i).at(yind);

			jacobianMatrix.at(1, 1) += dn.at(i, 1) * x;
			jacobianMatrix.at(2, 1) += dn.at(i, 1) * y;
			jacobianMatrix.at(1, 2) += dn.at(i, 2) * x;
			jacobianMatrix.at(2, 2) += dn.at(i, 2) * y;
		}
	}

	bool FEInterpolation2d::inside(const FloatArray& lcoords) const
	{
		FEM_ERROR("Not implemented.")
	}
} // end namespace fem
