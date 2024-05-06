/*****************************************************************//**
 * \file   fei1d_line_quad.h
 * \brief
 *
 * \author Leizs
 * \date   November 2023
 *********************************************************************/

#ifndef fei1dquad_h
#define fei1dquad_h

#include "feinterpol1d.h"

namespace fem
{
	/**
	 * Class representing a 1d linear isoparametric interpolation.
	 */
	class FEM_EXPORT FEI1dQuad : public FEInterpolation1d
	{
	protected:
		int cindx;
	public:
		FEI1dQuad(int coordIndx) : FEInterpolation1d(2), cindx(coordIndx) { }

		double giveLength() const override;

		static FloatArrayF<3> evalN(double ksi);
		std::pair<double, FloatMatrixF<1, 3>> evaldNdx(double ksi) const;

		void evalN(FloatArray& answer, const FloatArray& lcoords) const override;
		double evaldNdx(FloatMatrix& answer, const FloatArray& lcoords) const override;
		void local2global(FloatArray& answer, const FloatArray& lcoords) const override;
		int  global2local(FloatArray& answer, const FloatArray& lcoords) const override;
		double giveTransformationJacobian(const FloatArray& lcoords) const override;		
	};
} // end namespace oofem
#endif
