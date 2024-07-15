/*****************************************************************//**
 * \file   fei1d_line_hermite.h
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
	 * Class representing a 1d Hermitian cubic isoparametric interpolation.
	 * @author Mikael Öhman
	 */
	class FEM_EXPORT FEI1dHermite : public FEInterpolation1d
	{
	protected:
		int cindx;

	public:
		FEI1dHermite(int coordIndx) : FEInterpolation1d(2), cindx(coordIndx) { }

		double giveLength() const override;

		std::pair<double, FloatMatrixF<1, 4>> evaldNdx(double ksi) const;
		FloatMatrixF<1, 4> evald2Ndx2(double ksi) const;

		void evalN(FloatArray& answer, const FloatArray& lcoords) const override;
		double evaldNdx(FloatMatrix& answer, const FloatArray& lcoords) const override;
		void evald2Ndx2(FloatMatrix& answer, const FloatArray& lcoords) const override;
		void local2global(FloatArray& answer, const FloatArray& lcoords) const override;
		int  global2local(FloatArray& answer, const FloatArray& lcoords) const override;
		double giveTransformationJacobian(const FloatArray& lcoords) const override;
	};
} // end namespace oofem
#endif
