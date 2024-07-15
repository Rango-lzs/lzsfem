/*****************************************************************//**
 * \file   fei2d_line_quadratic.h
 * \brief
 *
 * \author Leizs
 * \date   November 2023
 *********************************************************************/

#ifndef fei2dlinequad_h
#define fei2dlinequad_h

#include "feinterpol2d.h"

namespace fem
{
	/**
	 * Class representing a 2d line quadratic interpolation.
	 * @todo Some more routines to be implemented here.
	 * @author Mikael Öhman
	      .
	      .       .2
	      .     .
	      .   .
	      . .1
	      . . . . . . .
	 */

	//二维线单元，一般在局部坐标系下插值
	class FEM_EXPORT FEI2dLineQuad : public FEInterpolation2d
	{
	public:
		FEI2dLineQuad(int ind1, int ind2) : FEInterpolation2d(2, ind1, ind2) { }

		double giveArea() const override { return 0.0; }

		void local2global(FloatArray& answer, const FloatArray& lcoords) const override;
		int global2local(FloatArray& answer, const FloatArray& gcoords) const override;

		// "Bulk"
		void evalN(FloatArray& answer, const FloatArray& lcoords) const override;
		double evaldNdx(FloatMatrix& answer, const FloatArray& lcoords) const override;
		void evaldNdxi(FloatMatrix& answer, const FloatArray& lcoords) const override;
		void giveJacobianMatrixAt(FloatMatrix& jacobianMatrix, const FloatArray& lcoords) const override;
		double giveTransformationJacobian(const FloatArray& lcoords) const override;

	protected:
		double edgeComputeLength(const IntArray& edgeNodes) const;
	};
} // end namespace oofem
#endif // fei2dlinequad_h
