/*****************************************************************//**
 * \file   fei2d_line_linear.h
 * \brief
 *
 * \author Leizs
 * \date   November 2023
 *********************************************************************/

#ifndef fei2dlinelin_h
#define fei2dlinelin_h

#include "feinterpol2d.h"

namespace fem
{
	/**
	 * Class representing a 2d line with linear interpolation.
	 * @todo{Some more routines to be implemented here}
	 * @author Mikael Öhman
	 */
	class FEM_EXPORT FEI2dLineLin : public FEInterpolation2d
	{
	public:
		FEI2dLineLin(int ind1, int ind2) : FEInterpolation2d(1, ind1, ind2) { }

		double giveArea() const override { return 0.0; }

		void local2global(FloatArray& answer, const FloatArray& lcoords) const override;
		int global2local(FloatArray& answer, const FloatArray& gcoords) const override;

		// "Bulk"
		static FloatArrayF<2> evalN(double xi);

		void evalN(FloatArray& answer, const FloatArray& lcoords) const override;
		double evaldNdx(FloatMatrix& answer, const FloatArray& lcoords) const override;
		void evaldNdxi(FloatMatrix& answer, const FloatArray& lcoords) const override;
		double giveTransformationJacobian(const FloatArray& lcoords) const override;

	protected:
		double edgeComputeLength(const IntArray& edgeNodes) const;
	};
} // end namespace oofem
#endif // fei2dlinelin_h
