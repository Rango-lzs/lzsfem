/*****************************************************************//**
 * \file   fei2d_line_hermite.h
 * \brief
 *
 * \author Leizs
 * \date   November 2023
 *********************************************************************/

#ifndef fei2dlinehermite_h
#define fei2dlinehermite_h

#include "feinterpol2d.h"

namespace fem
{
	/**
	 * Class representing a 2d line with Hermitian interpolation.
	 * The order used is cubic, quadratic, cubic, quadratic.
	 * The functions that need geometric information, a linear interpolation is assumed (for geometry). This means functions such as evaldNdx.
	 * @author Mikael Öhman
	 */
	class FEM_EXPORT FEI2dLineHermite : public FEInterpolation2d
	{
	public:
		FEI2dLineHermite(int ind1, int ind2) : FEInterpolation2d(1, ind1, ind2) { }
	
		double giveArea() const override { return 0.0; }
		double giveLength() const;

		void local2global(FloatArray& answer, const FloatArray& lcoords) const override;
		int global2local(FloatArray& answer, const FloatArray& gcoords) const override;

		// "Bulk"
		void evalN(FloatArray& answer, const FloatArray& lcoords) const override;
		double evaldNdx(FloatMatrix& answer, const FloatArray& lcoords) const override;
		double giveTransformationJacobian(const FloatArray& lcoords) const override;


		void edgeEvaldNds(FloatArray& answer, int iedge, const FloatArray& lcoords) const;
		void edgeEvald2Nds2(FloatArray& answer, int iedge, const FloatArray& lcoords) const;
		double edgeEvalNormal(FloatArray& normal, int iedge, const FloatArray& lcoords) const;

	};
} // end namespace fem
#endif // fei2dlinehermite_h
