/*****************************************************************//**
 * \file   fei2d_quad_linear.h
 * \brief
 *
 * \author Leizs
 * \date   October 2023
 *********************************************************************/

#ifndef fei2dquadlin_h
#define fei2dquadlin_h

#include "feinterpol2d.h"

namespace fem
{
	/**
	 * Class representing a 2d isoparametric linear interpolation based on natural coordinates
	 * for quadrilateral elements.
	 */
	class FEM_EXPORT FEI2dQuadLin : public FEInterpolation2d
	{
	public:
		FEI2dQuadLin(int ind1, int ind2) : FEInterpolation2d(1, ind1, ind2) { }

		integrationDomain giveIntegrationDomain() const override { return _Square; }

		//Element_Geometry_Type giveGeometryType() const override { return EGT_quad_1; }
		
		double giveArea() const override;

		// Bulk
		static FloatArrayF<4> evalN(const FloatArrayF<2>& lcoords);

		static FloatMatrixF<2, 4> evaldNdxi(const FloatArrayF<2>& lcoords);

		std::pair<double, FloatMatrixF<2, 4>> evaldNdx(const FloatArrayF<2>& lcoords) const;

		void evalN(FloatArray& answer, const FloatArray& lcoords) const override;

		void evaldNdxi(FloatMatrix& answer, const FloatArray& lcoords) const override;

		double evaldNdx(FloatMatrix& answer, const FloatArray& lcoords) const override;

		void local2global(FloatArray& answer, const FloatArray& lcoords) const override;

		int global2local(FloatArray& answer, const FloatArray& lcoords) const override;
		
		bool inside(const FloatArray& lcoords) const override;

		std::unique_ptr<IntegrationRule> giveIntegrationRule(int order) const override;
	};

	/**
	 * Class representing a 2d isoparametric linear interpolation based on natural coordinates
	 * for quadrilateral elements in axisymmetric setting.
	 */
	class FEM_EXPORT FEI2dQuadLinAxi : public FEI2dQuadLin
	{
	public:
		FEI2dQuadLinAxi(int ind1, int ind2) : FEI2dQuadLin(ind1, ind2) { }

		double giveTransformationJacobian(const FloatArray& lcoords) const override;		
	};

} // end namespace fem
#endif // fei2dquadlin_h
