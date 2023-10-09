#ifndef feinterpol1d_h
#define feinterpol1d_h

#include "fe_interpol.h"
#include <stdexcept>

namespace fem
{
	/**
	 * Class representing a general abstraction for finite element interpolation class.
	 */
	class FEM_EXPORT FEInterpolation1d : public FEInterpolation
	{
	public:
		FEInterpolation1d(int o) : FEInterpolation(o) { }
		//int giveNsd() const override { return 1; }

		//FloatArray giveParametricCenter() const override { return {0.}; }

		IntArray boundaryGiveNodes(int boundary) const override;
		void boundaryEvalN(FloatArray& answer, int boundary, const FloatArray& lcoords) const override;
		double boundaryEvalNormal(FloatArray& answer, int boundary, const FloatArray& lcoords) const  override;
		double boundaryGiveTransformationJacobian(int boundary, const FloatArray& lcoords) const override;
		void boundaryLocal2Global(FloatArray& answer, int boundary, const FloatArray& lcoords) const override;

		/**@name Surface interpolation services */
		//@{
		void boundarySurfaceEvalN(FloatArray& answer, int isurf, const FloatArray& lcoords) const override
		{
			FEM_ERROR("Functions not supported for this interpolator.");
		}
		void boundarySurfaceEvaldNdx(FloatMatrix& answer, int isurf, const FloatArray& lcoords) const override
		{
			FEM_ERROR("Functions not supported for this interpolator.");
		}
		double boundarySurfaceEvalNormal(FloatArray& answer, int isurf, const FloatArray& lcoords) const override
		{
			FEM_ERROR("Functions not supported for this interpolator.");
		}
		void boundarySurfaceLocal2global(FloatArray& answer, int isurf, const FloatArray& lcoords) const override
		{
			FEM_ERROR("Functions not supported for this interpolator.");
		}
		double boundarySurfaceGiveTransformationJacobian(int isurf, const FloatArray& lcoords) const override
		{
			FEM_ERROR("Functions not supported for this interpolator.");
		}
		IntArray boundarySurfaceGiveNodes(int boundary) const override
		{
			throw std::runtime_error("Functions not supported for this interpolator.");
		}
		//@}

		/**
		 * Computes the exact length.
		 * @param cellgeo Cell geometry for the element.
		 * @return Length of geometry.
		 */
		virtual double giveLength(const FEICellGeometry& cellgeo) const
		{
			FEM_ERROR("Not implemented in subclass.");
			return 0;
		}

		std::unique_ptr<IntegrationRule> giveIntegrationRule(int order) const override;
		std::unique_ptr<IntegrationRule> giveBoundaryIntegrationRule(int order, int boundary) const override;
		std::unique_ptr<IntegrationRule> giveBoundaryEdgeIntegrationRule(int order, int boundary) const override;
	};
} // end namespace fem
#endif // feinterpol1d_h
