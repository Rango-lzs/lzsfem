#ifndef feinterpol3d_h
#define feinterpol3d_h

#include "fe_interpol.h"

namespace fem
{
	class FEICellGeometry;
	/**
	 * Class representing a general abstraction for surface finite element interpolation class.
	 */
	class FEM_EXPORT FEInterpolation3d : public FEInterpolation
	{
	public:
		FEInterpolation3d(int o) : FEInterpolation(o) { }
		int giveNsd() const { return 3; }

		/**
		 * Computes the exact volume.
		 * @param cellgeo Cell geometry for the element.
		 * @return Volume of geometry.
		 */
		virtual double giveVolume(const FEICellGeometry& cellgeo) const;

		/**@name Edge interpolation services */
		//@{
		/**
		 * Evaluates the array of edge interpolation functions (shape functions) at given point.
		 * @param answer Contains resulting array of evaluated interpolation functions.
		 * @param iedge Edge number.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 */
		virtual void edgeEvalN(FloatArray& answer, int iedge, const FloatArray& lcoords, const FEICellGeometry& cellgeo) const = 0;
		/**
		 * Evaluates the matrix of derivatives of edge interpolation functions (shape functions) at given point.
		 * These derivatives are in global coordinate system (where the nodal coordinates are defined)
		 * @param answer Contains resulting matrix of derivatives, the member at i,j position contains value of dNj/dxi.
		 * @param iedge Determines the edge number.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 */
		virtual void edgeEvaldNdx(FloatMatrix& answer, int iedge, const FloatArray& lcoords, const FEICellGeometry& cellgeo) const = 0;
		/**
		 * Evaluates the matrix of derivatives of edge interpolation functions (shape functions) at given point.
		 * These derivatives are in local (parent) coordinate system
		 * @param answer Contains resulting matrix of derivatives, the member at i,j position contains value of dN_j/dxi_i.
		 * @param iedge Determines the edge number.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 */
		virtual void edgeEvaldNdxi(FloatArray& answer, int iedge, const FloatArray& lcoords, const FEICellGeometry& cellgeo) const;
		/**
		 * Evaluates edge global coordinates from given local ones.
		 * These derivatives are in global coordinate system (where the nodal coordinates are defined).
		 * @param answer Contains resulting global coordinates.
		 * @param iedge Determines edge number.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 */
		virtual void edgeLocal2global(FloatArray& answer, int iedge, const FloatArray& lcoords, const FEICellGeometry& cellgeo) const = 0;
		/**
		 * Evaluates the edge jacobian of transformation between local and global coordinates.
		 * @param iedge Determines edge number.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 * @return Determinant of the transformation.
		 */
		virtual double edgeGiveTransformationJacobian(int iedge, const FloatArray& lcoords, const FEICellGeometry& cellgeo) const = 0;

		virtual IntArray computeLocalEdgeMapping(int iedge) const = 0;
		IntArray computeEdgeMapping(const IntArray& elemNodes, int iedge) const;
		//@}

		/**@name Surface interpolation services */
		//@{
		/**
		 * Evaluates the array of edge interpolation functions (shape functions) at given point.
		 * @param answer Contains resulting array of evaluated interpolation functions.
		 * @param isurf Surface number.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 */
		virtual void surfaceEvalN(FloatArray& answer, int isurf, const FloatArray& lcoords, const FEICellGeometry& cellgeo) const = 0;
		/**
		 * Evaluates the matrix of derivatives of edge interpolation functions (shape functions) at given point.
		 * These derivatives are in global coordinate system (where the nodal coordinates are defined).
		 * @param answer Contains resulting matrix of derivatives, the member at i,j position contains value of dNj/dxi.
		 * @param isurf Determines the surface number.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 */
		virtual void surfaceEvaldNdx(FloatMatrix& answer, int isurf, const FloatArray& lcoords, const FEICellGeometry& cellgeo) const;
		/**
		 * Evaluates the normal out of the surface at given point.
		 * @param answer Contains resulting normal vector.
		 * @param isurf Determines the surface number.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 * @return Surface mapping jacobian.
		 */
		virtual double surfaceEvalNormal(FloatArray& answer, int isurf, const FloatArray& lcoords, const FEICellGeometry& cellgeo) const;

		/**
		 * Evaluates edge global coordinates from given local ones.
		 * These derivatives are in global coordinate system (where the nodal coordinates are defined).
		 * @param answer Contains resulting global coordinates.
		 * @param isurf Determines the surface number.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 */
		virtual void surfaceLocal2global(FloatArray& answer, int isurf, const FloatArray& lcoords, const FEICellGeometry& cellgeo) const = 0;
		/**
		 * Evaluates the edge jacobian of transformation between local and global coordinates.
		 * @param isurf Determines the surface number.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 * @return Determinant of the transformation.
		 */
		virtual double surfaceGiveTransformationJacobian(int isurf, const FloatArray& lcoords, const FEICellGeometry& cellgeo) const = 0;

		virtual IntArray computeLocalSurfaceMapping(int isurf) const = 0;
		IntArray computeSurfaceMapping(const IntArray& elemNodes, int isurf) const;
		//@}
	};
} // end namespace fem
#endif // feinterpol3d_h
