/*****************************************************************//**
 * \file   fe_interpol.h
 * \brief
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef feinterpol_h
#define feinterpol_h

#include "error.h"
#include "inputrecord.h"
#include "intarray.h"
#include "node.h"
#include "element.h"

namespace fem 
{
	class Element;
	class FloatArray;
	class FloatMatrix;
	class IntArray;
	class IntegrationRule;

	template <std::size_t N> class FloatArrayF;
	template <std::size_t N, std::size_t M> class FloatMatrixF;

	/**
	 * Class representing a general abstraction for finite element interpolation class.
	 * The boundary functions denote the (numbered) region that have 1 spatial dimension (i.e. edges) or 2 spatial dimensions.
	 */
	class FEM_EXPORT FEInterpolation
	{
	protected:
		int order = 0;

	public:
		FEInterpolation(int o) : order(o) { }
		virtual ~FEInterpolation() = default;
		
		/**
		 * Returns the interpolation order.
		 */
		int giveInterpolationOrder() const { return order; }

		/**
		 * Evaluates the array of interpolation functions (shape functions) at given point.
		 * @param answer Contains resulting array of evaluated interpolation functions.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 */
		virtual void evalN(FloatArray& answer, const FloatArray& lcoords) const = 0;
		/**
		 * Evaluates the matrix of derivatives of interpolation functions (shape functions) at given point.
		 * These derivatives are in global coordinate system (where the nodal coordinates are defined)
		 * @param answer Contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 * @return Determinant of the Jacobian.
		 */
		virtual double evaldNdx(FloatMatrix& answer, const FloatArray& lcoords) const = 0;
		/**
		 * Evaluates the matrix of second derivatives of interpolation functions (shape functions) at given point.
		 * These derivatives are in global coordinate system (where the nodal coordinates are defined)
		 * @param answer Contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 */
		virtual void evald2Ndx2(FloatMatrix& answer, const FloatArray& lcoords) const
		{
			FEM_ERROR("not implemented");
		}
		/**
		 * Evaluates the matrix of derivatives of interpolation functions (shape functions) at given point.
		 * These derivatives are wrt local (parent) coordinate system
		 * @param answer Contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxij.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 */
		virtual void evaldNdxi(FloatMatrix& answer, const FloatArray& lcoords) const
		{
			FEM_ERROR("not implemented");
		}
		/**
		 * Returns a matrix containing the local coordinates for each node corresponding to the interpolation
		 */
		virtual void giveLocalNodeCoords(FloatMatrix& answer) const
		{
			FEM_ERROR("FEInterpolation::giveLocalNodeCoords: not implemented");
		}
		/**
		 * Evaluates global coordinates from given local ones.
		 * @param answer Contains resulting global coordinates.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 */
		virtual void local2global(FloatArray& answer, const FloatArray& lcoords) const = 0;
		/**
		 * Evaluates local coordinates from given global ones.
		 * If local coordinates cannot be found (generate elements, or point far outside geometry,
		 * then the center coordinate will be used as a last resort, and the return value will be zero.
		 * @param answer Contains evaluated local coordinates.
		 * @param gcoords Array containing global coordinates.
		 * @param cellgeo Underlying cell geometry.
		 * @return Nonzero is returned if point is within the element geometry, zero otherwise.
		 */
		virtual int global2local(FloatArray& answer, const FloatArray& gcoords) const = 0;
		/**
		 * Evaluates the determinant of the transformation.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 * @return Determinant of the transformation.
		 */
		virtual double giveTransformationJacobian(const FloatArray& lcoords) const;
		/**
		 * Gives the jacobian matrix at the local coordinates.
		 * @param jacobianMatrix The requested matrix.
		 * @param lcoords Local coordinates.
		 * @param cellgeo Element geometry.
		 */
		virtual void giveJacobianMatrixAt(FloatMatrix& jacobianMatrix, const FloatArray& lcoords) const
		{
			FEM_ERROR("Not overloaded.");
		}

		/**
		 * Sets up a suitable integration rule for numerical integrating over volume.
		 * The required polynomial order for the determinant of the jacobian is added automatically.
		 * @param order Polynomial order of integrand (should NOT including determinant of jacobian).
		 */
		virtual std::unique_ptr<IntegrationRule> giveIntegrationRule(int order) const;
		//@}

		/** @name Edge boundary functions.
		 * Provide interpolation services for boundary edges (entity of dimension 1)
		 */
		 //@{
		 /**
		  * Evaluates the basis functions on the requested boundary.
		  * Only basis functions that are nonzero anywhere on the boundary are given. Ordering can be obtained from giveBoundaryNodes.
		  * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
		  * @param answer Basis functions Array to be filled with the boundary nodes.
		  * @param boundary Boundary number.
		  * @param lcoords The local coordinates (on the boundary local coordinate system).
		  * @param cellgeo Underlying cell geometry.
		  * @todo
		  */
		virtual void boundaryEdgeEvalN(FloatArray& answer, int boundary, const FloatArray& lcoords) const = 0;
		/**
		 * Evaluates the determinant of the transformation Jacobian on the requested boundary.
		 * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
		 * @param boundary Boundary number.
		 * @param lcoords The local coordinates (on the boundary local coordinate system).
		 * @param cellgeo Underlying cell geometry.
		 * @return The determinant of the boundary transformation Jacobian.
		 */
		virtual double boundaryEdgeGiveTransformationJacobian(int boundary, const FloatArray& lcoords) const = 0;
		/**
		 * Maps the local boundary coordinates to global.
		 * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
		 * @param answer Global coordinates.
		 * @param boundary Boundary number.
		 * @param lcoords The local coordinates (on the boundary local coordinate system).
		 * @param cellgeo Underlying cell geometry.
		 */
		virtual void boundaryEdgeLocal2Global(FloatArray& answer, int boundary, const FloatArray& lcoords) const = 0;
		
		/**
		 * Sets up a suitable integration rule for integrating over the requested boundary.
		 * The required polynomial order for the determinant of the jacobian is added automatically.
		 * @param order Polynomial order of the integrand (should NOT including determinant of jacobian).
		 * @param boundary Boundary number.
		 */
		virtual std::unique_ptr<IntegrationRule> giveBoundaryEdgeIntegrationRule(int order, int boundary) const;
		/**
		 * Gives the boundary nodes for requested boundary number.
		 * @param answer Array to be filled with the boundary nodes.
		 * @param boundary Boundary number.
		 */
		virtual IntArray boundaryEdgeGiveNodes(int boundary) const = 0;
		//@}

		/**@name Surface interpolation services
		 * Provide interpolation services for boundary edges (entities of dimension 2)
		 */
		 //@{
		 /**
		  * Evaluates the array of edge interpolation functions (shape functions) at given point.
		  * @param answer Contains resulting array of evaluated interpolation functions.
		  * @param isurf Surface number.
		  * @param lcoords Array containing (local) coordinates.
		  * @param cellgeo Underlying cell geometry.
		  */
		virtual void boundarySurfaceEvalN(FloatArray& answer, int isurf, const FloatArray& lcoords) const = 0;
		/**
		 * Evaluates the matrix of derivatives of edge interpolation functions (shape functions) at given point.
		 * These derivatives are in global coordinate system (where the nodal coordinates are defined).
		 * @param answer Contains resulting matrix of derivatives, the member at i,j position contains value of dNj/dxi.
		 * @param isurf Determines the surface number.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 */
		virtual void boundarySurfaceEvaldNdx(FloatMatrix& answer, int isurf, const FloatArray& lcoords) const = 0;
		/**
		 * Evaluates the normal out of the surface at given point.
		 * @param answer Contains resulting normal vector.
		 * @param isurf Determines the surface number.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 * @return Surface mapping jacobian.
		 */
		virtual double boundarySurfaceEvalNormal(FloatArray& answer, int isurf, const FloatArray& lcoords) const = 0;

		/**
		 * Evaluates edge global coordinates from given local ones.
		 * These derivatives are in global coordinate system (where the nodal coordinates are defined).
		 * @param answer Contains resulting global coordinates.
		 * @param isurf Determines the surface number.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 */
		virtual void boundarySurfaceLocal2global(FloatArray& answer, int isurf, const FloatArray& lcoords) const = 0;
		/**
		 * Evaluates the edge jacobian of transformation between local and global coordinates.
		 * @param isurf Determines the surface number.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 * @return Determinant of the transformation.
		 */
		virtual double boundarySurfaceGiveTransformationJacobian(int isurf, const FloatArray& lcoords) const = 0;
		
		/**
		 * Sets up a suitable integration rule for integrating over the requested boundary.
		 * The required polynomial order for the determinant of the jacobian is added automatically.
		 * @param order Polynomial order of the integrand (should NOT including determinant of jacobian).
		 * @param boundary Boundary number.
		 */
		virtual std::unique_ptr<IntegrationRule> giveBoundarySurfaceIntegrationRule(int order, int boundary) const;
		/**
		 * Gives the boundary nodes for requested boundary number.
		 * @param answer Array to be filled with the boundary nodes.
		 * @param boundary Boundary number.
		 */
		virtual IntArray boundarySurfaceGiveNodes(int boundary) const = 0;
		//@}

		/** @name General boundary interpolation functions.
		 * Provide interpolation servises for boundary entities with one dimension lower than the receiver interpolation.
		 * Typically these are mapped to boundaryEdge and boundarySurface methods depending on dimension.
		 *
		 */
		 //@{
		 /**
		  * Gives the boundary nodes for requested boundary number.
		  * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
		  * @param answer Array to be filled with the boundary nodes.
		  * @param boundary Boundary number.
		  */
		virtual IntArray boundaryGiveNodes(int boundary) const = 0;
		/**
		 * Evaluates the basis functions on the requested boundary.
		 * Only basis functions that are nonzero anywhere on the boundary are given. Ordering can be obtained from giveBoundaryNodes.
		 * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
		 * @param answer Basis functions Array to be filled with the boundary nodes.
		 * @param boundary Boundary number.
		 * @param lcoords The local coordinates (on the boundary local coordinate system).
		 * @param cellgeo Underlying cell geometry.
		 */
		virtual void boundaryEvalN(FloatArray& answer, int boundary, const FloatArray& lcoords) const = 0;
		/**
		 * Evaluates the normal on the requested boundary.
		 * @param answer The evaluated normal.
		 * @param boundary Boundary number.
		 * @param lcoords The local coordinates (on the boundary local coordinate system).
		 * @param cellgeo Underlying cell geometry.
		 * @return The boundary transformation Jacobian.
		 */
		virtual double boundaryEvalNormal(FloatArray& answer, int boundary, const FloatArray& lcoords) const = 0;
		/**
		 * Evaluates the determinant of the transformation Jacobian on the requested boundary.
		 * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
		 * @param boundary Boundary number.
		 * @param lcoords The local coordinates (on the boundary local coordinate system).
		 * @param cellgeo Underlying cell geometry.
		 * @return The determinant of the boundary transformation Jacobian.
		 */
		virtual double boundaryGiveTransformationJacobian(int boundary, const FloatArray& lcoords) const = 0;
		/**
		 * Maps the local boundary coordinates to global.
		 * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
		 * @param answer Global coordinates.
		 * @param boundary Boundary number.
		 * @param lcoords The local coordinates (on the boundary local coordinate system).
		 * @param cellgeo Underlying cell geometry.
		 */
		virtual void boundaryLocal2Global(FloatArray& answer, int boundary, const FloatArray& lcoords) const = 0;
		/**
		 * Computes the integral @f$ \int_S n \cdot x \mathrm{d}s @f$.
		 * @param boundary Boundary number.
		 * @param cellgeo Underlying cell geometry.
		 * @return Evaluated integral.
		 */
		virtual double evalNXIntegral(int boundary) const
		{
			FEM_ERROR("Not implemented");
			return 0.;
		}
		
		/**
		 * Sets up a suitable integration rule for integrating over the requested boundary.
		 * The required polynomial order for the determinant of the jacobian is added automatically.
		 * @param order Polynomial order of the integrand (should NOT including determinant of jacobian).
		 * @param boundary Boundary number.
		 */
		virtual std::unique_ptr<IntegrationRule> giveBoundaryIntegrationRule(int order, int boundary) const;
		//@}

		std::string errorInfo(const char* func) const { return func; } ///@todo Class name?
	};
} // end namespace fem
#endif // feinterpol_h
