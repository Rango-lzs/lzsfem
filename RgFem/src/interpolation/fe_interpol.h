/*****************************************************************//**
 * \file   fe_interpol.h
 * \brief
 *  单元插值基类
 *  1、
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef feinterpol_h
#define feinterpol_h

#include "femcore/fem_export.h"
#include "logger/error.h"
#include <memory>


namespace fem 
{
	class Element;
	class FloatArray;
	class FloatMatrix;
	class IntArray;
	class IntegrationRule;
	class FEIElementGeometry;

	/**
	 *
	 * 计算雅克比矩阵, 或者自然坐标到物理坐标的映射，需要知道单元形状信息, 此模块只计算插值函数相关信息
	 */
	class FEM_EXPORT FEInterpolation
	{
	private:		
		//the element to interpolated
		std::unique_ptr<FEIElementGeometry> m_elemGeom;

	public:
		FEInterpolation(std::unique_ptr<FEIElementGeometry> elemGeom);

		virtual ~FEInterpolation() = default;
		
		std::unique_ptr<FEIElementGeometry>& giveElemGeomety() const;

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
		virtual void evaldNdxi(FloatMatrix& answer, const FloatArray& lcoords) const = 0;

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
		virtual void local2global(const FloatArray& lcoords, FloatArray& gcoords) const = 0;
		/**
		 * Evaluates local coordinates from given global ones.
		 * If local coordinates cannot be found (generate elements, or point far outside geometry,
		 * then the center coordinate will be used as a last resort, and the return value will be zero.
		 * @param answer Contains evaluated local coordinates.
		 * @param gcoords Array containing global coordinates.
		 * @param cellgeo Underlying cell geometry.
		 * @return Nonzero is returned if point is within the element geometry, zero otherwise.
		 */
		virtual int global2local(const FloatArray& gcoords, FloatArray& lcoords) const = 0;
		
		/**
		 * Evaluates the determinant of the transformation.
		 * @param lcoords Array containing (local) coordinates.
		 * @param cellgeo Underlying cell geometry.
		 * @return Determinant of the transformation.
		 */
		virtual double calcJacobianDet(const FloatArray& lcoords) const;
		
		/**
		 * Gives the jacobian matrix at the local coordinates.
		 * @param jacobianMatrix The requested matrix.
		 * @param lcoords Local coordinates.
		 * @param cellgeo Element geometry.
		 */
		virtual void calcJacobianMatrixAt(FloatMatrix& jacobianMatrix, const FloatArray& lcoords) const
		{
			FEM_ERROR("Not overloaded.");
		}
	};
} // end namespace fem
#endif // feinterpol_h
