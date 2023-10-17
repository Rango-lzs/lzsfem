/*****************************************************************//**
 * \file   feinterpol2d.h
 * \brief  
 * 
 * \author Leizs
 * \date   October 2023
 *********************************************************************/

#ifndef feinterpol2d_h
#define feinterpol2d_h

#include "fe_interpol.h"
#include "mathfem.h"

namespace fem
{
	/**
	 * Class representing a general abstraction for surface finite element interpolation class.
	 */
	class FEM_EXPORT FEInterpolation2d : public FEInterpolation
	{
	protected:
		int xind, yind;

	public:
		FEInterpolation2d(int o, int ind1, int ind2) : FEInterpolation(o), xind(ind1), yind(ind2) { }

		int giveNsd() const override { return 2; }

		/**
		 * Computes the exact area.
		 * @param cellgeo Cell geometry for the element.
		 * @return Area of geometry.
		 */
		virtual double giveArea() const;

		/**
		 * Returns a characteristic length of the geometry, typically a diagonal or edge length.
		 * @param cellgeo Underlying cell geometry.
		 * @return Square root of area.
		 */
		virtual double giveCharacteristicLength() const { return sqrt(this->giveArea()); }

		/**
		 * Default implementation using Newton's method to find the local coordinates.
		 * Can be overloaded if desired.
		 */
		int global2local(FloatArray& answer, const FloatArray& gcoords) const override;

		void giveJacobianMatrixAt(FloatMatrix& jacobianMatrix, const FloatArray& lcoords) const override;

		virtual bool inside(const FloatArray& lcoords) const;
	};
} // end namespace fem
#endif // feinterpol2d_h
