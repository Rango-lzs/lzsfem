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

		int giveNsd() const override { return 1; }

		/**
		 * Computes the exact length.
		 * @param cellgeo Cell geometry for the element.
		 * @return Length of geometry.
		 */
		virtual double giveLength() const = 0;
	};
} // end namespace fem
#endif // feinterpol1d_h
