
#ifndef fmode_h
#define fmode_h

namespace fem
{
	/**
	 * Type representing the type of formulation (total or updated) of non-linear computation.
	 */
	enum fMode {
		UNKNOWN = 0, ///< Unknown.
		TL = 1,      ///< Total Lagrange.
		AL = 2,      ///< Updated Lagrange.
	};
} // end namespace fem
#endif // fmode_h
