#ifndef timediscretizationtype_h
#define timediscretizationtype_h

namespace fem
{
	/// Time discretization used by transient solvers.
	enum TimeDiscretizationType {
		TD_Unspecified = -1,     ///< Unspecified
		TD_ThreePointBackward = 0,     ///< Three-point Backward Euler method
		TD_TwoPointBackward = 1,     ///< Two-point Backward Euler method
		TD_Newmark = 2,     ///< Newmark-beta method
		TD_Wilson = 3,     ///< Wilson-theta method
		TD_Explicit = 4,     ///< Central difference
	};
} // end namespace fem
#endif // timediscretizationtype_h
