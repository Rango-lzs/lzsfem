
#ifndef contextoutputmode_h
#define contextoutputmode_h

namespace fem
{
	enum ContextOutputMode {
		COM_NoContext,     ///< No context.
		COM_Always,        ///< Enable for post-processing.
		COM_Required,      ///< If required (for backtracking computation).
		COM_UserDefined,   ///< Input attribute of domain (each n-th step).
	};
} // end namespace oofem
#endif // contextoutputmode_h
