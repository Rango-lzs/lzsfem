
#ifndef geneigvalsolvertype_h
#define geneigvalsolvertype_h

namespace fem
{
	/**
	 * Types of general eigenvalue solvers.
	 */
	enum GenEigvalSolverType {
		GES_SubspaceIt,
		GES_InverseIt,
		GES_SLEPc
	};
} // end namespace fem
#endif // geneigvalsolvertype_h
