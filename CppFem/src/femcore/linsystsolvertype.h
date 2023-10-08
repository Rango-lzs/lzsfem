
#ifndef linsystsolvertype_h
#define linsystsolvertype_h

namespace fem
{
	/**
	 * The values of this type should be related not to specific solvers,
	 * but more to specific packages that provide linear solver interface
	 * (possibly with many solver types) and are represented by a class
	 * derived from SparseLinearSystemNM.
	 * The selection of particular solver from package should be done using keywords,
	 * related to particular package.
	 */
	enum LinSystSolverType {
		ST_Direct = 0,
		ST_IML = 1,
		ST_Spooles = 2,
		ST_Petsc = 3,
		ST_DSS = 4,
		ST_Feti = 5,
		ST_MKLPardiso = 6,
		ST_SuperLU_MT = 7,
		ST_PardisoProjectOrg = 8 // experimental
	};
} // end namespace fem
#endif // linsystsolvertype_h
