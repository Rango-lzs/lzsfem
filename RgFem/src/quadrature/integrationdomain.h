
#ifndef integrationdomain_h
#define integrationdomain_h

namespace fem
{
	/**
	 * Used by integrator class to supply
	 * integration points for proper domain to be
	 * integrated (Area,Volume and its shape)
	 */
	enum integrationDomain {
		_UnknownIntegrationDomain,
		_Point,
		_Line,
		_Triangle,
		_Square,
		_Cube,
		_Tetrahedra,
		_Wedge,
		_Embedded2dLine,
		_3dDegShell,
	};
} // end namespace fem
#endif // integrationdomain_h
