
#ifndef crosssectextension_h
#define crosssectextension_h

namespace fem
{
	/// Type representing cross section extension for run time testing.
	enum CrossSectExtension {
		CS_StructuralCapability, ///< Structural capability.
		CS_StructuralInterfaceCapability, ///< Structural interface capability.
		CS_HeatCapability, ///< Heat capability.
		CS_LatticeStructuralCapability, ///< Structural lattice capability.
	};
} // end namespace fem
#endif // crosssectextension_h
