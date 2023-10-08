
#ifndef bctype_h
#define bctype_h

namespace fem
{
	/// Type representing the type of bc.
	enum bcType {
		UnknownBT,      ///< Unknown.
		DirichletBT,    ///< Prescribed value.
		TransmissionBC, ///< Neumann type (prescribed flux).
		ConvectionBC,   ///< Newton type - transfer coefficient
		SlipWithFriction,
		PenetrationWithResistance,
		OutFlowBC,
		RadiationBC     ///< Stefan-Boltzmann law.
	};
} // end namespace fem
#endif // bctype_h
