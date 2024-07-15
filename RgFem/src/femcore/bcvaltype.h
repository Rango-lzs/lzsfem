#ifndef bcvaltype_h
#define bcvaltype_h

namespace fem
{
	/// Type determining the type of general boundary condition.
	enum bcValType {
		UnknownBVT,
		TemperatureBVT,
		ForceLoadBVT,
		PressureBVT,
		HumidityBVT,
		VelocityBVT,
		DisplacementBVT,
		EigenstrainBVT,
		ReinforceBVT,
	};
} // end namespace fem
#endif // bcvaltype_h
