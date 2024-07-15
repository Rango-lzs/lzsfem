#ifndef varscaletype_h
#define varscaletype_h

namespace fem
{
	/// Type determining the scale corresponding to particular variable.
	enum VarScaleType {
		VST_Length,
		VST_Velocity,
		VST_Time,
		VST_Density,
		VST_Pressure,
		VST_Force,
		VST_Viscosity,
		VST_ReynoldsNumber
	};
} // end namespace fem
#endif // varscaletype_h
