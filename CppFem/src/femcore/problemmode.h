
#ifndef problemmode_h
#define problemmode_h

namespace fem
{
	enum problemMode {
		_processor,
		_postProcessor
	};

	/// Corresponds to macro- and micro-problem in multiscale simulations.
	enum problemScale {
		macroScale,
		microScale
	};
} // end namespace oofem
#endif // problemmode_h
