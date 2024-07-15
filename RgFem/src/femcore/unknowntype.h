#ifndef unknowntype_h
#define unknowntype_h
//#include "internalstatetype.h"
#include "enumitem.h"

namespace fem
{
#define UnknownType_DEF \
    ENUM_ITEM_WITH_VALUE(DisplacementVector, 1) \
    ENUM_ITEM_WITH_VALUE(GeneralizedDisplacementVector, 2) \
    ENUM_ITEM_WITH_VALUE(FluxVector, 3)                    \
    ENUM_ITEM_WITH_VALUE(VelocityVector, 4)                \
    ENUM_ITEM_WITH_VALUE(PressureVector, 5)                \
    ENUM_ITEM_WITH_VALUE(Temperature, 6)                   \
    ENUM_ITEM_WITH_VALUE(Humidity, 7)                      \
    ENUM_ITEM_WITH_VALUE(EigenVector, 8)                   \
    ENUM_ITEM_WITH_VALUE(DirectorField, 15) /* Vector field */ \
    ENUM_ITEM_WITH_VALUE(DeplanationFunction, 16)          \
    ENUM_ITEM_WITH_VALUE(MacroSlipVector, 17) \
    ENUM_ITEM_WITH_VALUE(ResidualForce, 18)
	/**
	 * Type representing particular unknown (its physical meaning).
	 */
	enum UnknownType {
		UnknownType_DEF
	};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h

	const char* __UnknownTypeToString(UnknownType _value);
} // end namespace fem
#endif // unknowntype_h
