#ifndef matesponsemode_h
#define matesponsemode_h

#include "enumitem.h"

namespace fem
{
#define MatResponseMode_DEF \
    ENUM_ITEM(TangentStiffness) \
    ENUM_ITEM(SecantStiffness) \
    ENUM_ITEM(ElasticStiffness)                                           \
    ENUM_ITEM(Conductivity)  /* element level conductivity matrix */ \
    ENUM_ITEM(Conductivity_ww) /* material level conductivity submatrix */ \
    ENUM_ITEM(Conductivity_hh) /* material level conductivity submatrix */ \
    ENUM_ITEM(Conductivity_hw) /* material level conductivity submatrix */ \
    ENUM_ITEM(Conductivity_wh) /* material level conductivity submatrix */ \
    ENUM_ITEM(Capacity)                                                   \
    ENUM_ITEM(Capacity_ww) /* material level capacity submatrix */ \
    ENUM_ITEM(Capacity_hh) /* material level capacity submatrix */ \
    ENUM_ITEM(Capacity_hw) /* material level capacity submatrix */ \
    ENUM_ITEM(Capacity_wh) /* material level capacity submatrix */ \
    ENUM_ITEM(IntSource)                                                  \
    ENUM_ITEM(IntSource_ww) /* material level internal source submatrix - water source */ \
    ENUM_ITEM(IntSource_hh) /*  - heat source */ \
    ENUM_ITEM(IntSource_hw) /*  - heat source dependency on water content change */ \
    ENUM_ITEM(IntSource_wh) /*  - water source dependency on temperature change */

	/**
	 * Describes the character of characteristic material matrix.
	 */
	enum MatResponseMode {
		MatResponseMode_DEF
	};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


	const char* __MatResponseModeToString(MatResponseMode _value);
} // end namespace fem
#endif // matesponsemode_h
