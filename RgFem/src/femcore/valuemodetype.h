/*****************************************************************//**
 * \file   valuemodetype.h
 * \brief
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef valuemodetype_h
#define valuemodetype_h

#include "enumitem.h"

namespace fem
{
	//
	// following mode determines the mode of particular unknown
	// which can be requested on DOF.
	// particular DOF contain for example displacement type unknown,
	// but we can request total value, increment of value or velocity or acceleration of
	// this unknown. This has been done mainly in  order to  improve runtime checking
	// ad Dof level.
	//
	// see also isUnknownModeIncrementalMode() function (cltypes.C)
	// when adding new ValueModeType mode.
	/*
	  Note:
	  VM_Total:   total value evaluated at the end of solution (time) step
	  VM_TotalIntrinsic: total value evaluated at intrinsic time
	*/


#define ValueModeType_DEF \
    ENUM_ITEM_WITH_VALUE(VM_Unknown, 0) \
    ENUM_ITEM_WITH_VALUE(VM_Total, 1)              \
    ENUM_ITEM_WITH_VALUE(VM_Velocity, 2)           \
    ENUM_ITEM_WITH_VALUE(VM_Acceleration, 3)       \
    ENUM_ITEM_WITH_VALUE(VM_Incremental, 4)        \
    ENUM_ITEM_WITH_VALUE(VM_RhsTotal, 5)           \
    ENUM_ITEM_WITH_VALUE(VM_RhsIncremental, 6)     \
    ENUM_ITEM_WITH_VALUE(VM_RhsInitial, 7)         \
    ENUM_ITEM_WITH_VALUE(VM_Intermediate, 8)       \
    ENUM_ITEM_WITH_VALUE(VM_TotalIntrinsic, 9)     \
    ENUM_ITEM_WITH_VALUE(VM_Residual, 99)



	/**
	 * Type representing the mode of UnknownType or CharType, or similar types.
	 * Afore mentioned types usually describes the physical meaning of
	 * value and ValueModeType provides the further necessary classification. For example "DisplacementVector"
	 * value can be further classified to be total displacement (TotalMode) or  velocity of
	 * displacement (VelocityMode) an so on.
	 */
	enum ValueModeType {
		ValueModeType_DEF
	};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


	const char* __ValueModeTypeToString(ValueModeType _value);
} // end namespace oofem
#endif // valuemodetype_h
