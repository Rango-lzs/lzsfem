#include "fecore_type.h"
#include <assert.h>

int fecore_data_size(FEDataType type)
{
	switch (type)
	{
	case FE_INVALID_TYPE: return 0; break;
	case FE_DOUBLE: return fecoreType<double>::size(); break;
	case FE_VEC2D : return fecoreType<Vector2d >::size(); break;
	case FE_Vector3d : return fecoreType<Vector3d >::size(); break;
	case FE_MAT3D : return fecoreType<Matrix3d >::size(); break;
	case FE_MAT3DS: return fecoreType<Matrix3ds>::size(); break;
	default:
		assert(false);
	}

	return 0;
};
