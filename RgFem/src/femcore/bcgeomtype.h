#ifndef bcgeomtype_h
#define bcgeomtype_h

namespace fem
{
	/// Type representing the geometric character of loading.
	enum bcGeomType {
		UnknownBGT,     ///< Unknown type.
		NodalLoadBGT,   ///< Concentrated nodal load.
		BodyLoadBGT,    ///< Distributed body load.
		EdgeLoadBGT,    ///< Distributed edge load.
		SurfaceLoadBGT, ///< Distributed surface load.
		PointLoadBGT,   ///< Concentrated point load (placed anywhere).
		GravityPressureBGT, ///<Pressure due to distributed body load.
	};
} // end namespace fem
#endif // bcgeomtype_h
