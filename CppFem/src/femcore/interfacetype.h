
#ifndef interfacetype_h
#define interfacetype_h

namespace fem
{
	/**
	 * Enumerative type, used to identify interface type.
	 * @see Interface More details.
	 */
	enum InterfaceType {
		UnknownInterfaceType,

		LayeredCrossSectionInterfaceType,
		FiberedCrossSectionInterfaceType,

		ZZNodalRecoveryModelInterfaceType,
		NodalAveragingRecoveryModelInterfaceType,
		SPRNodalRecoveryModelInterfaceType,

		ZZErrorEstimatorInterfaceType,
		HuertaErrorEstimatorInterfaceType,
		Huerta1dErrorEstimatorInterfaceType, // experimental

		SpatialLocalizerInterfaceType,

		EIPrimaryUnknownMapperInterfaceType,
		EIPrimaryFieldInterfaceType,

		NonlocalMaterialStatusExtensionInterfaceType,
		GradientDamageMaterialExtensionInterfaceType,
		GradientDamageMaterialStatusExtensionInterfaceType,

		NonlocalMaterialExtensionInterfaceType,
		NonlocalMaterialStiffnessInterfaceType,
		MaterialModelMapperInterfaceType,
		RandomMaterialStatusExtensionInterfaceType,

		HydrationModelInterfaceType,
		HydrationModelStatusInterfaceType,

		LEPlicElementInterfaceType,
		LevelSetPCSElementInterfaceType,

		XfemElementInterfaceType,
		VTKXMLExportModuleElementInterfaceType,
		FailureModuleElementInterfaceType,

		Beam3dSubsoilElementInterfaceType,
		Beam3dSubsoilMaterialInterfaceType,

		QCMaterialExtensionInterfaceType,

		MixedPressureMaterialExtensionInterfaceType
	};
} // end namespace fem
#endif // interfacetype_h
