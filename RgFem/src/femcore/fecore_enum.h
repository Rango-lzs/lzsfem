/*********************************************************************
 * \file   fecore_enum.h
 * \brief  
 * 
 * \author Leizs
 * \date   February 2025
 *********************************************************************/

#pragma once

//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//! This lists the super-class id's that can be used to register new classes
//! with the kernel. It effectively defines the base class that a class
//! is derived from.
enum SUPER_CLASS_ID {
	FEINVALID_ID,					// an invalid ID
	FEOBJECT_ID,					// derived from FEObjectBase (TODO: work in progress)
	FETASK_ID,                   	// derived from FECoreTask
	FESOLVER_ID,                 	// derived from FESolver
	FEMATERIAL_ID,               	// derived from FEMaterial
	FEMATERIALPROP_ID,				// derived from FEMaterialProperty
	FEDISCRETEMATERIAL_ID,			// derived from FEDiscreteMaterial
	FELOAD_ID,               	    // derived from FEModelLoad
	FENLCONSTRAINT_ID,           	// derived from FENLConstraint
	FEPLOTDATA_ID,               	// derived from FEPlotData
	FEANALYSIS_ID,               	// derived from FEAnalysis
	FESURFACEINTERFACE_ID, 			// derived from FESurfaceInterface
	FELOGNODEDATA_ID,            	// derived from FELogNodeData
	FELOGFACEDATA_ID,            	// derived from FELogFaceData
	FELOGELEMDATA_ID,            	// derived from FELogElemData
	FELOGOBJECTDATA_ID,            	// derived from FELogObjectData
	FELOGDOMAINDATA_ID,            	// derived from FELogDomainData
	FELOGNLCONSTRAINTDATA_ID,      	// derived from FELogNLConstraintData
	FELOGSURFACEDATA_ID,      		// derived from FELogSurfaceData
	FELOGMODELDATA_ID,            	// derived from FEModelLogData
	FEBC_ID,						// derived from FEBoundaryCondition
	FEGLOBALDATA_ID,				// derived from FEGlobalData
	FECALLBACK_ID,					// derived from FECallBack
	FESOLIDDOMAIN_ID,				// derived from FESolidDomain
	FESHELLDOMAIN_ID,				// derived from FEShellDomain
	FEBEAMDOMAIN_ID,				// derived from FEBeamDomain
	FEDISCRETEDOMAIN_ID,			// derived from FEDiscreteDomain
	FEDOMAIN2D_ID,					// derived from FEDomain2D
	FESURFACE_ID,					// derived from FESurface
	FEIC_ID,						// derived from FEInitialCondition
	FEMESHDATAGENERATOR_ID,			// derived from FEMeshDataGenerator
	FELOADCONTROLLER_ID,			// derived from FELoadContoller
	FEMODEL_ID,						// derived from FEModel (TODO: work in progress)
	FESCALARVALUATOR_ID,			// derived from FEScalarValuator
	FEVector3dVALUATOR_ID,				// derived from FEVectorValuator
	FEMAT3DVALUATOR_ID,				// derived from FEMAT3DValuator
	FEMAT3DSVALUATOR_ID,			// derived from FEMAT3DSValuator
	FEFUNCTION1D_ID,				// derived from FEFunction1D
	FELINEARSOLVER_ID,				// derived from LinearSolver
	FEMESHADAPTOR_ID,				// derived from FEMeshAdaptor
	FEMESHADAPTORCRITERION_ID,		// derived from FEMeshAdaptorCriterion
	FENEWTONSTRATEGY_ID,			// derived from FENewtonStrategy
	FETIMECONTROLLER_ID,			// derived from FETimeStepController
	FEEIGENSOLVER_ID,				// derived from EigenSolver
	FEDATARECORD_ID,				// derived from DataRecord
	FECLASS_ID,						// derived from FECoreClass
};

//-----------------------------------------------------------------------------
// Plot level sets the frequency of writes to the plot file.
enum FE_Plot_Level {
	FE_PLOT_NEVER,			// don't output anything
	FE_PLOT_MAJOR_ITRS,		// only output major iterations (i.e. converged time steps)
	FE_PLOT_MINOR_ITRS,		// output minor iterations (i.e. every Newton iteration)
	FE_PLOT_MUST_POINTS,	// output only on must-points
	FE_PLOT_FINAL,			// only output final converged state
	FE_PLOT_AUGMENTATIONS,	// plot state before augmentations
	FE_PLOT_STEP_FINAL,		// output the final step of a step
	FE_PLOT_USER1			// plot will only happen on CB_USER1 callback
};

//-----------------------------------------------------------------------------
// Plot hint
enum FE_Plot_Hint {
	FE_PLOT_NO_HINT = 0,
	FE_PLOT_APPEND = 1		// don't close plot file after run
};

//-----------------------------------------------------------------------------
// Output level sets the frequency of data output is written to the log or data files.
enum FE_Output_Level {
	FE_OUTPUT_NEVER,
	FE_OUTPUT_MAJOR_ITRS,
	FE_OUTPUT_MINOR_ITRS,
	FE_OUTPUT_MUST_POINTS,
	FE_OUTPUT_FINAL
};

//-----------------------------------------------------------------------------
//! Domain classes
//! The domain class defines the general catergory of element types
#define	FE_DOMAIN_SOLID		1
#define	FE_DOMAIN_SHELL		2
#define	FE_DOMAIN_BEAM		3
#define	FE_DOMAIN_SURFACE	4
#define	FE_DOMAIN_DISCRETE	5
#define	FE_DOMAIN_2D		6
#define FE_DOMAIN_EDGE		7

// --- data types ---
enum Var_Type { 
	PLT_FLOAT,		// scalar             : single fp
	PLT_VEC3F,		// 3D vector          : 3 fps
	PLT_MAT3FS,		// symm 2o tensor     : 6 fps
	PLT_MAT3FD,		// diagonal 2o tensor : 3 fps
	PLT_TENS4FS,	// symm 4o tensor     : 21 fps
	PLT_MAT3F,		// 2o tensor          : 9 fps
	PLT_ARRAY,		// variable array (see dictionary for size)
	PLT_ARRAY_VEC3F	// array of vec3f (see dictionary for size)
};

// --- storage format ---
// FMT_NODE : one value stored for each node of a region
// FMT_ITEM : one value stored for each item (e.g. element) of a region
// FMT_MULT : one value for each node of each item of a region
// FMT_REGION: one value per region (surface, domain)
enum Storage_Fmt { FMT_NODE, FMT_ITEM, FMT_MULT, FMT_REGION, FMT_MATPOINTS };

//-----------------------------------------------------------------------------
enum FEDataType {
	FE_INVALID_TYPE,
	FE_DOUBLE,
	FE_VEC2D,
	FE_VEC3D,
	FE_MAT3D,
	FE_MAT3DS
};


//-----------------------------------------------------------------------------
//! Different Matrix types. This is used when requesting a sparse Matrix format
//! from a linear solver. 
//! \sa LinearSolver::CreateSparseMatrix.
enum MatrixType {
	REAL_UNSYMMETRIC,			// non-symmetric 
	REAL_SYMMETRIC,				// symmetric (not necessarily positive definite)
	REAL_SYMM_STRUCTURE			// structurally symmetric
};
