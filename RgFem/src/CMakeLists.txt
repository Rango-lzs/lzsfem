set(TARGET_NAME RgFem)
project(${TARGET_NAME} LANGUAGES CXX)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(${TARGET_NAME})
target_compile_definitions(RgFem PRIVATE RgFem_Module)

target_include_directories(${TARGET_NAME} PRIVATE ${SOURCE_DIR})

#Add CLI11 Library, a head only comand parse library
set(CLI11_INCLUDE_DIR "${THIRD_PARTY_ROOT}/cli11-src/include")
target_include_directories(${TARGET_NAME} PRIVATE ${CLI11_INCLUDE_DIR})
message(status "Third party root is ****${THIRD_PARTY_ROOT}")

#Add Eigen3 Library
#set(Eigen_Root "D:/Program Files/eigen-3.4.0/eigen-3.4.0")
#target_include_directories(${TARGET_NAME} PRIVATE ${Eigen_Root})

#set(Eigen3_INCLUDE_DIR "D:/ThirdParty/eigen3") //at home

message(status "cmake module****${CMAKE_MODULE_PATH}")

message(status "cmake source****${CMAKE_SOURCE_DIR}")

# ģ·
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmakes")

# ҵEigen
find_package(Eigen3 REQUIRED)
 
# ҵEigen⣬Ŀ¼ӵĿ
if(Eigen3_FOUND)
  message(status "Eigen3 found****")
  include_directories(${Eigen3_INCLUDE_DIR})
endif()

set(PROJECT_SOURCES "")

# add app sources
#set(APP_SOURCES "")
#add_subdirectory(app)
#source_group("app" FILES ${APP_SOURCES})
#message(status "resource files is ${APP_SOURCES}")


# add basicio sources
#set(BASICIO_SOURCES "")
#add_subdirectory(basicio)
#source_group("basicio" FILES ${BASICIO_SOURCES})

#add app sources
set(APP_HEADERS
    app/AppUtils.h
    app/CmdOptions.h
    app/FEAppConfig.h
    app/RgFemApp.h
)
set(APP_SOURCES
    app/CmdOptions.cpp
    app/FEAppConfig.cpp
    app/AppUtils.cpp
    app/RgFemApp.cpp
    app/main.cpp
)


set(PROJECT_SOURCES
    ${PROJECT_SOURCES}
    ${APP_HEADERS}
    ${APP_SOURCES}
)
source_group("app" FILES ${APP_HEADERS} ${APP_SOURCES})


# add femcore sources
set(FEMCORE_HEADERS
    femcore/Callback.h
    femcore/FEModelObserver.h
    femcore/DOFS.h

    femcore/RgDofSchema.h
    femcore/RgDofSchema.cpp

    femcore/Equation/RgEquationNumbering.h
    femcore/Equation/RgEquationNumbering.cpp

    #femcore/DataRecord.h
    #femcore/Domain/DomainDataRecord.h
    #femcore/Domain/Rg3FieldElasticShellDomain.h
    #femcore/Domain/Rg3FieldElasticSolidDomain.h
    femcore/Domain/RgAssembler.h
    #femcore/Domain/RgBeamDomain.h
    femcore/Domain/RgDomain.h
    #femcore/Domain/RgDomain2D.h
    femcore/Domain/RgDomainList.h
    femcore/Domain/RgDomainParameter.h
    femcore/Domain/RgElasticAssembler.h
    #femcore/Domain/RgElasticDomain.h
    #femcore/Domain/RgElasticSolidDomain.h
    #femcore/Domain/RgElasticShellDomain.h
    #femcore/Domain/RgElasticTrussDomain.h
    #femcore/Domain/RgLinearTrussDomain.h
    femcore/Domain/RgShellDomain.h
    femcore/Domain/RgSolidDomain.h
    #femcore/Domain/RgSolidDomainFactory.h
    #femcore/Domain/RgTrussDomain.h
    #femcore/DomainDataRecord.h
    #femcore/ElementDataRecord.h
    #femcore/FE3FieldElasticSolidDomain.h

    #femcore/AnalysisStep/RgAnalysisStep.h
    #femcore/AnalysisStep/RgAnalysisStep.cpp
    #femcore/FEAnalysis/FEAnalysis.h
    #femcore/FEAnalysis/FEAnalysis.cpp

    femcore/AnalysisStep/RgAnalysis.h
    femcore/AnalysisStep/RgAnalysis.cpp

    femcore/AnalysisStep/AnalysisStep.h

    femcore/AnalysisStep/StaticStep.h
    femcore/AnalysisStep/StaticStep.cpp
    femcore/AnalysisStep/DynamicImplicitStep.cpp

    femcore/AnalysisStep/StepControl.h
    
    


    femcore/BoundaryCondition/RgBoundaryCondition.h
    femcore/BoundaryCondition/RgBoundaryCondition.cpp
    femcore/Load/RgLoad.h
    femcore/Load/RgLoad.cpp

    femcore/RgLoadController.h
    femcore/RgLoadController.cpp

    femcore/FEBodyLoad.h
    femcore/FEBoundaryCondition.h
    #femcore/FECore.h
    #femcore/FECoreClass.h
    #femcore/FECoreFactory.h
    #femcore/FECoreKernel.h
    femcore/FECoreTask.h
    femcore/FEDataMap.h
    femcore/FEDomainMap.h
    femcore/FEDiscreteSet.h
    femcore/FEDofList.h
    femcore/FEException.h
    femcore/FEFacetSet.h
    femcore/FEFunction1D.h
    femcore/FEGlobalData.h
    #femcore/FEGlobalMatrix.h
    femcore/FEGlobalVector.h
    femcore/FEInitialCondition.h
    femcore/FEItemList.h
    femcore/FELinearConstraintManager.h
    femcore/FELoadController.h
    femcore/FELoadCurve.h
    femcore/FEMat3dValuator.h
    femcore/FEMat3dsValuator.h
    femcore/FEMesh.h
    femcore/FEModel.h
    femcore/FEModelComponent.h
    #femcore/FEModelDataRecord.h
    femcore/FEModelLoad.h
    femcore/FEModelParam.h
    femcore/FEModule.h
    femcore/FENLConstraint.h
    femcore/FENode.h
    femcore/FENodeDataMap.h
    #femcore/FENodeElemList.h
    femcore/FENodeList.h
    femcore/FENodeNodeList.h
    femcore/FENodeSet.h
    femcore/FEObjectBase.h
    femcore/FEParam.h
    femcore/FEParamObject.h
    femcore/FEParamValidator.h
    femcore/FEParameterList.h
    #femcore/FEPlotDataStore.h
    femcore/FEPointFunction.h
    femcore/FEProperty.h
    femcore/FEPropertyT.h
    femcore/FERgModel.h
    femcore/FEScalarValuator.h
    #femcore/FESolidAnalysis.h
    femcore/FEStepComponent.h
    femcore/FESurfaceBC.h
    femcore/FESurfaceLoad.h
    femcore/FESurfaceMap.h
    femcore/FESurfacePair.h
    femcore/FESurfacePairConstraint.h
    femcore/FESurfacePairConstraintNL.h
    #femcore/FESurfaceToSurfaceMap.h
    #femcore/FETimeInfo.h
    femcore/FEValuator.h
    femcore/FEVec3dValuator.h
    femcore/FSPath.h
    #femcore/FaceDataRecord.h

    femcore/Matrix/CSRMatrix.h
    femcore/Matrix/CompactMatrix.h
    femcore/Matrix/CompactSymmMatrix.h
    femcore/Matrix/CompactUnSymmMatrix.h
    femcore/Matrix/DenseMatrix.h
    femcore/Matrix/FEGlobalMatrix.h
    femcore/Matrix/JFNKMatrix.h
    femcore/Matrix/MMatrix.h
    femcore/Matrix/MatrixOperator.h
    femcore/Matrix/MatrixProfile.h
    femcore/Matrix/SkylineMatrix.h
    femcore/Matrix/SparseMatrix.h  

    femcore/MatrixSolver/AccelerateSparseSolver.h
    femcore/MatrixSolver/BIPNSolver.h
    femcore/MatrixSolver/BiCGStabSolver.h
    femcore/MatrixSolver/BlockSolver.h
    femcore/MatrixSolver/BoomerAMGSolver.h
    femcore/MatrixSolver/FEASTEigenSolver.h
    femcore/MatrixSolver/FGMRESSolver.h
    femcore/MatrixSolver/HypreGMRESsolver.h
    femcore/MatrixSolver/MKLDSSolver.h
    femcore/MatrixSolver/PardisoProjectSolver.h
    femcore/MatrixSolver/PardisoSolver.h
    femcore/MatrixSolver/RCICGSolver.h
    femcore/MatrixSolver/SchurSolver.h
    femcore/MatrixSolver/StrategySolver.h
    femcore/MatrixSolver/TestSolver.h

    femcore/MetaObject.h
    femcore/MetaObjectDemo.h
    femcore/PointCurve.h
    femcore/RTTI/MetaClass.h
    femcore/RTTI/MetaClassStore.h
    femcore/RTTI/RTTIMacroDefine.h

    femcore/Solver/BFGSSolver.h
    #femcore/Solver/EigenSolver.h
    #femcore/Solver/FEBroydenStrategy.h
    #femcore/Solver/FECGSolidSolver.h
    #femcore/Solver/FEExplicitSolidSolver.h
    #femcore/Solver/FEFullNewtonStrategy.h
    #femcore/Solver/FEGlobalData.h
    femcore/Solver/FELineSearch.h
    femcore/Solver/FELinearSolver.h
    #femcore/Solver/FEModifiedNewtonStrategy.cpp
    #femcore/Solver/FENewtonSolver.cpp
    femcore/Solver/FENewtonStrategy.h
    femcore/Solver/FENodeReorder.h
    #femcore/Solver/FERigidSolver.h
    #femcore/Solver/FESolidSolver2.h
    #femcore/Solver/FESolver.h
    #femcore/Solver/JFNKStrategy.h
    femcore/Solver/LUSolver.h
    femcore/Solver/LinearSolver.h
    femcore/Solver/SkylineSolver.h
    femcore/Solver/FELevelStructure.h

    femcore/NewtonSolver/SolverBase.h
    femcore/NewtonSolver/SolverBase.cpp
    femcore/NewtonSolver/StaticSolver.h
    femcore/NewtonSolver/StaticSolver.cpp
    femcore/NewtonSolver/NewtonSolver.h
    femcore/NewtonSolver/NewtonSolver.cpp


    femcore/TimeStep/FETimeInfo.h
    femcore/TimeStep/FETimeStepController.h
    femcore/Timer.h
    femcore/fecore_api.h
    femcore/fecore_debug.h
    femcore/fecore_debug_t.h
    femcore/fecore_enum.h
    femcore/fecore_type.h
    femcore/fem_export.h
    femcore/sys.h
    femcore/FELinearSystem.h
    femcore/FESolidLinearSystem.h
    femcore/tools.h
    femcore/FEResidualVector.h
    femcore/units.h
    femcore/FELinearConstraint.h
    femcore/FENodalLoad.h
    femcore/FEDataArray.h
    femcore/FEPressureLoad.h
    femcore/FENodalBC.h
    femcore/FEFixedBC.h
    femcore/FEFixedDisplacement.h
    femcore/FEFixedRotation.h
    femcore/FETransform.h
    femcore/FEPrescribedDOF.h
    femcore/FEPrescribedBC.h
    femcore/FEBodyForce.h
    femcore/FEConstValueVec3.h
    femcore/MFunctions.h
    femcore/MTypes.h
    femcore/MMath.h
    femcore/MObj2String.h
    femcore/MItem.h
    femcore/MathObject.h
    femcore/MatrixTools.h
    femcore/BSpline.h
    femcore/eig3.h
    femcore/sys.h

    femcore/FEPrescribedDisplacement.h

    femcore/FESolidModule.h
    femcore/RgElementTraitsStore.h
    femcore/RgElementShapeStore.h
)

set(FEMCORE_SOURCES
    femcore/Callback.cpp
    femcore/DOFS.cpp
    #femcore/DataRecord.cpp

    #femcore/Domain/DomainDataRecord.cpp
    femcore/Domain/RgAssembler.cpp
    #femcore/Domain/RgBeamDomain.cpp
    femcore/Domain/RgDomain.cpp
    femcore/Domain/RgDomainList.cpp
    femcore/Domain/RgDomainParameter.cpp
    femcore/Domain/RgElasticAssembler.cpp
    #femcore/Domain/RgElasticDomain.cpp
    #femcore/Domain/RgElasticSolidDomain.cpp
    #femcore/Domain/RgElasticShellDomain.cpp
    #femcore/Domain/RgElasticTrussDomain.cpp
    #femcore/Domain/RgLinearTrussDomain.cpp
    femcore/Domain/RgShellDomain.cpp
    femcore/Domain/RgSolidDomain.cpp
    #femcore/Domain/RgSolidDomainFactory.cpp
    #femcore/Domain/RgTrussDomain.cpp

    #femcore/FEAnalysis/FEAnalysis.cpp
    femcore/FEBodyLoad.cpp
    femcore/FEBoundaryCondition.cpp
    #femcore/FECore.cpp
    #femcore/FECoreClass.cpp
    #femcore/FECoreFactory.cpp
    #femcore/FECoreKernel.cpp
    femcore/FECoreTask.cpp
    femcore/FEDataMap.cpp
    femcore/FEDomainMap.cpp
    femcore/FEDiscreteSet.cpp
    femcore/FEDofList.cpp
    femcore/FEException.cpp
    femcore/FEFacetSet.cpp
    femcore/FEFunction1D.cpp
    femcore/FEGlobalData.cpp
    #femcore/FEGlobalMatrix.cpp
    femcore/FEGlobalVector.cpp
    femcore/FEInitialCondition.cpp
    femcore/FEItemList.cpp
    femcore/FELinearConstraintManager.cpp
    femcore/FELoadController.cpp
    femcore/FELoadCurve.cpp
    femcore/FEMat3dValuator.cpp
    femcore/FEMat3dsValuator.cpp
    femcore/FEMesh.cpp
    femcore/FEModel.cpp
    #femcore/FEModelDataRecord.cpp
    femcore/FEModelLoad.cpp
    femcore/FEModelParam.cpp
    femcore/FEModule.cpp
    femcore/FENLConstraint.cpp
    femcore/FENode.cpp
    femcore/FENodeDataMap.cpp
    femcore/FENodeElemList.cpp
    femcore/FENodeList.cpp
    femcore/FENodeNodeList.cpp
    femcore/FENodeSet.cpp
    femcore/FEObjectBase.cpp
    femcore/FEParam.cpp
    femcore/FEParamValidator.cpp
    femcore/FEParameterList.cpp
    femcore/FEPointFunction.cpp
    femcore/FEProperty.cpp
    femcore/FERgModel.cpp
    femcore/FEScalarValuator.cpp
    #femcore/FESolidAnalysis.cpp
    femcore/FEStepComponent.cpp
    femcore/FESurface.cpp
    femcore/FESurfaceBC.cpp
    femcore/FESurfaceLoad.cpp
    femcore/FESurfaceMap.cpp
    femcore/FESurfacePair.cpp
    femcore/FESurfacePairConstraint.cpp
    femcore/FESurfacePairConstraintNL.cpp
    #femcore/FESurfaceToSurfaceMap.cpp
    #femcore/FETimeInfo.cpp
    femcore/FEVec3dValuator.cpp
    femcore/FSPath.cpp
    #femcore/FaceDataRecord.cpp
    femcore/Matrix/CSRMatrix.cpp
    femcore/Matrix/CompactMatrix.cpp
    femcore/Matrix/CompactSymmMatrix.cpp
    femcore/Matrix/CompactUnSymmMatrix.cpp
    femcore/Matrix/DenseMatrix.cpp
    femcore/Matrix/FEGlobalMatrix.cpp
    femcore/Matrix/JFNKMatrix.cpp
    femcore/Matrix/MMatrix.cpp
    femcore/Matrix/MatrixProfile.cpp
    femcore/Matrix/SkylineMatrix.cpp
    femcore/Matrix/SparseMatrix.cpp  
    

    femcore/MetaObject.cpp
    femcore/MetaObjectDemo.cpp
    femcore/PointCurve.cpp
    femcore/RTTI/MetaClass.cpp
    femcore/RTTI/MetaClassStore.cpp
    femcore/Solver/BFGSSolver.cpp
    #femcore/Solver/EigenSolver.cpp
    #femcore/Solver/FEBroydenStrategy.cpp
    #femcore/Solver/FECGSolidSolver.cpp
    #femcore/Solver/FEExplicitSolidSolver.cpp
    femcore/Solver/FEFullNewtonStrategy.cpp
    #femcore/Solver/FEGlobalData.cpp
    femcore/Solver/FELineSearch.cpp
    femcore/Solver/FELinearSolver.cpp
    #femcore/Solver/FEModifiedNewtonStrategy.cpp
    #femcore/Solver/FENewtonSolver.cpp
    femcore/Solver/FENewtonStrategy.cpp
    femcore/Solver/FENodeReorder.cpp

    #femcore/Solver/FERigidSolver.cpp
    #femcore/Solver/FESolidSolver2.cpp
    #femcore/Solver/FESolver.cpp
    #femcore/Solver/JFNKStrategy.cpp
    femcore/Solver/LUSolver.cpp
    femcore/Solver/LinearSolver.cpp
    femcore/Solver/SkylineSolver.cpp
    femcore/Solver/FELevelStructure.cpp

    femcore/TimeStep/FETimeInfo.cpp
    femcore/TimeStep/FETimeStepController.cpp
    femcore/Timer.cpp
    femcore/fecore_debug.cpp
    femcore/fecore_type.cpp
    femcore/FELinearSystem.cpp
    femcore/FESolidLinearSystem.cpp
    femcore/tools.cpp
    femcore/FEResidualVector.cpp
    femcore/FELinearConstraint.cpp
    femcore/FENodalLoad.cpp
    femcore/FEDataArray.cpp
    femcore/FEPressureLoad.cpp
    femcore/FENodalBC.cpp
    femcore/FEFixedBC.cpp
    femcore/FEFixedDisplacement.cpp
    femcore/FEFixedRotation.cpp
    femcore/FETransform.cpp
    femcore/FEPrescribedDOF.cpp
    femcore/FEPrescribedBC.cpp
    femcore/FEBodyForce.cpp
    femcore/FEConstValueVec3.cpp
    femcore/MathObject.cpp
    femcore/BSpline.cpp
    femcore/eig3.cpp

    femcore/FEModelComponent.cpp
    femcore/FEValuator.cpp
    femcore/MItem.cpp
    femcore/MFunctions.cpp

    femcore/FEPrescribedDisplacement.cpp
    femcore/FESolidModule.cpp

    femcore/colsol.cpp


    femcore/RgElementTraitsStore.cpp
    femcore/RgElementShapeStore.cpp
)


set(PROJECT_SOURCES
    ${PROJECT_SOURCES}
    ${FEMCORE_HEADERS}
    ${FEMCORE_SOURCES}
)
source_group("femcore" FILES ${FEMCORE_HEADERS} ${FEMCORE_SOURCES})

# add datastructure sources
set(DATA_STRUCTURE_HEADERS
    datastructure/FEBoundingBox.h  
    datastructure/MathUtils.h
    datastructure/Matrix2d.h
    datastructure/Matrix3d.h
    datastructure/Matrix3d.hpp
    datastructure/MatrixOperator.h
    datastructure/MatrixProfile.h
    datastructure/Vector2d.h
    datastructure/Vector3d.h
    datastructure/mathalg.h
    datastructure/Matrix.h
    datastructure/quatd.h
    datastructure/tens3d.h
    datastructure/tens3d.hpp
    datastructure/tens3dls.hpp
    datastructure/tens3drs.hpp
    datastructure/tens3ds.hpp
    datastructure/tens4d.h
    datastructure/tens4d.hpp
    datastructure/tens4dmm.hpp
    datastructure/tens4dms.hpp
    datastructure/tens4ds.hpp
    datastructure/tensor_base.h
    datastructure/vector_operator.h
    datastructure/table.h
    datastructure/Matrix3dSpecial.h
)

set(DATA_STRUCTURE_SOURCES
    datastructure/Matrix3d.cpp
    datastructure/mathalg.cpp
    datastructure/Matrix.cpp
    datastructure/quatd.cpp
    datastructure/tens3d.cpp
    datastructure/tens4d.cpp
    datastructure/vector_operator.cpp
    datastructure/Matrix3dSpecial.cpp
)

set(PROJECT_SOURCES
    ${PROJECT_SOURCES}
    ${DATA_STRUCTURE_HEADERS}
    ${DATA_STRUCTURE_SOURCES}
)
source_group("datastructure" FILES ${DATA_STRUCTURE_HEADERS} ${DATA_STRUCTURE_SOURCES})

# add element sources from subdirectory
add_subdirectory(elements)
set(PROJECT_SOURCES
    ${PROJECT_SOURCES}
    ${ELEMENT_HEADERS}
    ${ELEMENT_SOURCES}

)
source_group("elements" FILES ${ELEMENT_HEADERS} ${ELEMENT_SOURCES})



set(BASICIO_HEADERS
    basicio/DataRecord.h
    basicio/DataStore.h
    basicio/DomainDataRecord.h
    basicio/DumpFile.h
    basicio/DumpMemStream.h
    basicio/DumpStream.h
    basicio/ElementDataRecord.h
    basicio/FEModelDataRecord.h
    basicio/FEPlotData.h
    basicio/FEPlotDataStore.h
    basicio/FaceDataRecord.h
    basicio/NLConstraintDataRecord.h
    basicio/NodeDataRecord.h
    basicio/ObjectDataRecord.h
    basicio/SurfaceDataRecord.h
    basicio/FEDataExport.h
    #basicio/FEBioMechPlot.h
)

set(BASICIO_SOURCES
    basicio/DataRecord.cpp
    basicio/DataStore.cpp
    basicio/DomainDataRecord.cpp
    basicio/DumpFile.cpp
    basicio/DumpMemStream.cpp
    basicio/DumpStream.cpp
    basicio/ElementDataRecord.cpp
    basicio/FEModelDataRecord.cpp
    basicio/FEPlotData.cpp
    basicio/FEPlotDataStore.cpp
    basicio/FaceDataRecord.cpp
    basicio/NLConstraintDataRecord.cpp
    basicio/NodeDataRecord.cpp
    basicio/ObjectDataRecord.cpp
    basicio/SurfaceDataRecord.cpp
    basicio/FEDataExport.cpp
    #basicio/FEBioMechPlot.cpp
)
set(PROJECT_SOURCES
    ${PROJECT_SOURCES}
    ${BASICIO_HEADERS}
    ${BASICIO_SOURCES}
)
source_group("basicio" FILES ${BASICIO_HEADERS} ${BASICIO_SOURCES})


#add input sources
set(INPUT_HEADERS
    input/FileReader.h
    input/XML/XMLReader.h
    input/XML/XMLWriter.h
    #input/abaqus/AbaqusExport.h
    input/abaqus/AbaqusImport.h
    input/abaqus/AbaqusBoundaryParser.h
    input/abaqus/AbaqusLoadParser.h
    input/abaqus/AbaqusMaterialConverter.h
    #input/abaqus/AbaqusModel.h
    #input/abaqus/AbaqusParser.h
    #input/febio/FEBModel.h
    #input/febio/FEBioBoundarySection.h
    #input/febio/FEBioBoundarySection3.h
    #input/febio/FEBioCodeSection.h
    #input/febio/FEBioConstraintsSection.h
    #input/febio/FEBioContactSection.h
    #input/febio/FEBioControlSection.h
    #input/febio/FEBioControlSection4.h
    #input/febio/FEBioDiscreteSection.h
    #input/febio/FEBioGeometrySection.h
    #input/febio/FEBioGlobalsSection.h
    #input/febio/FEBioImport.h
    #input/febio/FEBioIncludeSection.h
    #input/febio/FEBioInitialSection.h
    #input/febio/FEBioInitialSection3.h
    #input/febio/FEBioLoadDataSection.h
    #input/febio/FEBioLoadsSection.h
    #input/febio/FEBioMaterialSection.h
    #input/febio/FEBioMeshAdaptorSection.h
    #input/febio/FEBioMeshDataSection.h
    #input/febio/FEBioMeshDomainsSection4.h
    #input/febio/FEBioMeshSection.h
    #input/febio/FEBioMeshSection4.h
    #input/febio/FEBioModelBuilder.h
    #input/febio/FEBioModuleSection.h
    #input/febio/FEBioOutputSection.h
    #input/febio/FEBioRigidSection.h
    #input/febio/FEBioRigidSection4.h
    #input/febio/FEBioStepSection.h
    #input/febio/FEBioStepSection4.h
    #input/febio/FEModelBuilder.h
    #input/febio/FEPlotDataStore.h
    #input/febio/FERestartImport.h
    #input/febio/FileImport.h
    #input/febio/febioxml_api.h
    #input/febio/stdafx.h
    #input/febio/xmltool.h
   #input/oofem/abaqus_input_parser.h
   #input/oofem/datareader.h
   #input/oofem/dynamic_data_reader.h
   #input/oofem/dynamic_input_record.h
   #input/oofem/inputrecord.h
   #input/oofem/text_data_reader.h
   #input/oofem/text_input_record.h
   #input/oofem/tokenizer.h
)

set(INPUT_SOURCES
    input/XML/XMLReader.cpp
    input/XML/XMLWriter.cpp
    #input/abaqus/AbaqusExport.cpp
    input/FileReader.cpp
    input/abaqus/AbaqusImport.cpp
    input/abaqus/AbaqusBoundaryParser.cpp
    input/abaqus/AbaqusLoadParser.cpp
    input/abaqus/AbaqusMaterialConverter.cpp
    #input/abaqus/AbaqusModel.cpp
    #input/febio/FEBModel.cpp
    #input/febio/FEBioBoundarySection.cpp
    #input/febio/FEBioBoundarySection3.cpp
    #input/febio/FEBioCodeSection.cpp
    #input/febio/FEBioConstraintsSection.cpp
    #input/febio/FEBioContactSection.cpp
    #input/febio/FEBioControlSection.cpp
    #input/febio/FEBioControlSection4.cpp
    #input/febio/FEBioDiscreteSection.cpp
    #input/febio/FEBioGeometrySection.cpp
    #input/febio/FEBioGlobalsSection.cpp
    #input/febio/FEBioImport.cpp
    #input/febio/FEBioIncludeSection.cpp
    #input/febio/FEBioInitialSection.cpp
    #input/febio/FEBioInitialSection3.cpp
    #input/febio/FEBioLoadDataSection3.cpp
    #input/febio/FEBioLoadDataSection.cpp
    #input/febio/FEBioLoadsSection.cpp
    #input/febio/FEBioMaterialSection.cpp
    #input/febio/FEBioMeshAdaptorSection.cpp
    #input/febio/FEBioMeshDataSection.cpp
    #input/febio/FEBioMeshDataSection4.cpp
    #input/febio/FEBioMeshDomainsSection4.cpp
    #input/febio/FEBioMeshSection.cpp
    #input/febio/FEBioMeshSection4.cpp
    #input/febio/FEBioModelBuilder.cpp
    #input/febio/FEBioModuleSection.cpp
    #input/febio/FEBioOutputSection.cpp
    #input/febio/FEBioRigidSection.cpp
    #input/febio/FEBioRigidSection4.cpp
    #input/febio/FEBioStepSection.cpp
    #input/febio/FEBioStepSection4.cpp
    #input/febio/FEModelBuilder.cpp
    #input/febio/FERestartImport.cpp
    #input/febio/FileImport.cpp
    #input/febio/xmltool.cpp
    #input/oofem/dynamic_data_reader.cpp
    #input/oofem/dynamic_input_record.cpp
    #input/oofem/inputrecord.cpp
    #input/oofem/text_data_reader.cpp
    #input/oofem/text_input_record.cpp
    #input/oofem/tokenizer.cpp
)

set(PROJECT_SOURCES
    ${PROJECT_SOURCES}
    ${INPUT_HEADERS}
    ${INPUT_SOURCES}
)
source_group("input" FILES ${INPUT_HEADERS} ${INPUT_SOURCES})


#add logger sources
set(LOGGER_HEADERS
    logger/LogFileStream.h
    logger/LogStream.h
    logger/Logfile.h
    #logger/error.h
    logger/log.h
    #logger/logger.h
    logger/verbose.h
)

set(LOGGER_SOURCES
    logger/LogFileStream.cpp
    logger/LogStream.cpp
    logger/Logfile.cpp
    #logger/error.cpp
    logger/log.cpp
    #logger/logger.cpp
)
set(PROJECT_SOURCES
    ${PROJECT_SOURCES}
    ${LOGGER_HEADERS}
    ${LOGGER_SOURCES}
)
source_group("logger" FILES ${LOGGER_HEADERS} ${LOGGER_SOURCES})


#add material sources
add_subdirectory(materials)
set(PROJECT_SOURCES
    ${PROJECT_SOURCES}
    ${MATERIAL_HEADERS}
    ${MATERIAL_SOURCES}
)


source_group("materials\\LargeDeformation" FILES ${LARGE_DEF_HEADERS} ${LARGE_DEF_SOURCES})
source_group("materials\\SmallDeformation" FILES ${SMALL_DEF_HEADERS} ${SMALL_DEF_SOURCES})
source_group("materials" FILES ${MATERIAL_HEADERS} ${MATERIAL_SOURCES})


#add plotfile sources
set(PLOT_FILE_HEADERS
    plotfile/FEBioPlotFile.h
    plotfile/PlotFile.h
    plotfile/PltArchive.h
    plotfile/FEVTKPlotFile.h
    plotfile/VtuWriter.h
)

set(PLOT_FILE_SOURCES
    plotfile/FEBioPlotFile.cpp
    plotfile/PlotFile.cpp
    plotfile/PltArchive.cpp
    plotfile/FEVTKPlotFile.cpp
    plotfile/VtuWriter.cpp
)
set(PROJECT_SOURCES
    ${PROJECT_SOURCES}
    ${PLOT_FILE_HEADERS}
    ${PLOT_FILE_SOURCES}
)
source_group("plotfile" FILES ${PLOT_FILE_HEADERS} ${PLOT_FILE_SOURCES})

set(VTU_HEADERS
	../../ThirdParty/paraviewo-src/src/paraviewo/ParaviewWriter.hpp
	../../ThirdParty/paraviewo-src/src/paraviewo/VTMWriter.hpp
	../../ThirdParty/paraviewo-src/src/paraviewo/HDF5VTUWriter.hpp
	../../ThirdParty/paraviewo-src/src/paraviewo/VTUWriter.hpp
	../../ThirdParty/paraviewo-src/src/paraviewo/PVDWriter.hpp
	../../ThirdParty/paraviewo-src/src/paraviewo/base64Layer.hpp
	../../ThirdParty/paraviewo-src/src/paraviewo/base64Layer.cpp
)

set(VTU_SOURCES
	../../ThirdParty/paraviewo-src/src/paraviewo/VTMWriter.cpp
	../../ThirdParty/paraviewo-src/src/paraviewo/HDF5VTUWriter.cpp
	../../ThirdParty/paraviewo-src/src/paraviewo/VTUWriter.cpp
	../../ThirdParty/paraviewo-src/src/paraviewo/PVDWriter.cpp
	../../ThirdParty/paraviewo-src/src/paraviewo/base64Layer.cpp
)

set(PROJECT_SOURCES
    ${PROJECT_SOURCES}
)
source_group("vtuwriter" FILES ${VTU_HEADERS} ${VTU_SOURCES})


target_sources(${TARGET_NAME} PRIVATE ${PROJECT_SOURCES})

#file(GLOB rootfiles *.cpp *.h)
#source_group("" FILES ${rootfiles})
#target_sources(${TARGET_NAME} PRIVATE ${rootfiles})
#
#sub_dir_list(SUB_DIRS ${SOURCE_DIR})
#
#foreach(sub_dir ${SUB_DIRS})
#    file(GLOB subdir_files ${sub_dir}/*.cpp ${sub_dir}/*.h)
#    string(REPLACE "/" "\\" sub_dir ${sub_dir})
#    source_group(${sub_dir} FILES ${subdir_files})
#    target_sources(${TARGET_NAME} PRIVATE ${subdir_files})
#endforeach()