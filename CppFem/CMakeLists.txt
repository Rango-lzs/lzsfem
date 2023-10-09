cmake_minimum_required(VERSION 3.0)
project(fem_test C CXX)

get_filename_component(ROOT_DIR ${CMAKE_CURRENT_LIST_FILE} DIRECTORY)

#Setting the outputs
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ROOT_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${ROOT_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ROOT_DIR}/lib)

#adding the sources manually

set(algebar_src
	src/algebra/floatmatrix.cpp
	src/algebra/floatarray.cpp
	src/algebra/dictionary.cpp
	src/algebra/scalarfunction.cpp
	src/algebra/mathfem.cpp
)

set(element_src 
	src/elements/element.cpp
	src/elements/Quad2.cpp
	src/elements/Bar1d.cpp
	#src/elements/Node.cpp
)

set(node_src
	src/node/node.cpp
)

set(meshing_src
	src/meshing/Mesh.cpp
	src/meshing/NaiveMesh.cpp
	src/meshing/CustomMesh.cpp
)

set(input_src
	src/input/inputrecord.cpp
	src/input/text_input_record.cpp
	src/input/dynamic_input_record.cpp
	src/input/text_data_reader.cpp
	src/input/dynamic_data_reader.cpp
	src/input/tokenizer.cpp
)

set(interpolation_src
	src/interpolation/Interpolation.cpp
	src/interpolation/InterQuad2d.cpp
	src/interpolation/fe_interpol.cpp
)

set(materials_src
	src/materials/ElasticMaterial.cpp
	src/materials/PlaneStress.cpp
	src/materials/PlaneStrain.cpp
	src/materials/material.cpp
	src/materials/VonMises.cpp
	src/materials/VonMisesPlaneStress.cpp
	src/materials/VonMisesPlaneStrain.cpp
	src/materials/ElastoPlasticMaterial.cpp
)

set(solver_src
	src/solver/Jacobi.cpp
	src/solver/LinearIterativeSolver.cpp
	src/solver/descent/DescentMethod.cpp
	src/solver/descent/GradientDescent.cpp
	src/solver/descent/ConjugateGradientDescent.cpp
	src/solver/descent/BiCGStab.cpp
	src/solver/NLSolver.cpp
	src/solver/LSolver.cpp
	src/solver/NLSolverCG.cpp
	src/solver/NLSolverBiCGStab.cpp
	src/solver/ArcLength.cpp
)

set(tensor_src

	src/tensors/Tensor.cpp
	src/tensors/Stress.cpp
)

set(femcore_src

	src/femcore/Function.h
	src/femcore/TestFunction.h
	src/femcore/Discretization.cpp
	src/femcore/Discretization.cpp
	src/femcore/ImplAssembly.cpp
	src/femcore/domain.cpp
	src/femcore/timer.cpp
	src/femcore/timestep.cpp
	src/femcore/fem_factory.cpp
	src/femcore/fem_application.cpp

	src/femcore/field.cpp
	src/femcore/field_manager.cpp
	src/femcore/engngm.cpp
	src/femcore/femcmpnn.cpp
)

set(physics_src
	src/physics/BoundaryConditions.cpp
	src/physics/MechanicalBoundaryConditions.cpp
	src/physics/Force.cpp
)

set(logger_src
	src/logger/logger.cpp
	src/logger/error.cpp
)

set(quadrature_src
	src/quadrature/integrationrule.cpp
)

set(SOURCES 
	${algebra_src}
	${meshing_src}				
	${element_src}
	${input_src}
	${interpolation_src}
	${materials_src}
	${solver_src}
	${tensor_src}
	${femcore_src}
	${physics_src}
	${logger_src}
	${node_src}
	${quadrature_src}
	#src/plotter/gnuplot_i.cpp
	#src/plotter/Plotter.cpp											
)

set(element_inc 
	src/elements/element.h
	src/elements/Quad2.h
	src/elements/Bar1d.h
	#src/elements/Node.h
)

set(node_inc
	src/node/node.h
)

set(algebra_inc
	src/algebra/Vector.h
	src/algebra/Matrix.h
	src/algebra/floatmatrix.h
	src/algebra/floatmatrixf.h
	src/algebra/floatarray.h
	src/algebra/floatarrayf.h
	src/algebra/intarray.h
	src/algebra/intarrayf.h
	src/algebra/dictionary.h
	src/algebra/range.h
	src/algebra/pair.h
	src/algebra/scalarfunction.h
	src/algebra/sparsemtrxtype.h
	src/algebra/mathfem.h
)

set(meshing_inc
	src/meshing/Mesh.h
	src/meshing/NaiveMesh.h
	src/meshing/CustomMesh.h
)

set(input_inc
	src/input/datareader.h
	src/input/inputrecord.h
	src/input/text_input_record.h
	src/input/dynamic_input_record.h
	src/input/text_data_reader.h
	src/input/dynamic_data_reader.h
	src/input/tokenizer.h
)

set(interpolation_inc
	src/interpolation/Interpolation.h
	src/interpolation/InterQuad2d.h
	src/interpolation/fe_interpol.h
)

set(materials_inc
	src/materials/ElasticMaterial.h
	src/materials/PlaneStress.h
	src/materials/PlaneStrain.h
	src/materials/material.h
	src/materials/VonMises.h
	src/materials/VonMisesPlaneStress.h
	src/materials/VonMisesPlaneStrain.h
	src/materials/ElastoPlasticMaterial.h
)

set(solver_inc
	src/solver/Jacobi.h
	src/solver/LinearIterativeSolver.h
	src/solver/descent/DescentMethod.h
	src/solver/descent/GradientDescent.h
	src/solver/descent/ConjugateGradientDescent.h
	src/solver/descent/BiCGStab.h
	src/solver/NLSolver.h
	src/solver/LSolver.h
	src/solver/NLSolverCG.h
	src/solver/NLSolverBiCGStab.h
	src/solver/ArcLength.h

	src/solver/linsystsolvertype.h
	src/solver/geneigvalsolvertype.h
)

set(tensor_inc

	src/tensors/Tensor.h
	src/tensors/Stress.h
)

set(femcore_inc

	src/femcore/Function.h
	src/femcore/TestFunction.h
	src/femcore/Discretization.h
	src/femcore/Discretization.h
	src/femcore/ImplAssembly.h
	src/femcore/domain.h
	src/femcore/timer.h
	src/femcore/timestep.h

	src/femcore/fem_export.h
	src/femcore/domaintype.h
	src/femcore/enumitem.h
	src/femcore/fem_factory.h
	src/femcore/fem_application.h

	src/femcore/contextmode.h
	src/femcore/contextoutputmode.h
	src/femcore/contextioresulttype.h

	src/femcore/field.h
	src/femcore/field_manager.h
	src/femcore/valuemodetype.h
	src/femcore/doftype.h
	src/femcore/dofiditem.h
	src/femcore/dof.h

	src/femcore/outputmanager.h
	src/femcore/engngm.h
	src/femcore/femcmpnn.h
	
)

set(physics_inc
	src/physics/BoundaryConditions.h
	src/physics/MechanicalBoundaryConditions.h
	src/physics/Force.h
)

set(logger_inc
	src/logger/logger.h
	src/logger/error.h
	src/logger/verbose.h
)

set(quadrature_inc
	src/quadrature/integrationrule.h
)

set(include_dir 
	src/algebra
	src/meshing
	src/interpolation
	src/elements
	src/materials
	src/solver
	src/tensors
	src/femcore
	src/physics
	src/input
	src/output
	src/logger
	src/node
	src/quadrature
)

set(INCLUDES
	${algebra_inc}
	${meshing_inc}				
	${element_inc}
	${input_inc}
	${interpolation_inc}
	${materials_inc}
	${solver_inc}
	${tensor_inc}
	${femcore_inc}
	${physics_inc}
	${logger_inc}
	${node_inc}
	${quadrature_inc}
)

#location of include directories
source_group(src/algebra FILES ${algebra_src})
source_group(src/meshing FILES ${meshing_src})
source_group(src/elements FILES ${element_src})
source_group(src/input  FILES ${input_src})
source_group(src/interpolation FILES ${interpolation_src})
source_group(src/materials FILES ${materials_src})
source_group(src/solver FILES ${solver_src})
source_group(src/tensor FILES ${tensor_src})
source_group(src/femcore FILES ${femcore_src})
source_group(src/physics FILES ${physics_src})

source_group(inc/algebra FILES ${algebra_inc})
source_group(inc/meshing FILES ${meshing_inc})
source_group(inc/elements FILES ${element_inc})
source_group(inc/input FILES ${input_inc})
source_group(inc/interpolation FILES ${interpolation_inc})
source_group(inc/materials FILES ${materials_inc})
source_group(inc/solver FILES ${solver_inc})
source_group(inc/tensor FILES ${tensor_inc})
source_group(inc/femcore FILES ${femcore_inc})
source_group(inc/physics FILES ${physics_inc})

add_definitions(-DFEM_MODULE)


include_directories(FEM ${include_dir})
add_library(FEM ${SOURCES} ${INCLUDES})

#creation of project
#add_executable(fem_test src/vonMisesTest.cpp)
#add_executable(fem_test src/plot_test.cpp)
#add_executable(fem_test src/2dSquare.cpp)
#add_executable(fem_test src/vector_test.cpp)
#add_executable(fem_test src/examples/vonMisesSquare.cpp)
#add_executable(fem_test src/examples/meshSets.cpp)
#add_executable(fem_test src/examples/nonLinear.cpp)
add_executable(fem_test src/examples/arcLengthTest.cpp src/examples/meshSets.cpp src/examples/nonLinear.cpp
src/examples/VonMisesSquare.cpp)
target_link_libraries(fem_test FEM)
target_include_directories(fem_test PUBLIC src/ src/plotter/)