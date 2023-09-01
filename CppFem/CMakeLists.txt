cmake_minimum_required(VERSION 3.0)
project(fem_test C CXX)

get_filename_component(ROOT_DIR ${CMAKE_CURRENT_LIST_FILE} DIRECTORY)

#Setting the outputs
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ROOT_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${ROOT_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ROOT_DIR}/lib)

#adding the sources manually

set(element_src 
	src/elements/Element.cpp
	src/elements/Quad2.cpp
	src/elements/Bar1d.cpp
	src/elements/Node.cpp
	src/elements/Rod.cpp
)

set(meshing_src
	src/meshing/Mesh.cpp
	src/meshing/NaiveMesh.cpp
	src/meshing/CustomMesh.cpp
)

set(interpolation_src
	src/interpolation/Interpolation.cpp
	src/interpolation/InterQuad2d.cpp
)

set(materials_src
	src/materials/ElasticMaterial.cpp
	src/materials/PlaneStress.cpp
	src/materials/PlaneStrain.cpp
	src/materials/Material.cpp
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

set(fem_src

	src/fem/Function.h
	src/fem/TestFunction.h
	src/fem/Discretization.cpp
	src/fem/Discretization.cpp
	src/fem/ImplAssembly.cpp
)

set(physics_src
	src/physics/BoundaryConditions.cpp
	src/physics/MechanicalBoundaryConditions.cpp
	src/physics/Force.cpp
)

set(SOURCES 
	${algebra_src}
	${meshing_src}				
	${element_src}
	${interpolation_src}
	${materials_src}
	${solver_src}
	${tensor_src}
	${fem_src}
	${physics_src}			
	#src/plotter/gnuplot_i.cpp
	#src/plotter/Plotter.cpp											
)

set(element_inc 
	src/elements/Element.h
	src/elements/Quad2.h
	src/elements/Bar1d.h
	src/elements/Node.h
	src/elements/Rod.h
)

set(algebra_inc
	src/algebra/Vector.h
	src/algebra/Matrix.h
)

set(meshing_inc
	src/meshing/Mesh.h
	src/meshing/NaiveMesh.h
	src/meshing/CustomMesh.h
)

set(interpolation_inc
	src/interpolation/Interpolation.h
	src/interpolation/InterQuad2d.h
)

set(materials_inc
	src/materials/ElasticMaterial.h
	src/materials/PlaneStress.h
	src/materials/PlaneStrain.h
	src/materials/Material.h
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
)

set(tensor_inc

	src/tensors/Tensor.h
	src/tensors/Stress.h
)

set(fem_inc

	src/fem/Function.h
	src/fem/TestFunction.h
	src/fem/Discretization.h
	src/fem/Discretization.h
	src/fem/ImplAssembly.h
)

set(physics_inc
	src/physics/BoundaryConditions.h
	src/physics/MechanicalBoundaryConditions.h
	src/physics/Force.h
)

set(include_dir 
	src/algebra
	src/meshing
	src/interpolation
	src/elements
	src/materials
	src/solver
	src/tensors
	src/fem
	src/physics
)

set(INCLUDES
	${algebra_inc}
	${meshing_inc}				
	${element_inc}
	${interpolation_inc}
	${materials_inc}
	${solver_inc}
	${tensor_inc}
	${fem_inc}
	${physics_inc}
)

#location of include directories
source_group(src/algebra FILES ${algebra_src})
source_group(src/meshing FILES ${meshing_src})
source_group(src/elements FILES ${element_src})
source_group(src/interpolation FILES ${interpolation_src})
source_group(src/materials FILES ${materials_src})
source_group(src/solver FILES ${solver_src})
source_group(src/tensor FILES ${tensor_src})
source_group(src/fem FILES ${fem_src})
source_group(src/physics FILES ${physics_src})

source_group(inc/algebra FILES ${algebra_inc})
source_group(inc/meshing FILES ${meshing_inc})
source_group(inc/elements FILES ${element_inc})
source_group(inc/interpolation FILES ${interpolation_inc})
source_group(inc/materials FILES ${materials_inc})
source_group(inc/solver FILES ${solver_inc})
source_group(inc/tensor FILES ${tensor_inc})
source_group(inc/fem FILES ${fem_inc})
source_group(inc/physics FILES ${physics_inc})


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