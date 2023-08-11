cmake_minimum_required(VERSION 3.0)
project(fem_test C CXX)

get_filename_component(ROOT_DIR ${CMAKE_CURRENT_LIST_FILE} DIRECTORY)

#Setting the outputs
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ROOT_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${ROOT_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ROOT_DIR}/lib)

#adding the sources manually
set(SOURCES src/algebra/Vector.h
			src/algebra/Matrix.h
			
			src/meshing/Mesh.cpp
			src/meshing/NaiveMesh.cpp
			src/meshing/CustomMesh.cpp

			src/interpolation/Interpolation.cpp
			src/interpolation/InterQuad2d.cpp
			
			src/elements/Element.cpp
			src/elements/Quad2.cpp
			src/elements/Node.cpp
			
			src/materials/ElasticMaterial.cpp
			src/materials/PlaneStress.cpp
			src/materials/PlaneStrain.cpp
			src/materials/Material.cpp
			src/materials/VonMises.cpp
			src/materials/VonMisesPlaneStress.cpp
			src/materials/VonMisesPlaneStrain.cpp
			src/materials/ElastoPlasticMaterial.cpp
			
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
			
			src/tensors/Tensor.cpp
			src/tensors/Stress.cpp
			
			src/fem/Function.h
			src/fem/TestFunction.h
			src/fem/GaussIntegration.cpp
			src/fem/Discretization.cpp
			src/fem/Discretization.cpp
			src/fem/ImplAssembly.cpp
			
			#src/plotter/gnuplot_i.cpp
			#src/plotter/Plotter.cpp
			
			src/physics/BoundaryConditions.cpp
			src/physics/MechanicalBoundaryConditions.cpp
			src/physics/Force.cpp						
)

set(INCLUDES src/algebra/Vector.h
			src/algebra/Matrix.h
			
			src/meshing/Mesh.h
			src/meshing/NaiveMesh.h
			src/meshing/CustomMesh.h

			src/interpolation/Interpolation.h
			src/interpolation/InterQuad2d.h
			
			src/elements/Element.h
			src/elements/Quad2.h
			src/elements/Node.h
			
			src/materials/ElasticMaterial.h
			src/materials/PlaneStress.h
			src/materials/PlaneStrain.h
			src/materials/Material.h
			src/materials/VonMises.h
			src/materials/VonMisesPlaneStress.h
			src/materials/VonMisesPlaneStrain.h
			src/materials/ElastoPlasticMaterial.h
			
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
			
			src/tensors/Tensor.h
			src/tensors/Stress.h
			
			src/fem/Function.h
			src/fem/TestFunction.h
			src/fem/GaussIntegration.h
			src/fem/Discretization.h
			src/fem/Discretization.h
			src/fem/ImplAssembly.h
			
			#src/plotter/gnuplot_i.h
			#src/plotter/Plotter.h
			
			src/physics/BoundaryConditions.h
			src/physics/MechanicalBoundaryConditions.h
			src/physics/Force.h						
)

#location of include directories

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