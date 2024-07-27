cmake_minimum_required(VERSION 3.1)

project(libfemp
	LANGUAGES CXX
	VERSION 0.1.0
	)

#add_subdirectory(libs/io)
#add_subdirectory(libs/mesh)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(libfemp_INCLUDES
	../include/libfemp/Analysis.hpp
	../include/libfemp/AnalysisResult.hpp
	../include/libfemp/Element.hpp
	../include/libfemp/ElementGroup.hpp
	../include/libfemp/ElementResults/ElementResults.hpp
	../include/libfemp/ElementResults/ElementResultsFactory.hpp
	../include/libfemp/ElementResults/ResultsRanges.hpp
	../include/libfemp/ElementResults/Strains.hpp
	../include/libfemp/ElementResults/Stresses.hpp
	../include/libfemp/elements/all.hpp
	../include/libfemp/elements/BaseElement.hpp
	../include/libfemp/elements/ElementFactory.hpp
	../include/libfemp/elements/Hexahedron20.hpp
	../include/libfemp/elements/Hexahedron27.hpp
	../include/libfemp/elements/Hexahedron8.hpp
	../include/libfemp/elements/HexahedronFamily.hpp
	../include/libfemp/elements/Prism15.hpp
	../include/libfemp/elements/Prism18.hpp
	../include/libfemp/elements/Prism6.hpp
	../include/libfemp/elements/PrismFamily.hpp
	../include/libfemp/elements/Tetrahedron10.hpp
	../include/libfemp/elements/Tetrahedron4.hpp
	../include/libfemp/elements/TetrahedronFamily.hpp
	../include/libfemp/Equation.hpp
	../include/libfemp/FemException.hpp
	../include/libfemp/GlobalDof.hpp
	../include/libfemp/GlobalDofMap.hpp
	../include/libfemp/Group.hpp
	../include/libfemp/io/export/JsonExporter.hpp
	../include/libfemp/io/import/ModelImporterFactory.hpp
	../include/libfemp/io/import/Parser.hpp
	../include/libfemp/LinearAnalysis.hpp
	../include/libfemp/LoadPattern.hpp
	../include/libfemp/loads/DomainLoad.hpp
	../include/libfemp/loads/NodalLoad.hpp
	../include/libfemp/loads/Quadrangle4.hpp
	../include/libfemp/loads/Quadrangle8.hpp
	../include/libfemp/loads/Quadrangle9.hpp
	../include/libfemp/loads/QuadrangleFamily.hpp
	../include/libfemp/loads/SurfaceLoad.hpp
	../include/libfemp/loads/Triangle10.hpp
	../include/libfemp/loads/Triangle3.hpp
	../include/libfemp/loads/Triangle6.hpp
	../include/libfemp/loads/TriangleFamily.hpp
	../include/libfemp/Material.hpp
	../include/libfemp/Model.hpp
	../include/libfemp/NodalDisplacement.hpp
	../include/libfemp/Node.hpp
	../include/libfemp/NodeGroup.hpp
	../include/libfemp/NodeRestrictions.hpp
	../include/libfemp/Point2D.hpp
	../include/libfemp/Point3D.hpp
	../include/libfemp/solvers/CGSolver.hpp
	../include/libfemp/solvers/CholeskySolver.hpp
	../include/libfemp/solvers/DenseCholeskySolver.hpp
	../include/libfemp/solvers/Solver.hpp
	../include/libfemp/solvers/UmfpackSolver.hpp
	../include/libfemp/Surface.hpp
	../include/libfemp/SurfaceLoadOperators/ConstantLoad.hpp
	../include/libfemp/SurfaceLoadOperators/SurfaceLoadOperator.hpp
	../include/libfemp/SurfaceLoadOperators/SurfaceNormalLoad.hpp
	)

set(libfemp_HEADERS
	io/import/FemJsonParser.hpp
	io/import/MshParser.hpp
	)

set(libfemp_SOURCES
	Analysis.cpp
	AnalysisResult.cpp
	Element.cpp
	ElementGroup.cpp
	ElementResults/ElementResults.cpp
	ElementResults/ElementResultsFactory.cpp
	elements/BaseElement.cpp
	elements/ElementFactory.cpp
	elements/Hexahedron20.cpp
	elements/Hexahedron27.cpp
	elements/Hexahedron8.cpp
	elements/HexahedronFamily.cpp
	elements/Prism15.cpp
	elements/Prism18.cpp
	elements/Prism6.cpp
	elements/PrismFamily.cpp
	elements/Tetrahedron10.cpp
	elements/Tetrahedron4.cpp
	elements/TetrahedronFamily.cpp
	Equation.cpp
	FemException.cpp
	GlobalDof.cpp
	GlobalDofMap.cpp
	Group.cpp
	io/export/JsonExporter.cpp
	io/import/FemJsonParser.cpp
	io/import/ModelImporterFactory.cpp
	io/import/MshParser.cpp
	io/import/Parser.cpp
	LinearAnalysis.cpp
	LoadPattern.cpp
	loads/DomainLoad.cpp
	loads/NodalLoad.cpp
	loads/Quadrangle4.cpp
	loads/Quadrangle8.cpp
	loads/Quadrangle9.cpp
	loads/QuadrangleFamily.cpp
	loads/SurfaceLoad.cpp
	loads/Triangle10.cpp
	loads/Triangle3.cpp
	loads/Triangle6.cpp
	loads/TriangleFamily.cpp
	Material.cpp
	Model.cpp
	NodalDisplacement.cpp
	Node.cpp
	NodeGroup.cpp
	NodeRestrictions.cpp
	Point2D.cpp
	Point3D.cpp
	Surface.cpp
	SurfaceLoadOperators/ConstantLoad.cpp
	SurfaceLoadOperators/SurfaceNormalLoad.cpp
	)

add_library(libfemp STATIC
	${libfemp_INCLUDES}
	${libfemp_HEADERS}
	${libfemp_SOURCES}
	)

target_link_libraries(libfemp
	PUBLIC
	Boost::headers
	#Eigen3::Eigen
	libla
	libquadrature::Object
	PRIVATE
	#libfemp_io::Object
	#libfemp_mesh::Object
	)

set_target_properties(libfemp PROPERTIES
	CXX_STANDARD 14
	CXX_EXTENSIONS OFF
	)

target_include_directories(libfemp
	PUBLIC 
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
	$<INSTALL_INTERFACE:include>
	)

# 添加头文件目录
set(RapidJson_Root "D:/Lzs/PhysicEngine/FEM_Project/femp")
include_directories(${RapidJson_Root})

add_library(libfemp::Static ALIAS libfemp)
