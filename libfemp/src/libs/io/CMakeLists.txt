cmake_minimum_required(VERSION 3.10)

set(libfemp_io_INCLUDES
	include/libfemp/io/import/Parser.hpp
	)

set(libfemp_io_HEADERS
    src/export/GmshExporter.hpp
    )

set(libfemp_io_SOURCES
    src/export/GmshExporter.cpp
    src/import/Parser.cpp
    )

add_library(libfemp_io OBJECT
	${libfemp_io_INCLUDES}
	${libfemp_io_HEADERS}
	${libfemp_io_SOURCES}
	)

target_link_libraries(libfemp_io
	PRIVATE
    libfemp_mesh::Object
    )

set_target_properties(libfemp_io PROPERTIES
	CXX_STANDARD 14
	CXX_EXTENSIONS OFF
	)

target_include_directories(libfemp_io
	PUBLIC 
	${CMAKE_CURRENT_SOURCE_DIR}/include
	PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/src
	)

add_library(libfemp_io::Object ALIAS libfemp_io)
