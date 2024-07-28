cmake_minimum_required(VERSION 3.1)

project(libfemp
	LANGUAGES CXX
	VERSION 0.1.0
	)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)


list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

list(APPEND CMAKE_PREFIX_PATH "D:/Program Files/eigen-3.4.0/eigen-3.4.0")
list(APPEND CMAKE_PREFIX_PATH "D:/Program Files/boost_1_79_0")
list(APPEND CMAKE_PREFIX_PATH "C:/local/boost_1_79_0")

option(USES_RE2C "Uses re2c to regenerate the lexer files")
option(BUILD_TESTING "Determines whether tests should be configured" OFF )
option(LIBFEMP_BUILD_TESTING "Determines whether libfemp's tests should be configured" OFF)

if(USES_RE2C)
	message(STATUS "Adding re2c custom rules")
	find_program(RE2C re2c REQUIRED)

	add_custom_command( OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/src/io/import/FemJsonParser.c++
		DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/io/import/FemJsonParser.cpp.re2c ${CMAKE_CURRENT_SOURCE_DIR}/src/io/import/FemJsonParser.hpp
		COMMAND ${RE2C}
		ARGS -i -T -o ${CMAKE_CURRENT_SOURCE_DIR}/src/io/import/FemJsonParser.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/io/import/FemJsonParser.cpp.re2c
		COMMENT "Generating re2c files"
		)
else(USES_RE2C)
	message(STATUS "Skipping re2c custom rules")
endif()

#set(Eigen3_DIR  "D:/ThirdParty/eigen3")
#find_package(Eigen3)
find_package(Boost)
#find_package(suitesparse MODULE REQUIRED)

include_directories("D:/Program Files/eigen-3.4.0/eigen-3.4.0")
include_directories("D:/ThirdParty/eigen3")

add_subdirectory(libs)
add_subdirectory(src)

add_subdirectory(example)

if(BUILD_TESTING AND LIBFEMP_BUILD_TESTING)
	find_package(GTest CONFIG REQUIRED)

	enable_testing()
	add_subdirectory(tests)
endif()
