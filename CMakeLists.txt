cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(RangoFem)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "CMakePredefinedTargets")

foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/$<CONFIG>)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/$<CONFIG>)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/$<CONFIG>)
endforeach()

macro(sub_dir_list result cur_dir)
    file(GLOB_RECURSE child_dir LIST_DIRECTORIES true RELATIVE ${cur_dir} ${cur_dir}/*)
    foreach(dir ${child_dir})
        if(IS_DIRECTORY ${cur_dir}/${dir})          
            list(APPEND ${result} ${dir})
            #message(status " ********* " result)
        endif()
    endforeach()   
endmacro(sub_dir_list)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "cmakes/")
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmakes)

set(THIRD_PARTY_ROOT "D:/Lzs/PhysicEngine/FEM_Project/lzsfem/ThirdParty") #office
set(THIRD_PARTY_ROOT "D:/ThirdParty") #home

option(ENABLE_UNIT_TESTS "Enable unit tests" OFF)

add_subdirectory(RgFem)
#add_subdirectory(libfemp)

if(ENABLE_UNIT_TESTS)
    set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT Tests)
    add_subdirectory(tests)
endif()

