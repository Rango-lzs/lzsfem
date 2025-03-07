cmake_minimum_required(version 3.2)
project(RgFem_Tests)

set(RG_INCLUDE_DIR "")
include_directories(${RG_INCLUDE_DIR})

#set(SOURCE_FILES main.cpp src/divider_tests.cpp)
file(GLOB SRC_FILES *.cpp *.h src/*.cpp src/*.h)

add_executable(rg_tests ${SRC_FILES})
target_link_libraries(rg_tests RgFem gtest)
install(TARGETS rg_tests DESTINATION ${RG_INSTALL_BIN_DIR})

