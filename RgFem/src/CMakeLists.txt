set(TARGET_NAME RgFem)
project(${TARGET_NAME} LANGUAGES CXX)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
add_library(${TARGET_NAME} SHARED "")
target_compile_definitions(RgFem PRIVATE RgFem_Module)

target_include_directories(${TARGET_NAME} PRIVATE ${SOURCE_DIR})



set(Eigen_Root "D:/Program Files/eigen-3.4.0/eigen-3.4.0")
target_include_directories(${TARGET_NAME} PRIVATE ${Eigen_Root})

#set(Eigen3_INCLUDE_DIR "D:/ThirdParty/eigen3")

message(status "flzs****${CMAKE_MODULE_PATH}")

# 找到Eigen库
#find_package(Eigen3 3.3 REQUIRED)
 
# 如果找到Eigen库，则将其包含目录添加到你的项目中
if(Eigen3_FOUND)
  message(status "Eigen3 found****")
  include_directories(${Eigen3_INCLUDE_DIR})
endif()

file(GLOB rootfiles *.cpp *.h)
source_group("" FILES ${rootfiles})
target_sources(${TARGET_NAME} PRIVATE ${rootfiles})


sub_dir_list(SUB_DIRS ${SOURCE_DIR})
 

foreach(sub_dir ${SUB_DIRS})
    file(GLOB subdir_files ${sub_dir}/*.cpp ${sub_dir}/*.h)
    string(REPLACE "/" "\\" sub_dir ${sub_dir})
    source_group(${sub_dir} FILES ${subdir_files})
    target_sources(${TARGET_NAME} PRIVATE ${subdir_files})
endforeach()