project(example)

file(GLOB SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_executable(example
	${SRC_FILES}
)

target_link_libraries(example
		libfemp::Static
	)

target_include_directories(example
		PUBLIC
		${PROJECT_SOURCE_DIR}
		"E:/Lzs/femp/external/libfemp/src"
)


