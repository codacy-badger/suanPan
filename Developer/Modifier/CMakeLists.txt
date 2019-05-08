cmake_minimum_required(VERSION 3.0)

project(ModifierExample)

add_definitions(-DSUANPAN_DLL)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

if(COMPILER_IDENTIFIER MATCHES "vs")
link_libraries(libopenblas)
else()
link_libraries(openblas)
endif()

add_library(${PROJECT_NAME} SHARED
    "ModifierExample.cpp"
    "../../Domain/Node.cpp"
    "../../Domain/Tag.cpp"
    "../../Element/Element.cpp"
    "../../Element/Modifier/Modifier.cpp"
    "../../Toolbox/debug.cpp"
    )
