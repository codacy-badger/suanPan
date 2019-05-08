cmake_minimum_required(VERSION 3.0)

project(ModuleBundle)

add_definitions(-DSUANPAN_DLL)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

include_directories(../Element)
include_directories(../Material)
include_directories(../Section)

if(COMPILER_IDENTIFIER MATCHES "vs")
link_libraries(libopenblas)
else()
link_libraries(openblas)
endif()

add_library(${PROJECT_NAME} SHARED
    "../../Domain/Node.cpp"
    "../../Domain/Tag.cpp"
    "../../Element/Element.cpp"
    "../../Material/Material.cpp"
    "../../Section/Section.cpp"
    "../../Toolbox/debug.cpp"
    "../Element/ElementExample.cpp"
    "../Material/MaterialExample.cpp"
    "../Section/SectionExample.cpp"
    )
