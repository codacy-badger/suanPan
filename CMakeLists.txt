cmake_minimum_required(VERSION 3.0)

if(${CMAKE_VERSION} VERSION_GREATER "3.13.0")
    cmake_policy(SET CMP0076 OLD)
endif()

project(suanPan C CXX)

include(CheckLanguage)
check_language(Fortran)
if(CMAKE_Fortran_COMPILER)
    set(FORTRAN_STATUS TRUE)
    enable_language(Fortran)
else()
    set(FORTRAN_STATUS FALSE)
endif()

set(ROOT .)

include(${ROOT}/.config.cmake)

link_directories(${PROJECT_BINARY_DIR}/Libs)

add_executable(${PROJECT_NAME}
        "suanPan.h"
        "suanPan.cpp"
        "Constraint/CMakeLists.txt"
        "Database/CMakeLists.txt"
        "Domain/CMakeLists.txt"
        "Load/CMakeLists.txt"
        "Recorder/CMakeLists.txt"
        "Step/CMakeLists.txt"
        "Toolbox/CMakeLists.txt")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
    target_sources(${PROJECT_NAME} PRIVATE "Resource/suanPan_gcc.rc")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC") # MSVC COMPILER
    target_sources(${PROJECT_NAME} PRIVATE "Resource/suanPan.rc")
endif()

add_subdirectory(Converger)
add_subdirectory(Element)
add_subdirectory(Material)
add_subdirectory(Section)
add_subdirectory(Solver)

add_subdirectory(Constraint)
add_subdirectory(Database)
add_subdirectory(Domain)
add_subdirectory(Load)
add_subdirectory(Recorder)
add_subdirectory(Step)
add_subdirectory(Toolbox)

if(BUILD_DLL_EXAMPLE)
    add_subdirectory(Developer/Element)
    add_subdirectory(Developer/Material)
    add_subdirectory(Developer/Section)
    add_subdirectory(Developer/Modifier)
    add_subdirectory(Developer/ModuleBundle)
endif()

if(FORTRAN_STATUS)
    add_subdirectory(Toolbox/amd-src)
    add_subdirectory(Toolbox/lapack-ext)
    add_subdirectory(Toolbox/arpack-src)
endif()

add_subdirectory(Toolbox/superlu-src)

target_link_libraries(${PROJECT_NAME}
    Converger
    Element
    Material
    Section
    Solver
    superlu)
	
if((NOT FORTRAN_STATUS) AND (COMPILER_IDENTIFIER MATCHES "vs"))
    target_link_libraries(${PROJECT_NAME} libamd libarpack libspmm)
else()
    target_link_libraries(${PROJECT_NAME} amd arpack spmm)
endif()

if(COMPILER_IDENTIFIER MATCHES "vs")
    target_link_libraries(${PROJECT_NAME} libmumpsd libopenblas)
else()
    target_link_libraries(${PROJECT_NAME} mumpsd openblas)
endif()
