cmake_minimum_required(VERSION 3.14)
project(sctp_sab CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Find Gurobi
# Set GUROBI_HOME environment variable or adjust the path below
if(DEFINED ENV{GUROBI_HOME})
    set(GUROBI_HOME $ENV{GUROBI_HOME})
else()
    # Common default paths
    file(GLOB GUROBI_SEARCH_PATHS
        /opt/gurobi*/linux64
        /Library/gurobi*/macos_universal2
        $ENV{HOME}/gurobi*/linux64
    )
    if(GUROBI_SEARCH_PATHS)
        list(GET GUROBI_SEARCH_PATHS -1 GUROBI_HOME)
    endif()
endif()

if(GUROBI_HOME)
    message(STATUS "Gurobi home: ${GUROBI_HOME}")
    include_directories(${GUROBI_HOME}/include)
    link_directories(${GUROBI_HOME}/lib)

    # Find the versioned gurobi library
    file(GLOB GUROBI_LIB_FILES "${GUROBI_HOME}/lib/libgurobi*.so")
    foreach(lib ${GUROBI_LIB_FILES})
        get_filename_component(lib_name ${lib} NAME_WE)
        string(REGEX REPLACE "^lib" "" lib_name ${lib_name})
        if(NOT lib_name MATCHES "gurobi_c\\+\\+")
            set(GUROBI_C_LIB ${lib_name})
        endif()
    endforeach()
    if(NOT GUROBI_C_LIB)
        set(GUROBI_C_LIB "gurobi110")  # fallback
    endif()
    message(STATUS "Gurobi C library: ${GUROBI_C_LIB}")
else()
    message(FATAL_ERROR "Gurobi not found. Set GUROBI_HOME environment variable.")
endif()

# Source files
set(SOURCES
    src/main.cpp
    src/read.cpp
    src/gp.cpp
    src/network.cpp
    src/instance.cpp
)

add_executable(sctp_sab ${SOURCES})

target_include_directories(sctp_sab PRIVATE include)
target_link_libraries(sctp_sab
    Eigen3::Eigen
    gurobi_c++
    ${GUROBI_C_LIB}
)
