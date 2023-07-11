
cmake_minimum_required (VERSION 3.14)
project (KMC
    VERSION 0.1.0
    DESCRIPTION "A C++ project using CMake"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS On)

find_package(OpenMP REQUIRED)

if(DEBUG)
  set(CMAKE_BUILD_TYPE "Debug")
else()
  set(CMAKE_BUILD_TYPE "Release")
endif()

## Add external dependencies here, e.g.:

## The directory for the library
## Change the directory name from "KMC" to "KMC_func" in case of confusion
add_subdirectory(KMC_func)

## Only build documentation if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    ## Docs only available if this is the main project
    find_package(Doxygen QUIET)
    if(Doxygen_FOUND)
        add_subdirectory(docs)
    else()
        message(STATUS "Doxygen not found, not building docs")
    endif()
    ## The example executable is only available for main project
    #add_subdirectory(examples)
endif()

## Allow testing for main and subprojects if -DTESTS=TRUE
if(TESTS)
  if(NOT DEBUG)
    message("!! Diagnostic testing is not run when compiled in release mode.")
  endif()
    enable_testing()
    add_subdirectory(tests)
endif()
