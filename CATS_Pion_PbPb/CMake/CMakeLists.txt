cmake_minimum_required(VERSION 3.0)
project(CMake_CATS)

# ================= PATHS ================= #
# destination of your tutorial project (where main.cpp / Basics.cpp etc. are)
set(PROJECT_DESTINATION "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_Pion_PbPb")

# CATS installation prefix (the folder that contains bin/cats-config)
set(CATS_PATH "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install")

# ROOT prefix (the folder that contains bin/root-config)
set(ROOT_PATH "/usr/local/Cellar/root/6.36.04_1")

# GSL (prefer brew /opt path; avoid Cellar version pin)
set(GSL_INCLUDE "/usr/local/opt/gsl/include")
set(GSL_INCLUDE2 "/usr/local/opt/gsl/include/gsl")
set(GSL_LIB "/usr/local/opt/gsl/lib")

# OpenMP (libomp from brew; avoid Cellar version pin)
set(OMP_INCLUDE "/usr/local/opt/libomp/include")
set(OMP_LIB "/usr/local/opt/libomp/lib")

# ================= OUTPUT ================= #
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_DESTINATION}/CMake/bin)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ================= EXECUTABLE ================= #
add_executable(${PROJECT_NAME}
  ${PROJECT_DESTINATION}/main.cpp
  ${PROJECT_DESTINATION}/Basics.cpp
  ${PROJECT_DESTINATION}/ExtendedCk.cpp
)

# ================= ROOT ================= #
find_program(ROOT_CONFIG_EXECUTABLE root-config HINTS
  "${ROOT_PATH}/bin"
  ENV PATH
)

if(NOT ROOT_CONFIG_EXECUTABLE)
  message(FATAL_ERROR "root-config not found. Either set ROOT_PATH correctly or source thisroot.sh so root-config is in PATH.")
endif()

execute_process(COMMAND bash -c "${ROOT_CONFIG_EXECUTABLE} --cflags" OUTPUT_VARIABLE ROOT_CXXFLAGS)
execute_process(COMMAND bash -c "${ROOT_CONFIG_EXECUTABLE} --libs"   OUTPUT_VARIABLE ROOT_LIBS)
execute_process(COMMAND bash -c "${ROOT_CONFIG_EXECUTABLE} --incdir" OUTPUT_VARIABLE ROOT_INCLUDE)

string(REGEX REPLACE "\n$" "" ROOT_CXXFLAGS "${ROOT_CXXFLAGS}")
string(REGEX REPLACE "\n$" "" ROOT_LIBS     "${ROOT_LIBS}")
string(REGEX REPLACE "\n$" "" ROOT_INCLUDE  "${ROOT_INCLUDE}")

# ================= CATS ================= #
execute_process(COMMAND bash -c "${CATS_PATH}/bin/cats-config --incdir" OUTPUT_VARIABLE CATS_INCLUDE)
execute_process(COMMAND bash -c "${CATS_PATH}/bin/cats-config --libdir" OUTPUT_VARIABLE CATS_LIBDIR)
string(REGEX REPLACE "\n$" "" CATS_INCLUDE "${CATS_INCLUDE}")
string(REGEX REPLACE "\n$" "" CATS_LIBDIR  "${CATS_LIBDIR}")

# ================= COMPILE FLAGS ================= #
# IMPORTANT: -lomp is LINK flag, NOT compile flag
set(CFLAGS " -O2 -std=c++20 -pthread ${ROOT_CXXFLAGS}")
string(APPEND CFLAGS " -Xpreprocessor -fopenmp")

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "${CFLAGS}")

# ================= INCLUDE DIRECTORIES ================= #
include_directories(
  ${PROJECT_DESTINATION}
  ${PROJECT_DESTINATION}/CATSpotentials
  ${CATS_INCLUDE}
  ${GSL_INCLUDE}
  ${GSL_INCLUDE2}
  ${OMP_INCLUDE}
  ${ROOT_INCLUDE}
)

# ================= LINK LIBRARIES ================= #
target_link_directories(${PROJECT_NAME} PRIVATE
  ${CATS_LIBDIR}
  ${GSL_LIB}
  ${OMP_LIB}
)

target_link_libraries(${PROJECT_NAME} PRIVATE
  CATS
  ${ROOT_LIBS}
  gsl gslcblas
  omp
)

# ================= RECORD CONFIG ================= #
file(WRITE  CMakeDLM.txt "${PROJECT_NAME}\n")
file(APPEND CMakeDLM.txt "${PROJECT_DESTINATION}\n")
file(APPEND CMakeDLM.txt "${ROOT_PATH}\n")
file(APPEND CMakeDLM.txt "${GSL_INCLUDE}\n")
file(APPEND CMakeDLM.txt "${GSL_LIB}\n")
file(APPEND CMakeDLM.txt "${OMP_INCLUDE}\n")
file(APPEND CMakeDLM.txt "${CATS_PATH}\n")
file(APPEND CMakeDLM.txt "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}\n")
