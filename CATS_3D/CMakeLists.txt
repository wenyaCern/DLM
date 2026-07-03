cmake_minimum_required(VERSION 3.16)
project(run_CATS3D)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CATS path
set(CATS_PATH "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install")

# ROOT
find_program(ROOT_CONFIG_EXECUTABLE root-config HINTS
  "/usr/local/bin"
  "/usr/local/Cellar/root/6.36.04_1/bin"
  ENV PATH
)

if(NOT ROOT_CONFIG_EXECUTABLE)
  message(FATAL_ERROR "root-config not found")
endif()

execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --cflags OUTPUT_VARIABLE ROOT_CXXFLAGS)
execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --libs   OUTPUT_VARIABLE ROOT_LIBS)
string(STRIP "${ROOT_CXXFLAGS}" ROOT_CXXFLAGS)
string(STRIP "${ROOT_LIBS}" ROOT_LIBS)
separate_arguments(ROOT_CXXFLAGS)
separate_arguments(ROOT_LIBS)

# CATS
execute_process(COMMAND ${CATS_PATH}/bin/cats-config --incdir OUTPUT_VARIABLE CATS_INCLUDE)
execute_process(COMMAND ${CATS_PATH}/bin/cats-config --libdir OUTPUT_VARIABLE CATS_LIBDIR)
string(STRIP "${CATS_INCLUDE}" CATS_INCLUDE)
string(STRIP "${CATS_LIBDIR}" CATS_LIBDIR)

function(configure_cats3d_target target_name)
  target_compile_options(${target_name} PRIVATE ${ROOT_CXXFLAGS})
  target_include_directories(${target_name} PRIVATE
    ${CATS_INCLUDE}
  )
  target_link_directories(${target_name} PRIVATE
    ${CATS_LIBDIR}
  )
  target_link_libraries(${target_name} PRIVATE
    CATS
    Minuit
    ${ROOT_LIBS}
  )
endfunction()

add_library(FitCATS3D SHARED
  Basics3D.cpp
  FitCATS3D.cpp
)
configure_cats3d_target(FitCATS3D)

add_library(FitCATS3DResonance SHARED
  Basics3D.cpp
  FitCATS3D.cpp
  CATS3DResonanceSourceBuilder.cpp
)
target_compile_definitions(FitCATS3DResonance PRIVATE FITCATS3D_ENABLE_RESONANCE)
configure_cats3d_target(FitCATS3DResonance)

add_executable(run_CATS3D
  main.cpp
)
configure_cats3d_target(run_CATS3D)

add_executable(fit_CATS3D
  fit_CATS3D_main.cpp
)
configure_cats3d_target(fit_CATS3D)
target_link_libraries(fit_CATS3D PRIVATE FitCATS3D)

add_executable(fit_CATS3D_resonance
  fit_CATS3D_resonance_main.cpp
)
configure_cats3d_target(fit_CATS3D_resonance)
target_link_libraries(fit_CATS3D_resonance PRIVATE FitCATS3DResonance)

add_executable(source_compare_CATS3D
  source_compare_CATS3D_main.cpp
)
configure_cats3d_target(source_compare_CATS3D)
target_link_libraries(source_compare_CATS3D PRIVATE FitCATS3DResonance)

add_executable(r20_check_CATS3D
  r20_check_CATS3D_main.cpp
)
configure_cats3d_target(r20_check_CATS3D)
target_link_libraries(r20_check_CATS3D PRIVATE FitCATS3D)
