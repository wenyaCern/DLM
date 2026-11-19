cmake_minimum_required(VERSION 3.16)
project(run_CATS3D)

get_filename_component(CATS3D_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CATS installation can be overridden on the farm with -DCATS_PATH=...
set(CATS_PATH "${CATS3D_ROOT}/../install" CACHE PATH "CATS installation prefix")

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
    ${CATS3D_ROOT}/fit/shared
    ${CATS3D_ROOT}/fit/common
    ${CATS3D_ROOT}
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
  ${CATS3D_ROOT}/fit/shared/Basics3D.cpp
  ${CATS3D_ROOT}/fit/shared/AdaptivePipiKernel.cpp
  ${CATS3D_ROOT}/fit/shared/CATS3DFastKernelCache.cpp
  ${CATS3D_ROOT}/fit/common/FitCATS3D.cpp
)
configure_cats3d_target(FitCATS3D)

add_executable(check_adaptive_kernel ${CATS3D_ROOT}/kernel_generation/check_adaptive_kernel.cpp ${CATS3D_ROOT}/fit/shared/AdaptivePipiKernel.cpp)
target_include_directories(check_adaptive_kernel PRIVATE ${CATS3D_ROOT})
configure_cats3d_target(check_adaptive_kernel)

add_library(FitCATS3DResonance SHARED
  ${CATS3D_ROOT}/fit/shared/Basics3D.cpp
  ${CATS3D_ROOT}/fit/shared/AdaptivePipiKernel.cpp
  ${CATS3D_ROOT}/fit/shared/CATS3DFastKernelCache.cpp
  ${CATS3D_ROOT}/fit/common/FitCATS3D.cpp
  ${CATS3D_ROOT}/fit/resonance/CATS3DResonanceSourceBuilder.cpp
)
target_compile_definitions(FitCATS3DResonance PRIVATE FITCATS3D_ENABLE_RESONANCE)
configure_cats3d_target(FitCATS3DResonance)

add_executable(run_CATS3D
  ${CATS3D_ROOT}/direct_cf/main.cpp
)
configure_cats3d_target(run_CATS3D)

add_executable(fit_CATS3D
  ${CATS3D_ROOT}/fit/gaussian/fit_CATS3D_main.cpp
)
configure_cats3d_target(fit_CATS3D)
target_link_libraries(fit_CATS3D PRIVATE FitCATS3D)

add_executable(fit_CATS3D_oo
  ${CATS3D_ROOT}/fit/nuclearStructureDatafit/fit_CATS3D_oo_main.cpp
)
configure_cats3d_target(fit_CATS3D_oo)
target_link_libraries(fit_CATS3D_oo PRIVATE FitCATS3D)

add_executable(fit_CATS3D_derived_full
  ${CATS3D_ROOT}/fit/gaussian/fit_CATS3D_derived_full_main.cpp
)
configure_cats3d_target(fit_CATS3D_derived_full)
target_link_libraries(fit_CATS3D_derived_full PRIVATE FitCATS3D)

add_executable(compare_cats3d_momentum_bins
  ${CATS3D_ROOT}/other/studies/compare_cats3d_momentum_bins.cpp
)
configure_cats3d_target(compare_cats3d_momentum_bins)
target_link_libraries(compare_cats3d_momentum_bins PRIVATE FitCATS3D)

add_executable(compare_cats3d_bin_effect_vs_radius
  ${CATS3D_ROOT}/other/studies/compare_cats3d_bin_effect_vs_radius.cpp
)
configure_cats3d_target(compare_cats3d_bin_effect_vs_radius)
target_link_libraries(compare_cats3d_bin_effect_vs_radius PRIVATE FitCATS3D)

add_executable(build_local_source_kernel_CATS3D
  ${CATS3D_ROOT}/kernel_generation/build_local_source_kernel_CATS3D_main.cpp
)
configure_cats3d_target(build_local_source_kernel_CATS3D)
target_link_libraries(build_local_source_kernel_CATS3D PRIVATE FitCATS3D)

add_executable(diagnose_cats3d_source_range
  ${CATS3D_ROOT}/other/studies/diagnose_cats3d_source_range.cpp
)
configure_cats3d_target(diagnose_cats3d_source_range)
target_link_libraries(diagnose_cats3d_source_range PRIVATE FitCATS3D)

add_executable(draw_approx_r40_fixedR_slices
  ${CATS3D_ROOT}/other/studies/draw_approx_r40_fixedR_slices.cpp
)
configure_cats3d_target(draw_approx_r40_fixedR_slices)
target_link_libraries(draw_approx_r40_fixedR_slices PRIVATE FitCATS3D)

add_executable(compare_momentum5_vs10_fixedR
  ${CATS3D_ROOT}/other/studies/compare_momentum5_vs10_fixedR.cpp
)
configure_cats3d_target(compare_momentum5_vs10_fixedR)
target_link_libraries(compare_momentum5_vs10_fixedR PRIVATE FitCATS3D)

add_executable(compare_spatial40_vs60_fixedR
  ${CATS3D_ROOT}/other/studies/compare_spatial40_vs60_fixedR.cpp
)
configure_cats3d_target(compare_spatial40_vs60_fixedR)
target_link_libraries(compare_spatial40_vs60_fixedR PRIVATE FitCATS3D)

add_executable(scan_radius_lambda_peak_shoulder
  ${CATS3D_ROOT}/other/studies/scan_radius_lambda_peak_shoulder.cpp
)
configure_cats3d_target(scan_radius_lambda_peak_shoulder)
target_link_libraries(scan_radius_lambda_peak_shoulder PRIVATE FitCATS3D)

add_executable(calculate_cats3d_r40_q5_octant
  ${CATS3D_ROOT}/direct_cf/calculate_cats3d_r40_q5_octant.cpp
)
configure_cats3d_target(calculate_cats3d_r40_q5_octant)
target_link_libraries(calculate_cats3d_r40_q5_octant PRIVATE FitCATS3D)

add_executable(build_adaptive_pipi_coulomb_kernel
  ${CATS3D_ROOT}/kernel_generation/build_adaptive_pipi_coulomb_kernel.cpp
)
configure_cats3d_target(build_adaptive_pipi_coulomb_kernel)

add_executable(fit_CATS3D_resonance
  ${CATS3D_ROOT}/fit/resonance/fit_CATS3D_resonance_main.cpp
)
configure_cats3d_target(fit_CATS3D_resonance)
target_link_libraries(fit_CATS3D_resonance PRIVATE FitCATS3DResonance)

add_executable(fit_resonance6291_mt0
  ${CATS3D_ROOT}/fit/resonance/fit_resonance6291_mt0_main.cpp
)
configure_cats3d_target(fit_resonance6291_mt0)
target_link_libraries(fit_resonance6291_mt0 PRIVATE FitCATS3DResonance)

add_executable(source_compare_CATS3D
  ${CATS3D_ROOT}/direct_cf/source_compare_CATS3D_main.cpp
)
configure_cats3d_target(source_compare_CATS3D)
target_link_libraries(source_compare_CATS3D PRIVATE FitCATS3DResonance)

add_executable(r20_check_CATS3D
  ${CATS3D_ROOT}/direct_cf/r20_check_CATS3D_main.cpp
)
configure_cats3d_target(r20_check_CATS3D)
target_link_libraries(r20_check_CATS3D PRIVATE FitCATS3D)

add_executable(kernel_cache_check_CATS3D
  ${CATS3D_ROOT}/kernel_generation/kernel_cache_check_CATS3D_main.cpp
)
configure_cats3d_target(kernel_cache_check_CATS3D)
target_link_libraries(kernel_cache_check_CATS3D PRIVATE FitCATS3D)

add_executable(build_data_kernel_CATS3D
  ${CATS3D_ROOT}/kernel_generation/build_data_kernel_CATS3D_main.cpp
)
configure_cats3d_target(build_data_kernel_CATS3D)
target_link_libraries(build_data_kernel_CATS3D PRIVATE FitCATS3D)

add_executable(convolve_kernel_gauss_CATS3D
  ${CATS3D_ROOT}/kernel_generation/convolve_kernel_gauss_CATS3D_main.cpp
)
configure_cats3d_target(convolve_kernel_gauss_CATS3D)
target_link_libraries(convolve_kernel_gauss_CATS3D PRIVATE FitCATS3D)

add_executable(scan_radius_external_kernel_CATS3D
  ${CATS3D_ROOT}/other/studies/scan_radius_external_kernel_CATS3D_main.cpp
)
configure_cats3d_target(scan_radius_external_kernel_CATS3D)
target_link_libraries(scan_radius_external_kernel_CATS3D PRIVATE FitCATS3D)
