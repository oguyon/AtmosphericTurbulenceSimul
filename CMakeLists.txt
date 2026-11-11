# =====================================================================
# milkatmturb - Atmospheric Turbulence Simulation Plugin for milk
# =====================================================================

cmake_minimum_required(VERSION 3.16)

if(POLICY CMP0104)
    cmake_policy(SET CMP0104 NEW)
endif()

set(LIBNAME "milkatmturb")
set(SRCNAME "atmturb")

message(" SRCNAME = ${SRCNAME} -> LIBNAME = ${LIBNAME}")

# Project definition
project(lib_${LIBNAME}_project C)

# Set C standard: adopt C23 by default if supported by the compiler, otherwise fall back to C11
if(NOT CMAKE_C_STANDARD)
    include(CheckCCompilerFlag)
    check_c_compiler_flag("-std=gnu23" COMPILER_SUPPORTS_GNU23)
    if(COMPILER_SUPPORTS_GNU23)
        set(CMAKE_C_STANDARD 23)
    else()
        check_c_compiler_flag("-std=c23" COMPILER_SUPPORTS_C23)
        if(COMPILER_SUPPORTS_C23)
            set(CMAKE_C_STANDARD 23)
        else()
            set(CMAKE_C_STANDARD 11)
        endif()
    endif()
endif()
set(CMAKE_C_STANDARD_REQUIRED ON)
message(STATUS "milkatmturb: Using C standard C${CMAKE_C_STANDARD}")

# Default to Release build type for maximum optimization (-O3)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Find OpenMP (required for multi-threaded turbulence simulation)
find_package(OpenMP REQUIRED)

# Find FFTW3 (for Fourier domain phase screens)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED fftw3)
pkg_check_modules(FFTWF REQUIRED fftw3f)
find_library(FFTWF_OMP_LIB NAMES fftw3f_omp fftw3f_threads)
if(FFTWF_OMP_LIB)
    list(APPEND FFTWF_LIBRARIES ${FFTWF_OMP_LIB})
endif()
find_library(FFTW_OMP_LIB NAMES fftw3_omp fftw3_threads)
if(FFTW_OMP_LIB)
    list(APPEND FFTW_LIBRARIES ${FFTW_OMP_LIB})
endif()

# Collect source files
set(SOURCEFILES
    src/milkatmturb.c
    src/AtmosphereModel/AtmosphereModel.c
    src/AtmosphereModel/atmmod_ria_table.c
    src/AtmosphereModel/atmmod_gas_mixture.c
    src/AtmosphereModel/atmmod_standard_model.c
    src/AtmosphereModel/atmmod_refraction.c
    src/AtmosphereModel/atmmod_config_parser.c
    src/AtmosphericTurbulence/nrlmsise-00.20131225/nrlmsise-00.c
    src/AtmosphericTurbulence/nrlmsise-00.20131225/nrlmsise-00_data.c
    src/AtmosphericTurbulence/AtmosphericTurbulence.c
    src/AtmosphericTurbulence/atmturb_config.c
    src/AtmosphericTurbulence/atmturb_profile.c
    src/AtmosphericTurbulence/atmturb_geometry.c
    src/AtmosphericTurbulence/atmturb_screens.c
    src/AtmosphericTurbulence/atmturb_wind.c
    src/AtmosphericTurbulence/atmturb_hvturb.c
    src/AtmosphericTurbulence/atmturb_cube_ops.c
    src/AtmosphericTurbulence/atmturb_psf_metrics.c
    src/AtmosphericTurbulence/atmturb_psf_sim.c
    src/AtmosphericTurbulence/atmturb_linpred_full.c
    src/AtmosphericTurbulence/atmturb_linpred_pixel.c
    src/AtmosphericTurbulence/atmturb_linpred_test.c
    src/AtmosphericTurbulence/atmturb_wfs_sim.c
    src/AtmosphericTurbulence/atmturb_wfs_render.c
    src/AtmosphericTurbulence/atmturb_wfs_render_cuda.c
    src/AtmosphericTurbulence/atmturb_wfs_cuda_common.c
    src/AtmosphericTurbulence/atmturb_wfs_stream_cuda.c
    src/AtmosphericTurbulence/atmturb_wfs_stream.c
    src/AtmosphericTurbulence/atmturb_superlayer.c
    src/AtmosphericTurbulence/atmturb_superlayer_nodes.c
    src/AtmosphericTurbulence/atmturb_fresnel.c
    src/AtmosphericTurbulence/atmturb_fresnel_simd_scalar.c
    src/AtmosphericTurbulence/atmturb_fresnel_simd_avx2.c
    src/AtmosphericTurbulence/atmturb_fresnel_simd_dispatch.c
    src/AtmosphericTurbulence/atmturb_rytov.c
    src/AtmosphericTurbulence/atmturb_rytov_kernels.c
    src/AtmosphericTurbulence/atmturb_rytov_simd_scalar.c
    src/AtmosphericTurbulence/atmturb_rytov_simd_avx2.c
    src/AtmosphericTurbulence/atmturb_rytov_simd_dispatch.c
    src/AtmosphericTurbulence/atmturb_lowfreq.c
    src/AtmosphericTurbulence/atmturb_rolling.c
    src/AtmosphericTurbulence/atmturb_simd_scalar.c
    src/AtmosphericTurbulence/atmturb_simd_avx2.c
    src/AtmosphericTurbulence/atmturb_simd_avx2_ops.c
    src/AtmosphericTurbulence/atmturb_simd_avx512.c
    src/AtmosphericTurbulence/atmturb_simd_dispatch.c
    src/AtmosphericTurbulence/atmturb_mkwfs_FPS.c
    src/AtmosphericTurbulence/atmturb_mkhvturb_FPS.c
    src/AtmosphericTurbulence/atmturb_mkmastert_FPS.c
    src/AtmosphericTurbulence/atmturb_mkvonkarman_FPS.c
    src/AtmosphericTurbulence/atmturb_aoloop_FPS.c
    src/AtmosphericTurbulence/atmturb_cli.c
)

set(INCLUDEFILES
    src/milkatmturb.h
    src/atmturb_c23_compat.h
    src/AtmosphereModel/AtmosphereModel.h
    src/AtmosphericTurbulence/AtmosphericTurbulence.h
    src/AtmosphericTurbulence/atmturb_profile.h
    src/AtmosphericTurbulence/atmturb_geometry.h
    src/AtmosphericTurbulence/atmturb_superlayer.h
    src/AtmosphericTurbulence/atmturb_fresnel.h
    src/AtmosphericTurbulence/atmturb_fresnel_simd.h
    src/AtmosphericTurbulence/atmturb_rytov.h
    src/AtmosphericTurbulence/atmturb_rytov_simd.h
    src/AtmosphericTurbulence/atmturb_lowfreq.h
    src/AtmosphericTurbulence/atmturb_rolling.h
    src/AtmosphericTurbulence/atmturb_simd.h
    src/AtmosphericTurbulence/atmturb_wfs_render.h
    src/AtmosphericTurbulence/atmturb_wfs_render_cuda.h
    src/AtmosphericTurbulence/atmturb_wfs_cuda_common.h
    src/AtmosphericTurbulence/atmturb_wfs_stream.h
    src/AtmosphericTurbulence/atmturb_mkwfs_FPS.h
    src/AtmosphericTurbulence/atmturb_mkhvturb_FPS.h
    src/AtmosphericTurbulence/atmturb_mkmastert_FPS.h
    src/AtmosphericTurbulence/atmturb_mkvonkarman_FPS.h
    src/AtmosphericTurbulence/atmturb_aoloop_FPS.h
    src/AtmosphericTurbulence/atmturb_cli.h
)

file(GLOB SCRIPTS "src/AtmosphericTurbulence/scripts/*")

if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    set_source_files_properties(
        src/AtmosphericTurbulence/atmturb_simd_avx2.c
        src/AtmosphericTurbulence/atmturb_simd_avx2_ops.c
        src/AtmosphericTurbulence/atmturb_rytov_simd_avx2.c
        src/AtmosphericTurbulence/atmturb_fresnel_simd_avx2.c
        PROPERTIES COMPILE_FLAGS "-mavx2 -mfma")
    set_source_files_properties(src/AtmosphericTurbulence/atmturb_simd_avx512.c
        PROPERTIES COMPILE_FLAGS "-mavx512f -mavx512dq -mavx512bw -mavx512vl -mfma")
endif()

# CUDA GPU Acceleration
if(HAVE_CUDA OR USE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
            if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24")
                set(CMAKE_CUDA_ARCHITECTURES "native")
            else()
                set(CMAKE_CUDA_ARCHITECTURES "70;75;80;86")
            endif()
        endif()
        enable_language(CUDA OPTIONAL)
        find_package(CUDAToolkit QUIET)
        if(CUDAToolkit_FOUND)
            set(USE_CUDA_ACCEL TRUE)
            add_compile_definitions(HAVE_CUDA)
            list(APPEND SOURCEFILES
                src/AtmosphericTurbulence/atmturb_cuda.cu
                src/AtmosphericTurbulence/atmturb_cuda_piston.cu
                src/AtmosphericTurbulence/atmturb_cuda_rytov.cu
                src/AtmosphericTurbulence/atmturb_cuda_rytov_extrude.cu
                src/AtmosphericTurbulence/atmturb_cuda_rytov_filter.cu
                src/AtmosphericTurbulence/atmturb_cuda_fresnel_kernels.cu
                src/AtmosphericTurbulence/atmturb_cuda_fresnel.cu
            )
            list(APPEND INCLUDEFILES
                src/AtmosphericTurbulence/atmturb_cuda.h
                src/AtmosphericTurbulence/atmturb_cuda_piston.h
                src/AtmosphericTurbulence/atmturb_cuda_rytov.h
                src/AtmosphericTurbulence/atmturb_cuda_fresnel_kernels.h
                src/AtmosphericTurbulence/atmturb_cuda_fresnel.h
            )
            message(STATUS "milkatmturb: CUDA GPU acceleration ENABLED")
        endif()
    endif()
endif()

set(_NRL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence/nrlmsise-00.20131225")
set_source_files_properties(
    src/AtmosphericTurbulence/nrlmsise-00.20131225/nrlmsise-00.c
    PROPERTIES COMPILE_FLAGS "-Wno-absolute-value"
)

if(COMMAND milk_add_runnable_module)
    # In-tree build under milk (feat/bash-engine)
    milk_add_runnable_module(${LIBNAME}
        SOURCES ${SOURCEFILES}
        LINK_LIBS
            milkOpticsMaterials
            milkWFpropagate
            milkCOREMODiofits
            milkfft
            milkpsf
            milklinalgebra
            milkimagegen
            milkimagefilter
            milkimagebasic
            milkstatistic
            milklinoptimtools
            OpenMP::OpenMP_C
            ${FFTW_LIBRARIES}
            ${FFTWF_LIBRARIES}
    )
    milk_runnable_module_usage(${LIBNAME}
        INCLUDES
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphereModel>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence>
            $<BUILD_INTERFACE:${_NRL_DIR}>
            $<INSTALL_INTERFACE:include>
            $<INSTALL_INTERFACE:include/${SRCNAME}>
    )
    if(USE_CUDA_ACCEL)
        target_link_libraries(${LIBNAME} PRIVATE CUDA::cudart CUDA::cufft)
    endif()
    milk_add_fpsexec_stubs(${LIBNAME})

    install(FILES ${INCLUDEFILES} DESTINATION include/${SRCNAME})
    install(PROGRAMS ${SCRIPTS} DESTINATION bin)
else()
    # Out-of-tree / standalone build against milk (feat/bash-engine)
    set(MILK_SOURCE_CANDIDATES
        ${MILK_SOURCE_DIR}
        $ENV{MILK_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/../../src-frameworkdev/milk
        ${CMAKE_CURRENT_SOURCE_DIR}/../milk
        ${CMAKE_CURRENT_SOURCE_DIR}/../milk-framework-dev
        $ENV{MILK_ROOT}
        /usr/local/milk
    )
    set(MILK_SOURCE_FOUND "")
    foreach(cand ${MILK_SOURCE_CANDIDATES})
        if(cand AND EXISTS "${cand}/src/coremods")
            get_filename_component(MILK_SOURCE_FOUND "${cand}" ABSOLUTE)
            message(STATUS "Found MILK source tree: ${MILK_SOURCE_FOUND}")
            break()
        endif()
    endforeach()

    set(_root_cands ${MILK_ROOT} $ENV{MILK_ROOT} $ENV{MILK_INSTALLDIR} ${CMAKE_INSTALL_PREFIX})
    if(MILK_SOURCE_FOUND)
        file(GLOB _ver_inst "${MILK_SOURCE_FOUND}/_install/milk-*")
        list(APPEND _root_cands ${_ver_inst} "${MILK_SOURCE_FOUND}/_install")
    endif()
    list(APPEND _root_cands /usr/local/milk /usr/local)
    foreach(rcand ${_root_cands})
        if(rcand AND EXISTS "${rcand}/include/libcomputeunit/milk_runnable.h")
            get_filename_component(MILK_ROOT "${rcand}" ABSOLUTE)
            break()
        endif()
    endforeach()

    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND MILK_ROOT)
        set(CMAKE_INSTALL_PREFIX "${MILK_ROOT}" CACHE PATH "Install prefix" FORCE)
    endif()
    message(STATUS "milkatmturb: MILK_ROOT=${MILK_ROOT}, INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}")

    set(CMAKE_INSTALL_RPATH
        "$ORIGIN;${CMAKE_INSTALL_PREFIX}/lib;${MILK_ROOT}/lib;/usr/local/milk/lib;/usr/local/lib")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

    find_path(MILK_INCLUDE_DIR
        NAMES libcomputeunit/milk_runnable.h CommandLineInterface/CLIcore.h
        HINTS
            ${MILK_ROOT}/include
            $ENV{MILK_ROOT}/include
            ${CMAKE_INSTALL_PREFIX}/include
            /usr/local/milk/include
            /usr/local/include
        NO_DEFAULT_PATH
    )
    if(NOT MILK_INCLUDE_DIR)
        find_path(MILK_INCLUDE_DIR
            NAMES libcomputeunit/milk_runnable.h CommandLineInterface/CLIcore.h
            HINTS /usr/local/milk/include /usr/local/include
        )
    endif()

    set(MILK_STANDALONE_INCLUDES
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphereModel
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence/nrlmsise-00.20131225
        ${CMAKE_CURRENT_SOURCE_DIR}/../OpticsMaterials
        ${CMAKE_CURRENT_SOURCE_DIR}/../WFpropagate
    )

    if(MILK_INCLUDE_DIR)
        list(APPEND MILK_STANDALONE_INCLUDES
            ${MILK_INCLUDE_DIR}
            ${MILK_INCLUDE_DIR}/CommandLineInterface
            ${MILK_INCLUDE_DIR}/CommandLineInterface/CLIcore
            ${MILK_INCLUDE_DIR}/libmilkcommon
            ${MILK_INCLUDE_DIR}/libmilkdata
            ${MILK_INCLUDE_DIR}/libfps
            ${MILK_INCLUDE_DIR}/libprocessinfo
            ${MILK_INCLUDE_DIR}/libcomputeunit
            ${MILK_INCLUDE_DIR}/ImageStreamIO
        )
    endif()

    if(MILK_SOURCE_FOUND)
        list(APPEND MILK_STANDALONE_INCLUDES
            ${MILK_SOURCE_FOUND}/src
            ${MILK_SOURCE_FOUND}/src/coremods
            ${MILK_SOURCE_FOUND}/plugins/milk-extra-src
            ${MILK_SOURCE_FOUND}/plugins/OpticsMaterials
            ${MILK_SOURCE_FOUND}/plugins/WFpropagate
            ${MILK_SOURCE_FOUND}/src/cli
            ${MILK_SOURCE_FOUND}/src/cli/CLIcore
            ${MILK_SOURCE_FOUND}/src/engine
            ${MILK_SOURCE_FOUND}/src/engine/libfps
            ${MILK_SOURCE_FOUND}/src/engine/libmilkcommon
            ${MILK_SOURCE_FOUND}/src/engine/libprocessinfo
            ${MILK_SOURCE_FOUND}/src/engine/libmilkdata
            ${MILK_SOURCE_FOUND}/src/engine/libcomputeunit
            ${MILK_SOURCE_FOUND}/src/engine/ImageStreamIO
        )
        if(EXISTS "${MILK_SOURCE_FOUND}/_build/config")
            list(APPEND MILK_STANDALONE_INCLUDES "${MILK_SOURCE_FOUND}/_build/config")
        endif()
    endif()

    add_library(${LIBNAME} SHARED ${SOURCEFILES})
    set_target_properties(${LIBNAME} PROPERTIES
        LINKER_LANGUAGE C
        INTERPROCEDURAL_OPTIMIZATION OFF
        INSTALL_RPATH "${CMAKE_INSTALL_RPATH}"
    )
    target_compile_definitions(${LIBNAME} PRIVATE MILK_NO_CLI)
    target_compile_options(${LIBNAME} PRIVATE
        $<$<COMPILE_LANGUAGE:C>:-Wall -Wextra -Wno-unused-const-variable -Wno-stringop-truncation>
    )
    target_include_directories(${LIBNAME} SYSTEM PRIVATE ${MILK_STANDALONE_INCLUDES})
    target_include_directories(${LIBNAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphereModel>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence>
        $<BUILD_INTERFACE:${_NRL_DIR}>
        $<INSTALL_INTERFACE:include>
        $<INSTALL_INTERFACE:include/${SRCNAME}>
    )

    set(MILK_LIB_NAMES
        milkcomputeunit
        milkfps
        milkprocessinfo
        milkdata
        ImageStreamIO
        milkOpticsMaterials
        milkWFpropagate
        milkCOREMODmemory
        milkCOREMODarith
        milkCOREMODiofits
        milkCOREMODtools
        milkfft
        milkpsf
        milklinalgebra
        milkimagegen
        milkimagefilter
        milkimagebasic
        milkstatistic
        milklinoptimtools
    )

    set(MILK_SEARCH_HINTS
        ${MILK_ROOT}/lib
        ${CMAKE_INSTALL_PREFIX}/lib
        ${CMAKE_CURRENT_SOURCE_DIR}/../OpticsMaterials/_build
        ${CMAKE_CURRENT_SOURCE_DIR}/../WFpropagate/_build
    )
    if(MILK_SOURCE_FOUND)
        list(APPEND MILK_SEARCH_HINTS ${MILK_SOURCE_FOUND}/_build/lib)
    endif()
    list(APPEND MILK_SEARCH_HINTS /usr/local/milk/lib /usr/local/lib)

    set(MILK_LIBRARIES "")
    foreach(lib ${MILK_LIB_NAMES})
        find_library(MILK_${lib}_LIB NAMES ${lib} HINTS ${MILK_SEARCH_HINTS} NO_DEFAULT_PATH)
        if(NOT MILK_${lib}_LIB)
            find_library(MILK_${lib}_LIB NAMES ${lib} HINTS ${MILK_SEARCH_HINTS})
        endif()
        if(MILK_${lib}_LIB)
            list(APPEND MILK_LIBRARIES ${MILK_${lib}_LIB})
        else()
            message(WARNING
                "MILK library '${lib}' not found in search paths; falling back to -l${lib}")
            list(APPEND MILK_LIBRARIES ${lib})
        endif()
    endforeach()

    target_link_libraries(${LIBNAME} PUBLIC
        ${MILK_LIBRARIES}
        OpenMP::OpenMP_C
        ${FFTW_LIBRARIES}
        ${FFTWF_LIBRARIES}
        m
    )

    if(USE_CUDA_ACCEL)
        target_link_libraries(${LIBNAME} PRIVATE CUDA::cudart CUDA::cufft)
    endif()

    install(TARGETS ${LIBNAME} DESTINATION lib)
    install(FILES ${INCLUDEFILES} DESTINATION include/${SRCNAME})
    install(PROGRAMS ${SCRIPTS} DESTINATION bin)

    # ── Standalone fpsexec stubs parsed from EXPORT_FPSEXEC(...) ───────
    set(_FPSEXEC_STUB_TPL [=[
#include "milk_skinny_runnable.h"

extern const MILK_RUNNABLE @UNIT@;

int main(int argc, char *argv[])
{
    MILK_EMBED_BUILD_TAG();
    return MILK_FPSEXEC_MAIN(argc, argv, &@UNIT@);
}
]=])

    string(CONCAT _FPSEXEC_RE
        "EXPORT_FPSEXEC[ \t]*\\([ \t]*([A-Za-z0-9_]+)[ \t]*,[ \t]*"
        "\"([^\"]+)\"[ \t]*,[ \t]*([A-Za-z_][A-Za-z0-9_]*)[ \t]*\\)"
    )

    foreach(_rel_src IN LISTS SOURCEFILES)
        set(_abs_src "${CMAKE_CURRENT_SOURCE_DIR}/${_rel_src}")
        if(NOT EXISTS "${_abs_src}")
            continue()
        endif()
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_abs_src}")
        file(READ "${_abs_src}" _src_text)
        string(REGEX MATCHALL "(^|\n)[ \t]*EXPORT_FPSEXEC[ \t]*\\([^\n]*\\)"
               _fpsexec_hits "${_src_text}")
        foreach(_hit IN LISTS _fpsexec_hits)
            if(_hit MATCHES "${_FPSEXEC_RE}")
                set(_exe "${CMAKE_MATCH_1}-fpsexec-${CMAKE_MATCH_2}")
                set(_unit "${CMAKE_MATCH_3}")
                string(REPLACE "@UNIT@" "${_unit}" _stub_src "${_FPSEXEC_STUB_TPL}")
                set(_stub_file "${CMAKE_CURRENT_BINARY_DIR}/fpsexec_stubs/${_exe}.c")
                file(GENERATE OUTPUT "${_stub_file}" CONTENT "${_stub_src}")
                add_executable(${_exe} "${_stub_file}")
                target_compile_definitions(${_exe} PRIVATE MILK_NO_CLI MILK_BUILD_BINNAME="${_exe}")
                target_include_directories(${_exe} SYSTEM PRIVATE ${MILK_STANDALONE_INCLUDES})
                target_link_libraries(${_exe} PRIVATE
                    ${LIBNAME}
                    ${MILK_LIBRARIES}
                    m
                    rt
                    pthread
                    -Wl,--allow-shlib-undefined
                )
                set_target_properties(${_exe} PROPERTIES
                    INSTALL_RPATH "$ORIGIN/../lib;${CMAKE_INSTALL_RPATH}"
                )
                install(TARGETS ${_exe} DESTINATION bin)
            endif()
        endforeach()
    endforeach()

    # Refresh milk module catalog after install when milk-cli is present
    install(CODE "
        set(_milk_cli \"\${CMAKE_INSTALL_PREFIX}/bin/milk-cli\")
        if(EXISTS \"\${_milk_cli}\")
            message(STATUS \"Refreshing milk module catalog via \${_milk_cli} --catalog\")
            execute_process(COMMAND \"\${_milk_cli}\" --catalog RESULT_VARIABLE _cat_rc)
        endif()
    ")

    # ── Statistical test validation executable ─────────────────────────
    add_executable(atmturb-validate-stats
        tests/atmturb_validate_stats.c
        tests/atmturb_validate_math.c
    )
    target_compile_definitions(atmturb-validate-stats PRIVATE MILK_NO_CLI)
    target_include_directories(atmturb-validate-stats PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
        ${CMAKE_CURRENT_SOURCE_DIR}/tests
    )
    target_include_directories(atmturb-validate-stats SYSTEM PRIVATE ${MILK_STANDALONE_INCLUDES})
    target_link_libraries(atmturb-validate-stats PRIVATE
        ${LIBNAME}
        ${MILK_LIBRARIES}
        ${FFTWF_LIBRARIES}
        ${FFTW_LIBRARIES}
        cfitsio
        m
    )

    # ── Performance benchmark executable ──────────────────────────────
    add_executable(atmturb-benchmark-perf
        tests/benchmark_perf.c
    )
    target_compile_definitions(atmturb-benchmark-perf PRIVATE MILK_NO_CLI)
    target_include_directories(atmturb-benchmark-perf PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
        ${CMAKE_CURRENT_SOURCE_DIR}/tests
    )
    target_include_directories(atmturb-benchmark-perf SYSTEM PRIVATE ${MILK_STANDALONE_INCLUDES})
    target_link_libraries(atmturb-benchmark-perf PRIVATE
        ${LIBNAME}
        ${MILK_LIBRARIES}
        ${FFTWF_LIBRARIES}
        ${FFTW_LIBRARIES}
        OpenMP::OpenMP_C
        m
    )
    if(USE_CUDA_ACCEL)
        target_link_libraries(atmturb-benchmark-perf PRIVATE CUDA::cudart CUDA::cufft)
    endif()

    # ── Streaming performance benchmark executable ────────────────────
    add_executable(atmturb-benchmark-stream
        tests/benchmark_streaming.c
    )
    target_compile_definitions(atmturb-benchmark-stream PRIVATE MILK_NO_CLI)
    target_include_directories(atmturb-benchmark-stream PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
        ${CMAKE_CURRENT_SOURCE_DIR}/tests
    )
    target_include_directories(atmturb-benchmark-stream SYSTEM PRIVATE ${MILK_STANDALONE_INCLUDES})
    target_link_libraries(atmturb-benchmark-stream PRIVATE
        ${LIBNAME}
        ${MILK_LIBRARIES}
        ${FFTWF_LIBRARIES}
        ${FFTW_LIBRARIES}
        OpenMP::OpenMP_C
        m
    )
    if(USE_CUDA_ACCEL)
        target_link_libraries(atmturb-benchmark-stream PRIVATE CUDA::cudart CUDA::cufft)
    endif()
endif()
