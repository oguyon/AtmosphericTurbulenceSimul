# =====================================================================
# milkatmturb - Atmospheric Turbulence Simulation Plugin for milk
# =====================================================================

cmake_minimum_required(VERSION 3.16)

set(LIBNAME "milkatmturb")
set(SRCNAME "atmturb")

message(" SRCNAME = ${SRCNAME} -> LIBNAME = ${LIBNAME}")

# Project definition
project(lib_${LIBNAME}_project C)

# Set C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Default to Release build type for maximum optimization (-O3)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# RPATH Configuration (preserve runtime search path in installed binaries)
set(CMAKE_INSTALL_RPATH
    "${CMAKE_INSTALL_PREFIX}/lib;/usr/local/milk/lib;/usr/local/milk-1.03.00/lib;/usr/local/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

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
    src/AtmosphericTurbulence/atmturb_wfs_stream.c
    src/AtmosphericTurbulence/atmturb_fresnel.c
    src/AtmosphericTurbulence/atmturb_lowfreq.c
    src/AtmosphericTurbulence/atmturb_rolling.c
    src/AtmosphericTurbulence/atmturb_simd_scalar.c
    src/AtmosphericTurbulence/atmturb_simd_avx2.c
    src/AtmosphericTurbulence/atmturb_simd_avx512.c
    src/AtmosphericTurbulence/atmturb_simd_dispatch.c
    src/AtmosphericTurbulence/atmturb_mkwfs_FPS.c
    src/AtmosphericTurbulence/atmturb_mkhvturb_FPS.c
    src/AtmosphericTurbulence/atmturb_mkmastert_FPS.c
    src/AtmosphericTurbulence/atmturb_mkvonkarman_FPS.c
    src/AtmosphericTurbulence/atmturb_aoloop_FPS.c
)

set(INCLUDEFILES
    src/milkatmturb.h
    src/AtmosphereModel/AtmosphereModel.h
    src/AtmosphericTurbulence/AtmosphericTurbulence.h
    src/AtmosphericTurbulence/atmturb_profile.h
    src/AtmosphericTurbulence/atmturb_geometry.h
    src/AtmosphericTurbulence/atmturb_fresnel.h
    src/AtmosphericTurbulence/atmturb_lowfreq.h
    src/AtmosphericTurbulence/atmturb_rolling.h
    src/AtmosphericTurbulence/atmturb_simd.h
    src/AtmosphericTurbulence/atmturb_wfs_stream.h
)

file(GLOB SCRIPTS "src/AtmosphericTurbulence/scripts/*")

if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    set_source_files_properties(src/AtmosphericTurbulence/atmturb_simd_avx2.c
        PROPERTIES COMPILE_FLAGS "-mavx2 -mfma")
    set_source_files_properties(src/AtmosphericTurbulence/atmturb_simd_avx512.c
        PROPERTIES COMPILE_FLAGS "-mavx512f -mavx512dq -mavx512bw -mavx512vl -mfma")
endif()

# CUDA GPU Acceleration
if(HAVE_CUDA OR USE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA OPTIONAL)
        find_package(CUDAToolkit QUIET)
        if(CUDAToolkit_FOUND)
            set(USE_CUDA_ACCEL TRUE)
            add_compile_definitions(HAVE_CUDA)
            list(APPEND SOURCEFILES
                src/AtmosphericTurbulence/atmturb_cuda.cu
            )
            list(APPEND INCLUDEFILES
                src/AtmosphericTurbulence/atmturb_cuda.h
            )
            message(STATUS "milkatmturb: CUDA GPU acceleration ENABLED")
        endif()
    endif()
endif()

add_library(${LIBNAME} SHARED ${SOURCEFILES})

set_target_properties(${LIBNAME} PROPERTIES
    LINKER_LANGUAGE C
    INTERPROCEDURAL_OPTIMIZATION OFF
)

if(TARGET CLIcore)
    # In-tree build under milk framework
    set(MILK_STANDALONE_INCLUDES
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphereModel
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence/nrlmsise-00.20131225
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/coremods
        ${CMAKE_SOURCE_DIR}/plugins
        ${CMAKE_SOURCE_DIR}/plugins/milk-extra-src
        ${CMAKE_SOURCE_DIR}/plugins/OpticsMaterials
        ${CMAKE_SOURCE_DIR}/plugins/WFpropagate
        ${CMAKE_SOURCE_DIR}/src/cli
        ${CMAKE_SOURCE_DIR}/src/cli/CLIcore
        ${CMAKE_SOURCE_DIR}/src/engine
        ${CMAKE_SOURCE_DIR}/src/engine/libfps
        ${CMAKE_SOURCE_DIR}/src/engine/libmilkcommon
        ${CMAKE_SOURCE_DIR}/src/engine/libprocessinfo
        ${CMAKE_SOURCE_DIR}/src/engine/libmilkdata
        ${CMAKE_SOURCE_DIR}/src/engine/ImageStreamIO
        $<TARGET_PROPERTY:CLIcore,INTERFACE_INCLUDE_DIRECTORIES>
    )

    set(MILK_LIBRARIES "")

    target_include_directories(${LIBNAME} PRIVATE
        ${MILK_STANDALONE_INCLUDES}
    )

    target_link_libraries(${LIBNAME} PUBLIC
        milkOpticsMaterials
        milkWFpropagate
        CLIcore
        ImageStreamIO
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
        OpenMP::OpenMP_C
        ${FFTW_LIBRARIES}
        ${FFTWF_LIBRARIES}
        m
    )
else()
    # Out-of-tree / standalone build
    set(MILK_SOURCE_CANDIDATES
        ${MILK_SOURCE_DIR}
        $ENV{MILK_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/../milk-framework-dev
        $ENV{MILK_ROOT}
        ${CMAKE_CURRENT_SOURCE_DIR}/../milk
    )
    set(MILK_SOURCE_FOUND "")
    foreach(cand ${MILK_SOURCE_CANDIDATES})
        if(cand AND EXISTS "${cand}/src/coremods" AND EXISTS "${cand}/plugins/milk-extra-src")
            get_filename_component(MILK_SOURCE_FOUND "${cand}" ABSOLUTE)
            message(STATUS "Found MILK source tree: ${MILK_SOURCE_FOUND}")
            break()
        endif()
    endforeach()

    find_path(MILK_INCLUDE_DIR
        NAMES CommandLineInterface/CLIcore.h CLIcore.h
        HINTS
            ${MILK_ROOT}/include
            $ENV{MILK_ROOT}/include
            /usr/local/milk/include
            /usr/local/include
    )

    set(MILK_STANDALONE_INCLUDES
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphereModel
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence/nrlmsise-00.20131225
        ${CMAKE_CURRENT_SOURCE_DIR}/../OpticsMaterials
        ${CMAKE_CURRENT_SOURCE_DIR}/../WFpropagate
    )

    # Include compatibility shim for installed milk layout quirks
    set(MILK_SHIM_DIR "${CMAKE_CURRENT_BINARY_DIR}/milk_compat_include")
    file(MAKE_DIRECTORY "${MILK_SHIM_DIR}")
    if(MILK_INCLUDE_DIR)
        if(NOT EXISTS "${MILK_SHIM_DIR}/libprocessinfo")
            file(CREATE_LINK "${MILK_INCLUDE_DIR}" "${MILK_SHIM_DIR}/libprocessinfo" SYMBOLIC)
        endif()
        if(EXISTS "${MILK_INCLUDE_DIR}/libmilkcommon/milk_compiler.h"
           AND NOT EXISTS "${MILK_SHIM_DIR}/milk_compiler.h")
            file(CREATE_LINK "${MILK_INCLUDE_DIR}/libmilkcommon/milk_compiler.h"
                 "${MILK_SHIM_DIR}/milk_compiler.h" SYMBOLIC)
        endif()
    endif()
    list(APPEND MILK_STANDALONE_INCLUDES "${MILK_SHIM_DIR}")

    if(MILK_SOURCE_FOUND)
        list(APPEND MILK_STANDALONE_INCLUDES
            ${MILK_SOURCE_FOUND}/src
            ${MILK_SOURCE_FOUND}/src/coremods
            ${MILK_SOURCE_FOUND}/plugins/milk-extra-src
            ${MILK_SOURCE_FOUND}/plugins/OpticsMaterials
            ${MILK_SOURCE_FOUND}/plugins/WFpropagate
            ${MILK_SOURCE_FOUND}/plugins/milk-extra-src/OpticsMaterials
            ${MILK_SOURCE_FOUND}/plugins/milk-extra-src/WFpropagate
            ${MILK_SOURCE_FOUND}/src/cli
            ${MILK_SOURCE_FOUND}/src/cli/CLIcore
            ${MILK_SOURCE_FOUND}/src/engine
            ${MILK_SOURCE_FOUND}/src/engine/libfps
            ${MILK_SOURCE_FOUND}/src/engine/libmilkcommon
            ${MILK_SOURCE_FOUND}/src/engine/libprocessinfo
            ${MILK_SOURCE_FOUND}/src/engine/libmilkdata
            ${MILK_SOURCE_FOUND}/src/engine/ImageStreamIO
        )
        if(EXISTS "${MILK_SOURCE_FOUND}/_build/config")
            list(APPEND MILK_STANDALONE_INCLUDES "${MILK_SOURCE_FOUND}/_build/config")
        endif()
    endif()

    if(MILK_INCLUDE_DIR)
        list(APPEND MILK_STANDALONE_INCLUDES
            ${MILK_INCLUDE_DIR}
            ${MILK_INCLUDE_DIR}/CommandLineInterface
            ${MILK_INCLUDE_DIR}/CommandLineInterface/CLIcore
            ${MILK_INCLUDE_DIR}/libmilkcommon
            ${MILK_INCLUDE_DIR}/libfps
            ${MILK_INCLUDE_DIR}/libmilkdata
            ${MILK_INCLUDE_DIR}/ImageStreamIO
        )
    endif()

    target_include_directories(${LIBNAME} PRIVATE ${MILK_STANDALONE_INCLUDES})

    set(MILK_LIB_NAMES
        milkOpticsMaterials
        milkWFpropagate
        CLIcore
        ImageStreamIO
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
        ${CMAKE_CURRENT_SOURCE_DIR}/../OpticsMaterials/_build
        ${CMAKE_CURRENT_SOURCE_DIR}/../WFpropagate/_build
    )
    if(MILK_SOURCE_FOUND)
        list(APPEND MILK_SEARCH_HINTS
            ${MILK_SOURCE_FOUND}/_build/lib
            ${MILK_SOURCE_FOUND}/_build/plugins/OpticsMaterials
            ${MILK_SOURCE_FOUND}/_build/plugins/WFpropagate
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/OpticsMaterials
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/WFpropagate
            ${MILK_SOURCE_FOUND}/_build/src/cli/CLIcore
            ${MILK_SOURCE_FOUND}/_build/src/engine/ImageStreamIO
            ${MILK_SOURCE_FOUND}/_build/src/coremods/COREMOD_memory
            ${MILK_SOURCE_FOUND}/_build/src/coremods/COREMOD_arith
            ${MILK_SOURCE_FOUND}/_build/src/coremods/COREMOD_iofits
            ${MILK_SOURCE_FOUND}/_build/src/coremods/COREMOD_tools
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/fft
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/psf
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/linalgebra
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/image_gen
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/image_filter
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/image_basic
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/statistic
            ${MILK_SOURCE_FOUND}/_build/plugins/milk-extra-src/linopt_imtools
        )
    endif()
    list(APPEND MILK_SEARCH_HINTS
        ${MILK_ROOT}/lib
        $ENV{MILK_ROOT}/lib
        /usr/local/milk/lib
        /usr/local/lib
    )


    set(MILK_LIBRARIES "")
    foreach(lib ${MILK_LIB_NAMES})
        find_library(MILK_${lib}_LIB NAMES ${lib} HINTS ${MILK_SEARCH_HINTS})
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
endif()

if(USE_CUDA_ACCEL)
    target_link_libraries(${LIBNAME} PRIVATE CUDA::cudart CUDA::cufft)
endif()

# Export and install shared library
install(
    TARGETS ${LIBNAME}
    EXPORT milkTargets
    DESTINATION lib
)

install(FILES ${INCLUDEFILES} DESTINATION include/${SRCNAME})
install(PROGRAMS ${SCRIPTS} DESTINATION bin)

# ── Standalone executables (FPS V2 pattern) ───────────────────────
if(NOT COMMAND add_milk_standalone AND MILK_SOURCE_FOUND
   AND EXISTS "${MILK_SOURCE_FOUND}/cmake/milk_standalone.cmake")
    if(EXISTS "${MILK_SOURCE_FOUND}/cmake/milk_extensions.cmake")
        include("${MILK_SOURCE_FOUND}/cmake/milk_extensions.cmake")
    endif()
    include("${MILK_SOURCE_FOUND}/cmake/milk_standalone.cmake")
endif()

if(COMMAND add_milk_standalone)
    link_directories(/usr/local/milk/lib ${MILK_SEARCH_HINTS})
    add_milk_standalone(atmturb-mkwfs src/AtmosphericTurbulence/atmturb_mkwfs_FPS.c)
    target_include_directories(milk-fpsexec-atmturb-mkwfs PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
        ${MILK_STANDALONE_INCLUDES}
    )
    target_link_libraries(milk-fpsexec-atmturb-mkwfs PRIVATE ${LIBNAME} ${MILK_LIBRARIES})

    add_milk_standalone(atmturb-mkhvturb src/AtmosphericTurbulence/atmturb_mkhvturb_FPS.c)
    target_include_directories(milk-fpsexec-atmturb-mkhvturb PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
        ${MILK_STANDALONE_INCLUDES}
    )
    target_link_libraries(milk-fpsexec-atmturb-mkhvturb PRIVATE ${LIBNAME} ${MILK_LIBRARIES})

    add_milk_standalone(atmturb-mkmastert src/AtmosphericTurbulence/atmturb_mkmastert_FPS.c)
    target_include_directories(milk-fpsexec-atmturb-mkmastert PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
        ${MILK_STANDALONE_INCLUDES}
    )
    target_link_libraries(milk-fpsexec-atmturb-mkmastert PRIVATE ${LIBNAME} ${MILK_LIBRARIES})

    add_milk_standalone(atmturb-mkvonkarman src/AtmosphericTurbulence/atmturb_mkvonkarman_FPS.c)
    target_include_directories(milk-fpsexec-atmturb-mkvonkarman PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
        ${MILK_STANDALONE_INCLUDES}
    )
    target_link_libraries(milk-fpsexec-atmturb-mkvonkarman PRIVATE ${LIBNAME} ${MILK_LIBRARIES})

    add_milk_standalone(atmturb-aoloop src/AtmosphericTurbulence/atmturb_aoloop_FPS.c)
    target_include_directories(milk-fpsexec-atmturb-aoloop PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
        ${MILK_STANDALONE_INCLUDES}
    )
    target_link_libraries(milk-fpsexec-atmturb-aoloop PRIVATE ${LIBNAME} ${MILK_LIBRARIES})
endif()

# ── Statistical test validation executable ─────────────────────────
add_executable(atmturb-validate-stats
    tests/atmturb_validate_stats.c
    tests/atmturb_validate_math.c
)
target_include_directories(atmturb-validate-stats PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
    ${CMAKE_CURRENT_SOURCE_DIR}/tests
    ${MILK_STANDALONE_INCLUDES}
)
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
target_include_directories(atmturb-benchmark-perf PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
    ${CMAKE_CURRENT_SOURCE_DIR}/tests
    ${MILK_STANDALONE_INCLUDES}
)
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

