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

# Find OpenMP (required for multi-threaded turbulence simulation)
find_package(OpenMP REQUIRED)

# Find GSL (GNU Scientific Library, for atmospheric modeling)
find_package(GSL REQUIRED)

# Find FFTW3 (for Fourier domain phase screens)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED fftw3)
pkg_check_modules(FFTWF REQUIRED fftw3f)

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
    src/AtmosphericTurbulence/atmturb_simd_scalar.c
    src/AtmosphericTurbulence/atmturb_simd_avx2.c
    src/AtmosphericTurbulence/atmturb_simd_avx512.c
    src/AtmosphericTurbulence/atmturb_simd_dispatch.c
    src/AtmosphericTurbulence/atmturb_mkwfs_FPS.c
    src/AtmosphericTurbulence/atmturb_mkhvturb_FPS.c
    src/AtmosphericTurbulence/atmturb_mkmastert_FPS.c
    src/AtmosphericTurbulence/atmturb_mkvonkarman_FPS.c
    src/OpticsMaterials/OpticsMaterials.c
    src/OpticsMaterials/optmat_solids.c
    src/OpticsMaterials/optmat_gases.c
    src/OpticsMaterials/optmat_pmgi.c
    src/OpticsMaterials/optmat_pmma.c
    src/WFpropagate/WFpropagate.c
    src/WFpropagate/wfprop_fresnel.c
    src/WFpropagate/wfprop_fresnel_FPS.c
    src/WFpropagate/wfprop_cube.c
    src/WFpropagate/wfprop_lyot.c
)

set(INCLUDEFILES
    src/milkatmturb.h
    src/AtmosphereModel/AtmosphereModel.h
    src/AtmosphericTurbulence/AtmosphericTurbulence.h
    src/AtmosphericTurbulence/atmturb_simd.h
    src/OpticsMaterials/OpticsMaterials.h
    src/WFpropagate/WFpropagate.h
)

file(GLOB SCRIPTS "src/AtmosphericTurbulence/scripts/*")

if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    set_source_files_properties(src/AtmosphericTurbulence/atmturb_simd_avx2.c
        PROPERTIES COMPILE_FLAGS "-mavx2 -mfma")
    set_source_files_properties(src/AtmosphericTurbulence/atmturb_simd_avx512.c
        PROPERTIES COMPILE_FLAGS "-mavx512f -mavx512dq -mavx512bw -mfma")
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
                src/WFpropagate/wfprop_fresnel_cuda.cu
                src/AtmosphericTurbulence/atmturb_cuda.cu
            )
            list(APPEND INCLUDEFILES
                src/WFpropagate/wfprop_fresnel_cuda.h
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
    target_include_directories(${LIBNAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphereModel
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence/nrlmsise-00.20131225
        ${CMAKE_CURRENT_SOURCE_DIR}/src/OpticsMaterials
        ${CMAKE_CURRENT_SOURCE_DIR}/src/WFpropagate
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/coremods
        ${CMAKE_SOURCE_DIR}/plugins/milk-extra-src
        $<TARGET_PROPERTY:CLIcore,INTERFACE_INCLUDE_DIRECTORIES>
    )

    target_link_libraries(${LIBNAME} PUBLIC
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
        GSL::gsl
        GSL::gslcblas
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
        ${CMAKE_CURRENT_SOURCE_DIR}/src/OpticsMaterials
        ${CMAKE_CURRENT_SOURCE_DIR}/src/WFpropagate
    )

    # Include compatibility shim for installed milk layout quirks
    set(MILK_SHIM_DIR "${CMAKE_CURRENT_BINARY_DIR}/milk_compat_include")
    file(MAKE_DIRECTORY "${MILK_SHIM_DIR}")
    if(MILK_INCLUDE_DIR)
        if(NOT EXISTS "${MILK_SHIM_DIR}/libprocessinfo")
            file(CREATE_LINK "${MILK_INCLUDE_DIR}" "${MILK_SHIM_DIR}/libprocessinfo" SYMBOLIC)
        endif()
        if(EXISTS "${MILK_INCLUDE_DIR}/libmilkcommon/milk_compiler.h" AND NOT EXISTS "${MILK_SHIM_DIR}/milk_compiler.h")
            file(CREATE_LINK "${MILK_INCLUDE_DIR}/libmilkcommon/milk_compiler.h" "${MILK_SHIM_DIR}/milk_compiler.h" SYMBOLIC)
        endif()
    endif()
    list(APPEND MILK_STANDALONE_INCLUDES "${MILK_SHIM_DIR}")

    if(MILK_SOURCE_FOUND)
        list(APPEND MILK_STANDALONE_INCLUDES
            ${MILK_SOURCE_FOUND}/src
            ${MILK_SOURCE_FOUND}/src/coremods
            ${MILK_SOURCE_FOUND}/plugins/milk-extra-src
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
        ${MILK_ROOT}/lib
        $ENV{MILK_ROOT}/lib
        /usr/local/milk/lib
        /usr/local/lib
    )
    if(MILK_SOURCE_FOUND)
        list(APPEND MILK_SEARCH_HINTS
            ${MILK_SOURCE_FOUND}/_build/lib
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

    set(MILK_LIBRARIES "")
    foreach(lib ${MILK_LIB_NAMES})
        find_library(MILK_${lib}_LIB NAMES ${lib} HINTS ${MILK_SEARCH_HINTS})
        if(MILK_${lib}_LIB)
            list(APPEND MILK_LIBRARIES ${MILK_${lib}_LIB})
        else()
            message(WARNING "MILK library '${lib}' not found in search paths; falling back to -l${lib}")
            list(APPEND MILK_LIBRARIES ${lib})
        endif()
    endforeach()

    target_link_libraries(${LIBNAME} PUBLIC
        ${MILK_LIBRARIES}
        OpenMP::OpenMP_C
        GSL::gsl
        GSL::gslcblas
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
if(COMMAND add_milk_standalone)
    add_milk_standalone(atmturb-mkwfs src/AtmosphericTurbulence/atmturb_mkwfs_FPS.c)
    target_include_directories(milk-fpsexec-atmturb-mkwfs PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
    )
    target_link_libraries(milk-fpsexec-atmturb-mkwfs PRIVATE ${LIBNAME})

    add_milk_standalone(atmturb-mkhvturb src/AtmosphericTurbulence/atmturb_mkhvturb_FPS.c)
    target_include_directories(milk-fpsexec-atmturb-mkhvturb PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
    )
    target_link_libraries(milk-fpsexec-atmturb-mkhvturb PRIVATE ${LIBNAME})

    add_milk_standalone(atmturb-mkmastert src/AtmosphericTurbulence/atmturb_mkmastert_FPS.c)
    target_include_directories(milk-fpsexec-atmturb-mkmastert PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
    )
    target_link_libraries(milk-fpsexec-atmturb-mkmastert PRIVATE ${LIBNAME})

    add_milk_standalone(atmturb-mkvonkarman src/AtmosphericTurbulence/atmturb_mkvonkarman_FPS.c)
    target_include_directories(milk-fpsexec-atmturb-mkvonkarman PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
    )
    target_link_libraries(milk-fpsexec-atmturb-mkvonkarman PRIVATE ${LIBNAME})

    add_milk_standalone(wfprop-fresnel src/WFpropagate/wfprop_fresnel_FPS.c)
    target_include_directories(milk-fpsexec-wfprop-fresnel PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/WFpropagate
    )
    target_link_libraries(milk-fpsexec-wfprop-fresnel PRIVATE ${LIBNAME})
endif()
