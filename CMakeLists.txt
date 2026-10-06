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
    src/AtmosphericTurbulence/AtmosphericTurbulence.c
    src/AtmosphericTurbulence/atmturb_mkwfs_FPS.c
    src/AtmosphericTurbulence/atmturb_mkhvturb_FPS.c
    src/OpticsMaterials/OpticsMaterials.c
    src/WFpropagate/WFpropagate.c
)

set(INCLUDEFILES
    src/milkatmturb.h
    src/AtmosphereModel/AtmosphereModel.h
    src/AtmosphericTurbulence/AtmosphericTurbulence.h
    src/OpticsMaterials/OpticsMaterials.h
    src/WFpropagate/WFpropagate.h
)

file(GLOB SCRIPTS "src/AtmosphericTurbulence/scripts/*")

# Build shared library for milk CLI
add_library(${LIBNAME} SHARED ${SOURCEFILES})

target_include_directories(${LIBNAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphereModel
    ${CMAKE_CURRENT_SOURCE_DIR}/src/AtmosphericTurbulence
    ${CMAKE_CURRENT_SOURCE_DIR}/src/OpticsMaterials
    ${CMAKE_CURRENT_SOURCE_DIR}/src/WFpropagate
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/coremods
    ${CMAKE_SOURCE_DIR}/plugins/milk-extra-src
)

target_include_directories(${LIBNAME} PRIVATE
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
endif()
