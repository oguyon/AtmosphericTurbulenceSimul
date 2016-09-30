bin_PROGRAMS = AtmosphericTurbulenceSimul
AtmosphericTurbulenceSimul_SOURCES = CLIcore.c CLIcore.h initmodules.c initmodules.h calc_bison.y calc_flex.l calc.h
AtmosphericTurbulenceSimul_CPPFLAGS = -I@abs_top_srcdir@/src -DCONFIGDIR=\"@abs_top_srcdir@/config\" -DSOURCEDIR=\"@abs_top_srcdir@/src\" -DABSSRCTOPDIR=\"@abs_top_srcdir@\" -fopenmp

ACLOCAL_AMFLAGS = -I m4

## For Bison and FLEX
AM_YFLAGS = -d
AM_LFLAGS = -i

noinst_HEADERS = calc.h
BUILT_SOURCES = calc_bison.h

SUBDIRS =
AtmosphericTurbulenceSimul_LDADD =



SUBDIRS += cudacomp
AtmosphericTurbulenceSimul_LDADD += cudacomp/libcudacomp.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/cudacomp

SUBDIRS += AtmosphericTurbulence
AtmosphericTurbulenceSimul_LDADD += AtmosphericTurbulence/libatmosphericturbulence.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/AtmosphericTurbulence

SUBDIRS += AtmosphereModel
AtmosphericTurbulenceSimul_LDADD += AtmosphereModel/libatmospheremodel.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/AtmosphereModel

SUBDIRS += psf
AtmosphericTurbulenceSimul_LDADD += psf/libpsf.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/psf

SUBDIRS += image_format
AtmosphericTurbulenceSimul_LDADD += image_format/libimageformat.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/image_format

SUBDIRS += OpticsMaterials
AtmosphericTurbulenceSimul_LDADD += OpticsMaterials/libopticsmaterials.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/OpticsMaterials

SUBDIRS += WFpropagate
AtmosphericTurbulenceSimul_LDADD += WFpropagate/libwfpropagate.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/WFpropagate

SUBDIRS += image_basic
AtmosphericTurbulenceSimul_LDADD += image_basic/libimagebasic.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/image_basic

SUBDIRS += image_filter
AtmosphericTurbulenceSimul_LDADD += image_filter/libimagefilter.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/image_filter

SUBDIRS += kdtree
AtmosphericTurbulenceSimul_LDADD += kdtree/libkdtree.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/kdtree

SUBDIRS += image_gen
AtmosphericTurbulenceSimul_LDADD += image_gen/libimagegen.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/image_gen

SUBDIRS += linopt_imtools
AtmosphericTurbulenceSimul_LDADD += linopt_imtools/liblinoptimtools.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/linopt_imtools

SUBDIRS += statistic
AtmosphericTurbulenceSimul_LDADD += statistic/libstatistic.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/statistic

SUBDIRS += fft
AtmosphericTurbulenceSimul_LDADD += fft/libfft.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/fft

SUBDIRS += info
AtmosphericTurbulenceSimul_LDADD += info/libinfo.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/info

SUBDIRS += COREMOD_arith
AtmosphericTurbulenceSimul_LDADD += COREMOD_arith/libcoremodarith.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_arith

SUBDIRS += COREMOD_iofits
AtmosphericTurbulenceSimul_LDADD += COREMOD_iofits/libcoremodiofits.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_iofits

SUBDIRS += COREMOD_memory
AtmosphericTurbulenceSimul_LDADD += COREMOD_memory/libcoremodmemory.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_memory

SUBDIRS += COREMOD_tools
AtmosphericTurbulenceSimul_LDADD += COREMOD_tools/libcoremodtools.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/COREMOD_tools

SUBDIRS += 00CORE
AtmosphericTurbulenceSimul_LDADD += 00CORE/lib00core.la
AtmosphericTurbulenceSimul_CPPFLAGS += -I@abs_top_srcdir@/src/00CORE

