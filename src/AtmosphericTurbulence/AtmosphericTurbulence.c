// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    AtmosphericTurbulence.c
 * @brief   Atmospheric turbulence module lifecycle and CLI command registration
 */

#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

extern DATA data;

static int make_AtmosphericTurbulence_wavefront_series_cli(void)
{
    if (CLI_checkarg(1, 1) + CLI_checkarg(2, 2) == 0)
    {
        make_AtmosphericTurbulence_wavefront_series(data.cmdargtoken[1].val.numf,
                                                   data.cmdargtoken[2].val.numl);
        return 0;
    }
    return 1;
}

static int make_AtmosphericTurbulence_vonKarmanWind_cli(void)
{
    if (CLI_checkarg(1, 2) + CLI_checkarg(2, 1) + CLI_checkarg(3, 1) +
        CLI_checkarg(4, 1) + CLI_checkarg(5, 2) + CLI_checkarg(6, 3) == 0)
    {
        make_AtmosphericTurbulence_vonKarmanWind(
            data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numf,
            data.cmdargtoken[3].val.numf, data.cmdargtoken[4].val.numf,
            data.cmdargtoken[5].val.numl, data.cmdargtoken[6].val.string);
        return 0;
    }
    return 1;
}

static int AtmosphericTurbulence_mkmastert_cli(void)
{
    if (CLI_checkarg(1, 3) + CLI_checkarg(2, 3) + CLI_checkarg(3, 2) +
        CLI_checkarg(4, 1) + CLI_checkarg(5, 1) == 0)
    {
        make_master_turbulence_screen(
            data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.numl, (float)data.cmdargtoken[4].val.numf,
            (float)data.cmdargtoken[5].val.numf, 0);
        return 0;
    }
    return 1;
}

static int AtmosphericTurbulence_makeHV_CN2prof_cli(void)
{
    if (CLI_checkarg(1, 1) + CLI_checkarg(2, 1) + CLI_checkarg(3, 1) +
        CLI_checkarg(4, 2) + CLI_checkarg(5, 3) == 0)
    {
        AtmosphericTurbulence_makeHV_CN2prof(
            data.cmdargtoken[1].val.numf, data.cmdargtoken[2].val.numf,
            data.cmdargtoken[3].val.numf, data.cmdargtoken[4].val.numl,
            data.cmdargtoken[5].val.string);
        return 0;
    }
    return 1;
}

static int AtmosphericTurbulence_measure_wavefront_series_expoframes_cli(void)
{
    if (CLI_checkarg(1, 1) + CLI_checkarg(2, 3) == 0)
    {
        measure_wavefront_series_expoframes(data.cmdargtoken[1].val.numf,
                                            data.cmdargtoken[2].val.string);
        return 0;
    }
    return 1;
}

static int AtmosphericTurbulence_mkTestTTseq_cli(void)
{
    if (CLI_checkarg(1, 1) + CLI_checkarg(2, 2) + CLI_checkarg(3, 2) +
        CLI_checkarg(4, 1) + CLI_checkarg(5, 2) + CLI_checkarg(6, 1) +
        CLI_checkarg(7, 2) == 0)
    {
        AtmosphericTurbulence_mkTestTTseq(
            data.cmdargtoken[1].val.numf, data.cmdargtoken[2].val.numl,
            data.cmdargtoken[3].val.numl, data.cmdargtoken[4].val.numf,
            data.cmdargtoken[5].val.numl, data.cmdargtoken[6].val.numf,
            data.cmdargtoken[7].val.numl);
        return 0;
    }
    return 1;
}

static int AtmosphericTurbulence_Build_LinPredictor_Full_cli(void)
{
    if (CLI_checkarg(1, 4) + CLI_checkarg(2, 4) + CLI_checkarg(3, 2) +
        CLI_checkarg(4, 1) + CLI_checkarg(5, 1) + CLI_checkarg(6, 1) == 0)
    {
        AtmosphericTurbulence_Build_LinPredictor_Full(
            data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.numl, data.cmdargtoken[4].val.numf,
            data.cmdargtoken[5].val.numf, data.cmdargtoken[6].val.numf);
        return 0;
    }
    return 1;
}

static int AtmosphericTurbulence_Apply_LinPredictor_Full_cli(void)
{
    if (CLI_checkarg(1, 2) + CLI_checkarg(2, 4) + CLI_checkarg(3, 4) +
        CLI_checkarg(4, 2) + CLI_checkarg(5, 1) + CLI_checkarg(6, 3) +
        CLI_checkarg(7, 3) == 0)
    {
        AtmosphericTurbulence_Apply_LinPredictor_Full(
            data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.string, data.cmdargtoken[4].val.numl,
            data.cmdargtoken[5].val.numf, data.cmdargtoken[6].val.string,
            data.cmdargtoken[7].val.string);
        return 0;
    }
    return 1;
}

static int AtmosphericTurbulence_LinPredictor_filt_2DKernelExtract_cli(void)
{
    if (CLI_checkarg(1, 4) + CLI_checkarg(2, 4) + CLI_checkarg(3, 2) +
        CLI_checkarg(4, 3) == 0)
    {
        AtmosphericTurbulence_LinPredictor_filt_2DKernelExtract(
            data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.numl, data.cmdargtoken[4].val.string);
        return 0;
    }
    return 1;
}

static int AtmosphericTurbulence_LinPredictor_filt_Expand_cli(void)
{
    if (CLI_checkarg(1, 4) + CLI_checkarg(2, 4) == 0)
    {
        AtmosphericTurbulence_LinPredictor_filt_Expand(
            data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string);
        return 0;
    }
    return 1;
}

static int AtmosphericTurbulence_Build_LinPredictor_cli(void)
{
    if (CLI_checkarg(1, 2) + CLI_checkarg(2, 1) + CLI_checkarg(3, 2) +
        CLI_checkarg(4, 2) + CLI_checkarg(5, 2) + CLI_checkarg(6, 2) +
        CLI_checkarg(7, 2) + CLI_checkarg(8, 1) == 0)
    {
        AtmosphericTurbulence_Build_LinPredictor(
            data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numf,
            data.cmdargtoken[3].val.numl, data.cmdargtoken[4].val.numl,
            data.cmdargtoken[5].val.numl, data.cmdargtoken[6].val.numl,
            data.cmdargtoken[7].val.numl, data.cmdargtoken[8].val.numf);
        return 0;
    }
    return 1;
}

static int AtmosphericTurbulence_psfCubeContrast_cli(void)
{
    if (CLI_checkarg(1, 4) + CLI_checkarg(2, 4) + CLI_checkarg(3, 3) == 0)
    {
        AtmosphericTurbulence_psfCubeContrast(
            data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.string);
        return 0;
    }
    return 1;
}

static int AtmosphericTurbulence_Test_LinPredictor_cli(void)
{
    if (CLI_checkarg(1, 2) + CLI_checkarg(2, 1) + CLI_checkarg(3, 4) +
        CLI_checkarg(4, 2) + CLI_checkarg(5, 2) + CLI_checkarg(6, 2) +
        CLI_checkarg(7, 1) == 0)
    {
        AtmosphericTurbulence_Test_LinPredictor(
            data.cmdargtoken[1].val.numl, data.cmdargtoken[2].val.numf,
            data.cmdargtoken[3].val.string, data.cmdargtoken[4].val.numl,
            data.cmdargtoken[5].val.numl, data.cmdargtoken[6].val.numl,
            data.cmdargtoken[7].val.numf);
        return 0;
    }
    return 1;
}

/**
 * init_AtmosphericTurbulence - Initialize module and register CLI commands
 *
 * Return: 0 on success.
 */
int init_AtmosphericTurbulence(void)
{
    RegisterCLIcommand("mkwfs", __FILE__, make_AtmosphericTurbulence_wavefront_series_cli,
                       "make wavefront series",
                       "<wavelength [nm]> <precision 0=single, 1=double>",
                       "mkwfs 1650.0 1",
                       "int make_AtmosphericTurbulence_wavefront_series(float slambdaum, long WFprecision)");

    RegisterCLIcommand("mkvonKarmanWind", __FILE__, make_AtmosphericTurbulence_vonKarmanWind_cli,
                       "make vonKarman wind model",
                       "<pixsize> <pixscale [m/pix]> <sigma windspeed [m/s]> <scale [m]> <size [long]> <output name>",
                       "mkvonKarmanWind 8192 0.1 20.0 50.0 512 vKmodel",
                       "long make_AtmosphericTurbulence_vonKarmanWind(long vKsize, float pixscale, float sigmawind, float Lwind, long size, char *IDout_name)");

    RegisterCLIcommand("mkmastert", __FILE__, AtmosphericTurbulence_mkmastert_cli,
                       "make 2 master phase screens",
                       "<screen0> <screen1> <size> <outerscale> <innerscale>",
                       "mkmastert scr0 scr1 2048 50.0 2.0",
                       "int make_master_turbulence_screen(char *ID_name1, char *ID_name2, long size, float outercale, float innercale, long WFprecision)");

    RegisterCLIcommand("mkHVturbprof", __FILE__, AtmosphericTurbulence_makeHV_CN2prof_cli,
                       "make Hufnagel-Valley turbulence profile",
                       "<high wind speed [m/s]> <r0 [m]> <site alt [m]> <NBlayers> <output file>",
                       "mkHVturbprof 21.0 0.15 4200 100 turbHV.prof",
                       "int AtmosphericTurbulence_makeHV_CN2prof(double wspeed, double r0, double sitealt, long NBlayer, char *outfile)");

    RegisterCLIcommand("atmturbmeasexpo", __FILE__, AtmosphericTurbulence_measure_wavefront_series_expoframes_cli,
                       "Measure long exposure time PSF from wavefront series",
                       "<etime [s]> <out name>",
                       "atmturbmeasexpo 1.0 outpsf",
                       "int measure_wavefront_series_expoframes(float etime, char *outfile)");

    RegisterCLIcommand("atmturbmktestTTs", __FILE__, AtmosphericTurbulence_mkTestTTseq_cli,
                       "make test TT sequence",
                       "<dt [s]> <number of pts per block> <number of blocks> <measurement noise> <accelerometer mode> <accelerometer noise> <mode>",
                       "atmturbmktestTTs 0.001 1000 10 0.1 0 0.0 0",
                       "int AtmosphericTurbulence_mkTestTTseq(double dt, long NBpts, long NBblocks, double measnoise, int ACCnmode, double ACCnoise, int MODE)");

    RegisterCLIcommand("atmturbwfpredictf", __FILE__, AtmosphericTurbulence_Build_LinPredictor_Full_cli,
                       "build full linear predictor from wavefront series",
                       "<input WF series (cube)> <mask image> <predictor order> <predictor time lag> <SVD eps> <RegLambda>",
                       "atmturbwfpredictf wfin wfmask 20 3.5 0.001 0.0",
                       "int AtmosphericTurbulence_Build_LinPredictor_Full(char *WFin_name, char *WFmask_name, int PForder, float PFlag, double SVDeps, double Rlambda)");

    RegisterCLIcommand("atmturbwfpapply", __FILE__, AtmosphericTurbulence_Apply_LinPredictor_Full_cli,
                       "Apply full linear predictor from wavefront series",
                       "<mode> <input WF series (cube)> <mask image> <predictor order> <predictor time lag> <predicted future values> <measured future values>",
                       "atmturbwfpapply 0 wfin wfmask 20 3.5 outp outf",
                       "int AtmosphericTurbulence_Apply_LinPredictor_Full(int MODE, char *WFin_name, char *WFmask_name, int PForder, float PFlag, char *WFoutp_name, char *WFoutf_name)");

    RegisterCLIcommand("atmturbwfp2Dkern", __FILE__, AtmosphericTurbulence_LinPredictor_filt_2DKernelExtract_cli,
                       "collapse WF predictor into 2D kernel",
                       "<input WF filter (cube)> <mask image> <kernel radius> <output kernel name>",
                       "atmturbwfp2Dkern wfpfilt wfmask 20 wfpkern",
                       "long AtmosphericTurbulence_LinPredictor_filt_2DKernelExtract(char *IDfilt_name, char *IDmask_name, long krad, char *IDkern_name)");

    RegisterCLIcommand("atmturbwfpexp", __FILE__, AtmosphericTurbulence_LinPredictor_filt_Expand_cli,
                       "Expand 3D filter cube into pixel-based 3D cube filters",
                       "<input WF filter (cube)> <mask image>",
                       "atmturbwfpexp wfpfilt wfmask ",
                       "long AtmosphericTurbulence_LinPredictor_filt_Expand(char *IDfilt_name, char *IDmask_name)");

    RegisterCLIcommand("atmturbwfpredict", __FILE__, AtmosphericTurbulence_Build_LinPredictor_cli,
                       "build linear predictor from wavefront series",
                       "<number steps input> <noise level [rad]> <predictor z size> <predictor xy radius> <lambda [um]>",
                       "atmturbwfpredict 1000 0.01 5 5 64 64 1.65",
                       "int AtmosphericTurbulence_Build_LinPredictor(long NB_WFstep, double WFphaNoise, long WFP_NBstep, long WFP_xyrad, long WFPiipix, long WFPjjpix, float slambdaum)");

    RegisterCLIcommand("atmturbmkpsfcc", __FILE__, AtmosphericTurbulence_psfCubeContrast_cli,
                       "measure contrast performance of WF cube",
                       "<input WF cube> <mask> <output psf cube>",
                       "atmturbmkpsfcc wfc mask psfc",
                       "long AtmosphericTurbulence_psfCubeContrast(char *IDwfc_name, char *IDmask_name, char *IDpsfc_name)");

    RegisterCLIcommand("atmturbwfptest", __FILE__, AtmosphericTurbulence_Test_LinPredictor_cli,
                       "Test linear predictor on wavefront series",
                       "<number steps input> <noise level [rad]> <predictor name> <lag> <iipix> <jjpix>",
                       "atmturbwfptest 1000 0.01 wfpfilt 1 32 54",
                       "int AtmosphericTurbulence_Test_LinPredictor(long NB_WFstep, double WFphaNoise, char *IDWFPfilt_name, long WFPlag, long WFPiipix, long WFPjjpix)");

    return 0;
}
