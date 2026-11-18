// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_cli.c
 * @brief   MILK_RUNNABLE CLI command wrappers for AtmosphereModel and AtmosphericTurbulence
 *
 * Defines interactive CLI commands `atmturb.mkatm`, `atmturb.measr0` (exporting
 * calculator variable `r0pix`), `atmturb.linpred_build`, `atmturb.linpred_apply`,
 * and `atmturb.contract`.
 */

#include <stdint.h>
#include <stdio.h>

#include "milk_runnable.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "AtmosphereModel/AtmosphereModel.h"
#include "AtmosphericTurbulence/AtmosphericTurbulence.h"
#include "atmturb_cli.h"
#include "atmturb_compat.h"

/* 1. mkatm: build NRLMSISE-00 atmosphere model and refraction paths */
static char   *mkatm_conf;
static double *mkatm_slambda;

static CLICMDARGDEF farg_mkatm[] = {
    { CLIARG_STR, ".conffile", "Atmosphere configuration file", "WFsim.conf",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT | FPFLAG_CLI_OPTIONAL),
      (void **) &mkatm_conf, NULL },
    { CLIARG_FLOAT64, ".slambda", "Secondary wavelength [um or m]", "1.65",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT | FPFLAG_CLI_OPTIONAL),
      (void **) &mkatm_slambda, NULL },
};

/**
 * cli_atmturb_mkatm - Build standard atmosphere profile and refraction paths from config
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE on error.
 */
static errno_t cli_atmturb_mkatm(void)
{
    double lam = *mkatm_slambda;
    if (lam > 20.0)
    {
        lam *= 1e-9;
    }
    else if (lam > 1e-3)
    {
        lam *= 1e-6;
    }

    if (AtmosphereModel_Create_from_CONF(mkatm_conf, (float) lam) != 0)
    {
        return RETURN_FAILURE;
    }
    return RETURN_SUCCESS;
}

MILK_RUNNABLE_RT_DEFINE(atmturb_mkatm);
const MILK_RUNNABLE runnable_atmturb_mkatm = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = {
        .cmdkey           = "mkatm",
        .description      = "Build NRLMSISE-00 atmosphere model and refraction paths",
        .description_long = "Reads site meteorology from .conffile, generates vertical "
                            "14-species NRLMSISE-00 profile `atm.txt`, writes site dispersion "
                            "`RindexSite.txt`, and traces curved refraction paths at 0.55 um "
                            "and .slambda."
    },
    .srcfile  = __FILE__,
    .farg     = farg_mkatm,
    .nbarg    = (int) (sizeof(farg_mkatm) / sizeof(farg_mkatm[0])),
    .cli_func = cli_atmturb_mkatm,
    .rt       = &atmturb_mkatm_rt,
};

/* 2. measr0: measure effective Fried parameter r0 [pixels] on a 2D/3D phase image */
static char *measr0_imname;

static CLICMDARGDEF farg_measr0[] = {
    { CLIARG_IMG, ".imname", "Input 2D/3D float phase image or stream", "outarraypha",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
      (void **) &measr0_imname, NULL },
};

/**
 * cli_atmturb_measr0 - Estimate Fried parameter r0 in pixels from phase structure function
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE on invalid image.
 */
static errno_t cli_atmturb_measr0(void)
{
    imageID id = image_ID(measr0_imname);
    if (id < 0 || dcimg[id].used != 1)
    {
        printf("ERROR: Image \"%s\" not found.\n", measr0_imname);
        return RETURN_FAILURE;
    }

    uint8_t naxis = dcimg[id].md[0].naxis;
    long    nx    = (long) dcimg[id].md[0].size[0];
    long    ny    = (naxis >= 2) ? (long) dcimg[id].md[0].size[1] : 0;
    long    nz    = (naxis >= 3) ? (long) dcimg[id].md[0].size[2] : 1;

    if (naxis < 2 || nx != ny || nx < 2 || dcimg[id].md[0].datatype != _DATATYPE_FLOAT)
    {
        printf("ERROR: Image \"%s\" must be square 2D/3D float32.\n", measr0_imname);
        return RETURN_FAILURE;
    }

    double sum_r0   = 0.0;
    long   slice_px = nx * ny;
    for (long z = 0; z < nz; z++)
    {
        sum_r0 += atmturb_measure_r0_pix(&dcimg[id].array.F[z * slice_px], nx);
    }

    double r0_pix = sum_r0 / (double) nz;
    create_variable_ID("r0pix", r0_pix);
    printf("[milkatmturb] %s: r0_pix = %.4f (%ld frame%s)\n",
           measr0_imname, r0_pix, nz, (nz == 1) ? "" : "s");
    fflush(stdout);
    return RETURN_SUCCESS;
}

MILK_RUNNABLE_RT_DEFINE(atmturb_measr0);
const MILK_RUNNABLE runnable_atmturb_measr0 = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = {
        .cmdkey           = "measr0",
        .description      = "Estimate Fried parameter r0 [pix] from 2D/3D phase image",
        .description_long = "Computes the mean lag-1 phase structure function D_phi(1 px) "
                            "across a square 2D or 3D float32 phase image/cube and exports "
                            "the inferred Kolmogorov Fried parameter (6.88 / D_phi)^(3/5) to "
                            "calculator variable `r0pix`.",
        .outvars          = "r0pix"
    },
    .srcfile  = __FILE__,
    .farg     = farg_measr0,
    .nbarg    = (int) (sizeof(farg_measr0) / sizeof(farg_measr0[0])),
    .cli_func = cli_atmturb_measr0,
    .rt       = &atmturb_measr0_rt,
};

/* 3. linpred_build: build full-aperture AR linear predictor matrix */
static char    *lpbuild_wfin;
static char    *lpbuild_mask;
static int64_t *lpbuild_order;
static double  *lpbuild_lag;
static double  *lpbuild_svdeps;
static double  *lpbuild_reglambda;

static CLICMDARGDEF farg_lpbuild[] = {
    { CLIARG_IMG, ".wfin", "Input 3D wavefront phase cube", "outarraypha",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &lpbuild_wfin, NULL },
    { CLIARG_IMG, ".wfmask", "2D pupil mask image", "wfmask",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &lpbuild_mask, NULL },
    { CLIARG_INT64, ".order", "Autoregressive history order", "3",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &lpbuild_order, NULL },
    { CLIARG_FLOAT64, ".lag", "Prediction horizon lag [frames]", "2.0",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &lpbuild_lag, NULL },
    { CLIARG_FLOAT64, ".svdeps", "SVD cutoff tolerance", "0.01",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT | FPFLAG_CLI_OPTIONAL),
      (void **) &lpbuild_svdeps, NULL },
    { CLIARG_FLOAT64, ".reglambda", "Tikhonov regularization parameter", "0.0",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT | FPFLAG_CLI_OPTIONAL),
      (void **) &lpbuild_reglambda, NULL },
};

/**
 * cli_atmturb_linpred_build - Train full-aperture linear wavefront predictor
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE on error.
 */
static errno_t cli_atmturb_linpred_build(void)
{
    int rc = AtmosphericTurbulence_Build_LinPredictor_Full(
        lpbuild_wfin, lpbuild_mask, (int) *lpbuild_order, (float) *lpbuild_lag,
        *lpbuild_svdeps, *lpbuild_reglambda);
    return (rc == 0) ? RETURN_SUCCESS : RETURN_FAILURE;
}

MILK_RUNNABLE_RT_DEFINE(atmturb_linpred_build);
const MILK_RUNNABLE runnable_atmturb_linpred_build = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = {
        .cmdkey           = "linpred_build",
        .description      = "Train full-aperture autoregressive linear wavefront predictor",
        .description_long = "Builds empirical spatiotemporal covariance matrices over active "
                            "pupil pixels in .wfmask and solves for the regularized SVD "
                            "spatiotemporal prediction filter `outPF`."
    },
    .srcfile  = __FILE__,
    .farg     = farg_lpbuild,
    .nbarg    = (int) (sizeof(farg_lpbuild) / sizeof(farg_lpbuild[0])),
    .cli_func = cli_atmturb_linpred_build,
    .rt       = &atmturb_linpred_build_rt,
};

/* 4. linpred_apply: apply trained full-aperture linear predictor */
static int64_t *lpapply_mode;
static char    *lpapply_wfin;
static char    *lpapply_mask;
static int64_t *lpapply_order;
static double  *lpapply_lag;
static char    *lpapply_outp;
static char    *lpapply_outf;

static CLICMDARGDEF farg_lpapply[] = {
    { CLIARG_INT64, ".mode", "Predictor mode flag", "0",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &lpapply_mode, NULL },
    { CLIARG_IMG, ".wfin", "Input 3D wavefront phase cube", "outarraypha",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &lpapply_wfin, NULL },
    { CLIARG_IMG, ".wfmask", "2D pupil mask image", "wfmask",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &lpapply_mask, NULL },
    { CLIARG_INT64, ".order", "Autoregressive history order", "3",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &lpapply_order, NULL },
    { CLIARG_FLOAT64, ".lag", "Prediction horizon lag [frames]", "2.0",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &lpapply_lag, NULL },
    { CLIARG_STR_NOT_IMG, ".wfoutp", "Output predicted wavefront cube name", "wfoutp",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT | FPFLAG_CLI_OPTIONAL),
      (void **) &lpapply_outp, NULL },
    { CLIARG_STR_NOT_IMG, ".wfoutf", "Output residual error cube name", "wfoutf",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT | FPFLAG_CLI_OPTIONAL),
      (void **) &lpapply_outf, NULL },
};

/**
 * cli_atmturb_linpred_apply - Apply full-aperture linear wavefront predictor
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE on error.
 */
static errno_t cli_atmturb_linpred_apply(void)
{
    int rc = AtmosphericTurbulence_Apply_LinPredictor_Full(
        (int) *lpapply_mode, lpapply_wfin, lpapply_mask, (int) *lpapply_order,
        (float) *lpapply_lag, lpapply_outp, lpapply_outf);
    return (rc == 0) ? RETURN_SUCCESS : RETURN_FAILURE;
}

MILK_RUNNABLE_RT_DEFINE(atmturb_linpred_apply);
const MILK_RUNNABLE runnable_atmturb_linpred_apply = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = {
        .cmdkey           = "linpred_apply",
        .description      = "Apply full-aperture linear wavefront predictor to cube",
        .description_long = "Applies the trained `outPF` autoregressive linear predictor to "
                            "input phase cube .wfin, writing predicted frames to .wfoutp and "
                            "prediction residuals to .wfoutf."
    },
    .srcfile  = __FILE__,
    .farg     = farg_lpapply,
    .nbarg    = (int) (sizeof(farg_lpapply) / sizeof(farg_lpapply[0])),
    .cli_func = cli_atmturb_linpred_apply,
    .rt       = &atmturb_linpred_apply_rt,
};

/* 5. contract: spatially decimate wavefront FITS cube(s) by 2^factor */
static char    *contract_inp;
static char    *contract_outp;
static int64_t *contract_factor;
static char    *contract_ina;
static char    *contract_outa;

static CLICMDARGDEF farg_contract[] = {
    { CLIARG_STR, ".inp", "Input phase FITS cube filename", "outarraypha.fits",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &contract_inp, NULL },
    { CLIARG_STR, ".outp", "Output phase FITS cube filename", "outarraypha_bin.fits",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &contract_outp, NULL },
    { CLIARG_INT64, ".factor", "Power-of-2 decimation exponent", "1",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT | FPFLAG_CLI_OPTIONAL),
      (void **) &contract_factor, NULL },
    { CLIARG_STR, ".ina", "Optional input amplitude FITS filename", "",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT | FPFLAG_CLI_OPTIONAL),
      (void **) &contract_ina, NULL },
    { CLIARG_STR, ".outa", "Optional output amplitude FITS filename", "",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT | FPFLAG_CLI_OPTIONAL),
      (void **) &contract_outa, NULL },
};

/**
 * cli_atmturb_contract - Spatially decimate phase or complex wavefront FITS cubes
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE on error.
 */
static errno_t cli_atmturb_contract(void)
{
    int rc;
    if (contract_ina[0] != '\0' && contract_outa[0] != '\0')
    {
        rc = contract_wavefront_cube(contract_ina, contract_inp, contract_outa, contract_outp,
                                     (int) *contract_factor);
    }
    else
    {
        rc = contract_wavefront_cube_phaseonly(contract_inp, contract_outp,
                                               (int) *contract_factor);
    }
    return (rc == 0) ? RETURN_SUCCESS : RETURN_FAILURE;
}

MILK_RUNNABLE_RT_DEFINE(atmturb_contract);
const MILK_RUNNABLE runnable_atmturb_contract = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = {
        .cmdkey           = "contract",
        .description      = "Spatially decimate wavefront FITS cube(s) by 2^factor",
        .description_long = "Bins 3D phase-only (or complex amplitude + phase when .ina/.outa "
                            "are provided) FITS wavefront cubes spatially by 2^(.factor) x "
                            "2^(.factor) per slice."
    },
    .srcfile  = __FILE__,
    .farg     = farg_contract,
    .nbarg    = (int) (sizeof(farg_contract) / sizeof(farg_contract[0])),
    .cli_func = cli_atmturb_contract,
    .rt       = &atmturb_contract_rt,
};
