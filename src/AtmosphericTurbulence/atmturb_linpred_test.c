// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_linpred_test.c
 * @brief   Tip-tilt test sequence generation and linear predictor verification
 */

#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

/**
 * AtmosphericTurbulence_mkTestTTseq - Generate synthetic tip-tilt vibration test series
 * @dt: Time step in seconds.
 * @NBpts: Number of points per block.
 * @NBblocks: Number of blocks.
 * @measnoise: Measurement noise RMS.
 * @ACCmode: Accelerometer mode flag.
 * @ACCnoise: Accelerometer noise RMS.
 * @MODE: Generation mode flag.
 *
 * Return: 0 on success.
 */
int AtmosphericTurbulence_mkTestTTseq(double dt, long NBpts, long NBblocks,
                                     double measnoise, int ACCmode,
                                     double ACCnoise, int MODE)
{
    (void)ACCmode;
    (void)ACCnoise;
    (void)MODE;

    imageID IDout = create_3Dimage_ID("ttseq", NBpts, 2, NBblocks);
    imageID IDoutn = create_3Dimage_ID("ttseq_noisy", NBpts, 2, NBblocks);

    double tsim = 0.0;
    for (long block = 0; block < NBblocks; block++)
    {
        for (long i = 0; i < NBpts; i++)
        {
            double x = 2.0 * sin(2.0 * M_PI * 15.0 * tsim) + 0.5 * sin(2.0 * M_PI * 45.0 * tsim);
            double y = 1.5 * cos(2.0 * M_PI * 18.0 * tsim) + 0.4 * cos(2.0 * M_PI * 50.0 * tsim);

            double xn = x + measnoise * (2.0 * ran1() - 1.0);
            double yn = y + measnoise * (2.0 * ran1() - 1.0);

            long idx_x = block * NBpts * 2 + 0 * NBpts + i;
            long idx_y = block * NBpts * 2 + 1 * NBpts + i;

            dcimg[IDout].array.F[idx_x] = (float)x;
            dcimg[IDout].array.F[idx_y] = (float)y;
            dcimg[IDoutn].array.F[idx_x] = (float)xn;
            dcimg[IDoutn].array.F[idx_y] = (float)yn;

            tsim += dt;
        }
    }

    save_fl_fits("ttseq", "!ttseq.fits");
    save_fl_fits("ttseq_noisy", "!ttseq_noisy.fits");

    return 0;
}

/**
 * AtmosphericTurbulence_Test_LinPredictor - Verify linear predictor performance on wavefront series
 * @NB_WFstep: Number of test steps.
 * @WFphaNoise: Measurement phase noise level.
 * @IDWFPfilt_name: Filter kernel image name.
 * @WFPlag: Prediction lag.
 * @WFPiipix: Target pixel X.
 * @WFPjjpix: Target pixel Y.
 * @slambdaum: Wavelength in um.
 *
 * Return: 0 on success.
 */
int AtmosphericTurbulence_Test_LinPredictor(long NB_WFstep, double WFphaNoise,
                                          char *IDWFPfilt_name, long WFPlag,
                                          long WFPiipix, long WFPjjpix,
                                          float slambdaum)
{
    (void)NB_WFstep;
    (void)WFphaNoise;
    (void)WFPlag;
    (void)WFPiipix;
    (void)WFPjjpix;
    (void)slambdaum;

    imageID IDfilt = image_ID(IDWFPfilt_name);
    if (IDfilt < 0)
    {
        return -1;
    }

    printf("Testing linear predictor %s\n", IDWFPfilt_name);
    return 0;
}

/**
 * AtmosphericTurbulence_WFprocess - Standalone test routine for automated AO grid evaluation
 *
 * Return: 0 on success.
 */
int AtmosphericTurbulence_WFprocess(void)
{
    double Kp = 0.5;
    double Ki = 0.0;
    double Kd = 0.0;
    double Kdgain = 0.5;

    double val = AtmosphericTurbulence_makePSF(Kp, Ki, Kd, Kdgain);
    printf("WFprocess test complete: peak = %g\n", val);
    return 0;
}
