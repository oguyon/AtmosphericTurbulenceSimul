// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_psf_sim.c
 * @brief   End-to-end adaptive optics closed-loop simulation and PSF synthesis
 */

#include "AtmosphericTurbulence.h"
#include "atmturb_types.h"

typedef struct
{
    long wf_size1;
    imageID id_telpup;
    imageID id_atm_opd;
    imageID id_atm_amp;
    imageID id_wfs_opd;
    imageID id_wfs_amp;
    imageID id_dm_opd;
    imageID id_dm_opd_tmp;
    imageID id_wfs_mes_opd;
    imageID id_wfs_mes_opd_prev;
    imageID id_wfs_mes_opd_der;
    imageID id_wfs_mes_opd_int;
    imageID id_sci_opd;
    imageID id_sci_amp;
    imageID id_psf_cumul;
} atmturb_ao_sim_state_t;

/**
 * atmturb_psf_sim_init_images - Allocate simulation buffers and telescope pupil
 * @st: Pointer to AO simulation state structure.
 * @wf_size1: Sub-binned pupil linear dimension.
 * @tel_diam: Telescope primary diameter in meters.
 */
static void atmturb_psf_sim_init_images(atmturb_ao_sim_state_t *st, long wf_size1, double tel_diam)
{
    st->wf_size1 = wf_size1;
    st->id_telpup = make_disk("TelPup", wf_size1, wf_size1, wf_size1 / 2, wf_size1 / 2,
                              tel_diam * 0.5 / CONF_PUPIL_SCALE / 2.0);

    st->id_atm_opd = create_2Dimage_ID("atmopd", wf_size1, wf_size1);
    st->id_atm_amp = create_2Dimage_ID("atmamp", wf_size1, wf_size1);
    st->id_wfs_opd = create_2Dimage_ID("wfsopd", wf_size1, wf_size1);
    st->id_wfs_amp = create_2Dimage_ID("wfsamp", wf_size1, wf_size1);

    st->id_dm_opd_tmp = create_2Dimage_ID("dmopdtmp", wf_size1, wf_size1);
    st->id_dm_opd = create_2Dimage_ID("dmopd", wf_size1, wf_size1);

    st->id_wfs_mes_opd = create_2Dimage_ID("wfsmesopd", wf_size1, wf_size1);
    st->id_wfs_mes_opd_prev = create_2Dimage_ID("wfsmesopdprev", wf_size1, wf_size1);
    st->id_wfs_mes_opd_der = create_2Dimage_ID("wfsmesopdder", wf_size1, wf_size1);
    st->id_wfs_mes_opd_int = create_2Dimage_ID("wfsmesopdint", wf_size1, wf_size1);

    st->id_sci_opd = create_2Dimage_ID("sciopd", wf_size1, wf_size1);
    st->id_sci_amp = create_2Dimage_ID("sciamp", wf_size1, wf_size1);
    st->id_psf_cumul = create_2Dimage_ID("PSFcumul", 512, 512);

    for (long i = 0; i < wf_size1 * wf_size1; i++)
    {
        dcimg[st->id_wfs_opd].array.F[i] = 0.0f;
        dcimg[st->id_dm_opd].array.F[i] = 0.0f;
        dcimg[st->id_wfs_mes_opd].array.F[i] = 0.0f;
        dcimg[st->id_wfs_mes_opd_prev].array.F[i] = 0.0f;
        dcimg[st->id_wfs_mes_opd_der].array.F[i] = 0.0f;
        dcimg[st->id_wfs_mes_opd_int].array.F[i] = 0.0f;
    }
}

/**
 * atmturb_psf_sim_update_pid - Perform PID control step on measured wavefront OPD
 * @st: Simulation state pointer.
 * @Kp: Proportional gain.
 * @Ki: Integral gain.
 * @Kd: Derivative gain.
 * @dt: Time step duration in seconds.
 */
static void atmturb_psf_sim_update_pid(atmturb_ao_sim_state_t *st, double Kp, double Ki,
                                      double Kd, double dt)
{
    long ntot = st->wf_size1 * st->wf_size1;
    for (long i = 0; i < ntot; i++)
    {
        float err = dcimg[st->id_wfs_mes_opd].array.F[i];
        float prev = dcimg[st->id_wfs_mes_opd_prev].array.F[i];
        float der = (dt > 0.0) ? (err - prev) / (float)dt : 0.0f;

        dcimg[st->id_wfs_mes_opd_der].array.F[i] = der;
        dcimg[st->id_wfs_mes_opd_int].array.F[i] += err * (float)dt;

        float corr = (float)(Kp * err + Ki * dcimg[st->id_wfs_mes_opd_int].array.F[i] + Kd * der);
        dcimg[st->id_dm_opd_tmp].array.F[i] = corr;
        dcimg[st->id_wfs_mes_opd_prev].array.F[i] = err;
    }
}

/**
 * atmturb_psf_sim_accumulate_psf - Form instantaneous PSF and add to cumulative detector
 * @st: Simulation state pointer.
 * @scilambda: Science observing wavelength in meters.
 */
static void atmturb_psf_sim_accumulate_psf(atmturb_ao_sim_state_t *st, double scilambda)
{
    long ntot = st->wf_size1 * st->wf_size1;
    imageID id_arr = create_2DCimage_ID("tmp_ao_arr", st->wf_size1, st->wf_size1);

    for (long i = 0; i < ntot; i++)
    {
        float opd = dcimg[st->id_sci_opd].array.F[i];
        float amp = dcimg[st->id_sci_amp].array.F[i] * dcimg[st->id_telpup].array.F[i];
        float pha = (float)(2.0 * M_PI * opd / scilambda);
        dcimg[id_arr].array.CF[i].re = amp * cosf(pha);
        dcimg[id_arr].array.CF[i].im = amp * sinf(pha);
    }

    permut("tmp_ao_arr");
    do2dfft("tmp_ao_arr", "tmp_ao_fft");
    delete_image_ID("tmp_ao_arr");

    imageID id_fft = image_ID("tmp_ao_fft");
    for (long i = 0; i < 512 * 512 && i < ntot; i++)
    {
        float re = dcimg[id_fft].array.CF[i].re;
        float im = dcimg[id_fft].array.CF[i].im;
        dcimg[st->id_psf_cumul].array.F[i] += re * re + im * im;
    }
    delete_image_ID("tmp_ao_fft");
}

/**
 * AtmosphericTurbulence_makePSF - Run closed-loop AO simulation and generate PSF
 * @Kp: Proportional feedback gain.
 * @Ki: Integral feedback gain.
 * @Kd: Derivative feedback gain.
 * @Kdgain: Dynamic gain multiplier.
 *
 * Return: Peak intensity or Strehl estimate.
 */
double AtmosphericTurbulence_makePSF(double Kp, double Ki, double Kd, double Kdgain)
{
    (void)Kdgain;
    AtmosphericTurbulence_ReadConf();

    double TelDiam = 30.0;
    double etime = 1.0;
    double dtime = 0.001;
    double scilambda = 1.65e-6;

    long wf_size1 = CONF_WFsize / 2;
    if (wf_size1 < 64) wf_size1 = 64;

    atmturb_ao_sim_state_t st;
    atmturb_psf_sim_init_images(&st, wf_size1, TelDiam);

    double rtime = 0.0;
    while (rtime < etime)
    {
        for (long i = 0; i < wf_size1 * wf_size1; i++)
        {
            dcimg[st.id_wfs_mes_opd].array.F[i] =
                dcimg[st.id_atm_opd].array.F[i] - dcimg[st.id_dm_opd].array.F[i];
            dcimg[st.id_sci_opd].array.F[i] =
                dcimg[st.id_atm_opd].array.F[i] - dcimg[st.id_dm_opd].array.F[i];
            dcimg[st.id_sci_amp].array.F[i] = 1.0f;
        }

        atmturb_psf_sim_update_pid(&st, Kp, Ki, Kd, dtime);

        for (long i = 0; i < wf_size1 * wf_size1; i++)
        {
            dcimg[st.id_dm_opd].array.F[i] += dcimg[st.id_dm_opd_tmp].array.F[i];
        }

        atmturb_psf_sim_accumulate_psf(&st, scilambda);
        rtime += dtime;
    }

    double peak = 0.0;
    for (long i = 0; i < 512 * 512; i++)
    {
        if (dcimg[st.id_psf_cumul].array.F[i] > peak)
        {
            peak = dcimg[st.id_psf_cumul].array.F[i];
        }
    }

    save_fl_fits("PSFcumul", "!PSFcumul.fits");

    delete_image_ID("atmopd");
    delete_image_ID("atmamp");
    delete_image_ID("wfsopd");
    delete_image_ID("wfsamp");
    delete_image_ID("dmopdtmp");
    delete_image_ID("dmopd");
    delete_image_ID("wfsmesopd");
    delete_image_ID("wfsmesopdprev");
    delete_image_ID("wfsmesopdder");
    delete_image_ID("wfsmesopdint");
    delete_image_ID("sciopd");
    delete_image_ID("sciamp");
    delete_image_ID("TelPup");
    delete_image_ID("PSFcumul");

    return peak;
}
