// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmmod_standard_model.c
 * @brief   Standard atmospheric profile generation, saving, and loading
 */

#include "AtmosphereModel.h"
#include "atmmod_types.h"
#include "nrlmsise-00.h"

/**
 * AtmosphereModel_H2O_Saturation - Water vapor saturation pressure
 * @T: Temperature in Kelvin.
 *
 * Return: Saturation vapor pressure in Pascals using IAPWS-95 formulation.
 */
double AtmosphereModel_H2O_Saturation(double T)
{
    const double Tc = 647.096;      // critical temperature [K]
    const double Pc = 22064000.0;   // critical pressure [Pa]
    const double C1 = -7.85951783;
    const double C2 = 1.84408259;
    const double C3 = -11.7866497;
    const double C4 = 22.6807411;
    const double C5 = -15.9618719;
    const double C6 = 1.80122502;
    const double Tn = 273.16;       // triple point [K]
    const double a0 = -13.928169;
    const double a1 = 34.707823;
    const double Pn = 611.657;

    if (T > 273.0)
    {
        double ups = 1.0 - T / Tc;
        return Pc * exp(Tc / T * (C1 * ups + C2 * pow(ups, 1.5) + C3 * pow(ups, 3.0) +
                                  C4 * pow(ups, 3.5) + C5 * pow(ups, 4.0) + C6 * pow(ups, 7.5)));
    }

    double ups = T / Tn;
    return Pn * exp(a0 * (1.0 - pow(ups, -1.5)) + a1 * (1.0 - pow(ups, -1.25)));
}

/**
 * AtmosphereModel_save_stdAtmModel - Save standard atmosphere profile to file
 * @fname: Output filename.
 *
 * Return: 0 on success, -1 on failure.
 */
int AtmosphereModel_save_stdAtmModel(char *fname)
{
    FILE *fp = fopen(fname, "w");
    if (fp == NULL)
    {
        return -1;
    }

    fprintf(fp, "#  1:alt[m]  2:denstot[part/cm3] 3:N2  4:O2  5:Ar  6:H2O  7:CO2  8:Ne  "
                "9:He  10:CH4  11:Kr  12:H2  13:N  14:O  15:H    16:density[g/cm3]  "
                "17:temperature[K] 18:pressure[stdatm]  19:RH\n");

    for (long i = 0; i < ATMMOD_NB_BINS; i++)
    {
        fprintf(fp, "%6.0f   %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  "
                    "%15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g  %15.8g   "
                    "%15.8g    %8.3lf   %12.10lf  %7.5f\n",
                10.0 * i, denstot[i], densN2[i], densO2[i], densAr[i], densH2O[i],
                densCO2[i], densNe[i], densHe[i], densCH4[i], densKr[i], densH2[i],
                densO3[i], densN[i], densO[i], densH[i], density[i], temperature[i],
                pressure[i], RH[i]);
    }
    fclose(fp);
    return 0;
}

/**
 * AtmosphereModel_load_stdAtmModel - Load atmosphere profile from file
 * @fname: Input filename.
 *
 * Return: 0 on success, -1 on failure.
 */
int AtmosphereModel_load_stdAtmModel(char *fname)
{
    printf("Loading atmosphere model \"%s\"\n", fname);
    fflush(stdout);

    FILE *fp = fopen(fname, "r");
    if (fp == NULL)
    {
        printf("ERROR: file \"%s\" missing\n", fname);
        return -1;
    }

    char line[500];
    if (fgets(line, sizeof(line), fp) == NULL)
    {
        fclose(fp);
        return -1;
    }

    for (long i = 0; i < ATMMOD_NB_BINS; i++)
    {
        float v[20];
        if (fscanf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                   &v[0], &v[1], &v[2], &v[3], &v[4], &v[5], &v[6], &v[7], &v[8], &v[9],
                   &v[10], &v[11], &v[12], &v[13], &v[14], &v[15], &v[16], &v[17], &v[18], &v[19]) == 20)
        {
            denstot[i] = v[1];
            densN2[i] = v[2];
            densO2[i] = v[3];
            densAr[i] = v[4];
            densH2O[i] = v[5];
            densCO2[i] = v[6];
            densNe[i] = v[7];
            densHe[i] = v[8];
            densCH4[i] = v[9];
            densKr[i] = v[10];
            densH2[i] = v[11];
            densO3[i] = v[12];
            densN[i] = v[13];
            densO[i] = v[14];
            densH[i] = v[15];
            temperature[i] = v[17];
            pressure[i] = v[18];
            RH[i] = v[19];
        }
    }
    fclose(fp);
    return 0;
}

/**
 * atmmod_edlen_refractivity - Compute fallback air refractivity via Edlen formula
 * @alt: Altitude above sea level in meters.
 * @lambda: Optical wavelength in meters.
 *
 * Return: Refractivity (n - 1) scaled by exponential scale height.
 */
static float atmmod_edlen_refractivity(float alt, float lambda)
{
    double lambda_um = (double)lambda * 1e6;
    if (lambda_um <= 0.01)
    {
        lambda_um = 0.55;
    }
    double s2 = 1.0 / (lambda_um * lambda_um);
    double n_minus_1_1e6 = 287.6155 + 1.62887 * s2 + 0.01360 * s2 * s2;
    double scale = exp(-(double)alt / 8400.0);
    return (float)(n_minus_1_1e6 * 1e-6 * scale);
}

/**
 * AtmosphereModel_stdAtmModel_N - Refractive index minus 1 at specified altitude
 * @alt: Altitude above sea level in meters.
 * @lambda: Optical wavelength in meters.
 * @mode: Verbosity flag (1 for debug print, 0 for silent).
 *
 * Return: Refractivity (n - 1).
 */
float AtmosphereModel_stdAtmModel_N(float alt, float lambda, int mode)
{
    if (densN2 == NULL)
    {
        return atmmod_edlen_refractivity(alt, lambda);
    }

    long i = (long)(alt / 10.0);
    if (i > 9998)
    {
        i = 9998;
    }
    float ifrac = 1.0f * alt / 10.0f - i;
    if (ifrac < 0.0f)
    {
        ifrac = 0.0f;
    }
    if (ifrac > 1.0f)
    {
        ifrac = 1.0f;
    }

    double w0 = 1.0 - ifrac;
    double w1 = ifrac;

    double dN2 = w0 * densN2[i] + w1 * densN2[i + 1];
    double dO2 = w0 * densO2[i] + w1 * densO2[i + 1];
    double dAr = w0 * densAr[i] + w1 * densAr[i + 1];
    double dH2O = w0 * densH2O[i] + w1 * densH2O[i + 1];
    double dCO2 = w0 * densCO2[i] + w1 * densCO2[i + 1];
    double dNe = w0 * densNe[i] + w1 * densNe[i + 1];
    double dHe = w0 * densHe[i] + w1 * densHe[i + 1];
    double dCH4 = w0 * densCH4[i] + w1 * densCH4[i + 1];
    double dKr = w0 * densKr[i] + w1 * densKr[i + 1];
    double dH2 = w0 * densH2[i] + w1 * densH2[i + 1];
    double dO3 = w0 * densO3[i] + w1 * densO3[i + 1];
    double dN = w0 * densN[i] + w1 * densN[i + 1];
    double dO = w0 * densO[i] + w1 * densO[i + 1];
    double dH = w0 * densH[i] + w1 * densH[i + 1];

    v_ABSCOEFF = 0.0;
    float val = (float)AirMixture_N(lambda, dN2, dO2, dAr, dH2O, dCO2, dNe, dHe,
                                   dCH4, dKr, dH2, dO3, dN, dO, dH);

    if (mode == 1)
    {
        printf("\nalt = %f m\ni = %ld\n\n", alt, i);
        printf("     N2  =  %5.3f\n     O2  =  %5.3f\n\n",
               dN2 * 1.0e6 / ATMMOD_LOSCHMIDT, dO2 * 1.0e6 / ATMMOD_LOSCHMIDT);
    }

    return val;
}

/**
 * atmmod_setup_nrlmsise_inputs - Initialize NRLMSISE-00 inputs and flags
 * @input: Array of 10000 NRLMSISE input structs.
 * @flags: NRLMSISE control flags struct.
 */
static void atmmod_setup_nrlmsise_inputs(struct nrlmsise_input *input,
                                        struct nrlmsise_flags *flags)
{
    flags->switches[0] = 0;
    for (int i = 1; i < 24; i++)
    {
        flags->switches[i] = 1;
    }

    int sec = (int)(3600.0 * (TimeLocalSolarTime - SiteLong / 15.0));
    if (sec < 0)
    {
        sec += 3600 * 24;
    }

    for (int i = 0; i < ATMMOD_NB_BINS; i++)
    {
        input[i].doy = TimeDayOfYear;
        input[i].year = 0;
        input[i].alt = 0.01 * i;
        input[i].g_lat = SiteLat;
        input[i].g_long = SiteLong;
        input[i].lst = TimeLocalSolarTime;
        input[i].sec = sec;
        input[i].f107A = 150;
        input[i].f107 = 150;
        input[i].ap = 4;
    }
}

/**
 * atmmod_match_site_conditions - Align model altitude and pressure to site
 * @input: Pointer to first NRLMSISE input.
 * @flags: NRLMSISE control flags.
 * @output: Output struct for site query.
 * @deltah_out: Computed altitude offset in meters.
 * @press_coeff_out: Computed pressure scaling factor.
 */
static void atmmod_match_site_conditions(struct nrlmsise_input *input,
                                         struct nrlmsise_flags *flags,
                                         struct nrlmsise_output *output,
                                         double *deltah_out, double *press_coeff_out)
{
    double h = SiteAlt;
    input->alt = 0.001 * h;
    gtd7(input, flags, output);
    double dens_ne = output->d[2] * 2.328e-5;
    double d_tot = output->d[0] + output->d[1] + output->d[2] + output->d[3] +
                   output->d[4] + output->d[6] + output->d[7] + output->d[8] + dens_ne;
    double press = d_tot * 1.0e6 / ATMMOD_LOSCHMIDT * (output->t[1] / 273.15);

    if (SiteTPauto == 0)
    {
        for (int k = 0; k < 10; k++)
        {
            input->alt = 0.001 * h;
            gtd7(input, flags, output);
            dens_ne = output->d[2] * 2.328e-5;
            d_tot = output->d[0] + output->d[1] + output->d[2] + output->d[3] +
                    output->d[4] + output->d[6] + output->d[7] + output->d[8] + dens_ne;
            press = d_tot * 1.0e6 / ATMMOD_LOSCHMIDT * (output->t[1] / 273.15);
            h += 1000.0 * (output->t[1] - SiteTemp) / 6.5;
        }
        *deltah_out = h - SiteAlt;
        *press_coeff_out = SitePress / press;
    }
    else
    {
        *deltah_out = 0.0;
        *press_coeff_out = 1.0;
    }
}

/**
 * atmmod_compute_raw_profiles - Query NRLMSISE for all altitude bins
 * @fp: File pointer to dump raw profile.
 * @input: Array of NRLMSISE inputs.
 * @flags: NRLMSISE control flags.
 * @output: Array of NRLMSISE outputs.
 * @deltah: Altitude offset in meters.
 * @press_coeff: Pressure scaling factor.
 * @TotPart0: Output array for total particles per bin.
 * @TotPart_out: Cumulative particles above site.
 * @TotPart1_out: Cumulative H2O background particles.
 * @TotPart2_out: Cumulative H2O exponential mixing particles.
 */
static void atmmod_compute_raw_profiles(FILE *fp, struct nrlmsise_input *input,
                                        struct nrlmsise_flags *flags,
                                        struct nrlmsise_output *output,
                                        double deltah, double press_coeff,
                                        double *TotPart0, double *TotPart_out,
                                        double *TotPart1_out, double *TotPart2_out)
{
    double tot_part = 0.0, tot_part1 = 0.0, tot_part2 = 0.0;

    for (long i = 0; i < ATMMOD_NB_BINS; i++)
    {
        double h = 10.0 * i + deltah;
        input[i].alt = 0.001 * h;
        gtd7(&input[i], flags, &output[i]);
        fprintf(fp, "%6.0f %12f %12f %12f %12f %12f %12f %12f %.8g %5f\n",
                input[i].alt * 1000.0, output[i].d[2], output[i].d[3], output[i].d[4],
                output[i].d[6], output[i].d[0], output[i].d[1] + output[i].d[8],
                output[i].d[7], output[i].d[5], output[i].t[1]);

        for (int k = 0; k < 9; k++)
        {
            output[i].d[k] *= press_coeff;
        }

        TotPart0[i] = output[i].d[0] + output[i].d[1] + output[i].d[2] + output[i].d[3] +
                      output[i].d[4] + output[i].d[6] + output[i].d[7] + output[i].d[8];

        densNe[i] = output[i].d[2] * 2.33e-5;
        TotPart0[i] += densNe[i];

        densCH4[i] = (h < 45000.0) ? TotPart0[i] * 2e-6 * (1.0 - h / 45000.0) : 0.0f;
        TotPart0[i] += densCH4[i];

        densKr[i] = output[i].d[2] * 1.46e-6;
        TotPart0[i] += densKr[i];

        densH2[i] = output[i].d[2] * 7.04e-7;
        TotPart0[i] += densH2[i];

        densO3[i] = (float)(300.0 * 2.69e20 / (4250.0 * sqrt(2.0 * M_PI)) *
                            exp(-0.5 * pow((h - 25000.0) / 4250.0, 2.0)) * 1.0e-6);

        densCO2[i] = (h < 70000.0) ? (float)(CO2_ppm * 1e-6 * TotPart0[i])
                                   : (float)((CO2_ppm - 0.007 * (h - 70000.0)) * 1e-6 * TotPart0[i]);
        output[i].d[3] -= densCO2[i];

        if (i == 0)
        {
            dens0 = (float)TotPart0[i];
        }
        TotPart0[i] *= 1000.0;
        temperature[i] = (float)output[i].t[1];

        if (h > SiteAlt)
        {
            tot_part += TotPart0[i];
            tot_part1 += 2.5e-6 * TotPart0[i] * exp(-3.0 * pow(h / 100000.0, 4.0));
            tot_part2 += exp(-h / SitePWSH) * TotPart0[i] * exp(-3.0 * pow(h / 100000.0, 4.0));
        }
        densH2O[i] = (float)(2.5e-6 * TotPart0[i] / 1000.0);
    }

    *TotPart_out = tot_part;
    *TotPart1_out = tot_part1;
    *TotPart2_out = tot_part2;
}

/**
 * atmmod_adjust_water_vapor - Scale water vapor profile based on site humidity
 * @TotPart0: Total particles per bin array.
 * @TotPart1: Background H2O particles.
 * @TotPart2: Exponential H2O particles.
 */
static void atmmod_adjust_water_vapor(const double *TotPart0, double TotPart1, double TotPart2)
{
    if (SiteH2OMethod == 1)
    {
        double X = (SiteTPW * 0.1) * 6.0221413e23 / 18.0;
        alpha1H2O = (float)((X - TotPart1) / TotPart2);
        if (alpha1H2O < 0.0f)
        {
            printf("ERROR: total precipitable water value is too low\n");
            exit(0);
        }
        for (long i = 0; i < ATMMOD_NB_BINS; i++)
        {
            double h = 10.0 * i;
            densH2O[i] += (float)(alpha1H2O * exp(-h / SitePWSH) * (TotPart0[i] / 1000.0) *
                                  exp(-3.0 * pow(h / 100000.0, 4.0)));
        }
    }
    else if (SiteH2OMethod == 2 || SiteH2OMethod == 3)
    {
        long i0 = (long)(SiteAlt / 10.0);
        float ifrac = SiteAlt / 10.0f - i0;
        double H2OSat = AtmosphereModel_H2O_Saturation((1.0 - ifrac) * temperature[i0] +
                                                       ifrac * temperature[i0 + 1]);
        double densH2O_site = (SiteRH / 100.0) * ATMMOD_LOSCHMIDT * 1e-6 * H2OSat / 101325.0;

        double pwsh = SitePWSH;
        if (SiteH2OMethod == 3)
        {
            double tpw0 = SiteTPW - (TotPart1 * 10.0 * 18.0 / 6.022e23);
            pwsh = 1000.0;
            double deltapwsh = 100.0;
            int dir = 1;
            for (long iter = 0; iter < 15; iter++)
            {
                double part2 = 0.0;
                for (long i = 0; i < ATMMOD_NB_BINS; i++)
                {
                    double h = 10.0 * i;
                    if (h > SiteAlt)
                    {
                        part2 += 1000.0 * densH2O_site * exp(-(h - SiteAlt) / pwsh) *
                                 exp(-3.0 * pow(h / 100000.0, 4.0));
                    }
                }
                double tpw = part2 * 10.0 * 18.0 / 6.022e23;
                int odir = dir;
                if (tpw > tpw0)
                {
                    pwsh -= deltapwsh;
                    dir = -1;
                }
                else
                {
                    pwsh += deltapwsh;
                    dir = 1;
                }
                if (odir != dir)
                {
                    deltapwsh *= 0.5;
                }
            }
        }

        for (long i = 0; i < ATMMOD_NB_BINS; i++)
        {
            double h = 10.0 * i;
            densH2O[i] += (float)(densH2O_site * exp(-(h - SiteAlt) / pwsh) *
                                  (TotPart0[i] / TotPart0[i0]) *
                                  exp(-3.0 * pow(h / 100000.0, 4.0)));
        }
    }
}

/**
 * atmmod_normalize_and_finalize - Normalize species fractions and compute thermodynamic fields
 * @output: Array of NRLMSISE outputs.
 */
static void atmmod_normalize_and_finalize(struct nrlmsise_output *output)
{
    for (long i = 0; i < ATMMOD_NB_BINS; i++)
    {
        denstot[i] = (float)(output[i].d[0] + output[i].d[1] + output[i].d[2] + output[i].d[3] +
                             output[i].d[4] + output[i].d[6] + output[i].d[7] + output[i].d[8] +
                             densCO2[i] + densNe[i]);
        double coeff = (denstot[i] + densH2O[i] + densO3[i]) / denstot[i];

        output[i].d[0] /= coeff;
        densHe[i] = (float)output[i].d[0];
        output[i].d[1] /= coeff;
        densO[i] = (float)output[i].d[1];
        output[i].d[2] /= coeff;
        densN2[i] = (float)output[i].d[2];
        output[i].d[3] /= coeff;
        densO2[i] = (float)output[i].d[3];
        output[i].d[4] /= coeff;
        densAr[i] = (float)output[i].d[4];
        output[i].d[6] /= coeff;
        densH[i] = (float)output[i].d[6];
        output[i].d[7] /= coeff;
        densN[i] = (float)output[i].d[7];
        output[i].d[8] /= coeff;
        densO[i] += (float)output[i].d[8];

        densH2O[i] /= (float)coeff;
        densCO2[i] /= (float)coeff;
        densNe[i] /= (float)coeff;

        densKr[i] = (float)(densN2[i] / 78.084 * 1.14e-6);
        densH2[i] = (float)(densN2[i] / 78.084 * 5.5e-7);

        density[i] = (float)output[i].d[5];
        temperature[i] = (float)output[i].t[1];
        pressure[i] = (float)(denstot[i] * 1.0e6 / ATMMOD_LOSCHMIDT * (output[i].t[1] / 273.15));
        RH[i] = (float)((densH2O[i] * 1e6 / ATMMOD_LOSCHMIDT) * 101325.0 /
                        AtmosphereModel_H2O_Saturation(temperature[i]));
    }
}

/**
 * AtmosphereModel_build_stdAtmModel - Generate profile from NRLMSISE-00
 * @fname: Output filename to save generated model.
 *
 * Return: 0 on success.
 */
int AtmosphereModel_build_stdAtmModel(char *fname)
{
    struct nrlmsise_output *output = malloc(sizeof(struct nrlmsise_output) * ATMMOD_NB_BINS);
    struct nrlmsise_input *input = malloc(sizeof(struct nrlmsise_input) * ATMMOD_NB_BINS);
    struct nrlmsise_flags flags;
    double *TotPart0 = malloc(sizeof(double) * ATMMOD_NB_BINS);

    atmmod_setup_nrlmsise_inputs(input, &flags);

    double deltah = 0.0;
    double press_coeff = 1.0;
    atmmod_match_site_conditions(&input[0], &flags, &output[0], &deltah, &press_coeff);

    FILE *fp = fopen(fname, "w");
    fprintf(fp, "#  alt[m]  N2  O2  Ar  H  He  O   N  density   Temperature\n");

    double tot_part = 0.0, tot_part1 = 0.0, tot_part2 = 0.0;
    atmmod_compute_raw_profiles(fp, input, &flags, output, deltah, press_coeff,
                               TotPart0, &tot_part, &tot_part1, &tot_part2);
    fclose(fp);

    atmmod_adjust_water_vapor(TotPart0, tot_part1, tot_part2);
    atmmod_normalize_and_finalize(output);

    free(TotPart0);
    free(input);
    free(output);

    AtmosphereModel_save_stdAtmModel(fname);
    return 0;
}
