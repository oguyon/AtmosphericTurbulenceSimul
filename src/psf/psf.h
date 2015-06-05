#ifndef _PSF_H
#define _PSF_H



int init_psf();


long PSF_makeChromatPSF(char *amp_name, char *pha_name, float coeff1, float coeff2, long NBstep, float ApoCoeff, char *out_name);

int PSF_finddiskcent(char *ID_name, float rad, float *result);

int PSF_finddiskcent_alone(char *ID_name, float rad);

int PSF_measurePhotocenter(char *ID_name);

float measure_enc_NRJ(char *ID_name, float xcenter, float ycenter, float fraction);

int measure_enc_NRJ1(char *ID_name, float xcenter, float ycenter, char *filename);

float measure_FWHM(char *ID_name, float xcenter, float ycenter, float step, long nb_step);

int center_PSF(char *ID_name, double *xcenter, double *ycenter, long box_size);

int fast_center_PSF(char *ID_name, double *xcenter, double *ycenter, long box_size);

int center_PSF_alone(char *ID_name);

int center_star(char *ID_in_name, double *x_star, double *y_star);

float get_sigma(char *ID_name, float x, float y, char *options);

float get_sigma_alone(char *ID_name);

int extract_psf(char *ID_name, char *out_name, long size);

long extract_psf_photcent(char *ID_name, char *out_name, long size);

int psf_variance(char *ID_out_m, char *ID_out_v, char *options);

int combine_2psf(char *ID_name, char *ID_name1, char *ID_name2, float radius, float index);

float psf_measure_SR(char *ID_name, float factor, float r1, float r2);

long PSF_coaddbest(char *IDcin_name, char *IDout_name, float r_pix);

#endif
