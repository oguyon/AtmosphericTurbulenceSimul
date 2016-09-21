#ifndef _AtmosphericTurbulence_H
#define _AtmosphericTurbulence_H


int init_AtmosphericTurbulence();

int AtmosphericTurbulence_change_configuration_file(char *fname);

long make_AtmosphericTurbulence_vonKarmanWind(long vKsize, float pixscale, float sigmawind, float Lwind, long size, char *IDout_name);

int make_master_turbulence_screen(char *ID_name1, char *ID_name2, long size, float outerscale, float innerscale);

int make_master_turbulence_screen_pow(char *ID_name1, char *ID_name2, long size, float power);

//int unwrap_phase_screen(char *ID_name);

int contract_wavefront_series(char *in_prefix, char *out_prefix, long NB_files);

int contract_wavefront_cube(char *ina_file, char *inp_file, char *outa_file, char *outp_file, int factor);

int contract_wavefront_cube_phaseonly(char *inp_file, char *outp_file, int factor);

int make_AtmosphericTurbulence_wavefront_series(float slambdaum);

int measure_wavefront_series(float factor);
 
int AtmosphericTurbulence_mkTestTTseq(double dt, long NBpts, long NBblocks, double measnoise, int ACCnmode, double ACCnoise, int MODE);

int AtmosphericTurbulence_Build_LinPredictor_Full(char *WFin_name, char *WFmask_name, int PForder, float PFlag, double SVDeps, double RegLambda);
int AtmosphericTurbulence_Apply_LinPredictor_Full(int MODE, char *WFin_name, char *WFmask_name, int PForder, float PFlag, char *WFoutp_name, char *WFoutf_name);
long AtmosphericTurbulence_LinPredictor_filt_2DKernelExtract(char *IDfilt_name, char *IDmask_name, long krad, char *IDkern_name);
long AtmosphericTurbulence_LinPredictor_filt_Expand(char *IDfilt_name, char *IDmask_name);

int AtmosphericTurbulence_Build_LinPredictor(long NB_WFstep, double WFphaNoise, long WFPlag, long WFP_NBstep, long WFP_xyrad, long WFPiipix, long WFPjjpix, float slambdaum);
long AtmosphericTurbulence_psfCubeContrast(char *IDwfc_name, char *IDmask_name, char *IDpsfc_name);
int AtmosphericTurbulence_Test_LinPredictor(long NB_WFstep, double WFphaNoise, char *IDWFPfilt_name, long WFPlag, long WFPiipix, long WFPjjpix, float slambdaum);

int measure_wavefront_series_expoframes(float etime, char *outfile);

int frame_select_PSF(char *logfile, long NBfiles, float frac);

int AtmosphericTurbulence_WFprocess();

int AtmosphericTurbulence_makeHV_CN2prof(double wspeed, double r0, double sitealt, long NBlayer, char *outfile);

#endif
