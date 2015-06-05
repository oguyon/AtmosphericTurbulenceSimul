#ifndef _AtmosphericTurbulence_H
#define _AtmosphericTurbulence_H


int init_AtmosphericTurbulence();

int AtmosphericTurbulence_change_configuration_file(char *fname);


int make_master_turbulence_screen(char *ID_name1, char *ID_name2, long size, float outerscale, float innerscale);

int make_master_turbulence_screen_pow(char *ID_name1, char *ID_name2, long size, float power);

//int unwrap_phase_screen(char *ID_name);

int contract_wavefront_series(char *in_prefix, char *out_prefix, long NB_files);

int contract_wavefront_cube(char *ina_file, char *inp_file, char *outa_file, char *outp_file, int factor);

int contract_wavefront_cube_phaseonly(char *inp_file, char *outp_file, int factor);

int make_AtmosphericTurbulence_wavefront_series();

int measure_wavefront_series(float factor);

int measure_wavefront_series_expoframes(float etime, char *outfile);

int frame_select_PSF(char *logfile, long NBfiles, float frac);

int AtmosphericTurbulence_WFprocess();

int AtmosphericTurbulence_makeHV_CN2prof(double wspeed, double r0, double sitealt, char *outfile);

#endif
