#ifndef _IMAGEFORMATMODULE_H
#define _IMAGEFORMATMODULE_H


int init_image_format();


int IMAGE_FORMAT_im_to_ASCII(char *IDname, char *foutname);

int IMAGE_FORMAT_FITS_to_ASCII(char *finname, char *foutname);

int image_writeBMP_auto(char *IDnameR, char *IDnameG, char *IDnameB, char *outname);

int CR2toFITS(char *fnameCR2, char *fnameFITS);

long IMAGE_FORMAT_FITS_to_ushortintbin_lock( char *IDname, char *fname);

long IMAGE_FORMAT_FITS_to_floatbin_lock(  char *IDname, char *fname);

long IMAGE_FORMAT_read_binary32f(char *fname, long xsize, long ysize, char *IDname);

int loadCR2toFITSRGB(char *fnameCR2, char *fnameFITSr, char *fnameFITSg, char *fnameFITSb);

int image_format_extract_RGGBchan(char *ID_name, char *IDoutR_name, char *IDoutG1_name, char *IDoutG2_name, char *IDoutB_name);

#endif

