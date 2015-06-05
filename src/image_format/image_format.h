#ifndef _IMAGEFORMATMODULE_H
#define _IMAGEFORMATMODULE_H

int IMAGE_FORMAT_im_to_ASCII(char *IDname, char *foutname);

int IMAGE_FORMAT_FITS_to_ASCII(char *finname, char *foutname);

int CR2toFITS(char *fnameCR2, char *fnameFITS);

long IMAGE_FORMAT_FITS_to_ushortintbin_lock( char *IDname, char *fname);

long IMAGE_FORMAT_FITS_to_floatbin_lock(  char *IDname, char *fname);

long IMAGE_FORMAT_read_binary32f(char *fname, long xsize, long ysize, char *IDname);

#endif

