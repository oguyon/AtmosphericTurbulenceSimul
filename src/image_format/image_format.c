#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <sys/file.h>

#include "CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "COREMOD_iofits/COREMOD_iofits.h"

#include "image_format/image_format.h"


//#include <modules/image_basic/image_basic.h>
#include "info/info.h"
#include "fft/fft.h"
//#include <modules/image_filter/image_filter.h>
#include "image_gen/image_gen.h"

#define SWAP(x,y)  tmp=(x);x=(y);y=tmp;

#define PI 3.14159265358979323846264338328

extern DATA data;

#define BMP_BIG_ENDIAN 0

static int CR2toFITS_NORM = 0; // 1 if FITS should be normalized to ISO = 1, exposure = 1 sec, and F/1.0
static float FLUXFACTOR = 1.0;

typedef struct {int rows; int cols; unsigned char* data;} sImage;

/* This pragma is necessary so that the data in the structures is aligned to 2-byte 
   boundaries.  Some different compilers have a different syntax for this line.  For
   example, if you're using cc on Solaris, the line should be #pragma pack(2).  
*/
#pragma pack(2)


/* Default data types.  Here, uint16 is an unsigned integer that has size 2 bytes (16 bits), 
   and uint32 is datatype that has size 4 bytes (32 bits).  You may need to change these 
   depending on your compiler. */
#define uint16 unsigned short
#define uint32 unsigned int

#define BI_RGB 0
#define BM 19778
#define BMP_FALSE 0
#define BMP_TRUE 1

typedef struct {
   uint16 bfType; 
   uint32 bfSize; 
   uint16 bfReserved1; 
   uint16 bfReserved2; 
   uint32 bfOffBits; 
} BITMAPFILEHEADER; 

typedef struct { 
   uint32 biSize;
   uint32 biWidth; 
   uint32 biHeight; 
   uint16 biPlanes; 
   uint16 biBitCount; 
   uint32 biCompression; 
   uint32 biSizeImage; 
   uint32 biXPelsPerMeter; 
   uint32 biYPelsPerMeter; 
   uint32 biClrUsed; 
   uint32 biClrImportant; 
} BITMAPINFOHEADER; 


typedef struct {
   unsigned char rgbBlue;
   unsigned char rgbGreen;
   unsigned char rgbRed;
   unsigned char rgbReserved;
} RGBQUAD;



// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//



//int IMAGE_FORMAT_2Dim_to_ASCII(char *IDname, char *fname)

int IMAGE_FORMAT_im_to_ASCII_cli()
{
  if(CLI_checkarg(1,4)+CLI_checkarg(2,3)==0)
    {
      IMAGE_FORMAT_im_to_ASCII(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string);
      return 0;
    }
  else
    return 1;
}

int CR2toFITS_cli()
{
  //  if(CLI_checkarg(1, 3)+CLI_checkarg(2, 3))
  CR2toFITS(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string);
  // else
  // return(0);
}


int IMAGE_FORMAT_FITS_to_ushortintbin_lock_cli()
{
  if(CLI_checkarg(1,4)+CLI_checkarg(2,3)==0)
    {
      IMAGE_FORMAT_FITS_to_ushortintbin_lock(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string);
      return 0;
    }
  else
    return 1;
}

int IMAGE_FORMAT_FITS_to_floatbin_lock_cli()
{
  if(CLI_checkarg(1,4)+CLI_checkarg(2,3)==0)
    {
      IMAGE_FORMAT_FITS_to_floatbin_lock(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string);
      return 0;
    }
  else
    return 1;
}


int IMAGE_FORMAT_read_binary32f_cli()
{
  if(CLI_checkarg(1,3)+CLI_checkarg(2,2)+CLI_checkarg(3,2)+CLI_checkarg(4,3)==0)
    {
      IMAGE_FORMAT_read_binary32f(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.numl, data.cmdargtoken[3].val.numl, data.cmdargtoken[4].val.string);
      return 0;
    }
  else
    return 1;
}


int IMAGE_FORMAT_extract_RGGBchan_cli()
{
  if(CLI_checkarg(1,4)+CLI_checkarg(2,3)+CLI_checkarg(3,3)+CLI_checkarg(4,3)+CLI_checkarg(5,3)==0)
    {
      image_format_extract_RGGBchan(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.string, data.cmdargtoken[4].val.string, data.cmdargtoken[5].val.string);
      return 0;
    }
  else
    return 1;
}

int IMAGE_FORMAT_loadCR2toFITSRGB_cli()
{
  if(CLI_checkarg(1,4)+CLI_checkarg(2,3)+CLI_checkarg(3,3)+CLI_checkarg(4,3)==0)
    {
      loadCR2toFITSRGB(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.string, data.cmdargtoken[4].val.string);
      return 0;
    }
  else
    return 1;
}

//int loadCR2toFITSRGB(char *fnameCR2, char *fnameFITSr, char *fnameFITSg, char *fnameFITSb)


int init_image_format()
{
  strcpy(data.module[data.NBmodule].name, __FILE__);
  strcpy(data.module[data.NBmodule].info, "conversion between image format, I/O");
  data.NBmodule++;
  
  strcpy(data.cmd[data.NBcmd].key,"im2ascii");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = IMAGE_FORMAT_im_to_ASCII_cli;
  strcpy(data.cmd[data.NBcmd].info,"convert image file to ASCII");
  strcpy(data.cmd[data.NBcmd].syntax,"<input image> <output ASCII file>");
  strcpy(data.cmd[data.NBcmd].example,"im2ascii im im.txt");
  strcpy(data.cmd[data.NBcmd].Ccall,"int IMAGE_FORMAT_im_to_ASCII(char *IDname, char *fname)");
  data.NBcmd++;
  
  strcpy(data.cmd[data.NBcmd].key,"cr2tofits");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = CR2toFITS_cli;
  strcpy(data.cmd[data.NBcmd].info,"convert cr2 file to fits");
  strcpy(data.cmd[data.NBcmd].syntax,"<input CR2 file> <output FITS file>");
  strcpy(data.cmd[data.NBcmd].example,"cr2tofits im01.CR2 im01.fits");
  strcpy(data.cmd[data.NBcmd].Ccall,"int CR2toFITS(char *fnameCR2, char *fnameFITS)");
  data.NBcmd++;
 
  strcpy(data.cmd[data.NBcmd].key,"writeushortintlock");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = IMAGE_FORMAT_FITS_to_ushortintbin_lock_cli;
  strcpy(data.cmd[data.NBcmd].info,"write unsigned short int with file locking");
  strcpy(data.cmd[data.NBcmd].syntax,"str1 is image, str2 is binary file");
  strcpy(data.cmd[data.NBcmd].example,"writeushortintlock im im.bin");
  strcpy(data.cmd[data.NBcmd].Ccall,"long IMAGE_FORMAT_FITS_to_ushortintbin_lock( char *IDname, char *fname)");
  data.NBcmd++;

  strcpy(data.cmd[data.NBcmd].key,"writefloatlock");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = IMAGE_FORMAT_FITS_to_floatbin_lock_cli;
  strcpy(data.cmd[data.NBcmd].info,"write float with file locking");
  strcpy(data.cmd[data.NBcmd].syntax,"str1 is image, str2 is binary file");
  strcpy(data.cmd[data.NBcmd].example,"writefloatlock im im.bin");
  strcpy(data.cmd[data.NBcmd].Ccall,"long IMAGE_FORMAT_FITS_to_floatbin_lock( char *IDname, char *fname)");
  data.NBcmd++;

  strcpy(data.cmd[data.NBcmd].key,"readb32fim");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = IMAGE_FORMAT_read_binary32f_cli;
  strcpy(data.cmd[data.NBcmd].info,"read 32-bit float RAW image");
  strcpy(data.cmd[data.NBcmd].syntax,"<bin file> <xsize> <ysize> <output image>");
  strcpy(data.cmd[data.NBcmd].example,"readb32fim im.bin xsize ysize im");
  strcpy(data.cmd[data.NBcmd].Ccall,"long IMAGE_FORMAT_read_binary32f(char *fname, long xsize, long ysize, char *IDname)");
  data.NBcmd++;

  strcpy(data.cmd[data.NBcmd].key,"extractRGGBchan");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = IMAGE_FORMAT_extract_RGGBchan_cli;
  strcpy(data.cmd[data.NBcmd].info,"extract RGGB channels from color image");
  strcpy(data.cmd[data.NBcmd].syntax,"<input image> <imR> <imG1> <imG2> <imB>");
  strcpy(data.cmd[data.NBcmd].example,"extractRGGBchan im imR imG1 imG2 imB");
  strcpy(data.cmd[data.NBcmd].Ccall,"int image_format_extract_RGGBchan(char *ID_name, char *IDoutR_name, char *IDoutG1_name, char *IDoutG2_name, char *IDoutB_name)");
  data.NBcmd++;

  strcpy(data.cmd[data.NBcmd].key,"loadcr2torgb");
  strcpy(data.cmd[data.NBcmd].module,__FILE__);
  data.cmd[data.NBcmd].fp = IMAGE_FORMAT_loadCR2toFITSRGB_cli;
  strcpy(data.cmd[data.NBcmd].info,"load CR2 file into R G B images");
  strcpy(data.cmd[data.NBcmd].syntax,"<input image> <imR> <imG> <imB>");
  strcpy(data.cmd[data.NBcmd].example,"loadcr2torgb im imR imG imB");
  strcpy(data.cmd[data.NBcmd].Ccall,"loadCR2toFITSRGB(char *fnameCR2, char *fnameFITSr, char *fnameFITSg, char *fnameFITSb)");
  data.NBcmd++;



 // add atexit functions here

  return 0;

}



int IMAGE_FORMAT_im_to_ASCII(char *IDname, char *foutname)
{
    long ii;
    long k;
    long ID;
    FILE *fpout;
    long naxis;
    long *coord;
    long npix;
    int kOK;
	
    ID = image_ID(IDname);
    naxis = data.image[ID].md[0].naxis;
    coord = (long*) malloc(sizeof(long)*naxis);
    npix = 1;
    for(k=0; k<naxis; k++)
    {
        npix *= data.image[ID].md[0].size[k];
        coord[k] = 0;
    }

    fpout = fopen(foutname, "w");

    for(ii=0; ii<npix; ii++)
    {
        for(k=0; k<naxis; k++)
            fprintf(fpout, "%4ld ", coord[k]);
        switch ( data.image[ID].md[0].atype ) {
        case CHAR:
            fprintf(fpout, " %5d\n", data.image[ID].array.C[ii]);
            break;
        case INT:
            fprintf(fpout, " %5d\n", data.image[ID].array.I[ii]);
            break;
        case FLOAT:
            fprintf(fpout, " %f\n", data.image[ID].array.F[ii]);
            break;
        case DOUBLE:
            fprintf(fpout, " %lf\n", data.image[ID].array.D[ii]);
            break;
        }
        coord[0]++;

        k = 0;
        kOK = 0;
        while((kOK==0)&&(k<naxis))
        {
            if(coord[k]==data.image[ID].md[0].size[k])
                {
					coord[k] = 0;
					coord[k+1]++;
				}
            else
                kOK = 1;
            k++;
        }
    }
    fclose(fpout);

    return 0;
}




/* This function is for byte swapping on big endian systems */
uint16 setUint16(uint16 x)
{
	if (BMP_BIG_ENDIAN)
		return (x & 0x00FF) << 8 | (x & 0xFF00) >> 8;
	else 
		return x;
}

/* This function is for byte swapping on big endian systems */
uint32 setUint32(uint32 x)
{
	if (BMP_BIG_ENDIAN)
		return (x & 0x000000FF) << 24 | (x & 0x0000FF00) << 8 | (x & 0x00FF0000) >> 8 | (x & 0xFF000000) >> 24;
	else 
		return x;
}

/*	This function writes out a 24-bit Windows bitmap file that is readable by Microsoft Paint.  
	The image data is a 1D array of (r, g, b) triples, where individual (r, g, b) values can 
	each take on values between 0 and 255, inclusive.
  
   The input to the function is:
	char *filename:					A string representing the filename that will be written
	uint32 width:					The width, in pixels, of the bitmap
	uint32 height:					The height, in pixels, of the bitmap
	unsigned char *image:				The image data, where each pixel is 3 unsigned chars (r, g, b)

   Written by Greg Slabaugh (slabaugh@ece.gatech.edu), 10/19/00
*/
uint32 write24BitBmpFile(char *filename, uint32 width, uint32 height, unsigned char *image)
{
	BITMAPINFOHEADER bmpInfoHeader;
	BITMAPFILEHEADER bmpFileHeader;
	FILE *filep;
	uint32 row, column;
	uint32 extrabytes, bytesize;
	unsigned char *paddedImage = NULL, *paddedImagePtr, *imagePtr;

	extrabytes = (4 - (width * 3) % 4) % 4;

        /* This is the size of the padded bitmap */
        bytesize = (width * 3 + extrabytes) * height;

        /* Fill the bitmap file header structure */
        bmpFileHeader.bfType = setUint16(BM);   /* Bitmap header */
        bmpFileHeader.bfSize = setUint32(0);      /* This can be 0 for BI_RGB bitmaps */
        bmpFileHeader.bfReserved1 = setUint16(0);
        bmpFileHeader.bfReserved2 = setUint16(0);
        bmpFileHeader.bfOffBits = setUint32(sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER));

        /* Fill the bitmap info structure */
        bmpInfoHeader.biSize = setUint32(sizeof(BITMAPINFOHEADER));
        bmpInfoHeader.biWidth = setUint32(width);
        bmpInfoHeader.biHeight = setUint32(height);
        bmpInfoHeader.biPlanes = setUint16(1);
        bmpInfoHeader.biBitCount = setUint16(24);            /* 24 - bit bitmap */
        bmpInfoHeader.biCompression = setUint32(BI_RGB);
        bmpInfoHeader.biSizeImage = setUint32(bytesize);     /* includes padding for 4 byte alignment */
        bmpInfoHeader.biXPelsPerMeter = setUint32(0);
        bmpInfoHeader.biYPelsPerMeter = setUint32(0);
        bmpInfoHeader.biClrUsed = setUint32(0);
        bmpInfoHeader.biClrImportant = setUint32(0);


        /* Open file */
        if ((filep = fopen(filename, "wb")) == NULL) {
                printf("Error opening file %s\n", filename);
                return BMP_FALSE;
        }

        /* Write bmp file header */
        if (fwrite(&bmpFileHeader, 1, sizeof(BITMAPFILEHEADER), filep) < sizeof(BITMAPFILEHEADER)) {
                printf("Error writing bitmap file header\n");
                fclose(filep);
                return BMP_FALSE;
        }

        /* Write bmp info header */
        if (fwrite(&bmpInfoHeader, 1, sizeof(BITMAPINFOHEADER), filep) < sizeof(BITMAPINFOHEADER)) {
                printf("Error writing bitmap info header\n");
                fclose(filep);
                return BMP_FALSE;
        }
     
        
	/* Allocate memory for some temporary storage */
	paddedImage = (unsigned char *)calloc(sizeof(unsigned char), bytesize);
	if (paddedImage == NULL) {
		printf("Error allocating memory \n");
		fclose(filep);
		return BMP_FALSE;
	}

	/* This code does three things.  First, it flips the image data upside down, as the .bmp
	format requires an upside down image.  Second, it pads the image data with extrabytes 
	number of bytes so that the width in bytes of the image data that is written to the
	file is a multiple of 4.  Finally, it swaps (r, g, b) for (b, g, r).  This is another
	quirk of the .bmp file format. */
	
	for (row = 0; row < height; row++) {
		imagePtr = image + (height - 1 - row) * width * 3;
		paddedImagePtr = paddedImage + row * (width * 3 + extrabytes);
		for (column = 0; column < width; column++) {
			*paddedImagePtr = *(imagePtr + 2);
			*(paddedImagePtr + 1) = *(imagePtr + 1);
			*(paddedImagePtr + 2) = *imagePtr;
			imagePtr += 3;
			paddedImagePtr += 3;
		}
	}

	/* Write bmp data */
	if (fwrite(paddedImage, 1, bytesize, filep) < bytesize) {
		printf("Error writing bitmap data\n");
		free(paddedImage);
		fclose(filep);
		return BMP_FALSE;
	}

	/* Close file */
	fclose(filep);
	free(paddedImage);
	return BMP_TRUE;
}


int image_writeBMP_auto(char *IDnameR, char *IDnameG, char *IDnameB, char *outname)
{
  long IDR,IDG,IDB;
  uint32 width;
  uint32 height;
  unsigned char *array;
  uint32 ii,jj;
  double minr,ming,minb,maxr,maxg,maxb;
  

  minr=img_min(IDnameR);
  ming=img_min(IDnameG);
  minb=img_min(IDnameB);
  
  maxr=img_max(IDnameR);
  maxg=img_max(IDnameG);
  maxb=img_max(IDnameB);
  
  IDR=image_ID(IDnameR);
  IDG=image_ID(IDnameG);
  IDB=image_ID(IDnameB);
  width = (uint32) data.image[IDR].md[0].size[0];
  height = (uint32) data.image[IDR].md[0].size[1];
  array = (unsigned char*) malloc(sizeof(unsigned char)*width*height*3);

  for(ii=0;ii<width;ii++)
    for(jj=0;jj<height;jj++)
      {
	array[(jj*width+ii)*3] = (unsigned char) ((data.image[IDR].array.F[(height-jj-1)*width+ii]-minr)*(255.0/(maxr-minr)));
	array[(jj*width+ii)*3+1] = (unsigned char) ((data.image[IDG].array.F[(height-jj-1)*width+ii]-ming)*(255.0/(maxg-ming)));
	array[(jj*width+ii)*3+2] = (unsigned char) ((data.image[IDB].array.F[(height-jj-1)*width+ii]-minb)*(255.0/(maxb-minb)));
      }
  write24BitBmpFile(outname,width,height,array);
  free(array);

  return(0);
}

int image_writeBMP(char *IDnameR, char *IDnameG, char *IDnameB, char *outname)
{
  long IDR,IDG,IDB;
  uint32 width;
  uint32 height;
  unsigned char *array;
  uint32 ii,jj;
  
  IDR=image_ID(IDnameR);
  IDG=image_ID(IDnameG);
  IDB=image_ID(IDnameB);
  width = (uint32) data.image[IDR].md[0].size[0];
  height = (uint32) data.image[IDR].md[0].size[1];
  array = (unsigned char*) malloc(sizeof(unsigned char)*width*height*3);

  for(ii=0;ii<width;ii++)
    for(jj=0;jj<height;jj++)
      {
	array[(jj*width+ii)*3] = (unsigned char) (data.image[IDR].array.F[(height-jj-1)*width+ii]);
	array[(jj*width+ii)*3+1] = (unsigned char) (data.image[IDG].array.F[(height-jj-1)*width+ii]);
	array[(jj*width+ii)*3+2] = (unsigned char) (data.image[IDB].array.F[(height-jj-1)*width+ii]);
      }
  write24BitBmpFile(outname,width,height,array);
  free(array);

  return(0);
}

long getImageInfo(FILE* inputFile, long offset, int numberOfChars)
{
  unsigned char			*ptrC;
  long				value=0L;
  int				i;
  unsigned char			dummy;
  int r;

  dummy = '0';
  ptrC = &dummy;

  fseek(inputFile, offset, SEEK_SET);

  for(i=1; i<=numberOfChars; i++)
  {
    r = fread(ptrC, sizeof(char), 1, inputFile);
    /* calculate value based on adding bytes */
    value = (long)(value + (*ptrC)*(pow(256, (i-1))));
  }

  return(value);
}

// ASCII format:
// ii jj value
// one line per pixel
long read_ASCIIimage(char *filename, char *ID_name, long xsize, long ysize)
{
  long ID;
  FILE *fp;
  long ii,jj;
  float value;

  ID = create_2Dimage_ID(ID_name,xsize,ysize);
   
  fp = fopen(filename,"r");
  if(fp==NULL)
    {
      fprintf(stderr,"ERROR: cannot open file \"%s\"\n",filename);
    }
  else
    {
      while((fscanf(fp,"%ld %ld %f\n",&ii,&jj,&value))==3)
	if((ii>-1)&&(ii<xsize)&&(jj>-1)&&(jj<ysize))
	  data.image[ID].array.F[jj*xsize+ii] = value;
      fclose(fp);
    }

  return(ID);
}

// ASCII format:
// value
long read_ASCIIimage1(char *filename, char *ID_name, long xsize, long ysize)
{
  long ID;
  FILE *fp;
  long ii,jj;
  double value;

  
  ID = create_2Dimage_ID(ID_name,xsize,ysize);
  
  fp = fopen(filename,"r");
  if(fp==NULL)
    {
      fprintf(stderr,"ERROR: cannot open file \"%s\"\n",filename);
    }
  else
    {
      for(ii=0;ii<xsize;ii++)
	for(jj=0;jj<ysize;jj++)
	  {
	    if(fscanf(fp,"%lf", &value)==1)
	      data.image[ID].array.F[jj*xsize+ii] = value;
	    else
	      {
		printERROR(__FILE__,__func__,__LINE__,"read error");
		exit(0);
	      }
	  }
      fclose(fp);
    }

  return(ID);
}



int read_BMPimage(char* filename, char *IDname_R, char *IDname_G, char *IDname_B)
{
  FILE				*bmpInput, *rasterOutput;
  sImage			originalImage;
  unsigned char			someChar;
  unsigned char			*pChar;
  long				fileSize;
  int				nColors;
  int				r, c;
  unsigned int BlueValue,RedValue,GreenValue;
  long IDR,IDG,IDB;

  /*--------INITIALIZE POINTER----------*/
  someChar = '0';
  pChar = &someChar;

  printf("Reading file %s\n", filename);

  /*----DECLARE INPUT AND OUTPUT FILES----*/
  bmpInput = fopen(filename, "rb");
  rasterOutput = fopen("data24.txt", "w");

  fseek(bmpInput, 0L, SEEK_END);

  /*-----GET BMP INFO-----*/
  originalImage.cols = (int)getImageInfo(bmpInput, 18, 4)+1;
  originalImage.rows = (int)getImageInfo(bmpInput, 22, 4);
  fileSize = getImageInfo(bmpInput, 2, 4);
  nColors = getImageInfo(bmpInput, 46, 4);

  /*----PRINT BMP INFO TO SCREEN-----*/
  printf("Width: %d\n", originalImage.cols);
  printf("Height: %d\n", originalImage.rows);
  printf("File size: %ld\n", fileSize);
  printf("Bits/pixel: %ld\n", getImageInfo(bmpInput, 28, 4));
  printf("No. colors: %d\n", nColors);


  IDR = create_2Dimage_ID(IDname_R,(long) originalImage.cols,(long) originalImage.rows);
  IDG = create_2Dimage_ID(IDname_G,(long) originalImage.cols,(long) originalImage.rows);
  IDB = create_2Dimage_ID(IDname_B,(long) originalImage.cols,(long) originalImage.rows);

  /*----FOR 24-BIT BMP, THERE IS NO COLOR TABLE-----*/
  fseek(bmpInput, 54, SEEK_SET);

  /*-----------READ RASTER DATA-----------*/
  for(r=0; r<=originalImage.rows-1; r++)
  {
    for(c=0; c<=originalImage.cols-1; c++)
    {
      
      /*----READ FIRST BYTE TO GET BLUE VALUE-----*/
      r = fread(pChar, sizeof(char), 1, bmpInput);
      BlueValue = *pChar;

      /*-----READ NEXT BYTE TO GET GREEN VALUE-----*/
      r = fread(pChar, sizeof(char), 1, bmpInput);
      GreenValue = *pChar;

      /*-----READ NEXT BYTE TO GET RED VALUE-----*/
      r = fread(pChar, sizeof(char), 1, bmpInput);
      RedValue = *pChar;

      /*---------WRITE TO FILES ---------*/
      /*fprintf(rasterOutput, "(%d %d) = \tRed \t%d", r, c, RedValue);
	fprintf(rasterOutput, "\tGreen \t%d \tBlue \t%d\n", GreenValue, BlueValue);*/
      data.image[IDR].array.F[r*originalImage.cols+c] = 1.0*RedValue;
      data.image[IDG].array.F[r*originalImage.cols+c] = 1.0*GreenValue;
      data.image[IDB].array.F[r*originalImage.cols+c] = 1.0*BlueValue;
    }
  }

  fclose(bmpInput);
  fclose(rasterOutput);

  return 0;
}

// reads PGM images (16 bit only)
// written to read output of "dcraw -t 0 -D -4 xxx.CR2" into FITS
int read_PGMimage(char *fname, char *ID_name)
{
    FILE *fp;
    char line1[100];
    long xsize,ysize;
    long maxval;
    long ID;
    double val;
    long ii,jj;
    int r;

    if((fp=fopen(fname,"r"))==NULL)
    {
        fprintf(stderr,"ERROR: cannot open file \"%s\"\n",fname);
    }
    else
    {
        r = fscanf(fp,"%s",line1);
        if(strcmp(line1,"P5")!=0)
            fprintf(stderr,"ERROR: File is not PGM image\n");
        else
        {
            r = fscanf(fp,"%ld %ld",&xsize,&ysize);
            printf("PGM image size: %ld x %ld\n",xsize,ysize);
            r = fscanf(fp,"%ld",&maxval);
            if(maxval!=65535)
                fprintf(stderr,"Not 16-bit image. Cannot read\n");
            else
            {
                printf("Reading PGM image\n");
                ID = create_2Dimage_ID(ID_name,xsize,ysize);
                fgetc(fp);
                for(jj=0; jj<ysize; jj++)
                {
                    for(ii=0; ii<xsize; ii++)
                    {
                        val = 256.0*((int) fgetc(fp)) + 1.0*((int) fgetc(fp));
                        data.image[ID].array.F[(ysize-jj-1)*xsize+ii] = val;
                    }
                }
            }
        }
        fclose(fp);
    }

    return(0);
}



// assumes dcraw is installed
int CR2toFITS(char *fnameCR2, char *fnameFITS)
{
    char command[200];
    FILE *fp;

    float iso;
    float shutter;
    float aperture;
    long ID;
    long xsize,ysize;
    long ii;
    int r;

    sprintf(command,"dcraw -t 0 -D -4 -c %s > _tmppgm.pgm",fnameCR2);
    r = system(command);


    read_PGMimage("_tmppgm.pgm","tmpfits1");
    r = system("rm _tmppgm.pgm");

    if(CR2toFITS_NORM==1)
    {
        sprintf(command,"dcraw -i -v %s | grep \"ISO speed\"| awk '{print $3}' > iso_tmp.txt",fnameCR2);
        r = system(command);
        fp = fopen("iso_tmp.txt","r");
        r = fscanf(fp,"%f\n",&iso);
        fclose(fp);
        r = system("rm iso_tmp.txt");
        printf("iso = %f\n",iso);

        sprintf(command,"dcraw -i -v %s | grep \"Shutter\"| awk '{print $2}' > shutter_tmp.txt",fnameCR2);
        r = system(command);
        fp = fopen("shutter_tmp.txt","r");
        r = fscanf(fp,"%f\n",&shutter);
        fclose(fp);
        r = system("rm shutter_tmp.txt");
        printf("shutter = %f\n",shutter);

        sprintf(command,"dcraw -i -v %s | grep \"Aperture\"| awk '{print $2}' > aperture_tmp.txt",fnameCR2);
        r = system(command);
        fp = fopen("aperture_tmp.txt","r");
        r = fscanf(fp,"f/%f\n",&aperture);
        fclose(fp);
        r = system("rm aperture_tmp.txt");
        printf("aperture = %f\n",aperture);

        ID = image_ID("tmpfits1");
        xsize = data.image[ID].md[0].size[0];
        ysize = data.image[ID].md[0].size[1];

        for(ii=0; ii<xsize*ysize; ii++)
            data.image[ID].array.F[ii] /= (shutter*aperture*aperture*iso);
    }

    save_fl_fits("tmpfits1",fnameFITS);
    delete_image_ID("tmpfits1");

    return(0);
}




// assumes dcraw is installed
long loadCR2(char *fnameCR2, char *IDname)
{
  char command[200];
  FILE *fp;

  float iso;
  float shutter;
  float aperture;
  long ID;
  long xsize,ysize;
  long ii;
  int r;

  sprintf(command,"dcraw -t 0 -D -4 -c %s > _tmppgm.pgm",fnameCR2);
  r = system(command);

  read_PGMimage("_tmppgm.pgm",IDname);
  r = system("rm _tmppgm.pgm");

  return(ID);
}


// load all images matching strfilter + .CR2
// return number of images converted
// FITS image name = CR2 image name with .CR2 -> .fits
long CR2toFITS_strfilter(char *strfilter)
{
  long i;
  long cnt = 0;
  char command[200];
  char fname[200];
  char fname1[200];
  FILE *fp;
  int r;

  sprintf(command,"ls %s.CR2 > flist.tmp\n",strfilter);
  r = system(command);
  
  fp = fopen("flist.tmp","r");
  while(fgets(fname,200,fp)!=NULL)
    {
      fname[strlen(fname)-1] = '\0';
      strncpy(fname1,fname,strlen(fname)-4);
      fname1[strlen(fname)-4] = '.';
      fname1[strlen(fname)-3] = 'f';
      fname1[strlen(fname)-2] = 'i';
      fname1[strlen(fname)-1] = 't';
      fname1[strlen(fname)] = 's';
      fname1[strlen(fname)+1] = '\0';

      CR2toFITS(fname,fname1);
      printf("File %s  -> file %s\n",fname,fname1);
      cnt++;
    }

  fclose(fp);
  r = system("rm flist.tmp");

  printf("%ld files converted\n",cnt);

  return(cnt);
}

//
// separates a single RGB image into its 4 channels
// output written in im_r, im_g1, im_g2 and im_b
//
int image_format_extract_RGGBchan(char *ID_name, char *IDoutR_name, char *IDoutG1_name, char *IDoutG2_name, char *IDoutB_name)
{
  long ID;
  long Xsize, Ysize;
  long IDr, IDg1, IDg2, IDb;
  long xsize2,ysize2;
  long ii,jj,ii1,jj1;
  int RGBmode = 0;
  long ID00, ID01, ID10, ID11;

  ID = image_ID(ID_name);
  Xsize = data.image[ID].md[0].size[0];
  Ysize = data.image[ID].md[0].size[1];

  printf("ID = %ld\n", ID);
  printf("size = %ld %ld\n", Xsize, Ysize);


  if((Xsize == 4770)&&(Ysize == 3178))
    RGBmode = 1;
  if((Xsize == 5202)&&(Ysize == 3465))
    RGBmode = 2;



  if(RGBmode == 0)
    {
      printERROR(__FILE__,__func__,__LINE__,"Unknown RGB image mode\n");
      exit(0);
    }

  xsize2 = Xsize/2;
  ysize2 = Ysize/2;

  printf("Creating color channel images, %ld x %ld\n", xsize2, ysize2);
  fflush(stdout);


  IDr = create_2Dimage_ID(IDoutR_name, xsize2, ysize2);
  IDg1 = create_2Dimage_ID(IDoutG1_name, xsize2, ysize2);
  IDg2 = create_2Dimage_ID(IDoutG2_name, xsize2, ysize2);
  IDb = create_2Dimage_ID(IDoutB_name, xsize2, ysize2);
  
  printf("STEP 2\n");
  fflush(stdout);

  if(RGBmode==1) // GBRG
    {
      ID00 = IDg1;
      ID10 = IDb;
      ID01 = IDr;
      ID11 = IDg2;
    }
  
  if(RGBmode==2)
    {
      ID00 = IDr;
      ID10 = IDg1;
      ID01 = IDg2;
      ID11 = IDb;
    }

  
  for(ii=0;ii<xsize2;ii++)
    for(jj=0;jj<ysize2;jj++)
      {
	ii1 = 2*ii;
	jj1 = 2*jj;
	data.image[ID01].array.F[jj*xsize2+ii] = data.image[ID].array.F[(jj1+1)*Xsize+ii1];
	data.image[ID00].array.F[jj*xsize2+ii] = data.image[ID].array.F[jj1*Xsize+ii1];
	data.image[ID11].array.F[jj*xsize2+ii] = data.image[ID].array.F[(jj1+1)*Xsize+(ii1+1)];
	data.image[ID10].array.F[jj*xsize2+ii] = data.image[ID].array.F[jj1*Xsize+(ii1+1)];
      }
  
  return(0);
}

//
// assembles 4 channels into a single image (inverse operation of routine above)
//
int image_format_reconstruct_from_RGGBchan(char *IDr_name, char *IDg1_name, char *IDg2_name, char *IDb_name, char *IDout_name)
{
  long ID;
  long IDr, IDg1, IDg2, IDb;
  long xsize1, ysize1, xsize2, ysize2;
  long ii1,jj1;
  int RGBmode = 0;
  long ID00, ID01, ID10, ID11;


  IDr = image_ID(IDr_name);
  IDg1 = image_ID(IDg1_name);
  IDg2 = image_ID(IDg2_name);
  IDb = image_ID(IDb_name);
  xsize1 = data.image[IDr].md[0].size[0];
  ysize1 = data.image[IDr].md[0].size[1];

  xsize2 = 2*xsize1;
  ysize2 = 2*ysize1;

  if((xsize2 == 4770)&&(ysize2 == 3178))
    RGBmode = 1;
  if((xsize2 == 5202)&&(ysize2 == 3465))
    RGBmode = 2;

  if(RGBmode == 0)
    {
      printERROR(__FILE__,__func__,__LINE__,"Unknown RGB image mode\n");
      exit(0);
    }

  if(RGBmode==1) // GBRG
    {
      ID00 = IDg1;
      ID10 = IDb;
      ID01 = IDr;
      ID11 = IDg2;
    }
  
  if(RGBmode==2)
    {
      ID00 = IDr;
      ID10 = IDg1;
      ID01 = IDg2;
      ID11 = IDb;
    }

  ID = create_2Dimage_ID(IDout_name,xsize2,ysize2);
  
  for(ii1=0;ii1<xsize1;ii1++)
    for(jj1=0;jj1<ysize1;jj1++)
      {
	data.image[ID].array.F[(2*jj1+1)*xsize2+2*ii1] = data.image[ID01].array.F[jj1*xsize1+ii1];
	data.image[ID].array.F[2*jj1*xsize2+2*ii1] = data.image[ID00].array.F[jj1*xsize1+ii1];
	data.image[ID].array.F[(2*jj1+1)*xsize2+(2*ii1+1)] = data.image[ID11].array.F[jj1*xsize1+ii1];
	data.image[ID].array.F[2*jj1*xsize2+(2*ii1+1)] = data.image[ID10].array.F[jj1*xsize1+ii1];
      }

  return(ID);
}



// convers a single raw bayer FITS frame into RGB FITS 
// uses "bias", "badpix" and "flat" if they exist
// output is imr, img, imb
// this is a simple interpolation routine
// IMPORTANT: input will be modified
// Sampling factor : 0=full resolution (slow), 1=half resolution (fast), 2=quarter resolution (very fast)
// Fast mode does not reject bad pixels
int convert_rawbayerFITStorgbFITS_simple(char *ID_name, char *ID_name_r, char *ID_name_g, char *ID_name_b, int SamplFactor)
{
  long ID;
  long Xsize,Ysize;
  long IDr,IDg,IDb,IDrc,IDgc,IDbc,IDbp;
  long IDbadpix;
  long IDflat;
  long IDdark;
  long IDbias;
  long ii,jj,ii1,jj1,ii2,jj2,iistart,iiend,jjstart,jjend,dii,djj;
  double v1,v2,v,vc,tmp1;
  long cnt;
  double coeff;
  long ID00, ID01, ID10, ID11;
  long ID00c, ID01c, ID10c, ID11c;
  double eps = 1.0e-8;
  int RGBmode = 0;

  int FastMode = 0;

  if(variable_ID("_RGBfast")!=-1)
    FastMode = 1;

  ID = image_ID(ID_name);
  Xsize = data.image[ID].md[0].size[0];
  Ysize = data.image[ID].md[0].size[1];

  printf("X Y  = %ld %ld\n", Xsize, Ysize);
  


  if((Xsize == 4290)&&(Ysize == 2856))
    RGBmode = 1;
  if((Xsize == 4770)&&(Ysize == 3178))
    RGBmode = 1;
  if((Xsize == 5202)&&(Ysize == 3465))
    RGBmode = 2;

  if(RGBmode == 0)
    {
      printERROR(__FILE__,__func__,__LINE__,"WARNING: Unknown RGB image mode\n");
      exit(0);
    }


  printf("FAST MODE = %d\n", FastMode);
  printf("RGBmode   = %d\n", RGBmode);
  //exit(0);

  if(FastMode==0)
    {
  // bias
  IDbias = image_ID("bias");
  if(IDbias == -1)
    {
      IDbias = create_2Dimage_ID("bias",Xsize,Ysize);
      for(ii=0;ii<Xsize*Ysize;ii++)
	data.image[IDbias].array.F[ii] = 0.0;
    }

  // dark
  IDdark = image_ID("dark");
  if(IDdark == -1)
    {
      IDdark = create_2Dimage_ID("dark",Xsize,Ysize);
      for(ii=0;ii<Xsize*Ysize;ii++)
	data.image[IDdark].array.F[ii] = 0.0;
    }

  // bad pixel map
  IDbadpix = image_ID("badpix");
  if(IDbadpix == -1)
    {
      IDbadpix = create_2Dimage_ID("badpix",Xsize,Ysize);
      for(ii=0;ii<Xsize*Ysize;ii++)
	data.image[IDbadpix].array.F[ii] = 0.0;
    }

  copy_image_ID("badpix", "badpix1", 0);
  IDbp = image_ID("badpix1");
  
  // flat field
  IDflat = image_ID("flat");
  if(IDflat == -1)
    {
      IDflat = create_2Dimage_ID("flat",Xsize,Ysize);
      for(ii=0;ii<Xsize*Ysize;ii++)
	data.image[IDflat].array.F[ii] = 1.0;
      //      arith_image_cstadd_inplace("flat",1.0);
    }

  


  // remove bias 
  if(IDbias!=-1)
    {
      for(ii=0;ii<Xsize;ii++)
	for(jj=0;jj<Ysize;jj++)
	  data.image[ID].array.F[jj*Xsize+ii] -= data.image[IDbias].array.F[jj*Xsize+ii];
    }
  // remove dark 
  if(IDdark!=-1)
    {
      for(ii=0;ii<Xsize;ii++)
	for(jj=0;jj<Ysize;jj++)
	  data.image[ID].array.F[jj*Xsize+ii] -= data.image[IDdark].array.F[jj*Xsize+ii];
    }

  
    
  // remove obvious isolated hot pixels
  cnt = 0;
  for(ii=0;ii<Xsize;ii++)
    for(jj=0;jj<Ysize;jj++)
      {
	v1 = data.image[ID].array.F[jj*Xsize+ii];
	iistart = ii-2;
	iiend = ii+2;
	if(iistart<0)
	  iistart = 0;
	if(iiend>Xsize-1)
	  iiend = Xsize-1;
	jjstart = jj-2;
	jjend = jj+2;
	if(jjstart<0)
	  jjstart = 0;
	if(jjend>Ysize-1)
	  jjend = Ysize-1;
	v2 = 0.0;
	for(ii1=iistart;ii1<iiend;ii1++)
	  for(jj1=jjstart;jj1<jjend;jj1++)
	    if((ii1!=ii)||(jj1!=jj))
	      {
		tmp1 = data.image[ID].array.F[jj1*Xsize+ii1];
		if(tmp1>v2)
		  v2 = tmp1;
	      }
	if(v1>4.0*v2+500.0)
	  {
	    data.image[ID].array.F[jj*Xsize+ii] = v2;
	    //		data.image[IDbp].array.F[jj*Xsize+ii] = 1.0;
	    cnt ++;
	  }
      }
  printf("%ld hot pixels removed\n",cnt);
   
 
  for(ii=0;ii<Xsize;ii++)
    for(jj=0;jj<Ysize;jj++)
      data.image[ID].array.F[jj*Xsize+ii] *= FLUXFACTOR;
    }
  


  switch ( SamplFactor ) {

  case 0 :

    if(image_ID(ID_name_r)!=-1)
      delete_image_ID(ID_name_r);
    IDr = create_2Dimage_ID(ID_name_r,Xsize,Ysize);
    IDrc = create_2Dimage_ID("imrc",Xsize,Ysize);
    
    if(image_ID(ID_name_g)!=-1)
      delete_image_ID(ID_name_g);
    IDg = create_2Dimage_ID(ID_name_g,Xsize,Ysize);
    IDgc = create_2Dimage_ID("imgc",Xsize,Ysize);
    
    if(image_ID(ID_name_b)!=-1)
      delete_image_ID(ID_name_b);
    IDb = create_2Dimage_ID(ID_name_b,Xsize,Ysize);
    IDbc = create_2Dimage_ID("imbc",Xsize,Ysize);
 
    if(RGBmode==1) // GBRG
      {
	ID00 = IDg;
	ID00c = IDgc;
	
	ID10 = IDb;
	ID10c = IDbc;
	
	ID01 = IDr;
	ID01c = IDrc;
	
	ID11 = IDg;
	ID11c = IDgc;
      }
    
    if(RGBmode==2)
      {
	ID00 = IDr;
	ID00c = IDrc;
	
	ID10 = IDg;
	ID10c = IDgc;
	
	ID01 = IDg;
	ID01c = IDgc;
	
	ID11 = IDb;
	ID11c = IDbc;
      }

    if(FastMode==0)
      {
	for(ii1=0;ii1<Xsize/2;ii1++)
	  for(jj1=0;jj1<Ysize/2;jj1++)
	    {
	      ii = ii1*2;
	      jj = jj1*2;
	      
	      ii2 = ii;
	      jj2 = jj+1;
	      data.image[ID01].array.F[jj2*Xsize+ii2] = data.image[ID].array.F[jj2*Xsize+ii2]/data.image[IDflat].array.F[jj2*Xsize+ii2];
	      data.image[ID01c].array.F[jj2*Xsize+ii2] = 1.0-data.image[IDbp].array.F[jj2*Xsize+ii2];
	      
	      ii2 = ii+1;
	      jj2 = jj+1;
	      data.image[ID11].array.F[jj2*Xsize+ii2] = data.image[ID].array.F[jj2*Xsize+ii2]/data.image[IDflat].array.F[jj2*Xsize+ii2];
	      data.image[ID11c].array.F[jj2*Xsize+ii2] = 1.0-data.image[IDbp].array.F[jj2*Xsize+ii2];
	  
	      ii2 = ii;
	      jj2 = jj;
	      data.image[ID00].array.F[jj2*Xsize+ii2] = data.image[ID].array.F[jj2*Xsize+ii2]/data.image[IDflat].array.F[jj2*Xsize+ii2];
	      data.image[ID00c].array.F[jj2*Xsize+ii2] = 1.0-data.image[IDbp].array.F[jj2*Xsize+ii2];
	      
	      ii2 = ii+1;
	      jj2 = jj;
	      data.image[ID10].array.F[jj2*Xsize+ii2] = data.image[ID].array.F[jj2*Xsize+ii2]/data.image[IDflat].array.F[jj2*Xsize+ii2];
	      data.image[ID10c].array.F[jj2*Xsize+ii2] = 1.0-data.image[IDbp].array.F[jj2*Xsize+ii2];
	    }
      	
	for(ii=0;ii<Xsize;ii++)
	  for(jj=0;jj<Ysize;jj++)
	    {
	      if(data.image[IDrc].array.F[jj*Xsize+ii]<0.5)
		{
		  v = 0.0;
		  vc = 0.0;
		  for(dii=-2;dii<3;dii++)
		    for(djj=-2;djj<3;djj++)
		      {
			ii1 = ii+dii;
			jj1 = jj+djj;
			if((ii1>-1)&&(jj1>-1)&&(ii1<Xsize)&&(jj1<Ysize))
			  if((dii!=0)||(djj!=0))
			    {
			      if(data.image[IDrc].array.F[jj1*Xsize+ii1]>0.5)
				{
				  coeff = exp(-5.0*(dii*dii+djj*djj));
				  vc += coeff;
				  v += data.image[IDr].array.F[jj1*Xsize+ii1]*coeff;
				}
			    }
		      }
		  data.image[IDr].array.F[jj*Xsize+ii] = v/vc;
		}
	      
	      if(data.image[IDgc].array.F[jj*Xsize+ii]<0.5)
		{
		  v = 0.0;
		  vc = 0.0;
		  for(dii=-2;dii<3;dii++)
		    for(djj=-2;djj<3;djj++)
		      {
			ii1 = ii+dii;
			jj1 = jj+djj;
			if((ii1>-1)&&(jj1>-1)&&(ii1<Xsize)&&(jj1<Ysize))
			  if((dii!=0)||(djj!=0))
			    {
			      if(data.image[IDgc].array.F[jj1*Xsize+ii1]>0.5)
				{
				  coeff = exp(-5.0*(dii*dii+djj*djj));
				  vc += coeff;
				  v += data.image[IDg].array.F[jj1*Xsize+ii1]*coeff;
				}
			    }
		      }
		  data.image[IDg].array.F[jj*Xsize+ii] = v/vc;
		}
	      
	      if(data.image[IDbc].array.F[jj*Xsize+ii]<0.5)
		{
		  v = 0.0;
		  vc = 0.0;
		  for(dii=-2;dii<3;dii++)
		    for(djj=-2;djj<3;djj++)
		      {
			ii1 = ii+dii;
			jj1 = jj+djj;
			if((ii1>-1)&&(jj1>-1)&&(ii1<Xsize)&&(jj1<Ysize))
			  if((dii!=0)||(djj!=0))
			    {
			      if(data.image[IDbc].array.F[jj1*Xsize+ii1]>0.5)
				{
				  coeff = exp(-5.0*(dii*dii+djj*djj));
				  vc += coeff;
				  v += data.image[IDb].array.F[jj1*Xsize+ii1]*coeff;
				}
			    }
		      }
		  data.image[IDb].array.F[jj*Xsize+ii] = v/vc;
		}
	    }
      }
    else
      {
	if(RGBmode==1) // GBRG
	  {
	    // G
	    for(ii1=0;ii1<Xsize/2;ii1++)
	      for(jj1=0;jj1<Ysize/2;jj1++)
		{
		  ii = ii1*2;
		  jj = jj1*2;
		 
		  
		  ii2 = ii;
		  jj2 = jj;
		  data.image[IDg].array.F[jj2*Xsize+ii2] = data.image[ID].array.F[jj2*Xsize+ii2];
		  ii2 = ii+1;
		  jj2 = jj+1;
		  data.image[IDg].array.F[jj2*Xsize+ii2] = data.image[ID].array.F[jj2*Xsize+ii2];
		}
	    // replace blue pixels
	    for(ii1=0;ii1<Xsize/2-1;ii1++)
	      for(jj1=1;jj1<Ysize/2;jj1++)
		{
		  ii = ii1*2;
		  jj = jj1*2;
		 		  
		  ii2 = ii+1;
		  jj2 = jj;
		  data.image[IDg].array.F[jj2*Xsize+ii2] = 0.25*(data.image[ID].array.F[jj2*Xsize+(ii2-1)]+data.image[ID].array.F[jj2*Xsize+(ii2+1)]+data.image[ID].array.F[(jj2+1)*Xsize+ii2]+data.image[ID].array.F[(jj2-1)*Xsize+ii2]);
		}
	    // replace red pixels
	    for(ii1=1;ii1<Xsize/2;ii1++)
	      for(jj1=0;jj1<Ysize/2-1;jj1++)
		{
		  ii = ii1*2;
		  jj = jj1*2;

		  ii2 = ii;
		  jj2 = jj+1;
		  data.image[IDg].array.F[jj2*Xsize+ii2] = 0.25*(data.image[ID].array.F[jj2*Xsize+(ii2-1)]+data.image[ID].array.F[jj2*Xsize+(ii2+1)]+data.image[ID].array.F[(jj2+1)*Xsize+ii2]+data.image[ID].array.F[(jj2-1)*Xsize+ii2]);
		}
	    


	    // R
	    for(ii1=0;ii1<Xsize/2;ii1++)
	      for(jj1=0;jj1<Ysize/2;jj1++)
		{
		   ii = ii1*2;
		   jj = jj1*2;
		   ii2 = ii;
		   jj2 = jj+1;
		   data.image[IDr].array.F[jj2*Xsize+ii2] = data.image[ID].array.F[jj2*Xsize+ii2];
		}
	    // replace g1 pixels
	    for(ii1=0;ii1<Xsize/2;ii1++)
	      for(jj1=1;jj1<Ysize/2;jj1++)
		{
		  ii = ii1*2;
		  jj = jj1*2;
		 		  
		  ii2 = ii;
		  jj2 = jj;
		  data.image[IDr].array.F[jj2*Xsize+ii2] = 0.5*(data.image[ID].array.F[(jj2-1)*Xsize+ii2] + data.image[ID].array.F[(jj2+1)*Xsize+ii2]);		  
		}
	    // replace g2 pixels
	    for(ii1=0;ii1<Xsize/2-1;ii1++)
	      for(jj1=0;jj1<Ysize/2;jj1++)
		{
		  ii = ii1*2;
		  jj = jj1*2;
		 		  
		  ii2 = ii+1;
		  jj2 = jj+1;
		  data.image[IDr].array.F[jj2*Xsize+ii2] = 0.5*(data.image[ID].array.F[jj2*Xsize+(ii2-1)] + data.image[ID].array.F[jj2*Xsize+(ii2+1)]);		  
		}
	    // replace b pixels
	    for(ii1=0;ii1<Xsize/2-1;ii1++)
	      for(jj1=1;jj1<Ysize/2;jj1++)
		{
		  ii = ii1*2;
		  jj = jj1*2;
		 		  
		  ii2 = ii+1;
		  jj2 = jj;
		  data.image[IDr].array.F[jj2*Xsize+ii2] = 0.25*(data.image[ID].array.F[(jj2-1)*Xsize+(ii2-1)] + data.image[ID].array.F[(jj2-1)*Xsize+(ii2+1)]+data.image[ID].array.F[(jj2+1)*Xsize+(ii2-1)] + data.image[ID].array.F[(jj2+1)*Xsize+(ii2+1)]);		  
		}	 
   

	    // B
	    for(ii1=0;ii1<Xsize/2;ii1++)
	      for(jj1=0;jj1<Ysize/2;jj1++)
		{
		   ii = ii1*2;
		   jj = jj1*2;
		   ii2 = ii+1;
		   jj2 = jj;
		   data.image[IDb].array.F[jj2*Xsize+ii2] = data.image[ID].array.F[jj2*Xsize+ii2];
		}

	    // replace g2 pixels
	    for(ii1=0;ii1<Xsize/2;ii1++)
	      for(jj1=0;jj1<Ysize/2-1;jj1++)
		{
		  ii = ii1*2;
		  jj = jj1*2;
		 		  
		  ii2 = ii+1;
		  jj2 = jj+1;
		  data.image[IDb].array.F[jj2*Xsize+ii2] = 0.5*(data.image[ID].array.F[(jj2-1)*Xsize+ii2] + data.image[ID].array.F[(jj2+1)*Xsize+ii2]);		  
		}
	    // replace g1 pixels
	    for(ii1=1;ii1<Xsize/2;ii1++)
	      for(jj1=0;jj1<Ysize/2;jj1++)
		{
		  ii = ii1*2;
		  jj = jj1*2;
		 		  
		  ii2 = ii;
		  jj2 = jj;
		  data.image[IDb].array.F[jj2*Xsize+ii2] = 0.5*(data.image[ID].array.F[jj2*Xsize+(ii2-1)] + data.image[ID].array.F[jj2*Xsize+(ii2+1)]);		  
		}
	    // replace r pixels
	    for(ii1=1;ii1<Xsize/2;ii1++)
	      for(jj1=0;jj1<Ysize/2-1;jj1++)
		{
		  ii = ii1*2;
		  jj = jj1*2;
		 		  
		  ii2 = ii;
		  jj2 = jj+1;
		  data.image[IDb].array.F[jj2*Xsize+ii2] = 0.25*(data.image[ID].array.F[(jj2-1)*Xsize+(ii2-1)] + data.image[ID].array.F[(jj2-1)*Xsize+(ii2+1)]+data.image[ID].array.F[(jj2+1)*Xsize+(ii2-1)] + data.image[ID].array.F[(jj2+1)*Xsize+(ii2+1)]);		  
		}	 

	  }
      }



    //  delete_image_ID("badpix1");
    
    delete_image_ID("imrc");
    delete_image_ID("imgc");
    delete_image_ID("imbc");
    //  delete_image_ID("imraw");
    break;
    
  case 1:
    if(image_ID(ID_name_r)!=-1)
      delete_image_ID(ID_name_r);
    IDr = create_2Dimage_ID(ID_name_r,Xsize/2,Ysize/2);
    IDrc = create_2Dimage_ID("imrc",Xsize/2,Ysize/2);
    
    if(image_ID(ID_name_g)!=-1)
      delete_image_ID(ID_name_g);
    IDg = create_2Dimage_ID(ID_name_g,Xsize/2,Ysize/2);
    IDgc = create_2Dimage_ID("imgc",Xsize/2,Ysize/2);
    
    if(image_ID(ID_name_b)!=-1)
      delete_image_ID(ID_name_b);
    IDb = create_2Dimage_ID(ID_name_b,Xsize/2,Ysize/2);
    IDbc = create_2Dimage_ID("imbc",Xsize/2,Ysize/2);
 
    if(RGBmode==1) // GBRG
      {
	ID00 = IDg;
	ID00c = IDgc;
	
	ID10 = IDb;
	ID10c = IDbc;
	
	ID01 = IDr;
	ID01c = IDrc;
	
	ID11 = IDg;
	ID11c = IDgc;
      }
    
    if(RGBmode==2)
      {
	ID00 = IDr;
	ID00c = IDrc;
	
	ID10 = IDg;
	ID10c = IDgc;
	
	ID01 = IDg;
	ID01c = IDgc;
	
	ID11 = IDb;
	ID11c = IDbc;
      }

    for(ii1=0;ii1<Xsize/2;ii1++)
      for(jj1=0;jj1<Ysize/2;jj1++)
	{
	  ii = ii1*2;
	  jj = jj1*2;
	  
	  ii2 = ii;
	  jj2 = jj+1;
	  data.image[ID01].array.F[jj1*Xsize/2+ii1] += data.image[ID].array.F[jj2*Xsize+ii2]/data.image[IDflat].array.F[jj2*Xsize+ii2];
	  data.image[ID01c].array.F[jj1*Xsize/2+ii1] += 1.0-data.image[IDbp].array.F[jj2*Xsize+ii2];
	  
	  ii2 = ii+1;
	  jj2 = jj+1;
	  data.image[ID11].array.F[jj1*Xsize/2+ii1] += data.image[ID].array.F[jj2*Xsize+ii2]/data.image[IDflat].array.F[jj2*Xsize+ii2];
	  data.image[ID11c].array.F[jj1*Xsize/2+ii1] += 1.0-data.image[IDbp].array.F[jj2*Xsize+ii2];
	  
	  ii2 = ii;
	  jj2 = jj;
	  data.image[ID00].array.F[jj1*Xsize/2+ii1] += data.image[ID].array.F[jj2*Xsize+ii2]/data.image[IDflat].array.F[jj2*Xsize+ii2];
	  data.image[ID00c].array.F[jj1*Xsize/2+ii1] += 1.0-data.image[IDbp].array.F[jj2*Xsize+ii2];
	  
	  ii2 = ii+1;
	  jj2 = jj;
	  data.image[ID10].array.F[jj1*Xsize/2+ii1] += data.image[ID].array.F[jj2*Xsize+ii2]/data.image[IDflat].array.F[jj2*Xsize+ii2];
	  data.image[ID10c].array.F[jj1*Xsize/2+ii1] += 1.0-data.image[IDbp].array.F[jj2*Xsize+ii2];	  
	
	  data.image[IDr].array.F[jj1*Xsize/2+ii1] /= data.image[IDrc].array.F[jj1*Xsize/2+ii1]+eps;
	  data.image[IDg].array.F[jj1*Xsize/2+ii1] /= data.image[IDgc].array.F[jj1*Xsize/2+ii1]+eps;
	  data.image[IDb].array.F[jj1*Xsize/2+ii1] /= data.image[IDbc].array.F[jj1*Xsize/2+ii1]+eps;
	}
    
    delete_image_ID("imrc");
    delete_image_ID("imgc");
    delete_image_ID("imbc");

    break;
 
  }


  return(0);
}










// assumes dcraw is installed
int loadCR2toFITSRGB(char *fnameCR2, char *fnameFITSr, char *fnameFITSg, char *fnameFITSb)
{
  FILE *fp;
  char command[200];
  float iso;
  float shutter;
  float aperture;
  long ID;
  long xsize,ysize;
  long ii;
  int r;


  sprintf(command,"dcraw -t 0 -D -4 -c %s > _tmppgm.pgm",fnameCR2);
  r = system(command);
  read_PGMimage("_tmppgm.pgm","tmpfits1");
//  r = system("rm _tmppgm.pgm");
  

  
  if(CR2toFITS_NORM==1)
    {
      sprintf(command,"dcraw -i -v %s | grep \"ISO speed\"| awk '{print $3}' > iso_tmp.txt",fnameCR2);
      r = system(command);
      fp = fopen("iso_tmp.txt","r");
      r = fscanf(fp,"%f\n",&iso);
      fclose(fp);
      r = system("rm iso_tmp.txt");
      printf("iso = %f\n",iso);
      
      sprintf(command,"dcraw -i -v %s | grep \"Shutter\"| awk '{print $2}' > shutter_tmp.txt",fnameCR2);
      r = system(command);
      fp = fopen("shutter_tmp.txt","r");
      r = fscanf(fp,"%f\n",&shutter);
      fclose(fp);
      r = system("rm shutter_tmp.txt");
      printf("shutter = %f\n",shutter);
      
      sprintf(command,"dcraw -i -v %s | grep \"Aperture\"| awk '{print $2}' > aperture_tmp.txt",fnameCR2);
      r = system(command);
      fp = fopen("aperture_tmp.txt","r");
      r = fscanf(fp,"f/%f\n",&aperture);
      fclose(fp);
      r = system("rm aperture_tmp.txt");
      printf("aperture = %f\n",aperture);

      ID = image_ID("tmpfits1");
      xsize = data.image[ID].md[0].size[0];
      ysize = data.image[ID].md[0].size[1];
      
      FLUXFACTOR = aperture*aperture/(shutter*iso);
    }
  else
    FLUXFACTOR = 1.0;

  printf("FLUXFACTOR = %g\n" ,FLUXFACTOR);

  if(variable_ID("RGBfullres")==-1)
    convert_rawbayerFITStorgbFITS_simple("tmpfits1", fnameFITSr, fnameFITSg, fnameFITSb, 1);
  else
    convert_rawbayerFITStorgbFITS_simple("tmpfits1", fnameFITSr, fnameFITSg, fnameFITSb, 0);

  delete_image_ID("tmpfits1");

  FLUXFACTOR = 1.0;

  return(0);
}






// convers a single raw bayer FITS frame into RGB FITS 
// uses "badpix" and "flat" if they exist
// output is imr, img, imb
// this is a more fancy algorithm that tries to match the data with an model with min information content
// works well for fields with lots of stars
/*
int convert_rawbayerFITStorgbFITS_mininfocontent(char *ID_name)
{
  long ID;
  long Xsize,Ysize;
  long PSFsamplingzoomfactor = 4;
  long PSFsize = 1024;
  long NBsource = 50;
    
  ID = image_ID(ID_name);
  Xsize = data.image[ID].md[0].size[0];
  Ysize = data.image[ID].md[0].size[1];
  
  

  return(0);
}
*/



int CR2tomov()
{
  char configfile[200];
  long ID,IDr,IDg,IDb;
  long ii,i;
  long cnt = 0;
  long cntmax;
  char command[200];
  char fname[200];
  char fnamer[200];
  char fnameg[200];
  char fnameb[200];
  char fnamestat[200];
  char fnameoutr[200];
  char fnameoutg[200];
  char fnameoutb[200];
  char fnamejpg[200];
  FILE *fp;
  FILE *fp1;
  
  long xsize,ysize;
  
  long IDrtot;
  long IDgtot;
  long IDbtot;

  double tot;

  double alpha = 0.7;
  double vr,vg,vb,v0r,v0g,v0b,v1r,v1g,v1b,tmp1,v0,r0,g0,b0;
  double tmpr,tmpg,tmpb,tmpr1,tmpg1,tmpb1;
  
  // conversion from CR2 to FITS RGB
  int CR2toFITSrgb;
  int CR2TOFITSRGB_FORCE;
  long maxnbFITSfiles;
  int binfact;


  // conversion from FITS RGB to JPEG
  int FITStoJPEG;
  double MINLEVEL;
  double MAXLEVEL;

  int MAXLEVEL_AUTO;
  double MAXLEVEL_AUTO_FLOOR;
  double MAXLEVEL_AUTO_CEIL;

  double MAX_PERC01_COEFF;
  double MAX_PERC05_COEFF;
  double MAX_PERC10_COEFF;
  double MAX_PERC20_COEFF;
  double MAX_PERC50_COEFF;
  double MAX_PERC80_COEFF;
  double MAX_PERC90_COEFF;
  double MAX_PERC95_COEFF;
  double MAX_PERC99_COEFF;
  double MAX_PERC995_COEFF;
  double MAX_PERC998_COEFF;
  double MAX_PERC999_COEFF;

  double RGBM_RR;
  double RGBM_RG;
  double RGBM_RB;
  double RGBM_GR;
  double RGBM_GG;
  double RGBM_GB;
  double RGBM_BR;
  double RGBM_BG;
  double RGBM_BB;
  double LUMR,LUMG,LUMB; // luminance vector
  double COLORSAT;

  double ALPHA;




  double vp01,vp05,vp10,vp20,vp50,vp80,vp90,vp95,vp99,vp995,vp998,vp999;
  double vp01r,vp05r,vp10r,vp20r,vp50r,vp80r,vp90r,vp95r,vp99r,vp995r,vp998r,vp999r;
  double vp01g,vp05g,vp10g,vp20g,vp50g,vp80g,vp90g,vp95g,vp99g,vp995g,vp998g,vp999g;
  double vp01b,vp05b,vp10b,vp20b,vp50b,vp80b,vp90b,vp95b,vp99b,vp995b,vp998b,vp999b;
  double *maxlevel;
  double *maxlevel1;
  double *array;
  double value,valuecnt;
  long boxsize;
  double sigma;
  long j,jstart,jend,j1;

  int NLCONV;
  double  NLCONV_OFFSET;
  double  NLCONV_LIMIT;
  double  NLCONV_FACT;
  double NLCONV_POW;
  double NLCONV_SIGMA;
  long IDr1,IDg1,IDb1,IDrp,IDgp,IDbp;

  long SKIP, SKIPcnt;
  long SKIP_FITStoJPEG, SKIPcnt_FITStoJPEG;

  int MKim;

  int r;

  sprintf(configfile,"cr2tojpegconf.txt");

  CR2toFITSrgb = read_config_parameter_int(configfile,"CR2TOFITSRGB");
  CR2TOFITSRGB_FORCE = read_config_parameter_int(configfile,"CR2TOFITSRGB_FORCE");
  maxnbFITSfiles = read_config_parameter_long(configfile,"CR2TOFITS_MAXNBFILE");
  binfact = read_config_parameter_int(configfile,"CR2TOFITSBIN");

  maxlevel = (double*) malloc(sizeof(double)*maxnbFITSfiles);

  FITStoJPEG = read_config_parameter_int(configfile,"FITStoJPEG");
  MINLEVEL = read_config_parameter_float(configfile,"MINLEVEL");
  MAXLEVEL = read_config_parameter_float(configfile,"MAXLEVEL");

  MAXLEVEL_AUTO = read_config_parameter_int(configfile,"MAXLEVEL_AUTO");
  MAXLEVEL_AUTO_FLOOR = read_config_parameter_float(configfile,"MAXLEVEL_AUTO_FLOOR");

  if(read_config_parameter_exists(configfile,"MAXLEVEL_AUTO_CEIL")==1)
    MAXLEVEL_AUTO_CEIL = read_config_parameter_float(configfile,"MAXLEVEL_AUTO_CEIL");
  else
    MAXLEVEL_AUTO_CEIL = 100000.0;

  MAX_PERC01_COEFF = read_config_parameter_float(configfile,"MAX_PERC01_COEFF");
  MAX_PERC05_COEFF = read_config_parameter_float(configfile,"MAX_PERC05_COEFF");
  MAX_PERC10_COEFF = read_config_parameter_float(configfile,"MAX_PERC10_COEFF");
  MAX_PERC20_COEFF = read_config_parameter_float(configfile,"MAX_PERC20_COEFF");
  MAX_PERC50_COEFF = read_config_parameter_float(configfile,"MAX_PERC50_COEFF");
  MAX_PERC80_COEFF = read_config_parameter_float(configfile,"MAX_PERC80_COEFF");
  MAX_PERC90_COEFF = read_config_parameter_float(configfile,"MAX_PERC90_COEFF");
  MAX_PERC90_COEFF = read_config_parameter_float(configfile,"MAX_PERC95_COEFF");
  MAX_PERC99_COEFF = read_config_parameter_float(configfile,"MAX_PERC99_COEFF");
  MAX_PERC995_COEFF = read_config_parameter_float(configfile,"MAX_PERC995_COEFF");
  MAX_PERC998_COEFF = read_config_parameter_float(configfile,"MAX_PERC998_COEFF");
  MAX_PERC999_COEFF = read_config_parameter_float(configfile,"MAX_PERC999_COEFF");


  RGBM_RR = read_config_parameter_float(configfile,"RGBM_RR");
  RGBM_RG = read_config_parameter_float(configfile,"RGBM_RG");
  RGBM_RB = read_config_parameter_float(configfile,"RGBM_RB");
  RGBM_GR = read_config_parameter_float(configfile,"RGBM_GR");
  RGBM_GG = read_config_parameter_float(configfile,"RGBM_GG");
  RGBM_GB = read_config_parameter_float(configfile,"RGBM_GB");
  RGBM_BR = read_config_parameter_float(configfile,"RGBM_BR");
  RGBM_BG = read_config_parameter_float(configfile,"RGBM_BG");
  RGBM_BB = read_config_parameter_float(configfile,"RGBM_BB");
  LUMR = read_config_parameter_float(configfile,"LUMR");
  LUMG = read_config_parameter_float(configfile,"LUMG");
  LUMB = read_config_parameter_float(configfile,"LUMB");
  COLORSAT = read_config_parameter_float(configfile,"COLORSAT");


  NLCONV = 0;
  if(read_config_parameter_exists(configfile,"NLCONV")==1)
    {
      NLCONV = read_config_parameter_int(configfile,"NLCONV");
      NLCONV_OFFSET = read_config_parameter_float(configfile,"NLCONV_OFFSET");
      NLCONV_LIMIT = read_config_parameter_float(configfile,"NLCONV_LIMIT");
      NLCONV_FACT = read_config_parameter_float(configfile,"NLCONV_FACT");
      NLCONV_POW = read_config_parameter_float(configfile,"NLCONV_POW");
      NLCONV_SIGMA = read_config_parameter_float(configfile,"NLCONV_SIGMA");
    }

  ALPHA = read_config_parameter_float(configfile,"ALPHA");
  
  

  SKIP = 0;

  ID = variable_ID("SKIP");
  if(ID!=1)
    SKIP = (long) (data.variable[ID].value.f + 0.1);
  printf("SKIP = %ld\n",SKIP);



  if(CR2toFITSrgb==1)
    {
      load_fits("bias.fits", "bias", 1);
      load_fits("dark.fits", "dark", 1);
      load_fits("badpix.fits", "badpix", 1);
      load_fits("flat.fits", "flat", 1);
      
      sprintf(command,"ls ./CR2/*.CR2 > flist.tmp\n");
      r = system(command);
      
      fp = fopen("flist.tmp","r");
      SKIPcnt = 0;
      while((fgets(fname,200,fp)!=NULL)&&(cnt<maxnbFITSfiles))
	{
	  sprintf(fnamestat,"./FITS/imgstats.%05ld.txt",cnt);
	  sprintf(fnameoutr,"./FITS/imr%05ld.fits",cnt);
	  sprintf(fnameoutg,"./FITS/img%05ld.fits",cnt);
	  sprintf(fnameoutb,"./FITS/imb%05ld.fits",cnt);

	  MKim = 0;
	  if((file_exists(fnameoutr)==1)&&(file_exists(fnameoutg)==1)&&(file_exists(fnameoutb)==1)&&(CR2TOFITSRGB_FORCE==0))
	    printf("Files %s %s %s exist, no need to recreate\n",fnameoutr,fnameoutg,fnameoutb);
	  else
	    {
	      if(SKIPcnt==0)
		{
		  MKim = 1;
		  printf("[%ld] working on file %s\n",cnt,fname);
		  fname[strlen(fname)-1] = '\0';
		  loadCR2toFITSRGB(fname,"imr","img","imb");
		  /*		  if(binfact!=1)
		    {
		      basic_contract("imr","imrc",binfact,binfact);
		      delete_image_ID("imr");
		      chname_image_ID("imrc","imr");
		      basic_contract("img","imgc",binfact,binfact);
		      delete_image_ID("img");
		      chname_image_ID("imgc","img");
		      basic_contract("imb","imbc",binfact,binfact);
		      delete_image_ID("imb");
		      chname_image_ID("imbc","imb");
		      }*/
		  ID = image_ID("imr");
		  xsize = data.image[ID].md[0].size[0];
		  ysize = data.image[ID].md[0].size[1];
		  
		  IDrtot = image_ID("imrtot");
		  if(IDrtot==-1)
		    {
		      IDrtot = create_2Dimage_ID("imrtot",xsize,ysize);
		      IDgtot = create_2Dimage_ID("imgtot",xsize,ysize);
		      IDbtot = create_2Dimage_ID("imbtot",xsize,ysize);
		    }
		  
		  IDr = image_ID("imr");
		  IDg = image_ID("img");
		  IDb = image_ID("imb");
		  
		  for(ii=0;ii<xsize*ysize;ii++)
		    {
		      data.image[IDr].array.F[ii] /= binfact*binfact;
		      data.image[IDg].array.F[ii] /= binfact*binfact;
		      data.image[IDb].array.F[ii] /= binfact*binfact;
		      
		      data.image[IDrtot].array.F[ii] += data.image[IDr].array.F[ii];
		      data.image[IDgtot].array.F[ii] += data.image[IDg].array.F[ii];
		      data.image[IDbtot].array.F[ii] += data.image[IDb].array.F[ii];
		    }
		  save_fl_fits("imrtot","!imrtot.fits");
		  save_fl_fits("imgtot","!imgtot.fits");
		  save_fl_fits("imbtot","!imbtot.fits");
		  
		  sprintf(fnameoutr,"!./FITS/imr%05ld.fits",cnt);
		  save_fl_fits("imr",fnameoutr);	      
		  sprintf(fnameoutg,"!./FITS/img%05ld.fits",cnt);
		  save_fl_fits("img",fnameoutg);
		  sprintf(fnameoutb,"!./FITS/imb%05ld.fits",cnt);
		  save_fl_fits("imb",fnameoutb);
		}
	    }




	  if(((MKim == 1)||(file_exists(fnamestat)==0))&&(SKIPcnt==0))
	    {
	      printf("[%ld] working on file %s (statistics)\n",cnt,fname);
	      if(MKim == 0)
		{
		  sprintf(fnameoutr,"./FITS/imr%05ld.fits",cnt);
		  sprintf(fnameoutg,"./FITS/img%05ld.fits",cnt);
		  sprintf(fnameoutb,"./FITS/imb%05ld.fits",cnt);		  
		  load_fits(fnameoutr, "imr", 1);
		  load_fits(fnameoutg, "img", 1);
		  load_fits(fnameoutb, "imb", 1);
		}

	      info_image_stats("imr","");
	      ID = variable_ID("vp01");
	      vp01r = data.variable[ID].value.f;
	      ID = variable_ID("vp05");
	      vp05r = data.variable[ID].value.f;
	      ID = variable_ID("vp10");
	      vp10r = data.variable[ID].value.f;
	      ID = variable_ID("vp20");
	      vp20r = data.variable[ID].value.f;
	      ID = variable_ID("vp50");
	      vp50r = data.variable[ID].value.f;
	      ID = variable_ID("vp80");
	      vp80r = data.variable[ID].value.f;
	      ID = variable_ID("vp90");
	      vp90r = data.variable[ID].value.f;
	      ID = variable_ID("vp95");
	      vp95r = data.variable[ID].value.f;
	      ID = variable_ID("vp99");
	      vp99r = data.variable[ID].value.f;
	      ID = variable_ID("vp995");
	      vp995r = data.variable[ID].value.f;
	      ID = variable_ID("vp998");
	      vp998r = data.variable[ID].value.f;
	      ID = variable_ID("vp999");
	      vp999r = data.variable[ID].value.f;
	      delete_image_ID("imr");	      	      
	      
	      info_image_stats("img","");
	      ID = variable_ID("vp01");
	      vp01g = data.variable[ID].value.f;
	      ID = variable_ID("vp05");
	      vp05g = data.variable[ID].value.f;
	      ID = variable_ID("vp10");
	      vp10g = data.variable[ID].value.f;
	      ID = variable_ID("vp20");
	      vp20g = data.variable[ID].value.f;
	      ID = variable_ID("vp50");
	      vp50g = data.variable[ID].value.f;
	      ID = variable_ID("vp80");
	      vp80g = data.variable[ID].value.f;
	      ID = variable_ID("vp90");
	      vp90g = data.variable[ID].value.f;
	      ID = variable_ID("vp95");
	      vp95g = data.variable[ID].value.f;
	      ID = variable_ID("vp99");
	      vp99g = data.variable[ID].value.f;
	      ID = variable_ID("vp995");
	      vp995g = data.variable[ID].value.f;
	      ID = variable_ID("vp998");
	      vp998g = data.variable[ID].value.f;
	      ID = variable_ID("vp999");
	      vp999g = data.variable[ID].value.f;
	      delete_image_ID("img");
	      
	      info_image_stats("imb","");
	      ID = variable_ID("vp01");
	      vp01b = data.variable[ID].value.f;
	      ID = variable_ID("vp05");
	      vp05b = data.variable[ID].value.f;
	      ID = variable_ID("vp10");
	      vp10b = data.variable[ID].value.f;
	      ID = variable_ID("vp20");
	      vp20b = data.variable[ID].value.f;
	      ID = variable_ID("vp50");
	      vp50b = data.variable[ID].value.f;
	      ID = variable_ID("vp80");
	      vp80b = data.variable[ID].value.f;
	      ID = variable_ID("vp90");
	      vp90b = data.variable[ID].value.f;
	      ID = variable_ID("vp95");
	      vp95b = data.variable[ID].value.f;
	      ID = variable_ID("vp99");
	      vp99b = data.variable[ID].value.f;
	      ID = variable_ID("vp995");
	      vp995b = data.variable[ID].value.f;
	      ID = variable_ID("vp998");
	      vp998b = data.variable[ID].value.f;
	      ID = variable_ID("vp999");
	      vp999b = data.variable[ID].value.f;
	      delete_image_ID("imb");
	      
	      
	      fp1 = fopen(fnamestat,"w");
	      fprintf(fp1,"%05ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",cnt,vp01r,vp05r,vp10r,vp20r,vp50r,vp80r,vp90r,vp95r,vp99r,vp995r,vp998r,vp999r,vp01g,vp05g,vp10g,vp20g,vp50g,vp80g,vp90g,vp95g,vp99g,vp995g,vp998g,vp999g,vp01b,vp05b,vp10b,vp20b,vp50b,vp80b,vp90b,vp95b,vp99b,vp995b,vp998b,vp999b);
	      fclose(fp1);
	    
	      if(MKim == 0)
		{
		  delete_image_ID("imr");
		  delete_image_ID("img");
		  delete_image_ID("imb");
		}
	    }
	     
	  if(MKim == 1)
	    {
	      delete_image_ID("imr");
	      delete_image_ID("img");
	      delete_image_ID("imb");
	    }	     
	     
	  
	  SKIPcnt++;
          if(SKIPcnt>SKIP-1)
            SKIPcnt = 0;	      
	
	  cnt++;	    
	}
      fclose(fp);
      r = system("rm flist.tmp");
      
      printf("%ld images processed\n",cnt);  
    }

  r = system("rm imgstats.txt");
  r = system("cat ./FITS/imgstats.*.txt > imgstats.txt");
  
  





  if(MAXLEVEL_AUTO == 1)
    {
      fp = fopen("imgstats.txt","r");
      while(fscanf(fp,"%05ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&cnt,&vp01r,&vp05r,&vp10r,&vp20r,&vp50r,&vp80r,&vp90r,&vp95r,&vp99r,&vp995r,&vp998r,&vp999r,&vp01g,&vp05g,&vp10g,&vp20g,&vp50g,&vp80g,&vp90g,&vp95g,&vp99g,&vp995g,&vp998g,&vp999g,&vp01b,&vp05b,&vp10b,&vp20b,&vp50b,&vp80b,&vp90b,&vp95b,&vp99b,&vp995b,&vp998b,&vp999b)==37)
	{
	  vp01 = ( vp01r + vp01g + vp01b ) / 3.0;
	  vp05 = ( vp05r + vp05g + vp05b ) / 3.0;
	  vp10 = ( vp10r + vp10g + vp10b ) / 3.0;
	  vp20 = ( vp20r + vp20g + vp20b ) / 3.0;
	  vp50 = ( vp50r + vp50g + vp50b ) / 3.0;
	  vp80 = ( vp80r + vp80g + vp80b ) / 3.0;
	  vp90 = ( vp90r + vp90g + vp90b ) / 3.0;
	  vp95 = ( vp95r + vp95g + vp95b ) / 3.0;
	  vp99 = ( vp99r + vp99g + vp99b ) / 3.0;
	  vp995 = ( vp995r + vp995g + vp995b ) / 3.0;
	  vp998 = ( vp998r + vp998g + vp998b ) / 3.0;
	  vp999 = ( vp999r + vp999g + vp999b ) / 3.0;
	  if(cnt < maxnbFITSfiles)
	    {
	      maxlevel[cnt] = vp01*MAX_PERC01_COEFF + vp05*MAX_PERC05_COEFF + vp10*MAX_PERC10_COEFF + vp20*MAX_PERC20_COEFF + vp50*MAX_PERC50_COEFF + vp80*MAX_PERC80_COEFF + vp90*MAX_PERC90_COEFF + vp95*MAX_PERC95_COEFF + vp99*MAX_PERC99_COEFF + vp995*MAX_PERC995_COEFF + vp998*MAX_PERC998_COEFF + vp999*MAX_PERC999_COEFF;
	      printf("[%ld %g]   ",cnt,maxlevel[cnt]);
	      maxlevel[cnt] = sqrt(maxlevel[cnt]*maxlevel[cnt]+MAXLEVEL_AUTO_FLOOR*MAXLEVEL_AUTO_FLOOR);
	      if(maxlevel[cnt] > MAXLEVEL_AUTO_CEIL)
		maxlevel[cnt] = MAXLEVEL_AUTO_CEIL;
	      printf("%ld -> %g\n",cnt, maxlevel[cnt]);
	    }
	}
      fclose(fp);

      if(0)
	{
	  // smooth the maxlevel in time
	  // scheme employed is running median/average
	  cntmax = maxnbFITSfiles;
	  if(cntmax > cnt+1)
	    cntmax = cnt+1;
	  
	  printf("cntmax = %ld\n",cntmax);
	  boxsize = 100;
	  if(boxsize>0.1*cntmax)
	    boxsize = (long) (0.1*cntmax);
	  sigma = 0.2*boxsize;
	  
	  if(boxsize==0)
	    boxsize = 1;
	  printf("boxsize = %ld\n",boxsize);
	  
	  
	  array = (double*) malloc(sizeof(double)*(2*boxsize+1));
	  maxlevel1 = (double*) malloc(sizeof(double)*cntmax);
	  for(i=0;i<cntmax;i++)
	    {
	      jstart = i-boxsize;
	      jend = i+boxsize+1;
	      /*	  jcent = 0;
		
	      while(jstart<0)
	      {
	      jstart++;
	      jend++;
	      jcent++;
	      }
	      while(jend>cntmax-1)
	      {
	      jstart--;
	      jend--;
	      jcent--;
	      }
	      */
	      
	      for(j=0;j<2*boxsize+1;j++)	   
		{
		  j1 = j+jstart;
		  if(j1<0)
		    j1 = 0;
		  if(j1>cntmax-1)
		    j1 = cntmax-1;

		  array[j] = maxlevel[j1];
		}
	      
	      quick_sort(array,2*boxsize+1);
	      
	      value = 0.0;
	      valuecnt = 0.0;
	      for(ii=0;ii<2*boxsize+1;ii++)
		{
		  tmp1 = 1.0*(ii-boxsize);
		  value += log10(array[ii])*exp(-tmp1*tmp1/sigma/sigma);
		  valuecnt += exp(-tmp1*tmp1/sigma/sigma);
		}
	      
	      maxlevel1[i] = pow(10.0,value/valuecnt);

	    }            
	  free(array);
	  
	  fp = fopen("maxlevel.log","w");
	  for(i=0;i<cntmax;i++)
	    {
	      printf("%ld MAXLEVEL : %g ---> %g\n",i,maxlevel[i],maxlevel1[i]);
	      fprintf(fp,"%ld %g %g\n",i,maxlevel[i],maxlevel1[i]);
	      maxlevel[i] = maxlevel1[i];		  
	    }
	  fclose(fp);
	  free(maxlevel1);      
	}
    }





  if(FITStoJPEG == 1)
    {
      printf("FITS to JPEG\n");

      SKIP_FITStoJPEG = 0;
      
      ID = variable_ID("SKIP_FITStoJPEG");
      if(ID!=1)
	SKIP_FITStoJPEG = (long) (data.variable[ID].value.f + 0.1);
      printf("SKIP FITS to JPEG = %ld\n",SKIP_FITStoJPEG);


      SKIPcnt_FITStoJPEG = 0;

      for(i=0;i<maxnbFITSfiles;i++)
	{
	  sprintf(fnamejpg,"./JPEG/im%05ld.jpg",i);
	  if(file_exists(fnamejpg)==1)
	    {
	      printf("Files %s exists, no need to recreate\n",fnamejpg);
	    }
	  else
	    {
	      sprintf(fnamer,"./FITS/imr%05ld.fits",i);
	      if(file_exists(fnamer)==1)
		{
		  if(SKIPcnt_FITStoJPEG==0)
		    {
		      printf("file %s exists\n",fnamer);
		      
		      sprintf(fnamer,"./FITS/imr%05ld.f.fits",i);
		      if(file_exists(fnamer)==1)
			IDr = load_fits(fnamer, "imr", 1);
		      else
			{
			  sprintf(fnamer,"./FITS/imr%05ld.fits",i);
			  IDr = load_fits(fnamer, "imr", 1);
			}

		      sprintf(fnameg,"./FITS/img%05ld.f.fits",i);
		      if(file_exists(fnameg)==1)
			IDg = load_fits(fnameg, "img", 1);
		      else
			{
			  sprintf(fnameg,"./FITS/img%05ld.fits",i);
			  IDg = load_fits(fnameg, "img", 1);
			}
		      
		      sprintf(fnameb,"./FITS/imb%05ld.f.fits",i);
		      if(file_exists(fnameb)==1)
			IDb = load_fits(fnameb, "imb", 1);
		      else
			{
			  sprintf(fnameb,"./FITS/imb%05ld.fits",i);
			  IDb = load_fits(fnameb, "imb", 1);
			}
		      

		      xsize = data.image[IDr].md[0].size[0];
		      ysize = data.image[IDr].md[0].size[1];
		      
		      if(MAXLEVEL_AUTO == 1)
			MAXLEVEL = maxlevel[i];
		      
		      
		      for(ii=0;ii<xsize*ysize;ii++)
			{
			  r0 = data.image[IDr].array.F[ii];
			  g0 = data.image[IDg].array.F[ii];
			  b0 = data.image[IDb].array.F[ii];
			  
			  r0 = (r0-MINLEVEL)/(MAXLEVEL-MINLEVEL);
			  g0 = (g0-MINLEVEL)/(MAXLEVEL-MINLEVEL);
			  b0 = (b0-MINLEVEL)/(MAXLEVEL-MINLEVEL);
			  
			  tmpr = r0*RGBM_RR + g0*RGBM_RG + b0*RGBM_RB;
			  tmpg = r0*RGBM_GR + g0*RGBM_GG + b0*RGBM_GB;
			  tmpb = r0*RGBM_BR + g0*RGBM_BG + b0*RGBM_BB;
			  
			  
			  tmpr1 = tmpr*((1.0-COLORSAT)*LUMR+COLORSAT) + tmpg*((1.0-COLORSAT)*LUMG) + tmpb*((1.0-COLORSAT)*LUMB);
			  tmpg1 = tmpr*((1.0-COLORSAT)*LUMR) + tmpg*((1.0-COLORSAT)*LUMG+COLORSAT) + tmpb*((1.0-COLORSAT)*LUMB);
			  tmpb1 = tmpr*((1.0-COLORSAT)*LUMR) + tmpg*((1.0-COLORSAT)*LUMG) + tmpb*((1.0-COLORSAT)*LUMB+COLORSAT);
			  
			  
			  data.image[IDr].array.F[ii] = tmpr1;
			  data.image[IDg].array.F[ii] = tmpg1;
			  data.image[IDb].array.F[ii] = tmpb1;
			}
		      
		      
		      for(ii=0;ii<xsize*ysize;ii++)
			{
			  vr = data.image[IDr].array.F[ii];
			  vg = data.image[IDg].array.F[ii];
			  vb = data.image[IDb].array.F[ii];
			  
			  if(vr<0.0)
			    vr = 0.0;
			  if(vg<0.0)
			    vg = 0.0;
			  if(vb<0.0)
			    vb = 0.0;
			}
		      
		      // non-linear convolution
		      
		      /*	      if(NLCONV==1)
			{
			  printf("NLCONV_OFFSET = %f\n",NLCONV_OFFSET);
			  printf("NLCONV_LIMIT = %f\n",NLCONV_LIMIT);
			  printf("NLCONV_FACT = %f\n",NLCONV_FACT);
			  printf("NLCONV_POW = %f\n",NLCONV_POW);
			  printf("NLCONV_SIGMA = %f\n",NLCONV_SIGMA);
			  
			  IDrp = create_2Dimage_ID("imrp",xsize,ysize);
			  IDgp = create_2Dimage_ID("imgp",xsize,ysize);
			  IDbp = create_2Dimage_ID("imbp",xsize,ysize);
			  
			  for(ii=0;ii<xsize*ysize;ii++)
			    {			 
			      if(data.image[IDr].array.F[ii]>NLCONV_LIMIT)
				data.image[IDr].array.F[ii] = NLCONV_LIMIT;
			      if(data.image[IDg].array.F[ii]>NLCONV_LIMIT)
				data.image[IDg].array.F[ii] = NLCONV_LIMIT;
			      if(data.image[IDb].array.F[ii]>NLCONV_LIMIT)
				data.image[IDb].array.F[ii] = NLCONV_LIMIT;
			      
			      
			      if(data.image[IDr].array.F[ii]>NLCONV_OFFSET)
				data.image[IDrp].array.F[ii] = pow((data.image[IDr].array.F[ii]-NLCONV_OFFSET)*NLCONV_FACT,NLCONV_POW);
			      if(data.image[IDg].array.F[ii]>NLCONV_OFFSET)
				data.image[IDgp].array.F[ii] = pow((data.image[IDg].array.F[ii]-NLCONV_OFFSET)*NLCONV_FACT,NLCONV_POW);
			      if(data.image[IDb].array.F[ii]>NLCONV_OFFSET)
				data.image[IDbp].array.F[ii] = pow((data.image[IDb].array.F[ii]-NLCONV_OFFSET)*NLCONV_FACT,NLCONV_POW);
			    }		      
			  make_gauss("kerg",xsize,ysize,NLCONV_SIGMA,1.0);
			  tot = arith_image_total("kerg");
			  arith_image_cstmult_inplace("kerg",1.0/tot);
			  
			  fconvolve_padd("imrp","kerg",(long) (10.0*NLCONV_SIGMA+2.0),"imr_c");
			  fconvolve_padd("imgp","kerg",(long) (10.0*NLCONV_SIGMA+2.0),"img_c");
			  fconvolve_padd("imbp","kerg",(long) (10.0*NLCONV_SIGMA+2.0),"imb_c");
			  delete_image_ID("kerg");
			  delete_image_ID("imrp");
			  delete_image_ID("imgp");
			  delete_image_ID("imbp");
			  IDr1 = image_ID("imr_c");
			  IDg1 = image_ID("img_c");
			  IDb1 = image_ID("imb_c");
			  for(ii=0;ii<xsize*ysize;ii++)
			    {			 
			      data.image[IDr].array.F[ii] += data.image[IDr1].array.F[ii];
			      data.image[IDg].array.F[ii] += data.image[IDg1].array.F[ii];
			      data.image[IDb].array.F[ii] += data.image[IDb1].array.F[ii];
			    }		      
			  delete_image_ID("imr_c");
			  delete_image_ID("img_c");
			  delete_image_ID("imb_c");
			  }*/
		      
		      for(ii=0;ii<xsize*ysize;ii++)
			{
			  vr = data.image[IDr].array.F[ii];
			  vg = data.image[IDg].array.F[ii];
			  vb = data.image[IDb].array.F[ii];
			  
			  if(vr<0.0)
			    vr = 0.0;
			  if(vg<0.0)
			    vg = 0.0;
			  if(vb<0.0)
			    vb = 0.0;
			  
			  if(vr>1.0)
			    vr = 1.0;
			  
			  if(vg>1.0)
			    vg = 1.0;
			  
			  if(vb>1.0)
			    vb = 1.0;
			  
			  vr = 255.0*pow(vr,ALPHA);
			  vg = 255.0*pow(vg,ALPHA);
			  vb = 255.0*pow(vb,ALPHA);
			  
			  data.image[IDr].array.F[ii] = vr;
			  data.image[IDg].array.F[ii] = vg;
			  data.image[IDb].array.F[ii] = vb;
			}
		      
		      
		      image_writeBMP("imr","img","imb","imrgb.bmp");
				      delete_image_ID("imr");
		      delete_image_ID("img");
		      delete_image_ID("imb");
		      //		  sprintf(fnamejpg,"./JPEG/im%05ld.jpg",i);
		      sprintf(command,"bmptoppm imrgb.bmp | ppmtojpeg --quality 95 > _tmpjpeg.jpg; mv _tmpjpeg.jpg %s",fnamejpg);
		      r = system(command);
		      r = system("rm imrgb.bmp");
		    }
		  SKIPcnt_FITStoJPEG++;
		  if(SKIPcnt_FITStoJPEG>SKIP_FITStoJPEG-1)
		    SKIPcnt_FITStoJPEG = 0;	      	      
		}
	    }      
	}
    }
  free(maxlevel);

  return(0);
}


// ron: readout noise in ADU
// gain: in e-/ADU
// alpha: 0 - 1, sets quantization noise at alpha x overall noise
// bias: image bias in ADU
long IMAGE_FORMAT_requantize(char *IDin_name, char *IDout_name, double alpha, double ron, double gain, double bias)
{
  long IDin, IDout;
  long ii;
  long xsize,ysize;
  double value;

  IDin = image_ID(IDin_name);
  xsize = data.image[IDin].md[0].size[0];
  ysize = data.image[IDin].md[0].size[1];
  
  IDout = create_2Dimage_ID(IDout_name,xsize,ysize);
  for(ii=0;ii<xsize*ysize;ii++)
    {
      value = data.image[IDin].array.F[ii];
      value = value - bias;
      if(value < 0.0)
	value = value/(alpha*ron);
      else
	value = 2.0/alpha*sqrt(gain)*(sqrt(gain*ron*ron+value)-sqrt(gain)*ron);
      data.image[IDout].array.F[ii] = value+0.5;
    }

  return(IDout);
}

long IMAGE_FORMAT_dequantize(char *IDin_name, char *IDout_name, double alpha, double ron, double gain, double bias)
{
  long IDin, IDout;
  long ii;
  long xsize,ysize;
  double value;

  IDin = image_ID(IDin_name);
  xsize = data.image[IDin].md[0].size[0];
  ysize = data.image[IDin].md[0].size[1];
  
  IDout = create_2Dimage_ID(IDout_name,xsize,ysize);
  for(ii=0;ii<xsize*ysize;ii++)
    {
      value = data.image[IDin].array.F[ii];
      if(value < 0.0)
	value = value*alpha*ron + bias;
      else
	{
	  value = alpha/2.0*value/sqrt(gain)+ron*sqrt(gain);
	  value = value*value;
	  value = value - gain*ron*ron + bias;
	}
      data.image[IDout].array.F[ii] = value;
    }

  return(IDout);
}

long IMAGE_FORMAT_read_binary16(char *fname, long xsize, long ysize, char *IDname)
{
  FILE *fp;
  char *buffer;
  unsigned long fileLen;
  long i, ii, jj;
  long ID;
  long v1;
  int r;
 
  //Open file
  fp = fopen(fname, "rb");
  if (!fp)
    {
      fprintf(stderr, "Unable to open file %s", fname);
      return(0);
    }
	
  //Get file length
  fseek(fp, 0, SEEK_END);
  fileLen=ftell(fp);
  fseek(fp, 0, SEEK_SET);
  
  //Allocate memory
  buffer=(char *)malloc(fileLen+1);
  if (!buffer)
    {
      fprintf(stderr, "Memory error!");
      fclose(fp);
      return(0);
    }

  //Read file contents into buffer
  r = fread(buffer, fileLen, 1, fp);
  fclose(fp);

  ID = create_2Dimage_ID(IDname, xsize, ysize);

  i = 0;
  for(jj=0;jj<ysize;jj++)
    for(ii=0;ii<xsize;ii++)
      {
	
	if(i<fileLen+1)
	  v1 = (long) (((unsigned char *)buffer)[i]) +  (long) (256*((unsigned char *)buffer)[i+1]);
	data.image[ID].array.F[jj*xsize+ii] = (float) v1;
	i += 2;
      }
  
  free(buffer);
  
  return(0);
}


long IMAGE_FORMAT_read_binary32f(char *fname, long xsize, long ysize, char *IDname)
{
  FILE *fp;
  float *buffer;
  unsigned long fileLen;
  long i, ii, jj;
  long ID;
  long v1;
  int r;

  //Open file
  fp = fopen(fname, "rb");
  if (!fp)
    {
      fprintf(stderr, "Unable to open file %s", fname);
      return(0);
    }
	
  //Get file length
  fseek(fp, 0, SEEK_END);
  fileLen=ftell(fp);
  fseek(fp, 0, SEEK_SET);
  
  //Allocate memory
  buffer = (float *)malloc(fileLen+1);
  if (!buffer)
    {
      fprintf(stderr, "Memory error!");
      fclose(fp);
      return(0);
    }

  //Read file contents into buffer
  r = fread(buffer, fileLen, 1, fp);
  fclose(fp);

  ID = create_2Dimage_ID(IDname, xsize, ysize);

  i = 0;
   for(jj=0;jj<ysize;jj++)
     for(ii=0;ii<xsize;ii++)
       {
	data.image[ID].array.F[jj*xsize+ii] = buffer[i];
	i++;
       }
  
  free(buffer);


  return(0);
}





long IMAGE_FORMAT_FITS_to_ushortintbin_lock( char *IDname, char *fname )
{
  long ID;
  long xsize, ysize;
  long ii;
  int fd;
  unsigned short int *valarray;
  int r;

  ID = image_ID(IDname);
  xsize = data.image[ID].md[0].size[0];
  ysize = data.image[ID].md[0].size[1];

  valarray = (unsigned short int*) malloc(sizeof(unsigned short int)*xsize*ysize);

  if(data.image[ID].md[0].atype == FLOAT)
    {
      printf("float -> unsigned short int array\n");
      for(ii=0;ii<xsize*ysize;ii++)
	valarray[ii] = (unsigned short int) data.image[ID].array.F[ii];
    }
  if(data.image[ID].md[0].atype == DOUBLE)
    {
      printf("double -> unsigned short int array\n");
      for(ii=0;ii<xsize*ysize;ii++)
	valarray[ii] = (unsigned short int) data.image[ID].array.D[ii];
    }

  fd = open(fname, O_RDWR | O_CREAT, S_IRUSR|S_IWUSR);
  flock(fd, LOCK_EX);
  if( fd < 0 )
    {
      perror( "Error opening file" );
      printf( "Error opening file \"%s\": %s\n", fname, strerror( errno ) );
    }
  r = write(fd, valarray, sizeof(unsigned short int)*xsize*ysize);     
  flock(fd, LOCK_UN);
  close(fd);


  free(valarray);

  return(0);
} 

long IMAGE_FORMAT_FITS_to_floatbin_lock(  char *IDname, char *fname )
{
  long ID;
  long xsize, ysize;
  long ii;
  int fd;
  float *valarray;
  int r;

  ID = image_ID(IDname);
  xsize = data.image[ID].md[0].size[0];
  ysize = data.image[ID].md[0].size[1];

  valarray = (float*) malloc(sizeof(float)*xsize*ysize);

  if(data.image[ID].md[0].atype == FLOAT)
    {
      printf("WRITING float array\n");
      for(ii=0;ii<xsize*ysize;ii++)
	valarray[ii] = data.image[ID].array.F[ii];
    }
  if(data.image[ID].md[0].atype == DOUBLE)
    {
      printf("WRITING double array\n");
      for(ii=0;ii<xsize*ysize;ii++)
	valarray[ii] = (float) data.image[ID].array.D[ii];
    }

  fd = open(fname, O_RDWR | O_CREAT, S_IRUSR|S_IWUSR);
  flock(fd, LOCK_EX);
  if( fd < 0 )
    printf( "Error opening file: %s\n", strerror( errno ) );
    
  r = write(fd, valarray, sizeof(float)*xsize*ysize);     
  //  for(ii=0;ii<xsize*ysize;ii++)
  //  printf("[%ld %f] ", ii, valarray[ii]);

  flock(fd, LOCK_UN);
  close(fd);


  free(valarray);

  return(0);
} 
