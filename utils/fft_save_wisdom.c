#include <stdio.h>
#include "define.inc"
#if FFT == _FFTW3_
#include <fftw3.h>
#elif FFT == _FFTW_
extern void fftw_export_wisdom_to_file(FILE *);
extern void fftw_import_wisdom_from_file(FILE *);
#endif

void save_wisdom_to_filename(char * filename){
#ifdef FFT
  FILE * file;
  file = fopen(filename, "w");
  fftw_export_wisdom_to_file(file);
  fclose(file);
#endif
  

}
void read_wisdom_from_filename(char * filename){
#ifdef FFT
  FILE * file;
	/*printf("HERE\n");*/
  file = fopen(filename, "r");
	/* Do nothing if the file doesn't exist*/
	if (file == NULL) return;
	/*printf("HERE2\n");*/
  fftw_import_wisdom_from_file(file);
	/*printf("HERE3\n");*/
  fclose(file);
	/*printf("HERE4\n");*/
#endif
  

}
