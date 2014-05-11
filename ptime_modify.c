// To do:
// 1. Must deallocate memory
// 2. Not getting identical results with paas (# of decimal places? Must check this)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ptimeLib.h"
#include <cpgplot.h>
#include "readPfits.h"
#include "fitsio.h"

// Program to modify templates

void writePaasModelFile(tmplStruct *tmpl,char *file);
void updatePfits(tmplStruct *tmpl,char *file,char *fname_out);
void findMinMax(float *y,int n,float *min,float *max);

int main(int argc,char *argv[])
{
  tmplStruct tmpl;
  int i;
  int nbin = 1024;
  char fname[128];
  char paasModelFile[128];
  char pfitsFile[128];
  char outputFits[128];
  int paasModel=0;
  int pfitsLoad = 1;

  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
      else if (strcmp(argv[i],"-pf")==0)
	{
	  strcpy(pfitsFile,argv[++i]);
	  strcpy(outputFits,argv[++i]);
	  pfitsLoad=1;
	}
      else if (strcmp(argv[i],"-paasModel")==0)
	{
	  paasModel=1;
	  strcpy(paasModelFile,argv[++i]);
	}
    }
  
  initialiseTemplate(&tmpl);
  printf("Reading template\n");
  readTemplate(fname,&tmpl); 
  printf("Complete reading template\n");
  
  if (paasModel==1){
    writePaasModelFile(&tmpl,paasModelFile);
  }
  if (pfitsLoad==1){
    updatePfits(&tmpl,pfitsFile,outputFits);
  }
  
  // Must deallocate the template memory
}

void updatePfits(tmplStruct *tmpl,char *fname,char *fname_out)
{
  pheader *header;
  fitsfile *fp_in;
  fitsfile *fp_out;
  int status=0;
  int colnum=0;
  int i,j,k;
  int nbin = 2048; // MUST READ FROM FILE
  int nchan = tmpl->nchan;
  int npol = tmpl->channel[0].nstokes;
  float dat_scl[nchan*npol];
  float dat_offs[nchan*npol];
  float fy[nbin*nchan];
  float ty[nbin];
  char fname_out_mod[1024];
  int ii=1;
  float min,max,scale,offs;

  sprintf(fname_out_mod,"!%s",fname_out);

  header = (pheader *)malloc(sizeof(pheader));

  fp_in = openFitsFile(fname);

  // Copy the input file to the output file
  /* Create the output file */
  printf("Creating file\n");
  if ( !fits_create_file(&fp_out, fname_out_mod, &status) )
    {
      /* Copy every HDU until we get an error */
      while( !fits_movabs_hdu(fp_in, ii++, NULL, &status) )
	fits_copy_hdu(fp_in, fp_out, 0, &status);
 
      /* Reset status after normal error */
      if (status == END_OF_FILE) status = 0;
    }
  fits_close_file(fp_in, &status);
  fits_close_file(fp_out, &status);
  printf("Finished creating file\n");

  fp_in = openFitsFile_readWrite(fname_out);
  loadPrimaryHeader(fp_out,header);
  // Update the data values in this pfits file

  fits_movnam_hdu(fp_out,BINARY_TBL,(char *)"SUBINT",0,&status);
  if (status) {fits_report_error(stdout,status); exit(1);}
  // Should do checks to ensure that the template will fit in the profile
  fits_get_colnum(fp_out, CASEINSEN, (char *)"DATA", &colnum, &status);  
  if (status) {fits_report_error(stdout,status); exit(1);}

  // Must set sensible scaling and dat_offs
  for (i=0;i<tmpl->channel[0].nstokes;i++){
    for (j=0;j<tmpl->nchan;j++){
      for (k=0;k<nbin;k++){
	ty[k] = (float)evaluateTemplateChannel(tmpl,(double)(k+0.5)/(double)nbin,j,i,0);
      }
      findMinMax(ty,nbin,&min,&max);
      scale = (max-min)/2.0/16384.0;
      offs = max - 16384.0*scale;
      for (k=0;k<nbin;k++){
	fy[j*nbin+k] = (ty[k]-offs)/scale;
      }
      dat_scl[j+i*nchan] = scale;
      dat_offs[j+i*nchan] = offs;
    }
    fits_write_col(fp_out,TFLOAT,colnum,1,1,nbin*nchan,fy,&status);
    if (status) {fits_report_error(stdout,status); exit(1);}
  }
  fits_get_colnum(fp_out, CASEINSEN, (char *)"DAT_OFFS", &colnum, &status);
  fits_write_col(fp_out,TFLOAT,colnum,1,1,nchan*npol,dat_offs,&status);
  fits_get_colnum(fp_out, CASEINSEN, (char *)"DAT_SCL", &colnum, &status);
  fits_write_col(fp_out,TFLOAT,colnum,1,1,nchan*npol,dat_scl,&status);
  closeFitsFile(fp_out);
  free(header);
}
void findMinMax(float *y,int n,float *min,float *max)
{
  int i;
  *min = *max = y[0];
  for (i=0;i<n;i++)
    {
      if (y[i] > *max) *max = y[i];
      if (y[i] < *min) *min = y[i];
    }

}

void writePaasModelFile(tmplStruct *tmpl,char *file)
{
  FILE *fout;

  if (!(fout = fopen(file,"w")))
    {
      printf("Unable to open output file >%s<\n",file);
    }
  else
    {
      int i;
      for (i=0;i<tmpl->channel[0].pol[0].nComp;i++)
	fprintf(fout,"%.5f %.5f %.5f\n",
		tmpl->channel[0].pol[0].comp[i].centroid,
		tmpl->channel[0].pol[0].comp[i].concentration,
		tmpl->channel[0].pol[0].comp[i].height);
      printf("Completed writing to: %s\n",file);
      fclose(fout);
    }
}
