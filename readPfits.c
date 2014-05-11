#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fitsio.h"
#include "readPfits.h"

// Routines for reading PSRFITS files

void closeFitsFile(fitsfile *fp)
{
  int status=0;
  fits_close_file(fp,&status);
  fits_report_error(stderr,status);
  if (status)
    {
      printf("Error closing file\n");
      exit(1);
    }

}

fitsfile * openFitsFile(char *fname)
{
  fitsfile *fp;
  int status=0;
  fits_open_file(&fp,fname,READONLY,&status);
  fits_report_error(stderr,status);
  if (status)
    {
      printf("Error opening file >%s<\n",fname);
      exit(1);
    }
  return fp;
}

fitsfile * openFitsFile_readWrite(char *fname)
{
  fitsfile *fp;
  int status=0;
  fits_open_file(&fp,fname,READWRITE,&status);
  fits_report_error(stderr,status);
  if (status)
    {
      printf("Error opening file >%s<\n",fname);
      exit(1);
    }
  return fp;
}

void loadPrimaryHeader(fitsfile *fp,pheader *phead)
{
  int status=0;
  int nkey=-1;
  int morekeys=-1;
  int i;
  char keyname[128],val[128],comment[128];

  fits_get_hdrspace(fp,&nkey,&morekeys,&status);

  phead->nhead = nkey;
  // Allocate memory
  phead->keyname = (char **)malloc(sizeof(char *)*nkey);
  phead->val = (char **)malloc(sizeof(char *)*nkey);
  phead->comment = (char **)malloc(sizeof(char *)*nkey);
  for (i=0;i<nkey;i++)
    {
      phead->keyname[i] = (char *)malloc(sizeof(char)*128);
      phead->val[i] = (char *)malloc(sizeof(char)*128);
      phead->comment[i] = (char *)malloc(sizeof(char)*128);
    }

  // Complete allocating memory

  for (i=1;i<=nkey;i++)
    {
      fits_read_keyn(fp,i+1,phead->keyname[i-1],phead->val[i-1],phead->comment[i-1],&status);
      if (strcmp(phead->keyname[i-1],"OBSFREQ")==0)
	sscanf(phead->val[i-1],"%f",&(phead->freq));
      else if (strcmp(phead->keyname[i-1],"STT_IMJD")==0)
	sscanf(phead->val[i-1],"%d",&(phead->imjd));
      else if (strcmp(phead->keyname[i-1],"STT_SMJD")==0)
	sscanf(phead->val[i-1],"%f",&(phead->smjd));
      else if (strcmp(phead->keyname[i-1],"STT_OFFS")==0)
	sscanf(phead->val[i-1],"%f",&(phead->stt_offs));
      else if (strcmp(phead->keyname[i-1],"OBSBW")==0)
	sscanf(phead->val[i-1],"%f",&(phead->bw));
    }
  // Read specific parameters
  fits_read_key(fp,TSTRING,(char *)"OBS_MODE",phead->obsMode,NULL,&status);
  fits_read_key(fp,TSTRING,(char *)"SRC_NAME",phead->source,NULL,&status);
  fits_read_key(fp,TSTRING,(char *)"TELESCOP",phead->telescope,NULL,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }

  // Now load information from the subintegration table
  fits_movnam_hdu(fp,BINARY_TBL,(char *)"SUBINT",1,&status);
  if (status)
    {
      printf("No subintegration table\n");
      status=0;
    }
  else
    {
      fits_read_key(fp,TINT,(char *)"NAXIS2",&(phead->nsub),NULL,&status);
      if (status)
	{
	  printf("Reading naxis2\n");
	  fits_report_error(stderr,status);
	  exit(1);
	}     
      fits_read_key(fp,TINT,(char *)"NCHAN",&(phead->nchan),NULL,&status);
      if (status)
	{
	  printf("Reading nchan\n");
	  fits_report_error(stderr,status);
	  exit(1);
	}
      
      /*      fits_read_key(fp,TFLOAT,(char *)"ZERO_OFF",&(phead->zeroOff),NULL,&status);
      if (status)
	{
	  printf("Reading zero_off\n");
	  fits_report_error(stderr,status);
	  phead->zeroOff = 0;
	  status=0;
	  printf("Complete reading zero_off\n");
	  }  */
      printf("COMMENTED OUT READING ZERO_OFF AND NBITS IN PFITS.C -- PUT BACK -- PROBLEM WITH SOME FILES\n"); 
      /*      fits_read_key(fp,TINT,(char *)"NBITS",&(phead->nbits),NULL,&status);
      if (status)
	{
	  printf("Reading nbits\n");
	  printf("Have error\n");
	  fits_report_error(stderr,status);
	  exit(1);
	  }*/
      fits_read_key(fp,TINT,(char *)"NPOL",&(phead->npol),NULL,&status);
      if (status)
	{
	  printf("Reading npol\n");
	  fits_report_error(stderr,status);
	  exit(1);
	}
      
      fits_read_key(fp,TINT,(char *)"NSBLK",&(phead->nsblk),NULL,&status);
      if (status)
	{
	  printf("Reading nsblk\n");
	  fits_report_error(stderr,status);
	  exit(1);
	}

      fits_read_key(fp,TINT,(char *)"NBIN",&(phead->nbin),NULL,&status);
      if (status)
	{
	  printf("Reading nbin\n");
	  fits_report_error(stderr,status);
	  exit(1);
	}

      //      printf("nbin = %d (%d)\n",phead->nbin,status);
      fits_read_key(fp,TFLOAT,(char *)"CHAN_BW",&(phead->chanbw),NULL,&status);
      if (phead->chanbw < 0 && phead->bw > 0)
	phead->bw*=-1;
      
      fits_read_key(fp,TFLOAT,(char *)"TBIN",&(phead->tsamp),NULL,&status);
      
    }
  fits_movnam_hdu(fp,BINARY_TBL,(char *)"PSRPARAM",1,&status);
  if (status)
    {
      printf("No PSRPARM table\n");
      status=0;
    }
  else
    {
      int len,i,colnum;
      char **line,str1[1024],str2[1024];
      char nval[128]="UNKNOWN";
      int anynul=0;
      float tt;
      fits_read_key(fp,TINT,(char *)"NAXIS2",&len,NULL,&status);

      fits_get_colnum(fp,CASEINSEN,(char *)"PARAM",&colnum,&status);
      if (status) {
	printf("Unable to find data in the psrparam table in FITS file\n");
	exit(1);
      }

      line = (char **)malloc(sizeof(char *));
      line[0] = (char *)malloc(sizeof(char)*1024); 

      for (i=0;i<len;i++)
	{
	  fits_read_col_str(fp,colnum,i+1,1,1,nval,line,&anynul,&status);
	  if (sscanf(line[0],"%s %s",str1,str2)==2)
	    {
	      if (strcasecmp(str1,"DM")==0)
		sscanf(str2,"%f",&(phead->dm));
	      if (strcasecmp(str1,"F0")==0)
		{
		  sscanf(str2,"%f",&tt);
		  phead->period = 1.0/tt;
		}
	    }
	  //	  printf("Read: %s\n",line[0]);
	}
      //      printf("Lenght = %d\n",len);
  free(line[0]);
  free(line);

    }
    
}

void allocateObsMemory(ptime_observation *obs,pheader *phead)
{
  int nbin,nchan,npol;
  int i,j;

  nchan = phead->nchan;
  npol  = phead->npol;
  nbin  = phead->nbin;

  obs->chan = (ptime_chan *)malloc(sizeof(ptime_chan)*nchan);
  obs->nchan = nchan;
  for (i=0;i<nchan;i++)
    {
      obs->chan[i].pol = (ptime_pol *)malloc(sizeof(ptime_pol)*npol);
      obs->chan[i].npol = npol;
      for (j=0;j<npol;j++)
	{
	  obs->chan[i].pol[j].val = (float *)malloc(sizeof(float)*nbin);
	  obs->chan[i].pol[j].nbin = nbin;
	}
    }
}

void deallocateMemory(ptime_observation *obs)
{
  int i,j;
  for (i=0;i<obs->nchan;i++)
    {
      for (j=0;j<obs->chan[i].npol;j++)
	{
	  
	  free(obs->chan[i].pol[j].val);
	}
      free(obs->chan[i].pol);
    }
  free(obs->chan);
  free(obs);
}

void readDatFreq(ptime_observation *obs,double *datFreq,fitsfile *fp,int nchan)
{
  int status=0;
  int colnum;
  int nval =0;
  int initflag=0;

  fits_movnam_hdu(fp,BINARY_TBL,(char *)"SUBINT",1,&status);
  if (status) {
    printf("Unable to move to subint table in FITS file\n");
    exit(1);
  }
  fits_get_colnum(fp,CASEINSEN,(char *)"DAT_FREQ",&colnum,&status);
  if (status) {
    printf("Unable to find DAT_FREQ in the subint table in FITS file\n");
    exit(1);
  }
  fits_read_col_dbl(fp,colnum,1,1,nchan,nval,datFreq,&initflag,&status);

}
void readSubintOffs(ptime_observation *obs,double *offs_sub,fitsfile *fp)
{
  int status=0;
  int nsub = 1;
  int colnum;
  int nval =0;
  int initflag=0;

  fits_movnam_hdu(fp,BINARY_TBL,(char *)"SUBINT",1,&status);
  if (status) {
    printf("Unable to move to subint table in FITS file\n");
    exit(1);
  }
  fits_get_colnum(fp,CASEINSEN,(char *)"OFFS_SUB",&colnum,&status);
  if (status) {
    printf("Unable to find OFFS_SUB in the subint table in FITS file\n");
    exit(1);
  }
  fits_read_col_dbl(fp,colnum,1,1,nsub,nval,offs_sub,&initflag,&status);

}

void readData(ptime_observation *obs,pheader *phead,fitsfile *fp)
{
  int status=0;
  int i,j,k,l;
  int initflag=0;
  int nval=0;
  int colnum;
  int nchan,nbin,npol;
  float ty[phead->nbin];
  float datScl[phead->nchan*phead->npol];
  float datOffs[phead->nchan*phead->npol];
  nchan = phead->nchan;
  nbin = phead->nbin;
  npol = phead->npol;

  //
  fits_movnam_hdu(fp,BINARY_TBL,(char *)"SUBINT",1,&status);
  if (status) {
    printf("Unable to move to subint table in FITS file\n");
    exit(1);
  }
  fits_get_colnum(fp,CASEINSEN,(char *)"DAT_SCL",&colnum,&status);
  if (status) {
    printf("Unable to find DAT_SCL in the subint table in FITS file\n");
    exit(1);
  }
  fits_read_col_flt(fp,colnum,1,1,phead->nchan*phead->npol,nval,datScl,&initflag,&status);


  fits_get_colnum(fp,CASEINSEN,(char *)"DAT_OFFS",&colnum,&status);
  if (status) {
    printf("Unable to find DAT_OFFS in the subint table in FITS file\n");
    exit(1);
  }
  fits_read_col_flt(fp,colnum,1,1,phead->nchan*phead->npol,nval,datOffs,&initflag,&status);

  fits_get_colnum(fp,CASEINSEN,(char *)"DATA",&colnum,&status);  
  if (status) {
    printf("Unable to find data in the subint table in FITS file\n");
    exit(1);
  }

  for (j=0;j<npol;j++)
    {
      for (i=0;i<nchan;i++)
	{
	  fits_read_col_flt(fp,colnum,1,j*(nchan*nbin)+i*nbin+1,nbin,nval,ty,&initflag,&status);
	  for (k=0;k<nbin;k++)
	    {
	      obs->chan[i].pol[j].val[k] = (ty[k]+datOffs[j*nchan+i])*datScl[j*nchan+i];
	    }
	}
    }
}

void removeBaseline(ptime_observation *obs,pheader *phead,int baselineType,float baselineFrac)
{
  int nchan,npol,nbin;
  int i,j,k,k0,k1,it;
  double bl,bl_best,bl2,bl2_best;
  int setbl=0;
  int best_k0=0;
  int best_k1=0;
  int nc=0;
  int best_nc=0;

  nchan = phead->nchan;
  nbin = phead->nbin;
  npol = phead->npol;

  //  for (i=0;i<nbin;i++)
  //    printf("blcalc: %d %g\n",i,obs->chan[0].pol[0].val[i]);


  for (i=0;i<nchan;i++)
    {
      for (j=0;j<npol;j++)
	{
	  if (j==0) {
	    setbl=0;
	    
	    for (k0=0;k0<nbin;k0++)
	      {
		for (it = 0; it < baselineFrac*nbin; it++)
		  {
		    bl=0;
		    bl2=0;
		    nc=0;
		    k1 = k0+it;
		    if (k1 > nbin) k1 = k1-nbin;
		    
		    if (k1 > k0) {
		      for (k=k0;k<k1;k++)
			{
			  bl += obs->chan[i].pol[j].val[k];
			  bl2 += pow(obs->chan[i].pol[j].val[k],2);
			  nc++;
			}		    
		    }
		    else if (k1 < k0) {
		      for (k=k0;k<nbin;k++)
			{
			  bl += obs->chan[i].pol[j].val[k];
			  bl2 += pow(obs->chan[i].pol[j].val[k],2);
			  nc++;
			}
		      for (k=0;k<k1;k++)
			{
			  bl += obs->chan[i].pol[j].val[k];
			  bl2 += pow(obs->chan[i].pol[j].val[k],2);
			  nc++;
			}
		    }
		  }
		if (setbl==0) {bl_best = bl; setbl=1; bl2_best = bl2; best_nc=nc;}
		else {
		  if (bl_best > bl) {bl_best = bl; best_k0=k0; best_k1=k1; bl2_best = bl2; best_nc=nc;}
		}
	      }
	  } else {
	    bl = 0;
	    bl2 = 0;
	    for (it = 0; it < baselineFrac*nbin; it++)
	      {
		bl=0;
		bl2=0;
		nc=0;
		k1 = best_k0+it;
		if (k1 > nbin) k1 = k1-nbin;
		
		if (k1 > best_k0) {
		  for (k=best_k0;k<k1;k++)
		    {
		      bl += obs->chan[i].pol[j].val[k];
		      bl2 += pow(obs->chan[i].pol[j].val[k],2);
		      nc++;
		    }
		}
		else if (k1 < best_k0) {
		  for (k=k0;k<nbin;k++)
		    {
		      bl += obs->chan[i].pol[j].val[k];
		      bl2 += pow(obs->chan[i].pol[j].val[k],2);
		      nc++;
		    }
		  for (k=0;k<k1;k++)
		    {
		      bl += obs->chan[i].pol[j].val[k];
		      bl2 += pow(obs->chan[i].pol[j].val[k],2);
		      nc++;
		    }
		}
	      }
	    bl_best = bl;
	    bl2_best= bl2;
	    best_nc = nc;
	  }


	  obs->chan[i].pol[j].baselineVal = bl_best/(double)best_nc;
	  obs->chan[i].pol[j].baseline_b0 = best_k0;
	  obs->chan[i].pol[j].baseline_b1 = best_k1;

	  {
	    // Calculate standard deviation
	    double sdev=0;
	    nc=0;
	    printf("Calculating sdev %d %d\n",best_k0,best_k1);
	    if (best_k1 > best_k0) {
	      for (k=best_k0;k<best_k1;k++)
		{
		  sdev += pow(obs->chan[i].pol[j].val[k]-bl_best/(double)best_nc,2);
		  nc++;
		}
	    }
	    else if (best_k1 < best_k0) {
	      for (k=best_k0;k<nbin;k++)
		{
		  sdev += pow(obs->chan[i].pol[j].val[k]-bl_best/(double)best_nc,2);
		  printf("Processing %g %g\n",obs->chan[i].pol[j].val[k],bl_best/(double)best_nc);
		  nc++;
		}
	      printf("Here %g %d\n",sdev,nc);
	      for (k=0;k<best_k1;k++)
		{
		  sdev += pow(obs->chan[i].pol[j].val[k]-bl_best/(double)best_nc,2);
		  nc++;
		}
	      printf("Here2 %g %d\n",sdev,nc);
	    }
	    
	    obs->chan[i].pol[j].sdev = sqrt(sdev/(double)nc);
	  }
	  printf("Using %g %g %g %g %d %d\n",bl2_best,bl_best,obs->chan[i].pol[j].baselineVal,obs->chan[i].pol[j].sdev,best_k0,best_k1);

	  for (k=0;k<nbin;k++)
	    (obs->chan[i].pol[j].val[k])-=(bl_best/(double)best_nc);

	}
    }
}
