/*=================================================================
 *
 * PERCWALK.C	Sample .MEX file corresponding to PERCWALK.M
 *	        Return random walk on percolation cluster
 *
 * The calling syntax is:
 *
 *		[yp,nstep] = percwalk(y,nw,nsaw)
 *
 *  Saw mode is currently not implemented
 *
 *  You may also want to look at the corresponding M-code, percwalk.m.
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 2004 Anders Malthe-Sorenssen
 * Physics of Geological Processes, University of Oslo, Norway
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	Y_IN		prhs[0]
#define	NW_IN		prhs[1]
#define	NSAW_IN		prhs[2]

/* Output Arguments */

#define	YP_OUT		plhs[0]
#define	NSTEP_OUT	plhs[1]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define PI 3.14159265

static void percwalk(
		     double	yp[],
		     double	y[],
		     double	nw[],
		     double	*nstep,
		     double	*nsaw,
		     unsigned int m,
		     unsigned int n,
		     unsigned int nwalk
		     )
{
  double	r1,r2,nn,c,pos,sstep;
  int nchange,niter,neighb,idir,nnr;
  int ntot,i,ii,iix,iiy,ix,iy,nnsaw,npos,step,nstop,nwalk1,nc;
  int dir[8];
  int neighblist[8];
  double yy;
  int ixmin,iymin,rmin,xc,yc,idx,idy,idd;
  
  dir[0] =  1;
  dir[1] =  0;
  dir[2] = -1;
  dir[3] =  0;
  dir[4] =  0;
  dir[5] =  1;
  dir[6] =  0;
  dir[7] = -1;

  ntot = m*n;
  nn = *nsaw;
  nnsaw = (int)nn;

/*  printf("m,n = %i,%i\n",m,n);
  printf("nsaw = %i\n",nnsaw);
*/	
  /* First, find random starting position */
  xc = m/2;
  yc = n/2;
  rmin = m*n;
  ixmin = 0;
  iymin = 0;
  for (ix=0;ix<m;ix++) {
      for (iy=0;iy<n;iy++) {
          i = iy*m + ix;
          c = y[i];
          if (c>0.1) {
              idx = ix - xc;
              idy = iy - yc;
              idd = idx*idx + idy*idy;
              if (idd<rmin) {
                  rmin = idd;
                  ixmin = ix;
                  iymin = iy;
              }
          }
      }
  }
  ix = ixmin;
  iy = iymin;
  
  /*c = nw[0]; Get first random number from string, use this for initial position 
  pos = ntot*c;
  npos = floor(pos);
  iy = npos/m;
  ix = npos - m*iy; */

/*  printf("pos,npos = %lf,%i\n",pos,npos);
  printf("ix,iy = %i,%i\n",ix,iy);
*/
  yp[0] = ix+1;
  yp[1] = iy+1;
  nstop = 0;
  step = 0;
  nwalk1 = nwalk - 1;
  i = iy*m + ix;
  c = y[i];
/*  printf("z(ix,iy) = %f\n",c);*/
  if (c>0.1) {
    nstop = 0;
/*    printf("c0 = %lf\n",c);*/
  }
  else {
/*    printf("c1 = %lf\n",c);*/
    nstop = 1;
  }
/*  printf("nstop , step, nwalk1 = %i,%i,%i\n",nstop,step,nwalk1);*/
  while ((nstop==0)&&(step<nwalk1)) {
    /* Check if there is a way out */
    neighb = 0;
    for (idir=0;idir<4;idir++) {
      iix = ix + dir[idir*2];
      iiy = iy + dir[idir*2+1];
      /* Check if inside */
      if ((iix>=0)&&(iix<m)&&(iiy>=0)&&(iiy<n)) {
      	ii = iiy*m + iix;
      	if (y[ii]>0.1) {
      		neighblist[neighb]	 = idir;
      		neighb++;
      	}
      }		
    }		
    /*      printf("neighb = %i\n",neighb);*/
    if (neighb>0) { /* Possible way out of this site */
      step++;
      c = nw[step];
      c = c*neighb;
      nc = floor(c);
      if (nc>3) nc = 3;
      if (nc<0) nc = 0;
      /*	printf("nc = %i\n",nc);*/
      idir = neighblist[nc];
      ix = ix + dir[idir*2];
      iy = iy + dir[idir*2+1];
      /*	printf("ix,iy = %i,%i\n",ix,iy);
			printf("step = %i\n",step);*/
      yp[2*step]   = ix+1;
      yp[2*step+1] = iy+1;
    }
    else { /* No way out - stop walker */
      nstop = 1;
    }
  }
	
  sstep = step + 1.0;
    
  *nstep = sstep;
/*  printf("Step = %i\n",step);*/
  return;
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
  double *yp; 
  double *nsaw,*y,*nw,*nstep;
  unsigned int m,n,nwalk,nwalk2;
    
  /* Check for proper number of arguments */
    
  if (nrhs != 3) { 
    mexErrMsgTxt("Three input arguments required."); 
  } else if (nlhs > 2) {
    mexErrMsgTxt("Too many output arguments."); 
  } 
    
  /* Find dimensions of Y. */
    
  m = mxGetM(Y_IN); 
  n = mxGetN(Y_IN);

  /* Find dimensions of NW. */
  nwalk = mxGetM(NW_IN); 
    
  /* Create a matrix for the return argument */ 
  nwalk2 = 2;
  YP_OUT = mxCreateDoubleMatrix(nwalk2, nwalk, mxREAL); 
  NSTEP_OUT = mxCreateDoubleMatrix(1, 1, mxREAL); 
    
  /* Assign pointers to the various parameters */ 
  yp = mxGetPr(YP_OUT);
  nstep = mxGetPr(NSTEP_OUT);
    
  y = mxGetPr(Y_IN);
  nw = mxGetPr(NW_IN); 
  nsaw = mxGetPr(NSAW_IN); 
 
  /* Do the actual computations in a subroutine */
  percwalk(yp,y,nw,nstep,nsaw,m,n,nwalk); 
  return;

}
