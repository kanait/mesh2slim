#ifndef JACOBI_H 
#define JACOBI_H 

#include <math.h>
#include <stdlib.h>
#define NRANSI
#include "nrutil.h"
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

#define NR_END 1
#define FREE_ARG char*

class Jacobi{
public:

static void jacobi(double **a, int n, double d[], double **v, int *nrot)
{
  int j,iq,ip,i;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
  
  b=vector(1,n);
  z=vector(1,n);
  for (ip=1;ip<=n;ip++) {
    for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=1;ip<=n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++)
        sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      free_vector(z,1,n);
      free_vector(b,1,n);
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++) {
        g=100.0*fabs(a[ip][iq]);
        if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
            && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
          a[ip][iq]=0.0;
        else if (fabs(a[ip][iq]) > tresh) {
          h=d[iq]-d[ip];
          if ((double)(fabs(h)+g) == (double)fabs(h))
            t=(a[ip][iq])/h;
          else {
            theta=0.5*h/(a[ip][iq]);
            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/sqrt(1+t*t);
          s=t*c;
          tau=s/(1.0+c);
          h=t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq]=0.0;
          for (j=1;j<=ip-1;j++) {
            ROTATE(a,j,ip,j,iq)
            }
          for (j=ip+1;j<=iq-1;j++) {
            ROTATE(a,ip,j,j,iq)
            }
          for (j=iq+1;j<=n;j++) {
            ROTATE(a,ip,j,iq,j)
            }
          for (j=1;j<=n;j++) {
            ROTATE(v,j,ip,j,iq)
            }
          ++(*nrot);
        }
      }
    }
    for (ip=1;ip<=n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  nrerror("Too many iterations in routine jacobi");
}

static void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  //exit(1);
}

static double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;
  
  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}


static void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}



};
#undef ROTATE
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 9z!+!1(t+%. */

#endif