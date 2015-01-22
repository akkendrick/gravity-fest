/* findfaultsM.c  -  find coordinates of fault nodes for interpolation */
/* usage:   findfaultsM coordfile fltfile outfile  */
/* format of coord file:
 numnp
 node#     xcoord  ycoord  zcoord
 ...
  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX(x,y)   (((x) < (y)) ? (y) : (x))
#define MIN(x,y)   (((x) < (y)) ? (x) : (y))

int main(
     int argc ,      /* number of input arguments */
     char *argv[]    /* input argument strings */
    )
{
 int  i , j , k , l , n , numnp , numel , nfnode , master , ndum1 , ndum2 ;
 int  i1 , i2 , i3 , n1 , n2, n3 , n4 , n5 , n6 , n7 , n8 , side , upper ;
 int  nn1 , nn2 , nn3 , nn4 , newel , f1 , f2 , f3 , mat , countup , countlow ;

 double  *x , *y , *z , *r , xx , yy , zz , rlevel , rcentroid , rtol ;
 double  bx,by,bz,sx,sy,sz,slip;
 int  *ien  ;

 FILE    *coordfile , *fltfile , *outfile ;
 
 coordfile = fopen(argv[1],"r") ;
 fltfile = fopen(argv[2],"r") ;
 outfile = fopen(argv[3],"w") ;

 fscanf(coordfile,"%d",&numnp) ;
 x = (double *)calloc(numnp, sizeof(double)) ;
 y = (double *)calloc(numnp, sizeof(double)) ;
 z = (double *)calloc(numnp, sizeof(double)) ;
 for(i=0 ; i<numnp ; i++)
    {
     fscanf(coordfile,"%d%lf%lf%lf",&n , &xx , &yy , &zz) ;
     x[n-1]=xx/1000.0 ;
     y[n-1]=yy/1000.0 ;
     z[n-1]=zz/1000.0 ;
    }

 fscanf(fltfile,"%d",&nfnode) ;

 for(i=0 ; i<nfnode ; i++)
    {
     fscanf(fltfile,"%d%d%lf%lf%lf%lf%lf%lf%d%lf",&n , &ndum1 , &bx , &by , &bz ,
            &sx , &sy , &sz , &ndum2 , &slip) ;
     xx = x[n-1] ;
     yy = y[n-1] ;
     zz = z[n-1] ;
     fprintf(outfile,"%d   %g   %g   %g\n", n , xx , yy , zz) ;
    }

  
 fclose(coordfile) ;
 fclose(fltfile) ;
 fclose(outfile) ;
}