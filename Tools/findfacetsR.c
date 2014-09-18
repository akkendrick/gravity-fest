/* findfacetsR.c  -  find [radial] surface or horizon facets from coords and ien */
/* usage:   findfacetsR infile outfile rvalue rtol */
/* rvalue is the floating point radius of the spherical surface to test for */
/* rtol is the allowable deviation from rvalue  */
/* format of input file:
 numnp
 node#     xcoord  ycoord  zcoord
 ...
 
 numel
 el#   0   mat#   n1  n2  n3  n4
 ...
 
 (this written for 3-d tet elements only)
 */

#include <stdio.h>
#include <math.h>
#define MAX(x,y)   (((x) < (y)) ? (y) : (x))
#define MIN(x,y)   (((x) < (y)) ? (x) : (y))

main(
     int argc ,      /* number of input arguments */
     char *argv[]    /* input argument strings */
    )
{
 int  i , j , k , l , n , numnp , numel , nfnode , master ;
 int  i1 , i2 , i3 , n1 , n2, n3 , n4 , n5 , n6 , n7 , n8 , side , upper ;
 int  nn1 , nn2 , nn3 , nn4 , newel , f1 , f2 , f3 , mat , countup , countlow ;

 double  *x , *y , *z , *r , xx , yy , zz , rlevel , rcentroid , rtol ;
 int  *ien  ;

 FILE    *infile , *outfile ;
 
 infile = fopen(argv[1],"r") ;
 outfile = fopen(argv[2],"w") ;
 sscanf(argv[3],"%lf",&rlevel) ;
 sscanf(argv[4],"%lf",&rtol) ;

 fscanf(infile,"%d",&numnp) ;
 x = (double *)calloc(numnp, sizeof(double)) ;
 y = (double *)calloc(numnp, sizeof(double)) ;
 z = (double *)calloc(numnp, sizeof(double)) ;
 r = (double *)calloc(numnp, sizeof(double)) ;
 for(i=0 ; i<numnp ; i++)
    {
     fscanf(infile,"%d%lf%lf%lf",&n , &xx , &yy , &zz) ;
     x[n-1]=xx ;
     y[n-1]=yy ;
     z[n-1]=zz ;
     r[n-1] = sqrt(xx*xx+yy*yy+zz*zz) ;
    }

 countup = 0 ;
 countlow = 0 ;
 fscanf(infile,"%d",&numel) ;
 /*
 ien = (int *)calloc(4*numel, sizeof(int)) ;
 */

 for(i=0 ; i<numel ; i++)
    {
     fscanf(infile,"%d%d%d%d%d%d%d",&n, &k, &mat, &n1 , &n2 , &n3 , &n4) ;
   /*
     *(ien+4*i) = n1 ;
     *(ien+4*i+1) = n2 ;
     *(ien+4*i+2) = n3 ;
     *(ien+4*i+3) = n4 ;
   */
     if(fabs(r[n1-1]-rlevel) < rtol && fabs(r[n2-1]-rlevel) < rtol &&
        fabs(r[n3-1]-rlevel) < rtol) side = 4 ;
     else if(fabs(r[n1-1]-rlevel) < rtol && fabs(r[n2-1]-rlevel) < rtol &&
        fabs(r[n4-1]-rlevel) < rtol) side = 3 ;
     else if(fabs(r[n1-1]-rlevel) < rtol && fabs(r[n3-1]-rlevel) < rtol &&
        fabs(r[n4-1]-rlevel) < rtol) side = 2 ;
     else if(fabs(r[n2-1]-rlevel) < rtol && fabs(r[n3-1]-rlevel) < rtol &&
        fabs(r[n4-1]-rlevel) < rtol) side = 1 ;
     else side = 0 ;
     rcentroid = (r[n1-1]+r[n2-1]+r[n3-1]+r[n4-1])/4.0 ;
     if(rcentroid >= rlevel) upper = 0 ;
     else  upper = 1 ;

     if(side && upper)
       {
        countup++ ;
        fprintf(outfile,"%d   %d      0.0    0.0   1.0(upper)\n", n , side) ;
       }
     else if(side && !upper)
       {
        countlow++ ;
        fprintf(outfile,"%d   %d      0.0    0.0   1.0(lower)\n", n , side) ;
       }
    }
  fprintf(outfile,"upper facets: %d\n",countup) ;
  fprintf(outfile,"lower facets: %d\n",countlow) ;
    
 
 fclose(infile) ;
 fclose(outfile) ;
}