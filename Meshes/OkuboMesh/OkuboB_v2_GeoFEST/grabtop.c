/* grabtop.c  -  get top node displacement */
/* usage:   grabtop GFoutfile nodelist selectout */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(
     int argc ,      /* number of input arguments */
     char *argv[]    /* input argument strings */
    )
{
 int  numnp,numel,i,node,nel,ndum,n1,n2,n3,n4,offset,mat ;
 double  xx,yy,zz,sx,sy,sz,vx,vy,vz,dt,st00,st11,st22,st01,st10,st02,st20,st12,st21 ;
 double  *x , *y , *z , *dx , *dy , *dz ;
 char  buf[80] ;
 FILE    *infile , *listfile , *outfile ;
 
 if(argc != 4)
    {
     printf("Wrong number of arguments.\n") ;
     exit(0) ;
    }

 infile = fopen(argv[1],"r") ;
 listfile = fopen(argv[2],"r") ;
 outfile = fopen(argv[3],"w") ;
 
/* ************************************************** */
fscanf(infile, "%d" , &numnp) ;
x = (double *) calloc(numnp,sizeof(double)) ;
y = (double *) calloc(numnp,sizeof(double)) ;
z = (double *) calloc(numnp,sizeof(double)) ;
dx = (double *) calloc(numnp,sizeof(double)) ;
dy = (double *) calloc(numnp,sizeof(double)) ;
dz = (double *) calloc(numnp,sizeof(double)) ;

for(i=0;i<numnp;i++)
   {
    fscanf(infile, "%s%d%lf%lf%lf%lf%lf%lf%lf%lf%lf" , buf ,&node , &xx , &yy , &zz , &sx , &sy , &sz , &vx , &vy , &vz) ;
    *(x + node -1) = xx ;
    *(y + node -1) = yy ;
    *(z + node -1) = zz ;
    *(dx + node -1) = sx ;
    *(dy + node -1) = sy ;
    *(dz + node -1) = sz ;
   }

fscanf(listfile, "%d" , &numnp) ;
for(i=0;i<numnp;i++)
   {
    fscanf(listfile, "%d" ,&node ) ;
    xx = *(x + node -1) ;
    yy = *(y + node -1) ;
    zz = *(z + node -1) ;
    sx = *(dx + node -1) ;
    sy = *(dy + node -1) ;
    sz = *(dz + node -1) ;
    fprintf(outfile,"%g   %g   %g      %g   %g   %g\n",xx,yy,zz,sx,sy,sz) ;
   }


fprintf(outfile,"\n") ;

 fclose(infile) ;
 fclose(listfile) ;
 fclose(outfile) ;
 printf("Processing complete.\n") ;

}
