/* gft2vtk.c  -  convert GeoFEST output to vtk legacy file */
/* usage:   gft2vtk coordfile ienfile dispfile strfile outfile.vtk */
/* version: Rupert Deese 2013 -- edits being made to restore program to
   functionality. */

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
 char  buf[80] ;
 FILE    *coordfile , *ienfile , *dispfile , *strfile , *outfile ;
 
 if(argc != 6)
    {
     printf("Wrong number of arguments.\n") ;
     printf("usage:   gft2vtk coordfile ienfile dispfile strfile outfile.vtu\n") ;
     exit(0) ;
    }

 coordfile = fopen(argv[1],"r") ;
 ienfile = fopen(argv[2],"r") ;
 dispfile = fopen(argv[3],"r") ;
 strfile = fopen(argv[4],"r") ;
 outfile = fopen(argv[5],"w") ;
 
/* ************************************************** 
    COORDFILE: Expected to have a first line with the number
    of nodes; all subsequent lines are node # and coords. */

fprintf(outfile,"# vtk DataFile Version 2.0\n") ;
fprintf(outfile,"FEM Grid Data\n") ;
fprintf(outfile,"ASCII\n") ;
fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n") ;
fscanf(coordfile, "%d" , &numnp) ;
fprintf(outfile,"POINTS %d float\n",numnp) ;
for(i=0;i<numnp;i++)
   {
    fscanf(coordfile, "%d%lf%lf%lf" , &node , &xx , &yy , &zz) ;
    fprintf(outfile,"%g   %g   %g\n",xx,yy,zz) ;
   }
fprintf(outfile,"\n") ;
printf("past coordfile");

/* ************************************************** 
    IENFILE: Expected to have a first line with the number of 
    elements; all subsequent lines are el #, dummy, mat. #,
    4 node connectivity. */

fscanf(ienfile, "%d" , &numel) ;
ndum=numel*5 ;
fprintf(outfile,"CELLS %d %d\n",numel,ndum) ;
for(i=0;i<numel;i++)
   {
    fscanf(ienfile, "%d%d%d%d%d%d%d" , &nel , &ndum , &mat , &n1 , &n2 , &n3 , &n4) ;
    n1-- ;
    n2-- ;
    n3-- ;
    n4-- ;
    fprintf(outfile,"4   %d   %d   %d   %d\n",n1,n2,n3,n4) ;
   }
fprintf(outfile,"\n") ;
fprintf(outfile,"CELL_TYPES %d\n",numel) ;
for(i=0;i<numel;i++)
   {
    fprintf(outfile,"10\n") ;
   }
fprintf(outfile,"\n") ;
printf("past ienfile");


/* **************************************************
    DISPFILE: First line containing the timestep; subsequent
    lines for each displacement. */

fprintf(outfile,"POINT_DATA %d\n",numnp) ;
fprintf(outfile,"VECTORS displacement float\n",numnp) ;
fscanf(dispfile,"%lf" , &dt) ;
printf("Time step=%g\n",dt) ;
for(i=0;i<numnp;i++)
   {
    fscanf(dispfile, "%s%d%lf%lf%lf%lf%lf%lf%lf%lf%lf" , buf , &node , &xx , &yy , &zz ,
            &sx , &sy , &sz , &vx , &vy , &vz) ;
    fprintf(outfile,"%g   %g   %g\n",sx,sy,sz) ;
   }
fprintf(outfile,"\n") ;
rewind(dispfile) ;

fprintf(outfile,"VECTORS velocity float\n",numnp) ;
fscanf(dispfile,"%g" , &dt) ;
for(i=0;i<numnp;i++)
   {
    fscanf(dispfile, "%s%d%lf%lf%lf%lf%lf%lf%lf%lf%lf" , buf , &node , &xx , &yy , &zz ,
            &sx , &sy , &sz , &vx , &vy , &vz) ;
    vx /= dt ;
    vy /= dt ;
    vz /= dt ;
    fprintf(outfile,"%g   %g   %g\n",vx,vy,vz) ;
   }
fprintf(outfile,"\n") ;


/* **************************************************
    STRFILE: No first line, just one line for every 
    stress. */

fprintf(outfile,"CELL_DATA %d\n",numel) ;
fprintf(outfile,"TENSORS stress float\n") ;
for(i=0;i<numel;i++)
   {
    fscanf(strfile, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%d" , &xx , &yy , &zz ,
            &st00 , &st11 , &st22 , &st01 , &st02 , &st12 , &nel) ;
    st10=st01 ;
    st20=st02 ;
    st21=st12 ;
    fprintf(outfile,"%g   %g   %g\n",st00,st01,st02) ;
    fprintf(outfile,"%g   %g   %g\n",st10,st11,st12) ;
    fprintf(outfile,"%g   %g   %g\n\n",st20,st21,st22) ;
   }


/* ************************************************** */

 fclose(coordfile) ;
 fclose(ienfile) ;
 fclose(dispfile) ;
 fclose(strfile) ;
 fclose(outfile) ;
 printf("Processing complete.\n") ;

}