/* gft2grav_vtk.c  -  convert GeoFEST output to vtk legacy file */
/* This version also does surface polygons for visualizing gravity changes */
/* usage:   gft2gvtk coordfile ienfile buoyfile dispfile strfile gravfile outfile.vtk */
/* version: GL 6/2014  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(
     int argc ,      /* number of input arguments */
     char *argv[]    /* input argument strings */
    )
{
 int  numnp,numel,i,node,nel,ndum,n1,n2,n3,n4,offset,mat,nfacet,numall,nf ;
 int *zienhold ;
 double  xx,yy,zz,sx,sy,sz,vx,vy,vz,dt,st00,st11,st22,st01,st10,st02,st20,st12,st21,
        xorg,yorg,zorg,rhog,dgx,dgy,dgz,dgnet,dN,Pzx,Pzy,Pzz ;
 char  buf[80] ;
 FILE    *coordfile , *ienfile , *dispfile , *strfile , *outfile , *buoyfile , *gravfile ;
 
 if(argc != 8)
    {
     printf("Wrong number of arguments.\n") ;
     printf("usage:   gft2gvtk coordfile ienfile buoyfile dispfile strfile gravfile outfile.vtu\n") ;
     exit(0) ;
    }

 coordfile = fopen(argv[1],"r") ;
 ienfile = fopen(argv[2],"r") ;
 buoyfile = fopen(argv[3],"r") ;
 dispfile = fopen(argv[4],"r") ;
 strfile = fopen(argv[5],"r") ;
 gravfile = fopen(argv[6],"r") ;
 outfile = fopen(argv[7],"w") ;
 
/* ************************************************** 
    COORDFILE: Expected to have a first line with the number
    of nodes; all subsequent lines are node # and coords. */

fprintf(outfile,"# vtk DataFile Version 2.3\n") ;
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
    4 node connectivity.

     BUOYFILE: Expected to have a first line with number of facets and
     three floats specifying center/upward of gravity and one more
     float specifying rho*g; all subsequent lines are el #, face #.  */

fscanf(ienfile, "%d" , &numel) ;
fscanf(buoyfile, "%d%lf%lf%lf%lf" , &nfacet , &xorg , &yorg , &zorg , &rhog) ;
if(nfacet < 0) {nfacet = -nfacet ;}

zienhold = (int *) calloc(4*numel,sizeof(int)) ;

ndum=numel*5 + nfacet*4 ;
numall = numel + nfacet ;
fprintf(outfile,"CELLS %d %d\n",numall,ndum) ;
for(i=0;i<numel;i++)
   {
    fscanf(ienfile, "%d%d%d%d%d%d%d" , &nel , &ndum , &mat , &n1 , &n2 , &n3 , &n4) ;
    n1-- ;
    n2-- ;
    n3-- ;
    n4-- ;
    fprintf(outfile,"4   %d   %d   %d   %d\n",n1,n2,n3,n4) ;
    *(zienhold+4*i+0) = n1 ;
    *(zienhold+4*i+1) = n2 ;
    *(zienhold+4*i+2) = n3 ;
    *(zienhold+4*i+3) = n4 ;
   }
for(i=0;i<nfacet;i++)
   {
    fscanf(buoyfile, "%d%d" , &nel , &nf) ;
    if(nf==1)
       {
        n1 = *(zienhold+4*(nel-1)+1) ;
        n2 = *(zienhold+4*(nel-1)+2) ;
        n3 = *(zienhold+4*(nel-1)+3) ;
       }
    else if(nf==2)
       {
        n1 = *(zienhold+4*(nel-1)+0) ;
        n2 = *(zienhold+4*(nel-1)+2) ;
        n3 = *(zienhold+4*(nel-1)+3) ;
       }
    else if(nf==3)
       {
        n1 = *(zienhold+4*(nel-1)+0) ;
        n2 = *(zienhold+4*(nel-1)+1) ;
        n3 = *(zienhold+4*(nel-1)+3) ;
       }
    else if(nf==4)
       {
        n1 = *(zienhold+4*(nel-1)+0) ;
        n2 = *(zienhold+4*(nel-1)+1) ;
        n3 = *(zienhold+4*(nel-1)+2) ;
       }
    fprintf(outfile,"3   %d   %d   %d\n",n1,n2,n3) ;
   }
fprintf(outfile,"\n") ;
fprintf(outfile,"CELL_TYPES %d\n",numall) ;
for(i=0;i<numel;i++)
   {
    fprintf(outfile,"10\n") ;
   }
for(i=0;i<nfacet;i++)
   {
    fprintf(outfile,"5\n") ;
   }
fprintf(outfile,"\n") ;
printf("past ienfile");


/* **************************************************
    DISPFILE: First line containing the timestep; subsequent
    lines for each displacement. */

fprintf(outfile,"POINT_DATA %d\n",numnp) ;
fprintf(outfile,"VECTORS displacement float\n") ;
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

fprintf(outfile,"VECTORS velocity float\n") ;
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
    GRAVFILE: contains X Y Z  dgx dgy dgz  dgnet dN  Pzx Pzy Pzz 
    only use dgnet and dN for now */

fprintf(outfile,"CELL_DATA %d\n",numel) ;
fprintf(outfile,"SCALARS gravity float 2\n") ;
for(i=0;i<numel;i++)
   {
    fprintf(outfile,"0   0\n") ;
   }
for(i=0;i<nfacet;i++)
   {
    fscanf(gravfile, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf" , &xx , &yy , &zz ,
            &dgx , &dgy , &dgz , &dgnet , &dN , &Pzx , &Pzy , &Pzz) ;
    fprintf(outfile,"%g   %g\n",dgnet,dN) ;
   }


/* ************************************************** */

/* **************************************************/
/* Leave out for now...
    STRFILE: No first line, just one line for every 
    stress. 

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

**** */
/* ************************************************** */

 fclose(coordfile) ;
 fclose(ienfile) ;
 fclose(buoyfile) ;
 fclose(dispfile) ;
 fclose(strfile) ;
 fclose(gravfile) ;
 fclose(outfile) ;
 printf("Processing complete.\n") ;

}