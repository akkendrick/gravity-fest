/* fixedge.c  -  replace tapered fault edge nodes */
/* usage:   fixedge fltdataold edgelist fltdatanew */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(
     int argc ,      /* number of input arguments */
     char *argv[]    /* input argument strings */
    )
{
 int  numnp,numel,i,j,node,nel,ndum,n1,n2,n3,n4,offset,mat,flag,tally ;
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
tally = 0 ;

for(i=0;i<numnp;i++)
   {
    fscanf(infile, "%d%d%lf%lf%lf%lf%lf%lf%d%lf" , &node,&ndum,&xx,&yy,&zz,&sx,&sy,&sz,&n1,&st00) ;

    fscanf(listfile, "%d" , &numel) ;
    flag = 0 ;
    for(j=0;j<numel;j++)
       {
        fscanf(listfile, "%d%d%lf%lf%lf%lf%lf%lf%d%lf" , &nel,&ndum,&xx,&yy,&zz,&sx,&sy,&sz,&n1,&st00) ;
        if(nel == node)  flag = 1 ;
       }

    if(flag)
       {
        fprintf(outfile,"%d  1  0.0 0.0 1.0   -1.0 0.0 0.0   1 2.5\n",node) ;
        tally++ ;
       }
    
    else
       fprintf(outfile,"%d  1  0.0 0.0 1.0   -1.0 0.0 0.0   1 5.0\n",node) ;
    
    rewind(listfile) ;
   }
   
 fclose(infile) ;
 fclose(listfile) ;
 fclose(outfile) ;
 printf("# of changed entries = %d; Processing complete.\n",tally) ;

}
