/*
***                         File solver.c                         ***
***                        GeoFEST version 6.0
*** Copyright (c) 2010, California Institute of Technology        ***
*** U.S.Sponsorship under NASA Contract NAS7-1407 is acknowledged ***
***
*** This software is designated for public release under JPL Task ***
*** Order Number NMO710991 and may be publicly released through
*** license with the Open Channel Foundation
***
*** This file contains modules for the sparse matrix storage,
*** substructure reduction or other solvers of the matrix system:
***
***          == solver ==
***          == addstiff ==
***          == addfor ==
***          == estiffprod ==
***          == reorder ==
***          == genrcm ==
***          == rcm ==
***          == degree ==
***          == fnroot ==
***          == rootls ==
***          == colht ==
***          == profile_diag ==
***          == factor ==
***          == full_back ==
***          == pcg_loop ==
***          == converged ==
***          == put_soln ==
***/

#define EXTERN extern
#include "solver.h"
#include "utility.h"

void elgrp_loop(
                int    task    /* which element-level task to perform */
               );

int compare_integers (const void *a, const void *b);

#define ES(i,j) *(es + nee*(j) + (i))

static int   first_time = TRUE ;

int compare_integers (const void *a, const void *b) {
   return (int) (*(int *)a - *(int *)b);
}



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   solver   ***************
*/

/*  INTERFACE:   */

     void
     solver(
            int  code     /* which solver option to use */
           )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine solver
** solver performs single right-hand-side matrix solution,
** currently using factorization and back substitution on 
** single processor OR using PCG iteration.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

    {
     switch(code)
        {
         case ITER:
           pcg_loop();
           completion() ;
           break;
           
         case FACBACK:
         case BACK:
          if( code == FACBACK )
           factor(stiff.matrix , stiff.diag , fe_sys.neq ) ;
     
           full_back(stiff.matrix , force.full_rhs , stiff.diag ,
               fe_sys.neq) ;     
           break;
        }
     printf("Out of solver; entering put_soln()...\n");
     

     put_soln( &force , global.id_pointer ) ;
    }

/************************** end of solver ************************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   addstiff   ***************
*/

/*  INTERFACE:   */

   void
   addstiff(
            PROFILE  *a ,  /* ptr to assembled profile stiffness storage */
            real  *es ,        /* ptr to element stiffness storage */
            ELEMENT_INFO  *info ,   /* ptr to element info struct */
            ELEMENT_DATA  *el_pt    /* ptr to current element struct */
           )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine addstiff
** addstiff adds element stiffness to global profile array,
** or to element storage if using PCG solver
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/
     
      {
       int  i_eq , j_eq , i , j , iadd , nee ;
       real  *stiff , *M ;

       nee = info->nee ;
       
       switch(fe_sys.solver)
          {
           case DIRECT:
           for (  j = 0 ; j < nee ; j++ )
              {
               if( (j_eq = (int) el_pt->lm[j]) != -1 )
                 {
                  for( i = 0 ; i <= j ; i++ )
                      {
                       if( (i_eq = (int) el_pt->lm[i]) != -1 )
                         {
                          if(j_eq >= i_eq)
                              {
                               iadd = (int) *((a->diag) + j_eq ) -
                                      j_eq + i_eq ;
                               *((a->matrix)+iadd) +=  ES(i,j) ;
                               }
                            else
                               {
                                iadd = (int) *((a->diag) + i_eq ) -
                                                i_eq + j_eq ;
                                *((a->matrix)+iadd) +=  ES(i,j) ;
                               } 
                         }
                     }
                  }
               }
           break ;
           
           case PCG:
           stiff = el_pt->stiff ;
           M = pcg.precond ;
           for (  j = 0 ; j < nee ; j++ )
               {
                  for( i = 0 ; i <= j ; i++ )
                     {
                      *(stiff + nee*j + i) = ES(i,j) ; /* store element stiffness */
                      
                      if( i == j && ( (i_eq = (int) el_pt->lm[i]) != -1 ) )
                         {
                          /* this is a diagonal element -- use in preconditioner */

                          *(M+i_eq) += ES(i,j) ;

                          /* Replaced by direct value of partial contributions
                          *            to be inverted later....
                          *if(*(M+i_eq) == ZERO)
                          *   {
                          *    *(M+i_eq) = ONE/ES(i,j) ;
                          *   }
                          *else
                          *   {
                          *    *(M+i_eq) = ONE/((ONE/ *(M+i_eq)) + ES(i,j)) ;
                          *   }
                          */

                         }
                     }
               }
           break ;
    
          } /* end switch */
               
                 
   }

/************************** end of addstiff **********************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   addfor   ***************
*/

/*  INTERFACE:   */

     void
     addfor(
            real  *dest ,  /* destination assembled rhs vector */
            real  *rhs ,        /* source element rhs vector */
            ELEMENT_INFO  *info ,   /* ptr to element info struct */
            ELEMENT_DATA  *el_pt    /* ptr to current element struct */
           )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine addfor
** addfor assembles the global r.h.s. vector.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      COUNTER      i    ;
      int          k    ;

      for ( i = 0 ; i < info->nee ; i++ )
         {
          k = el_pt->lm[i] ;
          if( k != (-1) )
             {
              *(dest + k) += rhs[i] ;
             }
         }
     }

/************************** end of addfor ************************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   estiffprod   ***************
*/

/*  INTERFACE:   */

 void
 estiffprod(
            GROUP   *grp_ptr    /* pointer to current element group */
           )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine estiffprod
** calculates the product of a given vector with stiffness in element storage
**    t = A * d
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

  {
   int     numel , nel , nee , neq , i , j , i_eq , j_eq ;
   real    *t , *d , *stiff, *piton, *crampon ;
   int *llm;
   ELEMENT_DATA   *el_pt ,*el_end;

   numel = (grp_ptr->el_info)->numel ;
   nee = (grp_ptr->el_info)->nee ;
   neq = fe_sys.neq ;
   t = pcg.temp ; /* destination vector */
   d = pcg.d ;  /* multiplicand vector */
   el_pt = grp_ptr->el_data;
   el_end = el_pt + numel; 
   for (  ; el_pt < el_end ; el_pt++ )
      {
       stiff = el_pt->stiff ;
       llm = el_pt->lm;
       /*
       * Use loops for any case on boundary or exceptional equation count.
       */
       if(el_pt->is_bc == 1 || nee != 12)
          {
          for (  j = 0 ; j < nee ; j++ )
             {
              if( (j_eq = llm[j]) != -1 )
                 {
                  for( i = 0 ; i < j ; i++ )
                     {
                      if( (i_eq =  llm[i]) != -1 )
                         {
                          *(t+i_eq) += *(d+j_eq) * *(stiff + nee*j + i) ;
                          *(t+j_eq) += *(d+i_eq) * *(stiff + nee*j + i) ;
                         }
                     } /* i loop */
                  *(t+j_eq) += *(d+j_eq) * *(stiff + nee*j + j) ;
                 }
             } /* j loop */
          } else {
                     /* Low-operation in-line commands for special case
                     *   of no boundary, nee==12
                     *
                     *  Following the older code, the assumptions here are
                     *  (a) stiff is allocated as a square, 12x12
                     *  (b) only the upper "half" gets used (i <= j)
                     *  (c) matrix is to be conceived as running down columns,
                     *     then across rows.  So "i" marks row.
                     *
                     *  So the goal is to fill the relevant parts of (neq-length)  t 
                     *  (pcg.temp, the destination)
                     *  with matrix-vector product stiff times d (ie, A*d)
                     *
                     *  Since t is accessed via the element's lm array, it won't
                     * have a convenient pattern.  So aim to finish each position
                     * before moving on.  
                     *
                     * "d" is accessed equally inconveniently (by lm), but we can store
                     * in d0..d11: local space that won't page out of cache, maybe registers
                     * if compiler is clever. 
                     *
                     * To access a full row of A without reading lower half,
                     * symmetry indicates we want to take items in column down to
                     * the diagonal, then move across row (see use of piton, crampon below).
                     */

                  real d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11;

                  
                  *(t+llm[0]) += (d0 = *(d+llm[0])) * *(stiff)
                              +  (d1 = *(d+llm[1])) * *(stiff +  nee) 
                              +  (d2 = *(d+llm[2])) * *(stiff + nee*2) 
                              +  (d3 = *(d+llm[3])) * *(stiff + nee*3) 
                              +  (d4 = *(d+llm[4])) * *(stiff + nee*4) 
                              +  (d5 = *(d+llm[5])) * *(stiff + nee*5) 
                              +  (d6 = *(d+llm[6])) * *(stiff + nee*6) 
                              +  (d7 = *(d+llm[7])) * *(stiff + nee*7) 
                              +  (d8 = *(d+llm[8])) * *(stiff + nee*8) 
                              +  (d9 = *(d+llm[9])) * *(stiff + nee*9) 
                              +  (d10 = *(d+llm[10])) * *(stiff + nee*10) 
                              +  (d11 = *(d+llm[11])) * *(stiff + nee*11); 

                    /* 
                    *   piton, crampon: pointers to landmark locations within the 
                    *   stiffness triangular matrix.
                    *
                    *   piton marks the top of the current column.
                    *
                    *   crampon marks the left edge of the current row.
                    *
                    *   Next item takes two items from column, then 10 from row.
                    *   Each succedding one takes one more from column, one less from row .
                    */

                  piton = stiff + nee;
                  crampon = stiff + 1;
                  *(t+llm[1]) += d0 * *piton
                              + d1 * *(piton + 1) 
                              + d2 * *(crampon + nee*2 ) 
                              + d3 * *(crampon + nee*3 )
                              + d4 * *(crampon + nee*4 )
                              + d5 * *(crampon + nee*5 )
                              + d6 * *(crampon + nee*6 )
                              + d7 * *(crampon + nee*7 )
                              + d8 * *(crampon + nee*8 )
                              + d9 * *(crampon + nee*9 )
                              + d10 * *(crampon + nee*10 )
                              + d11 * *(crampon + nee*11 );
                  piton += nee;
                  crampon++;
                  *(t+llm[2]) += d0 * *piton
                              + d1 * *(piton + 1) 
                              + d2 * *(piton + 2) 
                              + d3 * *(crampon + nee*3 )
                              + d4 * *(crampon + nee*4 )
                              + d5 * *(crampon + nee*5 )
                              + d6 * *(crampon + nee*6 )
                              + d7 * *(crampon + nee*7 )
                              + d8 * *(crampon + nee*8 )
                              + d9 * *(crampon + nee*9 )
                              + d10 * *(crampon + nee*10 )
                              + d11 * *(crampon + nee*11 );
                  piton += nee;
                  crampon++;
                  *(t+llm[3]) += d0 * *piton
                              + d1 * *(piton + 1) 
                              + d2 * *(piton + 2) 
                              + d3 * *(piton + 3)
                              + d4 * *(crampon + nee*4 )
                              + d5 * *(crampon + nee*5 )
                              + d6 * *(crampon + nee*6 )
                              + d7 * *(crampon + nee*7 )
                              + d8 * *(crampon + nee*8 )
                              + d9 * *(crampon + nee*9 )
                              + d10 * *(crampon + nee*10 )
                              + d11 * *(crampon + nee*11 );
                  piton += nee;
                  crampon++;
                  *(t+llm[4]) += d0 * *piton
                              + d1 * *(piton + 1) 
                              + d2 * *(piton + 2) 
                              + d3 * *(piton + 3)
                              + d4 * *(piton + 4)
                              + d5 * *(crampon + nee*5 )
                              + d6 * *(crampon + nee*6 )
                              + d7 * *(crampon + nee*7 )
                              + d8 * *(crampon + nee*8 )
                              + d9 * *(crampon + nee*9 )
                              + d10 * *(crampon + nee*10 )
                              + d11 * *(crampon + nee*11 );
                  piton += nee;
                  crampon++;
                  *(t+llm[5]) += d0 * *piton
                              + d1 * *(piton + 1) 
                              + d2 * *(piton + 2) 
                              + d3 * *(piton + 3)
                              + d4 * *(piton + 4)
                              + d5 * *(piton + 5)
                              + d6 * *(crampon + nee*6 )
                              + d7 * *(crampon + nee*7 )
                              + d8 * *(crampon + nee*8 )
                              + d9 * *(crampon + nee*9 )
                              + d10 * *(crampon + nee*10 )
                              + d11 * *(crampon + nee*11 );
                  piton += nee;
                  crampon++;
                  *(t+llm[6]) += d0 * *piton
                              + d1 * *(piton + 1) 
                              + d2 * *(piton + 2) 
                              + d3 * *(piton + 3)
                              + d4 * *(piton + 4)
                              + d5 * *(piton + 5)
                              + d6 * *(piton + 6)
                              + d7 * *(crampon + nee*7 )
                              + d8 * *(crampon + nee*8 )
                              + d9 * *(crampon + nee*9 )
                              + d10 * *(crampon + nee*10 )
                              + d11 * *(crampon + nee*11 );
                  piton += nee;
                  crampon++;
                  *(t+llm[7]) += d0 * *piton
                              + d1 * *(piton + 1) 
                              + d2 * *(piton + 2) 
                              + d3 * *(piton + 3)
                              + d4 * *(piton + 4)
                              + d5 * *(piton + 5)
                              + d6 * *(piton + 6)
                              + d7 * *(piton + 7)
                              + d8 * *(crampon + nee*8 )
                              + d9 * *(crampon + nee*9 )
                              + d10 * *(crampon + nee*10 )
                              + d11 * *(crampon + nee*11 );
                  piton += nee;
                  crampon++;
                  *(t+llm[8]) += d0 * *piton
                              + d1 * *(piton + 1) 
                              + d2 * *(piton + 2) 
                              + d3 * *(piton + 3)
                              + d4 * *(piton + 4)
                              + d5 * *(piton + 5)
                              + d6 * *(piton + 6)
                              + d7 * *(piton + 7)
                              + d8 * *(piton + 8)
                              + d9 * *(crampon + nee*9 )
                              + d10 * *(crampon + nee*10 )
                              + d11 * *(crampon + nee*11 );
                  piton += nee;
                  crampon++;
                  *(t+llm[9]) += d0 * *piton
                              + d1 * *(piton + 1) 
                              + d2 * *(piton + 2) 
                              + d3 * *(piton + 3)
                              + d4 * *(piton + 4)
                              + d5 * *(piton + 5)
                              + d6 * *(piton + 6)
                              + d7 * *(piton + 7)
                              + d8 * *(piton + 8)
                              + d9 * *(piton + 9)
                              + d10 * *(crampon + nee*10 )
                              + d11 * *(crampon + nee*11 );
                  piton += nee;
                  crampon++;
                  *(t+llm[10]) += d0 * *piton
                              + d1 * *(piton + 1) 
                              + d2 * *(piton + 2) 
                              + d3 * *(piton + 3)
                              + d4 * *(piton + 4)
                              + d5 * *(piton + 5)
                              + d6 * *(piton + 6)
                              + d7 * *(piton + 7)
                              + d8 * *(piton + 8)
                              + d9 * *(piton + 9)
                              + d10 * *(piton + 10)
                              + d11 * *(crampon + nee*11 );
                  piton += nee;
                  crampon++;
                  *(t+llm[11]) += d0 * *piton
                              + d1 * *(piton + 1) 
                              + d2 * *(piton + 2) 
                              + d3 * *(piton + 3)
                              + d4 * *(piton + 4)
                              + d5 * *(piton + 5)
                              + d6 * *(piton + 6)
                              + d7 * *(piton + 7)
                              + d8 * *(piton + 8)
                              + d9 * *(piton + 9)
                              + d10 * *(piton + 10)
                              + d11 * *(piton + 11);

          }
      } /* nel loop */

  }

/************************** end of estiffprod ********************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   reorder   ***************
*/

/*  INTERFACE:   */

    void
    reorder(
            GROUP   *grp_ptr    /* pointer to current element group */
           )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine reorder
** reorder uses ien information to build adjacency information and call
** permutation optimizing routines for minimizing matrix profile.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

{
   int *id;
   int index,
       *icol,
       *irow,
       *llist,
       *msk,
       *xls,
       *invper,
       nnodes = loc_sys.numnp,
       max_nz = MAX_CONNECT*nnodes,
       lfree,
       nel,
       idof,ndof,
       next_icol,
       num_col,
       total_nonzeros,
       pgnode,
       next_eq,
       this_link, next_link, ref_lnode, ref_gnode, aux_lnode, aux_gnode;

   ELEMENT_INFO  *info ;
   ELEMENT_DATA  *el_pt ;

/* FILE *tfp; */

   info = grp_ptr->el_info ;

/* use ien data to fill irow, icol arrays */

   global.permut = (int *) calloc(nnodes+1,sizeof(int));
   icol = (int *) calloc(max_nz,sizeof(int));
   irow = (int *) calloc(nnodes + 2,sizeof(int));
   llist = (int *) calloc(2*max_nz,sizeof(int));
   invper = (int *) calloc(nnodes+1,sizeof(int));
   msk = (int *) calloc(nnodes+1,sizeof(int));
   xls = (int *) calloc(nnodes+1,sizeof(int));


   lfree = nnodes + 1;

/* compile adjacency information from every element. . . */
   for( nel = 0 ;  nel < info->numel ; nel++ )
   {
      el_pt = grp_ptr->el_data + nel ;

/* . . . every pair of nodes for this element . . . */
/* NOTE CHECK FOR ZERO-BASED VS 1-BASED LISTS (Local, global nodes ) */

      for( ref_lnode = 0; ref_lnode < info->nen; ref_lnode++)
      {
         ref_gnode = el_pt->ien[ref_lnode];
         for( aux_lnode = 0; aux_lnode < info->nen; aux_lnode++)
         {
            aux_gnode = el_pt->ien[aux_lnode]; 

/* put this aux node into list for ref node.  Insert in node order.*/

            this_link = ref_gnode-1;
            next_link = llist[2*this_link + 1];

            while(next_link != 0 && llist[2*next_link+0] < aux_gnode)
            {
               this_link = next_link;
               next_link = llist[next_link*2 + 1];
            }

                /* note case llist[..) == aux_gnode 
                    requires no action, since link already exists.
                */

            if(next_link == 0 || llist[next_link*2+0] > aux_gnode)
            {
               if(lfree <= max_nz)
               {
                  llist[lfree*2+0] = aux_gnode;
                  llist[lfree*2 + 1] = next_link;
                  llist[this_link*2 + 1] = lfree++;
               }
               else
               {
                 printf("Bad assumptions on max_nz: rewrite alloc strategy?\n");
                  exit(0);
               } /* end if lfree . . .*/
            }    /* end if next_link . . .*/
         }       /* for aux_lnode . . . */
      }          /* for ref_lnode . . . */
   }             /* for nel. . . */

/* now we have the linked lists: build "icol" adjacency list, "irow" index */
   next_icol = 1;
   num_col = 0;
   for(ref_gnode = 1; ref_gnode <= nnodes; ref_gnode++)
   {
      this_link = llist[(ref_gnode-1)*2+1];
      irow[ref_gnode] = next_icol;
        
      while(this_link != 0)
      {
         icol[next_icol++]= llist[this_link*2+0];
         this_link = llist[2*this_link + 1];
         num_col++; 
      }
   }
   irow[nnodes+1] = next_icol;

   total_nonzeros  = irow[nnodes+1]-1;

/* call rcm optimizer */

   genrcm(nnodes,irow,icol,invper,msk,xls);

   for(ref_gnode = 1; ref_gnode <= nnodes; ref_gnode++)
      global.permut[invper[ref_gnode]] = ref_gnode;
/*
   for(ref_gnode = 1; ref_gnode <= nnodes; ref_gnode++)
      printf("gn = %d, pn = %d\n",ref_gnode,global.permut[ref_gnode]);
*/
/* write to a file to plot nonzero distribution */
/*
   tfp = fopen("nonzeros","w");
   for(ref_gnode = 1; ref_gnode < nnodes; ref_gnode++)
   {
      for(index = irow[ref_gnode]; index <= irow[ref_gnode+1]-1; index++)
      {
          fprintf(tfp,"%d\t%d\t%d\t%d\n",
              ref_gnode,icol[index],
              global.permut[ref_gnode],global.permut[icol[index]]);
      }
   }
   fclose(tfp);
*/
 
/* overwrite id array with permuted information */

   id = global.id_pointer ;

  /* now assign equations ...    (local equations start at zero)  */

     ndof = fe_sys.ndof;

     next_eq = 0 ;

     for( pgnode = 1 ; pgnode <=nnodes ; pgnode++ )
         {
          for(idof=0 ; idof<ndof ; idof++){
              id[ndof*(pgnode-1) + idof] = 
     ( *(global.eqtype + ndof*(invper[pgnode] - 1) + idof) != 0) ? next_eq++ : -1;
 /*             printf("pg = %d, idof = %d, id = %d\n",
                pgnode,idof,id[ndof*(pgnode-1) + idof]);
*/
           }

         }   /* end loop over permuted global node number */
      
     
/* overwrite lm arrays with permuted implied equations */

   for( nel = 0 ;  nel < info->numel ; nel++ )
   {
      el_pt = grp_ptr->el_data + nel ;
      index = 0 ;

      for( ref_lnode = 0 ; ref_lnode < info->nen ; ref_lnode++)
      {
         pgnode = global.permut[ el_pt->ien[ref_lnode] ] ;

         for( idof = 0 ; idof < ndof ; idof++){
            el_pt->lm[index++] = id[ndof*(pgnode-1) + idof];
         }
/*         printf("el = %d, ln = %d, ien = %d, lm0 = %d, lm1 = %d, lm2 = %d\n",
                    nel,ref_lnode,el_pt->ien[ref_lnode],
                                    el_pt->lm[ndof*ref_lnode + 0],
                                    el_pt->lm[ndof*ref_lnode + 1],
                                    el_pt->lm[ndof*ref_lnode + 2]);
*/
      }
   }
/* note that now the id array yields internal equations by permuted node,
   element ien yields the external (nonpermuted) node by local node,
   element lm yields the internal equations.
   If the element-wise list of permuted nodes is desired, must construct
   global.permut[el_pt->ien[lnode]].
*/
}
/************************** end of reorder ***********************************/ 


/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   genrcm   ***************
*/

/*  INTERFACE:   */

     void
     genrcm(
            int  neqns ,   /* NUMBER OF EQUATIONS */
            int  *xadj ,  /* ADJACENCY STRUCTURE OF THE GRAPH OF THE MATRIX */
            int  *adjncy ,/* ADJACENCY STRUCTURE OF THE GRAPH OF THE MATRIX */
            int  *perm ,   /* VECTOR THAT CONTAINS THE RCM ORDERING */
            int  *mask ,   /* MARK VARS NUMBERED IN ORDERING PROCESS */
            int  *xls   /* INDEX VECTOR FOR A LEVEL STRUCTURE */
           )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine genrcm
** genrcm computes the reverse cuthill mckee ordering for a general
** adjacency graph such as a sparse finite element matrix. 
**
** Original FORTRAN comments: (transcribed into C by JWP from
** p_reorder.f of PHOEBUS3_T3D):
*****************************************************************
*****************************************************************
**********   GENRCM ..... GENERAL REVERSE CUTHILL MCKEE   *******
*****************************************************************
*****************************************************************
**
**     PURPOSE - GENRCM FINDS THE REVERSE CUTHILL-MCKEE                     
**        ORDERING FOR A GENERAL GRAPH. FOR EACH CONNECTED                 
**        COMPONENT IN THE GRAPH, GENRCM OBTAINS THE ORDERING             
**        BY CALLING THE SUBROUTINE RCM.                                 
**
**     INPUT PARAMETERS -                                               
**        NEQNS - NUMBER OF EQUATIONS                                  
**        (XADJ, ADJNCY) - ARRAY PAIR CONTAINING THE ADJACENCY        
**               STRUCTURE OF THE GRAPH OF THE MATRIX.               
**
**     OUTPUT PARAMETER -                                         
**        PERM - VECTOR THAT CONTAINS THE RCM ORDERING.          
**
**     WORKING PARAMETERS -                                             
**        MASK - IS USED TO MARK VARIABLES THAT HAVE BEEN              
**               NUMBERED DURING THE ORDERING PROCESS. IT IS          
**               INITIALIZED TO 1, AND SET TO ZERO AS EACH NODE      
**               IS NUMBERED.                                       
**        XLS - THE INDEX VECTOR FOR A LEVEL STRUCTURE.  THE       
**               LEVEL STRUCTURE IS STORED IN THE CURRENTLY       
**               UNUSED SPACES IN THE PERMUTATION VECTOR PERM.   
**
**     PROGRAM SUBROUTINES -                                    
**        FNROOT, RCM.                                         
**
*****************************************************************
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

{
   int ccsize, i, nlvl,
       num = 1,
       root;

   for(i=1; i<=neqns; i++)
      mask[i] = 1;

/* for each masked connected component . . . */

   for(i=1; i<=neqns; i++)
   {
      if (mask[i] != 0) 
      {
         root = i;

/* first find a pseudo-peripheral node root.
   Note that the level structure found by 
   fnroot is stored starting at perm[num].
   Then rcm is called to order the component
   using root as the starting node.
*/

         fnroot( &root,xadj,adjncy,mask,&nlvl,xls,perm,num );

         rcm( root,xadj,adjncy,mask,perm,num,&ccsize,xls );

         num += ccsize;
      
         if (num > neqns) 
            break;
      }
   }
}
/************************** end of genrcm ************************************/ 


/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   rcm   ***************
*/

/*  INTERFACE:   */

        void
        rcm(
            int  root ,   /* NODE THAT DEFINES THE CONNECTED COMPONENT */
            int  *xadj ,  /* ADJACENCY STRUCTURE OF THE GRAPH OF THE MATRIX */
            int  *adjncy ,/* ADJACENCY STRUCTURE OF THE GRAPH OF THE MATRIX */
            int  *mask ,   /* MARK VARS NUMBERED IN ORDERING PROCESS */
            int  *perm ,   /* VECTOR THAT CONTAINS THE RCM ORDERING */
            int  num ,   /*  */
            int  *ccsize ,   /* SIZE OF THE CONNECTED COMPONENT */
            int  *deg   /* TEMPORARY VECTOR USED TO HOLD THE DEGREE */
           )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** rcm does connectivity analysis (jwp interpretation; see below).
** Original comments from Fortran code:
*****************************************************************
*****************************************************************
**********     RCM ..... REVERSE CUTHILL-MCKEE ORDERING   *******
*****************************************************************
*****************************************************************
**
**     PURPOSE - RCM NUMBERS A CONNECTED COMPONENT SPECIFIED BY             7.
**        MASK AND ROOT, USING THE RCM ALGORITHM.                           8.
**        THE NUMBERING IS TO BE STARTED AT THE NODE ROOT.                  9.
**
**     INPUT PARAMETERS -                                                  11.
**        ROOT - IS THE NODE THAT DEFINES THE CONNECTED                    12.
**               COMPONENT AND IT IS USED AS THE STARTING                  13.
**               NODE FOR THE RCM ORDERING.                                14.
**        (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR FOR                    15.
**               THE GRAPH.                                                16.
**
**     UPDATED PARAMETERS -                                                18.
**        MASK - ONLY THOSE NODES WITH NONZERO INPUT MASK                  19.
**               VALUES ARE CONSIDERED BY THE ROUTINE.  THE                20.
**               NODES NUMBERED BY RCM WILL HAVE THEIR                     21.
**               MASK VALUES SET TO ZERO.                                  22.
**
**     OUTPUT PARAMETERS -                                                 24.
**        PERM - WILL CONTAIN THE RCM ORDERING.                            25.
**        CCSIZE - IS THE SIZE OF THE CONNECTED COMPONENT                  26.
**               THAT HAS BEEN NUMBERED BY RCM.                            27.
**
**     WORKING PARAMETER -                                                 29.
**        DEG - IS A TEMPORARY VECTOR USED TO HOLD THE DEGREE              30.
**               OF THE NODES IN THE SECTION GRAPH SPECIFIED               31.
**               BY MASK AND ROOT.                                         32.
**
**     PROGRAM SUBROUTINES -                                               34.
**        DEGREE.                                                          35.
**
*****************************************************************
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

{
   int fnbr,i,j,jstop,jstrt,k,l,lbegin,lnbr,lperm,lvlend,nbr,node;

/* find the degrees of the nodes in the
   component specified by mask and root.
*/

   degree( root,xadj,adjncy,mask,deg,ccsize,perm,num );

/* watch out for pass-by-value:  if anything gets changed here, need
to pass address up, not value .. . */

   mask[root] = 0;
   if (*ccsize <= 1) return ;
   lvlend = 0;
   lnbr = 1;

   do
   {

/* lbegin and lvlend point to the beginning and
   the end of the current level respectively.
*/

      lbegin = lvlend + 1;
      lvlend = lnbr;

/* for each node in the current level. . . 
*/

      for(i=lbegin; i<=lvlend; i++)
      {
         node = perm[i+num-1];
         jstrt = xadj[node];
         jstop = xadj[node+1] - 1;

/* Find the unnumbered neighbors of node.
   fnbr and lnbr point the first and last
   unnumbered neighbors respectively of the current
   node in perm.
*/
         fnbr = lnbr + 1;
         for(j=jstrt; j<=jstop; j++)
         {
            nbr = adjncy[j];
            if( mask[nbr] == 0)
               continue;
            lnbr++;
            mask[nbr] = 0;
            perm[lnbr+num-1] = nbr;
         }
         if(fnbr >= lnbr) 
            continue;

/* Sort the neighbors of node in increasing 
   order by degree.  Linear insertion is used.
*/

         k = fnbr;
         do
         {
            l = k++;
            nbr = perm[k+num-1];
            while( l >= fnbr )  
            {
               lperm = perm[l+num-1];
               if( deg[lperm] <= deg[nbr] )
                  break;
               else
               {
                  perm[l+num] = lperm;
                  l--;
               }
            }
            perm[l+num] = nbr;
         } while(k < lnbr);
      }
   } while(lnbr > lvlend);

/* We now have the Cuthill McKee ordering.
   Reverse it below. . . 
*/

   k = *ccsize/2;
   l = *ccsize;
   for(i=1; i<=k; i++)
   {
      lperm             = perm[l + num - 1];
      perm[l + num - 1] = perm[i + num - 1];
      perm[l + num - 1] = lperm;
      l--;
   }
}
/************************** end of rcm ***************************************/ 
         

/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   degree   ***************
*/

/*  INTERFACE:   */

     void
     degree(
            int  root ,   /* NODE THAT DEFINES THE CONNECTED COMPONENT */
            int  *xadj ,  /* ADJACENCY STRUCTURE OF THE GRAPH OF THE MATRIX */
            int  *adjncy ,/* ADJACENCY STRUCTURE OF THE GRAPH OF THE MATRIX */
            int  *mask ,   /* SPECIFIES A SECTION SUBGRAPH */
            int  *deg  , /* ARRAY CONTAINING THE DEGREES OF THE NODES */
            int  *ccsize ,   /* SIZE OF THE COMPONENT SPECIFED */
            int  *ls ,   /* TEMPORARY VECTOR USED TO STORE THE NODES */
            int  num     /*  */
           )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine degree
** degree computes the graph steps to each node in the subgraph. 
** Original Comments from Fortran code:
*****************************************************************
*****************************************************************
**********     DEGREE ..... DEGREE IN MASKED COMPONENT   ********
*****************************************************************
*****************************************************************
**
**     PURPOSE - THIS ROUTINE COMPUTES THE DEGREES OF THE NODES             7.
**        IN THE CONNECTED COMPONENT SPECIFIED BY MASK AND ROOT.            8.
**        NODES FOR WHICH MASK IS ZERO ARE IGNORED.                         9.
**
**     INPUT PARAMETER -                                                   11.
**        ROOT - IS THE INPUT NODE THAT DEFINES THE COMPONENT.             12.
**        (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR.                       13.
**        MASK - SPECIFIES A SECTION SUBGRAPH.                             14.
**
**     OUTPUT PARAMETERS -                                                 16.
**        DEG - ARRAY CONTAINING THE DEGREES OF THE NODES IN               17.
**              THE COMPONENT.                                             18.
**        CCSIZE-SIZE OF THE COMPONENT SPECIFED BY MASK AND ROOT           19.
**
**     WORKING PARAMETER -                                                 21.
**        LS - A TEMPORARY VECTOR USED TO STORE THE NODES OF THE           22.
**               COMPONENT LEVEL BY LEVEL.                                 23.
**
*****************************************************************
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

{
   int i, ideg, j, jstop, jstrt,lbegin,lvlend,lvsize, nbr,node;

/* initialization . . .
   The array xadj is used as a temporary marker to 
   indicate which nodes have been considered so far.
*/

   ls[num] = root;
   xadj[root] *= -1;
   lvlend = 0;
   *ccsize = 1;

   do
   {

/* lbegin is the pointer to the beginning of the current
   level, and lvlend points to the end of this level.
*/

      lbegin = lvlend+1;
      lvlend = *ccsize;

/* Find the degrees of nodes in the current level,
   and at the same time, generate the next level.
*/

      for(i=lbegin; i<=lvlend; i++)
      {
         node = ls[i + num - 1];
         jstrt = -xadj[node];
         jstop = abs(xadj[node+1]) - 1;
         ideg = 0;
         if( jstop < jstrt) {
            deg[node] = ideg;
            continue;
          }
         for(j=jstrt; j<= jstop; j++)
         {
            nbr = adjncy[j];
            if( mask[nbr] == 0 ) continue;
            ideg++;
            if( xadj[nbr] < 0 ) continue;
            xadj[nbr] *= -1;
            (*ccsize)++;
            ls[*ccsize + num - 1] = nbr;
         }
         deg[node] = ideg;
      }

/* Compute the current level width.
   If it is nonzero, generate another level.
*/

      lvsize = *ccsize - lvlend;
   } while(lvsize > 0);

/* Reset xadj to its correct sign and return.
*/

   for(i=1; i<=*ccsize; i++)
   {
      node = ls[i + num - 1];
      xadj[node] *= -1;
   }
}
/************************** end of degree ************************************/ 


/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   fnroot   ***************
*/

/*  INTERFACE:   */

     void
     fnroot(
            int  *root ,   /* DEFINES THE COMP FOR PSEUDO-PERIPH NODE */
            int  *xadj ,  /* ADJACENCY STRUCTURE OF THE GRAPH OF THE MATRIX */
            int  *adjncy ,/* ADJACENCY STRUCTURE OF THE GRAPH OF THE MATRIX */
            int  *mask ,   /* SPECIFIES A SECTION SUBGRAPH */
            int  *nlvl ,  /* NUMBER OF LEVELS IN THE LEVEL STRUCTURE */
            int  *xls ,   /* LEVEL STRUCTURE ARRAY PAIR */
            int  *ls ,   /* LEVEL STRUCTURE ARRAY PAIR */
            int  num     /*  */
           )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine fnroot
** fnroot finds the pseudo-peripheral node for a given subgraph.
** Original comments from Fortran:
*****************************************************************
*****************************************************************
*********     FNROOT ..... FIND PSEUDO-PERIPHERAL NODE    *******
*****************************************************************
*****************************************************************
**
**    PURPOSE - FNROOT IMPLEMENTS A MODIFIED VERSION OF THE                 7.
**       SCHEME BY GIBBS, POOLE, AND STOCKMEYER TO FIND PSEUDO-             8.
**       PERIPHERAL NODES.  IT DETERMINES SUCH A NODE FOR THE               9.
**       SECTION SUBGRAPH SPECIFIED BY MASK AND ROOT.                      10.
**
**    INPUT PARAMETERS -                                                   12.
**       (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR FOR THE GRAPH.          13.
**       MASK - SPECIFIES A SECTION SUBGRAPH. NODES FOR WHICH              14.
**              MASK IS ZERO ARE IGNORED BY FNROOT.                        15.
**
**    UPDATED PARAMETER -                                                  17.
**       ROOT - ON INPUT, IT (ALONG WITH MASK) DEFINES THE                 18.
**              COMPONENT FOR WHICH A PSEUDO-PERIPHERAL NODE IS            19.
**              TO BE FOUND. ON OUTPUT, IT IS THE NODE OBTAINED.           20.
**
**    OUTPUT PARAMETERS -                                                  22.
**       NLVL - IS THE NUMBER OF LEVELS IN THE LEVEL STRUCTURE             23.
**              ROOTED AT THE NODE ROOT.                                   24.
**       (XLS,LS) - THE LEVEL STRUCTURE ARRAY PAIR CONTAINING              25.
**                  THE LEVEL STRUCTURE FOUND.                             26.
**
**    PROGRAM SUBROUTINES -                                                28.
**       ROOTLS.                                                           29.
**
*****************************************************************
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

{
   int ccsize, j, jstrt, k, kstop, kstrt, mindeg, nabor, ndeg, node, nunlvl;

/* determine the level structure rooted at root.
*/

   rootls(*root,xadj,adjncy,mask,nlvl,xls,ls,num);
   ccsize = xls[*nlvl+1] - 1;
   if(*nlvl == 1 || *nlvl == ccsize) return ;

/* Pick a node with minimum degree from the last level.
*/
   do
   {
      jstrt = xls[*nlvl];
      mindeg = ccsize;
      *root = ls[jstrt + num - 1];
      if( ccsize == jstrt )
      {
         rootls(*root,xadj,adjncy,mask,&nunlvl,xls,ls,num);
         if(nunlvl <= *nlvl) return ;
         *nlvl = nunlvl;
         if( *nlvl < ccsize) continue;
         return ;
      }
      for(j=jstrt; j<=ccsize; j++)
      {
         node = ls[j+num-1];
         ndeg = 0;
         kstrt = xadj[node];
         kstop = xadj[node+1] - 1;
         for(k=kstrt; k<= kstop; k++)
         {
             nabor = adjncy[k];
             if(mask[nabor] > 0) ndeg++;
         }
         if(ndeg >= mindeg) continue;
         *root = node;
         mindeg = ndeg;
      }

/* and generate its rooted level structure.
*/

      rootls(*root,xadj,adjncy,mask,&nunlvl,xls,ls,num);
      if(nunlvl <= *nlvl) return ;

      *nlvl = nunlvl;

   } while(*nlvl < ccsize);
}
/************************** end of fnroot ************************************/ 

      
/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   rootls   ***************
*/

/*  INTERFACE:   */

     void
     rootls(
            int  root ,   /* NODE AT WHICH LEVEL STRUCTURE IS ROOTED */
            int  *xadj ,  /* ADJACENCY STRUCTURE OF THE GRAPH OF THE MATRIX */
            int  *adjncy ,/* ADJACENCY STRUCTURE OF THE GRAPH OF THE MATRIX */
            int  *mask ,   /* SPECIFIES A SECTION SUBGRAPH */
            int  *nlvl ,  /* NUMBER OF LEVELS IN THE LEVEL STRUCTURE */
            int  *xls ,   /* LEVEL STRUCTURE ARRAY PAIR */
            int  *ls ,   /* LEVEL STRUCTURE ARRAY PAIR */
            int  num     /*  */
           )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine rootls
** rootls generates the level structure corresponding to "root".
** Original Fortran comments:
*****************************************************************
*****************************************************************
**********     ROOTLS ..... ROOTED LEVEL STRUCTURE      *********
*****************************************************************
*****************************************************************
**
**     PURPOSE - ROOTLS GENERATES THE LEVEL STRUCTURE ROOTED
**        AT THE INPUT NODE CALLED ROOT. ONLY THOSE NODES FOR
**        WHICH MASK IS NONZERO WILL BE CONSIDERED.
**
**     INPUT PARAMETERS -
**        ROOT - THE NODE AT WHICH THE LEVEL STRUCTURE IS TO
**               BE ROOTED.
**        (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR FOR THE
**               GIVEN GRAPH.
**        MASK - IS USED TO SPECIFY A SECTION SUBGRAPH. NODES
**               WITH MASK(I)=0 ARE IGNORED.
**
**     OUTPUT PARAMETERS -
**        NLVL - IS THE NUMBER OF LEVELS IN THE LEVEL STRUCTURE.
**        (XLS, LS) - ARRAY PAIR FOR THE ROOTED LEVEL STRUCTURE.
**
*****************************************************************
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

{
   int i,j,jstop,jstrt,lbegin,ccsize,lvlend,lvsize,nbr,node;

/* initialization . . .*/

   mask[root] = 0;
   ls[num] = root;
   *nlvl = 0;
   lvlend = 0;
   ccsize = 1;

   do
   {

/* lbegin is the pointer to the beginning of the current
   level, and lvlend points to the end of this level.
*/

      lbegin = lvlend + 1;
      lvlend = ccsize;
      (*nlvl)++;
      xls[*nlvl] = lbegin;

/* generate the next level by finding all the masked
   neighbors of nodes in the current level.
*/
      for(i=lbegin; i<= lvlend; i++)
      {
         node = ls[i+num-1];
         jstrt = xadj[node];
         jstop = xadj[node+1] - 1;
         if(jstop < jstrt) continue;
         for(j=jstrt; j<=jstop; j++)
         {
            nbr = adjncy[j];
            if(mask[nbr] == 0) continue;
            ccsize++;
            ls[ccsize+num-1] = nbr;
            mask[nbr] = 0;
         }
      }

/* Compute the current level width.
   If it is nonzero, generate the next level.
*/

      lvsize = ccsize - lvlend;
   } while (lvsize > 0);

   xls[*nlvl+1] = lvlend + 1;

/* Reset mask to one for the nodes in the level structure.
*/
   for(i=1; i<=ccsize; i++)
   {
      node = ls[i+num-1];
      mask[node] = 1;
   }
}
/************************** end of rootls ************************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   colht   ***************
*/

/*  INTERFACE:   */

      void
      colht(
            GROUP    *grp_ptr  ,  /* pointer to current element group */
            PROFILE  *a    /* ptr to assembled profile stiffness storage */
           )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine colht
** colht computes column heights of global array.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      ELEMENT_INFO  *info ;
      ELEMENT_DATA  *el_pt ;
      int    nel , min , i , num , m ;
 
/*    FILE *tfp;  */

      info = grp_ptr->el_info ;

      for( nel = 0 ;  nel < info->numel ; nel++ )   
          {
           min = 100000000 ;
           el_pt = grp_ptr->el_data + nel ;
           for ( i = 0 ; i < info->nee ; i++)
               {
                if( (num = (int) el_pt->lm[i] ) != -1 )
                   min = MIN(num,min) ;
               }

           for ( i = 0 ; i < info->nee ; i++ )
               {
                if( ( num = (int) el_pt->lm[i] ) != -1 )
                  if( ( m = num - min ) > (int) *(a->diag + num) )
                           *(a->diag + num) = (int) m ;
               }
            }

/* report colhts to file for debugging */
/*
        tfp = fopen("colhts","w");
        for(i=0;i<fe_sys.neq;i++)
          fprintf(tfp,"%d\t%d\n", i, a->diag[i]);
        fclose(tfp);
*/
      }

/************************** end of colht *************************************/ 


/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   profile_diag   ***************
*/

/*  INTERFACE:   */

 void
 profile_diag(
              PROFILE  *a  ,  /* ptr to assembled profile stiffness storage */
              int   profile_case    /* old flag for substructuring... */
             )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine profile_diag
** profile diag computes the diagonal addresses.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

    {
      int i, neq ;

      *(a->diag) = 0 ;
      neq       = a->neq  ;
      for( i=1 ; i < neq ; i++ ) 
          {
           if( profile_case == INTERIOR )
               *(a->diag + i) += *(a->diag + i - 1) + 1 ;
           else
               *(a->diag + i) = *(a->diag + i - 1) + (int ) i + 1 ;
          }
      a->size = (int) *(a->diag + neq - 1) + 1 ;
    }

/************************** end of profile_diag ******************************/ 


/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   factor   ***************
*/

/*  INTERFACE:   */
       void
       factor(
              real  *stiff ,         /* assembled profile matrix */
              int   *diag ,           /* array of diagonal addresses */
              int   number_of_eqs     /* number of equations */
             )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine factor
** factor performs profile-based factorization:
**                          t                       
**                     a=(u) * d * u (crout) 
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      int j_row , j , j_diag ,j_ht ,i_start , i_end ,
          k , i_diag  , i , i_ht , i_row ;
      real d ;

      j_row = (-1) ;
      for (j = 0 ; j < number_of_eqs; j++ )
          {
           j_diag = (int) *(diag + j) ;
           j_ht   = j_diag - j_row ;
           i_start = j - j_ht + 2  ;
           if( j_ht - 2 > 0 )
               {
                i_end = j - 1 ;
                k = j_row + 2 ;
                i_diag = (int) *(diag + i_start - 1) ;
                for ( i = i_start ; i <= i_end ; i++ )
                    {
                     i_row = i_diag ;
                     i_diag = (int) *(diag + i) ;
                     i_ht = MIN(i_diag - i_row - 1,i - i_start + 1) ;
                     if( i_ht > 0) 
             {
                       *(stiff + k) -= dot_real((stiff +  k - i_ht),
                                       (stiff + i_diag - i_ht),
                                        i_ht) ;
             }
                     k++ ;
                    }
                }
           if( j_ht - 1 > 0 )
              {
               i_row = j_row + 1 ;
               i_end = j_diag - 1 ;
               k = j - j_diag ;
               for ( i = i_row ; i <= i_end ; i++)
                   {
                    i_diag = (int) *(diag + i + k)   ;
                    d =  *(stiff + i) ;
                    *(stiff + i) /= *(stiff + i_diag) ;
                    *(stiff + j_diag) -= ( d * *(stiff + i) ) ;          
                   }
              }

           j_row = j_diag ;
          }
     }

/************************** end of factor ************************************/ 


/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   full_back   ***************
*/

/*  INTERFACE:   */

    void
    full_back(
              real  *stiff ,         /* factored profile matrix */
              real  *rhs ,         /* assembled rhs vector */
              int   *diag ,           /* array of diagonal addresses */
              int   number_of_eqs     /* number of equations */
             )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine full_back
** full_back performs all three steps of backsubstitution.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      real d ;
      int j_row , j , j_diag , i_start , 
          j_ht  , i , k ;

 /*  forward reduction . . .  */
  
      j_row = (-1) ;
      for ( j = 0 ; j < number_of_eqs; j++ )
          {
           j_diag = (int) *(diag + j) ;
           j_ht   =  j_diag - j_row ;
           i_start = j - j_ht + 2  ;
           if( j_ht - 1 > 0 )
       {
             *(rhs + j) 
                -= dot_real((stiff + j_row + 1) 
                ,(rhs + i_start - 1) , ( j_ht - 1 ) ) ;
       }
           j_row = j_diag ;
          }

  /*  diagonal scaling . . .  */
  
      for( j = 0 ; j < number_of_eqs ; j++)
         {
           j_diag = (int) *(diag + j) ;
           *(rhs + j ) /= *( stiff + j_diag) ;
         }

  /*  backsubstitution . . .  */
  
      j = number_of_eqs - 1 ;
      j_diag = (int) *(diag + j) ;
      FOREVER     
            {
             d = *(rhs + j) ;
             j-- ;
             if( j < 0 ) break ;
             j_row = (int) *(diag + j) ;
             if ( j_diag - j_row  >  1 )
                {
                 i_start = j - j_diag + j_row + 2 ;
                 k = j_row - i_start + 1 ;
                 for( i = i_start ; i <= j ; i++ )
                    *(rhs + i) -= ( *(stiff + i + k) * d ) ;
                 }
             j_diag = j_row ;
            }
     } 
/************************** end of full_back *********************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   pcg_loop   ***************
*/

/*  INTERFACE:   */

     void
     pcg_loop(  )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine pcg_loop
** pcg_loop performs the preconditioned conjugate gradient iteration
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

/*
- Modifications as of June 2007 -
No longer using intial guesses for CG algorithm; instead use zero vectors.
Also, use absolute norm size from first time step as convergence goal during
subsequent steps, rather than relative reduction.
*/
  {
   real  *M , *r , *d , *t , *x , *guess ;
   real  alpha , beta , rMr ;
   int   neq , k ;
   
   neq = loc_sys.neq ;
   /* set pointers */
   guess = force.last_result ;
   t = pcg.temp ;
   M = pcg.precond ;
   r = pcg.r ;
   d = pcg.d ;
   x = force.full_rhs ;
   
   /* intialize loop */
   k = -1 ;
   move_real(guess,d,neq) ; /* temporarily use "d" to hold x(0) */
   clear_real(t,neq) ;
   elgrp_loop( ESPROD ) ;   /* A x(0)  (result in "t") */
 /*  globalize( mesh, t, 0 );   */ /* update matrix-vector product result */
   vadd(x,-1.0,t,r,neq) ;   /* r = b - Ax */
   vouter(M,r,d,neq) ;      /* d = M r */
   move_real(guess,x,neq) ; /* overwrite "rhs" array with current soln */
   
   /* enter CG loop */
   
   do
      {
       k++ ;
       vouter(M,r,t,neq) ;         /* M r */
       rMr = dot_real(r,t,neq) ;    /* r M r */
     /*  rMr = dot_par(mesh,r,t,loc_sys.neq_owned) ; */  /* r M r */
       if(k==0) pcg.norm = rMr ;   /* save the starting norm */
       if(k==0 && time_data.time==ZERO)
          {
           time_data.step_started = 0 ;

           /* quick return if the starting RHS is zero (numerically ok) */
           if(rMr == ZERO ) 
               { 
                attempt = PCG_ZERO_RHS ;   
                return ; 
               }
          }
        if(isnan(rMr))
          { printf("Bad norm!\n") ; attempt = PCG_ZERO_RHS ;   return ; }

       clear_real(t,neq) ;
       elgrp_loop( ESPROD ) ;      /* A d  (result in "t") */
    /*   globalize( mesh, t, 0 );   */    /* update matrix-vector product result */
       alpha = rMr / dot_real(d,t,neq) ;  /* alpha = (r M r)/(d A d) */
     /*  alpha = rMr / dot_par(mesh,d,t,loc_sys.neq_owned) ; */ /* alpha = (r M r)/(d A d) */
       vadd(x,alpha,d,x,neq) ;    /* x <- x + alpha d */
       vadd(r,-alpha,t,r,neq) ;   /* r <- r - alpha A d */
       vouter(M,r,t,neq) ;        /* M r */
       beta = dot_real(r,t,neq) / rMr ;   /* beta = (r M r)new / (r M r)old */
      /* beta = dot_par(mesh,r,t,loc_sys.neq_owned) / rMr ;  */ /* beta = (r M r)new / (r M r)old */
       vadd(t,beta,d,d,neq) ;     /* d <- M r + beta d */
       *(pcg.hist+(k%CGMAX)) = rMr * beta ;
      } while(converged(k) != 1) ;
      
   /* iterations loop finished... 
      reinitialize and return back to solver() */

 /* No starting guesses any more ; use the zero vector */
 /* gl,cdn: fix to address ringing issue in solver convergence algorithm */

      printf("PCG is converged; clearing for next step...(neq=%d)\n",neq) ;
      clear_real(guess,neq) ;
      printf("Clearing done.\n") ;
  }
  
/************************** end of pcg_loop **********************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   converged   ***************
*/

/*  INTERFACE:   */

  int converged(
                int   k   /* iteration number */
               )

/*  RETURN VALUE:  flag indicating successful convergence (or not) */
/*  DESCRIPTION:   */
/*
** function converged()
** monitors convergence of CG solver from magnitude and history of 
** the residual norm
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

  {
   real  ratio , mean1 , mean2 , meanall , stddev ;
   int  i , done ;
  
   if(time_data.step_started)   /* Use established absolute norm goal */
      {
       if( *(pcg.hist+(k%CGMAX)) <= pcg.first_norm )
          {
           done = 1 ;
          }
       else
          done = 0 ;
      }

   else          /* First step; seek to reduce norm by relative CGTOL factor */
      {
       if( (ratio = *(pcg.hist+(k%CGMAX))/pcg.norm) <= CGTOL )
          {
           done = 1 ;
           if(time_data.time > ZERO)
              {
               time_data.step_started = 1 ;
               pcg.first_norm = *(pcg.hist+(k%CGMAX)) ;
               /* we now use first_norm to store the absolute goal */
              }
          }
       else
          done = 0 ;
      }
  
  
  
   if( done )
      {
           fprintf(cgconv_file,"time=%g converged in %d iterations\n",
                 time_data.time,k) ;
           fprintf(cgconv_file,"starting norm=%g , ending norm=%g\n",
                 pcg.norm,*(pcg.hist+(k%CGMAX)) ) ;
           fflush(cgconv_file) ;
       return(1) ;  /* converged within tolerance */
      }
   else if( k < CGMAX || k%CGMAX != 0 )
      {
       return(0) ;  /* not converged, but insufficient history to declare it stalled */
      }
   else
      {
       mean1 = ZERO ;
       mean2 = ZERO ;
       stddev = ZERO ;
       
       for(i=0;i<(CGMAX/2);i++)
          {
           mean1 += *(pcg.hist+i)/(((real)(CGMAX/2))*pcg.norm) ;
           mean2 += *(pcg.hist+i+(CGMAX/2))/(((real)(CGMAX/2))*pcg.norm) ;
          }
       meanall = (mean1+mean2)/2.0 ;
       for(i=0;i<CGMAX;i++)
          {
           stddev +=
                 ((*(pcg.hist+i)/pcg.norm) - meanall)*((*(pcg.hist+i)/pcg.norm) - meanall)/((real)CGMAX) ;
          }
       stddev = sqrt(stddev) ;
       if( mean1 - mean2 < (stddev/((real)(CGMAX/2))) )   /* progress in CGMAX iterations is less than approx sdom */
          {
              
               printf("*****WARNING***** conjugate gradient STALLED at t=%g\n",time_data.time) ;
               fprintf(cgconv_file,"time=%g stalled in %d iterations; norm goal= %g\n",
                     time_data.time,k,pcg.first_norm) ;
               fprintf(cgconv_file,"starting=%g , ending=%g\n",
                     pcg.norm,*(pcg.hist+(k%CGMAX)) ) ;
               fflush(cgconv_file) ;
              
           return(1) ;  /* declare it stalled and bail out */
          }
       else
          {
           return(0) ;  /* still making progress; keep going */
          }
      }
  }
  
/************************** end of converged *********************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   put_soln   ***************
*/

/*  INTERFACE:   */

       void
       put_soln(
                RHS_DATA  *solve , /* structure containing rhs arrays */
                int       *id      /* id array of eq numbers */
               )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine put_soln
** put_soln transfers solution to nodal storage.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      int      eq_number ;
      COUNTER  node , j , ndof , pgnode;

      ndof = fe_sys.ndof ;
      
      /* march thru the external (unpermuted) nodes */
      for(node = 0 ; node < loc_sys.numnp ; node++)
            {
             pgnode = global.permut[node+1];

             for( j = 0 ; j < ndof ; j++ )
                     {
                      eq_number = id[(pgnode - 1)*ndof + j] ;

                      if(  eq_number != (-1) )
                            {
                             *(global.del_displ + node*ndof + j) =
                                  *(solve->full_rhs + eq_number)  ;
                            }

                      else if( time_data.elastic )
                            {
                             *(global.del_displ + node*ndof + j) =
                                  *(global.forv + node*ndof + j)  ;
                            }

                      else if( time_data.fail != 1 )
                            {
                             *(global.del_displ + node*ndof + j) =
                                  *(global.forv + node*ndof + j) *
                                  time_data.dt  ;
                            }
                     }
            }
     }
/************************** end of put_soln **********************************/ 



