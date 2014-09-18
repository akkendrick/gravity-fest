/*
***                         File solver.h                         ***
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

#include <stdio.h>
#include <stdlib.h>
#include "mydefs.h"
#include "errcodes.h"
#include "finel.h"
#include "main.h"
#include <math.h>



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   solver   ***************
*/

/*  INTERFACE:   */

     void
     solver(
            int  code     /* which solver option to use */
           );

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
           );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine addstiff
** addstiff adds element stiffness to global profile array,
** or to element storage if using PCG solver
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



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
           );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine addfor
** addfor assembles the global r.h.s. vector.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   estiffprod   ***************
*/

/*  INTERFACE:   */

 void
 estiffprod(
            GROUP   *grp_ptr    /* pointer to current element group */
           );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine estiffprod
** calculates the product of a given vector with stiffness in element storage
**    t = A * d
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   reorder   ***************
*/

/*  INTERFACE:   */

    void
    reorder(
            GROUP   *grp_ptr    /* pointer to current element group */
           );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine reorder
** reorder uses ien information to build adjacency information and call
** permutation optimizing routines for minimizing matrix profile.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



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
           );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine genrcm
** genrcm computes the reverse cuthill mckee ordering for a general
** adjacency graph such as a sparse finite element matrix. 
*/
/*---------------------------------------------------------------------------*/



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
           );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** rcm does connectivity analysis (jwp interpretation; see below).
*/
/*---------------------------------------------------------------------------*/



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
           );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine degree
** degree computes the graph steps to each node in the subgraph. 
*/
/*---------------------------------------------------------------------------*/



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
           );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine fnroot
** fnroot finds the pseudo-peripheral node for a given subgraph.
*/
/*---------------------------------------------------------------------------*/



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
           );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine rootls
** rootls generates the level structure corresponding to "root".
*/
/*---------------------------------------------------------------------------*/



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
           );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine colht
** colht computes column heights of global array.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



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
             );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine profile_diag
** profile diag computes the diagonal addresses.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



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
             );

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
             );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine full_back
** full_back performs all three steps of backsubstitution.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   pcg_loop   ***************
*/

/*  INTERFACE:   */

     void
     pcg_loop(  );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine pcg_loop
** pcg_loop performs the preconditioned conjugate gradient iteration
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   converged   ***************
*/

/*  INTERFACE:   */

  int converged(
                int   k   /* iteration number */
               );

/*  RETURN VALUE:  flag indicating successful convergence (or not) */
/*  DESCRIPTION:   */
/*
** function converged()
** monitors convergence of CG solver from magnitude and history of 
** the residual norm
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



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
               );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine put_soln
** put_soln transfers solution to nodal storage.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/


