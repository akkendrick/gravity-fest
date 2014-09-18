/*
***                           File main.h                         ***
***                        GeoFEST version 6.0
*** Copyright (c) 2010, California Institute of Technology        ***
*** U.S.Sponsorship under NASA Contract NAS7-1407 is acknowledged ***
***
*** This software is designated for public release under JPL Task ***
*** Order Number NMO710991 and may be publicly released through
*** license with the Open Channel Foundation
***
*** This file contains the program main entry point, the main
*** task driver, and modules for driving high-level functions
*** and interactions with the operator:
***
***          == main == ++
***          == elastic == ++
***          == time_step == ++
***          == clear_stiff == ++
***          == elgrp_loop == ++
***          == node_load == ++
***          == accumulate == ++
***          == completion == ++
***          == wrt_save == ++
***          == vrestart == ++
***          == el_save == ++
***          == check_monitor == ++
***          == terminate == ++
***/

#include <stdio.h>
#include <stdlib.h>
#include "mydefs.h"
#include "finel.h"
#include "errcodes.h"
#include <math.h>
#include <time.h>


/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   elastic   ***************
*/

/*  INTERFACE:   */

void
elastic(int slip_all);

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine elastic
** elastic performs an elastic solution of the 
** finite element problem.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   time_step   ***************
*/

/*  INTERFACE:   */

void
time_step();

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine time_step
** time_step computes the visco-elastic time stepping solutions
** to the finite element problem.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   clear_stiff   ***************
*/

/*  INTERFACE:   */

void
clear_stiff(void);

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine clear_stiff
** clear_stiff nulls the stiffness matrix memory.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   elgrp_loop   ***************
*/

/*  INTERFACE:   */

void
     elgrp_loop(
                int  task      /* which element-level task to perform */
               );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine elgrp_loop
** elgrp_loop invokes a particular task that relies on 
** element structures and affects all the elements in 
** every element group.
*/
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   node_load   ***************
*/

/*  INTERFACE:   */

void
     node_load(
               real  time ,   /* current simulation time */
               int reform     /* flag indicating stiffness reform status */
              );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine node_load
** node_load applies node-based forces that may accumulate
** with time.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   accumulate   ***************
*/

/*  INTERFACE:   */

void
accumulate(void);

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine accumulate
** accumulate updates the total displacement based on the 
** single-step increment.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   completion   ***************
*/

/*  INTERFACE:   */

void
completion(void);

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine completion
** completion syncronizes multiple processors. 
** Keeps track of incurred errors in each processor
** Just a dummy call in sequential code
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   wrt_save   ***************
*/

/*  INTERFACE:   */
     void
     wrt_save();


/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine wrt_save
** wrt_save dumps state to disk for restarts.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   vrestart   ***************
*/

/*  INTERFACE:   */

void
vrestart();

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** vrestart restores state from previous dump and 
** continues simulation.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   el_save   ***************
*/

/*  INTERFACE:   */
     void
     el_save(
             GROUP   *grp_ptr ,       /* pointer to current element group */
             int     code             /* code specifying which action */
            );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine el_save
** el_save dumps or restores element data.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   check_monitor   ***************
*/

/*  INTERFACE:   */
     void
     check_monitor(
                   int   nstep    /* count of steps for backup */
                  );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine check_monitor
** check_monitor allows user run intervention via strings stored
** in the file monitor.fem.
*/
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*                                                                                                                                  
*************   ROUTINE:   terminate   ***************                                                                          
*/

/*  INTERFACE:   */

void
terminate();

 /*  RETURN VALUE:  -none- */
 /*  DESCRIPTION:   */
 /*                                                                                                                                  
 ** Routine terminate                                                                                                               
 ** Normal exit of the program from simulation flow                                                                        
 **/
 /*   EOP   */
 /*---------------------------------------------------------------------------*/
 
 
