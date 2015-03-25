/*
***                           File inphase.h                      ***
***                        GeoFEST version 6.0
*** Copyright (c) 2010, California Institute of Technology        ***
*** U.S.Sponsorship under NASA Contract NAS7-1407 is acknowledged ***
***
*** This software is designated for public release under JPL Task ***
*** Order Number NMO710991 and may be publicly released through
*** license with the Open Channel Foundation
***
*** This file contains subroutines which handle input, 
*** equation numbering, and memory allocation.
***
***          == input_phase ==
***          == matrix_alloc ==
***          == gen_number ==
***          == output_phase ==
***          == el_output ==
***          == locate_pt ==
***          == set_flow_from_string ==
***          == set_param_from_string ==
***/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "mydefs.h"
#include "errcodes.h"
#include "finel.h"
#include "output_phase.h"



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   input_phase   ***************
*/

/*  INTERFACE:   */

void
input_phase();

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine input_phase
** input_phase is the top-level input file reading function.  
** Most input-file records are read here or in routines 
** in generat.c that are called by this function.
** A few items are read by main.c as well.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   matrix_alloc   ***************
*/

/*  INTERFACE:   */

void
matrix_alloc();

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine matrix_alloc
** allocates space for the sparse stiffness matrix _or_
** the PCG solver arrays
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   gen_number   ***************
*/

/*  INTERFACE:   */

void
gen_number();

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine gen_number
** gen_number assigns equation numbers based on nodes and ndof
** Note: Constructed inverse id_pointer array to find global
**       1-based node numbers from equation numbers (in dot_par)
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   output_phase   ***************
*/

/*  INTERFACE:   */

   void
   output_phase(
                int    quake,  /* flag to indicate quake or not */
                int    redo    /* flag to repeat prior timestep */
               );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine output_phase
** output_phase writes out requested data at a quake 
** or scheduled time.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   el_output   ***************
*/

/*  INTERFACE:   */

      void
      el_output(
                GROUP		*grp_ptr/* pointer to current element group */
               );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine el_output
** el_output prints out element stress information.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   locate_pt   ***************
*/

/*  INTERFACE:   */

      void
      locate_pt(
                ELEMENT_DATA  *el_pt  ,  /* ptr to current element struct */
                ELEMENT_INFO  *info ,   /* ptr to element info struct */
                real   xpt[3][10]     /* return array of computed coords */
               );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine locate_pt
** locate_pt maps the integration points into physical space.
** Global coordinates are placed in xpt[j]
** Called by el_output element stress output function (above).
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   set_flow_from_string   ***************
*/

/*  INTERFACE:   */

      void
        set_flow_from_string(
            char in_string[]
                            );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine set_flow_from_string
** reads 1-line string (obtained from a file) for one flow_control name value pair
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   set_param_from_string   ***************
*/

/*  INTERFACE:   */

      void
        set_param_from_string(
            char in_string[]
                            );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine set_param_from_string
** reads 1-line string (obtained from a file) for one parameter name value pair
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/


