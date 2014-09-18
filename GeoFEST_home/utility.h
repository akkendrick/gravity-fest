/*
***                           File utility.h                      ***
***                        GeoFEST version 6.0
*** Copyright (c) 2010, California Institute of Technology        ***
*** U.S.Sponsorship under NASA Contract NAS7-1407 is acknowledged ***
***
*** This software is designated for public release under JPL Task ***
*** Order Number NMO710991 and may be publicly released through
*** license with the Open Channel Foundation
***
*** This file contains miscellaneous utility routines used throughout
*** the finite element program: 
***
***          == move_real ==
***          == clear_real ==
***          == dot_real ==
***          == incr_real ==
***          == vadd ==
***          == vouter ==
***          == squawk ==
***          == real_data_revcompar
***          == split_local_node_compare
***/


#include <stdio.h>
#include "mydefs.h"
#include "errcodes.h"
#include "finel.h"



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   move_real   ***************
*/

/*  INTERFACE:   */

      void
      move_real(
                real   *from , /* source array */
                real   *to ,   /* destination array */
                int    n       /* number of entries */
               );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine move_real
** move_real copies data to a new location.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   clear_real   ***************
*/

/*  INTERFACE:   */

     void
     clear_real(
                real   *array , /* array to be zapped */
                int    n       /* number of entries */
               );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine clear_real
** clear_real nulls a data area.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   dot_real   ***************
*/

/*  INTERFACE:   */

  real dot_real(
                real   *vect_1 , /* multiplicand array */
                real   *vect_2 , /* multiplicand array */
                int    n       /* number of entries */
               );

/*  RETURN VALUE:  value of the scalar product */
/*  DESCRIPTION:   */
/*
** Function dot_real
** dot_real forms dot product of two vectors. 
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   incr_real   ***************
*/

/*  INTERFACE:   */

      void
      incr_real(
                real   *from , /* source array */
                real   *to ,   /* destination array */
                int    n       /* number of entries */
               );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine incr_real
** incr_real increments the contents of array "to" with the contents of "from".
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   vadd   ***************
*/

/*  INTERFACE:   */

           void
           vadd(
                real   *v1 , /* first array */
                real   mult , /* scalar multiplier for second array */
                real   *v2 , /* second array */
                real   *dest , /* destination array */
                int    neq       /* number of entries */
               );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine vadd
** calculates the linear combination of vectors
**    dest = v1 + mult*v2
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   vouter   ***************
*/

/*  INTERFACE:   */

         void
         vouter(
                real   *v1 , /* first array */
                real   *v2 , /* second array */
                real   *dest , /* destination array */
                int    neq       /* number of entries */
               );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine vouter
** calculates the outer product of two vectors
** actually a misnomer -- it's really the element-by-element product
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   squawk   ***************
*/

/*  INTERFACE:   */

         void
         squawk(
                char   *message  /* string to be sent to stdout */
               );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine squawk
** Prints a message to stdout that is common to all processors only on
** processor zero
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   real_data_revcompar***************
*/

/*  INTERFACE:   */

         int
         real_data_revcompar(
               const void* p,
               const void* q
               );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine real_data_revcompar
** Compares two reals for reverse-numeric qsort
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   split_local_node_compare***************
*/

/*  INTERFACE:   */

         int
         split_local_node_compare(
               const void* p,
               const void* q
               );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine split_local_node_compare
** Compares two node numbers for numeric qsort
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/


