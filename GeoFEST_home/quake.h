/*
***                           File quake.h                        ***
***                        GeoFEST version 6.0
*** Copyright (c) 2010, California Institute of Technology        ***
*** U.S.Sponsorship under NASA Contract NAS7-1407 is acknowledged ***
***
*** This software is designated for public release under JPL Task ***
*** Order Number NMO710991 and may be publicly released through
*** license with the Open Channel Foundation
***
*** This file contains modules which perform active fault slip in
*** response to evolving elastic stresses.
***
***          == failcheck ==
***          == fail_loop ==
***/

#include <stdio.h>
#include <stdlib.h>
#include "mydefs.h"
#include "errcodes.h"
#include "finel.h"



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   failcheck   ***************
*/

/*  INTERFACE:   */

     void
     failcheck(
               GROUP   *grp_ptr    /* pointer to current element group */
              );

/*  RETURN VALUE:  n/a */
/*  DESCRIPTION:   */
/*
** Routine failcheck
** failcheck looks at the elements comprising the specified split node strand
** and determines whether a failure condition has been met, returning a flag
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   fail_loop   ***************
*/

/*  INTERFACE:   */

     void
     fail_loop(
                int     f_index         /* 0-based fault strand identifier */
               );

/*  RETURN VALUE:  n/a */
/*  DESCRIPTION:   */
/*
** Routine fail_loop
** fail_loop is the topmost controlling routine that for a specific split node
** strand that has been declared 'failed', iteratively applies varying amounts of
** slip with the objective of minimizing global strain energy
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/


