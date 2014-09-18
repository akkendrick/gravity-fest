/*
***                           File grav.h                       ***
***                        GeoFEST version 6.0
*** Copyright (c) 2010, California Institute of Technology        ***
*** U.S.Sponsorship under NASA Contract NAS7-1407 is acknowledged ***
***
*** This software is designated for public release under JPL Task ***
*** Order Number NMO710991 and may be publicly released through
*** license with the Open Channel Foundation
***
*** This file contains modules which perform calculations 
*** related to element stresses and strains:
***
***          == grav_calc ==
***          == dgrav_form ==
***          == dgrav_flux ==
***/


#include <stdio.h>
#include <stdlib.h>
#include "mydefs.h"
#include "errcodes.h"
#include "finel.h"
#include <math.h>
#define EXTERN extern


/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   grav_calc   ***************
*/

/*  INTERFACE:   */

      void
      grav_calc(
                GROUP		*grp_ptr , /* pointer to current element group */
                int         task
               );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine grav_calc
** grav_calc computes and outputs gravity change information.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/
/************************** end of grav_calc *********************************/ 




/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   dgrav_form   ***************
*/

/*  INTERFACE:   */

      void
      dgrav_form(
                 real  *dgrav ,
                 ELEMENT_INFO  *info ,
                 int  side ,
                 real  dmass ,
                 real  *midpt ,
                 ELEMENT_DATA  *bel_pt
                );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine dgrav_form
** dgrav_form performs the part of the gravity change calculation resulting from
** the net Eulerian flow of mass in/out of each volume element, the results
** being evaluated at the position of each surface facet that is a buoyancy element.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/
/************************** end of dgrav_form *********************************/ 




/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   dgrav_flux   ***************
*/

/*  INTERFACE:   */

      void
      dgrav_flux(
                 real  *dgrav ,
                 ELEMENT_INFO  *info ,
                 int  rside ,
                 int  dside ,
                 real delta_rho ,
                 real *upvec ,
                 real big_g ,
                 ELEMENT_DATA  *del_pt ,
                 ELEMENT_DATA  *rel_pt ,
                 real *fmid
                );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine dgrav_flux
** dgrav_flux performs the part of the gravity change calculation resulting from
** the Eulerian mass flow imbalance occurring at the interface between
** solid and vaccuum or between contrasting densities, the results
** being evaluated at the position of each surface facet that is a buoyancy element.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/
/************************** end of dgrav_flux *********************************/ 


