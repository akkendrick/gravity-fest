/*
***                           File strain.h                       ***
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
***          == form_stress ==
***          == form_beta ==
***          == form_dbar ==
***          == duostress ==
***          == find_max_strain ==
***          == real_cmp ==
***/

#include <stdio.h>
#include "mydefs.h"
#include "errcodes.h"
#include "finel.h"
#include <math.h>


/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   form_stress   ***************
*/

/*  INTERFACE:   */
     void
     form_stress(
                 GROUP   *grp_ptr ,   /* pointer to current element group */
                 int     code         /* code specifying elastic or visco */
                );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine form_stress
** form_stress computes the stress or time-step stress increment
** based on the FE solution.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   form_beta   ***************
*/

/*  INTERFACE:   */
     void
     form_beta(
               ELEMENT_DATA  *el_pt  , /* ptr to current element struct */
               ELEMENT_MAT   *mat_pt , /* ptr to el material prop struct */
               ELEMENT_INFO  *info     /* ptr to element info struct */
              );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine form_beta
** form_beta computes beta, the viscoplastic strain rate for
** this element and time step.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   form_dbar   ***************
*/

/*  INTERFACE:   */
     void
     form_dbar(
               ELEMENT_DATA  *el_pt ,  /* ptr to current element struct */
               ELEMENT_MAT   *mat_pt , /* ptr to el material prop struct */
               ELEMENT_INFO  *info     /* ptr to element info struct */
              );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine form_dbar
** form_dbar computes dbar, the VE single-step 
** constitutive matrix:
**
** dbar === (S + alpha delta_t beta')^-1
**
** but we compute the factored form, not the explicit inverse.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   duostress   ***************
*/

/*  INTERFACE:   */
     void
     duostress(
               ELEMENT_DATA  *el_pt ,  /* ptr to current element struct */
               ELEMENT_MAT   *mat_pt , /* ptr to el material prop struct */
               ELEMENT_INFO  *info     /* ptr to element info struct */
              );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine duostress
** calculates stresses associated with 2-node truss or spring type elements
** (special case of 4-node tetrahedral elements)
**
** 1 component for trusses; 3 components for springs
**
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   find_max_strain   ***************
*/

/*  INTERFACE:   */
     void
     find_max_strain(
                     GROUP   *grp_ptr    /* pointer to current element group */
                    );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine find_max_strain
** find_max_strain sorts the viscoplastic strains in all elements
** and stores an adjusted maximum value for use in time step computation
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   real_cmp   ***************
*/

/*  INTERFACE:   */
     int
     real_cmp(
              const void *first_arg ,
              const void *second_arg
             );

/*  RETURN VALUE:  -1, 0, or 1 */
/*  DESCRIPTION:   */
/*
** Routine real_cmp
** comparison function to sort real numbers with qsort
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/


