/*
***                           File stiff.h                        ***
***                        GeoFEST version 6.0
*** Copyright (c) 2010, California Institute of Technology        ***
*** U.S.Sponsorship under NASA Contract NAS7-1407 is acknowledged ***
***
*** This software is designated for public release under JPL Task ***
*** Order Number NMO710991 and may be publicly released through
*** license with the Open Channel Foundation
***
*** This file contains modules which construct the finite
*** element stiffness matrix and the right-hand-side "force" vector.
***
***          == form_stiff ==
***          == form_rhs ==
***          == form_bc ==
***          == form_slip ==
***          == adjust_bc ==
***          == lame_form ==
***          == duo_form ==
***          == force_form ==
***          == surf_form ==
***          == buoy_form ==
***          == shape ==
***          == dotsh ==
***          == adfldp ==
***          == p_shape ==
***          == form_equil ==
***/

#include <stdio.h>
#include <stdlib.h>
#include "mydefs.h"
#include "errcodes.h"
#include "finel.h"
#include <math.h>

/* Gaussian quadrature points (without scale factors) and weights */

#define WT1  0.78867513
#define WT2  0.21132487
#define RTET1 0.138196601125011
#define RTET2 0.585410196624969

/* 1/24 */
#define WTET 0.041666666666667
/* 4/3 */
#define WHEX 1.333333333333333


/* Arrays of integration constants common to stiff.c, strain.c */
static real
   x_2x2[4] = {-ONE  , ONE  , ONE  , -ONE } ,
   y_2x2[4] = {-ONE  ,-ONE  , ONE  ,  ONE } ,
   wt_2x2[4] = { ONE  , ONE  , ONE  ,  ONE} ,

   x_2x2x2[8] = {-ONE  , ONE  , ONE  , -ONE , -ONE  , ONE  , ONE , -ONE } ,
   y_2x2x2[8] = {-ONE  , -ONE  , ONE  , ONE , -ONE  , -ONE  , ONE , ONE } ,
   z_2x2x2[8] = {-ONE  , -ONE  , -ONE  , -ONE , ONE  , ONE  , ONE , ONE } ,
   wt_2x2x2[8] = {ONE  , ONE  , ONE  , ONE , ONE  , ONE  , ONE , ONE } ,

   x_3x3[9] =
    {-ONE  , ONE  , ONE  , -ONE , ZERO  ,  ONE , ZERO  , -ONE , ZERO} ,
   y_3x3[9] =
    {-ONE  ,-ONE  , ONE  ,  ONE , -ONE  , ZERO ,  ONE  , ZERO , ZERO} ,
   wt_3x3[9] =
    { 25.0/81.0 , 25.0/81.0 , 25.0/81.0 , 25.0/81.0 , 40.0/81.0 ,
      40.0/81.0 , 40.0/81.0 , 40.0/81.0 , 64.0/81.0               } ,

   x_tet[4] = {RTET1, RTET2, RTET1, RTET1} ,
   y_tet[4] = {RTET1, RTET1, RTET2, RTET1} ,
   z_tet[4] = {RTET1, RTET1, RTET1, RTET2} ,
   wt_tet[4] = {WTET,WTET,WTET,WTET} ,
   
   r_3edge[3] = {PT5  , PT5  , ZERO } ,
   s_3edge[3] = {ZERO  , PT5  , PT5 } ,
   t_3edge[3] = {PT5  , ZERO  , PT5 } ,
   wt_3edge[3] = { ONE/SIX  , ONE/SIX  , ONE/SIX } ;


static real
   x_tri[6] = {-ONE, ZERO, ZERO, ONE, ZERO, ZERO} ,
   y_tri[6] = {ZERO, -ONE, ZERO, ZERO, ONE, ZERO} ,
   z_tri[6] = {ZERO, ZERO, -ONE, ZERO, ZERO, ONE} ,
   wt_tri[6] = {WHEX,WHEX,WHEX,WHEX,WHEX,WHEX} ,
   ra[8] = {-0.5, 0.5, 0.5,-0.5, -0.5, 0.5,0.5,-0.5} ,
   sa[8] = {-0.5,-0.5, 0.5, 0.5, -0.5,-0.5,0.5, 0.5} ,
   ta[8] = {-0.5,-0.5,-0.5,-0.5,  0.5, 0.5,0.5, 0.5} ;



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   form_stiff   ***************
*/

/*  INTERFACE:   */

     void
     form_stiff(
                GROUP   *grp_ptr ,   /* pointer to current element group */
                int     code         /* code specifying elastic or visco */
               );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine form_stiff
** form_stiff computes the element-wise stiffness array for the elt group.
** "code" indicates if this is for an elastic problem or a VE step.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   form_rhs   ***************
*/

/*  INTERFACE:   */

     void
     form_rhs(
              GROUP   *grp_ptr ,   /* pointer to current element group */
              int     code         /* code specifying elastic or visco */
             );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine form_rhs
** form_rhs computes the right-hand-side vector for the FE problem.
** "code" indicates if this is for the elastic problem or a VE step.
** There are four stages to form_rhs, doing set-up and invocation of::
**    force_form (computes body force and VE strain-related force)
**    surf_form (computes rhs term for surface tractions)
**    buoy_form (computes rhs term for buoyancy at density contrast horizons)
**    form_slip (computes rhs term for split nodes)
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   form_bc   ***************
*/

/*  INTERFACE:   */

     void
     form_bc(
             GROUP   *grp_ptr ,   /* pointer to current element group */
             int     code         /* code specifying elastic or visco */
            );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine form_bc
** form_bc computes the boundary condition terms for stiffness and rhs
** due to imposed conditions
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   form_slip   ***************
*/

/*  INTERFACE:   */

     void
     form_slip(
               int  code ,     /* code specifying elastic or visco */
               int  cslip ,     /* code specifying continuous slip */
               real  *rhs ,   /* temp array to hold element rhs contrib */
               real  *stiff , /* temp array to hold el stiffness contrib */
               ELEMENT_INFO  *info ,   /* ptr to element info struct */
               int  nel ,       /* current element number */
               int  node ,      /* current node number */
               real  *slip ,   /* array containing fault slip amplitudes */
               real frac ,     /* fraction of total slip amount for quakes  */
               ELEMENT_DATA  *el_pt ,  /* ptr to current element struct */
               ELEMENT_MAT   *mat_pt  /* ptr to el material prop struct */
              );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine form_slip
** form_slip computes the split-nodes fault offsets and their influence 
** on the finite element stiffness and RHS.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   adjust_bc   ***************
*/

/*  INTERFACE:   */

     void
     adjust_bc(
               real  *disp ,   /* array containing fault slip amplitudes */
               real  *rhs ,   /* temp array to hold element rhs contrib */
               real  *es ,    /* temp array to hold el stiffness contrib */
               int   nee        /* number of element equations */
              );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine adjust_bc
** adjust_bc computes the necessary terms that modify the right-hand-side
** due to imposed displacements.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   lame_form   ***************
*/

/*  INTERFACE:   */

     void
     lame_form(
               int  code ,     /* code specifying elastic or visco */
               real  *es ,    /* temp array to hold el stiffness contrib */
               ELEMENT_INFO  *info ,   /* ptr to element info struct */
               ELEMENT_MAT   *mat_pt ,  /* ptr to el material prop struct */
               ELEMENT_DATA  *el_pt    /* ptr to current element struct */
              );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine lame_form
** lame_form computes the volumetric stiffness term for one element,
** that is based on the constitutive volume constants. 
** Results are stored in the element stiffness array.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   duo_form   ***************
*/

/*  INTERFACE:   */

     void
     duo_form(
               int  code ,     /* code specifying elastic or visco */
               real  *es ,    /* temp array to hold el stiffness contrib */
               ELEMENT_INFO  *info ,   /* ptr to element info struct */
               ELEMENT_MAT   *mat_pt ,  /* ptr to el material prop struct */
               ELEMENT_DATA  *el_pt    /* ptr to current element struct */
              );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine duo_form
** duo_form computes the stiffness term for a two-node element,
** with options for either a truss "roller" or an isotropic spring. 
** Results are stored in the element stiffness array.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   force_form   ***************
*/

/*  INTERFACE:   */

    void
    force_form(
               int  code ,     /* code specifying elastic or visco */
               real  *rhs ,   /* temp array to hold element rhs contrib */
               ELEMENT_INFO  *info ,   /* ptr to element info struct */
               ELEMENT_MAT   *mat_pt ,  /* ptr to el material prop struct */
               ELEMENT_DATA  *el_pt    /* ptr to current element struct */
              );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine force_form
** force_form computes the "dumb" gravity contribution to the element
** right-hand side term for the elastic case, and the strain-rate
** term for the VE case.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   surf_form   ***************
*/

/*  INTERFACE:   */

     void
     surf_form(
               int  code ,     /* code specifying elastic or visco */
               real  *rhs ,   /* temp array to hold element rhs contrib */
               ELEMENT_INFO  *info ,   /* ptr to element info struct */
               int  side ,    /* index number of the element side */
               real  *trac ,   /* array containg traction components */
               ELEMENT_DATA  *el_pt    /* ptr to current element struct */
              );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine surf_form
** surf_form computes the surface traction forcing term for the 
** right-hand-side in the finite element problem.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   buoy_form   ***************
*/

/*  INTERFACE:   */

     void
     buoy_form(
               int  code ,     /* code specifying radial or rectilinear */
               real  *rhs ,   /* temp array to hold element rhs contrib */
               ELEMENT_INFO  *info ,   /* ptr to element info struct */
               int  side ,    /* index number of the element side */
               real  *up_vec ,   /* array containing local "up" */
               real  rho_g ,   /* array containing density*g contrast */
               ELEMENT_DATA  *el_pt    /* ptr to current element struct */
              );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine buoy_form
** buoy_form computes the buoyancy traction term due to displacement
** of a density contrast in the finite element grid.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   shape   ***************
*/

/*  INTERFACE:   */

     void
     shape(
           ELEMENT_DATA  *el_pt  ,  /* ptr to current element struct */
           ELEMENT_INFO  *info ,   /* ptr to element info struct */
           real  *sh_pt ,   /* el storage to hold shape functions */
           real  *det_pt ,   /* el storage to hold determinant */
           int  type ,    /* element type number */
           int  nen ,    /* number of element nodes */
           int  nint ,   /* number of Gauss integration points */
           int  all     /* flag indicating we need derivatives also */
          );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine shape
** shape computes the shape functions and gradients for a single element.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   dotsh   ***************
*/

/*  INTERFACE:   */

real dotsh(
           real  *array ,  /* array of nodal quantities to interpolate */
           int  i ,     /* derivative component index number */
           int  j ,     /* dof index number of quantity component */
           real  *sh ,   /* el storage containing shape functions */
           ELEMENT_DATA  *el_pt  ,  /* ptr to current element struct */
           ELEMENT_INFO  *info ,   /* ptr to element info struct */
           int  nen ,    /* number of element nodes */
           int  ndim ,    /* number of dimensions in array */
           int  st_flag  ,    /* flag to invoke split node adjustment */
           int  task_code   /* elastic, stepping or quake  */
          );

/*  RETURN VALUE:  the element-interpolated sum value */
/*  DESCRIPTION:   */
/*
** Function dotsh
** dotsh computes the dot product of a set of shape-functions (or gradients)
** with a vector of values, hence computing an interpolated coordinate
** or function. If st_flag, x is assumed to represent displacement, and
** a correction for split-nodes is applied prior to the dot product.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   adfldp   ***************
*/

/*  INTERFACE:   */

real adfldp(
            int  j ,     /* dof index number of slip component */
            ELEMENT_DATA  *el_pt  ,  /* ptr to current element struct */
            int  node ,    /* number of current node */
            ELEMENT_INFO  *info  ,  /* ptr to element info struct */
            int  task_code   /* elastic, stepping or quake  */
           );

/*  RETURN VALUE:  the required additional nodal displacement */
/*  DESCRIPTION:   */
/*
** Function adfldp
** adfldp returns to dotsh() the added nodal displacements
** needed to account for split node fault slip in calculating strain
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   p_shape   ***************
*/

/*  INTERFACE:   */

     void
     p_shape(
             GROUP   *grp_ptr     /* pointer to current element group */
            );

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine p_shape
** p_shape computes the parent-space shape function for common element types.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   form_equil   ***************
*/

/*  INTERFACE:   */
     void
     form_equil(
                GROUP   *grp_ptr    /* pointer to current element group */
               );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine form_equil
** computes corrections necessary to preserve equilibrium 
** in stress increments calculated during viscous stepping
**
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/


