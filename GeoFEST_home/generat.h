/*
***                           File generat.h                      ***
***                        GeoFEST version 6.0
*** Copyright (c) 2010, California Institute of Technology        ***
*** U.S.Sponsorship under NASA Contract NAS7-1407 is acknowledged ***
***
*** This software is designated for public release under JPL Task ***
*** Order Number NMO710991 and may be publicly released through
*** license with the Open Channel Foundation
***
***
*** This file contains modules for generating and reading 
*** node and element data and related isoparametric tasks 
***
***          == gen_element ==
***          == next_element ==
***          == read_surf ==
***          == read_buoy ==
***          == read_slip ==
***          == splitfind ==
***          == compute_fterms ==
***          == load_element ==
***          == gen_map ==
***          == gen_real ==
***          == dot_sh ==
***          == set_param_from_string ==
***
***/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "mydefs.h"
#include "errcodes.h"
#include "finel.h"
#include <math.h>


/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   gen_element   ***************
*/

/*  INTERFACE:   */

    void
    gen_element(
                GROUP   *grp_ptr    /* pointer to current element group */
               );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine gen_element
** gen_element generates element storage and ien array
** It reads in element group information, individual elements,
** and surface traction records.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   next_element   ***************
*/

/*  INTERFACE:   */

   int  next_element(
                     int   *node_pointer , /* pointer to hold node */
                     int   *gen_pointer  , /* pointer to gen parameter */
                     FILE  *readfile
                    );

/*  RETURN VALUE:  number of the next node or object read in */
/*  DESCRIPTION:   */
/*
** Function next_element
** next_element reads the first part of a record in a double 
** null-terminated list, so the calling function may 
** catch the termination.
**
** The use of this function is not limited to finite elements, but is
** also employed for other lists, such as nodes and tractions.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   read_surf   ***************
*/

/*  INTERFACE:   */

     void
     read_surf(
               int  *tally ,  /* counter for records read in */
               int  *temp_list ,  /* storage for surf location data */
               real  *temp_trac ,  /* storage for surf traction data */
               int  type ,    /* element type number */
               int  numel ,  /* number of elements */
               int  ndof  ,  /* number of disp degrees of freedom */
               int  *more    /* more time dependent tractions to read in */
              );

/*  RETURN VALUE:   -none- */
/*  DESCRIPTION:   */
/*
** Routine read_surf
** read_surf reads in surface traction records.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   read_buoy   ***************
*/

/*  INTERFACE:   */

     void
     read_buoy(
               int  *tally ,  /* counter for records read in */
               BUOY_DATA *buoy_ptr ,  /* buoyancy data structure */
               int  type ,    /* element type number */
               int  numel ,  /* number of elements */
               int  ndof    /* number of disp degrees of freedom */
              );

/*  RETURN VALUE:   -none- */
/*  DESCRIPTION:   */
/*
** Routine read_buoy
** read_buoy reads in buoyancy surface records.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   read_slip   ***************
*/

/*  INTERFACE:   */

     void
     read_slip(
               SPLITNODE *fspl, /* structure of node, group, vectors */
               int  type ,    /* element type number */
               int  *ptr_nsplit ,  /* local number of split nodes */
               int  g_nsplit ,  /* global number of split nodes */
               int  ndof    /* number of disp degrees of freedom */
              );

/*  RETURN VALUE:   -none- */
/*  DESCRIPTION:   */
/*
** Routine read_slip
** read_slip reads split-node records for symmetric fault 
** slip accumulation. By examining the global set of
** split nodes, it retains and counts the local split nodes.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   splitfind   ***************
*/

/*  INTERFACE:   */

 int splitfind(
               int   qnode ,  /* query node to search for */
               SPLITNODE *fspl ,  /* split nodes */
               int  nsplit    /* number of split nodes */
              );

/*  RETURN VALUE:   index number of found split node entry; -1 if not found */
/*  DESCRIPTION:   */
/*
** function splitfind
** splitfind performs an efficient bisection search on the list of split nodes 
** to see if qnode is among them.  If so, it returns the index number of
** that entry.  Otherwise, returns -1.  Assumes that the fspl list is
** provided sorted in ascending node number order!
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   compute_fterms   ***************
*/

/*  INTERFACE:   */

 void compute_fterms(
               SPLITNODE *fspl,   /* local split nodes */
               GROUP   *grp_ptr,  /* pointer to current element group */
               int nsplit,        /* number split nodes this processor */
               int g_nsplit       /* global number split nodes */
              );

/*  DESCRIPTION:   */
/*
** function compute_fterms
** finds elements containing split nodes, computes component of slip
** for side element is on, stores in group's slip list (val, grp).
** Also sorts local fspl nodes by local node number, stores fspl records
** in global.splitn for pyramid processing.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   load_element   ***************
*/

/*  INTERFACE:   */

  void
  load_element(
               ELEMENT_DATA  *el_pt  ,  /* ptr to current element struct */
               ELEMENT_INFO  *info ,   /* ptr to element info struct */
               int      nel_glob ,      /* one-based global element number */
               int      ien[] ,     /* holding for ien assignments */
               int      mat       /* one-based material index number */
              );

/*  RETURN VALUE:   -none- */
/*  DESCRIPTION:   */
/*
** Routine load_element 
** load_element loads element arrays with ien and equation data.
** The information for the element has been previously read into a 
** temporary array. 
**
** This function checks for degeneracies of 2D elements,
** right-hand rule order of tet nodes (swapping when necessary),
** and fills in the ien and lm arrays for the element.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   gen_map   ***************
*/

/*  INTERFACE:   */

    void
    gen_map();

/*  RETURN VALUE:   -none- */
/*  DESCRIPTION:   */
/*
** Routine gen_map
** gen_map fills the nodal activity map array.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   gen_real   ***************
*/

/*  INTERFACE:   */

    void
    gen_real(
             real      *array ,  /* array of generated values */
             int        dim   ,  /* dimensionality of grid */
             FILE      *the_file  /* input file to read from */
            );

/*  RETURN VALUE:   -none- */
/*  DESCRIPTION:   */
/*
** Routine gen_real
** gen_real reads or generates floating point global array.
** Used to read or generate nodes, displacement bcs, or 
** velocity bcs.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   dot_sh   ***************
*/

/*  INTERFACE:   */

    real  dot_sh(
                 real      *shape  ,     /* shape function storage */
                 real      *variable  ,  /* interpolate variable storage */
                 int        n            /* number to sum over */
                );

/*  RETURN VALUE:   dot product value */
/*  DESCRIPTION:   */
/*
** Function dot_sh
** dot_sh performs the dot product of the vector of shape functions
** (may be gradients of shape functions)
** with the vector "variable"
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
        set_eldata_from_string(
            char in_string[] ,
            GROUP   *grp_ptr    /* pointer to current element group */
                            );

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine set_param_from_string
** reads 1-line string (obtained from a file) for one parameter name value pair
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/


