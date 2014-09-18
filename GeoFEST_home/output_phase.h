#ifndef	GUARD_OUTPUT_PHASE_H
#define	GUARD_OUTPUT_PHASE_H 1
/*
***                        File output_phase.h                         ***
***                        GeoFEST version 4.8.4
*** Copyright (c) 2010, California Institute of Technology        ***
*** U.S.Sponsorship under NASA Contract NAS7-1407 is acknowledged ***
***
*** This software is designated for public release under JPL Task ***
*** Order Number NMO710991 and may be publicly released through
*** license with the Open Channel Foundation
***
*** This file contains modules for parallel assembly and output
*** of node and element based results.
***
***          == out_node_data_init==
***          == out_node_data_destroy==
***          == out_node_data_destroy_list==
***          == out_node_data_fprintf==
***          == out_node_data_compar==
***          == out_node_data_sort==
***          == out_node_data_fprintf_list==
***          == out_element_data_init==
***          == out_element_data_destroy==
***          == out_element_data_destroy_list==
***          == out_element_data_fprintf==
***          == out_element_data_compar==
***          == out_element_data_sort==
***          == out_element_data_fprintf_list==
*/

/*
***GeoFEST version 4.8.4***
Early pre-numbered modifications:
G. A. Lyzenga -- unified parallel/serial version and bug fixes 11/17/03
G. A. Lyzenga -- split node enhancements 8/16/02
G. A. Lyzenga -- modification to correctly relax hydrostatic stresses 6/15/00
J. W. Parker -- 3D support and deleted hypercube branches 9/16/98
G. A. Lyzenga  -- last modified on 3/03/95

Modification history:
2.1: fixes bug in boundary velocities and changes
     the way stresses are output.
2.2: implements split node faults
2.3: adds pass-thru comments
3.0: 3-D tetrahedra supported, hypercube branches deleted
3.1: hydrostatic relaxation changes
4.0: name change from visco to GeoFEST; correct implementation of 2-D out-of-
     plane stresses as well as hydrostatic relaxation; culling of remaining 
     vestiges of hypercube stuff; one code without recompilation for both
     2-D and 3-D domains
4.1: Introduction of trilinear hexahedral bricks and PCG iterative solver option
4.2: Bug fixes; standard source prologues added (June 02)
4.3: Corrects inefficiency in split node algorithm and introduces new
     split node input format
4.4: Adds dummy calls to MPI and Pyramid to provide a single unified
     source for parallel and serial codes.  Also introduces a major
     algorithmic enhancement; stress equilibrium maintenance during visco
     elastic stepping for long-term stability.  Also more minor bug fixes.
4.5: Supports surface tractions and roller/spring duo options for tet elements.
     Provides enhanced conjugate gradient convergence monitoring. 
4.6: Supports multiple split node fault strands/histories; buoyancy reaction forces
4.7: Includes simulation flow control codes
4.8: Supports nonplanar faults, optional and default control codes using
     controls.fem file, parallel support for boundary velocity conditions,
     better iteration convergence.
4.8.4: Support for these new functions:
       - Radial planetary geometry for buoyancy
       - Tensile mode disclacement on split node faults
       - More intelligent adaptive time stepping for nonlinear visco problems
       Critical bug fixes for:
       - Surface forces and buoyancy for tet elements of left-handed chirality
       - Surface forces and buoyancy for hexahedral element facets

Additional CVS-generated notes pertaining to 4.3 through 4.5:
(CVS was adapted under CT project and used for v4.3 to 4.5;
note CVS imposes its own version numbers with the comments)
 $Log: main.c,v $

 Revision 1.11  2005/05/02 23:38:23  jwparker
 Charles says he fixed an inadvertent modification that did more graphic
 related output but stopped early.  So this should be the working 4.5g.

 Revision 1.10  2005/04/25 23:39:26  jwparker
 This version is the Milestone G version which ran on Cosmos, and also the one
 we must post on OpenChannel for any verification. Further work should compare
 this with the previous revision, which ran on Los Angeles and so was
 further developed towards automatic thresholding for AMR.  There is also
 a beta version of 4.5 plus buoyancy and fault physics advances I think, yet
 to be checked in.

 This version may only be used on parallel computers due to AMR hard-wiring.

 Revision 1.8  2004/04/05 21:15:23  jwparker
 This is the update version 4.5 with several improvements over the
 milestone F code. These include single source file for parallel
 and sequential use (supported by Sequential_Library dummy calls
 to mimic MPI functions), physics based correction for stress
 equilibirium in long-term viscous solutions,enhanced convergence
 monitoring, and support for surface tractions and roller/spring duo
 options for tet elements.

 Revision 1.4  2003/02/09 06:10:55  lyzenga
 *** empty log message ***

END  of version history
*/

#include <stdio.h>
#include <stdlib.h>

#define INLINE inline
#define	SPATIAL_DIMENSIONS 3
#define	DEGREES_OF_FREEDOM 3

typedef	struct Coordinates {
  double element[SPATIAL_DIMENSIONS];
  } Coordinates;

typedef	struct Displacement {
  double element[DEGREES_OF_FREEDOM];
  } Displacement;

typedef	struct deltaDisplacement {
  double element[DEGREES_OF_FREEDOM];
  } deltaDisplacement;

typedef struct out_node_data {
  size_t		node;
  size_t		group;
  Coordinates		coords;
  Displacement		displ;
  deltaDisplacement	del_displ;
  } OUT_NODE_DATA;

INLINE static
void out_node_data_init(
    OUT_NODE_DATA*	p,
    size_t		node,
    size_t		group,
    const double	coords[],
    const double	displ[],
    const double	del_displ[]
    ) {
  size_t  i ;
  
  p->node	= node;
  p->group	= group;
  for(i=0;i<SPATIAL_DIMENSIONS;i++)
     {
      if(i >= fe_sys.nsd)  p->coords.element[i]=0.0 ;
      else  p->coords.element[i]=coords[i] ;
     }
  for(i=0;i<DEGREES_OF_FREEDOM;i++)
     {
      if(i >= fe_sys.ndof)  {p->displ.element[i]=0.0 ; p->del_displ.element[i]=0.0 ;}
      else  {p->displ.element[i]=displ[i] ; p->del_displ.element[i]=del_displ[i] ;}
     }
  }

INLINE static
void out_node_data_destroy(
    const
    OUT_NODE_DATA*	p
    ) {
  }

INLINE static
void out_node_data_destroy_list(
    const
    OUT_NODE_DATA	p[],
    size_t              records
    ) {
  }

INLINE static
int out_node_data_fprintf(
    FILE*	stream,
    const	OUT_NODE_DATA* p
  ) {
  int characters = fprintf(stream, "node %d", (int)p->node);
  { const
    size_t	n = fe_sys.nsd;
    size_t	j = 0;
    for (j = 0; j < n; ++j)
      characters += fprintf(stream, "   %f", p->coords.element[j]);
    }
  { const
    size_t	n = fe_sys.ndof;
    size_t	j = 0;
    for (j = 0; j < n; ++j)
      characters += fprintf(stream, "   %g", p->displ.element[j]);
    }
  { const
    size_t	n = fe_sys.ndof;
    size_t	j = 0;
    for (j = 0; j < n; ++j)
      characters += fprintf(stream, "   %g", p->del_displ.element[j]);
    }
  characters += fprintf(stream, "\n");
  return characters;
  }

INLINE static
int out_node_data_compar(const void* p, const void* q) {
  return (*((const OUT_NODE_DATA**)p))->node
       - (*((const OUT_NODE_DATA**)q))->node;
  }

INLINE static const
OUT_NODE_DATA** out_node_data_sort(
    /* This function allocates and returns an array of pointers
     * which must be free'd by the calling program.			*/
    const
    OUT_NODE_DATA	q[],
    size_t		records
    ) {
  const
  OUT_NODE_DATA**	p = calloc(records, sizeof(OUT_NODE_DATA*));
  size_t		record = 0;
  for (record = 0; record < records; ++record)
    p[record] = &q[record];
  qsort((void*)p, records, sizeof(OUT_NODE_DATA*), out_node_data_compar);
  return p;
  }

INLINE static
int out_node_data_fprintf_list(
    FILE*		stream,
    const
    OUT_NODE_DATA	q[],
    size_t		records
    ) {
  const
  OUT_NODE_DATA**	p = out_node_data_sort(q, records);
  int			characters = 0;
  size_t		record = 0;
  for (record = 0; record < records; ++record)
    characters += out_node_data_fprintf(stream, p[record]);
  free((void*)p);
  return characters;
  }

typedef	struct Stress {
  double element[2*DEGREES_OF_FREEDOM];
  } Stress;

typedef struct out_element_data {
  size_t		element;
  size_t		group;
  Coordinates		coords;
  Stress		stress;
  } OUT_ELEMENT_DATA;

INLINE static
void out_element_data_init(
    OUT_ELEMENT_DATA*	p,
    size_t		element,
    size_t		group,
    size_t      str_size,
    const double	coords[],
    const double	stress[]
    ) {
  size_t  i ;

  p->element	= element;
  p->group	= group;
  for(i=0;i<SPATIAL_DIMENSIONS;i++)
     {
      if(i >= fe_sys.nsd)  p->coords.element[i]=0.0 ;
      else  p->coords.element[i]=coords[i] ;
     }
  for(i=0;i<2*DEGREES_OF_FREEDOM;i++)
     {
      if(i >= str_size)  p->stress.element[i]=0.0 ;
      else  p->stress.element[i]=stress[i] ;
     }
  }

INLINE static
void out_element_data_destroy(
    const
    OUT_ELEMENT_DATA*	p
    ) {
  }

INLINE static
void out_element_data_destroy_list(
    const
    OUT_ELEMENT_DATA	p[],
    size_t              records
    ) {
  }

INLINE static
int out_element_data_fprintf(
    FILE*	stream,
    const	OUT_ELEMENT_DATA* p,
    size_t  str_size
  ) {
  int characters = 0;
  { const
    size_t	n = fe_sys.nsd;
    size_t	j = 0;
    for (j = 0; j < n; ++j)
      characters += fprintf(stream, "   %f", p->coords.element[j]);
    }
  { const
    size_t	n = str_size;
    size_t	j = 0;
    for (j = 0; j < n; ++j)
      characters += fprintf(stream, "   %f", p->stress.element[j]);
    }
  characters += fprintf(stream, "   %d\n", (int)p->element);
  return characters;
  }

INLINE static
int out_element_data_compar(const void* p, const void* q) {
  /* This routine must be modified to sort on group number.		*/
  return (*((const OUT_ELEMENT_DATA**)p))->element
       - (*((const OUT_ELEMENT_DATA**)q))->element;
  }

INLINE static const
OUT_ELEMENT_DATA** out_element_data_sort(
    /* This function allocates and returns an array of pointers
     * which must be free'd by the calling program.			*/
    const
    OUT_ELEMENT_DATA	q[],
    size_t		records
    ) {
  const
  OUT_ELEMENT_DATA**	p = calloc(records, sizeof(OUT_ELEMENT_DATA*));
  size_t		record = 0;
  for (record = 0; record < records; ++record)
    p[record] = &q[record];
  qsort((void*)p, records, sizeof(OUT_ELEMENT_DATA*), out_element_data_compar);
  return p;
  }

INLINE static
int out_element_data_fprintf_list(
    FILE*		stream,
    const
    OUT_ELEMENT_DATA	q[],
    size_t		records,
    size_t      str_size
    ) {
  const
  OUT_ELEMENT_DATA**	p = out_element_data_sort(q, records);
  int			characters = 0;
  size_t		record = 0;
  for (record = 0; record < records; ++record)
    characters += out_element_data_fprintf(stream, p[record], str_size);
  free((void*)p);
  return characters;
  }

#undef	DEGREES_OF_FREEDOM
#undef	SPATIAL_DIMENSIONS
#undef	INLINE
#endif/*GUARD_OUTPUT_PHASE_H 1	*/

