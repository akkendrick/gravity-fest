#ifndef	GUARD_FINEL_H
#define	GUARD_FINEL_H 1
/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   INCLUDE:   finel.h   ***************
*/

/*  DESCRIPTION:
     global variables and definitions for GeoFEST finite element program
*/
/*   EOP   */
/*---------------------------------------------------------------------------*/
/*
***                           File finel.h                        ***
***                        GeoFEST version 6.0
*** Copyright (c) 2010, California Institute of Technology        ***
*** U.S.Sponsorship under NASA Contract NAS7-1407 is acknowledged ***
***
*** This software is designated for public release under JPL Task ***
*** Order Number NMO710991 and may be publicly released through
*** license with the Open Channel Foundation
***/


/* global variables and definitions for finite element program */

/*
#include<pamr_mesh.h>
*/
/*
#define REVNO "4.8"
*/

/* finite element standard defs for 2-d  or 3-d code */

/* define system types                  */
#define  INTERIOR     0
#define  SHARED       1

/* definitions specific to this version */
#define MAX_STRING_LENGTH 250
#define MAXEL   5000
#define MAX_GROUP  1
#define MAX_PROCS 32
#define MAX_NEIGHBOR 5
/* generous a priori bound on how many nodes connect to a node, on average */
/* to be used to size adjacency array icol(number of nonzeros). */
#define MAX_CONNECT 50

/* macro definitions */

#define GLOBALIZE \
        ((solve->neq_all == fe_sys.neq) ? (GLOBALIZE_F) : (GLOBALIZE_S))
#define ELASTIC   0
#define VISCO     1
#define QUAKEEVT     2
#define FACBACK   0
#define BACK      1
#define ITER      2
#define SUB       0
#define FULL      1
#define DIRECT      1
#define PCG      2
#define BILIN         1
#define SERENDIP      2
#define QUAKE         3
#define TET           4
#define TRUSS         5
#define TRILIN        6
#define GENERATE     1
#define OUTPUT       2
#define SHAPE        3
#define FORMS_ELAS   4
#define FORMS_STEP   5
#define RHS_ELAS     6
#define RHS_STEP     7
#define BC_ELAS      8
#define BC_STEP      9
#define COL_HT      10
#define E_STRESS    11
#define V_STRESS    12
#define DUMP        13
#define RESTORE     14
#define NEW_LM      15
#define FAIL_CHECK  16
#define SMOOTH_BEGIN    17
#define SMOOTH_END      18
#define REORDER         19
#define ESPROD         20
#define EQUIL         21
#define MAX_STRAIN     22
#define FORMS_QUAKE     23
#define RHS_QUAKE     24
#define Q_STRESS     25
#define GRAV_CALC     26
#define SETUP            1
#define DRYRUN           2
#define GLOBALIZE_S      3
#define GLOBALIZE_F      4
#define DOT_PROD         5
#define FOLD2            6
#define FOLD4            7
#define FOLD6            8
#define FOLD8            9
#define FOLDI            10
#define ROUNDTOL  1e-10
#define CGTOL      1.0e-18 
#define CGFILE     "cghist.txt"
#define GRAVFILE   "gravin.dat"
#define BASICFILE  "basic.dat"
#define BCCFILE    "bcc.dat"
#define COORDFILE  "coord.dat"
#define BCVFILE    "bcv.dat"
#define CTRLFILE   "controls.fem"
#define PRTFILE    "print.dat"
#define TIMEFILE   "time.dat"
#define MONFILE    "monitor.fem"
#define ELEMFILE    "eldata.dat"
#define SURFFILE   "surfdata.dat"
#define BUOYFILE   "buoydata.dat"
#define FLTFILE    "fltdata.dat"
#define CGMAX        1000
/* 7 reals represent the B and S vectors and the slip magnitude */
/* #define NSPLITATTR       9 */
#define NSPLITATTR       24
#define REFINE_THRESHOLD      1.99e-4 
#define SLIP_TOL      1.0e-2 
#define BIG_G      6.674e-11 

/*  global parameters  */

EXTERN  int             attempt ;
EXTERN  int             iam ;
EXTERN  int             nproc ;
EXTERN  int             num_node_terms ;
EXTERN  int             num_face_terms ;
EXTERN  int             num_element_terms ;
/* EXTERN  char            msg[MAX_STRING_LENGTH] ; */

/******************************************/  
struct split_node_info
    {
     int  node ;
     int  ngrp ;
     int  grp[8] ;
     real ivec[3];
     real svec[3];
     real slip[8];
    };

/* a fault node can belong to up to 8 strands */

typedef struct split_node_info SPLITNODE ;

struct flow_control_struct
    {
     int solve_initial_elastic,
         write_initial_elas_solns,

         refine_elas_count,

         write_refined_mesh_sms,
         write_refined_mesh_toptris,
         write_refined_mesh_geofest,

         solve_refined_elas,
         write_refined_elas_soln,
         
         solve_visco,

         write_initial_post_solns,
         refine_post_count,
         refine_post_pct_goal,
         write_mesh_refine_post,
         solve_refined_post,
         write_refined_post_soln,
     
         write_initial_slip_solns,
         refine_slip_count,
         refine_slip_pct_goal,
         write_mesh_refine_slip,
         solve_refined_slip, 
         write_refined_slip_soln,
         first_soln_slips_all;

     real refine_elas_pct_goal;
    };
typedef struct flow_control_struct FLOW_SIM ;
EXTERN FLOW_SIM flow_control;

/******************************************/
 struct system_info
    {
     int               numnp ;
     int                 nsd ;
     int                 neq ;
     int           neq_owned ;
     int                ndof ;
     int              nrates ;
     int         save_shape  ;
     int         solver     ;
    } ;
typedef struct system_info SYSTEM_INFO ;
EXTERN  SYSTEM_INFO   fe_sys ;
EXTERN  SYSTEM_INFO   loc_sys ;

/******************************************/  
 struct  time_info
      {
       int   ntime_grp             ;
       int   *nsteps               ;
       real *alpha                ;  
       real *delt                 ;
       int   nreform               ;
       int   nprints               ;
       real *prt_time             ;
       int   nprt_nodes            ;
       int   loc_nprt_nodes            ;
       int   *node_list            ;
       int   nprt_elem             ;
       int   loc_nprt_elem             ;
       int   *elem_list            ;
       int   restart               ;
       int   sav_state             ;
       int    quake_prt            ;
       int    fail                 ;
       int    gravcalc             ;
       int    elastic              ;
       real  alpha_delt           ;
       real  dt                   ;
       real  time                 ;
       real  slip_interval        ;
       real  next_slip            ;
       real  traction_time[6]     ;
       real  currtsuf             ;
       int   itsuf                ;
       int   numtsuf              ;
       int   do_slip              ;
       int    paused               ;
       int    nbackup              ;
       int  nsub ;
       real *maxstep ;
       real *endtime ;
       real max_strain ;
       real scale_ratio ;
       int nfgrps       ;
       int step_started ;
      } ;
typedef   struct time_info   TIME_INFO  ;
EXTERN  TIME_INFO    time_data  ;
 
/******************************************/  
 struct global_data
    {
     int       n_group ;
     real     *coords ;
     real     *displ  ;
     real     *del_displ ;
     real     *forv ;
     real     *rate ;
     real     *bctime ;
     real     *pnode ;        /* nodal pressures for smoothing */
     real     *yp ;       /* diagonalized matrix for smoothing */
     int     *g_nodes ;
     int     *id_pointer ;
     int     *owned_eqs;   /* eqn. numbers this processor owns */
     int     *eqtype ;
     real     *splitn ; /* place to store split-node attributes */
     real   elas_energy ;
     int    current_flt_poll ;
     int    current_flt_status ;
     int     *permut ; /* profile-optimizing permutation of global nodes */
     char  regularInput[MAX_STRING_LENGTH];
     char  extendedInput[MAX_STRING_LENGTH];
     char  regularOutput[MAX_STRING_LENGTH];
     char  cgconvOutput[MAX_STRING_LENGTH];
     char  mesh0Output[MAX_STRING_LENGTH];
     char  mesh1toptris[MAX_STRING_LENGTH];
     char  mesh1Output[MAX_STRING_LENGTH];
     char  gravInput[MAX_STRING_LENGTH];
     char  startfile[MAX_STRING_LENGTH]         ; /* start file name */
     char  savefile[MAX_STRING_LENGTH]          ; /* save file name */
     char inputDirectoryPath[MAX_STRING_LENGTH]; /* input directory path string */
     char outputDirectoryPath[MAX_STRING_LENGTH]; /* output directory path string */
    } ;
typedef struct global_data GLOBAL_DATA ;
EXTERN  GLOBAL_DATA  global ;

/******************************************/  
 struct element_data
    {
     int     nel    ;
     real    str_energy;
     real     *stiff ;
     real     *shape ;
     real     *deter   ;
     real     *stress  ;
     real     *bta     ;
     real     *dbar    ;
     real     *dilat   ;
     int      *is_split ;
     int     mat    ;
     int     is_bc  ;
     int     degen  ;
     int     nint  ;
     int     *ien   ;
     int     *lm    ;
    } ;
typedef struct element_data ELEMENT_DATA ;

/******************************************/

struct buoy_data
{
  int   numface ;
  int   tally ;
  real  upvec[3] ;
  real  rho_g ;
  real  little_g ;
  real  delta_rho ;
  real  *dgrav ;
  int   *buoy_list ;
  int   rad_flag ;
  int   grav_out_flag ;
} ;
typedef  struct buoy_data BUOY_DATA ;

/******************************************/
 struct  fltgrp_data
 {
   real  interval           ;
   real  first_event        ;
   real  fric_coeff         ;
   real  fail_limit         ;
   real  fail_quorum        ;
   int   prev_num           ;
   int   due_now            ;
   int   aperiodic          ;
   real *schedule           ;
   real q_amount            ;
 } ;
typedef   struct fltgrp_data   FLTGRP  ;
EXTERN  FLTGRP    *fltgrp_ptr  ;

/******************************************/  
 struct element_info
    {
     int      type        ;
     int      numel       ;
     int      g_numel     ;
     int      numat       ;
     int      numsuf      ;
     int      numbuoy     ;
     int      nsplit      ; /* number split nodes in this processor */
     int      g_nsplit    ; /* global number split nodes */
     int      nfterm      ; /* number of node && element terms in this processor */
     int      nen         ;
     int      nee         ;
     int      el_size     ;
     int      stress_size ;
     int      dbar_size   ;
     int      failed      ;
     int    *g_elem     ;
     int    *surf_list0 ;
     int    *surf_list1 ;
     int    *surf_list2 ;
     int    *surf_list3 ;
     int    *surf_list4 ;
     int    *surf_list5 ;
     real   *surf_trac0 ;
     real   *surf_trac1 ;
     real   *surf_trac2 ;
     real   *surf_trac3 ;
     real   *surf_trac4 ;
     real   *surf_trac5 ;
/*   int    *(surf_list[6])  ;  */
/*   real   *(surf_trac[6])  ;  */
     int    *slip_list  ;
     real    *slip_val   ;
     int    *slip_grp  ;
     real    *slip_normal   ;
     BUOY_DATA    *buoy  ;
     real    *strain_list  ;
    } ;
typedef struct element_info ELEMENT_INFO ;

/******************************************/  
 struct element_mat
    {
     real     *mu     ;
     real     *lambda ;
     real     *visc   ;
     real     *bforce ;
     real     *density ;
     real     *big_g  ;
     real     *little_g  ;
     int     *plastic ;
    } ;
typedef struct element_mat ELEMENT_MAT ;

/******************************************/  
 struct group_dir
    {
     ELEMENT_DATA     *el_data ;
     ELEMENT_INFO     *el_info ;
     ELEMENT_MAT      *el_mat  ;
     int              group_num ;
    } ;
typedef struct group_dir GROUP ;
EXTERN  GROUP  groups[MAX_GROUP] ;

/******************************************/  
 struct rhs_data
    {
     real     *full_rhs        ;
     real     *incr_rhs        ;
     real     *last_result     ;
     real     *first_result    ;
     real     *net_external    ;
    } ;
typedef struct rhs_data RHS_DATA ;
EXTERN  RHS_DATA  force ;

/******************************************/  
 struct profile_data
    {
     real     *matrix ;
     int       *diag   ;
     int       neq     ;
     int       size    ;
    } ;
typedef struct profile_data PROFILE ;
EXTERN  PROFILE  stiff ;

 struct cg_data
    {
     int       neq     ;
     real     *precond ;
     real     *r ;
     real     *d ;
     real     *temp ;
     real     *hist ;
     real     norm ;
     real     first_norm ;
    } ;
typedef struct cg_data CONGRAD ;
EXTERN  CONGRAD  pcg ;

/******************************************/  
typedef real *COLUMN ;
COLUMN  *matcol ;
/******************************************/  

EXTERN  FILE  *in_file  ;
EXTERN  FILE  *out_file ;
EXTERN  FILE  *bcc_file ;
EXTERN  FILE  *coord_file ;
EXTERN  FILE  *bcv_file ;
EXTERN  FILE  *dump_file ;
EXTERN  FILE  *cgconv_file ;
EXTERN  FILE  *newmesh_file ;
EXTERN  FILE  *grav_file ;
EXTERN  FILE  *elem_file ;
EXTERN  FILE  *surf_file ;
EXTERN  FILE  *buoy_file ;
EXTERN  FILE  *flt_file ;
EXTERN  FILE  *print_file ;
EXTERN  FILE  *time_file ;
EXTERN  int  bugflag ;

#endif	/*GUARD_FINEL_H	*/
