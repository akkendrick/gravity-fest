/*
***                           File inphase.c                      ***
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

#define EXTERN extern
#include "inphase.h"
#include "utility.h"
#include "output_phase.h"

void completion(void);
void gen_map();
void
gen_real(
		 real      *array ,  /* array of generated values */
		 int        dim   ,  /* dimensionality of grid */
		 FILE      *the_file  /* input file to read from */
		);
void
elgrp_loop(
           int	task	/* which element-level task to perform */
          );
void vrestart();
void
profile_diag(
	 PROFILE*	a,	/* ptr to assembled profile stiffness storage   */
	 int	profile_case	/* old flag for substructuring...		*/
            );
void
shape(
  ELEMENT_DATA*	el_pt,	/* ptr to current element struct		*/
  ELEMENT_INFO*	info,	/* ptr to element info struct			*/
  real*		sh_pt,	/* el storage to hold shape functions		*/
  real*		det_pt,	/* el storage to hold determinant		*/
  int		type,	/* element type number				*/
  int		nen,	/* number of element nodes			*/
  int		nint,	/* number of Gauss integration points		*/
  int		all	/* flag indicating we need derivatives also	*/
      );
      
real   dotsh() ;

static real
   sh_temp[256] , det_temp[9] ;

char  line[80] ;

/**********************************************************************/




/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   input_phase   ***************
*/

/*  INTERFACE:   */

void
input_phase()

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

    {
     FILE *c_file;
     COUNTER  i , n , onoff ;
     int   numnp , ndof , nrates , nsd ; 
     char  simtask[10] ;
     char *test;
     char default_string[MAX_STRING_LENGTH];
     char in_string[MAX_STRING_LENGTH];
     char controlfile_string[MAX_STRING_LENGTH];
     char            msg[MAX_STRING_LENGTH] ;
     
/* *************************** SIM CONTROL CODES ***************** */
/* 
*  set defaults : 
*  sequential case: skip refinement, just do elastic (and out), visco
*
* Notes on scheduled slip:
* default: scheduled faults slip.  May be overridden so all faults slip
*  at first elastic solution on coarse mesh, for refinement purposes. 
*  Note may be overridden by flow control line 
*  (eg, in  run directory file controls.fem) with line similar to
*  FIRST_SOLN_SLIPS_ALL 1
*/
    squawk("\nFlow control flags - internal defaults:\n\n");
    if(1)  /* sequential case  ON */
       {
        set_flow_from_string("ELASTIC1 1");
        set_flow_from_string("ELAS_OUT1 1");
        set_flow_from_string("FIRST_SOLN_SLIPS_ALL 0");
        set_flow_from_string("REFINE 0");
        set_flow_from_string("REFINE_ELASTIC_PERCENT_GOAL 0.00");
        set_flow_from_string("REFINE_OUT_SMS 0");
        set_flow_from_string("REFINE_OUT_TOPTRIS 0");
        set_flow_from_string("REFINE_OUT 0");
        set_flow_from_string("ELASTIC2 0");
        set_flow_from_string("ELAS_OUT2 0");
        set_flow_from_string("VISCO 1");
       }
    else /* parallel case OFF */
       {
        set_flow_from_string("ELASTIC1 1");
        set_flow_from_string("ELAS_OUT1 0");
        set_flow_from_string("FIRST_SOLN_SLIPS_ALL 0");
        set_flow_from_string("REFINE 1");
        set_flow_from_string("REFINE_ELASTIC_PERCENT_GOAL 5.78");
        set_flow_from_string("REFINE_OUT_SMS 1");
        set_flow_from_string("REFINE_OUT_TOPTRIS 1");
        set_flow_from_string("REFINE_OUT 0");
        set_flow_from_string("ELASTIC2 1");
        set_flow_from_string("ELAS_OUT2 1");
        set_flow_from_string("VISCO 1");
       }
    
    /* 
    * "while" loop to skip blank lines, since using fgets not fscanf  here.
    * This loop should handle any number of blank lines and lines starting with
    * non-digits. First line starting with a digit breaks out.
    */

    squawk("\nFlow control flags - regular input file:\n\n");

    while(TRUE)
       {
        test = fgets(in_string,MAX_STRING_LENGTH,in_file);

        if(test == NULL)  /* abnormal end of file or hard error */
           {attempt = FLOW_CODES; completion();}

        if(strlen(in_string) > 1) /* skip blank line -- strlen==1*/
           {
            if(isdigit(in_string[0]))
               break;

            /* not a digit or a blank line: process flow line, any order. */
            set_flow_from_string(in_string);
           }
       }
    /* 
       Check for an override in controls.fem file, if any
    */
     sprintf(msg,"%s/%s", global.inputDirectoryPath, CTRLFILE);
     if(NULL == (c_file = fopen( msg , "r")))
         squawk("No controls.fem file - continuing.\n\n"); 
     else
        {
         squawk("\nFlow control flags - **overrides** - file \"controls.fem\":\n\n");
         while(TRUE)
            {
             test = fgets(controlfile_string,MAX_STRING_LENGTH,c_file);
             if(test == NULL)  /* presumably this is end of file*/
                break;

             if(strlen(test) > 1)
                {
                 set_flow_from_string(controlfile_string);
                }
            }
        }

     squawk("\nFINAL control flags:\n\n");
     sprintf(msg,"ELASTIC1(solve initial elastic) => %d\n",
        flow_control.solve_initial_elastic); squawk(msg);
     sprintf(msg,"ELAS_OUT1(write initial elastic solutions) => %d\n",
        flow_control.write_initial_elas_solns); squawk(msg);
     sprintf(msg,"FIRST_SOLN_SLIPS_ALL(first solution slips all faults) => %d\n",
        flow_control.first_soln_slips_all); squawk(msg);
     sprintf(msg,"REFINE(refine elastic solution count) => %d\n",
        flow_control.refine_elas_count); squawk(msg);
     sprintf(msg,"REFINE_ELASTIC_PERCENT_GOAL(refine elastic per cent goal) => %lf\n",
        flow_control.refine_elas_pct_goal); squawk(msg);
     sprintf(msg,"REFINE_OUT_SMS(write refined mesh sms file) => %d\n",
        flow_control.write_refined_mesh_sms); squawk(msg);
     sprintf(msg,"REFINE_OUT_TOPTRIS(write refined mesh toptris file) => %d\n",
        flow_control.write_refined_mesh_toptris); squawk(msg);
     sprintf(msg,"REFINE_OUT(write refined mesh geofest file) => %d\n",
        flow_control.write_refined_mesh_geofest); squawk(msg);
     sprintf(msg,"ELASTIC2(solve refined elastic) => %d\n",
        flow_control.solve_refined_elas); squawk(msg);
     sprintf(msg,"ELAS_OUT2(write refined elastic solution) => %d\n",
        flow_control.write_refined_elas_soln); squawk(msg);
     sprintf(msg,"VISCO(solve visco) => %d\n",
        flow_control.solve_visco); squawk(msg);


/* *************************** GLOBAL INPUT PHASE ***************** */

    squawk("\nNow reading global parameters from basic input file:\n\n");
    
 /*  here we'll read in numnp, nsd, ndof, nrates, save_shape, solver, numgroups,
     nprintnodes, nprintels, ntimegroups, nreform, nbackup, nfltgroups
 */
    while(TRUE)
       {
        test = fgets(in_string,MAX_STRING_LENGTH,in_file);

        if(test == NULL)  /* abnormal end of file or hard error */
           {attempt = FLOW_CODES; completion();}

        if(strlen(in_string) > 1) /* skip blank line -- strlen==1*/
           {
            if(isdigit(in_string[0]))
               break;

            /* not a digit or a blank line: process flow line, any order. */
            set_param_from_string(in_string);
           }
       }


     sprintf(msg,"Global number of nodes = %d\n\n",fe_sys.numnp) ;
     squawk(msg);

     sprintf(msg,
        "nsd = %d  ;  ndof = %d  ;  nrates = %d  ;  save_shape = %d  ;  solver = %d\n\n",
        fe_sys.nsd , fe_sys.ndof , fe_sys.nrates , fe_sys.save_shape , fe_sys.solver ) ;
     squawk(msg);

        loc_sys.numnp = fe_sys.numnp ;

     numnp = loc_sys.numnp;
     nsd = fe_sys.nsd ;
     ndof    = fe_sys.ndof    ;
     nrates  = fe_sys.nrates  ;
     
     global.coords = (real * ) calloc(numnp*nsd,sizeof(real)) ;
     if (global.coords == NULL) attempt = CRD_MEM ;
     completion() ;
     
     global.displ = (real * ) calloc(numnp*ndof,sizeof(real)) ;
     if (global.displ == NULL) attempt = CRD_MEM ;
     completion() ;
     
     global.del_displ = (real * ) calloc(numnp*ndof,sizeof(real)) ;
     if (global.del_displ == NULL) attempt = CRD_MEM ;
     completion() ;
     
     global.forv = (real * ) calloc(numnp*ndof,sizeof(real)) ;
     if (global.forv == NULL) attempt = CRD_MEM ;
     completion() ;
     
     global.eqtype = (int * ) calloc(numnp*ndof,sizeof(int)) ;
     if (global.eqtype == NULL) attempt = CRD_MEM ;
     completion() ;
     
     global.splitn = (real * ) calloc(numnp*NSPLITATTR,sizeof(real)) ;
     if (global.splitn == NULL) attempt = CRD_MEM ;
     completion() ;
     
     if( nrates > 0 )
        {
         global.bctime = (real * ) calloc(nrates,sizeof(real)) ;
         global.rate = (real * ) calloc(numnp*ndof*nrates,
                                         sizeof(real)) ;
         if (global.rate == NULL) attempt = CRD_MEM ;
        }
     completion() ;
     
     gen_map() ;  /*  fill the decomp map array  */
     squawk("Global b.c. map generated...\n") ;
         
     gen_number() ;  /* obtain local node and equation numbering */
        
     squawk("Equation numbering complete...\n") ; 
    
      sprintf(msg,"%s/%s", global.inputDirectoryPath, COORDFILE);
      coord_file = fopen(msg, "r");
      if (NULL == coord_file) 
         {
          fprintf(stderr, "Couldn't open coord input file: %s!\n", msg);
          attempt = FILE_OPEN;
         }
      completion();

     gen_real( global.coords,nsd,coord_file ) ;  /*  fill the coord array  */
     squawk("Coordinate generation complete...\n") ;
     
      sprintf(msg,"%s/%s", global.inputDirectoryPath, BCVFILE);
      bcv_file = fopen(msg, "r");
      if (NULL == bcv_file) 
         {
          fprintf(stderr, "Couldn't open coord input file: %s!\n", msg);
          attempt = FILE_OPEN;
         }
      completion();

     gen_real( global.forv,ndof,bcv_file ) ; /*  fill the b. c. array  */
     squawk("Displacement boundary value generation complete...\n") ;
     
     for( n=0 ; n<nrates ; n++ )
        {
         fscanf(bcv_file,"%lf",global.bctime+n) ;
         sprintf(msg,"Reading in boundary velocity set #%d for time = %g...\n",
		n+1 , *(global.bctime+n) ) ;
         squawk(msg);

         gen_real( (global.rate + n*ndof*numnp),ndof,bcv_file ) ;
         /*  fill the velocity arrays  */
        }
     squawk("Velocity boundary rates generation complete...\n") ;
     

     
/* ************************ ELEMENT INPUT PHASE ***************** */

     if(global.n_group > MAX_GROUP)
         {
          squawk("Too many element groups\n") ;
          exit(1) ;
         }

     for(i=0 ; i<global.n_group ; i++)
         {
          groups[i].el_info =
             (ELEMENT_INFO *) calloc(1,sizeof(ELEMENT_INFO)) ;
          if (groups[i].el_info == NULL) attempt = ELI_MEM + i ;
          completion() ;

          groups[i].el_mat =
             (ELEMENT_MAT *) calloc(1,sizeof(ELEMENT_MAT)) ;
          if (groups[i].el_mat == NULL) attempt = MAT_MEM + i ;
          completion() ;
          
          groups[i].group_num = i+1 ;
         }
                  
     squawk("\nStarting element generation...\n") ;

     elgrp_loop( GENERATE ) ;   /* read in and generate element geometry */

     squawk("Doing shape functions...\n") ;

     elgrp_loop( SHAPE ) ;    /* calculate element shape functions
                                 for future use */

     squawk("Element generation complete.  Setting up linear algebra...\n") ;
     

     
/* ************************ EQUATION SETUP PHASE ***************** */

     matrix_alloc() ;      /*  set up storage and indexing for the 
                              skyline profile elimination solver OR
                              the preconditioned CG solver  */


     squawk("Equation setup completed.\n") ;


/* ************************ TIME INPUT PHASE ***************** */

     time_data.time = ZERO ;
     time_data.loc_nprt_nodes = ZERO ;
     time_data.loc_nprt_elem = ZERO ;
     
  /*   fscanf(in_file,"%d%d",&time_data.nprt_nodes,&time_data.nprt_elem) ; */
     sprintf(msg,"%s/%s", global.inputDirectoryPath, PRTFILE);
     if(NULL == (print_file = fopen( msg , "r")))
         squawk("No print instruction file - problem...\n\n"); 

     sprintf(msg,"number of printed nodes = %d\nnumber of printed elements = %d\n",
        time_data.nprt_nodes,time_data.nprt_elem) ;
     squawk(msg);

     if(time_data.nprt_nodes > 0)
        {
         time_data.node_list
            = ( int  * ) calloc(time_data.nprt_nodes,sizeof(int)) ;
         if (time_data.node_list == NULL) attempt = TGP_MEM ;
         completion() ;

         for( i = 0 ; i < time_data.nprt_nodes ; i++ )
            fscanf(print_file,"%d", time_data.node_list + i ) ;

        }

     if(time_data.nprt_elem > 0)
        {
         time_data.elem_list
            = ( int  * ) calloc(2*time_data.nprt_elem,sizeof(int)) ;
         if (time_data.elem_list == NULL) attempt = TGP_MEM ;
         completion() ;

         for( i = 0 ; i < time_data.nprt_elem ; i++ )
            fscanf(print_file,"%d%d", time_data.elem_list + 2*i
                   , time_data.elem_list + 2*i + 1 ) ;

                                    /* group # and elem. # */

        }

/*
     fscanf(in_file,"%d%d%d%d",&time_data.ntime_grp,&time_data.nreform,
                           &time_data.nbackup,&time_data.nfgrps) ;
*/
     sprintf(msg,"number of time groups = %d\nnumber of steps per reform = %d\n",
        time_data.ntime_grp,time_data.nreform) ;
     squawk(msg);
     sprintf(msg,"number of steps per backup = %d\nnumber of fault groups = %d\n",
        time_data.nbackup,time_data.nfgrps) ;
     squawk(msg);
     
     if(time_data.ntime_grp == 0)  /* this is an elastic-only run */
        {
         time_data.nprints = 1 ;
         time_data.prt_time
             = ( real  * ) calloc(time_data.nprints,sizeof(real)) ;
         *(time_data.prt_time) = ZERO ;
         fltgrp_ptr = (FLTGRP *) calloc(time_data.nfgrps,sizeof(FLTGRP)) ;
         squawk("No time data  generation (elastic only)...\n") ;
         return ;
        }
     
     /* the rest of the input required only for time-stepping runs... */
     sprintf(msg,"%s/%s", global.inputDirectoryPath, TIMEFILE);
     if(NULL == (time_file = fopen( msg , "r")))
         squawk("No time instruction file - problem...\n\n"); 
     
     time_data.maxstep
        = ( real  * ) calloc(time_data.ntime_grp,sizeof(real)) ;
     if (time_data.maxstep == NULL) attempt = TGP_MEM ;
     time_data.endtime
        = ( real  * ) calloc(time_data.ntime_grp,sizeof(real)) ;
     if (time_data.endtime == NULL) attempt = TGP_MEM ;
     time_data.nsteps
        = ( int  * ) calloc(time_data.ntime_grp,sizeof(int)) ;
     if (time_data.nsteps == NULL) attempt = TGP_MEM ;
     completion() ;
                 
     time_data.alpha
         = ( real  * ) calloc(time_data.ntime_grp,sizeof(real)) ;
     if (time_data.alpha == NULL) attempt = TGP_MEM ;
     completion() ;
                 
     time_data.delt 
         = ( real  * ) calloc(time_data.ntime_grp,sizeof(real)) ;
     if (time_data.delt == NULL) attempt = TGP_MEM ;
     completion() ;
                 
     for( i = 0 ; i < time_data.ntime_grp ; i++ )
        fscanf(time_file,"%lf%lf%lf",
             time_data.endtime + i, time_data.alpha + i , time_data.delt + i ) ;

     fscanf(time_file,"%d", &time_data.nprints ) ;
     if( time_data.nprints > 0 )
         {
          time_data.prt_time
             = ( real  * ) calloc(time_data.nprints,sizeof(real)) ;
          if (time_data.prt_time == NULL) attempt = TGP_MEM ;

          for( i = 0 ; i < time_data.nprints ; i++ )
            fscanf(time_file,"%lf", time_data.prt_time + i ) ;
         }
         
 /* now comes the input block for split node slip history groups...*/
 
     fltgrp_ptr = (FLTGRP *) calloc(time_data.nfgrps,sizeof(FLTGRP)) ;
     for(i=0 ; i<time_data.nfgrps ; i++)
        {
         fscanf(time_file,"%d",&n) ;
         fscanf(time_file,"%lf",&((fltgrp_ptr+n-1)->interval) ) ;
         fscanf(time_file,"%lf",&((fltgrp_ptr+n-1)->first_event) ) ;
 /* Three new parameters added to support active failing faults!!! */
 /*
     fric_coeff - Coulomb friction coefficient for this fault in calculating failure
     fail_limit - Coulomb stress (in appropriate units) for fault failure
     fail_quorum - net fraction (0.0 to 1.0) of fault area * stress exceeding
                   fail_limit in order to declare failure
 */
         fscanf(time_file,"%lf",&((fltgrp_ptr+n-1)->fric_coeff) ) ;
         fscanf(time_file,"%lf",&((fltgrp_ptr+n-1)->fail_limit) ) ;
         fscanf(time_file,"%lf",&((fltgrp_ptr+n-1)->fail_quorum) ) ;
         
         printf("Fault group #%d:\ninterval=%g   first time=%g   friction coef=%g   fail lim=%g   fail pct=%g\n\n",
                n,(fltgrp_ptr+n-1)->interval,(fltgrp_ptr+n-1)->first_event,(fltgrp_ptr+n-1)->fric_coeff,
                (fltgrp_ptr+n-1)->fail_limit,(fltgrp_ptr+n-1)->fail_quorum);
        }

     squawk("Time data  generation complete...\n") ;
     
     time_data.restart = 0 ;
     fscanf(in_file,"%s",line) ;
     if(strcmp(line , "NO_RESTART") != 0 )
        {
         sprintf(global.startfile,"%s",line) ;
         time_data.restart = 1 ;
         vrestart() ;
        }

     time_data.sav_state = 0 ;
     fscanf(in_file,"%s",line) ;
     if(strcmp(line , "NO_SAVE") != 0 )
        {
         sprintf(global.savefile,"%s",line) ;
         time_data.sav_state = 1 ;
        }

    } /* end function */ 

/************************** end of input_phase *******************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   matrix_alloc   ***************
*/

/*  INTERFACE:   */

void
matrix_alloc()

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine matrix_alloc
** allocates space for the sparse stiffness matrix _or_
** the PCG solver arrays
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

   {
    int   i ;
    
    switch(fe_sys.solver)
       {
        case DIRECT:
         stiff.neq = fe_sys.neq = loc_sys.neq ; /* only valid for single processor! */

         stiff.diag = (int  * ) calloc(fe_sys.neq,sizeof(int )) ;
         if (stiff.diag == NULL) attempt = DIAG_MEM ;
         completion() ;
         
         elgrp_loop( REORDER );
         elgrp_loop( COL_HT ) ;
         profile_diag( &stiff , INTERIOR ) ;
         
         printf("Total profile size is %d\n",stiff.size);

         stiff.matrix =
            (real * ) calloc(stiff.size,sizeof(real)) ;
         if (stiff.matrix == NULL) attempt = PROF_MEM ;
         completion() ;

         break;
/*---------------------------------------------------------------*/
         
        case PCG:
         pcg.neq = loc_sys.neq ;

         pcg.precond = (real  * ) calloc(loc_sys.neq,sizeof(real)) ;
         if (pcg.precond == NULL) attempt = PCG_MEM ;

         pcg.r = (real  * ) calloc(loc_sys.neq,sizeof(real)) ;
         if (pcg.r == NULL) attempt = PCG_MEM ;

         pcg.d = (real  * ) calloc(loc_sys.neq,sizeof(real)) ;
         if (pcg.d == NULL) attempt = PCG_MEM ;

         pcg.temp = (real  * ) calloc(loc_sys.neq,sizeof(real)) ;
         if (pcg.temp == NULL) attempt = PCG_MEM ;

         completion() ;
         
         pcg.hist = (real  * ) calloc(CGMAX,sizeof(real)) ;

         global.permut = (int *) calloc(loc_sys.numnp+1,sizeof(int));
         for(i=0;i<=loc_sys.numnp;i++)
           global.permut[i] = i ;  /* null permutation */
    /*   printf("Null node permutation allocated with size %d\n",loc_sys.numnp) ; */

         force.last_result = (real * ) calloc(loc_sys.neq,sizeof(real)) ;
         if (force.last_result == NULL) attempt = R_MEM ;

         force.first_result = (real * ) calloc(loc_sys.neq,sizeof(real)) ;
         if (force.first_result == NULL) attempt = R_MEM ;

         break;
       }
       
         force.full_rhs = (real * ) calloc(loc_sys.neq,sizeof(real)) ;
         if (force.full_rhs == NULL) attempt = R_MEM ;
         force.net_external = (real * ) calloc(loc_sys.neq,sizeof(real)) ;
         if (force.net_external == NULL) attempt = R_MEM ;
         completion() ;

   }

/************************** end of matrix_alloc ******************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   gen_number   ***************
*/

/*  INTERFACE:   */

void
gen_number()

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

   {
    COUNTER   i , j , k , index ;
    int   eqtype[3] ;
    int   ieq_local, ieq_owned ;
    int   *id_pointer ;
    int   neq_local , ndof , numnp ;
    char            msg[MAX_STRING_LENGTH] ;

     neq_local = 0 ;
     numnp = loc_sys.numnp ;
     ndof = fe_sys.ndof ; 

   
     for( i=0 ; i< numnp ; i++ )
         {
          index = i ;

          for(k=0 ; k<ndof ; k++)
             {
              eqtype[k] = *(global.eqtype+ndof*index+k) ;
             }
                      /* 1 = active , 0 = boundary  */

               for(k=0 ; k<ndof ; k++)
                  if(eqtype[k])  neq_local++ ;
              
         }   /* end loop over nodes */

     loc_sys.neq = neq_local ;
     squawk("Decomposed node and equation totals:\n") ;
     completion() ;
     printf ("Number of nodes = %d,  equations = %d\n",
	 loc_sys.numnp, loc_sys.neq) ;
     fflush(stdout) ;
     
  /* we now will assign equation numbers */
  
  /* first, allocate space for the id array  */
  /* second, allocate space for the inverse id array  */
     global.id_pointer =
        (int * ) calloc(numnp*ndof,sizeof(int)) ;
     if (global.id_pointer == NULL) attempt = ID_MEM ;
     global.owned_eqs =
        (int * ) calloc(loc_sys.neq,sizeof(int)) ;
     if (global.owned_eqs == NULL) attempt = ID_MEM ;
     completion() ;

     id_pointer = global.id_pointer ;
     
  /* now assign equations ...    (equations start at zero)  */
   
     ieq_local = 0 ;
     ieq_owned = 0 ;
     
     for(j=0 ; j<numnp*ndof ; j++) *(id_pointer + j) = (-1) ;
     
     for( i=0 ; i<numnp ; i++ )
         {
          /*index = pamr_my_node_gid( m, i+1) ;*/ /* return global index */
          index = i+1 ;

          for(k=0 ; k<ndof ; k++)
             {
              eqtype[k] = *(global.eqtype+ndof*i+k) ;
             }

          for(k=0 ; k<ndof ; k++)
             {
              if(eqtype[k])
                 {
                  *(id_pointer + ndof*(i) + k) = ieq_local ;
            /*
                  if (!pamr_my_partition_node( m, index ) ||
                       pamr_my_node_leader(m, index)) {
                       global.owned_eqs[ieq_owned] = ieq_local;
                       ieq_owned++ ;
                  }
             */
                   global.owned_eqs[ieq_owned] = ieq_local;
                   ieq_owned++ ;
                  ieq_local++ ;
                 }
             }
         }   /* end loop over global node number */
         
         
     loc_sys.neq_owned = ieq_owned;
/*     MPI_Allreduce((void *)&loc_sys.neq_owned,(void *)&fe_sys.neq_owned,1,
                     MPI_INT,MPI_SUM,MPI_COMM_WORLD); 
*/
     fe_sys.neq_owned = loc_sys.neq_owned ;
     sprintf(msg,"Global number of equations: %d\n",fe_sys.neq_owned); 
     squawk(msg);
        

   }

/************************** end of gen_number ********************************/ 



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
               )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine output_phase
** output_phase writes out requested data at a quake 
** or scheduled time.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

    {
     COUNTER  i , n , k , ntot , index , ndof , numnp , numnp_all , nsd ;
     int  mine , print , curr_node , count , prt_node ;
     static real  outbuf[9] ;
     real  time;
         
     int owned_nodes, out_num, tot_num, j, s_size, tmp_sum, indx, global_node;
     int *rcounts, *displs;

     ndof = fe_sys.ndof ;
     numnp = loc_sys.numnp ;
     numnp_all = fe_sys.numnp ;
     nsd = fe_sys.nsd ;
     time = time_data.time ;
     
/*  
   Here "time" is steadily increasing, vs.  the print times (prt_time[])
   Note that once used, prt_time[n] is set to -ONE.
   So for each item in prt_time[] we skip the -ONE's, 
   then check if "time" has met or exceeded a remaining prt_time;
   we count it if it's (ROUNDTOL) close, so 15.9999999999 matches 16.
   If so, set that prt_time to -ONE and "print" flag to 1.
*/

   
     print = 0 ;
     for(n=0 ; n<time_data.nprints ; n++) {
         if(time_data.prt_time[n] == (-ONE)){ continue; }
         else if (time > time_data.prt_time[n] - ROUNDTOL)
            {
             time_data.prt_time[n] = (-ONE);
             print = 1;
             break;
            }
         else
            {
             break;
            }
     }

     if (redo)
       print = 1;

     if( !(print || quake) )
        {
         printf("No output scheduled at time = %f ; (next at %f)\n",time,time_data.prt_time[n]) ;
         fflush(stdout);
         return ;
        }

	if(quake)
	    printf("Performing quake output at time = %f ...\n",time) ;
	else
	    printf("\n ***Performing scheduled output at time = %f ***\n\n",time) ;

    fflush(stdout);

	fprintf(out_file,
	"\n Global coordinates & displacements & delt displacements \n") ;
	fprintf(out_file,"\n Simulation time = %f ; step size = %g\n",
	    time,time_data.dt) ;
	if(quake)
	  fprintf(out_file,"***Fault Rupture*** at t=%g\n",time) ;
	fprintf(out_file,
 "\n          coordinates          displacements       del_displacements\n") ;
     
      if (time_data.nprt_nodes == (-1))
	      ntot = numnp_all ;
      else
          ntot = time_data.nprt_nodes ;
             
      for(i=0 ; i<ntot ; i++) {

          if(time_data.nprt_nodes == (-1))
             curr_node = i+1 ;
          else
             curr_node = *(time_data.node_list+i) ;

	     index = curr_node - 1 ;
	     prt_node = curr_node ;

         fprintf(out_file,"node %d    %g   %g   %g      %g   %g   %g      %g   %g   %g\n",
                  curr_node, *(global.coords + nsd*index), *(global.coords + nsd*index + 1),
                  *(global.coords + nsd*index + 2), *(global.displ + ndof*index),
                  *(global.displ + ndof*index + 1), *(global.displ + ndof*index + 2),
                  *(global.del_displ + ndof*index), *(global.del_displ + ndof*index + 1),
                  *(global.del_displ + ndof*index + 2) ) ;
      }        


     elgrp_loop(OUTPUT) ;
     
     if( time_data.gravcalc )
        {
         squawk("Calculating gravity changes...\n") ;
         elgrp_loop( GRAV_CALC ) ;
        }

     squawk("Output step finished.\n") ;

    }

/************************** end of output_phase ******************************/ 




/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   el_output   ***************
*/

/*  INTERFACE:   */

      void
      el_output(
                GROUP		*grp_ptr/* pointer to current element group */
               )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine el_output
** el_output prints out element stress information.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      int g_numel , nint , nel , use_list , n , str_size , inode1 , inode2 ,
	  ipt , j , ng , ntot , nsd , mine , prt_n , count , index ;
      real   xpt[3][10] , str_ctr[9] ;
      ELEMENT_DATA   *el_pt ;

      int out_num, tot_num, i, k, s_size, tmp_sum, indx ;
      int *rcounts, *displs;


          if(time_data.nprt_elem == -1)
             out_num =  (grp_ptr->el_info)->numel ;
          else
             out_num = time_data.nprt_elem ;
      
      



      nsd = fe_sys.nsd ;
      g_numel = (grp_ptr->el_info)->g_numel ;
      str_size = (grp_ptr->el_info)->stress_size ;

	 fprintf(out_file,
	   "\n Element stresses  -- element group #%d \n",grp_ptr->group_num) ;
	 fprintf(out_file,"\n Simulation time = %f ; step size = %g\n",
	   time_data.time,time_data.dt) ;
	 fprintf(out_file,
      "\n coordinates                          stresses . . .->      el.#\n") ;
     if(str_size == 6)
	 fprintf(out_file,
      "   X             Y            Z       sig_XX        sig_YY      sig_ZZ      sig_XY        sig_XZ      sig_YZ\n") ;
     
     else if(str_size == 4)
	 fprintf(out_file,
      "   X             Y       sig_XX        sig_YY      sig_XY        sig_ZZ\n") ;
     
	fflush(out_file) ;

      
      switch((grp_ptr->el_info)->type)
          {
           case BILIN:
           case SERENDIP:
           case TET:
           case TRILIN:
               if(time_data.nprt_elem == (-1))
                  {
                   ntot = g_numel ;
                   use_list = 0 ;
                  }
               else
                  {
                   ntot = time_data.nprt_elem ;
                   use_list = 1 ;
                  }
               for ( nel = 0, indx = 0 ; nel < ntot ; nel++)
                   {
                    if( use_list )
                       {
                        if( *(time_data.elem_list + 2*nel) !=
                            grp_ptr->group_num ) continue ;
                        ng = *(time_data.elem_list + 2*nel + 1) ;
                        n = ng ;
                       }
                    else
                       {
                        ng = nel+1 ;
                        n = ng ;
                       }
                    if( n == 0 ) continue ;
                    
                      index = n - 1 ;
                      prt_n = n ;

                      el_pt = grp_ptr->el_data + index ;

                      nint = el_pt->nint ;
                      
                      if((grp_ptr->el_info)->type == TET && el_pt->degen == 1)
                           /* this is a "duo" element - no integration points  */
                         {
                          for(j=0 ; j<str_size ; j++)
                             str_ctr[j] = *(el_pt -> stress + j) ;
                             
                          inode1 = *(el_pt->ien    ) - 1 ;  /* zero-based node index */
                          inode2 = *(el_pt->ien + 1) - 1 ;
                          for(j=0 ; j < nsd  ; j++)
                             str_ctr[str_size+j] =  PT5 * (*(global.coords+inode1*nsd+j) +
                                       *(global.coords+inode2*nsd+j) );
                          continue ;
                         }
                      
                      locate_pt( el_pt , grp_ptr->el_info , xpt ) ;

                      for(j=0 ; j<str_size ; j++)
                          str_ctr[j] = ZERO ;  /* midpoint stresses */

                      for(j=0 ; j < nsd ; j++)
                          xpt[j][9] = ZERO ;  /* midpoint location */

                      for (ipt = 0 ; ipt < nint ; ipt++)
                        {
                           for(j=0; j < nsd; j++)
                               xpt[j][9] += xpt[j][ipt] / ((real) nint) ;

                           for(j=0 ; j<str_size ; j++)
                               str_ctr[j] +=
                               *((el_pt -> stress) + ipt*str_size + j)  / ((real) nint) ;
                        }

                      for(j=0 ; j < nsd; j++)
   		          str_ctr[str_size+j] = xpt[j][9] ;

		  fprintf(out_file,"%g   %g   %g      %g   %g   %g   %g   %g   %g   %d\n",
		          xpt[0][9],xpt[1][9],xpt[2][9], str_ctr[0],str_ctr[1],str_ctr[2],
		          str_ctr[3],str_ctr[4],str_ctr[5],prt_n) ;
		}




          break;
          }
    }

/************************** end of el_output *********************************/ 



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
               )

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

  {
   int    j , ipt , save_shape , type , nint , nen , nsd , nsdp1 ;
   real  *shape_ptr , *det_ptr , *sh ;

   save_shape = fe_sys.save_shape ;
   nen = info->nen ;
   type = info->type ;
   nint = el_pt->nint ;
   nsd = fe_sys.nsd ;
   nsdp1 = nsd + 1 ;
   
   if( !save_shape )
        {
         shape_ptr = sh_temp ;
         det_ptr = det_temp ;
         shape( el_pt , info , shape_ptr , det_ptr , type ,
                nen , nint , 0 ) ;
        }
   else
        {
         shape_ptr = el_pt -> shape ;
        }

     for (ipt = 0 ; ipt < nint ; ipt++)
          {
           sh = shape_ptr + ipt*nen*nsdp1 ;

           for(j=0; j < nsd ; j++)
              xpt[j][ipt] = dotsh( global.coords ,
                  nsd , j , sh , el_pt , info , nen , nsd , 0 , ELASTIC) ;
          }

  }

/************************** end of locate_pt *********************************/ 




/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   set_flow_from_string   ***************
*/

/*  INTERFACE:   */

      void
        set_flow_from_string(
            char in_string[]
                            )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine set_flow_from_string
** reads 1-line string (obtained from a file) for one flow_control name value pair
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/
{
    char  simtask[MAX_STRING_LENGTH] ;
    char value[MAX_STRING_LENGTH];
    char            msg[MAX_STRING_LENGTH] ;
    int ivalue;
    real fvalue;

    sscanf(in_string,"%s%s",simtask,value);
    ivalue = atoi(value); /* ignored for float cases below */

    if(0  == strcmp(simtask, "ELASTIC1"))
       flow_control.solve_initial_elastic = ivalue ;
    else if(0  == strcmp(simtask, "ELAS_OUT1"))
       flow_control.write_initial_elas_solns = ivalue ;
    else if(0  == strcmp(simtask, "REFINE"))
       flow_control.refine_elas_count = ivalue ;
    else if(0  == strcmp(simtask, "REFINE_OUT_SMS"))
       flow_control.write_refined_mesh_sms = ivalue ;
    else if(0  == strcmp(simtask, "REFINE_OUT_TOPTRIS"))
       flow_control.write_refined_mesh_toptris= ivalue ;
    else if(0  == strcmp(simtask, "REFINE_OUT"))
       flow_control.write_refined_mesh_geofest= ivalue ;
    else if(0  == strcmp(simtask, "ELASTIC2"))
       flow_control.solve_refined_elas = ivalue ;
    else if(0  == strcmp(simtask, "ELAS_OUT2"))
       flow_control.write_refined_elas_soln = ivalue ;
    else if(0  == strcmp(simtask, "VISCO"))
       flow_control.solve_visco = ivalue ;
    else if(0  == strcmp(simtask, "REFINE_ELASTIC_PERCENT_GOAL"))
       {
       sscanf(value,"%lf",&flow_control.refine_elas_pct_goal);
       fvalue = flow_control.refine_elas_pct_goal;
       }
    else if(0 == strcmp(simtask, "FIRST_SOLN_SLIPS_ALL"))
       flow_control.first_soln_slips_all = ivalue ;
    else
       {
        attempt = FLOW_CODES; completion();
       }

    /*
    * Report that this value is set. 
    * Treat ints different from floats.
    */
    if(0  != strcmp(simtask, "REFINE_ELASTIC_PERCENT_GOAL"))
       {
        sprintf(msg,"Set: %s -> %d\n",simtask,ivalue);
       }
    else
       {
        sprintf(msg,"Set: %s -> %lf\n",simtask,fvalue);
       }
    squawk(msg);
}
/************************** end of set_flow_from_string ********************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   set_param_from_string   ***************
*/

/*  INTERFACE:   */

      void
        set_param_from_string(
            char in_string[]
                            )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine set_param_from_string
** reads 1-line string (obtained from a file) for one parameter name value pair
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/
{
    char  simtask[MAX_STRING_LENGTH] ;
    char value[MAX_STRING_LENGTH];
    char            msg[MAX_STRING_LENGTH] ;
    int ivalue;
    real fvalue;

    sscanf(in_string,"%s%s",simtask,value);
    ivalue = atoi(value); /* ignored for float cases below */

    if(0  == strcmp(simtask, "NUMNP"))
       {
        ivalue = atoi(value);
        fe_sys.numnp = ivalue ;
       }
    else if(0  == strcmp(simtask, "NSD"))
       {
        ivalue = atoi(value);
        fe_sys.nsd = ivalue ;
       }
    else if(0  == strcmp(simtask, "NDOF"))
       {
        ivalue = atoi(value);
        fe_sys.ndof = ivalue ;
       }
    else if(0  == strcmp(simtask, "NRATES"))
       {
        ivalue = atoi(value);
        fe_sys.nrates = ivalue ;
       }
    else if(0  == strcmp(simtask, "SAVE_SHAPE"))
       {
        ivalue = atoi(value);
        fe_sys.save_shape = ivalue ;
       }
    else if(0  == strcmp(simtask, "SOLVER"))
       {
        ivalue = atoi(value);
        fe_sys.solver = ivalue ;
       }
    else if(0  == strcmp(simtask, "NUMGROUPS"))
       {
        ivalue = atoi(value);
        global.n_group = ivalue ;
       }
    else if(0  == strcmp(simtask, "NPRINTNODES"))
       {
        ivalue = atoi(value);
        time_data.nprt_nodes = ivalue ;
       }
    else if(0  == strcmp(simtask, "NPRINTELS"))
       {
        ivalue = atoi(value);
        time_data.nprt_elem = ivalue ;
       }
    else if(0  == strcmp(simtask, "NTIMEGROUPS"))
       {
        ivalue = atoi(value);
        time_data.ntime_grp = ivalue ;
       }
    else if(0  == strcmp(simtask, "NREFORM"))
       {
        ivalue = atoi(value);
        time_data.nreform = ivalue ;
       }
    else if(0  == strcmp(simtask, "NBACKUP"))
       {
        ivalue = atoi(value);
        time_data.nbackup = ivalue ;
       }
    else if(0  == strcmp(simtask, "NFLTGROUPS"))
       {
        ivalue = atoi(value);
        time_data.nfgrps = ivalue ;
       }
    else
       {
        squawk("Input parameter not understood!\n");
        attempt = FLOW_CODES; completion();
       }

    /*
    * Report that this value is set. 
    */
        sprintf(msg,"Set: %s -> %d\n",simtask,ivalue);

    squawk(msg);
}
/************************** end of set_param_from_string ********************************/ 
