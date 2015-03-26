/*
***                           File main.c                         ***
***                        GeoFEST version 6.0
*** Copyright (c) 2010, California Institute of Technology        ***
*** U.S.Sponsorship under NASA Contract NAS7-1407 is acknowledged ***
***
*** This software is designated for public release under JPL Task ***
*** Order Number NMO710991 and may be publicly released through
*** license with the Open Channel Foundation
***
*** This file contains the program main entry point, the main
*** task driver, and modules for driving high-level functions
*** and interactions with the operator:
***
***          == main == ++
***          == elastic == ++
***          == time_step == ++
***          == clear_stiff == ++
***          == elgrp_loop == ++
***          == node_load == ++
***          == accumulate == ++
***          == completion == ++
***          == wrt_save == ++
***          == vrestart == ++
***          == el_save == ++
***          == check_monitor == ++
***          == terminate == ++
***/

#define EXTERN
#include "main.h"
#include "utility.h"
#include "solver.h"
#include "inphase.h"
#include "generat.h"
#include "stiff.h"
#include "strain.h"

void input_phase();
void grav_calc();
void fail_loop();
void failcheck();


/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   main   ***************
*/

/*  INTERFACE:   */
     int
     main(
          int argc ,      /* number of input arguments */
          char *argv[]    /* input argument strings */
         )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine main
** main prints a banner, 
** interprets up to two args as input and output files,
** scans two comment lines in the input, 
** and divides (and times) the processing
** into input, elastic solution, time-stepping, and output.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      int slip_all;     /* flag: all faults slip for elastic scratch */
      char line[MAX_STRING_LENGTH]; /* temp space for read-in strings*/
      char            msg[MAX_STRING_LENGTH] ;

      /* 
      *  Deliver greeting.
      */
      squawk("\n\n**********************************************************\n");
      sprintf(msg,"Welcome to the Geophysical Finite Element Simulation Tool (%s)\n", argv[0]); 
      squawk(msg);
      squawk("\n**********************************************************\n\n");

      /* 
      *  The first GeoFEST command line option is the input directory path.
      */
      if( 1 < argc)
          strcpy(global.inputDirectoryPath,argv[1]);
      else
          strcpy(global.inputDirectoryPath,".");

      /* 
      *  Set outpuDirectoryPath from 2nd arg or use default "." 
      */
      if(2<argc) 
          strcpy(global.outputDirectoryPath,argv[2]);
      else 
          strcpy(global.outputDirectoryPath,".");

      sprintf(msg,"InputDirectoryPath set to '%s'\n", global.inputDirectoryPath);
      squawk(msg);
      sprintf(msg,"OutputDirectoryPath set to '%s'\n", global.outputDirectoryPath);
      squawk(msg);

/* open basic input data file */

      sprintf(msg,"%s/%s", global.inputDirectoryPath, BASICFILE);
      in_file = fopen(msg, "r");
      if (NULL == in_file) 
         {
          fprintf(stderr, "Couldn't open input file: %s!\n", msg);
          attempt = FILE_OPEN;
         }
      completion();

      squawk("\nReading basic GeoFEST input data \n");
      
      /* 
      *  Read output file name into string line from the input file.
      */
      
      fscanf(in_file, "%s", line);
      sprintf(global.regularOutput, "%s/%s", global.outputDirectoryPath, line);

          out_file = fopen(global.regularOutput, "w");
          if (NULL == out_file) 
             {
              fprintf(stderr, "Couldn't open output file %s!\n", global.regularOutput);
              attempt = FILE_OPEN;
             }
          else 
             {
              fprintf(stdout, "Opening output file %s ...\n\n", global.regularOutput);
             }
      completion();

/* open bcc input file */

      sprintf(msg,"%s/%s", global.inputDirectoryPath, BCCFILE);
      bcc_file = fopen(msg, "r");
      if (NULL == bcc_file) 
         {
          fprintf(stderr, "Couldn't open bcc input file: %s!\n", msg);
          attempt = FILE_OPEN;
         }
      completion();

/* open coord input file */

      sprintf(msg,"%s/%s", global.inputDirectoryPath, COORDFILE);
      coord_file = fopen(msg, "r");
      if (NULL == coord_file) 
         {
          fprintf(stderr, "Couldn't open coord input file: %s!\n", msg);
          attempt = FILE_OPEN;
         }
      completion();

/* open boundary disp/vel input file */

      sprintf(msg,"%s/%s", global.inputDirectoryPath, BCVFILE);
      bcv_file = fopen(msg, "r");
      if (NULL == bcv_file) 
         {
          fprintf(stderr, "Couldn't open boundary disp/vel input file: %s!\n", msg);
          attempt = FILE_OPEN;
         }
      completion();

/* open element input data file */

      sprintf(msg,"%s/%s", global.inputDirectoryPath, ELEMFILE);
      elem_file = fopen(msg, "r");
      if (NULL == elem_file) 
         {
          fprintf(stderr, "Couldn't open element input file: %s!\n", msg);
          attempt = FILE_OPEN;
         }
      completion();

      /*
      *  Open cgfile for iteration history
      *  also open an aux input file containing grav data
      */
      sprintf(global.cgconvOutput, "%s/%s", global.outputDirectoryPath, CGFILE);
      sprintf(global.gravInput, "%s/%s", global.inputDirectoryPath, GRAVFILE);

          cgconv_file = fopen(global.cgconvOutput, "w");
          if (NULL == cgconv_file) 
             {
              fprintf(stderr, "Couldn't open convergence history file %s!\n",
                 global.cgconvOutput);
              attempt = FILE_OPEN;
             }
      completion();
      /* check for presence of a gravity info data file */
      
     if(NULL == (grav_file = fopen(global.gravInput, "r")))
          {
           time_data.gravcalc = 0 ;  /* no gravity input; turn off flag */
          squawk("\nNo gravity file detected; proceeding without gravity calculation...\n") ;
          }
      else
         {
          time_data.gravcalc = 1 ;  /* set flag and leave the file open to read later */
          squawk("\nGravity file detected and opened...\n") ;
         }

      sprintf(msg,"%s/%s", global.inputDirectoryPath, SURFFILE);
      surf_file = fopen(msg, "r");
      if(NULL == surf_file)
          {
          squawk("\nNo surface traction file detected; proceeding...\n") ;
          }

      sprintf(msg,"%s/%s", global.inputDirectoryPath, BUOYFILE);
      buoy_file = fopen(msg, "r");
      if(NULL == buoy_file)
          {
          squawk("\nNo buoyancy file detected; proceeding...\n") ;
          }

      sprintf(msg,"%s/%s", global.inputDirectoryPath, FLTFILE);
      flt_file = fopen(msg, "r");
      if(NULL == flt_file)
          {
          squawk("\nNo fault node file detected; proceeding...\n") ;
          }

      /*
      * Read and echo the pass-thru comment lines 
      */
      fscanf(in_file,"%[^*]*\n",line) ;
      squawk("-------------Input file begins with comments:-----------\n");
          printf("%s\n",line) ;
          fprintf(out_file,"%s\n",line) ;

      fscanf(in_file,"%[^*]*\n",line) ;
          printf("%s\n",line) ;
          fprintf(out_file,"%s\n",line) ;
      squawk("\n---------------(End Input File Comments)--------------\n");
      
      /*
      *  Read in remainder of input by call to input_phase
      */

      input_phase() ;

      squawk("\n");

      completion() ;

      if(flow_control.solve_initial_elastic == 0) terminate() ;
      /* No solution requested - stop here. *********************************/

      if( !time_data.restart ) /* there is no prior run or restart file */
         {

          slip_all = flow_control.first_soln_slips_all;

          /*
          * Perform elastic solution - call "elastic"
          */
          squawk("\nStarting Elastic Solution...\n") ;

          elastic(slip_all) ;


          /* Initial elastic result output is requested.**********************/
          if(flow_control.write_initial_elas_solns == 1  && flow_control.solve_visco == 0) 
             {
              output_phase(0, 0) ;
             }

         } /* end of no-restart case that calls elastic: "!time_data.restart" */

      if(flow_control.solve_visco == 0) terminate() ; 
      
      /* No time stepping requested - stop here. *****************************/

      /*
      * Begin Time Stepping.
      */

         time_step() ;
         
      /* 
      * We are done.  Save state if requested (input file lists 
      * something other than NO_SAVE for save file name )
      */
          
      if(time_data.sav_state)
         wrt_save() ;
          
      terminate() ;
     } 
     
/****************************** end of main **********************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   elastic   ***************
*/

/*  INTERFACE:   */

void
elastic(int slip_all)

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine elastic
** elastic performs an elastic solution of the 
** finite element problem.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/
    {
     int  i ;
     char            msg[MAX_STRING_LENGTH] ;
    
     squawk("Starting elastic solution...\n") ;

     time_data.elastic = 1 ;
     clear_real(global.displ,loc_sys.numnp*fe_sys.ndof) ;
     clear_real(global.del_displ,loc_sys.numnp*fe_sys.ndof) ;
     clear_real(force.full_rhs,loc_sys.neq) ;
     clear_real(force.net_external,loc_sys.neq) ;
     clear_stiff() ;
     squawk("Storage arrays cleared...\n") ;

     elgrp_loop( BC_ELAS ) ;
     squawk("BC's computed...\n") ;
     time_data.do_slip = 0 ; /* obsolete; replaced by fault slip groups below */
     sprintf(msg,"Number of fault slip groups = %d\n",time_data.nfgrps) ;
     squawk(msg);
     
     if( time_data.ntime_grp == 0 )
        {
         /* in elastic-only problems, slip all faults */
         squawk("No time stepping; slip all faults now...\n") ;
         for(i=0 ; i<time_data.nfgrps ; i++)
            {
             time_data.do_slip = 1 ;
             (fltgrp_ptr+i)->due_now = 1 ; /* slip this fault at t=0 */
             (fltgrp_ptr+i)->first_event = ZERO ;
            }
        }
     else if (slip_all == 1)
        {
         squawk("Slip all faults now...\n") ;
         for(i=0 ; i<time_data.nfgrps ; i++)
            {
             time_data.do_slip = 1 ;
             (fltgrp_ptr+i)->due_now = 1 ; /* slip this fault at t=0 */
            }
        }
     else
        {
     
        for(i=0 ; i<time_data.nfgrps ; i++)
           {
            (fltgrp_ptr+i)->prev_num = 0 ;
            
            if( (fltgrp_ptr+i)->first_event == ZERO && (fltgrp_ptr+i)->interval >= ZERO)
               {
                /* A fault with a first event time of zero will move during the elastic step;
                   it will then not move again until a full repeat interval has elapsed 
                */
                time_data.do_slip = 1 ;
                (fltgrp_ptr+i)->due_now = 1 ; /* slip this fault at t=0 */
                sprintf(msg,"Fault #%d slip occurring at t = 0 \n",i+1) ;
                squawk(msg);
               }
            else
               {
                /* a fault with a non-zero first event time will not move during the elastic
                   step, but will at the first event time and every repeat interval thereafter 
                */
                (fltgrp_ptr+i)->due_now = 0 ;
                (fltgrp_ptr+i)->prev_num = -1 ;
               }
           } /* end for nfgrps */
        } /* end of not elastic-only "else" case */

     elgrp_loop( FORMS_ELAS ) ;
     squawk("Stiffness computed...\n") ;
     elgrp_loop( RHS_ELAS ) ;
     squawk("RHS done...\n") ;

     /*  update pcg.precond (preconditioner) */
          
     if (PCG == fe_sys.solver)
        {
         /* Finish forming preconditioner */
             { int ii;
               for (ii = 0; ii < loc_sys.neq; ii++) 
                  {
                   if(*(pcg.precond+ii)==ZERO) printf("bad precond value at ieq=%d\n",ii) ;
                   if( isnan(*(pcg.precond+ii) ) ) printf("NaN precond value at ieq=%d\n",ii) ;
                          *(pcg.precond+ii) = ONE / *(pcg.precond+ii);
                  }
             }
        }

     node_load(ZERO,YES) ;

     squawk("\nCalling SOLVER\n") ;

     if (DIRECT == fe_sys.solver)
        solver( FACBACK ) ;
     else if (PCG == fe_sys.solver)
        solver( ITER ) ;
     squawk("\nOut of SOLVER; begin accumulate()\n") ;

     accumulate() ;
     
     squawk("\nCalculating elastic stresses\n") ;
     elgrp_loop( E_STRESS ) ;
     squawk("\nStresses complete\n") ;
          
     squawk("Elastic solution complete.\n") ;
    }
/*************************** end of elastic **********************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   time_step   ***************
*/

/*  INTERFACE:   */

void
time_step()

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine time_step
** time_step computes the visco-elastic time stepping solutions
** to the finite element problem.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/
    {
     COUNTER  i , eventnum ;
     int   step_count , cslip , fcase , anyfail ;
     char            msg[MAX_STRING_LENGTH] ;

     int new_step = FALSE ;
     real time = time_data.time ;
     int time_group = 0 ;


     completion() ;
     output_phase(0, 0) ;
     if (0 == time_data.ntime_grp)  return;
     time_data.elastic = 0 ;

     squawk("\nStarting time-stepping algorithm...\n") ;

     force.incr_rhs = (real * ) calloc(loc_sys.neq,sizeof(real)) ;
     if (NULL == force.incr_rhs) attempt = R_MEM ;
     completion() ;
     clear_real(force.full_rhs,loc_sys.neq) ;


     step_count = 0 ;

     for ( time_group = 0 ; time_group < time_data.ntime_grp ; time_group++)
        {
         real dt = *(time_data.delt + time_group);
         real alpha = *(time_data.alpha + time_group);
         real endtime = *(time_data.endtime + time_group);
         int step = -1;

         time_data.nsteps[time_group] = 0;
         time_data.nsub = 1;
         /*
         *  Step through the time steps in this time group.
         */
         while (time_data.time < endtime - ROUNDTOL)
            {
             int reform = NO;
             real tstep = dt/(real)time_data.nsub;
             ++time_data.nsteps[time_group];
             ++step;
             ++step_count;
             time += tstep;
             time_data.dt = tstep;
             time_data.alpha_delt = tstep*alpha;
             time_data.time = time;
             sprintf(msg,"Start of time step; t=%f , delta_t=%f \n",time,tstep) ;
             squawk(msg);
             
             time_data.do_slip = 0;
             time_data.fail = 0 ;
             fcase = 0 ; /* by default no split node activity occurring */
             
             /*
             * Test fault groups for triggered slip  at this time.
             */
             anyfail = 0 ;
             for(i=0 ; i<time_data.nfgrps ; i++)
                {
                 if((fltgrp_ptr+i)->interval < ZERO)  /* active "quake" slip enabled */
                   {
                    /* testing and slipping of active fault */
                    (fltgrp_ptr+i)->due_now = 0 ;
                    global.current_flt_poll = i+1 ; /* must be one-based */
                    global.current_flt_status = FALSE ;
                    elgrp_loop( FAIL_CHECK ) ;
                    if(global.current_flt_status==TRUE)
                       {
                        anyfail = 1 ;
                        time_data.do_slip = 1 ;
                        (fltgrp_ptr+i)->due_now = 1 ;
                        sprintf(msg,"Triggered slip on fault #%d at t = %f \n",i+1,time) ;
                        squawk(msg);
                        time_data.fail = 1 ; /* freeze BCs during quake slip */
                        fail_loop(i) ; /* perform the quake slip on this strand */
                        time_data.do_slip = 0 ; /* reset fault flag for non-quakes */
                        time_data.fail = 0 ;
                       }
                   }
                } /*end i loop */
                
                
             /*
             * Test fault groups for scheduled slip event at this time.
             */
             for(i=0 ; i<time_data.nfgrps ; i++)
                {
                 (fltgrp_ptr+i)->due_now = 0 ; /* initialize */
                 cslip = 0 ;
                 
                 if((fltgrp_ptr+i)->interval == ZERO) /* this is a continuously slipping fault */
                   {
                    if(time >= (fltgrp_ptr+i)->first_event - ROUNDTOL)
                       {
						eventnum = (fltgrp_ptr+i)->prev_num + 1 ;
						((fltgrp_ptr+i)->prev_num)++ ;
						cslip = 1 ;
						fcase = 1 ;
						time_data.do_slip = 1 ;
						(fltgrp_ptr+i)->due_now = 2 ; /* flagged as continuously slipping */
						sprintf(msg,"Steadily slipping fault #%d at t = %f \n",i+1,time) ;
						squawk(msg);
					   }
                   }

                 else if((fltgrp_ptr+i)->interval < ZERO)  /* active "quake" slip */
                   {
                    /* active fault already dealt with - do nothing here */
                   }

                 else /* regular scheduled slip */
                   {
					 eventnum = (int)(
					   ( time + ROUNDTOL - (fltgrp_ptr+i)->first_event ) /
						   (fltgrp_ptr+i)->interval
					   ) ;
                 if( eventnum > (fltgrp_ptr+i)->prev_num && eventnum >= 0 
                     && time >= (fltgrp_ptr+i)->first_event - ROUNDTOL )
						{
						 fcase = 2 ;
						 time_data.do_slip = 1 ;
						 (fltgrp_ptr+i)->due_now = 1 ;
						 ((fltgrp_ptr+i)->prev_num)++ ;
						 sprintf(msg,"Scheduled fault #%d slip occurring at t = %f \n",i+1,time) ;
						 squawk(msg);
						}

                   }

                } /*end i loop */

             if( time_data.do_slip || anyfail ||
                  ((alpha != ZERO) && (step % time_data.nreform == 0)) )
                {
                 clear_stiff() ;
                 elgrp_loop( FORMS_STEP ) ;
                 reform = YES ;
                }
             /*   
             *   Recent change: postpone this clearing until after buoyancy calc is done
             *
             *          clear_real(global.del_displ, loc_sys.numnp*fe_sys.ndof) ;
             */         

             new_step = FALSE ;

             /*
             * Prepare to do next time step: set up velocity boundaries array
             */
             if(fe_sys.nrates == 0) /* this case: no velocity boundary conditions */
                {
                 clear_real(force.incr_rhs,loc_sys.neq) ;
                 clear_real(global.forv , loc_sys.numnp*fe_sys.ndof) ;
                }
             else
                {int n = 0;
                 for(n=0 ; n<fe_sys.nrates ; n++)
                    {
                     real test = *(global.bctime + n) ;
                     if(ZERO <= test && test <= time)
                        {
                         new_step = TRUE ;
                         move_real(global.rate + n*loc_sys.numnp*fe_sys.ndof,
                                    global.forv , loc_sys.numnp*fe_sys.ndof) ;
                         *(global.bctime + n) = (-ONE) ;
                         sprintf(msg,"Beginning boundary rate #%d at t = %f \n",n+1,time) ;
						 squawk(msg);

                         break ;
                        }
                    }
                }
  
             if(reform || new_step)  /* only redo the rhs if the bc or stiffness have changed */
                {
                 clear_real(force.incr_rhs,loc_sys.neq) ;
                 elgrp_loop( BC_STEP ) ;
                }


             incr_real(force.incr_rhs , force.full_rhs , loc_sys.neq) ;
             /* need to preserve the equil correction already in full_rhs from last step */


            time_data.max_strain = ZERO ;
      
             elgrp_loop( RHS_STEP ) ;
             
             sprintf(msg,"Stiffness and element RHS complete; make preconditioner...\n") ;
             squawk(msg);
               
             /* Finish forming preconditioner */
             if(reform == YES && PCG == fe_sys.solver)
                { 
                 int ii;
                 for (ii = 0; ii < loc_sys.neq; ii++) 
                    {
                     *(pcg.precond+ii) = ONE / *(pcg.precond+ii);
                     if(*(pcg.precond+ii)==ZERO) printf("bad precond value at ieq=%d\n",ii) ;
                     if( isnan(*(pcg.precond+ii) ) ) printf("NaN precond value at ieq=%d\n",ii) ;
                    }
                 sprintf(msg,"(Finishing preconditioner) Simulation t=%f...\n",time) ;
                 squawk(msg);
                } 
    
             /*   now we add any nodal load forces  (GL)  */
             node_load(time,reform) ;

             clear_real(global.del_displ, loc_sys.numnp*fe_sys.ndof) ;

             /* 
             *  Unlabled block: find max strain, adjust time step's nsub
             */
                { int new_nsub = time_data.nsub;
                 real global_max_strain;

/* Task MAX_STRAIN sorts the element strains and then stores in  
    time_data.max_strain the value of strain in the element whose strain exceeds that 
    in all but STRAIN_PCT% of the elements (suggest 1%)
*/
                 elgrp_loop( MAX_STRAIN ) ;
	
                    global_max_strain = time_data.max_strain;

                 if (global_max_strain > ZERO)
                    {
                     /*
                     * The largest viscoplastic strain change in any non-Newtonian element
                     * should not be greater than 1/STRAIN_FAC of the current elastic strain.
                     * The 'ratio' factor leads to the new nsub that must subdivide the
                     * original default time step in order to meet this criterion.
                     * Earlier versions used a STRAIN_FAC of 5, but this may be too small.
                     */
                     
                     real ratio = global_max_strain*STRAIN_FAC*time_data.nsub;
                     double       exp = log(ratio)/log(TWO);
                     exp = (double)((int)exp);
                     new_nsub = MIN ((int) (0.01 + pow(2.0, exp + 1.0)), MAXSUB);
                     new_nsub = MAX (1, new_nsub);
                    }
       
                 else
                    {
                     new_nsub = time_data.nsub ;
                    }
       
                 if (new_nsub != time_data.nsub)
                    {
                     int  failed = 0 ;  /* temporarily not implemented */

                     if( !failed )
                          time -= tstep ;

                     time_data.time = time ;

                     sprintf(msg,"Changing nsub from %d to %d at time = %f\n", 
                             time_data.nsub, new_nsub, time) ;
                     squawk(msg);

                     time_data.nsub = new_nsub ;
                     step = (-1) ;
                     continue ;    /* roll back a time step */
                    }
                } /* end of unlabled block introducing new_nsub, global_max_strain */
                            
             sprintf(msg,"Entering solver...\n") ;
             squawk(msg);

             if (PCG == fe_sys.solver)
                {
                 solver( ITER ) ;
                }
             else if(reform)
                {
                 solver( FACBACK ) ;
                }
             else
                {
                 solver( BACK ) ;
                }
                               
             accumulate() ;

             elgrp_loop( V_STRESS ) ;
  
             /*
             * perform stress equilibrium correction 
             */
             clear_real(force.full_rhs,loc_sys.neq) ;
             elgrp_loop( EQUIL ) ;
             incr_real(force.net_external , force.full_rhs , loc_sys.neq) ;

             output_phase( 0, 0) ;

             check_monitor( step_count+1) ;
                  
            /* fflush( stdout ) ; */
            }     /* end loop stepping time within time group */
        }   /* end loop over time groups  */
    }
             
/*************************** end of time_step ********************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   clear_stiff   ***************
*/

/*  INTERFACE:   */

void
clear_stiff(void)

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine clear_stiff
** clear_stiff nulls the stiffness matrix memory.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

{
    switch(fe_sys.solver)
       {
        case DIRECT:
         clear_real(stiff.matrix,stiff.size) ;
         break;
         
        case PCG:
         clear_real(pcg.precond,pcg.neq) ;
         break;
       }

}
/*************************** end of clear_stiff ******************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   elgrp_loop   ***************
*/

/*  INTERFACE:   */

void
     elgrp_loop(
                int  task      /* which element-level task to perform */
               )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine elgrp_loop
** elgrp_loop invokes a particular task that relies on 
** element structures and affects all the elements in 
** every element group.
**
** Supported tasks are
** GENERATE: read or generate element data, fill ien, lm arrays
** OUTPUT: compute and write out element centroid stresses
** SHAPE: compute parent-space shape functions and gradients
** FORMS_ELAS: compute elastic stiffness matrix
** FORMS_STEP: compute VE single-step stiffness matrix
** FORMS_QUAKE: compute elastic stiffness matrix for triggered fault event
** RHS_ELAS: compute elastic right-hand side vector for FE problem
** RHS_STEP: compute VE right-hand side vector for FE problem
** RHS_QUAKE: compute elastic right-hand side vector for triggered fault event
** BC_ELAS:  compute boundary-condition terms for elastic FE problem
** BC_STEP: compute boundary-condition terms for VE FE single step
** COL_HT: compute the column heights for profile stiffness storage
** E_STRESS: compute the stresses from elastic solution
** V_STRESS: compute the stresses from a VE step solution
** Q_STRESS: compute the stresses from a triggered fault event
** FAIL_CHECK: check the stresses for fault failure condition
** ESPROD: perform element-wise siffness-vector product for CG solver
** DUMP: dump information needed for restarts
** RESTORE: read and restore solution vector and info for a restart
** REORDER: analyze adjacency, find profile-optimizing node permutation for direct solver
** EQUIL: compute equilibrium stress corrections for visco steps
** MAX_STRAIN: scan the domain for the largest viscoplastic strain rates
** GRAV_CALC: calculate changes in gravity potential and geoid from disp solution
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

  {
   int i = 0;
   real xnorm2;
   char taskstr[MAX_STRING_LENGTH];
   char            msg[MAX_STRING_LENGTH] ;
   
   for(i=0 ; i<global.n_group ; i++)
      {
       switch( task )
               {
                case GENERATE:
                    gen_element( &groups[i] ) ;
                    strcpy(taskstr,"GENERATE=gen_element");
                    break ;
                case OUTPUT:
                    el_output( &groups[i] ) ;
                    strcpy(taskstr,"OUTPUT=el_output");
                    completion() ;
                    break ;
                case SHAPE:
                    p_shape( &groups[i] ) ;
                    strcpy(taskstr,"SHAPE=p_shape");
                    completion() ;
                    break ;
                case FORMS_ELAS:
                    form_stiff( &groups[i] , ELASTIC ) ;
                    strcpy(taskstr,"FORMS_ELAS=form_stiff(..ELASTIC)");
                    completion() ;
                    break ;
                case FORMS_STEP:
                    form_stiff( &groups[i] , VISCO ) ;
                    strcpy(taskstr,"FORMS_STEP=form_stiff(..VISCO)");
                    completion() ;
                    break ;
                case FORMS_QUAKE:
                    form_stiff( &groups[i] , QUAKEEVT ) ;
                    strcpy(taskstr,"FORMS_STEP=form_stiff(..QUAKE)");
                    completion() ;
                    break ;
                case RHS_ELAS:
                    form_rhs( &groups[i] , ELASTIC ) ;
                    strcpy(taskstr,"RHS_ELAS=form_rhs(..ELASTIC)");
                    completion() ;
                    break ;
                case RHS_STEP:
                    form_rhs( &groups[i] , VISCO ) ;
                    strcpy(taskstr,"RHS_STEP=form_rhs(..VISCO)");
                    completion() ;
                    break ;
                case RHS_QUAKE:
                    form_rhs( &groups[i] , QUAKEEVT ) ;
                    strcpy(taskstr,"RHS_STEP=form_rhs(..QUAKE)");
                    completion() ;
                    break ;
                case BC_ELAS:
                    form_bc( &groups[i] , ELASTIC ) ;
                    strcpy(taskstr,"BC_ELAS=form_bc(..ELASTIC)");
                    completion() ;
                    break ;
                case BC_STEP:
                    form_bc( &groups[i] , VISCO ) ;
                    strcpy(taskstr,"BC_STEP=form_bc(..VISCO)");
                    completion() ;
                    break ;
                case COL_HT:
                    colht( &groups[i] , &stiff ) ;
                    strcpy(taskstr,"COL_HT=colht");
                    completion() ;
                    break ;
                case E_STRESS:
                    form_stress(&groups[i] , ELASTIC) ;
                    strcpy(taskstr,"E_STRESS=form_stress(..ELASTIC)");
                    completion() ;
                    break ;
                case V_STRESS:
                    form_stress(&groups[i] , VISCO) ;
                    strcpy(taskstr,"V_STRESS=form_stress(..VISCO)");
                    completion() ;
                    break ;
                case Q_STRESS:
                    form_stress(&groups[i] , QUAKEEVT) ;
                    strcpy(taskstr,"Q_STRESS=form_stress(..QUAKE)");
                    completion() ;
                    break ;
                case FAIL_CHECK:
                    failcheck(&groups[i]) ;
                    strcpy(taskstr,"FAIL_CHECK=failcheck()");
                    completion() ;
                    break ;
                case ESPROD:
                    estiffprod( &groups[i] ) ;
                    strcpy(taskstr,"ESPROD=estiffprod");
                    break ;
                case DUMP:
                    el_save(&groups[i] , DUMP) ;
                    strcpy(taskstr,"DUMP=el_save(..DUMP)");
                    break ;
                case RESTORE:
                    el_save(&groups[i] , RESTORE) ;
                    strcpy(taskstr,"RESTORE=el_save(..RESTORE)");
                    break ;
                case REORDER:
                    reorder( &groups[i] );
                    strcpy(taskstr,"REORDER=reorder");
                    completion() ;
                    break;
                case EQUIL:
                    form_equil( &groups[i] );
                    strcpy(taskstr,"EQUIL=form_equil");
                    break;
                case MAX_STRAIN:
                    find_max_strain( &groups[i] );
                    break;
                case GRAV_CALC:
                    grav_calc( &groups[i] , ELASTIC  );
                    break;
                default:
                    attempt = BAD_TASK ;
                    completion() ;
                    break ;
               }

      }  /* end loop */

      /* Skip some: ESPROD, SHAPE, *STRESS do not yield valid RHS information. */
      if(!(ESPROD == task || SHAPE == task || E_STRESS == task || FORMS_STEP == task || 
           MAX_STRAIN == task || V_STRESS == task || OUTPUT == task)){
         if(force.full_rhs)
            xnorm2 = dot_real(force.full_rhs,force.full_rhs,loc_sys.neq_owned) ;   
         else
            xnorm2 = 0.0;
 /*      sprintf(msg, "After task %d: %s rhs norm squared is............... %lf\n", task,taskstr, xnorm2); */
    /*     squawk(msg); */
      }

  }
/*************************** end of elgrp_loop *******************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   node_load   ***************
*/

/*  INTERFACE:   */

void
     node_load(
               real  time ,   /* current simulation time */
               int reform     /* flag indicating stiffness reform status */
              )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine node_load
** node_load applies node-based forces that may accumulate
** with time.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

    {
     real   val ;
     int  k,  id , new_step ;
     int pgnode; /* one-based permuted node numbers */
    
     int ndof = fe_sys.ndof;
     int numnp = loc_sys.numnp;

     if(ZERO == time)          /* the ELASTIC case... */
        {
         /*  put in free nodal forces, if any  */
         int i = 0;
         for(i=0 ; i<numnp ; i++)
            {
             for(k=0 ; k<ndof ; k++)
                {
                 pgnode = global.permut[i+1];
                 id = *(global.id_pointer+(pgnode-1)*ndof+k) ;
                 val = *(global.forv+i*ndof+k) ;
                 if(val != ZERO && id != (-1))
                    {
                     *(force.full_rhs + id) += val ; /* Due to movement of call */
                     *(force.net_external + id) -= val ;
                    }
                }
            }
        }
     
   else 
      { /* the TIME-STEPPING case... */
   
  
/*  these tasks need to be performed  post-globalize in parallel code (GL) */

       int i = 0;
       for(i=0 ; i<numnp ; i++)
          {
           int k = 0;
           for(k=0 ; k<ndof ; k++)
              {
               int pgnode = global.permut[i+1];
               int id = *(global.id_pointer+(pgnode-1)*ndof+k) ;
               real val = *(global.forv+i*ndof+k) * time_data.dt ;
               if(val != ZERO && id != (-1))
                  {
                   *(force.full_rhs + id) += val ; /* Due to movement of call */
                   *(force.net_external + id) -= val ;
                  }
              }
          }
       }
    }
/*************************** end of node_load ********************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   accumulate   ***************
*/

/*  INTERFACE:   */

void
accumulate(void)

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine accumulate
** accumulate updates the total displacement based on the 
** single-step increment.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

    {
     COUNTER   i=0 ;

     for(i=0 ; i<loc_sys.numnp*fe_sys.ndof ; i++)
        {
         *(global.displ + i) += *(global.del_displ + i) ;
        }     
    }

/*************************** end of accumulate *******************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   completion   ***************
*/

int attempt = 0;   /* global variable definition */

/*  INTERFACE:   */

void
completion(void)

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine completion
** completion syncronizes multiple processors. 
** Keeps track of incurred errors in each processor
** Just a dummy call in sequential code
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

    {
     int   all = 0 , contrib = attempt , ierror ;
     
        all = attempt ;
        
     if(all)
        {
         squawk("An error occurred:\n") ;
         printf("\tIAM: %d\t  error code = %d\n", iam, attempt) ;
         fflush(stdout);
         exit(EXIT_FAILURE) ;
        }

    }

/*************************** end of completion *******************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   wrt_save   ***************
*/

/*  INTERFACE:   */
     void
     wrt_save()


/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine wrt_save
** wrt_save dumps state to disk for restarts.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

    {
     static real  outbuf[4] ;
         
     COUNTER ndof = fe_sys.ndof;
     COUNTER numnp = fe_sys.numnp;
     real  time = time_data.time;
     
     if (NULL == (dump_file = fopen(global.savefile, "w")))
          {
           printf("Dump file couldn't be opened!\n") ;
           return ;
          }

     printf("Saving simulation at time = %f to disk...\n",time) ;
     fprintf(dump_file,"%f\n",time) ;
     
  { COUNTER i = 0;
     for(i=0 ; i<numnp ; i++)
                 {
                  COUNTER index = i;
                  COUNTER k = 0;
                  for(k=0 ; k<ndof ; k++)
                         {
                          outbuf[k] =
                             *(global.displ + ndof*index + k) ;
                         }
                            
                  fprintf(dump_file,"%d",i+1) ;
                  for(k=0 ; k<ndof ; k++)
                     {
                      fprintf(dump_file,"  %f",outbuf[k]) ;
                     }
                  fprintf(dump_file,"\n") ;
                 }        
    }
             fflush(dump_file) ;
             elgrp_loop( DUMP ) ;
             fclose(dump_file) ;
             printf("Simulation dump complete.\n") ;
    }

/*************************** end of wrt_save *********************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   vrestart   ***************
*/

/*  INTERFACE:   */

void
vrestart()

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** vrestart restores state from previous dump and 
** continues simulation.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

    {
     static real  outbuf[4] ;
         
     COUNTER ndof = fe_sys.ndof;
     COUNTER numnp = fe_sys.numnp;
     real  time = time_data.time;
     
     if (NULL == (dump_file = fopen(global.startfile, "r"))) 
          {
           printf("Restart file %s couldn't be opened!\n",global.startfile) ;
           exit(1) ;
          }

     fscanf(dump_file,"%lf\n",&time) ;
     time_data.time = time ;
     printf("Restarting simulation at time = %f ...\n",time) ;
     
  { COUNTER i = 0;
     for(i=0 ; i<numnp ; i++)
        {
         int n = 0;
         fscanf(dump_file,"%d",&n) ;
          { COUNTER k = 0;
                  for(k=0 ; k<ndof ; k++)
                     {
                      fscanf(dump_file,"%lf",outbuf+k) ;
                     }
       }
          { COUNTER index = i;
         COUNTER   k = 0;
            for(k=0 ; k<ndof ; k++)
     {
      *(global.displ + ndof*index + k) =
      outbuf[k] ;
     }
          }        
      }
    }

  elgrp_loop( RESTORE ) ;
  fclose(dump_file) ;
  printf("Simulation restoration complete.\n") ;
    }

/*************************** end of vrestart *********************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   el_save   ***************
*/

/*  INTERFACE:   */
     void
     el_save(
             GROUP   *grp_ptr ,       /* pointer to current element group */
             int     code             /* code specifying which action */
            )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine el_save
** el_save dumps or restores element data.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

{
 if (DUMP == code) 
     {
   int numel = (grp_ptr->el_info)->numel;
   int str_size = (grp_ptr->el_info)->stress_size;
      fflush(dump_file) ;
      
      switch((grp_ptr->el_info)->type)
          {
           case BILIN:
           case SERENDIP:
           case TET:
      { int  nel = 0;
               for ( nel = 0 ; nel < numel ; nel++)
                   {
        int  n = nel + 1;
        ELEMENT_DATA* el_pt = grp_ptr->el_data + n - 1;
        int  nint = el_pt->nint;
        int  ipt = 0;
                    fflush(dump_file) ;

                    for (ipt = 0 ; ipt < nint ; ipt++)
                       {
                     int  j = 0;
                        for(j=0 ; j<str_size ; j++)
                           fprintf( dump_file,"%f  ",
                            *((el_pt -> stress) + ipt*str_size + j) ) ;
                        fprintf(dump_file,"\n") ;
                       }
                   }
 }
               fflush(dump_file) ;
               break ;
               
          }
     }

else if(RESTORE == code)
     {
      int numel = (grp_ptr->el_info)->numel;
      int str_size = (grp_ptr->el_info)->stress_size;
      
      switch((grp_ptr->el_info)->type)
          {
           case BILIN:
           case SERENDIP:
      { int  nel = 0;
               for ( nel = 0 ; nel < numel ; nel++)
                   {
     int  n = nel + 1;
     ELEMENT_DATA* el_pt = grp_ptr->el_data + n - 1;
     int  nint = el_pt->nint;
     int  ipt = 0;
                    for (ipt = 0 ; ipt < nint ; ipt++)
                       {
                        int  j = 0;
                        for(j=0 ; j<str_size ; j++)
                           {
                         real str_temp;
                            fscanf( dump_file,"%lf",
                                    &str_temp ) ;
                            if(n != 0)
                               {
                                el_pt = grp_ptr->el_data + n - 1 ;
                                *((el_pt -> stress) + ipt*str_size + j) =
                                  str_temp ;
                               }
                           }
                       }
                   }
 }
               break ;
               
          }
     }

    }

/**************************** end of el_save *********************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   check_monitor   ***************
*/

/*  INTERFACE:   */
     void
     check_monitor(
                   int   nstep    /* count of steps for backup */
                  )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine check_monitor
** check_monitor allows user run intervention via strings stored
** in the file monitor.fem.
**
** Supported options:
** OK: take no action - continue running normally
** REPORT: report cuurent time step to stdout and continue running normally
** KILL: close files and abort
** SAVE_STOP: write state, then abort
** SAVE_GO:  write dumps as requested in input, and continue
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

    {

/* These defines are replaced by "enum" type immediately following.
#define OK  0
#define KILL  1
#define SAVE_STOP 2
#define SAVE_GO  3
#define REPORT  4
#define NCODES  5
*/
     typedef enum {
        OK,
        KILL,
        SAVE_STOP,
        SAVE_GO,
        REPORT,
        NCODES
        } CODE;

     char  *mon_codes[NCODES] =
             { "OK" , "KILL" , "SAVE_STOP" , "SAVE_GO" , "REPORT" } ;
       
     int  flag = 0 ;
     int  n = 0 ;
     char line[MAX_STRING_LENGTH];
     char            msg[MAX_STRING_LENGTH] ;
     
     if(nstep % time_data.nbackup == 0)
        {
         sprintf(msg,"Simulation running normally at time = %f ; step size = %f\n",
            time_data.time , time_data.dt) ;
         squawk(msg);
         flag = 1 ;
        }
     n = NCODES;
     sprintf(msg,"%s/%s", global.inputDirectoryPath, MONFILE);
     if(NULL == (dump_file = fopen( msg , "r")))
          {
           return ;
          }

     fscanf(dump_file,"%s",line) ;

        { int i = 0;
         for(i=0 ; i<NCODES ; i++)
            {
             if(0 == strcmp(line, mon_codes[i]))
                {
                 n = i ;
                 break ;
                }
            }
        }

     switch( n )
        {
         case OK:
            fclose(dump_file) ;
            if(flag && time_data.sav_state)
               {
                wrt_save() ;
                squawk("Run backed up and continuing...\n") ;
               }
            return ;
            break ;

         case REPORT:
            if(0 == iam)
               {
                rewind(dump_file) ;
                fprintf(dump_file,"OK\n") ;
                printf("Simulation running normally at time = %f ; step size = %f\n",
                       time_data.time , time_data.dt) ;
               }
            fclose(dump_file) ;
            return ;
            break ;

         case KILL:
            fclose(in_file) ;
            if(0 == iam)
                {
                 fclose(out_file) ;
                 rewind(dump_file) ;
                 fprintf(dump_file,"OK\n") ;
                 printf("CPU time = %f\nTime at Forced Termination: ",
                         (double)clock()/CLOCKS_PER_SEC);
                }
            fclose(dump_file) ;
            squawk("Run aborted immediately by monitor request.\n") ;
            exit(EXIT_SUCCESS) ;
            break ;

         case SAVE_STOP:
            if(0 == iam)
                {
                 rewind(dump_file) ;
                 fprintf(dump_file,"OK\n") ;
                }
             fclose(dump_file) ;
             if(time_data.sav_state)
                wrt_save() ;
             else
                {
                 sprintf(global.savefile,"save.def") ;
                 wrt_save() ;
                }
             fclose(in_file) ;
                fclose(out_file) ;
                squawk("Run saved and halted by monitor request.\n") ;
                exit(EXIT_SUCCESS) ;
                break ;

         case SAVE_GO:
                 rewind(dump_file) ;
                 fprintf(dump_file,"OK\n") ;
            fclose(dump_file) ;
            if(time_data.sav_state)
               wrt_save() ;
            else
               {
                sprintf(global.savefile,"save.def") ;
                wrt_save() ;
               }
            squawk("Run backed up and continuing...\n") ;
            return ;
            break ;

         default:
            fclose(dump_file) ;
            return ;
            break ;
        }
    }

/************************** end of check_monitor *****************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*                                                                                                                                  
*************   ROUTINE:   terminate   ***************                                                                          
*/

/*  INTERFACE:   */

void
terminate()

 /*  RETURN VALUE:  -none- */
 /*  DESCRIPTION:   */
 /*                                                                                                                                  
 ** Routine terminate                                                                                                               
 ** Normal exit of the program from simulation flow                                                                        
 **/
 /*   EOP   */
 /*---------------------------------------------------------------------------*/
{
       fclose(in_file) ;
       
       fclose(out_file) ;

       squawk("Normal program termination.\n") ;


    exit( 0 );

}
/*************************** end of terminate **********************************/
