/*
***                           File quake.c                        ***
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



#define EXTERN extern
#include "quake.h"
#include "stiff.h"
static real
    sh_temp[256] , det_temp[9] ;



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   failcheck   ***************
*/

/*  INTERFACE:   */

     void
     failcheck(
               GROUP   *grp_ptr    /* pointer to current element group */
              )

/*  RETURN VALUE:  n/a */
/*  DESCRIPTION:   */
/*
** Routine failcheck
** failcheck looks at the elements comprising the specified split node strand
** and determines whether a failure condition has been met, returning a flag
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      int    n , numel , nel , ndof , numsuf , nfterm , side , type , nint , ipt ,
             nen , nsd , nsdp1 , node , nee , i , numface , skip , cslip , rad_code ;
      int   *surf_list , *slip_list , *buoy_list , *slip_grp , ifltgrp , save_shape ,
             str_size ;
      real  *surf_val , *slip_val , *slip_normal , *trac , *slip , *up_ptr , rho_g ,
            *vnorm , slip_amt , sum_failstr , sum_comparstr , limit , fail_pct , fric ,
            *shape_ptr , *det_ptr , *wt , *sh , w_det , *str , px , py , pz , pn ,
            psq , ps , pc ;
      ELEMENT_DATA   *el_pt ;
      BUOY_DATA      *buoy_ptr ;


      save_shape = fe_sys.save_shape;
      nsd = fe_sys.nsd ;
      nsdp1 = nsd + 1 ;
      numel = (grp_ptr->el_info)->numel ;
      nfterm = (grp_ptr->el_info)->nfterm ;
      ndof = fe_sys.ndof ;
	  slip_list = (grp_ptr->el_info) -> slip_list ;
	  slip_val = (grp_ptr->el_info) -> slip_val ;
	  slip_grp = (grp_ptr->el_info) -> slip_grp ;
	  slip_normal = (grp_ptr->el_info) -> slip_normal ;
      str_size = (grp_ptr->el_info)->stress_size ;
      nen  = (grp_ptr->el_info)->nen ;
      type = (grp_ptr->el_info)->type ;

      sum_failstr = ZERO ;
      sum_comparstr = ZERO ;
      for(n=0 ; n<nfterm ; n++)
		   {
		    if( *(slip_grp+n) != global.current_flt_poll ) continue ; 
		          /* filter the split nodes belonging to the currently polled group */
				  
			ifltgrp = *(slip_grp+n) - 1 ;
			limit = (fltgrp_ptr+ifltgrp)->fail_limit ;
			fail_pct = (fltgrp_ptr+ifltgrp)->fail_quorum ;
			fric = (fltgrp_ptr+ifltgrp)->fric_coeff ;
			nel = (int) *(slip_list + 2*n) ;
			node = (int) *(slip_list + 2*n + 1) ;
			vnorm = slip_normal+ndof*n ;
			el_pt = grp_ptr->el_data + nel ;
			nint = el_pt->nint ;
		    switch( nint )
				 {
				  case 4:
					 wt = ((BILIN == type) ? wt_2x2 : wt_tet) ;
					 break ;
				  case 6:
					 wt = wt_tri ;
					 break ;
				  case 8:
					 wt = wt_2x2x2 ;
					 break ;
				  case 9:
					 wt = wt_3x3 ;
					 break ;
				 }

			  if( !save_shape )
				 {
				  shape_ptr = sh_temp ;
				  det_ptr = det_temp ;
				  shape( el_pt , grp_ptr->el_info , shape_ptr , det_ptr , type , nen , nint , 1 ) ;
				 }
			  else
				 {
				  shape_ptr = el_pt -> shape ;
				  det_ptr = el_pt -> deter ;
				 }

		   for (ipt = 0 ; ipt < nint ; ipt++)
			   {
				sh = shape_ptr + ipt*nen*nsdp1 ;
				w_det = wt[ipt] * det_ptr[ipt];
				str = (el_pt -> stress) + ipt*str_size ;
				
				px = str[0]*vnorm[0]+str[3]*vnorm[1]+str[4]*vnorm[2] ; /* traction components */
				py = str[3]*vnorm[0]+str[1]*vnorm[1]+str[5]*vnorm[2] ;
				pz = str[4]*vnorm[0]+str[5]*vnorm[1]+str[2]*vnorm[2] ;
				pn = px*vnorm[0]+py*vnorm[1]+pz*vnorm[2] ;   /* normal stress */
				psq = px*px+py*py+pz*pz ;
				ps = sqrt(psq-pn*pn) ;    /* shear stress */
				pc = ps + fric*pn ;    /* Coulomb failure stress */
/*
        printf("Shear stress=%g  Coulomb stress=%g  vnorm=(%g %g %g)  wdet=%g  lim=%g\n",
                ps , pc , vnorm[0] , vnorm[1] , vnorm[2] , w_det , limit) ;
*/

				sum_failstr += w_det*pc ;
				sum_comparstr += w_det*limit ;
			   }

		   }
        if(sum_failstr/sum_comparstr > fail_pct) global.current_flt_status = TRUE ;
        printf("Checking active fault #%d: current=%g  limit=%g\n",
               global.current_flt_poll,sum_failstr,sum_comparstr) ;


          
     }

/************************** end of failcheck ********************************/ 



 
/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   fail_loop   ***************
*/

/*  INTERFACE:   */

     void
     fail_loop(
                int     f_index         /* 0-based fault strand identifier */
               )

/*  RETURN VALUE:  n/a */
/*  DESCRIPTION:   */
/*
** Routine fail_loop
** fail_loop is the topmost controlling routine that for a specific split node
** strand that has been declared 'failed',  applies varying amounts of
** slip with the objective of minimizing global strain energy
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      int  reform , nextslot , lastslot ;
      real  slip_delta , energy1 , energy2 , energy3 , total_slip ;
      real  kink , min_slip , concav ;
      
     
      clear_stiff() ;
      clear_real(force.full_rhs,loc_sys.neq) ; /* we elect not to do equil correction during fail_loop() */
      elgrp_loop( FORMS_QUAKE ) ;  reform = YES ;

      slip_delta = PT5 ;
      total_slip = ZERO ;
      (fltgrp_ptr+f_index)->q_amount = ONE - slip_delta ;
      total_slip += (fltgrp_ptr+f_index)->q_amount ;
      elgrp_loop( RHS_QUAKE ) ;  /* first slip (guess-delta) */
               

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
		} 
    
      clear_real(global.del_displ, loc_sys.numnp*fe_sys.ndof) ;
                   
      if (PCG == fe_sys.solver)
		{
		 solver( ITER ) ;
		}
      else
		{
		 solver( FACBACK ) ;
		}
		
      accumulate() ;
      global.elas_energy = ZERO ;
      elgrp_loop( Q_STRESS ) ;
      energy1 = global.elas_energy ;

      clear_real(force.full_rhs,loc_sys.neq) ;

  /*  first exploratory slip complete...  */
  printf("\nslip delta=%g  applied=%g  slip total= %g  strain energy=%g\n\n",
  slip_delta,(fltgrp_ptr+f_index)->q_amount,total_slip,global.elas_energy) ;

      (fltgrp_ptr+f_index)->q_amount = slip_delta ;
      total_slip += (fltgrp_ptr+f_index)->q_amount ;
      elgrp_loop( RHS_QUAKE ) ;  /* 2nd slip (+delta) */
      clear_real(global.del_displ, loc_sys.numnp*fe_sys.ndof) ;
                   
      if (PCG == fe_sys.solver)
		{
		 solver( ITER ) ;
		}
      else
		{
		 solver( BACK ) ;
		}
		
      accumulate() ;
      global.elas_energy = ZERO ;
      elgrp_loop( Q_STRESS ) ;
      energy2 = global.elas_energy ;
      clear_real(force.full_rhs,loc_sys.neq) ;

  /*  2nd exploratory slip complete...  */
  printf("\nslip delta=%g  applied=%g  slip total= %g  strain energy=%g\n\n",
  slip_delta,(fltgrp_ptr+f_index)->q_amount,total_slip,global.elas_energy) ;

      (fltgrp_ptr+f_index)->q_amount = slip_delta ;
      total_slip += (fltgrp_ptr+f_index)->q_amount ;
      elgrp_loop( RHS_QUAKE ) ;  /* 3rd slip (+delta) */
      clear_real(global.del_displ, loc_sys.numnp*fe_sys.ndof) ;
                   
      if (PCG == fe_sys.solver)
		{
		 solver( ITER ) ;
		}
      else
		{
		 solver( BACK ) ;
		}
		
      accumulate() ;
      global.elas_energy = ZERO ;
      elgrp_loop( Q_STRESS ) ;
      energy3 = global.elas_energy ;
      clear_real(force.full_rhs,loc_sys.neq) ;

  /*  3rd exploratory slip complete... now estimate position of energy minimum...  */
  printf("\nslip delta=%g  applied=%g  slip total= %g  strain energy=%g\n\n",
  slip_delta,(fltgrp_ptr+f_index)->q_amount,total_slip,global.elas_energy) ;
  
      kink = (energy1 - energy2) / (energy2 - energy3) ;
      min_slip = PT5 - PT25*(THREE*kink - ONE)/(ONE - kink) ;
      concav = (energy2 - energy1) / (0.75 - min_slip) ;
      
      if( concav <= ZERO )
         {
          printf("\nWARNING - no energy minimum found for fault slip; aborting\n") ;
		  (fltgrp_ptr+f_index)->q_amount = -THREE*slip_delta;
		  total_slip += (fltgrp_ptr+f_index)->q_amount ;
		  elgrp_loop( RHS_QUAKE ) ;  /* return to zero total slip (-3*delta) */
		  clear_real(global.del_displ, loc_sys.numnp*fe_sys.ndof) ;
					   
		  if (PCG == fe_sys.solver)
			{
			 solver( ITER ) ;
			}
		  else
			{
			 solver( BACK ) ;
			}
			
		  accumulate() ;
		  global.elas_energy = ZERO ;
		  elgrp_loop( Q_STRESS ) ;
		  clear_real(force.full_rhs,loc_sys.neq) ;
		  printf("\nslip delta=%g  applied=%g  slip total= %g  strain energy=%g\n\n",
                 slip_delta,(fltgrp_ptr+f_index)->q_amount,total_slip,global.elas_energy) ;
          printf("\nNet slip results:\n                time   strand   net slip   energy\n") ;
          printf("EVENT*** %g   %d   %g   %g\n\n", time_data.time , f_index+1 , total_slip , global.elas_energy) ;
         }
    
      else if( min_slip <= ZERO )
         {
          printf("\nWARNING - non-positive slip found for fault; aborting\n") ;
		  (fltgrp_ptr+f_index)->q_amount = -THREE*slip_delta;
		  total_slip += (fltgrp_ptr+f_index)->q_amount ;
		  elgrp_loop( RHS_QUAKE ) ;  /* return to zero total slip (-3*delta) */
		  clear_real(global.del_displ, loc_sys.numnp*fe_sys.ndof) ;
					   
		  if (PCG == fe_sys.solver)
			{
			 solver( ITER ) ;
			}
		  else
			{
			 solver( BACK ) ;
			}
			
		  accumulate() ;
		  global.elas_energy = ZERO ;
		  elgrp_loop( Q_STRESS ) ;
		  clear_real(force.full_rhs,loc_sys.neq) ;
		  printf("\nslip delta=%g  applied=%g  slip total= %g  strain energy=%g\n\n",
                 slip_delta,(fltgrp_ptr+f_index)->q_amount,total_slip,global.elas_energy) ;
          printf("\nNet slip results:\n                time   strand   net slip   energy\n") ;
          printf("EVENT** %g   %d   %g   %g\n\n", time_data.time , f_index+1 , total_slip , global.elas_energy) ;
         }


      else if( min_slip > ONE )
         {
          printf("\n Slip maxed out for fault; applying max\n") ;
		  (fltgrp_ptr+f_index)->q_amount = -slip_delta;
		  total_slip += (fltgrp_ptr+f_index)->q_amount ;
		  elgrp_loop( RHS_QUAKE ) ;  /* return to max nominal slip (-delta) */
		  clear_real(global.del_displ, loc_sys.numnp*fe_sys.ndof) ;
					   
		  if (PCG == fe_sys.solver)
			{
			 solver( ITER ) ;
			}
		  else
			{
			 solver( BACK ) ;
			}
			
		  accumulate() ;
		  global.elas_energy = ZERO ;
		  elgrp_loop( Q_STRESS ) ;
		  clear_real(force.full_rhs,loc_sys.neq) ;
		  printf("\nslip delta=%g  applied=%g  slip total= %g  strain energy=%g\n\n",
                 slip_delta,(fltgrp_ptr+f_index)->q_amount,total_slip,global.elas_energy) ;
          printf("\nNet slip results:\n                time   strand   net slip   energy\n") ;
          printf("EVENT** %g   %d   %g   %g\n\n", time_data.time , f_index+1 , total_slip , global.elas_energy) ;
         }


      else
         {
		  (fltgrp_ptr+f_index)->q_amount = min_slip-THREE*slip_delta;
		  total_slip += (fltgrp_ptr+f_index)->q_amount ;
		  elgrp_loop( RHS_QUAKE ) ;  /* slip by amount to minimize energy (quadratic fit) */
		  clear_real(global.del_displ, loc_sys.numnp*fe_sys.ndof) ;
					   
		  if (PCG == fe_sys.solver)
			{
			 solver( ITER ) ;
			}
		  else
			{
			 solver( BACK ) ;
			}
			
		  accumulate() ;
		  global.elas_energy = ZERO ;
		  elgrp_loop( Q_STRESS ) ;
		  clear_real(force.full_rhs,loc_sys.neq) ;
		  printf("\nslip delta=%g  applied=%g  slip total= %g  strain energy=%g\n\n",
                 slip_delta,(fltgrp_ptr+f_index)->q_amount,total_slip,global.elas_energy) ;
          printf("\nNet slip results:\n                time   strand   net slip   energy\n") ;
          printf("EVENT* %g   %d   %g   %g\n\n", time_data.time , f_index+1 , total_slip , global.elas_energy) ;
         }
         
     }


/************************** end of fail_loop ********************************/ 

