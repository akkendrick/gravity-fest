/*
***                           File strain.c                       ***
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

#define EXTERN extern
#include "stiff.h"
#include "strain.h"


real dotsh(
  real*		array,	/* array of nodal quantities to interpolate	*/
  int		i,	/* derivative component index number		*/
  int		j,	/* dof index number of quantity component	*/
  real*		sh,	/* el storage containing shape functions	*/
  ELEMENT_DATA*	el_pt,	/* ptr to current element struct		*/
  ELEMENT_INFO*	info,	/* ptr to element info struct			*/
  int		nen,	/* number of element nodes			*/
  int		ndim,	/* number of dimensions in array		*/
  int		st_flag,	/* flag to invoke split node adjustment		*/
  int  task_code   /* elastic, stepping or quake  */
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
  
void
  full_back(
  real*	stiff,		/* factored profile matrix			*/
  real*	rhs,		/* assembled rhs vector				*/
  int*	diag,		/* array of diagonal addresses			*/
  int	number_of_eqs	/* number of equations				*/
  );
  
void
factor(
  real*	stiff,		/* assembled profile matrix			*/
  int*	diag,		/* array of diagonal addresses			*/
  int	number_of_eqs	/* number of equations				*/
  );
  

static
real	sh_temp[256], det_temp[9];

static
int	el_diag[24] = {0, 2, 5, 9, 14, 20, 27, 35, 44, 54, 65, 77,
	90, 104, 119, 135, 152, 170, 189, 209, 230, 252, 275, 299};




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
                )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine form_stress
** form_stress computes the stress or time-step stress increment
** based on the FE solution.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      ELEMENT_DATA  *el_pt ;
      ELEMENT_MAT   *mat_pt ;
      ELEMENT_INFO  *info ;
      real  mu , lambda , dil , e1 , e2 , e3 , e4 , e5 , e6, dstr_energy, w_det , contrib ,
            product , trace ;
      real  dstr[6] , *str , *shape_ptr , *det_ptr , *sh ,
             *dbar , *bta, *wt , *strain_list ;
      int    i , numel , ipt , nel , nsd , nsdp1 ,
             nint , nen , ndof , str_size , dbar_size ,
             save_shape , type , this_code ;
     
      nsd = fe_sys.nsd ;
      nsdp1 = nsd + 1 ;
      ndof = fe_sys.ndof ;
      numel = (grp_ptr->el_info)->numel ;
      nen = (grp_ptr->el_info)->nen ;
      type = (grp_ptr->el_info)->type ;
      str_size = (grp_ptr->el_info)->stress_size ;
      dbar_size = (grp_ptr->el_info)->dbar_size ;
      save_shape = fe_sys.save_shape;
      mat_pt = grp_ptr->el_mat ;
      info = grp_ptr->el_info ;
      strain_list = info->strain_list ;

      switch((grp_ptr->el_info)->type)
	{
	case BILIN:
	  wt = wt_2x2;
	  break;
	case SERENDIP:
	  wt = wt_2x2x2;
	  break;
	case TET:
	  wt = wt_tet;
	  break;
	case TRILIN:
	  wt = wt_3x3;
	  break;
	}

      switch((grp_ptr->el_info)->type)
          {
           case BILIN:
           case SERENDIP:
           case TET:
           case TRILIN:
               for ( nel = 0 ; nel < numel ; nel++)
                   {
                    el_pt = grp_ptr->el_data + nel ;
                    dstr_energy = 0.0e0;
                    nint = el_pt->nint ;
                    this_code = code ;
                    if( !(*(mat_pt->plastic + el_pt->mat)) && (this_code != QUAKEEVT))
                       this_code = ELASTIC ;
                    mu = *(mat_pt->mu + (el_pt->mat)) ;
                    lambda = *(mat_pt->lambda + (el_pt->mat)) ;
	            	shape_ptr = (save_shape)? el_pt->shape: sh_temp ;
	            	det_ptr = (save_shape)? el_pt->deter: det_temp ;
	            	
	     /* filter out any duo elements for special handling */
	                if(info->type == TET && el_pt->degen == 1)
	                   {
	                    duostress(el_pt , mat_pt , info) ;
	                    continue ;
	                   }

                    if( !save_shape )
                         {
                          shape( el_pt , grp_ptr->el_info , shape_ptr , det_ptr , type ,
                                 nen , nint , 1 ) ;
                         }

                      el_pt->str_energy = 0.0;
                      *(strain_list+nel) = ZERO ;   /* initialize */
                      for (ipt = 0 ; ipt < nint ; ipt++)
                           {
                            sh = shape_ptr + ipt*nen*nsdp1 ;
                            w_det = wt[ipt] * det_ptr[ipt];
                            str = (el_pt -> stress) + ipt*str_size ;
                            bta    = (el_pt->bta) + ipt*str_size    ;
                            dbar = (el_pt -> dbar) + ipt*dbar_size ;
/*                          dilat = (el_pt -> dilat) + ipt ;          */

                            switch( ndof )
                               {
                                case 1:
                      /* no pressures in the 1-dimensional case */
                                   e1 = dotsh( global.del_displ ,
                                           0 , 0 , sh , el_pt , info , nen , ndof , 1 , code ) ;
                                   e2 = dotsh( global.del_displ ,
                                           1 , 0 , sh , el_pt , info , nen , ndof , 1 , code ) ;

                                   if( ELASTIC == this_code )
                                      {
                                       dstr[0] = mu * e1 ;
                                       dstr[1] = mu * e2 ;
                                      }
                 
                                    else if( VISCO == this_code )
                                      {
                                       dstr[0] = e1 - *(bta) ;
                                       dstr[1] = e2 - *(bta+1) ;
                                       full_back(dbar,dstr,el_diag, str_size) ;
                                      }
                                  break ;
                                   
                                case 2:
                                   e1 = dotsh( global.del_displ ,
                                            0 , 0 , sh , el_pt , info , nen , ndof , 1 , code ) ;
                                   e2 = dotsh( global.del_displ ,
                                            1 , 1 , sh , el_pt , info , nen , ndof , 1 , code ) ;
                                   e3 = dotsh( global.del_displ ,
                                            1 , 0 , sh , el_pt , info , nen , ndof , 1 , code ) +
                                        dotsh( global.del_displ ,
                                            0 , 1 , sh , el_pt , info , nen , ndof , 1 , code ) ;
           /* not used             ed = dotsh( global.displ ,
                                            0 , 0 , sh , el_pt , info , nen , ndof , 1 ) +
                                        dotsh( global.displ ,
                                            1 , 1 , sh , el_pt , info , nen , ndof , 1 ) ;
                                   *dilat = ed * lambda ;
             ...    */
                                   if( ELASTIC == this_code )
                                      {
                                       dil = lambda * (e1 + e2) ;
                                       dstr[0] = dil + TWO * mu * e1 ;
                                       dstr[1] = dil + TWO * mu * e2 ;
                                       dstr[2] = dil ; /* out of plane component */
                                       dstr[3] = mu * e3 ;
                                      }
                 
                                    else if( VISCO == this_code )
                                      {
                                       dstr[0] = e1 - *(bta) ;
                                       dstr[1] = e2 - *(bta+1) ;
                                       dstr[2] = ZERO - *(bta+2) ; /* out of plane component */
                                       dstr[3] = e3 - *(bta+3) ;
                                       full_back(dbar,dstr,el_diag, str_size) ;
                                      }
                                   break ;
                                   
                                case 3:
                                   if(2 == nsd) 
                                   {
                                      e1 = dotsh( global.del_displ ,
                                              0 , 0 , sh , el_pt , info , nen , ndof , 1 , code ) ;
                                      e2 = dotsh( global.del_displ ,
                                              1 , 1 , sh , el_pt , info , nen , ndof , 1 , code ) ;
                                      e3 = ZERO ;
                                      
                                      e4 = dotsh( global.del_displ ,
                                              1 , 0 , sh , el_pt , info , nen , ndof , 1 , code ) +
                                           dotsh( global.del_displ ,
                                              0 , 1 , sh , el_pt , info , nen , ndof , 1 , code ) ;
                                      e5 = dotsh( global.del_displ ,
                                           0 , 2 , sh , el_pt , info , nen , ndof , 1 , code ) ;
                                      e6 = dotsh( global.del_displ ,
                                           1 , 2 , sh , el_pt , info , nen , ndof , 1 , code ) ;
                /* not used           ed = dotsh( global.displ ,
                                            0 , 0 , sh , el_pt , info , nen , ndof , 1 ) +
                                           dotsh( global.displ ,
                                            1 , 1 , sh , el_pt , info , nen , ndof , 1 ) ;
                                      *dilat = ed * lambda ;
                  ...    */
                                      if( ELASTIC == this_code )
                                         {
                                          dil = lambda * (e1 + e2) ;
                                          dstr[0] = dil + TWO * mu * e1 ;
                                          dstr[1] = dil + TWO * mu * e2 ;
                                          dstr[2] = dil ;
                                          dstr[3] = mu * e4 ;
                                          dstr[4] = mu * e5 ;
                                          dstr[5] = mu * e6 ;
                                         }
                    
                                       else if( VISCO == this_code )
                                         {
                                          dstr[0] = e1 - *(bta) ;
                                          dstr[1] = e2 - *(bta+1) ;
                                          dstr[2] = e3 - *(bta+2) ;
                                          dstr[3] = e4 - *(bta+3) ;
                                          dstr[4] = e5 - *(bta+4) ;
                                          dstr[5] = e6 - *(bta+5) ;
                                          full_back(dbar,dstr,el_diag, str_size) ;
                                          }
                                   } else if (3 == nsd){

                                      e1 = dotsh( global.del_displ ,
                                              0 , 0 , sh , el_pt , info , nen , ndof , 1 , code ) ;
                                      e2 = dotsh( global.del_displ ,
                                              1 , 1 , sh , el_pt , info , nen , ndof , 1 , code ) ;
                                      e3 = dotsh( global.del_displ ,
                                              2 , 2 , sh , el_pt , info , nen , ndof , 1 , code ) ;
                                      
                                      e4 = dotsh( global.del_displ ,
                                              1 , 0 , sh , el_pt , info , nen , ndof , 1 , code ) +
                                           dotsh( global.del_displ ,
                                              0 , 1 , sh , el_pt , info , nen , ndof , 1 , code ) ;
                                      e5 = dotsh( global.del_displ ,
                                              2 , 0 , sh , el_pt , info , nen , ndof , 1 , code ) +
                                           dotsh( global.del_displ ,
                                              0 , 2 , sh , el_pt , info , nen , ndof , 1 , code ) ;
                                      e6 = dotsh( global.del_displ ,
                                              2 , 1 , sh , el_pt , info , nen , ndof , 1 , code ) +
                                           dotsh( global.del_displ ,
                                              1,  2 , sh , el_pt , info , nen , ndof , 1 , code ) ;
/*
   In addition, we need to calculate the total (as opposed to incremental) pressure stress
   that would exist if the medium had zero rigidity; this is to be used as a goal toward
   which the viscoplastic corrections will drive the dilatational strain...
*/
      /*   but this is now not used...
                                      ed = dotsh( global.displ ,
                                              0 , 0 , sh , el_pt , info , nen , ndof , 1 ) +
                                           dotsh( global.displ ,
                                              1,  1 , sh , el_pt , info , nen , ndof , 1 ) +
                                           dotsh( global.displ ,
                                              2,  2 , sh , el_pt , info , nen , ndof , 1 ) ;
                                      *dilat = ed * lambda ;
                  ...    */
                                      if( ELASTIC == this_code || QUAKEEVT == this_code )
                                         {
                                          dil = lambda * (e1 + e2 + e3) ;
                                          dstr[0] = dil + TWO * mu * e1 ;
                                          dstr[1] = dil + TWO * mu * e2 ;
                                          dstr[2] = dil + TWO * mu * e3;
                                          dstr[3] = mu * e4 ;
                                          dstr[4] = mu * e5 ;
                                          dstr[5] = mu * e6 ;
/* experimental strain energy computation */
                                          dstr_energy = (e1*dstr[0] + e2*dstr[1] + e3*dstr[2] +
                                                         e4*dstr[3] + e5*dstr[4] + e6*dstr[5]);
/* this next correction (rhs) is the element volume. */
                                          dstr_energy *= w_det;
                                          el_pt->str_energy += dstr_energy;
                                         }
                    
                                       else if( VISCO == this_code )
                                         {
                                          dstr[0] = e1 - *(bta) ;
                                          dstr[1] = e2 - *(bta+1) ;
                                          dstr[2] = e3 - *(bta+2) ;
                                          dstr[3] = e4 - *(bta+3) ;
                                          dstr[4] = e5 - *(bta+4) ;
                                          dstr[5] = e6 - *(bta+5) ;
                                          full_back(dbar,dstr,el_diag, str_size) ;
                                          }
                                   }
                                   break ;
                                   
                               }
                            for(i=0 ; i<str_size ; i++)
                               {
                                *(str + i) += dstr[i] ;
                               }
                               
                            if( QUAKEEVT == this_code ) /* global strain energy tally */
                               {
                                product = str[0]*str[0]+str[1]*str[1]+str[2]*str[2]+str[3]*str[3]+
                                          str[4]*str[4]+str[5]*str[5] ;
                                trace = str[0]+str[1]+str[2] ;
                                contrib = w_det*( PT25*product/mu - 
                                  PT25*lambda*trace*trace/(mu*(THREE*lambda+TWO*mu)) ) ;
                                global.elas_energy += contrib ;
                               }
                            
                           }
                   }
               break ;
               
          }
    }

/************************** end of form_stress *******************************/ 



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
              )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine form_beta
** form_beta computes beta, the viscoplastic strain rate for
** this element and time step.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

       {
        real  *stress  ,  *bta  , *visc , *strain_list ;
        real  fn1  ,  sxx  , syy  ,  sxy  , sigma  , ev , ep , mu_el ,
               sxz , syz , s3 , s4 , s5 , del1 , del2 , del3 , szz ;
        int    ipt , ndof , nint , str_size , plastic , iel ;
	    real  step_strain ;
        
        ndof = fe_sys.ndof ;
        nint = el_pt->nint ;
        str_size = info->stress_size ;
        plastic = *(mat_pt->plastic + (el_pt->mat)) ;
        visc  = mat_pt->visc + 2*(el_pt->mat) ;
        mu_el   = *(mat_pt->mu + el_pt->mat) ;
        strain_list = info->strain_list ;
            
        fn1 = *(visc + 1) - ONE ;
        ev  = time_data.dt / (*visc) ;
        ep = ev ;

        for (ipt = 0 ; ipt < nint ; ipt++)
           {
            stress = (el_pt->stress) + ipt*str_size ;
            bta    = (el_pt->bta) + ipt*str_size    ;
/*          dilat  = (el_pt->dilat) + ipt           ;  */
    
            switch( ndof )
               {
                case 1:
                    sxz =  *(stress) ;
                    syz =  *(stress + 1 ) ;
                    
                    if ( plastic > 1 )
                       {
                        sigma = sxz * sxz + syz * syz ;
                        ev = pow(sigma,PT5*fn1) * time_data.dt / (*visc) ;
                       }
                    
                    *bta        = ev * sxz ;
                    *(bta + 1 ) = ev * syz ;
                    break ;
    
                case 2:
                    sxx =  *(stress) ;
                    syy =  *(stress + 1 ) ;
                    szz =  *(stress + 2 ) ;
                    del1  =  sxx - syy ;
                    del2  =  sxx - szz ;
                    del3  =  syy - szz ;
                    sxy =  *(stress + 3 ) ;
/*                  press = (sxx+syy+szz)/THREE - *dilat ;  */
                    
                    if ( plastic > 1 )
                       {
                        sigma = (del1*del1 + del2*del2 + del3*del3)/SIX
                                 + sxy*sxy ;
                        ev = pow(sigma,PT5*fn1) * time_data.dt / (*visc) ;
/*                      ep = pow(fabs(press),fn1) * time_data.dt / (*visc) ;  */
                       }
                    
                    *bta        = ev * ( sxx/THREE - (syy+szz)/SIX ) /*  + (ep * press)  */ ;
                    *(bta + 1 ) = ev * ( syy/THREE - (sxx+szz)/SIX ) /*  + (ep * press)  */ ;
                    *(bta + 2 ) = ev * ( szz/THREE - (sxx+syy)/SIX ) /*  + (ep * press)  */ ;
                    *(bta + 3 ) = ev * sxy ;
                    break ;
    
                case 3:
                    sxx =  *(stress) ;
                    syy =  *(stress + 1 ) ;
                    szz =  *(stress + 2 ) ;
                    s3  =  *(stress + 3 ) ;
                    s4  =  *(stress + 4 ) ;
                    s5  =  *(stress + 5 ) ;
                    del1  =  sxx - syy ;
                    del2  =  sxx - szz ;
                    del3  =  syy - szz ;
/* calculate the difference between the current pressure and the zero rigidity "goal" */
/*                  press = (sxx+syy+szz)/THREE - *dilat ;     */
                    
                    if ( plastic > 1 )
                       {
                        sigma = (del1*del1 + del2*del2 + del3*del3)/SIX
                                 + s3*s3 + s4*s4 + s5*s5 ;
                        ev = pow(sigma,PT5*fn1) * time_data.dt / (*visc) ;
/*                      ep = pow(fabs(press),fn1) * time_data.dt / (*visc) ;   */
                       }
                    
                    *bta        = ev * ( sxx/THREE - (syy+szz)/SIX ) /*  + (ep * press)  */ ;
                    *(bta + 1 ) = ev * ( syy/THREE - (sxx+szz)/SIX ) /*  + (ep * press)  */ ;
                    *(bta + 2 ) = ev * ( szz/THREE - (sxx+syy)/SIX ) /*  + (ep * press)  */ ;
                    *(bta + 3 ) = ev * s3 ;
                    *(bta + 4 ) = ev * s4 ;
                    *(bta + 5 ) = ev * s5 ;
                    break ;
               }

	       if (plastic > 1)
	       {
	       /* determine how large the step in strain is, in case the time step size
	          needs to be reduced (only needed for non-Newtonian problems) */
	          
            sigma = sqrt( sigma ) ;
		 switch (ndof)
		 {
		   case 1:  
		     step_strain = sqrt (*bta * *bta + *(bta+1) * *(bta+1)) ;
		     break ;
		   case 2:  
		     step_strain = sqrt (*bta * *bta + *(bta+1) * *(bta+1) 
					 + *(bta+2) * *(bta+2)) ;
		     break ;
		   case 3:  
		     step_strain = sqrt (*bta * *bta + *(bta+1) * *(bta+1)
				+ *(bta+2) * *(bta+2) + *(bta+3) * *(bta+3) 
				+ *(bta+4) * *(bta+4) + *(bta+5) * *(bta+5)) ;
		     break ;
		 }
                 step_strain *= (mu_el/sigma) ;
                 
         /* No longer perform the test here... instead store it  */
		/* time_data.max_strain = MAX(time_data.max_strain, step_strain) ; */
		
		         iel = el_pt->nel - 1 ;  /* zero-based index */
		         if( *(strain_list+iel) < step_strain )
		             *(strain_list+iel) = step_strain ;
		        
	       }

           }
       }
           
/************************** end of form_beta *********************************/ 



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
              )

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

       {
        real  *stress  ,  *visc , *dbar ;
        real  fn1  ,  sxx  , syy  ,  sxy  , sigma  , ev , 
               sxz , syz , s3 , s4 , s5 , del1 , del2 , del3 , szz , 
               a , t1 , t2 , t3 , sx , sy , sz , ev_sig , en ;
               
        real  lam_el , mu_el , lam_2mu , det ;
        real  c1 , c2 , c3 , lame1 , lame2 , mu_cube , cs ;
        int    ndof , nint , plastic , str_size , dbar_size ;
        int    ipt ;
        
        
        ndof = fe_sys.ndof ;
        nint = el_pt->nint ;
        plastic = *(mat_pt->plastic + el_pt->mat) ;
        str_size = info->stress_size ;
        dbar_size = info->dbar_size ;
        lam_el  = *(mat_pt->lambda + el_pt->mat) ;
        mu_el   = *(mat_pt->mu + el_pt->mat) ;
        lam_2mu = lam_el + TWO * mu_el ;
        visc  = (mat_pt->visc + 2*(el_pt->mat)) ;
            
        fn1 = *(visc + 1) - ONE ;
        ev  = time_data.alpha_delt / (*visc) ;
/*      ep = ev ;   */
        
        switch( ndof )
           {
            case 2:
              det = 8.0*pow(mu_el,FOUR) + 12.0*lam_el*pow(mu_el,THREE) ;
              c1 = ( FOUR * mu_el * mu_el * (mu_el+lam_el) ) / det ;
              c2 = ( -TWO * mu_el * mu_el * lam_el ) / det ;
              break ;
              
            case 3:
              lame1 = lam_2mu * lam_2mu  -  lam_el * lam_el  ;
              lame2 = lam_el * lam_el  -  lam_el * lam_2mu  ;
              mu_cube = pow(mu_el,THREE) ;
              cs = lam_2mu * lame1 + TWO * lam_el * lame2 ;
              det = cs * mu_cube ;
              c1 = lame1 * mu_cube / det ;
              c2 = lame2 * mu_cube / det ;
              c3 = cs * mu_el * mu_el / det ;
              break ;
           }

        for (ipt = 0 ; ipt < nint ; ipt++)
           {
            stress = (el_pt->stress) + ipt*str_size ;
            dbar = (el_pt -> dbar) + ipt*dbar_size ;
/*          dilat = (el_pt -> dilat) + ipt ;     */
   
            switch( ndof )
               {
                case 1:
                    if( 2 == plastic )
                       {
                        sxz =  *(stress) ;
                        syz =  *(stress + 1 ) ;
                        sigma = sxz * sxz + syz * syz ;
                        ev_sig = pow(sigma,PT5*fn1) * ev ;
                        a = fn1/sigma ;
                        *(dbar) = (ONE / mu_el) + ev_sig *
                                              (1 + a*sxz*sxz) ;
                        *(dbar+1) = ev_sig*a*sxz*syz ;
                        *(dbar+2) = (ONE / mu_el) + ev_sig *
                                              (1 + a*syz*syz) ;
                       }
                    else
                       {
                        *(dbar) = (ONE / mu_el) + ev ;
                        *(dbar+1) = ZERO ;
                        *(dbar+2) = (ONE / mu_el) + ev ;
                       }
                    break ;
    
                case 2:

                    if( 2 == plastic )
                       {
                        sxx =  *(stress) ;
                        syy =  *(stress + 1 ) ;
                        szz =  *(stress + 2 ) ;
                        sxy =  *(stress + 3 ) ;
                        del1  =  sxx - syy ;
                        del2  =  sxx - szz ;
                        del3  =  syy - szz ;
                        sigma = (del1*del1 + del2*del2 + del3*del3)/SIX
                                 + sxy*sxy ;
/*                      press = (sxx+syy+szz)/THREE - *dilat ;     */
                        ev_sig = pow(sigma,PT5*fn1) * ev ;
/*                      ep = pow(fabs(press),fn1) * *(visc+1) * ep ;   */
                        en = ev_sig * fn1 * PT25 / sigma ;
                        sx = (TWO*sxx - syy - szz) / THREE ;
                        sy = (TWO*syy - sxx - szz) / THREE ;
                        sz = (TWO*szz - sxx - syy) / THREE ;
                        t3 = TWO * sxy ;
                        *(dbar)       = c1 + ev_sig / THREE + en*sx*sx  /* + ep/THREE */ ;
                        *(dbar +  1 ) = c2 - ev_sig / SIX + en*sx*sy    /* + ep/THREE */ ;
                        *(dbar +  2 ) = c1 + ev_sig / THREE + en*sy*sy  /* + ep/THREE */ ;
                        *(dbar +  3 ) = c2 - ev_sig / SIX + en*sx*sz    /* + ep/THREE */ ;
                        *(dbar +  4 ) = c2 - ev_sig / SIX + en*sy*sz    /* + ep/THREE */ ;
                        *(dbar +  5 ) = c1 + ev_sig / THREE + en*sz*sz  /* + ep/THREE */ ;
                        *(dbar +  6 ) = en*sx*t3 ;
                        *(dbar +  7 ) = en*sy*t3 ;
                        *(dbar +  8 ) = en*sz*t3 ;
                        *(dbar +  9 ) = (ONE/mu_el) + ev_sig  + en*t3*t3 ;
                       }
                    else
                       {
                        *(dbar)       = c1 + ev / THREE  /* + ep/THREE */ ;
                        *(dbar +  1 ) = c2 - ev / SIX    /* + ep/THREE */ ;
                        *(dbar +  2 ) = c1 + ev / THREE  /* + ep/THREE */ ;
                        *(dbar +  3 ) = c2 - ev / SIX    /* + ep/THREE */ ;
                        *(dbar +  4 ) = c2 - ev / SIX    /* + ep/THREE */ ;
                        *(dbar +  5 ) = c1 + ev / THREE  /* + ep/THREE */ ;
                        *(dbar +  6 ) = ZERO ;
                        *(dbar +  7 ) = ZERO ;
                        *(dbar +  8 ) = ZERO ;
                        *(dbar +  9 ) = (ONE/mu_el) + ev ;
                       }
                    break ;
    
                case 3:

                    if( 2 == plastic )
                       {
                        sxx =  *(stress) ;
                        syy =  *(stress + 1 ) ;
                        szz =  *(stress + 2 ) ;
                        s3  =  *(stress + 3 ) ;
                        s4  =  *(stress + 4 ) ;
                        s5  =  *(stress + 5 ) ;
                        del1  =  sxx - syy ;
                        del2  =  sxx - szz ;
                        del3  =  syy - szz ;
                        sigma = (del1*del1 + del2*del2 + del3*del3)/SIX
                                 + s3*s3 + s4*s4 + s5*s5 ;
/* calculate the difference between the current pressure and the zero rigidity "goal" */
/*                      press = (sxx+syy+szz)/THREE - *dilat ;       */
                        ev_sig = pow(sigma,PT5*fn1) * ev ;
/*                      ep = pow(fabs(press),fn1) * *(visc+1) * ep ;   */
                        en = ev_sig * fn1 * PT25 / sigma ;
                        sx = (TWO*sxx - syy - szz) / THREE ;
                        sy = (TWO*syy - sxx - szz) / THREE ;
                        sz = (TWO*szz - sxx - syy) / THREE ;
                        t1 = TWO * s5 ;
                        t2 = TWO * s4 ;
                        t3 = TWO * s3 ;

                        *(dbar)       = c1 + ev_sig / THREE + en*sx*sx   /* + ep/THREE */ ;
                        *(dbar +  1 ) = c2 - ev_sig / SIX + en*sx*sy     /* + ep/THREE */ ;
                        *(dbar +  2 ) = c1 + ev_sig / THREE + en*sy*sy   /* + ep/THREE */ ;
                        *(dbar +  3 ) = c2 - ev_sig / SIX + en*sx*sz     /* + ep/THREE */ ;
                        *(dbar +  4 ) = c2 - ev_sig / SIX + en*sy*sz     /* + ep/THREE */ ;
                        *(dbar +  5 ) = c1 + ev_sig / THREE + en*sz*sz   /* + ep/THREE */ ;
                        *(dbar +  6 ) = en*sx*t3 ;
                        *(dbar +  7 ) = en*sy*t3 ;
                        *(dbar +  8 ) = en*sz*t3 ;
                        *(dbar +  9 ) = c3 + ev_sig  + en*t3*t3 ;
                        *(dbar + 10 ) = en*sx*t2 ;
                        *(dbar + 11 ) = en*sy*t2 ;
                        *(dbar + 12 ) = en*sz*t2 ;
                        *(dbar + 13 ) = en*t3*t2 ;
                        *(dbar + 14 ) = c3 + ev_sig  + en*t2*t2 ;
                        *(dbar + 15 ) = en*sx*t1 ;
                        *(dbar + 16 ) = en*sy*t1 ;
                        *(dbar + 17 ) = en*sz*t1 ;
                        *(dbar + 18 ) = en*t3*t1 ;
                        *(dbar + 19 ) = en*t1*t2 ;
                        *(dbar + 20 ) = c3 + ev_sig  + en*t1*t1 ;
                       }
                    else
                       {
                        *(dbar)       = c1 + ev / THREE   /* + ep/THREE */ ;
                        *(dbar +  1 ) = c2 - ev / SIX     /* + ep/THREE */ ;
                        *(dbar +  2 ) = c1 + ev / THREE   /* + ep/THREE */ ;
                        *(dbar +  3 ) = c2 - ev / SIX     /* + ep/THREE */ ;
                        *(dbar +  4 ) = c2 - ev / SIX     /* + ep/THREE */ ;
                        *(dbar +  5 ) = c1 + ev / THREE   /* + ep/THREE */ ;
                        *(dbar +  6 ) = ZERO ;
                        *(dbar +  7 ) = ZERO ;
                        *(dbar +  8 ) = ZERO ;
                        *(dbar +  9 ) = c3 + ev ;
                        *(dbar + 10 ) = ZERO ;
                        *(dbar + 11 ) = ZERO ;
                        *(dbar + 12 ) = ZERO ;
                        *(dbar + 13 ) = ZERO ;
                        *(dbar + 14 ) = c3 + ev ;
                        *(dbar + 15 ) = ZERO ;
                        *(dbar + 16 ) = ZERO ;
                        *(dbar + 17 ) = ZERO ;
                        *(dbar + 18 ) = ZERO ;
                        *(dbar + 19 ) = ZERO ;
                        *(dbar + 20 ) = c3 + ev ;
                       }
                    break ;
               }
            factor(dbar , el_diag , str_size) ;
           }
       }

/************************** end of form_dbar *********************************/ 



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
              )

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

       {
        int  ndof , option , inode1 , inode2 ;
        real  *str  ,  *cvect ;
        real  cx , cy , cz , ux1 , ux2 , uy1 , uy2 , uz1 , uz2 , kspring , len ;
        
        ndof = fe_sys.ndof ;
        len  = *(mat_pt->lambda + el_pt->mat) ;
        kspring   = *(mat_pt->mu + el_pt->mat) ;
        
        cvect = mat_pt->bforce + ndof*(el_pt->mat) ; /* grav is used as holding space for cvect */
        if( (int) ( *(mat_pt->visc + 2*(el_pt->mat)) ) ) option = 1 ;
        else option = 0 ;
         /* visc == zero for truss ; visc non-zero for spring */
         
        str = (el_pt -> stress) ;
         
         switch( option )
            {
             case 0:
                cx = *cvect ;
                cy = *(cvect+1) ;
                cz = *(cvect+2) ;
                inode1 = *(el_pt->ien    ) - 1 ;  /* zero-based node index */
                inode2 = *(el_pt->ien + 1) - 1 ;
                ux1 = *(global.del_displ+inode1*ndof) ;
                ux2 = *(global.del_displ+inode2*ndof) ;
                uy1 = *(global.del_displ+inode1*ndof+1) ;
                uy2 = *(global.del_displ+inode2*ndof+1) ;
                uz1 = *(global.del_displ+inode1*ndof+2) ;
                uz2 = *(global.del_displ+inode2*ndof+2) ;
    /* for a roller, there is only a single scalar stress component, str[0] */
                *str += (kspring/len)*(cx*(ux2-ux1)+cy*(uy2-uy1)+cz*(uz2-uz1)) ;
                break ;

             case 1:
                inode1 = *(el_pt->ien    ) - 1 ;  /* zero-based node index */
                inode2 = *(el_pt->ien + 1) - 1 ;
                ux1 = *(global.del_displ+inode1*ndof) ;
                ux2 = *(global.del_displ+inode2*ndof) ;
                uy1 = *(global.del_displ+inode1*ndof+1) ;
                uy2 = *(global.del_displ+inode2*ndof+1) ;
                uz1 = *(global.del_displ+inode1*ndof+2) ;
                uz2 = *(global.del_displ+inode2*ndof+2) ;
    /* for a spring, there are 3 stress components, str[0] thru str[2] */
                *str     += (kspring/len)*(ux2-ux1) ;
                *(str+1) += (kspring/len)*(uy2-uy1) ;
                *(str+2) += (kspring/len)*(uz2-uz1) ;
                break ;
            }
  /* for both kinds of duo, we'll use str[3] thru str[5] to hold the displacement offset */
        *(str+3)  +=  ux2-ux1 ;
        *(str+4)  +=  uy2-uy1 ;
        *(str+5)  +=  uz2-uz1 ;
       }

/************************** end of duostress *********************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   find_max_strain   ***************
*/

/*  INTERFACE:   */
     void
     find_max_strain(
                     GROUP   *grp_ptr    /* pointer to current element group */
                    )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine find_max_strain
** find_max_strain sorts the viscoplastic strains in all elements
** and stores an adjusted maximum value for use in time step computation
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      ELEMENT_INFO  *info ;
      real  *strain_list ;
      int    i , numel ;
     
      info = grp_ptr->el_info ;
      strain_list = info->strain_list ;
      numel = info->numel ;
      
      /* First sort the saved strains in ascending order  */
      qsort( strain_list , numel , sizeof(real) , real_cmp ) ;
      
      /* Now select the value that is STRAIN_PCT% from the top and store it  */
      i = (int) ( numel*(ONE - (STRAIN_PCT/100.0)) ) - 1 ;
      time_data.max_strain = *(strain_list + i) ;      

      }

/************************** end of find_max_strain *********************************/ 



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
             )

/*  RETURN VALUE:  -1, 0, or 1 */
/*  DESCRIPTION:   */
/*
** Routine real_cmp
** comparison function to sort real numbers with qsort
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      real  first  = *(real *)first_arg ;
      real  second = *(real *)second_arg ;
      if( first < second )
         {
          return -1 ;
         }
      else if( first == second )
         {
          return 0 ;
         }
      else
         {
          return 1 ;
         }
     }

/************************** end of real_cmp *********************************/ 



