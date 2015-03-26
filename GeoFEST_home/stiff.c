/*
***                           File stiff.c                        ***
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

#define EXTERN extern
#include "stiff.h"

#define SH(i,n) *(sh + nsdp1*(n) + (i))
#define ES(i,j) *(es + nee*(j) + (i))

#define PSH_BILIN(i,enode,point)\
 *(psh_b + 4*nsdp1*(enode) + nsdp1*(point) + (i))

#define PSH_SEREN(i,enode,point,rule)\
 *(psh_s + 18*nsdp1*(enode) + 9*nsdp1*(rule) + nsdp1*(point) + (i))

#define PSH_TET(i,enode,point)\
 *(psh_t + 4*nsdp1*(enode) + nsdp1*(point) + (i))

#define PSH_TRILIN6(i,enode,point)\
 *(psh_h6 + 6*nsdp1*(enode) + nsdp1*(point) + (i))

#define PSH_TRILIN8(i,enode,point)\
 *(psh_h8 + 8*nsdp1*(enode) + nsdp1*(point) + (i))


static real
  psh_b[48] ,
   /* parent space shape functions for 4-node bilinear quad */
  psh_s[432]   ,
   /* parent space shape functions for 8-node serendipity quad */
  psh_t[64]    ,
   /* parent space shape functions for 4-node tet */
  psh_h6[192]    ,
   /* parent space shape functions for 8-node hexahedron (6-pt rule) */
  psh_h8[256]    ;
   /* parent space shape functions for degen 8-node hexahedron (8-pt rule) */

static real
   el_stiff[576] , el_rhs[24] , temp[24] , sh_temp[256] , det_temp[9] ;

static int  el_diag[24] =
            { 0 , 2 , 5 , 9 , 14 , 20 , 27 , 35 , 44 , 54 , 65 , 77 ,
             90 , 104 , 119 , 135 , 152 , 170 , 189 , 209 , 230 ,
             252 , 275 , 299 } ;
             
void completion(void);

/*
real  dotsh();
real adfldp();
*/
void
addstiff(
  PROFILE*	a,	/* ptr to assembled profile stiffness storage	*/
  real*		es,	/* ptr to element stiffness storage		*/
  ELEMENT_INFO*	info,	/* ptr to element info struct			*/
  ELEMENT_DATA* el_pt	/* ptr to current element struct		*/
  );

void
addfor (
  real*		dest,	/* destination assembled rhs vector		*/
  real*		rhs,	/* source element rhs vector			*/
  ELEMENT_INFO*	info,	/* ptr to element info struct			*/
  ELEMENT_DATA*	el_pt	/* ptr to current element struct		*/
  );

void
form_dbar(
  ELEMENT_DATA* el_pt,	/* ptr to current element struct		*/
  ELEMENT_MAT*  mat_pt,	/* ptr to el material prop struct		*/
  ELEMENT_INFO* info	/* ptr to element info struct			*/
  );

void
full_back(
  real* stiff,		/* factored profile matrix			*/
  real* rhs,		/* assembled rhs vector				*/
  int*  diag,		/* array of diagonal addresses			*/
  int   number_of_eqs	/* number of equations				*/
  );

void
form_beta(
  ELEMENT_DATA*	el_pt,	/* ptr to current element struct		*/
  ELEMENT_MAT*	mat_pt,	/* ptr to el material prop struct		*/
  ELEMENT_INFO*	info	/* ptr to element info struct			*/
  );

void
move_real(
  real*	from,	/* source array						*/
  real*	to,	/* destination array					*/
  int	n	/* number of entries					*/
  );




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
               )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine form_stiff
** form_stiff computes the element-wise stiffness array for the elt group.
** "code" indicates if this is for an elastic problem or a VE step.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      int     numel , nel , nsd , nsdp1 ;

      ELEMENT_DATA   *el_pt ;


      nsd = fe_sys.nsd ;
      nsdp1 = nsd + 1 ;
      numel = (grp_ptr->el_info)->numel ;
      
      switch((grp_ptr->el_info)->type)
          {
           case TET:
           case BILIN:
           case SERENDIP:
           case TRILIN:
              for ( nel = 0 ; nel < numel ; nel++)
                   {
                    el_pt = grp_ptr->el_data + nel ;
                    if((grp_ptr->el_info)->type == TET && el_pt->degen == 1)
                      duo_form( code , el_stiff , grp_ptr->el_info ,
                              grp_ptr->el_mat, el_pt ) ;
                    else
                      lame_form( code , el_stiff , grp_ptr->el_info ,
                              grp_ptr->el_mat, el_pt ) ;

                    addstiff( &stiff , el_stiff ,
                              grp_ptr->el_info , el_pt ) ;
                   }
               break ;
               
          }
    }

/************************** end of form_stiff ********************************/ 



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
             )

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

     {
      int    n , numel , nel , ndof , numsuf , nfterm , side , itsuf , do_it ,
             nsd , nsdp1 , node , nee , i , numface , skip , cslip , rad_code ;
      int   *surf_list , *slip_list , *buoy_list , *slip_grp ;
      real  *surf_val , *slip_val , *trac , *slip , *up_ptr , rho_g , slip_amt ;
      ELEMENT_DATA   *el_pt ;
      BUOY_DATA      *buoy_ptr ;


      nsd = fe_sys.nsd ;
      nsdp1 = nsd + 1 ;
      numel = (grp_ptr->el_info)->numel ;
      numsuf = (grp_ptr->el_info)->numsuf ;
      nfterm = (grp_ptr->el_info)->nfterm ;
      nee    = (grp_ptr->el_info)->nee ;
      ndof = fe_sys.ndof ;

      switch((grp_ptr->el_info)->type)
          {
           case TET:
           case BILIN:
           case SERENDIP:
           case TRILIN:

               for ( nel = 0 ; nel < numel ; nel++)
                   {

                    el_pt = grp_ptr->el_data + nel ;

                    skip = 0 ;
                    if((grp_ptr->el_info)->type == TET && el_pt->degen == 1) skip = 1 ;
                        /* we don't want to do rhs calculations on duo elements! */
                    if(QUAKEEVT == code)  skip = 1 ;  /* no grav forces for quake events */
                    if(!skip)
                       {
                        force_form( code , el_rhs , grp_ptr->el_info ,
                              grp_ptr->el_mat, el_pt ) ;
                        addfor( force.full_rhs , el_rhs , grp_ptr->el_info ,
                            el_pt ) ;
                        if(ELASTIC == code) /* save forces for equilibrium correction */
                           {
                            addfor( force.net_external , el_rhs , grp_ptr->el_info ,
                               el_pt ) ;
                           }
                       } /* end if !skip */
                   }

                itsuf = time_data.itsuf ;      /* current surface traction index */
                do_it = 0 ;  /* don't do surface tractions unless it is time */
                if( time_data.elastic ) /* do  t=zero surface tractions */
                   {
                    surf_list = (grp_ptr->el_info) -> surf_list0 ;
                    surf_val = (grp_ptr->el_info) -> surf_trac0 ;
                    time_data.itsuf++ ;  /* increment surface traction index */
                    time_data.currtsuf = time_data.traction_time[time_data.itsuf] ;
                    printf("After elastic step, traction index = %d , currtsuf = %lg\n",time_data.itsuf,time_data.currtsuf) ;
                    printf("traction_time[0] = %lg\n",time_data.traction_time[0]) ;
                    printf("traction_time[1] = %lg\n",time_data.traction_time[1]) ;
                    printf("traction_time[2] = %lg\n",time_data.traction_time[2]) ;
                    printf("traction_time[3] = %lg\n",time_data.traction_time[3]) ;
                    printf("traction_time[4] = %lg\n",time_data.traction_time[4]) ;
                    printf("traction_time[5] = %lg\n",time_data.traction_time[5]) ;
                    do_it = 1 ;
                   }
                else if( time_data.time >= time_data.currtsuf && numsuf > 0) /* do  t>zero surface tractions */
                   {
                    printf("traction index = %d , currtsuf = %lg\n",itsuf,time_data.currtsuf) ;
/*                  surf_list = (grp_ptr->el_info) -> surf_list[itsuf] ;  */
/*                  surf_val = (grp_ptr->el_info) -> surf_trac[itsuf] ;  */
					switch( itsuf )
					   {
						case 0:
						  surf_list = (grp_ptr->el_info)->surf_list0 ;
						  surf_val = (grp_ptr->el_info)->surf_trac0 ;
						  break ;
		 
						 case 1:
						  surf_list = (grp_ptr->el_info)->surf_list1 ;
						  surf_val = (grp_ptr->el_info)->surf_trac1 ;
						  break ;
		 
						case 2:
						  surf_list = (grp_ptr->el_info)->surf_list2 ;
						  surf_val = (grp_ptr->el_info)->surf_trac2 ;
						  break ;
		 
						case 3:
						  surf_list = (grp_ptr->el_info)->surf_list3 ;
						  surf_val = (grp_ptr->el_info)->surf_trac3 ;
						  break ;
		 
						case 4:
						  surf_list = (grp_ptr->el_info)->surf_list4 ;
						  surf_val = (grp_ptr->el_info)->surf_trac4 ;
						  break ;
		 
						case 5:
						  surf_list = (grp_ptr->el_info)->surf_list5 ;
						  surf_val = (grp_ptr->el_info)->surf_trac5 ;
						  break ;           
					   }

                    time_data.itsuf++ ;  /* increment surface traction index */
                    time_data.currtsuf = time_data.traction_time[time_data.itsuf] ;
                    printf("Surface traction delta applied at t = %lg\n", time_data.time) ;
                    do_it = 1 ;
                   }

                if( do_it != 0 && numsuf > 0 ) /* do surface tractions */
                   {

/************************************************************************/
/*
surf_list[] is a numsuf-by-2-deep array of ints containing:
        [0] = global element number
        [1] = face number; 
    example: face 3 is the side comprised of nodes 1,2 and 4

surf_val[] is a numsuf-by-ndof-deep array of reals containing:
         [0] = x-component of surface traction
         [1] = y-component of surface traction
         [2] = z-component of surface traction
*/
/************************************************************************/

                    for(n=0 ; n<numsuf ; n++)
                       {
                        nel = *(surf_list + 2*n) ;
                        side = (int) *(surf_list + 2*n + 1) ;
                        trac = (surf_val + ndof*n) ;
                        el_pt = grp_ptr->el_data + nel ;
                        surf_form( code , el_rhs , grp_ptr->el_info ,
                                   side , trac , el_pt ) ;
                        addfor( force.full_rhs , el_rhs ,
                                grp_ptr->el_info , el_pt ) ;
                                
                    /* save forces for equilibrium correction */

                        addfor( force.net_external , el_rhs , grp_ptr->el_info ,
                            el_pt ) ;

                       }
                   }


                if( !time_data.elastic && (QUAKEEVT != code) ) /* do buoyancy tractions */
                   {
                    for(i=0 ; i < (grp_ptr->el_info)->numbuoy ; i++)
                       {
                        buoy_ptr = (grp_ptr->el_info)->buoy + i ;
                        buoy_list = buoy_ptr->buoy_list ;
						numface = buoy_ptr->numface ;
						rho_g = buoy_ptr->rho_g ;
						rad_code = buoy_ptr->rad_flag ;
		/*
						printf("Before...\n") ;
						printf("numface=%d\nrho_g=%g\nfirst el=%d\nfirst side=%d\n\n",
						       numface,rho_g,*buoy_list,*(buoy_list+1)) ;
		*/
						for(n=0 ; n<numface ; n++)
						   {
							nel = *(buoy_list + 2*n) ;
							side = (int) *(buoy_list + 2*n + 1) ;
							if (side <= 0) {
						printf("\tIAM:%d,BadSide:%d,nel:%d,numface=%d,first el=%d\nfirst side=%d\n\n",
						       iam, side, nel, numface,*buoy_list,*(buoy_list+1)) ;
							}
							up_ptr = buoy_ptr->upvec ;
							el_pt = grp_ptr->el_data + nel ;
							buoy_form( rad_code , el_rhs , grp_ptr->el_info ,
									   side , up_ptr , rho_g , el_pt ) ;
							addfor( force.full_rhs , el_rhs ,
									grp_ptr->el_info , el_pt ) ;
									
						/* save forces for equilibrium correction */
	
							addfor( force.net_external , el_rhs , grp_ptr->el_info ,
								el_pt ) ;
	
						   }
		/*
						printf("After...\n") ;
						printf("numface=%d\nrho_g=%g\nfirst el=%d\nfirst side=%d\n\n",
						       numface,rho_g,*buoy_list,*(buoy_list+1)) ;
		*/
   
                       }
                   }


                if( time_data.do_slip )
                   {
                    slip_list = (grp_ptr->el_info) -> slip_list ;
                    slip_val = (grp_ptr->el_info) -> slip_val ;
                    slip_grp = (grp_ptr->el_info) -> slip_grp ;
                   for(n=0 ; n<nfterm ; n++)
                       {
                        i = *(slip_grp+n) - 1 ;
                        slip_amt = (fltgrp_ptr+i)->q_amount ;
                        if((fltgrp_ptr+i)->due_now == 0 )
                           continue ; /* this fault group not moving at this time */
                           
                        cslip = 0 ;
                        if((fltgrp_ptr+i)->due_now == 2 )
                           cslip = 1 ; /* this fault group is continuously slipping */
                              
                        nel = (int) *(slip_list + 2*n) ;
                        node = (int) *(slip_list + 2*n + 1) ;
                        slip = (slip_val + ndof*n) ;
                        el_pt = grp_ptr->el_data + nel ;
                        form_slip( code , cslip , el_rhs , el_stiff , grp_ptr->el_info ,
                                   nel , node , slip , slip_amt , el_pt , grp_ptr->el_mat ) ;
                        addfor( force.full_rhs , el_rhs ,
                                grp_ptr->el_info , el_pt ) ;
                       }
                   }

               break ;
               
          }
    }

/************************** end of form_rhs **********************************/ 



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
            )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine form_bc
** form_bc computes the boundary condition terms for stiffness and rhs
** due to imposed conditions
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      int     numel , nen , nee , nel , i , j ,
              active , node , ndof , nsd , nsdp1 ;
      ELEMENT_DATA   *el_pt ;


      nsd = fe_sys.nsd ;
      nsdp1 = nsd + 1 ;
      numel = (grp_ptr->el_info)->numel ;
      nen = (grp_ptr->el_info)->nen ;
      nee = (grp_ptr->el_info)->nee ;
      ndof = fe_sys.ndof ;

      switch((grp_ptr->el_info)->type)
          {
           case TET:
           case BILIN:
           case SERENDIP:
           case TRILIN:
               for( nel = 0 ; nel < numel ; nel++)
                   {
                    el_pt = grp_ptr->el_data + nel ;
                    if( !(el_pt->is_bc) )  continue ;
                    
                    lame_form( code , el_stiff , grp_ptr->el_info ,
                               grp_ptr->el_mat, el_pt ) ;

                    for ( i = 0 ; i < nee ; i++ )
                         el_rhs[i] = ZERO ;
                        
                        active = 0 ;
                        for(j=0 ; j<nen ; j++)
                           {
                            node = *(el_pt->ien + j) - 1 ;
                            for(i=0 ; i<ndof ; i++)
                               {
                                if( *(el_pt->lm + ndof*j + i) != (-1) ||
                                    *(global.forv + node*ndof + i) ==
                                    ZERO )
                                     {
                                      *(temp + ndof*j + i) = ZERO ;
                                     }
                                else
                                     {
                                      *(temp + ndof*j + i) =
                                        *(global.forv + node*ndof + i) ;
                                      active = 1 ;
                                     }
                               }
                           }
                        if(VISCO == code)
                           {
                            for(i=0 ; i<nee ; i++)
                               *(temp + i) *= time_data.dt ;
                           }
                        if( active )
                           adjust_bc( temp , el_rhs , el_stiff , nee ) ;

                    if(VISCO == code)
                       addfor( force.incr_rhs , el_rhs , grp_ptr->el_info ,
                               el_pt ) ;
                    else
                       addfor( force.full_rhs , el_rhs , grp_ptr->el_info ,
                               el_pt ) ;
                   }

               break ;
               
          }
    }

/************************** end of form_bc ***********************************/ 



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
              )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine form_slip
** form_slip computes the split-nodes fault offsets and their influence 
** on the finite element stiffness and RHS.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      int  i , nee , ndof , nen , j , active ;
     
      nen     = info->nen ;
      nee     = info->nee ;
      ndof     = fe_sys.ndof ;
      
      for (  i = 0 ; i < nee ; i++ )
          *(rhs + i) = ZERO ;
      active = 0 ;
          
      lame_form( code , stiff , info , mat_pt , el_pt ) ;

      for(j=0 ; j<nen ; j++)
         {
          if( *(el_pt->ien + j) == node )
             {
              for(i=0 ; i<ndof ; i++)
                  {
                   *(temp + ndof*j + i) = *(slip + i) ;
                   if(cslip)
                     *(temp + ndof*j + i) *= time_data.dt ;
                       /* slip interpreted as a rate, and multiplied by time step */
                    if(QUAKEEVT == code)
                     *(temp + ndof*j + i) *= frac ;
                      /* slip modified by slip fraction in quake adjustment  */
                 }
              active = 1 ;
             }
          else
             {
              for(i=0 ; i<ndof ; i++)
                  {
                   *(temp + ndof*j + i) = ZERO ;
                  }
             }
          }

      if( active )
         adjust_bc( temp , rhs , stiff , nee ) ;

         }

/************************** end of form_slip *********************************/ 



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
              )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine adjust_bc
** adjust_bc computes the necessary terms that modify the right-hand-side
** due to imposed displacements.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

   {
    int i , j ;
    
    for(j=0 ; j<nee ; j++)
       {
        for(i=0 ; i<=j ; i++)
           {
            if(ZERO == disp[i]) continue ;
            
            rhs[j] -= ES(i,j) * disp[i] ;
           }
            
        if( j < (nee-1) )
           {
            for(i=(j+1) ; i<nee ; i++)
               {
                if(ZERO == disp[i]) continue ;
            
                rhs[j] -= ES(j,i) * disp[i] ;
               }
           }
       }
   }
/************************** end of adjust_bc *********************************/ 



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
              )

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

     {
      real  w_det , lam_2mu , lam_el  , mu_el ;
      real  db11 , db21 , db31 , db41 , db51 , db61 ,
             db12 , db22 , db32 , db42 , db52 , db62 ,
             db13 , db23 , db33 , db43 , db53 , db63 ;
      real  bcol[6] , *shape_ptr , *det_ptr , *sh , *dbar , *wt ;
      int    i , j , ipt , it2 , jt2 , it3 , jt3 ,
             nint , nen , nee , nee_sq , ndof , str_size , dbar_size ,
             save_shape , type , nsd , nsdp1 ;
             

      nsd = fe_sys.nsd ;
      nsdp1 = nsd + 1 ;
      nint = el_pt->nint ;
      lam_el  = *(mat_pt->lambda + el_pt->mat) ;
      mu_el   = *(mat_pt->mu + el_pt->mat) ;
      lam_2mu = lam_el + TWO * mu_el ;
      nen     = info->nen ;
      nee     = info->nee ;
      nee_sq = info->el_size ;
      type = info->type ;

      if(TET == type && 1 == el_pt->degen)
         return;

      ndof     = fe_sys.ndof ;
      str_size = info->stress_size ;
      dbar_size = info->dbar_size ;
      save_shape = fe_sys.save_shape;
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
      if( !(*(mat_pt->plastic + el_pt->mat)) )
         code = ELASTIC ;

      if( !save_shape )
         {
          shape_ptr = sh_temp ;
          det_ptr = det_temp ;
          shape( el_pt , info , shape_ptr , det_ptr , type , nen , nint , 1 ) ;
         }
      else
         {
          shape_ptr = el_pt -> shape ;
          det_ptr = el_pt -> deter ;
         }

      for (  i = 0 ; i < nee_sq ; i++ )
          *(es + i) = ZERO ;

      if(VISCO == code)
        {
         form_dbar(el_pt  ,  mat_pt  , info) ;
        }
          
      for (ipt = 0 ; ipt < nint ; ipt++)
           {
            sh = shape_ptr + ipt*nen*nsdp1 ;
            w_det = wt[ipt] * *(det_ptr + ipt) ;
            dbar = (el_pt -> dbar) + ipt*dbar_size ;
            
            switch( ndof )
               {
                case 1:
                       if(ELASTIC == code || QUAKEEVT == code)
                          {
                           for(j=0;j<nen;j++)
                              {
                               for(i=0;i<=j;i++)
                                  {
                                   ES(i,j) += w_det*mu_el*
                                     (SH(0,i)*SH(0,j)+SH(1,i)*SH(1,j)) ;
                                  }
                              }
                          }
                       else if(VISCO == code)
                          {
                           for(j=0;j<nen;j++)
                              {
                               bcol[0]  = w_det * SH(0,j) ;
                               bcol[1]  = w_det * SH(1,j) ;
                               full_back(dbar,bcol,el_diag,str_size) ;
                               db11 = bcol[0] ;
                               db21 = bcol[1] ;
                               for(i=0;i<=j;i++)
                                  {
                                   ES(i,j) +=
                                      SH(0,i)*db11+SH(1,i)*db21 ;
                                  }
                              }
                          }
                   break ;
                   
                case 2:
                       if(ELASTIC == code || QUAKEEVT == code)
                          {
                           for(j=0;j<nen;j++)
                              {
                               jt2 = j+j ;
                               db11 = w_det * lam_2mu * SH(0,j) ;
                               db12 = w_det * lam_el * SH(1,j) ;
                               db21 = w_det * lam_el * SH(0,j) ;
                               db22 = w_det * lam_2mu * SH(1,j) ;
                               db31 = w_det * lam_el * SH(0,j) ;
                               db32 = w_det * lam_el * SH(1,j) ;
                               db41 = w_det * mu_el * SH(1,j) ;
                               db42 = w_det * mu_el * SH(0,j) ;
                               for(i=0;i<=j;i++)
                                  {
                                   it2 = i+i ;
                                   ES(it2 , jt2)     +=
                                      SH(0,i)*db11 + SH(1,i)*db41 ;
                                   ES(it2 , jt2+1)   +=
                                      SH(1,i)*db42 + SH(0,i)*db12 ;
                                   ES(it2+1 , jt2)   +=
                                      SH(0,i)*db41 + SH(1,i)*db21 ;
                                   ES(it2+1 , jt2+1) +=
                                      SH(1,i)*db22 + SH(0,i)*db42 ;
                                  }
                              }
                          }

                       else if(VISCO == code)
                          {
                           for(j=0;j<nen;j++)
                              {
                               jt2 = j+j ;
                               bcol[0]  = w_det * SH(0,j) ;
                               bcol[1]  = ZERO ;
                               bcol[2]  = ZERO ;
                               bcol[3]  = w_det * SH(1,j) ;
                               full_back(dbar,bcol,el_diag,str_size) ;
                               db11 = bcol[0] ;
                               db21 = bcol[1] ;
                               db31 = bcol[2] ;
                               db41 = bcol[3] ;
                               bcol[0]  = ZERO ;
                               bcol[1]  = w_det * SH(1,j) ;
                               bcol[2]  = ZERO ;
                               bcol[3]  = w_det * SH(0,j) ;
                               full_back(dbar,bcol,el_diag,str_size) ;
                               db12 = bcol[0] ;
                               db22 = bcol[1] ;
                               db32 = bcol[2] ;
                               db42 = bcol[3] ;
                               for(i=0;i<=j;i++)
                                  {
                                   it2 = i+i ;
                                   ES(it2 , jt2)     +=
                                      SH(0,i)*db11 + SH(1,i)*db41 ;
                                   ES(it2 , jt2+1)   +=
                                      SH(1,i)*db42 + SH(0,i)*db12 ;
                                   ES(it2+1 , jt2)   +=
                                      SH(0,i)*db41 + SH(1,i)*db21 ;
                                   ES(it2+1 , jt2+1) +=
                                      SH(1,i)*db22 + SH(0,i)*db42 ;
                                  }
                              }
                          }
                   break ;
                   
                case 3:
                       if((ELASTIC == code || QUAKEEVT == code)  && 2 == nsd)
                          {
                           for(j=0;j<nen;j++)
                              {
                               jt3 = 3*j ;
                               db11 = w_det * lam_2mu * SH(0,j) ;
                               db12 = w_det * lam_el * SH(1,j) ;
                               db21 = w_det * lam_el * SH(0,j) ;
                               db22 = w_det * lam_2mu * SH(1,j) ;
                               db41 = w_det * mu_el * SH(1,j) ;
                               db42 = w_det * mu_el * SH(0,j) ;
                               for(i=0;i<=j;i++)
                                  {
                                   it3 = 3*i ;
                                   ES(it3 , jt3)     +=
                                      SH(0,i)*db11 + SH(1,i)*db41 ;
                                   ES(it3+1 , jt3)   +=
                                      SH(1,i)*db21 + SH(0,i)*db41 ;
                                   ES(it3+2 , jt3)    = ZERO ;
                                   ES(it3 , jt3+1)   +=
                                      SH(0,i)*db12 + SH(1,i)*db42 ;
                                   ES(it3+1 , jt3+1) +=
                                      SH(1,i)*db22 + SH(0,i)*db42 ;
                                   ES(it3+2 , jt3+1)  = ZERO ;
                                   ES(it3 , jt3+2)    = ZERO ;
                                   ES(it3+1 , jt3+2)  = ZERO ;
                                   ES(it3+2 , jt3+2) +=
                                      SH(0,i)*db42 + SH(1,i)*db41 ;
                                  }
                              }
                          }

                       else if(VISCO == code && 2 == nsd)
                          {
                           for(j=0;j<nen;j++)
                              {
                               jt3 = 3*j ;
                               bcol[0]  = w_det * SH(0,j) ;
                               bcol[1]  = ZERO ;
                               bcol[2]  = ZERO ;
                               bcol[3]  = w_det * SH(1,j) ;
                               bcol[4]  = ZERO ;
                               bcol[5]  = ZERO ;
                               full_back(dbar,bcol,el_diag,str_size) ;
                               db11 = bcol[0] ;
                               db21 = bcol[1] ;
                               db31 = bcol[2] ;
                               db41 = bcol[3] ;
                               db51 = bcol[4] ;
                               db61 = bcol[5] ;
                               bcol[0]  = ZERO ;
                               bcol[1]  = w_det * SH(1,j) ;
                               bcol[2]  = ZERO ;
                               bcol[3]  = w_det * SH(0,j) ;
                               bcol[4]  = ZERO ;
                               bcol[5]  = ZERO ;
                               full_back(dbar,bcol,el_diag,str_size) ;
                               db12 = bcol[0] ;
                               db22 = bcol[1] ;
                               db32 = bcol[2] ;
                               db42 = bcol[3] ;
                               db52 = bcol[4] ;
                               db62 = bcol[5] ;
                               bcol[0]  = ZERO ;
                               bcol[1]  = ZERO ;
                               bcol[2]  = ZERO ;
                               bcol[3]  = ZERO ;
                               bcol[4]  = w_det * SH(0,j) ;
                               bcol[5]  = w_det * SH(1,j) ;
                               full_back(dbar,bcol,el_diag,str_size) ;
                               db13 = bcol[0] ;
                               db23 = bcol[1] ;
                               db33 = bcol[2] ;
                               db43 = bcol[3] ;
                               db53 = bcol[4] ;
                               db63 = bcol[5] ;
                               for(i=0;i<=j;i++)
                                  {
                                   it3 = 3*i ;
                                   ES(it3 , jt3)     +=
                                      SH(0,i)*db11 + SH(1,i)*db41 ;
                                   ES(it3+1 , jt3)   +=
                                      SH(1,i)*db21 + SH(0,i)*db41 ;
                                   ES(it3+2 , jt3)    =
                                      SH(0,i)*db51 + SH(1,i)*db61 ;
                                   ES(it3 , jt3+1)   +=
                                      SH(0,i)*db12 + SH(1,i)*db42 ;
                                   ES(it3+1 , jt3+1) +=
                                      SH(1,i)*db22 + SH(0,i)*db42 ;
                                   ES(it3+2 , jt3+1)  =
                                      SH(0,i)*db52 + SH(1,i)*db62 ;
                                   ES(it3 , jt3+2)    =
                                      SH(0,i)*db13 + SH(1,i)*db43 ;
                                   ES(it3+1 , jt3+2)  =
                                      SH(1,i)*db23 + SH(0,i)*db43 ;
                                   ES(it3+2 , jt3+2) +=
                                      SH(0,i)*db53 + SH(1,i)*db63 ;
                                  }
                              }
                          }
                    if((ELASTIC == code || QUAKEEVT == code) && 3 == nsd)
                          {
                           for(j=0;j<nen;j++)
                              {
                               jt3 = 3*j ;
                               db11 = w_det * lam_2mu * SH(0,j) ;
                               db12 = w_det * lam_el * SH(1,j) ;
                               db13 = w_det * lam_el * SH(2,j) ;
                               db21 = w_det * lam_el * SH(0,j) ;
                               db22 = w_det * lam_2mu * SH(1,j) ;
                               db23 = w_det * lam_el * SH(2,j) ;
                               db31 = w_det * lam_el * SH(0,j) ;
                               db32 = w_det * lam_el * SH(1,j) ;
                               db33 = w_det * lam_2mu * SH(2,j) ;
                               db41 = w_det * mu_el * SH(1,j) ;
                               db42 = w_det * mu_el * SH(0,j) ;
                               db51 = w_det * mu_el * SH(2,j) ;
                               db53 = w_det * mu_el * SH(0,j) ;
                               db62 = w_det * mu_el * SH(2,j) ;
                               db63 = w_det * mu_el * SH(1,j) ;
                               for(i=0;i<=j;i++)
                                  {
                                   it3 = 3*i ;
                                   ES(it3 , jt3)     +=
                                    SH(0,i)*db11 + SH(1,i)*db41 + SH(2,i)*db51 ;
                                   ES(it3+1 , jt3)   +=
                                      SH(1,i)*db21 + SH(0,i)*db41 ;
                                   ES(it3+2 , jt3)    +=
                                      SH(2,i)*db31 + SH(0,i)*db51 ;
                                   ES(it3 , jt3+1)   +=
                                      SH(0,i)*db12 + SH(1,i)*db42 ;
                                   ES(it3+1 , jt3+1) +=
                                    SH(1,i)*db22 + SH(0,i)*db42 + SH(2,i)*db62 ;
                                   ES(it3+2 , jt3+1)  +=  
                                      SH(2,i)*db32 + SH(1,i)*db62 ;
                                   ES(it3 , jt3+2)    += 
                                      SH(0,i)*db13 + SH(2,i)*db53 ;
                                   ES(it3+1 , jt3+2)  +=
                                      SH(1,i)*db23 + SH(2,i)*db63 ;
                                   ES(it3+2 , jt3+2) +=
                                    SH(2,i)*db33 + SH(0,i)*db53 + SH(1,i)*db63;
                                  }
                              }
                          }

                       else if(VISCO == code && 3 == nsd)
                          {
                           for(j=0;j<nen;j++)
                              {
                               jt3 = 3*j ;
                               bcol[0]  = w_det * SH(0,j) ;
                               bcol[1]  = ZERO ;
                               bcol[2]  = ZERO ;
                               bcol[3]  = w_det * SH(1,j) ;
                               bcol[4]  = w_det * SH(2,j) ;
                               bcol[5]  = ZERO ;
                               full_back(dbar,bcol,el_diag,str_size) ;
                               db11 = bcol[0] ;
                               db21 = bcol[1] ;
                               db31 = bcol[2] ;
                               db41 = bcol[3] ;
                               db51 = bcol[4] ;
                               db61 = bcol[5] ;
                               bcol[0]  = ZERO ;
                               bcol[1]  = w_det * SH(1,j) ;
                               bcol[2]  = ZERO ;
                               bcol[3]  = w_det * SH(0,j) ;
                               bcol[4]  = ZERO ;
                               bcol[5]  = w_det * SH(2,j) ;
                               full_back(dbar,bcol,el_diag,str_size) ;
                               db12 = bcol[0] ;
                               db22 = bcol[1] ;
                               db32 = bcol[2] ;
                               db42 = bcol[3] ;
                               db52 = bcol[4] ;
                               db62 = bcol[5] ;
                               bcol[0]  = ZERO ;
                               bcol[1]  = ZERO ;
                               bcol[2]  = w_det * SH(2,j) ;
                               bcol[3]  = ZERO ;
                               bcol[4]  = w_det * SH(0,j) ;
                               bcol[5]  = w_det * SH(1,j) ;
                               full_back(dbar,bcol,el_diag,str_size) ;
                               db13 = bcol[0] ;
                               db23 = bcol[1] ;
                               db33 = bcol[2] ;
                               db43 = bcol[3] ;
                               db53 = bcol[4] ;
                               db63 = bcol[5] ;
                               for(i=0;i<=j;i++)
                                  {
                                   it3 = 3*i ;
                                   ES(it3 , jt3)     +=
                                    SH(0,i)*db11 + SH(1,i)*db41 + SH(2,i)*db51;
                                   ES(it3+1 , jt3)   +=
                                    SH(1,i)*db21 + SH(0,i)*db41 + SH(2,i)*db61;
                                   ES(it3+2 , jt3)    =
                                    SH(0,i)*db51 + SH(1,i)*db61 + SH(2,i)*db31;
                                   ES(it3 , jt3+1)   +=
                                    SH(0,i)*db12 + SH(1,i)*db42 + SH(2,i)*db52;
                                   ES(it3+1 , jt3+1) +=
                                    SH(1,i)*db22 + SH(0,i)*db42 + SH(2,i)*db62;
                                   ES(it3+2 , jt3+1)  =
                                    SH(0,i)*db52 + SH(1,i)*db62 + SH(2,i)*db32 ;
                                   ES(it3 , jt3+2)    =
                                    SH(0,i)*db13 + SH(1,i)*db43 + SH(2,i)*db53;
                                   ES(it3+1 , jt3+2)  =
                                    SH(1,i)*db23 + SH(0,i)*db43 + SH(2,i)*db63;
                                   ES(it3+2 , jt3+2) +=
                                    SH(0,i)*db53 + SH(1,i)*db63 + SH(2,i)*db33;
                                  }
                              }
                          }
                   break ;
                   
               }
                 
           }
  }  

/************************** end of lame_form *********************************/ 



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
              )

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
{
      real  bmat[6] , *cvect , len , mod ;
      int  nsd , nsdp1 , nen , nee , nee_sq , type , ndof , option ;
      int  i , j ;

      nsd = fe_sys.nsd ;
      nsdp1 = nsd + 1 ;
      len = *(mat_pt->lambda + el_pt->mat) ;
      mod   = *(mat_pt->mu + el_pt->mat) / (len*len) ;
      nen     = info->nen ;
      nee     = info->nee ;
      nee_sq = info->el_size ;
      type = info->type ;
      ndof     = fe_sys.ndof ;
      cvect = mat_pt->bforce + ndof*(el_pt->mat) ; /* grav is used as holding space for cvect */
      option = (int) ( *(mat_pt->visc + 2*(el_pt->mat)) ) ;
         /* visc == zero for truss ; visc non-zero for spring */

      for (  i = 0 ; i < nee_sq ; i++ )
          *(es + i) = ZERO ;
          
            switch( ndof )
               {
                case 1:
                       if(option == 0)  /* truss */
                          {
                           /* no such thing in 1-D */
                          }
                          
                       else  /* isospring */
                          {
                           ES(0 , 0) =   mod ;
                           ES(0 , 1) =  -mod ;
                           ES(1 , 0) =  -mod ;
                           ES(1 , 1) =   mod ;
                          }
                   break ;
                   
                case 2:
                       if(option == 0)  /* truss */
                          {
                           for(i=0;i<2;i++)
                              {
                               *(bmat+i) = -(*(cvect+i)) ;
                               *(bmat+i+2) = *(cvect+i) ;
                              }
                           for(j=0;j<nee;j++)
                              {
                               for(i=0;i<=j;i++)
                                  {
                                   ES(i , j) = *(bmat+i) * mod * *(bmat+j) ;
                                   ES(j , i) = ES(i , j) ;
                                  }
                               }
                          }
                          
                       else  /* isospring */
                          {
                           for(j=0;j<nee;j++)
                              {
                               for(i=0;i<=j;i++)
                                  {
                                   if(i == j)
                                     ES(i , j) =  mod ;
                                   else if(i == j-2)
                                     ES(i , j) =  -mod ;
                                   else
                                     ES(i , j) =  ZERO ;
                                   ES(j , i) = ES(i , j) ;
                                  }
                               }
                          }
                   break ;
                   
                case 3:
                       if(option == 0)  /* truss */
                          {
                           for(i=0;i<3;i++)
                              {
                               *(bmat+i) = -(*(cvect+i)) ;
                               *(bmat+i+3) = *(cvect+i) ;
                              }
                           for(j=0;j<nee;j++)
                              {
                               if( j > 5) continue ;
                                /* nodes 3 and 4 are unconnected and dormant */
                                
                               for(i=0;i<=j;i++)
                                  {
                                   ES(i , j) = *(bmat+i) * mod * *(bmat+j) ;
                                   if(i != j)
                                      ES(j , i) = ES(i , j) ;
                                  }
                               }
                              
                          }
                       else   /* isospring */
                          {
                           for(j=0;j<nee;j++)
                              {
                               if( j > 5) continue ;
                                /* nodes 3 and 4 are unconnected and dormant */
                                
                               for(i=0;i<=j;i++)
                                  {
                                   if(i == j)
                                     ES(i , j) =  mod ;
                                   else if(i == j-3)
                                     ES(i , j) =  -mod ;
                                   else
                                     ES(i , j) =  ZERO ;
                                   if(i != j)
                                      ES(j , i) = ES(i , j) ;
                                  }
                               }
                          }
                   break ;
                   
               }
}


/************************** end of duo_form *********************************/ 



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
              )

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

     {
      real    *grav , *shape_ptr , *det_ptr , *sh , *dbar , *wt ;
      real    w_det , btemp[6] ;
      int      i , j , ipt , active , jt2 , jt3 , nint , nen , nee ,
               ndof , str_size , dbar_size , save_shape , type , nsd , nsdp1 ;
               
      nsd = fe_sys.nsd ;
      nsdp1 = nsd + 1 ;
      ndof  = fe_sys.ndof ;
      nint = el_pt->nint ;
      grav  = mat_pt->bforce + ndof*(el_pt->mat) ;
      nen     = info->nen ;
      nee     = info->nee ;
      type = info->type ;

      if(TET == type && 1 == el_pt->degen)
         return;

      str_size = info->stress_size ;
      dbar_size = info->dbar_size ;
      save_shape = fe_sys.save_shape;
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
      

      if( !(*(mat_pt->plastic + el_pt->mat)) )
         code = ELASTIC ;

      for (  i = 0 ; i < nee ; i++ )
          *(rhs + i) = ZERO ;
          

      if(VISCO == code)
         form_beta(el_pt  ,  mat_pt  , info) ;
          
      if( (ELASTIC == code) && time_data.elastic)
         {
          active = 0 ;
          for(i=0 ; i<ndof ; i++)
            if( *(grav + i) != ZERO )  active = 1 ;
          if( !active ) return ;
         }
        
      if( !save_shape )
         {
          shape_ptr = sh_temp ;
          det_ptr = det_temp ;
          shape( el_pt , info , shape_ptr , det_ptr , type , nen , nint , 1 ) ;
         }
      else
         {
          shape_ptr = el_pt -> shape ;
          det_ptr = el_pt -> deter ;
         }

      for (ipt = 0 ; ipt < nint ; ipt++)
           {
            sh = shape_ptr + ipt*nen*nsdp1 ;
            w_det = wt[ipt] * *(det_ptr + ipt) ;
            dbar = (el_pt -> dbar) + ipt*dbar_size ;
            
            switch( ndof )
               {
                case 1:
                   if(ELASTIC == code)
                      {
                       for(j=0;j<nen;j++)
                          {
                           *(rhs + j) += w_det * (*grav) * SH(nsd,j) ;
                          }
                      }
 
                    else if(VISCO == code)
                      {
                       move_real( (el_pt->bta)+ipt*str_size ,
                                   btemp , str_size ) ;
                       full_back(dbar,btemp,el_diag,str_size) ;

                       for(j=0;j<nen;j++)
                          {
                           *(rhs + j) += w_det *
                             (btemp[0] * SH(0,j) + btemp[1] * SH(1,j)) ;
                          }
                      }
                  break ;
                   
                case 2:
                   if(ELASTIC == code)
                      {
                       for(j=0;j<nen;j++)
                          {
                           jt2 = j+j ;
                           *(rhs + jt2)     +=
                              w_det * (*grav) * SH(nsd,j) ;
                           *(rhs + jt2 + 1) +=
                              w_det * (*(grav + 1)) * SH(nsd,j) ;
                          }
                      }
 
                    else if(VISCO == code)
                      {
                       move_real( (el_pt->bta)+ipt*str_size ,
                                   btemp , str_size ) ;
                       full_back(dbar,btemp,el_diag,str_size) ;

                       for(j=0;j<nen;j++)
                          {
                           jt2 = j+j ;

                           *(rhs + jt2)     += w_det *
                             (btemp[0] * SH(0,j) + btemp[3] * SH(1,j)) ;
                           *(rhs + jt2 + 1) += w_det *
                             (btemp[1] * SH(1,j) + btemp[3] * SH(0,j)) ;
                          }
                      }
                   break ;
                   
                case 3:
                   if(ELASTIC == code)
                      {
                       for(j=0;j<nen;j++)
                          {
                           jt3 = 3*j ;
                           *(rhs + jt3)     +=
                              w_det * (*grav) * SH(nsd,j) ;
                           *(rhs + jt3 + 1) +=
                              w_det * (*(grav + 1)) * SH(nsd,j) ;
                           *(rhs + jt3 + 2) +=
                              w_det * (*(grav + 2)) * SH(nsd,j) ;
                          }
                      }
 
                    else if(VISCO == code)
                      {
                       move_real( (el_pt->bta)+ipt*str_size ,
                                   btemp , str_size ) ;
                       full_back(dbar,btemp,el_diag,str_size) ;

                       if(2 == nsd)
                       {
                       for(j=0;j<nen;j++)
                          {
                           jt3 = 3*j ;

                           *(rhs + jt3)     += w_det * 
                             (btemp[0] * SH(0,j) + btemp[3] * SH(1,j)) ;
                           *(rhs + jt3 + 1) += w_det *
                             (btemp[1] * SH(1,j) + btemp[3] * SH(0,j)) ;
                           *(rhs + jt3 + 2) += w_det *
                             (btemp[4] * SH(0,j) + btemp[5] * SH(1,j)) ;
                          }
                       } else if (3 == nsd)
                       {
                       for(j=0;j<nen;j++)
                          {
                           jt3 = 3*j ;

                           *(rhs + jt3)     += w_det * (btemp[4] * SH(2,j) +
                             btemp[0] * SH(0,j) + btemp[3] * SH(1,j)) ;
                           *(rhs + jt3 + 1) += w_det * (btemp[5] * SH(2,j) +
                             btemp[1] * SH(1,j) + btemp[3] * SH(0,j)) ;
                           *(rhs + jt3 + 2) += w_det * (btemp[2] * SH(2,j) +
                             btemp[4] * SH(0,j) + btemp[5] * SH(1,j)) ;
                          }
                       }
                      }
                   break ;
                   
               }
                 
           }
  }  
/************************** end of force_form ********************************/ 



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
              )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine surf_form
** surf_form computes the surface traction forcing term for the 
** right-hand-side in the finite element problem.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      int  i , nee , ndof , n1 , n2 , n3 , n4 , node1 , node2 , node3 , node4 , nsd , nsdp1 ;
      int  ipt , n , j ;
      real  xx1 , yy1 , zz1 , xx2 , yy2 , zz2 , xx3 , yy3 , zz3 , xx4 , yy4 , zz4 ,
            alength , dlength , elength , adotb , sum , xpt , ypt , hold , hold1 , w_det ,
            dx , dy , dz , dl , ax , ay , az , bx , by , bz , cx , cy , cz , area ,
            ex , ey , ez , fx , fy , fz , cdote , rpt , spt ;
      real  xplan[4][2]  ; /* 4 nodes by 2 components  */
      real  sh2all[48] , detall[4] ; /* 4 integration pts, 4 nodes , 3 SH's */
      real  xs[2][2] ;
      real  *sh2 , *det , *wt ;
     
      nsd   = fe_sys.nsd ;
      nsdp1 = nsd + 1 ;
      nee   = info->nee ;
      ndof  = fe_sys.ndof ;
      
      for ( i = 0 ; i < nee ; i++ )
          *(rhs + i) = ZERO ;
          
      
      if(2 == nsd)
         {
/* ***** Applies surface tractions to an element edge ***** */
/* ***** This version applies two-point surface tractions to a
         bilinear element edge.  DO NOT USE SERENDIP OPTION ***** */
         node1 = side - 1 ;
         node2 = side ;
         if( node2 == 4 )  node2 = 0 ;
         n1 = *(el_pt->ien + node1) - 1 ;
         n2 = *(el_pt->ien + node2) - 1 ;
         xx1 = *(global.coords + n1*nsd) ;
         yy1 = *(global.coords + n1*nsd + 1) ;
         xx2 = *(global.coords + n2*nsd) ;
         yy2 = *(global.coords + n2*nsd + 1) ;
         dx = xx2 - xx1 ;
         dy = yy2 - yy1 ;
         dl = sqrt(dx*dx + dy*dy) ;

/* changed to global traction reference, rather than normal-tangential
   formulation, and changed to single midpoint - implemented 10/19/01  */
         switch( ndof )
            {
             case 1:
                *(rhs + node1) =
                  dl * PT5 * (*trac) ;

                *(rhs + node2) =
                  dl * PT5 * (*trac) ;
                break ;
                   
             case 2:
                *(rhs + 2*node1    ) =
                  dl * PT5 * (*trac) ;
                *(rhs + 2*node1 + 1) =
                  dl * PT5 * (*(trac+1)) ;

                *(rhs + 2*node2    ) =
                  dl * PT5 * (*trac) ;
                *(rhs + 2*node2 + 1) =
                  dl * PT5 * (*(trac+1)) ;
                break ;
                   
             case 3:
                *(rhs + 3*node1    ) =
                  dl * PT5 * (*trac) ;
                *(rhs + 3*node1 + 1) =
                  dl * PT5 * (*(trac+1)) ;
                *(rhs + 3*node1 + 2) =
                  dl * PT5 * (*(trac+2)) ;

                *(rhs + 3*node2    ) =
                  dl * PT5 * (*trac) ;
                *(rhs + 3*node2 + 1) =
                  dl * PT5 * (*(trac+1)) ;
                *(rhs + 3*node2 + 2) =
                  dl * PT5 * (*(trac+2)) ;
                break ;
         }
      } 
      else if (3 == nsd && (TET == info->type))
         {
/************************************************************************/
/*
This version implements one-point surface tractions at the face
centroid of tetrahedral elements.  
GL 10/23/98

8/21/09 - adding modifications to do full 2D Gaussian integration
of surface tractions
*/
/************************************************************************/

  /* check and correct for tets that have been flopped in chirality */
         if(el_pt->degen == 2)
            {
             if(side == 3)  side = 4 ;
             else if(side == 4)  side = 3 ;
            }

         switch(side)
           {
            case 1:
              node1 = 2 ;
              node2 = 3 ;
              node3 = 4 ;
              break;
            case 2:
              node1 = 1 ;
              node2 = 3 ;
              node3 = 4 ;
              break;
            case 3:
              node1 = 1 ;
              node2 = 2 ;
              node3 = 4 ;
              break;
            case 4:
              node1 = 1 ;
              node2 = 2 ;
              node3 = 3 ;
              break;
         }
         n1 = *(el_pt->ien + node1 - 1) - 1 ;  /* look up global node #'s */
         n2 = *(el_pt->ien + node2 - 1) - 1 ;
         n3 = *(el_pt->ien + node3 - 1) - 1 ;
         xx1 = *(global.coords + n1*nsd) ;     /* then get the global coords */
         yy1 = *(global.coords + n1*nsd + 1) ;
         zz1 = *(global.coords + n1*nsd + 2) ;
         xx2 = *(global.coords + n2*nsd) ;
         yy2 = *(global.coords + n2*nsd + 1) ;
         zz2 = *(global.coords + n2*nsd + 2) ;
         xx3 = *(global.coords + n3*nsd) ;
         yy3 = *(global.coords + n3*nsd + 1) ;
         zz3 = *(global.coords + n3*nsd + 2) ;

         ax = xx2 - xx1 ;    /* compute the edge vectors... */
         ay = yy2 - yy1 ;
         az = zz2 - zz1 ;
         bx = xx3 - xx1 ;
         by = yy3 - yy1 ;
         bz = zz3 - zz1 ;
         
         /* now create local planar coordinates for shape function use */
         
         xplan[0][0] = ZERO ;  /* node 1 */
         xplan[0][1] = ZERO ;
         
         xplan[1][0] = alength = sqrt(ax*ax+ay*ay+az*az) ;  /* node 2 */
         xplan[1][1] = ZERO ;
         
         adotb = ax*bx+ay*by+az*bz ;
         cx = bx - ax*adotb/(alength*alength) ;
         cy = by - ay*adotb/(alength*alength) ;
         cz = bz - az*adotb/(alength*alength) ;
         
         xplan[2][0] = adotb/alength ;  /* node 3 */
         xplan[2][1] = sqrt(cx*cx+cy*cy+cz*cz) ;
                  


	for(ipt=0 ; ipt<3 ; ipt++)
	   {
	    sh2 = sh2all + ipt*3*3 ;
	    det = detall +ipt ;
	    
		   rpt = r_3edge[ipt] ;
		   spt = s_3edge[ipt] ;
		   wt = wt_3edge ;
	    
/* next generate the 2D shape functions */

	   *(sh2 + 0*3) = ONE ;
	   *(sh2 + 0*3 + 1) = ZERO ; 
	   *(sh2 + 0*3 + 2) = rpt ;

	   *(sh2 + 1*3) = ZERO ;
	   *(sh2 + 1*3 + 1) = ONE ; 
	   *(sh2 + 1*3 + 2) = spt ;

	   *(sh2 + 2*3) = -ONE ;
	   *(sh2 + 2*3 + 1) = -ONE ; 
	   *(sh2 + 2*3 + 2) = ONE - rpt - spt ;

          for(i=0 ; i<2 ; i++)
             {
              for(j=0 ; j<2 ; j++)
                 {
              
/* note that the first argument (i) counts the physical space coordinate x, 
while the second (j) counts the parent space coordinate Tse. 
Hence xs[j][i] = dX(j) / dTse(i).
*/

                sum = ZERO ;
                for(n=0 ; n<3 ; n++)
					{
					 sum += xplan[n][j] * *(sh2+3*n+i) ;
					}
                 xs[j][i] = sum ;
                 }
             }

           *det = xs[0][0]*xs[1][1] - xs[0][1]*xs[1][0] ;
           				  
           if(*det < ZERO)
              {
               printf("Negative face det: element=%d  ipt=%d\n",el_pt->nel,ipt) ;
               /* this should never happen ; just check in case */
              }
           
          for(i=0 ; i<3 ; i++)
             {
              hold = (xs[1][1] * *(sh2 + i*3) -
                      xs[1][0] * *(sh2 + i*3 + 1)) ;

              hold1 = (xs[0][0] * *(sh2 + i*3 + 1) -
                       xs[0][1] * *(sh2 + i*3)) ;

              *(sh2 + i*3)     = hold/ *det ;
              *(sh2 + i*3 + 1) = hold1/ *det ;
             }

        w_det = wt[ipt] * *det ;
        

		*(rhs + 3*(node1-1)    ) +=  w_det * (*trac) * *(sh2+(3*0)+2) ;
		*(rhs + 3*(node1-1) + 1) +=  w_det * (*(trac+1)) * *(sh2+(3*0)+2) ;
		*(rhs + 3*(node1-1) + 2) +=  w_det * (*(trac+2)) * *(sh2+(3*0)+2) ;

		*(rhs + 3*(node2-1)    ) +=  w_det * (*trac) * *(sh2+(3*1)+2) ;
		*(rhs + 3*(node2-1) + 1) +=  w_det * (*(trac+1)) * *(sh2+(3*1)+2) ;
		*(rhs + 3*(node2-1) + 2) +=  w_det * (*(trac+2)) * *(sh2+(3*1)+2) ;

		*(rhs + 3*(node3-1)    ) +=  w_det * (*trac) * *(sh2+(3*2)+2) ;
		*(rhs + 3*(node3-1) + 1) +=  w_det * (*(trac+1)) * *(sh2+(3*2)+2) ;
		*(rhs + 3*(node3-1) + 2) +=  w_det * (*(trac+2)) * *(sh2+(3*2)+2) ;


	   }

           }

      else if (3 == nsd && (TRILIN == info->type))
         {
/************************************************************************/
/*
This version implements flattened quad integrated tractions at the faces
 of brick elements.  
GL 10/02/09
*/
/************************************************************************/
         switch(side)
           {
            case 1:
              node1 = 1 ;
              node2 = 4 ;
              node3 = 3 ;
              node4 = 2 ;
              break;
            case 2:
              node1 = 5 ;
              node2 = 6 ;
              node3 = 7 ;
              node4 = 8 ;
              break;
            case 3:
              node1 = 1 ;
              node2 = 2 ;
              node3 = 6 ;
              node4 = 5 ;
              break;
            case 4:
              node1 = 2 ;
              node2 = 3 ;
              node3 = 7 ;
              node4 = 6 ;
              break;
            case 5:
              node1 = 3 ;
              node2 = 4 ;
              node3 = 8 ;
              node4 = 7 ;
              break;
            case 6:
              node1 = 1 ;
              node2 = 5 ;
              node3 = 8 ;
              node4 = 4 ;
              break;
         }
         
         n1 = *(el_pt->ien + node1 - 1) - 1 ;  /* look up global node #'s */
         n2 = *(el_pt->ien + node2 - 1) - 1 ;
         n3 = *(el_pt->ien + node3 - 1) - 1 ;
         n4 = *(el_pt->ien + node4 - 1) - 1 ;
         xx1 = *(global.coords + n1*nsd) ;     /* then get the global coords */
         yy1 = *(global.coords + n1*nsd + 1) ;
         zz1 = *(global.coords + n1*nsd + 2) ;
         xx2 = *(global.coords + n2*nsd) ;
         yy2 = *(global.coords + n2*nsd + 1) ;
         zz2 = *(global.coords + n2*nsd + 2) ;
         xx3 = *(global.coords + n3*nsd) ;
         yy3 = *(global.coords + n3*nsd + 1) ;
         zz3 = *(global.coords + n3*nsd + 2) ;
         xx4 = *(global.coords + n4*nsd) ;
         yy4 = *(global.coords + n4*nsd + 1) ;
         zz4 = *(global.coords + n4*nsd + 2) ;

         ax = xx2 - xx1 ;    /* compute the edge vectors... */
         ay = yy2 - yy1 ;
         az = zz2 - zz1 ;
         bx = xx4 - xx1 ;
         by = yy4 - yy1 ;
         bz = zz4 - zz1 ;
         cx = xx3 - xx1 ;
         cy = yy3 - yy1 ;
         cz = zz3 - zz1 ;
         
         /* now create local planar coordinates for shape function use */
         
         xplan[0][0] = ZERO ;  /* node 1 */
         xplan[0][1] = ZERO ;
         
         xplan[1][0] = alength = sqrt(ax*ax+ay*ay+az*az) ;  /* node 2 */
         xplan[1][1] = ZERO ;
         
         adotb = ax*bx+ay*by+az*bz ;
         
         dx = bx - ax*adotb/(alength*alength) ; /* d is a vector perp to a with b's projected length */
         dy = by - ay*adotb/(alength*alength) ;
         dz = bz - az*adotb/(alength*alength) ;
         
         xplan[3][0] = adotb/alength ;  /* node 4 */
         xplan[3][1] = dlength = sqrt(dx*dx+dy*dy+dz*dz) ;
         
         ex = ay*bz - az*by ;
         ey = az*bx - ax*bz ;
         ez = ax*by - ay*bx ;
         elength = sqrt(ex*ex+ey*ey+ez*ez) ;
         ex /= elength ;
         ey /= elength ;
         ez /= elength ;   /* e is a unit vector perp to the axb plane */
         
         cdote = cx*ex + cy*ey + cz*ez ;
         
         fx = cx - ex*cdote ; /* f is the c vector flattened into the axb plane */
         fy = cy - ey*cdote ;
         fz = cz - ez*cdote ;
         
         xplan[2][0] = (fx*ax+fy*ay+fz*az)/alength ;  /* node 3 */
         xplan[2][1] = (fx*dx+fy*dy+fz*dz)/dlength ;


	for(ipt=0 ; ipt<4 ; ipt++)
	   {
	    sh2 = sh2all + ipt*4*3 ;
	    det = detall +ipt ;
	    
		   xpt = x_2x2[ipt] ;
		   ypt = y_2x2[ipt] ;
		   xpt /= (real) sqrt(THREE) ;
		   ypt /= (real) sqrt(THREE) ;
		   wt = wt_2x2 ;
	    
/* next generate the 2D shape functions */
               *(sh2 + 0*3) = PT25 * (ypt - ONE) ;
               *(sh2 + 0*3 + 1) = PT25 * (xpt - ONE) ; 
               *(sh2 + 0*3 + 2) = PT25 * (ONE - xpt) * (ONE - ypt) ;

               *(sh2 + 1*3) = PT25 * (ONE - ypt) ;
               *(sh2 + 1*3 + 1) = (-PT25) * (ONE + xpt) ; 
               *(sh2 + 1*3 + 2) = PT25 * (ONE + xpt) * (ONE - ypt) ;

               *(sh2 + 2*3) = PT25 * (ONE + ypt) ;
               *(sh2 + 2*3 + 1) = PT25 * (ONE + xpt) ; 
               *(sh2 + 2*3 + 2) = PT25 * (ONE + xpt) * (ONE + ypt) ;

               *(sh2 + 3*3) = (-PT25) * (ONE + ypt) ;
               *(sh2 + 3*3 + 1) = PT25 * (ONE - xpt) ; 
               *(sh2 + 3*3 + 2) = PT25 * (ONE - xpt) * (ONE + ypt) ;

          for(i=0 ; i<2 ; i++)
             {
              for(j=0 ; j<2 ; j++)
                 {
              
/* note that the first argument (i) counts the physical space coordinate x, 
while the second (j) counts the parent space coordinate Tse. 
Hence xs[j][i] = dX(j) / dTse(i).
*/

                sum = ZERO ;
                for(n=0 ; n<4 ; n++)
					{
					 sum += xplan[n][j] * *(sh2+3*n+i) ;
					}
                 xs[j][i] = sum ;
                 }
             }

           *det = xs[0][0]*xs[1][1] - xs[0][1]*xs[1][0] ;
           
           if(*det < ZERO)
              {
               printf("Negative face det: element=%d  ipt=%d\n",el_pt->nel,ipt) ;
               /* this should never happen ; just check in case */
              }

          for(i=0 ; i<4 ; i++)
             {
              hold = (xs[1][1] * *(sh2 + i*3) -
                      xs[1][0] * *(sh2 + i*3 + 1)) ;

              hold1 = (xs[0][0] * *(sh2 + i*3 + 1) -
                       xs[0][1] * *(sh2 + i*3)) ;

              *(sh2 + i*3)     = hold/ *det ;
              *(sh2 + i*3 + 1) = hold1/ *det ;
             }

        w_det = wt[ipt] * *det ;
        

		*(rhs + 3*(node1-1)    ) +=  w_det * (*trac) * *(sh2+(3*0)+2) ;
		*(rhs + 3*(node1-1) + 1) +=  w_det * (*(trac+1)) * *(sh2+(3*0)+2) ;
		*(rhs + 3*(node1-1) + 2) +=  w_det * (*(trac+2)) * *(sh2+(3*0)+2) ;

		*(rhs + 3*(node2-1)    ) +=  w_det * (*trac) * *(sh2+(3*1)+2) ;
		*(rhs + 3*(node2-1) + 1) +=  w_det * (*(trac+1)) * *(sh2+(3*1)+2) ;
		*(rhs + 3*(node2-1) + 2) +=  w_det * (*(trac+2)) * *(sh2+(3*1)+2) ;

		*(rhs + 3*(node3-1)    ) +=  w_det * (*trac) * *(sh2+(3*2)+2) ;
		*(rhs + 3*(node3-1) + 1) +=  w_det * (*(trac+1)) * *(sh2+(3*2)+2) ;
		*(rhs + 3*(node3-1) + 2) +=  w_det * (*(trac+2)) * *(sh2+(3*2)+2) ;

		*(rhs + 3*(node4-1)    ) +=  w_det * (*trac) * *(sh2+(3*3)+2) ;
		*(rhs + 3*(node4-1) + 1) +=  w_det * (*(trac+1)) * *(sh2+(3*3)+2) ;
		*(rhs + 3*(node4-1) + 2) +=  w_det * (*(trac+2)) * *(sh2+(3*3)+2) ;

	   }


           }

  }  
/************************** end of surf_form *********************************/ 



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
              )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine buoy_form
** buoy_form computes the buoyancy traction term due to displacement
** of a density contrast in the finite element grid.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      int  i , nee , ndof , n1 , n2 , n3 , n4 , node1 , node2 , node3 , node4 , nsd , nsdp1 ;
      int  ipt , n , j ;
      real  xx1 , yy1 , zz1 , xx2 , yy2 , zz2 , xx3 , yy3 , zz3 ,  xx4 , yy4 , zz4 ,
            alength , dlength , elength , adotb , sum , xpt , ypt , hold , hold1 , w_det ,
            dx , dy , dz , dl , ax , ay , az , bx , by , bz , cx , cy , cz , area ,
            ux1 , uy1 , uz1 , ux2 , uy2 , uz2 , ux3 , uy3 , uz3 , avex , avey , avez ,
            ex , ey , ez , fx , fy , fz , cdote , ux4 , uy4 , uz4 , rpt , spt ,
            reactx , reacty , reactz , centrx , centry , centrz , local_up[3] , dist ;
      real  xplan[4][2]  ; /* 4 nodes by 2 components  */
      real  sh2all[48] , detall[4] ; /* 4 integration pts, 4 nodes , 3 SH's */
      real  xs[2][2] ;
      real  *sh2 , *det , *wt ;
     
      nsd   = fe_sys.nsd ;
      nsdp1 = nsd + 1 ;
      nee   = info->nee ;
      ndof  = fe_sys.ndof ;
      
      for ( i = 0 ; i < nee ; i++ )
          *(rhs + i) = ZERO ;
          
          
      if (3 == nsd && (TET == info->type))
         {

  /* check and correct for tets that have been flopped in chirality */
         if(el_pt->degen == 2)
            {
             if(side == 3)  side = 4 ;
             else if(side == 4)  side = 3 ;
            }

         switch(side)
           {
            case 1:
              node1 = 2 ;
              node2 = 3 ;
              node3 = 4 ;
              break;
            case 2:
              node1 = 1 ;
              node2 = 3 ;
              node3 = 4 ;
              break;
            case 3:
              node1 = 1 ;
              node2 = 2 ;
              node3 = 4 ;
              break;
            case 4:
              node1 = 1 ;
              node2 = 2 ;
              node3 = 3 ;
              break;
	   default:
              printf("IAM:%d, Bad Side:%d", iam, side);
              exit(EXIT_SUCCESS);
              break;
         }
         n1 = *(el_pt->ien + node1 - 1) - 1 ;  /* look up global node #'s */
         n2 = *(el_pt->ien + node2 - 1) - 1 ;
         n3 = *(el_pt->ien + node3 - 1) - 1 ;
         xx1 = *(global.coords + n1*nsd) ;     /* then get the global coords */
         yy1 = *(global.coords + n1*nsd + 1) ;
         zz1 = *(global.coords + n1*nsd + 2) ;
         xx2 = *(global.coords + n2*nsd) ;
         yy2 = *(global.coords + n2*nsd + 1) ;
         zz2 = *(global.coords + n2*nsd + 2) ;
         xx3 = *(global.coords + n3*nsd) ;
         yy3 = *(global.coords + n3*nsd + 1) ;
         zz3 = *(global.coords + n3*nsd + 2) ;

         ax = xx2 - xx1 ;    /* compute the edge vectors... */
         ay = yy2 - yy1 ;
         az = zz2 - zz1 ;
         bx = xx3 - xx1 ;
         by = yy3 - yy1 ;
         bz = zz3 - zz1 ;

         /* now create local planar coordinates for shape function use */
         
         xplan[0][0] = ZERO ;  /* node 1 */
         xplan[0][1] = ZERO ;
         
         xplan[1][0] = alength = sqrt(ax*ax+ay*ay+az*az) ;  /* node 2 */
         xplan[1][1] = ZERO ;
         
         adotb = ax*bx+ay*by+az*bz ;
         cx = bx - ax*adotb/(alength*alength) ;
         cy = by - ay*adotb/(alength*alength) ;
         cz = bz - az*adotb/(alength*alength) ;
         
         xplan[2][0] = adotb/alength ;  /* node 3 */
         xplan[2][1] = sqrt(cx*cx+cy*cy+cz*cz) ;
         
         /* now get the displacements of the face nodes... */
         
         ux1 = *(global.del_displ + n1*ndof) ;
         uy1 = *(global.del_displ + n1*ndof + 1) ;
         uz1 = *(global.del_displ + n1*ndof + 2) ;
         ux2 = *(global.del_displ + n2*ndof) ;
         uy2 = *(global.del_displ + n2*ndof + 1) ;
         uz2 = *(global.del_displ + n2*ndof + 2) ;
         ux3 = *(global.del_displ + n3*ndof) ;
         uy3 = *(global.del_displ + n3*ndof + 1) ;
         uz3 = *(global.del_displ + n3*ndof + 2) ;
         avex = (ux1 + ux2 + ux3)/THREE ;
         avey = (uy1 + uy2 + uy3)/THREE ;
         avez = (uz1 + uz2 + uz3)/THREE ;
         
         /* ...and calculate the reaction forces */
         
         if(code == 0) /* buoyancy in a rectilinear direction */
            {
			 reactx = -rho_g * (avex*up_vec[0]+avey*up_vec[1]+avez*up_vec[2])*up_vec[0] ;
			 reacty = -rho_g * (avex*up_vec[0]+avey*up_vec[1]+avez*up_vec[2])*up_vec[1] ;
			 reactz = -rho_g * (avex*up_vec[0]+avey*up_vec[1]+avez*up_vec[2])*up_vec[2] ;
            }
         else if(code == 1)  /* buoyancy spherically radial to an origin point */
            {
             centrx = (xx1 + xx2 + xx3)/THREE ;
             centry = (yy1 + yy2 + yy3)/THREE ;
             centrz = (zz1 + zz2 + zz3)/THREE ;
             local_up[0] = centrx - up_vec[0] ; /* up_vec here is actually the center location*/
             local_up[1] = centry - up_vec[1] ; /* of radial gravity in spherical cases */
             local_up[2] = centrz - up_vec[2] ;
             dist = sqrt(local_up[0]*local_up[0]+local_up[1]*local_up[1]
                          +local_up[2]*local_up[2]) ;
             local_up[0] /= dist ;
             local_up[1] /= dist ;
             local_up[2] /= dist ; /* now have a normalized unit vector locally "up" */
             
			 reactx = -rho_g * (avex*local_up[0]+avey*local_up[1]+avez*local_up[2])*local_up[0] ;
			 reacty = -rho_g * (avex*local_up[0]+avey*local_up[1]+avez*local_up[2])*local_up[1] ;
			 reactz = -rho_g * (avex*local_up[0]+avey*local_up[1]+avez*local_up[2])*local_up[2] ;
            }
            
     /* now do 2D Gaussian integration of the tractions */
            
	for(ipt=0 ; ipt<3 ; ipt++)
	   {
	    sh2 = sh2all + ipt*3*3 ;
	    det = detall +ipt ;
	    
		   rpt = r_3edge[ipt] ;
		   spt = s_3edge[ipt] ;
		   wt = wt_3edge ;
	    
/* next generate the 2D shape functions */

	   *(sh2 + 0*3) = ONE ;
	   *(sh2 + 0*3 + 1) = ZERO ; 
	   *(sh2 + 0*3 + 2) = rpt ;

	   *(sh2 + 1*3) = ZERO ;
	   *(sh2 + 1*3 + 1) = ONE ; 
	   *(sh2 + 1*3 + 2) = spt ;

	   *(sh2 + 2*3) = -ONE ;
	   *(sh2 + 2*3 + 1) = -ONE ; 
	   *(sh2 + 2*3 + 2) = ONE - rpt - spt ;

          for(i=0 ; i<2 ; i++)
             {
              for(j=0 ; j<2 ; j++)
                 {
              
/* note that the first argument (i) counts the physical space coordinate x, 
while the second (j) counts the parent space coordinate Tse. 
Hence xs[j][i] = dX(j) / dTse(i).
*/

                sum = ZERO ;
                for(n=0 ; n<3 ; n++)
					{
					 sum += xplan[n][j] * *(sh2+3*n+i) ;
					}
                 xs[j][i] = sum ;
                 }
             }

           *det = xs[0][0]*xs[1][1] - xs[0][1]*xs[1][0] ;
           
           if(*det < ZERO)
              {
               printf("Negative buoyface det: element=%d  ipt=%d\n",el_pt->nel,ipt) ;
               /* this should never happen ; just check in case */
              }
           
          for(i=0 ; i<3 ; i++)
             {
              hold = (xs[1][1] * *(sh2 + i*3) -
                      xs[1][0] * *(sh2 + i*3 + 1)) ;

              hold1 = (xs[0][0] * *(sh2 + i*3 + 1) -
                       xs[0][1] * *(sh2 + i*3)) ;

              *(sh2 + i*3)     = hold/ *det ;
              *(sh2 + i*3 + 1) = hold1/ *det ;
             }

        w_det = wt[ipt] * *det ;
        

		*(rhs + 3*(node1-1)    ) +=  w_det * reactx * *(sh2+(3*0)+2) ;
		*(rhs + 3*(node1-1) + 1) +=  w_det * reacty * *(sh2+(3*0)+2) ;
		*(rhs + 3*(node1-1) + 2) +=  w_det * reactz * *(sh2+(3*0)+2) ;

		*(rhs + 3*(node2-1)    ) +=  w_det * reactx * *(sh2+(3*1)+2) ;
		*(rhs + 3*(node2-1) + 1) +=  w_det * reacty * *(sh2+(3*1)+2) ;
		*(rhs + 3*(node2-1) + 2) +=  w_det * reactz * *(sh2+(3*1)+2) ;

		*(rhs + 3*(node3-1)    ) +=  w_det * reactx * *(sh2+(3*2)+2) ;
		*(rhs + 3*(node3-1) + 1) +=  w_det * reacty * *(sh2+(3*2)+2) ;
		*(rhs + 3*(node3-1) + 2) +=  w_det * reactz * *(sh2+(3*2)+2) ;


	   }
            
            
         
           }

      else if (3 == nsd && (TRILIN == info->type))
         {
/************************************************************************/
/*
This version implements flattened quad integrated tractions at the faces
 of brick elements.  
GL 10/16/09
*/
/************************************************************************/
         switch(side)
           {
            case 1:
              node1 = 1 ;
              node2 = 4 ;
              node3 = 3 ;
              node4 = 2 ;
              break;
            case 2:
              node1 = 5 ;
              node2 = 6 ;
              node3 = 7 ;
              node4 = 8 ;
              break;
            case 3:
              node1 = 1 ;
              node2 = 2 ;
              node3 = 6 ;
              node4 = 5 ;
              break;
            case 4:
              node1 = 2 ;
              node2 = 3 ;
              node3 = 7 ;
              node4 = 6 ;
              break;
            case 5:
              node1 = 3 ;
              node2 = 4 ;
              node3 = 8 ;
              node4 = 7 ;
              break;
            case 6:
              node1 = 1 ;
              node2 = 5 ;
              node3 = 8 ;
              node4 = 4 ;
              break;
         }
         
         n1 = *(el_pt->ien + node1 - 1) - 1 ;  /* look up global node #'s */
         n2 = *(el_pt->ien + node2 - 1) - 1 ;
         n3 = *(el_pt->ien + node3 - 1) - 1 ;
         n4 = *(el_pt->ien + node4 - 1) - 1 ;
         xx1 = *(global.coords + n1*nsd) ;     /* then get the global coords */
         yy1 = *(global.coords + n1*nsd + 1) ;
         zz1 = *(global.coords + n1*nsd + 2) ;
         xx2 = *(global.coords + n2*nsd) ;
         yy2 = *(global.coords + n2*nsd + 1) ;
         zz2 = *(global.coords + n2*nsd + 2) ;
         xx3 = *(global.coords + n3*nsd) ;
         yy3 = *(global.coords + n3*nsd + 1) ;
         zz3 = *(global.coords + n3*nsd + 2) ;
         xx4 = *(global.coords + n4*nsd) ;
         yy4 = *(global.coords + n4*nsd + 1) ;
         zz4 = *(global.coords + n4*nsd + 2) ;

         ax = xx2 - xx1 ;    /* compute the edge vectors... */
         ay = yy2 - yy1 ;
         az = zz2 - zz1 ;
         bx = xx4 - xx1 ;
         by = yy4 - yy1 ;
         bz = zz4 - zz1 ;
         cx = xx3 - xx1 ;
         cy = yy3 - yy1 ;
         cz = zz3 - zz1 ;
         

         /* now create local planar coordinates for shape function use */
         
         xplan[0][0] = ZERO ;  /* node 1 */
         xplan[0][1] = ZERO ;
         
         xplan[1][0] = alength = sqrt(ax*ax+ay*ay+az*az) ;  /* node 2 */
         xplan[1][1] = ZERO ;
         
         adotb = ax*bx+ay*by+az*bz ;
         
         dx = bx - ax*adotb/(alength*alength) ; /* d is a vector perp to a with b's projected length */
         dy = by - ay*adotb/(alength*alength) ;
         dz = bz - az*adotb/(alength*alength) ;
         
         xplan[3][0] = adotb/alength ;  /* node 4 */
         xplan[3][1] = dlength = sqrt(dx*dx+dy*dy+dz*dz) ;
         
         ex = ay*bz - az*by ;
         ey = az*bx - ax*bz ;
         ez = ax*by - ay*bx ;
         elength = sqrt(ex*ex+ey*ey+ez*ez) ;
         ex /= elength ;
         ey /= elength ;
         ez /= elength ;   /* e is a unit vector perp to the axb plane */
         
         cdote = cx*ex + cy*ey + cz*ez ;
         
         fx = cx - ex*cdote ; /* f is the c vector flattened into the axb plane */
         fy = cy - ey*cdote ;
         fz = cz - ez*cdote ;
         
         xplan[2][0] = (fx*ax+fy*ay+fz*az)/alength ;  /* node 3 */
         xplan[2][1] = (fx*dx+fy*dy+fz*dz)/dlength ;

         /* now get the displacements of the face nodes... */
         
         ux1 = *(global.del_displ + n1*ndof) ;
         uy1 = *(global.del_displ + n1*ndof + 1) ;
         uz1 = *(global.del_displ + n1*ndof + 2) ;
         ux2 = *(global.del_displ + n2*ndof) ;
         uy2 = *(global.del_displ + n2*ndof + 1) ;
         uz2 = *(global.del_displ + n2*ndof + 2) ;
         ux3 = *(global.del_displ + n3*ndof) ;
         uy3 = *(global.del_displ + n3*ndof + 1) ;
         uz3 = *(global.del_displ + n3*ndof + 2) ;
         ux4 = *(global.del_displ + n4*ndof) ;
         uy4 = *(global.del_displ + n4*ndof + 1) ;
         uz4 = *(global.del_displ + n4*ndof + 2) ;
         avex = (ux1 + ux2 + ux3 + ux4)/4.0 ;
         avey = (uy1 + uy2 + uy3 + uy4)/4.0 ;
         avez = (uz1 + uz2 + uz3 + uz4)/4.0 ;
         
         /* ...and calculate the reaction forces */
         
         if(code == 0) /* buoyancy in a rectilinear direction */
            {
			 reactx = -rho_g * (avex*up_vec[0]+avey*up_vec[1]+avez*up_vec[2])*up_vec[0] ;
			 reacty = -rho_g * (avex*up_vec[0]+avey*up_vec[1]+avez*up_vec[2])*up_vec[1] ;
			 reactz = -rho_g * (avex*up_vec[0]+avey*up_vec[1]+avez*up_vec[2])*up_vec[2] ;
            }
         else if(code == 1)  /* buoyancy spherically radial to an origin point */
            {
             centrx = (xx1 + xx2 + xx3 + xx4)/4.0 ;
             centry = (yy1 + yy2 + yy3 + yy4)/4.0 ;
             centrz = (zz1 + zz2 + zz3 + zz4)/4.0 ;
             local_up[0] = centrx - up_vec[0] ; /* up_vec here is actually the center location*/
             local_up[1] = centry - up_vec[1] ; /* of radial gravity in spherical cases */
             local_up[2] = centrz - up_vec[2] ;
             dist = sqrt(local_up[0]*local_up[0]+local_up[1]*local_up[1]
                          +local_up[2]*local_up[2]) ;
             local_up[0] /= dist ;
             local_up[1] /= dist ;
             local_up[2] /= dist ; /* now have a normalized unit vector locally "up" */
             
			 reactx = -rho_g * (avex*local_up[0]+avey*local_up[1]+avez*local_up[2])*local_up[0] ;
			 reacty = -rho_g * (avex*local_up[0]+avey*local_up[1]+avez*local_up[2])*local_up[1] ;
			 reactz = -rho_g * (avex*local_up[0]+avey*local_up[1]+avez*local_up[2])*local_up[2] ;
            }
            
     /* now do 2D Gaussian integration of the tractions */
            
	for(ipt=0 ; ipt<4 ; ipt++)
	   {
	    sh2 = sh2all + ipt*4*3 ;
	    det = detall +ipt ;
	    
		   xpt = x_2x2[ipt] ;
		   ypt = y_2x2[ipt] ;
		   xpt /= (real) sqrt(THREE) ;
		   ypt /= (real) sqrt(THREE) ;
		   wt = wt_2x2 ;
	    
/* next generate the 2D shape functions */
               *(sh2 + 0*3) = PT25 * (ypt - ONE) ;
               *(sh2 + 0*3 + 1) = PT25 * (xpt - ONE) ; 
               *(sh2 + 0*3 + 2) = PT25 * (ONE - xpt) * (ONE - ypt) ;

               *(sh2 + 1*3) = PT25 * (ONE - ypt) ;
               *(sh2 + 1*3 + 1) = (-PT25) * (ONE + xpt) ; 
               *(sh2 + 1*3 + 2) = PT25 * (ONE + xpt) * (ONE - ypt) ;

               *(sh2 + 2*3) = PT25 * (ONE + ypt) ;
               *(sh2 + 2*3 + 1) = PT25 * (ONE + xpt) ; 
               *(sh2 + 2*3 + 2) = PT25 * (ONE + xpt) * (ONE + ypt) ;

               *(sh2 + 3*3) = (-PT25) * (ONE + ypt) ;
               *(sh2 + 3*3 + 1) = PT25 * (ONE - xpt) ; 
               *(sh2 + 3*3 + 2) = PT25 * (ONE - xpt) * (ONE + ypt) ;

          for(i=0 ; i<2 ; i++)
             {
              for(j=0 ; j<2 ; j++)
                 {
              
/* note that the first argument (i) counts the physical space coordinate x, 
while the second (j) counts the parent space coordinate Tse. 
Hence xs[j][i] = dX(j) / dTse(i).
*/

                sum = ZERO ;
                for(n=0 ; n<4 ; n++)
					{
					 sum += xplan[n][j] * *(sh2+3*n+i) ;
					}
                 xs[j][i] = sum ;
                 }
             }

           *det = xs[0][0]*xs[1][1] - xs[0][1]*xs[1][0] ;
           
           if(*det < ZERO)
              {
               printf("Negative buoyface det: element=%d  ipt=%d\n",el_pt->nel,ipt) ;
               /* this should never happen ; just check in case */
              }
                      
          for(i=0 ; i<4 ; i++)
             {
              hold = (xs[1][1] * *(sh2 + i*3) -
                      xs[1][0] * *(sh2 + i*3 + 1)) ;

              hold1 = (xs[0][0] * *(sh2 + i*3 + 1) -
                       xs[0][1] * *(sh2 + i*3)) ;

              *(sh2 + i*3)     = hold/ *det ;
              *(sh2 + i*3 + 1) = hold1/ *det ;
             }

        w_det = wt[ipt] * *det ;
        

		*(rhs + 3*(node1-1)    ) +=  w_det * reactx * *(sh2+(3*0)+2) ;
		*(rhs + 3*(node1-1) + 1) +=  w_det * reacty * *(sh2+(3*0)+2) ;
		*(rhs + 3*(node1-1) + 2) +=  w_det * reactz * *(sh2+(3*0)+2) ;

		*(rhs + 3*(node2-1)    ) +=  w_det * reactx * *(sh2+(3*1)+2) ;
		*(rhs + 3*(node2-1) + 1) +=  w_det * reacty * *(sh2+(3*1)+2) ;
		*(rhs + 3*(node2-1) + 2) +=  w_det * reactz * *(sh2+(3*1)+2) ;

		*(rhs + 3*(node3-1)    ) +=  w_det * reactx * *(sh2+(3*2)+2) ;
		*(rhs + 3*(node3-1) + 1) +=  w_det * reacty * *(sh2+(3*2)+2) ;
		*(rhs + 3*(node3-1) + 2) +=  w_det * reactz * *(sh2+(3*2)+2) ;

		*(rhs + 3*(node4-1)    ) +=  w_det * reactx * *(sh2+(3*3)+2) ;
		*(rhs + 3*(node4-1) + 1) +=  w_det * reacty * *(sh2+(3*3)+2) ;
		*(rhs + 3*(node4-1) + 2) +=  w_det * reactz * *(sh2+(3*3)+2) ;

	   }
            
         
           }


  }  
/************************** end of buoy_form *********************************/ 



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
          )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine shape
** shape computes the shape functions and gradients for a single element.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

   {
    real  xs[3][3] , hold,hold1,hold2 ;
    real  *sh , *det ;
    int    i , j , ipt , rule , nsd , nsdp1 ;
    
    nsd = fe_sys.nsd ;
    nsdp1 = nsd + 1 ;
    for(ipt=0 ; ipt<nint ; ipt++)
       {
        sh = sh_pt + ipt*nen*nsdp1 ;
        det = det_pt + ipt ;
        
        switch( type )
           {
            case TET:
                  if(1 == el_pt->degen)  /* duo elt requires no shape nor determinant. */
                     {
                      continue;
                     }
                  for(i=0 ; i<nen ; i++)
                     {
                      if( all )
                         {
                          *(sh + i*nsdp1) = PSH_TET(0,i,ipt) ;
                          *(sh + i*nsdp1 + 1) = PSH_TET(1,i,ipt) ;
                          *(sh + i*nsdp1 + 2) = PSH_TET(2,i,ipt) ;
                         }
                      *(sh + i*nsdp1 + nsd) = PSH_TET(nsd,i,ipt) ;
                     }
    
               break ;

            case BILIN:
               for(i=0 ; i<nen ; i++)
                  {
                   if( all )
                      {
                       *(sh + i*nsdp1) = PSH_BILIN(0,i,ipt) ;
                       *(sh + i*nsdp1 + 1) = PSH_BILIN(1,i,ipt) ;
                      }
                   *(sh + i*nsdp1 + nsd) = PSH_BILIN(nsd,i,ipt) ;
                  }
               if(el_pt->degen)
                  {
                   for(j=0 ; j<nsdp1 ; j++)
                      {
                       *(sh + 2*nsdp1 + j) += *(sh + 3*nsdp1 + j) ;
                       *(sh + 3*nsdp1 + j) = ZERO ;
                      }
                  }
    
               break ;
    
            case TRILIN:
               for(i=0 ; i<nen ; i++)
                  {
                   if(6 == nint)
                      {
					   if( all )
						  {
						   *(sh + i*nsdp1) = PSH_TRILIN6(0,i,ipt) ;
						   *(sh + i*nsdp1 + 1) = PSH_TRILIN6(1,i,ipt) ;
						   *(sh + i*nsdp1 + 2) = PSH_TRILIN6(2,i,ipt) ;
						  }
					   *(sh + i*nsdp1 + nsd) = PSH_TRILIN6(nsd,i,ipt) ;
                      }
                   else if(8 == nint)
                      {
					   if( all )
						  {
						   *(sh + i*nsdp1) = PSH_TRILIN8(0,i,ipt) ;
						   *(sh + i*nsdp1 + 1) = PSH_TRILIN8(1,i,ipt) ;
						   *(sh + i*nsdp1 + 2) = PSH_TRILIN8(2,i,ipt) ;
						  }
					   *(sh + i*nsdp1 + nsd) = PSH_TRILIN8(nsd,i,ipt) ;
                      }
                  }
               if(1 == el_pt->degen)  /* triangular wedge */
                  {
                   for(j=0 ; j<nsdp1 ; j++)
                      {
                       *(sh + 2*nsdp1 + j) += *(sh + 3*nsdp1 + j) ;
                       *(sh + 3*nsdp1 + j) = ZERO ;
                       *(sh + 6*nsdp1 + j) += *(sh + 7*nsdp1 + j) ;
                       *(sh + 7*nsdp1 + j) = ZERO ;
                      }
                  }
               if(2 == el_pt->degen)  /* tetrahedron */
                  {
                   for(j=0 ; j<nsdp1 ; j++)
                      {
                       *(sh + 2*nsdp1 + j) += *(sh + 3*nsdp1 + j) ;
                       *(sh + 3*nsdp1 + j) = ZERO ;
                       *(sh + 4*nsdp1 + j) += 
                         (*(sh + 5*nsdp1 + j) + *(sh + 6*nsdp1 + j) + *(sh + 7*nsdp1 + j)) ;
                       *(sh + 5*nsdp1 + j) = ZERO ;
                       *(sh + 6*nsdp1 + j) = ZERO ;
                       *(sh + 7*nsdp1 + j) = ZERO ;
                      }
                  }
    
               break ;
    
            case SERENDIP:
               if(4 == nint) rule = 0 ;
               else if(9 == nint) rule = 1 ;
               for(i=0 ; i<nen ; i++)
                  {
                   if( all )
                      {
                       *(sh + i*nsdp1)       = PSH_SEREN(0,i,ipt,rule) ;
                       *(sh + i*nsdp1 + 1)   = PSH_SEREN(1,i,ipt,rule) ;
                      }
                   *(sh + i*nsdp1 + nsd) = PSH_SEREN(nsd,i,ipt,rule) ;
                  }
               if(el_pt->degen)
                  {
                   for(j=0 ; j<nsdp1 ; j++)
                      {
                       *(sh + 2*nsdp1 + j) += *(sh + 3*nsdp1 + j) +
                                              *(sh + 6*nsdp1 + j) ;
                       *(sh + 3*nsdp1 + j) = ZERO ;
                       *(sh + 6*nsdp1 + j) = ZERO ;
                      }
                  }
    
               break ;
           }
    
          if( !all )  continue ;
          
          for(i=0 ; i<nsd ; i++)
             {
              for(j=0 ; j<nsd ; j++)
/* note that the first argument counts the physical space coordinate x, 
while the second counts the parent space coordinate Tse. 
Hence xs[j][i] = dX(j) / dTse(i).
*/
                xs[j][i] = dotsh( global.coords , i , j ,
                                        sh , el_pt , info , nen , nsd , 0 , ELASTIC) ;
             }
    
          if(2 == nsd){

             *det = xs[0][0]*xs[1][1] - xs[0][1]*xs[1][0] ;

          } else if (3 == nsd) {

             *det = xs[0][0]*xs[1][1]*xs[2][2]
                   -xs[2][0]*xs[1][1]*xs[0][2]
                   +xs[0][1]*xs[1][2]*xs[2][0]
                   -xs[2][1]*xs[1][2]*xs[0][0]
                   +xs[0][2]*xs[1][0]*xs[2][1]
                   -xs[2][2]*xs[1][0]*xs[0][1];

          } else {
            attempt = DET_ERR;
            printf("NSD %d Not == 2 or 3!\n",nsd);
          }


          if( *det <= ZERO && el_pt->degen != 1)
             {
              attempt = DET_ERR ;
              printf("Bad determinant in element %d: %8.5g\n",
                      el_pt->nel,*det) ;

             }
/* jwp cutting out the determinant division: it doesn't work in 3d.
hence alter the 2d case below as well.
          else
             {
              for(i=0 ; i<nsd ; i++)
                 {
                  for(j=0 ; j<nsd ; j++)
                    xs[i][j] /= (*det) ;
                 }
             }
*/
/* 
C   compute gradient of shape fns w/respect to x
C                                           /     \ -1
C      dNj       dNj       dTse     dNj    |  dX   |
C     -----  =  -----  *  ------ = ----- * | ----- |
C      dX        dTse      dX       dTse   |  dTse |
C                                           \     /
C
*/
/* this code corrected as per JWP bug fix of Sept. 2000 */

          if(2 == nsd){

          for(i=0 ; i<nen ; i++)
             {
              hold = (xs[1][1]* *(sh + i*nsdp1) -
                     xs[1][0]* *(sh + i*nsdp1 + 1)) ;

              hold1 = (xs[0][0]* *(sh + i*nsdp1 + 1) -
                   xs[0][1]* *(sh + i*nsdp1)) ;

              *(sh + i*nsdp1) = hold/ *det ;
              *(sh + i*nsdp1 + 1) = hold1/ *det ;
             }

          } else if (3 == nsd){

          for(i=0 ; i<nen ; i++)
             {
              hold = (xs[1][1]*xs[2][2] - xs[1][2]*xs[2][1]) * *(sh + i*nsdp1) +
                ( xs[1][2]*xs[2][0] - xs[1][0]*xs[2][2]) * *(sh + i*nsdp1 + 1) +
                ( xs[1][0]*xs[2][1] - xs[1][1]*xs[2][0]) * *(sh + i*nsdp1 + 2) ;

              hold1= (xs[0][2]*xs[2][1] - xs[0][1]*xs[2][2]) * *(sh + i*nsdp1) +
                ( xs[0][0]*xs[2][2] - xs[0][2]*xs[2][0]) * *(sh + i*nsdp1 + 1) +
                ( xs[0][1]*xs[2][0] - xs[0][0]*xs[2][1]) * *(sh + i*nsdp1 + 2) ;

              hold2= (xs[0][1]*xs[1][2] - xs[0][2]*xs[1][1]) * *(sh + i*nsdp1) +
                ( xs[0][2]*xs[1][0] - xs[0][0]*xs[1][2]) * *(sh + i*nsdp1 + 1) +
                ( xs[0][0]*xs[1][1] - xs[0][1]*xs[1][0]) * *(sh + i*nsdp1 + 2) ;

              *(sh + i*nsdp1) = hold/ *det ;
              *(sh + i*nsdp1 + 1) = hold1/ *det ;
              *(sh + i*nsdp1 + 2) = hold2/ *det ;
             }

          } else {
            attempt = DET_ERR;
            printf("NSD %d Not == 2 or 3!\n",nsd);
          }
       }
   }
/************************** end of shape *************************************/ 



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
          )

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

   {
    real  sum = ZERO ;
    real  x ;
    int    n , node , nsd , nsdp1 ;
    
    nsd = fe_sys.nsd ;
    nsdp1 = nsd + 1 ;
    for(n=0 ; n<nen ; n++)
        {
         node = *(el_pt->ien + n) - 1 ;
         x = *(array + node*ndim + j) ;

         if( st_flag && time_data.do_slip && *(el_pt->is_split) )
            x += adfldp( j , el_pt , node , info , task_code ) ;
         
         sum += x * SH(i,n) ;
        }
    return( sum ) ;
   }
/************************** end of dotsh *************************************/ 



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
           )

/*  RETURN VALUE:  the required additional nodal displacement */
/*  DESCRIPTION:   */
/*
** Function adfldp
** adfldp returns to dotsh() the added nodal displacements
** needed to account for split node fault slip in calculating strain
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

   {
    int    i , ndof , total , index , nel , ifltgrp , active ;
    int    *slip_list , *slip_grp ;
    real   *slip_val ;
    
    ndof = fe_sys.ndof ;
    nel = el_pt->nel ;
    total = *(el_pt->is_split) ;
    slip_list = info -> slip_list ;
    slip_val = info -> slip_val ;
    slip_grp = info -> slip_grp ;
    for(i=0 ; i<total ; i++)
        {
         index = *(el_pt->is_split + i + 1) ;
         ifltgrp = *(slip_grp+index) - 1 ; /* zero based fault group */
         active = (fltgrp_ptr+ifltgrp)->due_now ;
         /* slip_list[2*index] is the zero-based element number.
            slip_list[2*index+1] is the 1-based node number.
            nel is (here) the one-based element (from el_pt->nel)
            node is the zero-based node number: calling dotsh created it as "ien" - 1
          */
         if( slip_list[2*index]   == (nel - 1)  &&  /* testing zero-based slip_list entry vs. one-based nel */
             slip_list[2*index+1] == (node + 1) &&  /* testing one-based slip_list entry vs. zero-based node*/
              active )
           {
            if(ELASTIC == task_code)
               return(  *(slip_val + ndof*index + j)  );
               /* simply add the fault slip  */

            else if(VISCO == task_code && active == 2)
               return(  *(slip_val + ndof*index + j) * time_data.dt  );
               /* interpret as slip RATE; multiply by dt  */

            else if(VISCO == task_code && active == 1)
               return(  *(slip_val + ndof*index + j)  );
               /* scheduled slip event, just add  */

            else if(QUAKEEVT == task_code)
               return(  *(slip_val + ndof*index + j) * (fltgrp_ptr+ifltgrp)->q_amount  );
               /* quake event; modify nominal slip by trial amount  */
           }
        }
    return( ZERO ) ;
   }
/************************** end of adfldp ************************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   p_shape   ***************
*/

/*  INTERFACE:   */

     void
     p_shape(
             GROUP   *grp_ptr     /* pointer to current element group */
            )

/*  RETURN VALUE:  - none - */
/*  DESCRIPTION:   */
/*
** Routine p_shape
** p_shape computes the parent-space shape function for common element types.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

   {
    int  ipt , rule , nen , nint , nel , numel , type , nsd , nsdp1 , n ;
    real  xpt , ypt, zpt ;
    ELEMENT_INFO  *info ;
    ELEMENT_DATA   *el_pt ;
    
/*  squawk("Doing parent functions...\n") ;  */

    nsd = fe_sys.nsd ;
    nsdp1 = nsd + 1 ;
    switch((grp_ptr->el_info)->type)
       {
        case QUAKE:
           return  ; /* no shape functions for quake elements */
           break ;

        case TET:
           for(ipt = 0; ipt < 4; ipt++)
              {
               xpt = x_tet[ipt] ;
               ypt = y_tet[ipt] ;
               zpt = z_tet[ipt] ;

/*
C
C  3D parent coordinates : (r,s,t)
C     define u = 1 - r - s - t
C  Shape functions and derivatives :
C
C      |               |  pN  |  pN  |  pN              s
C      |       N       | ---- | ---- | ----             |
C      |               |  pr  |  ps  |  pt              3
C -----|---------------|------|------|------            |
C   1  | 1 - r - s - t |  -1  |  -1  |  -1
C -----|---------------|------|------|------            |
C   2  |       r       |  1   |  0   |  0               |
C -----|---------------|------|------|------            |
C   3  |       s       |  0   |  1   |  0               1-------- 2 --r
C -----|---------------|------|------|------           /
C   4  |       t       |  0   |  0   |  1            /
C -----|---------------|------|------|------       4
C                                                /
C                                              t
*/
               PSH_TET(0,0,ipt) = -ONE ;
               PSH_TET(1,0,ipt) = -ONE ;
               PSH_TET(2,0,ipt) = -ONE ;
               PSH_TET(nsd,0,ipt) = ONE - xpt - ypt - zpt;

               PSH_TET(0,1,ipt) = ONE ;
               PSH_TET(1,1,ipt) = ZERO ;
               PSH_TET(2,1,ipt) = ZERO ;
               PSH_TET(nsd,1,ipt) = xpt ;

               PSH_TET(0,2,ipt) = ZERO  ;
               PSH_TET(1,2,ipt) = ONE;
               PSH_TET(2,2,ipt) = ZERO ;
               PSH_TET(nsd,2,ipt) = ypt ;

               PSH_TET(0,3,ipt) = ZERO  ;
               PSH_TET(1,3,ipt) = ZERO  ;
               PSH_TET(2,3,ipt) = ONE;
               PSH_TET(nsd,3,ipt) = zpt ;
              }

           break ;

        case BILIN:
           for(ipt=0 ; ipt<4 ; ipt++)
              {
               xpt = x_2x2[ipt] ;
               ypt = y_2x2[ipt] ;
               xpt /= (real) sqrt(THREE) ;
               ypt /= (real) sqrt(THREE) ;
               
               PSH_BILIN(0,0,ipt) = PT25 * (ypt - ONE) ;
               PSH_BILIN(1,0,ipt) = PT25 * (xpt - ONE) ; 
               PSH_BILIN(nsd,0,ipt) = PT25 * (ONE - xpt) * (ONE - ypt) ;

               PSH_BILIN(0,1,ipt) = PT25 * (ONE - ypt) ;
               PSH_BILIN(1,1,ipt) = (-PT25) * (ONE + xpt) ; 
               PSH_BILIN(nsd,1,ipt) = PT25 * (ONE + xpt) * (ONE - ypt) ;

               PSH_BILIN(0,2,ipt) = PT25 * (ONE + ypt) ;
               PSH_BILIN(1,2,ipt) = PT25 * (ONE + xpt) ; 
               PSH_BILIN(nsd,2,ipt) = PT25 * (ONE + xpt) * (ONE + ypt) ;

               PSH_BILIN(0,3,ipt) = (-PT25) * (ONE + ypt) ;
               PSH_BILIN(1,3,ipt) = PT25 * (ONE - xpt) ; 
               PSH_BILIN(nsd,3,ipt) = PT25 * (ONE - xpt) * (ONE + ypt) ;
              }

           break ;

        case TRILIN:
             for(ipt=0 ; ipt<6 ; ipt++)
              {
               for(n=0 ; n<8 ; n++)
                  {
                   PSH_TRILIN6(0,n,ipt) = ra[n]*(0.5+sa[n]*y_tri[ipt])*(0.5+ta[n]*z_tri[ipt]);
                   PSH_TRILIN6(1,n,ipt) = sa[n]*(0.5+ra[n]*x_tri[ipt])*(0.5+ta[n]*z_tri[ipt]);
                   PSH_TRILIN6(2,n,ipt) = ta[n]*(0.5+ra[n]*x_tri[ipt])*(0.5+sa[n]*y_tri[ipt]);
                   PSH_TRILIN6(nsd,n,ipt) =
                     (0.5+ra[n]*x_tri[ipt])*(0.5+sa[n]*y_tri[ipt])*(0.5+ta[n]*z_tri[ipt]);
                  }
              }

             for(ipt=0 ; ipt<8 ; ipt++)
              {
               xpt = x_2x2x2[ipt] ;
               ypt = y_2x2x2[ipt] ;
               zpt = z_2x2x2[ipt] ;
               xpt /= (real) sqrt(THREE) ;
               ypt /= (real) sqrt(THREE) ;
               zpt /= (real) sqrt(THREE) ;
               
               for(n=0 ; n<8 ; n++)
                  {
                   PSH_TRILIN8(0,n,ipt) = ra[n]*(0.5+sa[n]*ypt)*(0.5+ta[n]*zpt);
                   PSH_TRILIN8(1,n,ipt) = sa[n]*(0.5+ra[n]*xpt)*(0.5+ta[n]*zpt);
                   PSH_TRILIN8(2,n,ipt) = ta[n]*(0.5+ra[n]*xpt)*(0.5+sa[n]*ypt);
                   PSH_TRILIN8(nsd,n,ipt) =
                     (0.5+ra[n]*xpt)*(0.5+sa[n]*ypt)*(0.5+ta[n]*zpt);
                  }
              }

           break ;

        case SERENDIP:
           for(rule=0 ; rule<2 ; rule++)
               {
                for(ipt=0 ; ipt<9 ; ipt++)
                  {
                   if(rule == 0 && ipt >= 4) continue ;
                   
                   if(rule == 0)
                      {
                       xpt = x_2x2[ipt] ;
                       ypt = y_2x2[ipt] ;
                       xpt /= (real) sqrt(THREE) ;
                       ypt /= (real) sqrt(THREE) ;
                      }
                   else if(rule == 1)
                      {
                       xpt = x_3x3[ipt] ;
                       ypt = y_3x3[ipt] ;
                       xpt *= (real) sqrt(THREE/FIVE) ;
                       ypt *= (real) sqrt(THREE/FIVE) ;
                      }
               
                   PSH_SEREN(0,4,ipt,rule) = xpt * (ypt - ONE) ;
                   PSH_SEREN(1,4,ipt,rule) = PT5 * (xpt*xpt - ONE) ;
                   PSH_SEREN(nsd,4,ipt,rule) = PT5 * (ONE - xpt*xpt) *
                                               (ONE - ypt) ;
               
                   PSH_SEREN(0,5,ipt,rule) = PT5 * (ONE - ypt*ypt) ;
                   PSH_SEREN(1,5,ipt,rule) = (-ypt) * (ONE + xpt) ;
                   PSH_SEREN(nsd,5,ipt,rule) = PT5 * (ONE + xpt) *
                                               (ONE - ypt*ypt) ;
               
                   PSH_SEREN(0,6,ipt,rule) = (-xpt) * (ONE + ypt) ;
                   PSH_SEREN(1,6,ipt,rule) = PT5 * (ONE - xpt*xpt) ;
                   PSH_SEREN(nsd,6,ipt,rule) = PT5 * (ONE - xpt*xpt) *
                                               (ONE + ypt) ;
               
                   PSH_SEREN(0,7,ipt,rule) = PT5 * (ypt*ypt - ONE) ;
                   PSH_SEREN(1,7,ipt,rule) = ypt * (xpt - ONE) ;
                   PSH_SEREN(nsd,7,ipt,rule) = PT5 * (ONE - xpt) *
                                               (ONE - ypt*ypt) ;
               
                   PSH_SEREN(0,0,ipt,rule) = PT25 * (ypt - ONE) -
                         PT5 * (PSH_SEREN(0,7,ipt,rule) +
                         PSH_SEREN(0,4,ipt,rule)) ;
                   PSH_SEREN(1,0,ipt,rule) = PT25 * (xpt - ONE) -
                         PT5 * (PSH_SEREN(1,7,ipt,rule) +
                         PSH_SEREN(1,4,ipt,rule)) ;
                   PSH_SEREN(nsd,0,ipt,rule) = PT25 * (ONE - xpt) *
                                               (ONE - ypt) -
                         PT5 * (PSH_SEREN(nsd,7,ipt,rule) +
                         PSH_SEREN(nsd,4,ipt,rule)) ;
               
                   PSH_SEREN(0,1,ipt,rule) = PT25 * (ONE - ypt) -
                         PT5 * (PSH_SEREN(0,4,ipt,rule) +
                         PSH_SEREN(0,5,ipt,rule)) ;
                   PSH_SEREN(1,1,ipt,rule) = (-PT25) * (ONE + xpt) -
                         PT5 * (PSH_SEREN(1,4,ipt,rule) +
                         PSH_SEREN(1,5,ipt,rule)) ;
                   PSH_SEREN(nsd,1,ipt,rule) = PT25 * (ONE + xpt) *
                                               (ONE - ypt) -
                         PT5 * (PSH_SEREN(nsd,4,ipt,rule) +
                         PSH_SEREN(nsd,5,ipt,rule)) ;
               
                   PSH_SEREN(0,2,ipt,rule) = PT25 * (ONE + ypt) -
                         PT5 * (PSH_SEREN(0,5,ipt,rule) +
                         PSH_SEREN(0,6,ipt,rule)) ;
                   PSH_SEREN(1,2,ipt,rule) = PT25 * (ONE + xpt) -
                         PT5 * (PSH_SEREN(1,5,ipt,rule) +
                         PSH_SEREN(1,6,ipt,rule)) ;
                   PSH_SEREN(nsd,2,ipt,rule) = PT25 * (ONE + xpt) *
                                               (ONE + ypt) -
                         PT5 * (PSH_SEREN(nsd,5,ipt,rule) +
                         PSH_SEREN(nsd,6,ipt,rule)) ;
               
                   PSH_SEREN(0,3,ipt,rule) = (-PT25) * (ONE + ypt) -
                         PT5 * (PSH_SEREN(0,6,ipt,rule) +
                         PSH_SEREN(0,7,ipt,rule)) ;
                   PSH_SEREN(1,3,ipt,rule) = PT25 * (ONE - xpt) -
                         PT5 * (PSH_SEREN(1,6,ipt,rule) +
                         PSH_SEREN(1,7,ipt,rule)) ;
                   PSH_SEREN(nsd,3,ipt,rule) = PT25 * (ONE - xpt) *
                                               (ONE + ypt) -
                         PT5 * (PSH_SEREN(nsd,6,ipt,rule) +
                         PSH_SEREN(nsd,7,ipt,rule)) ;
               
                  }
               }

           break ;
       }

    if( fe_sys.save_shape )
       {
        info = grp_ptr->el_info ;
        nen     = info->nen ;
        numel     = info->numel ;
        type      = info->type ;

/*  printf("Allocating memory for shape functions...\n") ;  */

        for( nel = 0 ; nel < numel ; nel++)
            {
             el_pt = grp_ptr->el_data + nel ;
             nint = el_pt->nint ;
        
             el_pt->shape =
                  (real *) calloc(nsdp1*nen*nint,sizeof(real)) ;
             el_pt->deter = (real *) calloc(nint,sizeof(real)) ;

             if(el_pt->shape == NULL || el_pt->deter == NULL)
                attempt = SH_MEM ;
            }
        completion() ;

/*  printf("Doing daughter functions...\n") ;  */

        for( nel = 0 ; nel < numel ; nel++)
            {
             el_pt = grp_ptr->el_data + nel ;
             nint = el_pt->nint ;

             shape( el_pt , grp_ptr->el_info , el_pt->shape , el_pt->deter , type ,
                    nen , nint , 1 ) ;
            }
       }
       
   }
/************************** end of p_shape ***********************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   form_equil   ***************
*/

/*  INTERFACE:   */
     void
     form_equil(
                GROUP   *grp_ptr    /* pointer to current element group */
               )

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

     {
      int     numel , nel , nsd , nsdp1 , ndof , nee , nen , str_size ,
              save_shape , i , j , jt2 , jt3 , ipt , nint , type ;
      real   *wt , *shape_ptr , *det_ptr , *sh , *btemp ;
      real   w_det ;

      ELEMENT_DATA   *el_pt ;
      ELEMENT_INFO  *info ;
      ELEMENT_MAT   *mat_pt ;


      nsd = fe_sys.nsd ;
      nsdp1 = nsd + 1 ;
      numel = (grp_ptr->el_info)->numel ;
      mat_pt = grp_ptr->el_mat ;
      info = grp_ptr->el_info ;
      ndof  = fe_sys.ndof ;
      nee     = info->nee ;
      nen     = info->nen ;
      str_size = info->stress_size ;
      save_shape = fe_sys.save_shape ;
      type = (grp_ptr->el_info)->type ;
      
      switch(type)
          {
           case TET:
           case BILIN:
           case SERENDIP:
           case TRILIN:
              for ( nel = 0 ; nel < numel ; nel++)
                   {
                    el_pt = grp_ptr->el_data + nel ;
                    if((grp_ptr->el_info)->type == TET && el_pt->degen == 1) continue ;
                    /* we don't want to do rhs calculations on duo elements! */

                        
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

				    for (  i = 0 ; i < nee ; i++ )
					    el_rhs[i] = ZERO ;
          
					  if( !save_shape )
						 {
						  shape_ptr = sh_temp ;
						  det_ptr = det_temp ;
						  shape( el_pt , info , shape_ptr , det_ptr ,
						         type , nen , nint , 1 ) ;
						 }
					  else
						 {
						  shape_ptr = el_pt -> shape ;
						  det_ptr = el_pt -> deter ;
						 }

					  for (ipt = 0 ; ipt < nint ; ipt++)
						   {
							sh = shape_ptr + ipt*nen*nsdp1 ;
							w_det = wt[ipt] * *(det_ptr + ipt) ;
							/* btemp here points to stress */
							btemp = (el_pt->stress)+ipt*str_size ;
							
							switch( ndof )
							   {
								case 1:

								  break ;
								   
								case 2:

								   break ;
								   
								case 3:
								   if(2 == nsd)
								   {
								   for(j=0;j<nen;j++)
									  {
									   jt3 = 3*j ;
			
									   *(el_rhs + jt3)     -= w_det * 
										 (btemp[0] * SH(0,j) + btemp[3] * SH(1,j)) ;
									   *(el_rhs + jt3 + 1) -= w_det *
										 (btemp[1] * SH(1,j) + btemp[3] * SH(0,j)) ;
									   *(el_rhs + jt3 + 2) -= w_det *
										 (btemp[4] * SH(0,j) + btemp[5] * SH(1,j)) ;
									  }
								   } else if (3 == nsd)
								   {
								   for(j=0;j<nen;j++)
									  {
									   jt3 = 3*j ;
			
									   *(el_rhs + jt3)     -= w_det * (btemp[4] * SH(2,j) +
										 btemp[0] * SH(0,j) + btemp[3] * SH(1,j)) ;
									   *(el_rhs + jt3 + 1) -= w_det * (btemp[5] * SH(2,j) +
										 btemp[1] * SH(1,j) + btemp[3] * SH(0,j)) ;
									   *(el_rhs + jt3 + 2) -= w_det * (btemp[2] * SH(2,j) +
										 btemp[4] * SH(0,j) + btemp[5] * SH(1,j)) ;
									  }
								   }
								   break ;
								   
							   }

						   } /* end loop over integration points */

                    addfor( force.full_rhs , el_rhs , grp_ptr->el_info ,
                            el_pt ) ;

                   }
               break ;
               
          }
     }

/************************** end of form_equil *********************************/ 



