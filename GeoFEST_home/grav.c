/*
***                           File grav.c                       ***
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
***          == grav_calc ==
***          == dgrav_form ==
***          == dgrav_flux ==
***/


#define EXTERN extern
#include "grav.h"
#include "stiff.h"

static real
    sh_temp[256] , det_temp[9] ;




/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   grav_calc   ***************
*/

/*  INTERFACE:   */

      void
      grav_calc(
                GROUP		*grp_ptr , /* pointer to current element group */
                int         task
               )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine grav_calc
** grav_calc computes and outputs gravity change information.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

     {
      int g_numel , nint , nel , use_list , n , str_size , inode1 , inode2 ,
	  ipt , j , ng , ntot , nsd , nsdp1 , mine , prt_n , count , index , numface ,
	  numel , nen , bnel , side , dnumface , nn , dnel , dside , save_shape , ndof ;
      real   xpt[3][10] , str_ctr[9] , midpt[3] , fmid[3] , dvol , w_det ,
             ed , dmass , density , big_g , delta_rho , little_g , dgnet ;
      real   *shape_ptr , *det_ptr , *sh , *dgrav , *up_ptr , *wt , rho_g , dN ;
      ELEMENT_DATA   *el_pt , *bel_pt , *del_pt ;
      ELEMENT_MAT   *mat_pt ;
      BUOY_DATA      *buoy_ptr , *receive_ptr , *donor_ptr ;
      /* MPI_Status  status ; */

      int out_num, tot_num, i, ii, k, s_size, tmp_sum, indx , rad_code ,
          inode ;
      int  *buoy_list , *dbuoy_list ;
      int *rcounts, *displs;
      /* OUT_ELEMENT_DATA *s_data, *r_data; */


      nsd = fe_sys.nsd ;  nsdp1 = nsd + 1 ;
      ndof = fe_sys.ndof ;
      g_numel = (grp_ptr->el_info)->g_numel ;
      str_size = (grp_ptr->el_info)->stress_size ;
      numel = (grp_ptr->el_info)->numel ;
      nen     = (grp_ptr->el_info)->nen ;
      mat_pt = grp_ptr->el_mat ;
      save_shape = fe_sys.save_shape;

     if (0 == iam)
	{
	 fprintf(out_file,
	   "\n Gravity changes  -- element group #%d \n",grp_ptr->group_num) ;
	 fprintf(out_file,"\n Simulation time = %f ; step size = %g\n",
	   time_data.time,time_data.dt) ;
	 fprintf(out_file,
      "\n coordinates                          accelerations . . .->      geoid change & gradients\n") ;
     if(str_size == 6)
	 fprintf(out_file,
      "   X             Y            Z       dg_x        dg_y      dg_z      dg_net        dN    P_zx     P_zy     P_zz\n") ;     
	fflush(out_file) ;
	}

/*
The structure of this coding is based on first looping over each element and computing the 
dilatation for each one.  This need not be saved in element storage, as it is used immediately.
Then nested within the element loop is a loop over buoyancy facets, modeled after the coding
in form_rhs().  Within this sub-loop, each facet accumulates the gravity contributions from
the enclosing loop element as just computed.

At the completion of the buoy facet loop, the buoy structure storage should contain a preliminary
total of the delta_g components and the delta_potl values at the centroid of each facet.  What
remains is to compute the contributions from the boundary mass fluxes through each facet and
the final computation of the net g change and net geoid change.  For this last step, there is a
double nested loop over buoy facets.  The outer loop gathers the contributions from each of the 
inner loop facets (excluding self-contribution) and then before looping back, does the final 
g and geoid calculation for that facet.
*/
      
      switch((grp_ptr->el_info)->type)
          {
           case TET:
           case TRILIN:
           
    /* first initialize to zero the storage that will receive
        gravity contributions at the surface facets! */
        
		for (i=0 ; i < (grp_ptr->el_info)->numbuoy ; i++)  /* loop over facets */
			   {
				buoy_ptr = (grp_ptr->el_info)->buoy + i ;
				numface = buoy_ptr->numface ;
				
				for(n=0 ; n<numface ; n++)
				   {
					dgrav = buoy_ptr->dgrav + 7*n ; /* contains 3 components of dg and scalar dV and 3 gradients */
					for(ii=0 ; ii < 7 ; ii++)
						{
						 dgrav[ii] = ZERO ;
						}
				   }
			   }   /* numbuoy loop */
      /*   initialization complete    */
           
           
           
           for ( nel = 0 ; nel < numel ; nel++)
                   {
                    el_pt = grp_ptr->el_data + nel ;
                    big_g  = *(mat_pt->big_g + el_pt->mat) ;
                    little_g  = *(mat_pt->little_g + el_pt->mat) ;
                    density  = *(mat_pt->density + el_pt->mat) ;
                    nint = el_pt->nint ;
	            	shape_ptr = (save_shape)? el_pt->shape: sh_temp ;
	            	det_ptr = (save_shape)? el_pt->deter: det_temp ;
	            	
	     /* filter out any duo elements for special handling */
	                if((grp_ptr->el_info)->type == TET && el_pt->degen == 1)
	                   {
	                    continue ;
	                   }

			  switch( nint )
				 {
				  case 4:
					 wt = ((BILIN == (grp_ptr->el_info)->type) ? wt_2x2 : wt_tet) ;
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
                          shape( el_pt , grp_ptr->el_info , shape_ptr , det_ptr ,
                                (grp_ptr->el_info)->type , nen , nint , 1 ) ;
                         }

                      dvol = ZERO ;
                      
                      for (ipt = 0 ; ipt < nint ; ipt++)
                           {
                            sh = shape_ptr + ipt*nen*nsdp1 ;
                            w_det = wt[ipt] * det_ptr[ipt];
                    /*      dilat = (el_pt -> dilat) + ipt ; */

                            switch( ndof )
                               {
                                case 2:
                                   ed = dotsh( global.displ ,
                                            0 , 0 , sh , el_pt , (grp_ptr->el_info) , nen , ndof , 1 , task ) +
                                        dotsh( global.displ ,
                                            1 , 1 , sh , el_pt , (grp_ptr->el_info) , nen , ndof , 1 , task ) ;
                          /*       *dilat = ed * lambda ;  */

                                   break ;
                                   
                                case 3:
                                   if(2 == nsd) 
                                   {
                                      ed = dotsh( global.displ ,
                                            0 , 0 , sh , el_pt , (grp_ptr->el_info) , nen , ndof , 1 , task ) +
                                           dotsh( global.displ ,
                                            1 , 1 , sh , el_pt , (grp_ptr->el_info) , nen , ndof , 1 , task ) ;
                               /*     *dilat = ed * lambda ;  */

                                   } else if (3 == nsd){

                                      ed = dotsh( global.displ ,
                                              0 , 0 , sh , el_pt , (grp_ptr->el_info) , nen , ndof , 1 , task ) +
                                           dotsh( global.displ ,
                                              1,  1 , sh , el_pt , (grp_ptr->el_info) , nen , ndof , 1 , task ) +
                                           dotsh( global.displ ,
                                              2,  2 , sh , el_pt , (grp_ptr->el_info) , nen , ndof , 1 , task ) ;
                               /*     *dilat = ed * w_det ;     */
                                      dvol += ed * w_det ;

                                   }
                                   break ;
                                   
                               }
                                 
                           }  /* integration loop -ipt- */
                           
                      dmass = -dvol*density*big_g ; /* Eulerian mass change of element */

        /* Compute midpoint of contributing element */
					  midpt[0] = ZERO ;
					  midpt[1] = ZERO ;
					  midpt[2] = ZERO ;

                      for(i=0 ; i<nen ; i++)
                         {
                          inode = *(el_pt->ien + i) - 1 ;
                          midpt[0] += *(global.coords + inode*nsd) / ((real) nen) ;
                          midpt[1] += *(global.coords + inode*nsd + 1) / ((real) nen) ;
                          midpt[2] += *(global.coords + inode*nsd + 2) / ((real) nen) ;
                         }

                      for (i=0 ; i < (grp_ptr->el_info)->numbuoy ; i++)  /* loop over facets */
                           {
							buoy_ptr = (grp_ptr->el_info)->buoy + i ;
							buoy_list = buoy_ptr->buoy_list ;
							numface = buoy_ptr->numface ;
							rho_g = buoy_ptr->rho_g ;
							rad_code = buoy_ptr->rad_flag ;
                            
							for(n=0 ; n<numface ; n++)
							   {
								bnel = *(buoy_list + 2*n) ;
								side = (int) *(buoy_list + 2*n + 1) ;
								dgrav = buoy_ptr->dgrav + 7*n ; /* contain 3 components of dg and scalar dV plus 3 gradients */
								
								if (side <= 0) {
							printf("\tIAM:%d,BadSide:%d,nel:%d,numface=%d,first el=%d\nfirst side=%d\n\n",
								   iam, side, nel, numface,*buoy_list,*(buoy_list+1)) ;  }
								   
								up_ptr = buoy_ptr->upvec ;
								bel_pt = grp_ptr->el_data + bnel ;
								dgrav_form( dgrav , grp_ptr->el_info ,
										   side , dmass , midpt , bel_pt ) ;
		
							   }
                            
                           }   /* numbuoy loop */
                           
                           
                   } /* element loop -nel- */

    /*  final facet double-loop to finalize gravity contributions */
    
                        for (i=0 ; i < (grp_ptr->el_info)->numbuoy ; i++)  /* outer loop over facets */
                           {
							receive_ptr = (grp_ptr->el_info)->buoy + i ;
							buoy_list = receive_ptr->buoy_list ;
							numface = receive_ptr->numface ;
							rad_code = receive_ptr->rad_flag ;
                            
							for(n=0 ; n<numface ; n++)
							   {
								bnel = *(buoy_list + 2*n) ;
								el_pt = grp_ptr->el_data + bnel ;
								side = (int) *(buoy_list + 2*n + 1) ;
								dgrav = receive_ptr->dgrav + 7*n ; /* contain 3 components of dg and scalar dV plus 3 gradients */
								up_ptr = receive_ptr->upvec ;


                                for (ii=0 ; ii < (grp_ptr->el_info)->numbuoy ; ii++)  /* inner loop over facets */
                                   {
                                   	donor_ptr = (grp_ptr->el_info)->buoy + ii ;
							        dbuoy_list = donor_ptr->buoy_list ;
							        dnumface = donor_ptr->numface ;
							        delta_rho = donor_ptr->delta_rho ;

                                    for(nn=0 ; nn<dnumface ; nn++)
                                       {
                                        if(i==ii && n==nn) continue ;
                                        dnel = *(dbuoy_list + 2*nn) ;
                                        del_pt = grp_ptr->el_data + dnel ;
                                        dside = (int) *(dbuoy_list + 2*nn + 1) ;
                                        dgrav_flux( dgrav , grp_ptr->el_info ,
										   side , dside , delta_rho , up_ptr , big_g , del_pt , el_pt , fmid ) ;
                                       }
                                   }

           /* now collect the net results... */
           
       dgnet = little_g * ( sqrt((dgrav[0]/little_g - up_ptr[0])*(dgrav[0]/little_g - up_ptr[0]) +
                                 (dgrav[1]/little_g - up_ptr[1])*(dgrav[1]/little_g - up_ptr[1]) +
                                 (dgrav[2]/little_g - up_ptr[2])*(dgrav[2]/little_g - up_ptr[2])
                                ) - ONE ) ;  /* the change in the magnitude of local g */

       dN = - dgrav[3]/little_g ;   /* the change in geoid height along local vertical */
       
           /* print the gravity results... */
           
           
      fprintf(out_file, "%g  %g  %g      %g  %g  %g     %g     %g    %g     %g     %g\n",
              fmid[0],fmid[1],fmid[2],dgrav[0],dgrav[1],dgrav[2],dgnet,dN,dgrav[4],dgrav[5],dgrav[6]) ;     

		
							   }
                            
                           }


           break;
          }    /* -- end switch -- */
    }

/************************** end of grav_calc *********************************/ 




/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   dgrav_form   ***************
*/

/*  INTERFACE:   */

      void
      dgrav_form(
                 real  *dgrav ,
                 ELEMENT_INFO  *info ,
                 int  side ,
                 real  dmass ,
                 real  *midpt ,
                 ELEMENT_DATA  *bel_pt
                )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine dgrav_form
** dgrav_form performs the part of the gravity change calculation resulting from
** the net Eulerian flow of mass in/out of each volume element, the results
** being evaluated at the position of each surface facet that is a buoyancy element.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/



{
 int  nsd , ndof , node1 , node2 , node3 , node4 , n1 , n2 , n3 , n4 ;
 real  xx1 , xx2 , xx3 , xx4 , yy1 , yy2 , yy3 , yy4 , zz1 , zz2 , zz3 , zz4 , 
       xxave , yyave , zzave , dist ;
 real  radial[3] ;


      nsd   = fe_sys.nsd ;
      ndof  = fe_sys.ndof ;
      
          
      if (3 == nsd && (TET == info->type))
         {

  /* check and correct for tets that have been flopped in chirality */
         if(bel_pt->degen == 2)
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
         n1 = *(bel_pt->ien + node1 - 1) - 1 ;  /* look up global node #'s */
         n2 = *(bel_pt->ien + node2 - 1) - 1 ;
         n3 = *(bel_pt->ien + node3 - 1) - 1 ;
         xx1 = *(global.coords + n1*nsd) ;     /* then get the global coords */
         yy1 = *(global.coords + n1*nsd + 1) ;
         zz1 = *(global.coords + n1*nsd + 2) ;
         xx2 = *(global.coords + n2*nsd) ;
         yy2 = *(global.coords + n2*nsd + 1) ;
         zz2 = *(global.coords + n2*nsd + 2) ;
         xx3 = *(global.coords + n3*nsd) ;
         yy3 = *(global.coords + n3*nsd + 1) ;
         zz3 = *(global.coords + n3*nsd + 2) ;

         xxave = (xx1+xx2+xx3)/THREE ;    /* compute the face middle... */
         yyave = (yy1+yy2+yy3)/THREE ;
         zzave = (zz1+zz2+zz3)/THREE ;

           }

      else if (3 == nsd && (TRILIN == info->type))
         {

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
         
         n1 = *(bel_pt->ien + node1 - 1) - 1 ;  /* look up global node #'s */
         n2 = *(bel_pt->ien + node2 - 1) - 1 ;
         n3 = *(bel_pt->ien + node3 - 1) - 1 ;
         n4 = *(bel_pt->ien + node4 - 1) - 1 ;
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

         xxave = (xx1+xx2+xx3+xx4)/FOUR ;    /* compute the face middle... */
         yyave = (yy1+yy2+yy3+yy4)/FOUR ;
         zzave = (zz1+zz2+zz3+zz4)/FOUR ;

           }
           
      radial[0] = xxave - midpt[0] ;  /* field point - source point  */
      radial[1] = yyave - midpt[1] ;
      radial[2] = zzave - midpt[2] ;
      
      dist = sqrt(radial[0]*radial[0] + radial[1]*radial[1] + radial[2]*radial[2]) ;
      
      dgrav[0] += -dmass * radial[0] / (dist*dist*dist) ;
      dgrav[1] += -dmass * radial[1] / (dist*dist*dist) ;
      dgrav[2] += -dmass * radial[2] / (dist*dist*dist) ;
      dgrav[3] += -dmass/dist ;          /* scalar potential contribution */
      dgrav[4] += -THREE*dmass*radial[0]*radial[2] / (dist*dist*dist*dist*dist) ;  /* P_zx gradient */
      dgrav[5] += -THREE*dmass*radial[1]*radial[2] / (dist*dist*dist*dist*dist) ;  /* P_zy gradient */
      dgrav[6] += dmass*(ONE/(dist*dist*dist) - 
                       (THREE*radial[2]*radial[2])/(dist*dist*dist*dist*dist)) ;   /* P_zz gradient */

}

/************************** end of dgrav_form *********************************/ 

/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   dgrav_flux   ***************
*/

/*  INTERFACE:   */

      void
      dgrav_flux(
                 real  *dgrav ,
                 ELEMENT_INFO  *info ,
                 int  rside ,
                 int  dside ,
                 real delta_rho ,
                 real *upvec ,
                 real big_g ,
                 ELEMENT_DATA  *del_pt ,
                 ELEMENT_DATA  *rel_pt ,
                 real *fmid
                )

/*  RETURN VALUE:  -none- */
/*  DESCRIPTION:   */
/*
** Routine dgrav_flux
** dgrav_flux performs the part of the gravity change calculation resulting from
** the Eulerian mass flow imbalance occurring at the interface between
** solid and vaccuum or between contrasting densities, the results
** being evaluated at the position of each surface facet that is a buoyancy element.
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
            reactx , reacty , reactz , centrx , centry , centrz , local_up[3] , dist , dmass ,
            xdmid , ydmid , zdmid , xrmid , yrmid , zrmid , fluxden ;
      real  xplan[4][2]  ; /* 4 nodes by 2 components  */
      real  sh2all[48] , detall[4] ; /* 4 integration pts, 4 nodes , 3 SH's */
      real  xs[2][2] , radial[3] ;
      real  *sh2 , *det , *wt ;
     
      nsd   = fe_sys.nsd ;
      nsdp1 = nsd + 1 ;
      nee   = info->nee ;
      ndof  = fe_sys.ndof ;
      
      dmass = ZERO ;
          
          
      if (3 == nsd && (TET == info->type))
         {

  /* check and correct for tets that have been flopped in chirality */
         if(del_pt->degen == 2)
            {
             if(dside == 3)  dside = 4 ;
             else if(dside == 4)  dside = 3 ;
            }

         switch(dside)
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
              printf("IAM:%d, Bad Side:%d", iam, dside);
              exit(EXIT_SUCCESS);
              break;
         }
         n1 = *(del_pt->ien + node1 - 1) - 1 ;  /* look up global node #'s */
         n2 = *(del_pt->ien + node2 - 1) - 1 ;
         n3 = *(del_pt->ien + node3 - 1) - 1 ;
         xx1 = *(global.coords + n1*nsd) ;     /* then get the global coords */
         yy1 = *(global.coords + n1*nsd + 1) ;
         zz1 = *(global.coords + n1*nsd + 2) ;
         xx2 = *(global.coords + n2*nsd) ;
         yy2 = *(global.coords + n2*nsd + 1) ;
         zz2 = *(global.coords + n2*nsd + 2) ;
         xx3 = *(global.coords + n3*nsd) ;
         yy3 = *(global.coords + n3*nsd + 1) ;
         zz3 = *(global.coords + n3*nsd + 2) ;
         
         /* while we're at it, compute the face center position of the *donor* facet */
         
         xdmid = (xx1+xx2+xx3)/THREE ;
         ydmid = (yy1+yy2+yy3)/THREE ;
         zdmid = (zz1+zz2+zz3)/THREE ;


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
         
         ex = ay*bz - az*by ;
         ey = az*bx - ax*bz ;
         ez = ax*by - ay*bx ;
         elength = sqrt(ex*ex+ey*ey+ez*ez) ;
         ex /= elength ;
         ey /= elength ;
         ez /= elength ;   /* e is a unit vector perp to the axb plane */
         
         if(ex*upvec[0]+ey*upvec[1]+ez*upvec[2] < ZERO)    /* verify outward normal is "up" */
            {
             ex = -ex ;
             ey = -ey ;
             ez = -ez ;
            }              /* note that this revision only correctly uses rectilinear upvec, -not- radial */
         
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
         
         /* ...and calculate mass flux */
         
         fluxden = ex*avex + ey*avey + ez*avez ;
                     
     /* now do 2D Gaussian integration of the flux density */
            
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
               printf("Negative buoyface det: element=%d  ipt=%d\n",del_pt->nel,ipt) ;
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
        
        dmass += w_det * fluxden * delta_rho * big_g ;

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
         switch(dside)
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
         
         n1 = *(del_pt->ien + node1 - 1) - 1 ;  /* look up global node #'s */
         n2 = *(del_pt->ien + node2 - 1) - 1 ;
         n3 = *(del_pt->ien + node3 - 1) - 1 ;
         n4 = *(del_pt->ien + node4 - 1) - 1 ;
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

         /* while we're at it, compute the face center position of the *donor* facet */
         
         xdmid = (xx1+xx2+xx3+xx4)/FOUR ;
         ydmid = (yy1+yy2+yy3+yy4)/FOUR ;
         zdmid = (zz1+zz2+zz3+zz4)/FOUR ;

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
         
         /* ...and calculate mass flux */
         
         if(ex*upvec[0]+ey*upvec[1]+ez*upvec[2] < ZERO)    /* verify outward normal is "up" */
            {
             ex = -ex ;
             ey = -ey ;
             ez = -ez ;
            }              /* note that this revision only correctly uses rectilinear upvec, -not- radial */
         
         fluxden = ex*avex + ey*avey + ez*avez ;
            
     /* now do 2D Gaussian integration of the flux density */
            
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
               printf("Negative buoyface det: element=%d  ipt=%d\n",del_pt->nel,ipt) ;
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
        

        dmass += w_det * fluxden * delta_rho * big_g ;

	   }
            
         
           }
           
/*  calculate the grav/potential changes from dmass and locations  */

/* We now need the center locations of the *receiver* facets... */

      if (3 == nsd && (TET == info->type))
         {

  /* check and correct for tets that have been flopped in chirality */
         if(rel_pt->degen == 2)
            {
             if(rside == 3)  rside = 4 ;
             else if(rside == 4)  rside = 3 ;
            }

         switch(rside)
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
              printf("IAM:%d, Bad Side:%d", iam, rside);
              exit(EXIT_SUCCESS);
              break;
         }
         n1 = *(rel_pt->ien + node1 - 1) - 1 ;  /* look up global node #'s */
         n2 = *(rel_pt->ien + node2 - 1) - 1 ;
         n3 = *(rel_pt->ien + node3 - 1) - 1 ;
         xx1 = *(global.coords + n1*nsd) ;     /* then get the global coords */
         yy1 = *(global.coords + n1*nsd + 1) ;
         zz1 = *(global.coords + n1*nsd + 2) ;
         xx2 = *(global.coords + n2*nsd) ;
         yy2 = *(global.coords + n2*nsd + 1) ;
         zz2 = *(global.coords + n2*nsd + 2) ;
         xx3 = *(global.coords + n3*nsd) ;
         yy3 = *(global.coords + n3*nsd + 1) ;
         zz3 = *(global.coords + n3*nsd + 2) ;

         xrmid = (xx1+xx2+xx3)/THREE ;    /* compute the face middle... */
         yrmid = (yy1+yy2+yy3)/THREE ;
         zrmid = (zz1+zz2+zz3)/THREE ;

		  fmid[0] = xrmid ;    /* make these available to the calling routine... */
		  fmid[1] = yrmid ;
		  fmid[2] = zrmid ;
           }

      else if (3 == nsd && (TRILIN == info->type))
         {

         switch(rside)
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
         
         n1 = *(rel_pt->ien + node1 - 1) - 1 ;  /* look up global node #'s */
         n2 = *(rel_pt->ien + node2 - 1) - 1 ;
         n3 = *(rel_pt->ien + node3 - 1) - 1 ;
         n4 = *(rel_pt->ien + node4 - 1) - 1 ;
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

         xrmid = (xx1+xx2+xx3+xx4)/FOUR ;    /* compute the face middle... */
         yrmid = (yy1+yy2+yy3+yy4)/FOUR ;
         zrmid = (zz1+zz2+zz3+zz4)/FOUR ;

		  fmid[0] = xrmid ;    /* make these available to the calling routine... */
		  fmid[1] = yrmid ;
		  fmid[2] = zrmid ;

           }

      radial[0] = xrmid - xdmid ;  /* field point - source point  */
      radial[1] = yrmid - ydmid ;
      radial[2] = zrmid - zdmid ;
      dist = sqrt(radial[0]*radial[0] + radial[1]*radial[1] + radial[2]*radial[2]) ;
      
      dgrav[0] += -dmass * radial[0] / (dist*dist*dist) ;
      dgrav[1] += -dmass * radial[1] / (dist*dist*dist) ;
      dgrav[2] += -dmass * radial[2] / (dist*dist*dist) ;
      dgrav[3] += -dmass/dist ;          /* scalar potential contribution */
      dgrav[4] += -THREE*dmass*radial[0]*radial[2] / (dist*dist*dist*dist*dist) ;  /* P_zx gradient */
      dgrav[5] += -THREE*dmass*radial[1]*radial[2] / (dist*dist*dist*dist*dist) ;  /* P_zy gradient */
      dgrav[6] += dmass*(ONE/(dist*dist*dist) - 
                       (THREE*radial[2]*radial[2])/(dist*dist*dist*dist*dist)) ;   /* P_zz gradient */


  }

/************************** end of dgrav_flux *********************************/ 


