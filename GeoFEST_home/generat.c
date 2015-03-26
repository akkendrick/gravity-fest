/*
***                           File generat.c                      ***
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

#define EXTERN extern
#include "generat.h"
#include "utility.h"
void completion(void);

/* GC: global coordinate value of co (0,1,2 => x,y,z) for node nd */
#define GC(co,nd) (global.coords[(co) + (nd)*nsd])

/* XV:  tet basis vector, component co (0,1,2) for tet node en from node 0*/
#define XV(co,en) (GC(co,el_pt->ien[(en)] - 1) - GC(co,el_pt->ien[0] - 1))

static
     int      el_inc[3] ,
              nod_inc[3] ,
              el_end[3] ;
static
     int         current_element  ,
                 node             ,
                 gen_order        ,
                 nodes_gen[3]     ,
                 node_inc[3]      ;
static int    last_node ;
static int      *temp_list , *temp_where , *temp_grp ;
static real     *temp_trac , *temp_amount , *temp_normal ;
static
     real  temp[3][8]   ;

static
     real shape[8] = { ZERO , ZERO , ZERO , ZERO
              , ZERO , ZERO , ZERO , ZERO } ;



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   gen_element   ***************
*/

/*  INTERFACE:   */

    void
    gen_element(
                GROUP   *grp_ptr    /* pointer to current element group */
               )

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

   {
    COUNTER  i , j , k , m , nel ;
    int       type , numsuf , numbuoy , ndof , nsd , nsf , plastic ,
              numface , tally , nsplit, g_nsplit , more  ;
    int       nel_local , mat , numel , numat , nen , nee,
              g_numel , num_local , surf_done , isurf ;
    int       stress_size  , dbar_size , nint ; 
    int       x_index , y_index , x_offset , y_offset , next_el ;
    SPLITNODE *fspl;
    real      coeff , exponent , trac_time ;
    real      rho_g , upvec[3] ;
    ELEMENT_DATA   *el_pt ;
    BUOY_DATA      *buoy_ptr ;
    int       ien_current[8]  ;
    int       ien_temp[8] ;
    char *test;
    char in_string[MAX_STRING_LENGTH];
    char            msg[MAX_STRING_LENGTH] ;

     /*  begin with element constants . . .  */
     
      /*  here we'll read in numel, type, numat, numsuf, numbuoy, nsplit */
      
    while(TRUE)
       {
        test = fgets(in_string,MAX_STRING_LENGTH,elem_file);

        if(test == NULL)  /* abnormal end of file or hard error */
           {attempt = FLOW_CODES; completion();}

        if(strlen(in_string) > 1) /* skip blank line -- strlen==1*/
           {
            if(isdigit(in_string[0]))
               break;

            /* not a digit or a blank line: process flow line, any order. */
            set_eldata_from_string(in_string,grp_ptr);
           }
       }

     
     numel = (grp_ptr->el_info)->numel ;
     g_numel = numel ;
     type = (grp_ptr->el_info)->type ;
     numat = (grp_ptr->el_info)->numat ;
     numsuf = (grp_ptr->el_info)->numsuf ;
     numbuoy = (grp_ptr->el_info)->numbuoy ;
     g_nsplit = (grp_ptr->el_info)->g_nsplit ;
     (grp_ptr->el_info)->failed = 0 ;

     sprintf(msg,"Global number of elements = %d\n\n",(grp_ptr->el_info)->g_numel) ;
     squawk(msg); 
     sprintf(msg,"Global number of split nodes: %d\n", (grp_ptr->el_info)->g_nsplit);
     squawk(msg);


     grp_ptr->el_data =
        (ELEMENT_DATA *) calloc(numel,sizeof(ELEMENT_DATA)) ;
     if (NULL == grp_ptr->el_data) attempt = ELD_MEM ;

     completion() ;

     switch( type )
        {
         case BILIN:
         case TET:
            nen = (grp_ptr->el_info)->nen = 4 ;
            break ;

         case SERENDIP:
         case TRILIN:
            nen = (grp_ptr->el_info)->nen = 8 ;
            break ;
        }
        
     ndof = fe_sys.ndof ;
     nsd = fe_sys.nsd ;
     nee = (grp_ptr->el_info)->nee = nen*ndof ;
     (grp_ptr->el_info)->el_size = nee*nee ;
     
     if( type == BILIN || type == SERENDIP || type == TET || type == TRILIN )
        {
         switch( ndof )
            {
             case 1:
                stress_size = (grp_ptr->el_info)->stress_size = 2 ;
                dbar_size   = (grp_ptr->el_info)->dbar_size   = 3 ;
                break ;
             case 2:
                stress_size = (grp_ptr->el_info)->stress_size = 4 ;
                dbar_size   = (grp_ptr->el_info)->dbar_size   = 10 ;
                break ;
             case 3:
                stress_size = (grp_ptr->el_info)->stress_size = 6 ;
                dbar_size   = (grp_ptr->el_info)->dbar_size   = 21 ;
                break ;
            }
         nsf = ndof ;
        }
     
     /*  next get element material properties . . .  */

     (grp_ptr->el_mat)->mu     =
              (real *) calloc(numat,sizeof(real)) ;
     if (NULL == (grp_ptr->el_mat)->mu) attempt = MAT_MEM ;

     (grp_ptr->el_mat)->lambda =
              (real *) calloc(numat,sizeof(real)) ;
     if (NULL == (grp_ptr->el_mat)->lambda) attempt = MAT_MEM ;

     (grp_ptr->el_mat)->visc   =
            (real *) calloc(2*numat,sizeof(real)) ;
     if (NULL == (grp_ptr->el_mat)->visc) attempt = MAT_MEM ;

     (grp_ptr->el_mat)->bforce =
         (real *) calloc(ndof*numat,sizeof(real)) ;
     if (NULL == (grp_ptr->el_mat)->bforce) attempt = MAT_MEM ;
     
	 if( time_data.gravcalc )
		{
		 (grp_ptr->el_mat)->density =
			 (real *) calloc(numat,sizeof(real)) ;
		 if (NULL == (grp_ptr->el_mat)->density) attempt = MAT_MEM ;
/*  no longer used ***		 
		 (grp_ptr->el_mat)->big_g =
			 (real *) calloc(ndof*numat,sizeof(real)) ;
		 if (NULL == (grp_ptr->el_mat)->big_g) attempt = MAT_MEM ;
		 
		 (grp_ptr->el_mat)->little_g =
			 (real *) calloc(ndof*numat,sizeof(real)) ;
		 if (NULL == (grp_ptr->el_mat)->little_g) attempt = MAT_MEM ;
*** */
		}

     (grp_ptr->el_mat)->plastic=
              (int *) calloc(numat,sizeof(int)) ;
     if (NULL == (grp_ptr->el_mat)->plastic) attempt = MAT_MEM ;
     
     for(i=0 ; i<numat ; i++)
        {
         fscanf(elem_file,"%lf",(grp_ptr->el_mat)->mu+i) ;
         fscanf(elem_file,"%lf",(grp_ptr->el_mat)->lambda+i) ;
         fscanf(elem_file,"%lf%lf",&coeff,&exponent) ;
         *((grp_ptr->el_mat)->visc+2*i) = coeff ;
         *((grp_ptr->el_mat)->visc+2*i+1) = exponent ;

	 if( time_data.gravcalc )
		{
/*  no longer used ***		 
		 fscanf(grav_file,"%lf",(grp_ptr->el_mat)->big_g+i) ;
		 fscanf(grav_file,"%lf",(grp_ptr->el_mat)->little_g+i) ;
*** */
		 fscanf(grav_file,"%lf",(grp_ptr->el_mat)->density+i) ;
		}

         sprintf(msg,"Material with id: %d, %lg, %lg, %lg, %lg\n",
            i,
            grp_ptr->el_mat->mu[i], grp_ptr->el_mat->lambda[i], 
            grp_ptr->el_mat->visc[2*i], grp_ptr->el_mat->visc[2*i+1] );
         squawk(msg);

         /*
         *  Set plastic code: elastic->0, linear visco ->1, nonlinear visco->2
         */
         if(coeff == ZERO || exponent == ZERO)
            *((grp_ptr->el_mat)->plastic+i) = 0 ;
         else if(exponent == ONE)
            *((grp_ptr->el_mat)->plastic+i) = 1 ;
         else
            *((grp_ptr->el_mat)->plastic+i) = 2 ;
                
         /* printf("mat=%d  mu=%g  lambda=%g  coeff=%g  exp=%g\n",
         *          i+1,*((grp_ptr->el_mat)->mu+i),*((grp_ptr->el_mat)->lambda+i),
         *          coeff,exponent) ;
         * 
         *  printf("material #%d plasticity code = %d\n",
         *            i+1,*((grp_ptr->el_mat)->plastic+i) ) ;
         */

         for(j=0;j<ndof;j++)
            {
             fscanf(elem_file,"%lf",(grp_ptr->el_mat)->bforce+ndof*i+j) ;
            }
        } /* end for loop over numat */
    
     for( nel=0 ; nel<numel ; nel++ )
        {
         el_pt = grp_ptr->el_data + nel ;
          
         /*  now allocate memory for the ien and lm arrays... */

         el_pt->ien = (int *) calloc(nen,sizeof(int)) ;
         if (NULL == el_pt->ien) attempt = IEN_MEM ;

         el_pt->lm = (int *) calloc(nee,sizeof(int)) ;
         if (NULL == el_pt->lm) attempt = LM_MEM ;

        }

     completion() ;

     /*  now generate elements, filling up the arrays just created... */

     nel_local = 0 ;  /* keep a tally of the elements local to this proc */
   
     while (next_element( &current_element , &gen_order , elem_file ) !=  0 )
        {
          fscanf(elem_file,"%d",&mat);  
          for ( i = 0 ; i < nen ; ++i) 
             {
              fscanf(elem_file,"%d",&ien_current[i]);  
             }
             
          if ( gen_order != 0 )
             {
              squawk("Element generation error.\n"); 
              attempt = IEN_MEM ;
              completion() ;
             }
                

                 num_local = current_element ;
                 el_pt = grp_ptr->el_data + (num_local - 1) ;
                 load_element(el_pt,grp_ptr->el_info,
                              num_local,ien_current,mat) ;
                 nel_local++ ;


    
        } /* end while next element */

     if(numel != nel_local) 
        {
         attempt = NEL_ERR ;
        }
/*
     sprintf(msg,"size of int = %d , size of real = %d\n",
         (int)sizeof(int),(int)sizeof(real)) ;
     squawk(msg);
*/
     completion() ;

     for( nel=0 ; nel<numel ; nel++ )
        {
         el_pt = grp_ptr->el_data + nel ;
         m = el_pt->mat ;
         plastic = *( (grp_ptr->el_mat) -> plastic + m ) ;
         nint = el_pt->nint ;
          
         /*  next allocate memory for the stresses and visco arrays... */

         if( type == BILIN || type == SERENDIP || type == TET || type == TRILIN)
            {
             el_pt->stress =
                ( real * ) calloc( nint*stress_size , sizeof(real)) ;
             if (NULL == el_pt->stress) attempt = STRESS_MEM  ;
             
             el_pt->dilat =
                ( real * ) calloc( nint , sizeof(real)) ;
             if (NULL == el_pt->dilat) attempt = STRESS_MEM  ;

             el_pt->is_split =
                ( int * ) calloc( 8 , sizeof(int)) ;
             /* there should never be more than 7 split nodes in an element, I hope */  

             if (NULL == el_pt->is_split) attempt = STRESS_MEM  ;
             *(el_pt->is_split) = 0 ; /* initialize just to be sure */

             if( plastic )
                {
                 el_pt->bta =
                 ( real * ) calloc( nint*stress_size , sizeof(real)) ;
                 if (NULL == el_pt->bta) attempt = BTA_MEM  ;
            
                 el_pt->dbar =
                 ( real * ) calloc( nint*dbar_size , sizeof(real)) ;
                 if (NULL == el_pt->dbar) attempt = DBAR_MEM  ;

                }
            }

         /*  if PCG is being used, allocate memory for element stiffness storage */
         if(PCG == fe_sys.solver)
            {
             el_pt->stiff =
                ( real * ) calloc( (grp_ptr->el_info)->el_size , sizeof(real)) ;
             if (NULL == el_pt->stiff) attempt = ELSTIFF_MEM  ;
            }
        } /* end for over numel */
     completion() ;

     /*  now allocate memory and read in surface tractions... */
     if( numsuf != 0 )
        {
         squawk("Processing surface tractions. . .\n") ;
         surf_done = 0 ;
         isurf = 0 ;
         time_data.traction_time[0] = 0.0 ;
         time_data.traction_time[1] = 1.0e+99 ;
         time_data.traction_time[2] = 1.0e+99 ;
         time_data.traction_time[3] = 1.0e+99 ;
         time_data.traction_time[4] = 1.0e+99 ;
         time_data.traction_time[5] = 1.0e+99 ;
         time_data.currtsuf = 0.0 ;  /* current onset time */
         time_data.itsuf = 0 ;  /* current onset index */
                    printf("traction_time[0] = %lg\n",time_data.traction_time[0]) ;
                    printf("traction_time[1] = %lg\n",time_data.traction_time[1]) ;
                    printf("traction_time[2] = %lg\n",time_data.traction_time[2]) ;
                    printf("traction_time[3] = %lg\n",time_data.traction_time[3]) ;
                    printf("traction_time[4] = %lg\n",time_data.traction_time[4]) ;
                    printf("traction_time[5] = %lg\n",time_data.traction_time[5]) ;
      do{
         if(isurf > 0)
            {
             fscanf(surf_file,"%lf",&trac_time) ;
             time_data.traction_time[isurf] = trac_time ; /* store onset time */
             printf("Traction delta  #%d reading in...\n",isurf) ;
            }
         
         temp_list = (int *) calloc(numsuf,2*sizeof(int)) ;
         if (NULL == temp_list) attempt = TRAC_MEM ;
         temp_trac = (real *) calloc(numsuf,nsf*sizeof(real)) ;
         if (NULL == temp_trac) attempt = TRAC_MEM ;
         
         tally = 0 ;
         
         read_surf( &tally , temp_list , temp_trac ,
                    type , numel , ndof , &more ) ;

         if(tally != 0)
            {
             temp_list = (int *) realloc(temp_list , 2*sizeof(int)*tally) ;
             temp_trac = (real *) realloc(temp_trac , nsf*sizeof(real)*tally) ;
            }
         else  /* just a length-1 place holder */
            {
             temp_list = (int *) realloc(temp_list , 2*sizeof(int)) ;
             temp_trac = (real *) realloc(temp_trac , nsf*sizeof(real)) ;
            }
         if (NULL == temp_list) attempt = TRAC_MEM ;
         if (NULL == temp_trac) attempt = TRAC_MEM ;

         (grp_ptr->el_info)->numsuf = tally ; /* numsuf now contains the proc LOCAL number! */
         printf (" local numsuf = %d\n" , tally) ;
         fflush(stdout) ;
        
        switch( isurf )
           {
            case 0:
              (grp_ptr->el_info)->surf_list0 = temp_list ;
              (grp_ptr->el_info)->surf_trac0 = temp_trac ;
              break ;
         
             case 1:
              (grp_ptr->el_info)->surf_list1 = temp_list ;
              (grp_ptr->el_info)->surf_trac1 = temp_trac ;
              break ;
         
            case 2:
              (grp_ptr->el_info)->surf_list2 = temp_list ;
              (grp_ptr->el_info)->surf_trac2 = temp_trac ;
              break ;
         
            case 3:
              (grp_ptr->el_info)->surf_list3 = temp_list ;
              (grp_ptr->el_info)->surf_trac3 = temp_trac ;
              break ;
         
            case 4:
              (grp_ptr->el_info)->surf_list4 = temp_list ;
              (grp_ptr->el_info)->surf_trac4 = temp_trac ;
              break ;
         
            case 5:
              (grp_ptr->el_info)->surf_list5 = temp_list ;
              (grp_ptr->el_info)->surf_trac5 = temp_trac ;
              break ;           
           }
         
         surf_done = 1 ;
         if(more == 1)
            {
             ++isurf ;
             surf_done = 0 ;
            }
        } while(surf_done != 1) ;
        time_data.numtsuf = isurf + 1 ; /* number of traction epochs */
        time_data.itsuf = 0 ;  /* initialize */
        time_data.currtsuf = time_data.traction_time[time_data.itsuf] ;
        
         squawk("Surface tractions all read in . . .\n") ;
        }
     else
        {
         time_data.traction_time[0] = 0.0 ;
         time_data.traction_time[1] = 1.0e+99 ;
         time_data.traction_time[2] = 1.0e+99 ;
         time_data.traction_time[3] = 1.0e+99 ;
         time_data.traction_time[4] = 1.0e+99 ;
         time_data.traction_time[5] = 1.0e+99 ;
        squawk("No surface tractions to process . . .\n") ;
        }

     /*  now allocate memory and read in buoyancy surfaces... */

     if( numbuoy != 0 )
        {
         squawk("Processing buoyancy elements. . .\n") ;
          
         (grp_ptr->el_info)->buoy = (BUOY_DATA *) calloc(numbuoy,sizeof(BUOY_DATA)) ;
         if (NULL == (grp_ptr->el_info)->buoy) attempt = BUOY_MEM ;

         for(i=0 ; i<numbuoy ; i++)
            {
             buoy_ptr = (grp_ptr->el_info)->buoy + i ;
             fscanf(buoy_file,"%d%lf%lf%lf%lf",&numface,&upvec[0],&upvec[1],&upvec[2],&rho_g);
   /*
      If gravity calculation is TRUE then for each buoy surface, we need to read in
       delta_rho, the mass density contrast across the surface, as well as the local
       acceleration of gravity and a flag saying whether to output gravity changes
       on this surface.
   */
             if( time_data.gravcalc )
                {
                 fscanf(grav_file,"%lf%lf%d",
                   &(buoy_ptr->delta_rho),&(buoy_ptr->little_g),&(buoy_ptr->grav_out_flag)) ;
                }
             buoy_ptr->upvec[0] = upvec[0] ;
             buoy_ptr->upvec[1] = upvec[1] ;
             buoy_ptr->upvec[2] = upvec[2] ;
             buoy_ptr->rad_flag = 0 ;
             if(numface < 0)
                 {
                  buoy_ptr->rad_flag = 1 ; /* true - means buoyancy will be radial to a point */
                  numface = -numface ;  /* user specified by giving negative numface */
                  
      /* NOTE by GL 8/6/09 -- I think these radial buoyancy modifications will go through
         transparently to the parallel case, since the radial flag is set in all processors
         and a positive value of numface is used subsequently.  HOWEVER, I think that more
         attention will be required to make this work properly with refinement, since the
         new mesh currently will inherit a numface value that is always positive.  This should
         be an easy fix to make if/when desired.
      */
                 }
             buoy_ptr->numface = numface ;
             buoy_ptr->rho_g = rho_g ;
             /*
             * printf("rho_g=%g\nupvec=(%g,%g,%g)\nnumface=%d\n\n",
             *         rho_g,upvec[0],upvec[1],upvec[2],numface) ;
             */
             buoy_ptr->buoy_list = (int *) calloc(numface,2*sizeof(int)) ;
             if (NULL == buoy_ptr->buoy_list) attempt = BUOYLIST_MEM ;
             tally = 0 ;
             read_buoy( &tally , buoy_ptr , 
                    type , numel , ndof ) ;
             if(tally != 0)
                {
                 buoy_ptr->buoy_list = 
		    (int *) realloc(buoy_ptr->buoy_list , 2*sizeof(int)*tally) ;
                 }
             else  /* just a length-1 place holder */
		 {
		  buoy_ptr->buoy_list =
		     (int *) realloc(buoy_ptr->buoy_list , 2*sizeof(int)) ;
		 }
             if (NULL == buoy_ptr->buoy_list) attempt = BUOYLIST_MEM ;

             buoy_ptr->numface = tally ; /* the proc LOCAL number! */
             printf (" local numface = %d\n" , tally) ;
             fflush(stdout) ;

     /* allocate storage for gravity changes... */
     
             if( time_data.gravcalc && buoy_ptr->grav_out_flag )
                {
             buoy_ptr->dgrav = (real *) calloc(numface,7*sizeof(real)) ;
             if (NULL == buoy_ptr->dgrav) attempt = BUOYLIST_MEM ;
                }

            }
         
         squawk("Buoyancy elements all read in . . .\n") ;
        }
     else
         squawk("No buoyancy elements to process . . .\n") ;
         
   /*  allocate storage for nonlinear visco strain checking  */
      (grp_ptr->el_info)->strain_list = (real *) calloc(numel,sizeof(real)) ;


     if( 0 == g_nsplit )
        {
         sprintf(msg,
            "\nNo fault nodes; generation of element group #%d complete.\n",
             grp_ptr->group_num) ;
         squawk(msg);
         return ;
        }
     else
        {
         sprintf(msg,"Number of split nodes to read from input: %d\n",g_nsplit); 
         squawk(msg);
        }
     /*  Now allocate memory for and read in and process split node faults... 
     * Note we don't know how many local split nodes we will have - 
     *safe choice is to allocate space as though all belong to this node.
     *Note these are all freed by end of this subroutine.
     */

     fspl = (SPLITNODE *) calloc(g_nsplit,sizeof(SPLITNODE));
     if (NULL == fspl) attempt = SPLITNODE_MEM ;
     completion() ;
          
     read_slip( fspl, type , &nsplit, g_nsplit , ndof ) ;

     /*
     * read_slip finds, saves those split nodes that are in this processor's
     * partition.  So this is the first point where we know the number that
     * pertain to this processor. Save in the el_info structure.
     */
     (grp_ptr->el_info)->nsplit = nsplit ;  

     compute_fterms(fspl,grp_ptr,nsplit,g_nsplit) ; 

         printf("Number of split node element terms: %d\n",(grp_ptr->el_info)->nfterm);

     squawk("Split nodes all processed . . .\n") ;

     sprintf(msg,"Generation of element group #%d complete.\n\n",
             grp_ptr->group_num) ;
     squawk(msg);
   }

/************************** end of gen_element *******************************/ 



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
                    )

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

            {
             fscanf(readfile,"%d",node_pointer) ;  
             fscanf(readfile,"%d",gen_pointer) ;

             return(*node_pointer);
            }

/************************** end of next_element *******************************/ 



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
              )

/*  RETURN VALUE:   -none- */
/*  DESCRIPTION:   */
/*
** Routine read_surf
** read_surf reads in surface traction records.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

   {
    int  nel , side , j , n ;
    real  trac[3] ;
    
    switch( type )
       {
        case BILIN:
        case SERENDIP:
        case TET:
        case TRILIN:
          while (next_element( &nel , &side , surf_file ) !=  0 )
             {
              for(j=0;j<ndof;j++)
                fscanf(surf_file,"%lf",trac+j) ;
                
              
              n = nel - 1 ; /* zero-based local */

              *(temp_list+2*(*tally)) = n ;
              *(temp_list+2*(*tally)+1) = side ;
              
              for(j=0;j<ndof;j++)
                *(temp_trac+ndof*(*tally)+j) = trac[j] ;
                
              (*tally)++ ;  /* increment the counter */
             }
             *more = side ;
          break ;

        case QUAKE:
          break ;
       }
   }
/************************** end of read_surf *********************************/ 



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
              )

/*  RETURN VALUE:   -none- */
/*  DESCRIPTION:   */
/*
** Routine read_buoy
** read_buoy reads in buoyancy surface records.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

   {
    int  nel , side , j , n ;
    
    switch( type )
       {
        case BILIN:
        case SERENDIP:
        case TET:
        case TRILIN:
          while (next_element( &nel , &side , buoy_file ) !=  0 )
             {                
              
              n = nel - 1 ; /* zero-based local */

              *(buoy_ptr->buoy_list+2*(*tally)) = n ;
              *(buoy_ptr->buoy_list+2*(*tally)+1) = side ;
              
              (*tally)++ ;  /* increment the counter */
             }
          break ;

        case QUAKE:
          break ;
       }
   }
/************************** end of read_buoy *********************************/ 



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
              )

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

   {
    int   n, i, j, k, itemp, count;
    real ftemp ;
    SPLITNODE xspl; /* temporary space while we examine if in this partition */
    
    switch( type )
       {
        case BILIN:
        case SERENDIP:
        case TET:
        case TRILIN:
          count=0;
          for(n=0 ; n<g_nsplit ; n++)
             {
              /* 
              *  node # , # of groups , ivec[0] , ivec[1] , ivec[2] , 
              *      svec[0] , svec[1] , svec[2] , grp1 , slip1 , grp2 , slip2...
              */
              fscanf(flt_file,"%d%d%lf%lf%lf%lf%lf%lf",
                &xspl.node, &xspl.ngrp, &xspl.ivec[0], &xspl.ivec[1],&xspl.ivec[2],
                &xspl.svec[0],&xspl.svec[1],&xspl.svec[2]);
               for(i=0;i<xspl.ngrp;i++)
                  fscanf(flt_file,"%d%lf",&xspl.grp[i],&xspl.slip[i]) ;

              if(xspl.node >= 0)
                 {
                  fspl[count] = xspl; /* structure copy */
                  count++;
                 }
             }
          break ;

        case QUAKE:
          break ;
       }
    *ptr_nsplit = count;
    
  
   }

/************************** end of read_slip *********************************/ 



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
              )

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

{
 int  bot , mid , top ;
 
 if(nsplit < 1) return ( -1);
 bot = 0 ;
 top = nsplit-1 ;
 mid = top/2 ;
 
 if( !BETWEEN(fspl[bot].node, qnode,fspl[top].node ) )  return( -1 ) ;
 if(qnode == fspl[bot].node)  return( bot ) ;
 if(qnode == fspl[top].node)  return( top ) ;
 if(qnode == fspl[mid].node)  return( mid ) ;
 while((top-bot) > 2)
    {
     if( BETWEEN(fspl[bot].node,qnode,fspl[mid].node))
        {
         top = mid ;
         mid = bot + (top-bot)/2 ;
         if(qnode == fspl[mid].node) return( mid ) ;
        }
     else
        {
         bot = mid ;
         mid = bot + (top-bot)/2 ;
         if(qnode == fspl[mid].node )  return( mid ) ;
        }
    }
 return( -1 ) ;
}
/************************** end of splitfind *********************************/ 



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
              )

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

{
      COUNTER i,j,k, count, tally,nel,ii;

      int nfterm,this_node,multiplicity,tensile,index;

      int numel = (grp_ptr->el_info)->numel ,
          nen = (grp_ptr->el_info)->nen,
          type = (grp_ptr->el_info)->type,
          nsd = fe_sys.nsd ,
          ndof = fe_sys.ndof ;

      real      xctr , yctr , zctr , dot , mult , norm , move ;

      real      slip[8] , ivec[3] , svec[3] , cross[3] ;
      ELEMENT_DATA   *el_pt ;

      if(nsplit > 1){
         qsort((void*)fspl, (size_t) nsplit, sizeof(SPLITNODE), split_local_node_compare);
      }
      
                     
      /* 
      *  Save this in attributes area.
      *  Note these are pre-zeroed  (by calloc), in inphase.c
      *  So the values will be zero, except for those in the nsplit list 
      *  (in fspl i=0:nsplit-1)
      */

/* the block below is not yet finished - GL */
/*
      for(i=0; i < nsplit; i++) 
         {
          if (fspl[i].node < 0) continue; 
          global.splitn[NSPLITATTR*(fspl[i].node-1)]   = fspl[i].ivec[0];
          global.splitn[NSPLITATTR*(fspl[i].node-1)+1] = fspl[i].ivec[1];
          global.splitn[NSPLITATTR*(fspl[i].node-1)+2] = fspl[i].ivec[2];
          global.splitn[NSPLITATTR*(fspl[i].node-1)+3] = fspl[i].svec[0];
          global.splitn[NSPLITATTR*(fspl[i].node-1)+4] = fspl[i].svec[1];
          global.splitn[NSPLITATTR*(fspl[i].node-1)+5] = fspl[i].svec[2];
          global.splitn[NSPLITATTR*(fspl[i].node-1)+6] = fspl[i].slip[0];
          global.splitn[NSPLITATTR*(fspl[i].node-1)+7] = (real) 1;
          global.splitn[NSPLITATTR*(fspl[i].node-1)+8] = (real) fspl[i].grp[0];
         }
*/

      completion() ;
      
      /*
      * Dry run: just count the instances so we can allocate before filling.
      */
      
      tally = 0 ;
      for( nel=0 ; nel<numel ; nel++ )
         {
          el_pt = grp_ptr->el_data + nel ;
          for(i=0 ; i<nen ; i++)
             {
              /* careful not to double count degenerate quads! */
              if(type == BILIN && i==(nen-1) && el_pt->degen == 1) continue ;
              
              this_node = *(el_pt->ien + i) ;

              /* 
              *  j is the zero-based index number of the searched-for split node input entry
              */
              j = splitfind(this_node,fspl,nsplit) ;
              if(j != -1)
                 {
                  multiplicity = fspl[j].ngrp ;
                  tally += multiplicity ;
                 }
             }     /* loop i over element nodes */
         }         /* loop nel over elements */

      completion() ;

      (grp_ptr->el_info)->nfterm = nfterm = tally ;  /* total number of instances:nel, this_node */

      temp_where = (int *) calloc(nfterm,2*sizeof(int)) ;
      temp_amount = (real *) calloc(nfterm,ndof*sizeof(real)) ;
      temp_grp = (int *) calloc(nfterm,sizeof(int)) ;
      temp_normal = (real *) calloc(nfterm,ndof*sizeof(real)) ;

      if (temp_where == NULL || temp_amount == NULL || temp_grp == NULL) attempt = SPLIT_MEM  ;

      tally = 0 ;
      for( nel=0 ; nel<numel ; nel++ )
         {
          el_pt = grp_ptr->el_data + nel ;
          for(i=0 ; i<nen ; i++)
             {
             /* careful not to double count degenerate quads! */
              if(type == BILIN && i==(nen-1) && el_pt->degen == 1) continue ;
              
              this_node = *(el_pt->ien + i) ;

                  /* 
                  *  j is the zero-based index number of the searched-for split node input entry
                  */
                  j = splitfind(this_node,fspl,nsplit) ;
                  if(j != -1)
                     {
                      multiplicity = fspl[j].ngrp ;
                      for(ii=0;ii<multiplicity;ii++)
                         {
                          index = tally + ii ;
                          *(temp_where+2*index) = nel ;         /* zero-based el# */
                          *(temp_where+2*index+1) = this_node ; /* one-based node# */
                          *(temp_grp+index) = fspl[j].grp[ii] ;  /* strand group# */
                          count = (*(el_pt->is_split)) + 1 ;
                          *(el_pt->is_split) = count ;
                          *(el_pt->is_split + count) = index ;
                         }
                      
                      /* 
                      *  Flag and indices set for easy lookup of split nodes later during stress calc 
                      */

                      ivec[0] = fspl[j].ivec[0] ;
                      ivec[1] = fspl[j].ivec[1] ; /* ivec is the direction perp to slip */
                      ivec[2] = fspl[j].ivec[2] ;
                      svec[0] = fspl[j].svec[0] ;
                      svec[1] = fspl[j].svec[1] ; /* svec is the slip vector direction */
                      svec[2] = fspl[j].svec[2] ; ;
                      for(ii=0;ii<multiplicity;ii++)
                         {
                          slip[ii] = fspl[j].slip[ii]  ; /* slip is the full slip amplitude */
                         }
                     
   /*  7.2.09 GL - If ivec is given as zero, then the sense of fault motion is assumed tensile  */
   
                     
                      norm = (float)sqrt(ivec[0]*ivec[0] + ivec[1]*ivec[1] + ivec[2]*ivec[2]) ;
                      if(norm == ZERO) tensile = 1 ;
                      else             tensile = 0 ;

		      /* Renormalizing 9.25.2007 - jwp */
		      
		              if(!tensile)
		                {
                         ivec[0] = ivec[0] / norm ;
                         ivec[1] = ivec[1] / norm ;
                         ivec[2] = ivec[2] / norm ;
                        }

		      /*
                      *if((float)fabs(norm - ONE) > 1.0e-03)
                      *   {
                      *    squawk("ivec not normalized!\n");
                      *   printf("node #%d , splitnode #%d , instance #%d\n",
                      *          this_node,j+1,tally) ;
                      *   printf("ivec = [ %g %g %g ]\n",ivec[0],ivec[1],ivec[2]) ;
                      *   printf("|ivec| = %g\n", norm );
                      *   attempt = IVEC_NORM ;
                      *  }
		      */

                      norm = sqrt(svec[0]*svec[0] + svec[1]*svec[1] + svec[2]*svec[2]) ;

		      /* Renormalizing 9.25.2007 - jwp */
                      svec[0] = svec[0] / norm ;
                      svec[1] = svec[1] / norm ;
                      svec[2] = svec[2] / norm ;

		      /*
                      * if(fabs(norm - ONE) > 1.0e-03)
                      *  {
                      *   squawk("svec not normalized!\n");
                      *   printf("node #%d , splitnode #%d , instance #%d\n",
                      *          this_node,j+1,tally) ;
                      *   printf("svec = [ %g %g %g ]\n",svec[0],svec[1],svec[2]) ;
                      *   printf("|svec| = %g\n", norm );
                      *   attempt = SVEC_NORM ;
                      *  }
                      * if(fabs(ivec[0]*svec[0]+ivec[1]*svec[1]+ivec[2]*svec[2]) > 1.0e-03)
                      *  {
                      *   squawk("ivec and svec not perpendicular!\n");
                      *   attempt = ISVEC_PERP ;
                      *  }
		      */
                     
		           if(!tensile)
		             {
                      cross[0] = svec[1]*ivec[2] - svec[2]*ivec[1] ;
                      cross[1] = svec[2]*ivec[0] - svec[0]*ivec[2] ;
                      cross[2] = svec[0]*ivec[1] - svec[1]*ivec[0] ;
                     }
		           else
		             {
                      cross[0] = svec[0] ;
                      cross[1] = svec[1] ;
                      cross[2] = svec[2] ;
                     }
                      xctr = ZERO ;
                      yctr = ZERO ;
                      zctr = ZERO ;
                      for(k=0 ; k<nen ; k++)
                         {
                          xctr += *(global.coords + nsd * ( *(el_pt->ien + k) - 1 ) ) ;
                          yctr += *(global.coords + nsd * ( *(el_pt->ien + k) - 1 ) + 1 ) ;
                          if(nsd == 3)
                          zctr += *(global.coords + nsd * ( *(el_pt->ien + k) - 1 ) + 2 ) ;
                         }
                      xctr = xctr/(real)nen - *(global.coords + nsd * ( this_node - 1 ) ) ;
                      yctr = yctr/(real)nen - *(global.coords + nsd * ( this_node - 1 ) + 1 ) ;
                      if(nsd == 3)
                      zctr = zctr/(real)nen - *(global.coords + nsd * ( this_node - 1 ) + 2 ) ;
                      
                      dot = cross[0]*xctr + cross[1]*yctr + cross[2]*zctr ;

                      if( dot > ZERO )
                         mult = PT5 ;
                      else
                         mult = -PT5 ;
                         
                      norm = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]) ;
                      cross[0] = cross[0] / norm ;
                      cross[1] = cross[1] / norm ;
                      cross[2] = cross[2] / norm ;
                       
                         
                      for(ii=0;ii<multiplicity;ii++)
                         {
                          index = tally + ii ;
                          move = mult * slip[ii] ;
                          
                      switch( ndof )
                        {
                          case 1:
                              *(temp_amount+ndof*index) = svec[2]*move ;
                              break ;
                          case 2:
                              *(temp_amount+ndof*index) = svec[0]*move ;
                              *(temp_amount+ndof*index+1) = svec[1]*move ;
                              break ;
                          case 3:
                              *(temp_amount+ndof*index) = svec[0]*move ;
                              *(temp_amount+ndof*index+1) = svec[1]*move ;
                              *(temp_amount+ndof*index+2) = svec[2]*move ;
                              
                              *(temp_normal+ndof*index) = cross[0] ;
                              *(temp_normal+ndof*index+1) = cross[1] ;
                              *(temp_normal+ndof*index+2) = cross[2] ;
                              break ;
                         }
                        }
                      tally += multiplicity ;

                     }
                  
             }     /* loop i over element nodes */
         }         /* loop nel over elements */

      completion() ;
         
      (grp_ptr->el_info)->slip_list = temp_where ;
      (grp_ptr->el_info)->slip_val = temp_amount ;
      (grp_ptr->el_info)->slip_grp = temp_grp ;
      (grp_ptr->el_info)->slip_normal = temp_normal ;

      free( fspl ) ;
}
/************************** end of compute_fterms *********************************/ 



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
              )

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

    {
     COUNTER i , j , k , index , ieq , ndof , nsd ;
     int   node_number ;
     real s;
     int itmp;
     
     ndof = fe_sys.ndof ; 
     nsd = fe_sys.nsd ; 
     el_pt->nel =  nel_glob ;

/* load the ien array with one-based node numbers ...  */

     for(i=0 ; i<info->nen ; i++)
        {
         if( (*(el_pt->ien + i) = ien[i]) )
	      {
           continue ;
	      }
         else
	      {
           attempt = IEN_ERR*1000 + nel_glob ;
	      }
             /* Error if glob_to_local returns zero local node # */
        }
/* jwp: if tet, check triple product and reorder if necessary. TP should
be > 0 if tet is right-handed. Check for determinant >0 (later) at all
wt points .  
*/

     el_pt->degen = 0 ;

     if(info->nen == 4 && nsd == 3) {
        s = XV(0,1)*XV(1,2)*XV(2,3) - XV(2,1)*XV(1,2)*XV(0,3) +
            XV(1,1)*XV(2,2)*XV(0,3) - XV(0,1)*XV(2,2)*XV(1,3) +
            XV(2,1)*XV(0,2)*XV(1,3) - XV(1,1)*XV(0,2)*XV(2,3);

        if(s < 0.0 && mat > 0)  /* negative: switch if a real tet */
        {
           itmp          = el_pt->ien[3];
           el_pt->ien[3] = el_pt->ien[2];
           el_pt->ien[2] = itmp;
           el_pt->degen = 2 ;    /* flag this as a tet with flopped chirality */
        }
     }

     switch( info->type )
        {
         case BILIN:
             el_pt->nint = 4 ;
             if( ien[2] == ien[3] )
                el_pt->degen = 1 ;
             break ;

         case TET:
             el_pt->nint = 4 ;
             if(mat < 0)
                {
                 mat = -mat ;
                 el_pt->degen = 1 ;  /* this is actually a "duo" element */
                }
             break ;

         case TRILIN:
             el_pt->nint = 6 ;
             if( ien[2] == ien[3] && ien[6] == ien[7] )
                el_pt->degen = 1 ; /* triangular prism */
             if( ien[2] == ien[3] && ien[4] == ien[5] && ien[4] == ien[6] && ien[4] == ien[7] )
                el_pt->degen = 2 ; /* tetrahedron */
             if(el_pt->degen != 0)  el_pt->nint = 8 ;
             break ;

         case SERENDIP:
             el_pt->nint = 9 ;
             if( ien[2] == ien[3] && ien[2] == ien[6] )
                el_pt->degen = 1 ;
             break ;
        }
     el_pt->mat = mat-1 ;        /* mat is stored in zero-based form  */

/* and load the lm array with zero-based equation numbers  */

     el_pt->is_bc = 0 ;
     index = 0 ;
     for( j = 0 ; j < info->nen ; j++)
        {
         node_number = *(el_pt->ien + j) ;
         for( k = 0 ; k < ndof ; k++)
             {
              ieq = *(global.id_pointer + ndof*(node_number-1) + k) ;
              *(el_pt->lm + index) = ieq ;
              if( ieq == (-1) )
                 el_pt->is_bc = 1 ;
              index++ ;
             }
        }
    }

/************************** end of load_element ******************************/ 



/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   ROUTINE:   gen_map   ***************
*/

/*  INTERFACE:   */

    void
    gen_map()

/*  RETURN VALUE:   -none- */
/*  DESCRIPTION:   */
/*
** Routine gen_map
** gen_map fills the nodal activity map array.
**/
/*   EOP   */
/*---------------------------------------------------------------------------*/

   {
    int		ndof = fe_sys.ndof;
    int		current_node, node_indx;
    int	    next_node;
    int	gen_inc ;
    int	end_node ;
    COUNTER	k , i ;
    int	bctype[3];
    
     
    /* the map array contains space for 
    up to 3 d.o.f. type designators */

    squawk("Start reading boundary info . . .\n") ;
    
    /* by default, make all entries free unless overidden */
    for(i=0 ; i<loc_sys.numnp ; i++)
       {
        for(k=0 ; k<ndof ; k++)
          *(global.eqtype+ndof*i+k) = 1 ;
       }
     
  /* now read in specified designators */
  
    while (next_element( &current_node , &gen_order , bcc_file ) !=  0 )
      {
        for(k=0 ; k<ndof ; k++)
           fscanf(bcc_file,"%d",&bctype[k]) ;

            node_indx = current_node - 1 ; /* make it zero-based */
            for(k=0 ; k<ndof ; k++)
                *(global.eqtype+ndof*(node_indx)+k) = bctype[k] ;


    if(0 != gen_order)  /* non-functional old code! */
       {
              squawk("BCC generation error.\n"); 
              attempt = IEN_MEM ;
              completion() ;
       }
    
      }  /* end while */
     
     completion() ; 

    squawk("Boundary code info input complete.\n") ;
     
    }

/************************** end of gen_map ***********************************/ 



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
            )

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

     {

      COUNTER i , j  ;
      int    l_node , n ;
      real   value ;

      n = -1 ;
      while ( n != 0 )
            {
             fscanf(the_file,"%d",&n) ;
             if(0 != n)
               {
		        l_node = n ;
                for(i = 0 ; i < dim ; ++i)
                   {
                    fscanf(the_file,"%lf",&value) ;
                      *(array + (l_node-1)*dim + i) = value ;
                   } 
               }
            }
    }

/************************** end of gen_real **********************************/ 



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
                )

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

      {
       real sum = ZERO             ,
             *shape_end = shape + n ;

       while (shape < shape_end) sum += (*shape++)*(*variable++) ;

       return(sum) ;
      }

/************************** end of dot_sh ************************************/ 




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

    if(0  == strcmp(simtask, "NUMEL"))
       {
        ivalue = atoi(value);
        (grp_ptr->el_info)->numel = ivalue ;
        (grp_ptr->el_info)->g_numel = ivalue ;
       }
    else if(0  == strcmp(simtask, "TYPE"))
       {
        ivalue = atoi(value);
        (grp_ptr->el_info)->type = ivalue ;
       }
    else if(0  == strcmp(simtask, "NUMAT"))
       {
        ivalue = atoi(value);
        (grp_ptr->el_info)->numat = ivalue ;
       }
    else if(0  == strcmp(simtask, "NUMSUF"))
       {
        ivalue = atoi(value);
        (grp_ptr->el_info)->numsuf = ivalue ;
       }
    else if(0  == strcmp(simtask, "NUMBUOY"))
       {
        ivalue = atoi(value);
        (grp_ptr->el_info)->numbuoy = ivalue ;
       }
    else if(0  == strcmp(simtask, "NSPLIT"))
       {
        ivalue = atoi(value);
        (grp_ptr->el_info)->g_nsplit = ivalue ;
       }
    else
       {
        squawk("Element input parameter not understood!\n");
        attempt = FLOW_CODES; completion();
       }

    /*
    * Report that this value is set. 
    */
        sprintf(msg,"Set: %s -> %d\n",simtask,ivalue);

    squawk(msg);
}
/************************** end of set_param_from_string ********************************/ 
