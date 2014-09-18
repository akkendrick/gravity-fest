/*---------------------------------------------------------------------------*/
/*   BOP   */
/*
   *************   INCLUDE:   mydefs.h   ***************
*/

/*  DESCRIPTION:
     general purpose definitions for GeoFEST finite element program
*/
/*   EOP   */
/*---------------------------------------------------------------------------*/
/*
***                           File mydefs.h                       ***
***                        GeoFEST version 6.0
*** Copyright (c) 2010, California Institute of Technology        ***
*** U.S.Sponsorship under NASA Contract NAS7-1407 is acknowledged ***
***
*** This software is designated for public release under JPL Task ***
*** Order Number NMO710991 and may be publicly released through
*** license with the Open Channel Foundation
***/


/* general purpose definitions */
#define TRUE 1
#define FALSE 0
#define MAXSUB  (pow(2.0,20.0))
#define STRAIN_FAC  20.0
#define STRAIN_PCT  1.0
#define FOREVER    for (;;)
#define FAIL       1
#define SUCCEED    0
#define YES     1
#define NO      0
#define ABS(x)     (((x) < 0) ? -(x) : (x))
#define MAX(x,y)   (((x) < (y)) ? (y) : (x))
#define MIN(x,y)   (((x) < (y)) ? (x) : (y))
#define FSIGNUM(x)   (((x) >= 0.0) ? (1.0) : (-1.0))

#define real double
#define HZ 100

#define SET_TRUE(x)   x=1
#define SET_FALSE(x)  x=0
#define BETWEEN(x,y,z) (((x) <= (y) && (y) <= (z)) ? 1 : 0)
#define COUNTER    int
#define ZERO       0.0
#define PT25       0.25
#define PT5        0.5
#define ONE        1.0
#define TWO        2.0
#define THREE      3.0
#define FOUR       4.0
#define FIVE       5.0
#define SIX        6.0
