#include <Rdefines.h>
#include <Rinternals.h>
#include <R.h>

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>   
#include <algorithm>
#include <numeric>
#include "matalloc.h"
#include <omp.h>

using namespace std;

int             outofbounds_int;
bool            outofbounds_bool;
double          outofbounds_double;

// 160 by 98  works fine
#define SPP1CAPACITY     5400    // max number of first choke sp (larger than MAXNOINC * MAXHORIZON *NOSIZES= 20 * 5* 52 = 5200) 
#define SPP2CAPACITY     5400    // max number of second choke sp 
#define NOSPEC              5    // number of species used in the analysis
#define NOSIZES             5    // number of size classes used in the analysis
#define MAXNOINC           20    // max number of increments
#define MAXHORIZON         52    // number of seasons to fish (= time)

typedef unsigned UFINT;

typedef float  (*FTYPE)[SPP2CAPACITY]; //matalloc(sizeof(float), (void *)0, 2, kSpp1Cap, kSpp2Cap); /* make pointer called FTYPE for FF0Star and FF1 arrays (that only has states and stores V for best choices) */
typedef float  (*FCTYPE)[SPP2CAPACITY][3]; // (FCTYPE) matalloc(sizeof(float), (void*) 0,3, kSpp1Capacity, kSpp2Capacity,kNPatch) ;  /* make pointer called FCTYPE for FF0 array (that has states and patch, storing V for all patches/choices)_  */
typedef float  (*ITYPE)[SPP1CAPACITY][SPP2CAPACITY][3]; //  matalloc(sizeof(float), (void *)0, 4, MAXHORIZON, kSpp1Capacity,kSpp2Capacity,kNPatch);  /* make pointer called ITYPE optimal choice array*/
typedef double (*PTYPE)[MAXHORIZON][NOSPEC][NOSIZES][MAXNOINC]; /* make pointer called PTYPE for information on statistical distribution prey*/
typedef double (*ATYPE)[MAXHORIZON][NOSPEC][(NOSIZES*MAXNOINC) - 1];  /* make pointer called ATYPE for information on statistical distribution prey aggregated over sizes*/
typedef float  (*PITYPE)[NOSPEC][NOSIZES];                         // matalloc(sizeof(float), (void *)0, 3, MAXHORIZON, NOSPEC, NOSIZES);; /* make pointer called PITYPE for price information*/


/* random number generator*/
float ranf ()
{
  return ((float)rand()/RAND_MAX); //genereert random getal tussen 0 en 1
}

/* initialise utility function */
float utility (int aLndSpp1,  int aLndSpp1Quota, float aSpp1LndQuotaFine, 
	       int aLndSpp2,  int aLndSpp2Quota, float aSpp2LndQuotaFine)
{                                                              
  float Fine1 = 0;
  float Fine2 = 0;
  if ( aLndSpp1 >= aLndSpp1Quota)      Fine1  = (aLndSpp1 - aLndSpp1Quota)     * aSpp1LndQuotaFine;
  if ( aLndSpp2 >= aLndSpp2Quota)      Fine2  = (aLndSpp2 - aLndSpp2Quota)     * aSpp2LndQuotaFine;
  return ( - (Fine1 + Fine2));
}

/* Calc economic return per choice to include in  */  /* only the landings take part on the short gains  */
float shortTermGains (int aTime, int aNoInc, PTYPE aLndParms, int aPatch, PITYPE aSppPrice)
{      
 float mean0 = 0;
 for (int inc = 0; inc < aNoInc; inc++){
       for (int s = 0; s < NOSPEC; s++){
             for (int si = 0; si < NOSIZES; si++){
                   mean0  += aLndParms[aPatch][aTime][s][si][inc] * inc * aSppPrice[aTime][s][si];     /* NA; only LndParms */ 
        }
     }
  }
  return (mean0) ;
}

float shortTermCosts (int aTime, int aPatch, float aPriceEffort, int anEffortArray[][MAXHORIZON])
{      
 return ((anEffortArray[aPatch][aTime] * aPriceEffort) ) ;
 //aPrice Effort is (fuel use* fuel price) + gear maintenance
}

float FFF (int aLndSpp1,int aLndSpp2, int aNoInc, ATYPE aLndParmsAgg, int aTime, int aNPatch, double aShortTermEcon[], FCTYPE anFF0, FTYPE anFF0Star, FTYPE anFF1, ITYPE aProbChoice, int verbose)
{      
  if ((aLndSpp1 + aNoInc) >  SPP1CAPACITY) Rprintf("running out of FF1 capacity in backward calcs");
  
  double rhs;
  double vMax = -2E120;
  int theLndSpp1Val, theDisSpp1Val,theLndSpp2Val,theDisSpp2Val;
  
  for (int i = 0; i < aNPatch; i++){
    rhs = 0;
    int incL1, incL2;
    float xL1,  xL2; /* L for landings  */
    
    for ( incL1 = 0; incL1 < aNoInc; incL1++) {
      xL1= aLndParmsAgg[i][aTime][0][incL1];
      theLndSpp1Val = aLndSpp1 + incL1;
      for ( incL2 = 0; incL2 <  aNoInc; incL2++) {
	       xL2= aLndParmsAgg[i][aTime][1][incL2];
	       theLndSpp2Val = aLndSpp2 + incL2;
	       rhs += xL1 * xL2 * anFF1[theLndSpp1Val][theLndSpp2Val];  
      }                                  
    }                            
    
    rhs = rhs + aShortTermEcon[i] ; /* the mean should come from the economic funtion */ 
    anFF0[aLndSpp1][aLndSpp2][i] = rhs;
      
    if (rhs > vMax)  {
      vMax = rhs;
      anFF0Star[aLndSpp1][aLndSpp2] = vMax; 	/* Vmax is defined as a very small number */ 
    }						/* so the best choice will be higher than this Vmax */ 
  }						/* rhs is stored in Istar; best choice */ 
  return vMax;
}

/* initialise Simulation experiment (SEXP= S-expressie) */

SEXP SimulateF ( int aSimNumber, int aHorizon, ITYPE aProbChoice, int aNPatch, int aNoInc, vector< float > aSizeSppInc, PTYPE aLndParms, PTYPE aDisParms, int anEffortArray[][MAXHORIZON], int verbose){
  SEXP ReturnObject, SimDims2, SimDims3, choice, spp1Rand, spp2Rand, spp3Rand, spp4Rand,spp5Rand, 
    spp1Landings, spp2Landings, spp3Landings, spp4Landings, spp5Landings,spp1LndHold, spp2LndHold,
    spp3LndHold, spp4LndHold, spp5LndHold,  spp1Discards, spp2Discards, spp3Discards, spp4Discards, 
    spp5Discards,spp1DisHold, spp2DisHold, spp3DisHold, spp4DisHold, spp5DisHold,anEffort;
  
  PROTECT(ReturnObject      = NEW_OBJECT(MAKE_CLASS("Sim")));
  PROTECT(SimDims2          = allocVector(INTSXP,3));
  INTEGER(SimDims2)[0]      = 1; 
  INTEGER(SimDims2)[1]      = aSimNumber; 
  INTEGER(SimDims2)[2]      = aHorizon; 
  PROTECT(SimDims3          = allocVector(INTSXP,3));
  INTEGER(SimDims3)[0]      = NOSIZES; //number of size classes
  INTEGER(SimDims3)[1]      = aSimNumber; 
  INTEGER(SimDims3)[2]      = aHorizon; 

  PROTECT(choice         = allocArray(INTSXP,SimDims2));
  PROTECT(spp1Rand       = allocArray(REALSXP,SimDims3));
  PROTECT(spp2Rand       = allocArray(REALSXP,SimDims3));
  PROTECT(spp3Rand       = allocArray(REALSXP,SimDims3));
  PROTECT(spp4Rand       = allocArray(REALSXP,SimDims3));
  PROTECT(spp5Rand       = allocArray(REALSXP,SimDims3));
  PROTECT(spp1Landings   = allocArray(REALSXP,SimDims3)); //landings by time step
  PROTECT(spp2Landings   = allocArray(REALSXP,SimDims3));
  PROTECT(spp3Landings   = allocArray(REALSXP,SimDims3));
  PROTECT(spp4Landings   = allocArray(REALSXP,SimDims3));
  PROTECT(spp5Landings   = allocArray(REALSXP,SimDims3));
  PROTECT(spp1LndHold    = allocArray(REALSXP,SimDims3)); //cumulative sum of landings
  PROTECT(spp2LndHold    = allocArray(REALSXP,SimDims3));
  PROTECT(spp3LndHold    = allocArray(REALSXP,SimDims3));
  PROTECT(spp4LndHold    = allocArray(REALSXP,SimDims3));
  PROTECT(spp5LndHold    = allocArray(REALSXP,SimDims3));
  PROTECT(spp1Discards   = allocArray(REALSXP,SimDims3)); //discards by time step
  PROTECT(spp2Discards   = allocArray(REALSXP,SimDims3));
  PROTECT(spp3Discards   = allocArray(REALSXP,SimDims3));
  PROTECT(spp4Discards   = allocArray(REALSXP,SimDims3));
  PROTECT(spp5Discards   = allocArray(REALSXP,SimDims3));
  PROTECT(spp1DisHold    = allocArray(REALSXP,SimDims3)); //cumulative sum of discards
  PROTECT(spp2DisHold    = allocArray(REALSXP,SimDims3));
  PROTECT(spp3DisHold    = allocArray(REALSXP,SimDims3));
  PROTECT(spp4DisHold    = allocArray(REALSXP,SimDims3));
  PROTECT(spp5DisHold    = allocArray(REALSXP,SimDims3));
  PROTECT(anEffort       = allocArray(INTSXP,SimDims2));
  
  for (int s= 0; s < aSimNumber; s++)	 { /* go through the x simulations */   
     
    if (verbose==1 && s < 20) Rprintf ("ves %d ", s); R_FlushConsole();

    float Q;
    vector <int> aSpp1LndHold  (NOSIZES,0);
    vector <int> aSpp2LndHold  (NOSIZES,0);
    vector <int> aSpp3LndHold  (NOSIZES,0); 
    vector <int> aSpp4LndHold  (NOSIZES,0);
    vector <int> aSpp5LndHold  (NOSIZES,0); 
    vector <int> aSpp1DisHold  (NOSIZES,0);
    vector <int> aSpp2DisHold  (NOSIZES,0);
    vector <int> aSpp3DisHold  (NOSIZES,0); 
    vector <int> aSpp4DisHold  (NOSIZES,0);
    vector <int> aSpp5DisHold  (NOSIZES,0); 
    
    int Effort = 0;
    int aChoice;
    
    for (int t = 0; t < aHorizon; t++) { 
       
      // get patch from the array
      float Cprobl  = 0;
      float Cprobup = 0;

      aChoice = 0;  //default if it is not after first stick   
      Q = ranf();	/* set starting choice random number 0-1 */

      int totSpp1Hold =0;
      int totSpp2Hold =0;
      
      for (int si = 0; si < NOSIZES; si++) {
	totSpp1Hold += aSpp1LndHold[si];
        totSpp2Hold += aSpp2LndHold[si];
      }
      
      Cprobl = aProbChoice[t][totSpp1Hold][totSpp2Hold][0];

      for (int stick = 1; stick < aNPatch; stick++) {
	Cprobup = aProbChoice[t][totSpp1Hold][totSpp2Hold][stick] + Cprobl;

	if ((Q > Cprobl) && (Q <= Cprobup)){
	  aChoice = stick;
	}
	Cprobl = Cprobup;
      }

      if (verbose==1 && s < 20 ){
	Rprintf (" time %d ", t); 
	Rprintf ("spp1,2  %d,", totSpp1Hold); 	
	Rprintf ("%d ", totSpp2Hold); 	
	
	for (int ccc = 0; ccc < aNPatch; ccc++) {
	  Rprintf(" probs %f,",  aProbChoice[t][totSpp1Hold][totSpp2Hold][ccc]);
	}
	Rprintf (" choice %d \n", aChoice); R_FlushConsole();
      }
      
      INTEGER(choice)[s+ t*aSimNumber] = aChoice + 1;  /* +1 because choice in c++ starts at 0 */
      INTEGER(anEffort)[s+ t*aSimNumber]= anEffortArray[aChoice][t];
      Effort = Effort + INTEGER(anEffort)[s+ t*aSimNumber]; /*calculate fueluse*/
      
      /********************************************************************************/
      /* GET CONSEQUENCES OF CHOICE                                                   */
      /********************************************************************************/
      float Lprobl  = 0;
      float Lprobup = 0;
      float Dprobl  = 0;
      float Dprobup = 0;
      
      /* calculate spp1[1] landings and discards in the patch */        
      for (int si = 0; si < NOSIZES; si++){
	Q = ranf ();	/* set starting ITQ for first sp = random number 0-1 */
	REAL(spp1Rand)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = Q;
	Lprobl  = 0;
	Lprobup = aLndParms[aChoice][t][0][si][0];
	Dprobl  = 0;
	Dprobup = aDisParms[aChoice][t][0][si][0];
	if (Q <= Lprobup) REAL(spp1Landings)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = 0;
	if (Q <= Dprobup) REAL(spp1Discards)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = 0;
	Lprobl = Lprobup;
               Dprobl = Dprobup;
	       
               for (int stick = 1; stick < aNoInc; stick++) {
                   Lprobup = aLndParms[aChoice][t][0][si][stick] + Lprobl;
                   Dprobup = aDisParms[aChoice][t][0][si][stick] + Dprobl;
                   if ((Q > Lprobl) && (Q <= Lprobup)){
                      aSpp1LndHold[si] = aSpp1LndHold[si] + stick;
                      REAL(spp1Landings)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = stick * aSizeSppInc[0];
                   }
                   Lprobl = Lprobup;
                   if ((Q > Dprobl) && (Q <= Dprobup)){
                      aSpp1DisHold[si] = aSpp1DisHold[si] + stick;
                      REAL(spp1Discards)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = stick * aSizeSppInc[0];
                   }
                   Dprobl = Dprobup;
               }
          }
      
          /* calculate spp2  catch in the patch */
          for (int si = 0; si < NOSIZES; si++){
              Q = ranf ();	/* set starting ITQ for first sp = random number 0-1 */
              REAL(spp2Rand)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = Q;
              Lprobl  = 0;
              Lprobup = aLndParms[aChoice][t][1][si][0];
              Dprobl  = 0;
              Dprobup = aDisParms[aChoice][t][1][si][0];
              if (Q <= Lprobup) REAL(spp2Landings)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = 0;
              if (Q <= Dprobup) REAL(spp2Discards)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = 0;
              Lprobl = Lprobup;
              Dprobl = Dprobup;
            
              for (int stick = 1; stick < aNoInc; stick++) {
                  Lprobup = aLndParms[aChoice][t][1][si][stick] + Lprobl;
                  Dprobup = aDisParms[aChoice][t][1][si][stick] + Dprobl;
                  if ((Q > Lprobl) && (Q <= Lprobup)){
                     aSpp2LndHold[si] = aSpp2LndHold[si] + stick;
                     REAL(spp2Landings)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = stick * aSizeSppInc[1];
                  }
                  Lprobl = Lprobup;
                  if ((Q > Dprobl) && (Q <= Dprobup)){
                     aSpp2DisHold[si] = aSpp2DisHold[si] + stick;
                     REAL(spp2Discards)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = stick * aSizeSppInc[1];
                  }
                  Dprobl = Dprobup;
              }
          }
        
          /* calculate spp3  catch in the patch */
          for (int si = 0; si < NOSIZES; si++){
              Q = ranf ();	/* set starting ITQ for first sp = random number 0-1 */
              REAL(spp3Rand)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = Q;
              Lprobl  = 0;
              Lprobup = aLndParms[aChoice][t][2][si][0];
              Dprobl  = 0;
              Dprobup = aDisParms[aChoice][t][2][si][0];
              if (Q <= Lprobup) REAL(spp3Landings)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = 0;
              if (Q <= Dprobup) REAL(spp3Discards)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = 0;
              Lprobl = Lprobup;
              Dprobl = Dprobup;
            
              for (int stick = 1; stick < aNoInc; stick++) {
                  Lprobup = aLndParms[aChoice][t][2][si][stick] + Lprobl;
                  Dprobup = aDisParms[aChoice][t][2][si][stick] + Dprobl;
                  if ((Q > Lprobl) && (Q <= Lprobup)){
                     aSpp3LndHold[si] = aSpp3LndHold[si] + stick;
                     REAL(spp3Landings)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = stick * aSizeSppInc[2];
                  }
                  Lprobl = Lprobup;
                  if ((Q > Dprobl) && (Q <= Dprobup)){
                     aSpp3DisHold[si] = aSpp3DisHold[si] + stick;
                     REAL(spp3Discards)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = stick * aSizeSppInc[2];
                  }
                  Dprobl = Dprobup;
              }
          }
        
          /* calculate spp4  catch in the patch */
          for (int si = 0; si < NOSIZES; si++){
              Q = ranf ();	/* set starting ITQ for first sp = random number 0-1 */
              REAL(spp4Rand)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = Q;
              Lprobl  = 0;
              Lprobup = aLndParms[aChoice][t][3][si][0];
              Dprobl  = 0;
              Dprobup = aDisParms[aChoice][t][3][si][0];
              if (Q <= Lprobup) REAL(spp4Landings)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = 0;
              if (Q <= Dprobup) REAL(spp4Discards)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = 0;
              Lprobl = Lprobup;
              Dprobl = Dprobup;
         
              for (int stick = 1; stick < aNoInc; stick++) {
                  Lprobup = aLndParms[aChoice][t][3][si][stick] + Lprobl;
                  Dprobup = aDisParms[aChoice][t][3][si][stick] + Dprobl;
                  if ((Q > Lprobl) && (Q <= Lprobup)){
                     aSpp4LndHold[si] = aSpp4LndHold[si] + stick;
                     REAL(spp4Landings)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = stick * aSizeSppInc[3];
                  }
                  Lprobl = Lprobup;
                  if ((Q > Dprobl) && (Q <= Dprobup)){
                     aSpp4DisHold[si] = aSpp4DisHold[si] + stick;
                     REAL(spp4Discards)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = stick * aSizeSppInc[3];
                  }
                  Dprobl = Dprobup;
              }
          }
        
          /* calculate spp5  catch in the patch */
          for (int si = 0; si < NOSIZES; si++){
              Q = ranf ();	/* set starting ITQ for first sp = random number 0-1 */
              REAL(spp5Rand)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = Q;
              Lprobl  = 0;
              Lprobup = aLndParms[aChoice][t][4][si][0];
              Dprobl  = 0;
              Dprobup = aDisParms[aChoice][t][4][si][0];
              if (Q <= Lprobup) REAL(spp5Landings)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = 0;
              if (Q <= Dprobup) REAL(spp5Discards)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = 0;
              Lprobl = Lprobup;
              Dprobl = Dprobup;
            
              for (int stick = 1; stick < aNoInc; stick++) {
                  Lprobup = aLndParms[aChoice][t][4][si][stick] + Lprobl;
                  Dprobup = aDisParms[aChoice][t][4][si][stick] + Dprobl;
                  if ((Q > Lprobl) && (Q <= Lprobup)){
                     aSpp5LndHold[si] = aSpp5LndHold[si] + stick;
                     REAL(spp5Landings)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = stick * aSizeSppInc[4];
                  }
                  Lprobl = Lprobup;
                  if ((Q > Dprobl) && (Q <= Dprobup)){
                     aSpp5DisHold[si] = aSpp5DisHold[si] + stick;
                     REAL(spp5Discards)[si + s*NOSIZES + t*aSimNumber*NOSIZES] = stick * aSizeSppInc[4];
                  }
                  Dprobl = Dprobup;
              }
          }

          for (int si = 0; si < NOSIZES; si++){   
              REAL(spp1LndHold)[si + s*NOSIZES + t*aSimNumber*NOSIZES]  = aSpp1LndHold[si] * aSizeSppInc[0];
              REAL(spp2LndHold)[si + s*NOSIZES + t*aSimNumber*NOSIZES]  = aSpp2LndHold[si] * aSizeSppInc[1];
              REAL(spp3LndHold)[si + s*NOSIZES + t*aSimNumber*NOSIZES]  = aSpp3LndHold[si] * aSizeSppInc[2];
              REAL(spp4LndHold)[si + s*NOSIZES + t*aSimNumber*NOSIZES]  = aSpp4LndHold[si] * aSizeSppInc[3];
              REAL(spp5LndHold)[si + s*NOSIZES + t*aSimNumber*NOSIZES]  = aSpp5LndHold[si] * aSizeSppInc[4];
              REAL(spp1DisHold)[si + s*NOSIZES + t*aSimNumber*NOSIZES]  = aSpp1DisHold[si] * aSizeSppInc[0];
              REAL(spp2DisHold)[si + s*NOSIZES + t*aSimNumber*NOSIZES]  = aSpp2DisHold[si] * aSizeSppInc[1];
              REAL(spp3DisHold)[si + s*NOSIZES + t*aSimNumber*NOSIZES]  = aSpp3DisHold[si] * aSizeSppInc[2];
              REAL(spp4DisHold)[si + s*NOSIZES + t*aSimNumber*NOSIZES]  = aSpp4DisHold[si] * aSizeSppInc[3];
              REAL(spp5DisHold)[si + s*NOSIZES + t*aSimNumber*NOSIZES]  = aSpp5DisHold[si] * aSizeSppInc[4];
          }
      }
  }
  
  SET_SLOT(ReturnObject, install("choice"),  choice);
  SET_SLOT(ReturnObject, install("spp1Rand"),   spp1Rand);
  SET_SLOT(ReturnObject, install("spp2Rand"),   spp2Rand);
  SET_SLOT(ReturnObject, install("spp3Rand"),   spp3Rand);
  SET_SLOT(ReturnObject, install("spp4Rand"),   spp4Rand);
  SET_SLOT(ReturnObject, install("spp5Rand"),   spp5Rand);
  SET_SLOT(ReturnObject, install("spp1Landings"),   spp1Landings);
  SET_SLOT(ReturnObject, install("spp2Landings"),   spp2Landings);
  SET_SLOT(ReturnObject, install("spp3Landings"),   spp3Landings);
  SET_SLOT(ReturnObject, install("spp4Landings"),   spp4Landings);
  SET_SLOT(ReturnObject, install("spp5Landings"),   spp5Landings);
  SET_SLOT(ReturnObject, install("spp1LndHold"),  spp1LndHold);
  SET_SLOT(ReturnObject, install("spp2LndHold"),  spp2LndHold);
  SET_SLOT(ReturnObject, install("spp3LndHold"),  spp3LndHold);
  SET_SLOT(ReturnObject, install("spp4LndHold"),  spp4LndHold);
  SET_SLOT(ReturnObject, install("spp5LndHold"),  spp5LndHold);
  SET_SLOT(ReturnObject, install("spp1Discards"),   spp1Discards);
  SET_SLOT(ReturnObject, install("spp2Discards"),   spp2Discards);
  SET_SLOT(ReturnObject, install("spp3Discards"),   spp3Discards);
  SET_SLOT(ReturnObject, install("spp4Discards"),   spp4Discards);
  SET_SLOT(ReturnObject, install("spp5Discards"),   spp5Discards);
  SET_SLOT(ReturnObject, install("spp1DisHold"),  spp1DisHold);
  SET_SLOT(ReturnObject, install("spp2DisHold"),  spp2DisHold);
  SET_SLOT(ReturnObject, install("spp3DisHold"),  spp3DisHold);
  SET_SLOT(ReturnObject, install("spp4DisHold"),  spp4DisHold);
  SET_SLOT(ReturnObject, install("spp5DisHold"),  spp5DisHold);
  SET_SLOT(ReturnObject, install("effort"),  anEffort);
  
  return ReturnObject;
}


void printPatchDetails ( int aPatch, int aTime, int aNoInc, PTYPE aLndParms){

  Rprintf("Patch ");Rprintf("%d, ",aPatch), Rprintf(" time"); Rprintf("%d\n ",aTime);
  for (int s = 0; s < NOSPEC; s++){
    for (int si = 0; si < NOSIZES; si++){
      for (int inc = 0; inc < aNoInc; inc++){
        Rprintf("%1.4f ",aLndParms[aPatch][aTime][s][si][inc]);
      }
      Rprintf("\n");
    }
  }  
}


void nonZeroRanges ( int aHorizon, int aNoInc, int aNPatch, ATYPE aLndParmsAgg, int whatRangeLT[MAXHORIZON]){

  
  // estimate ranges for which we have nonzeros
  float whatRangeL[aNoInc] = {0};                      // !! remember that the aNoinc we have from the argument is already muliplied by number of size classes in function call
  
  for (int t = 0; t < aHorizon; t++){
    for (int inc = 0; inc < aNoInc; inc++){            // !! remember that the aNoinc we have from the argument is already muliplied by number of size classes in function call
      whatRangeL[inc] = 0;
    }
    for (int i = 0; i < aNPatch; i++){
      for (int s = 0; s < 2; s++){                     // only loop over sp0 and 1, because those are the only ones used in backward
	for (int inc = 0; inc < aNoInc; inc++){
	  whatRangeL[inc] += aLndParmsAgg[i][t][s][inc]; 
	}
      }
    }
    int inc = aNoInc;
    while ( whatRangeL[inc]<0.00000000001){ inc--; };
    whatRangeLT[t]= std::min(inc + 20, aNoInc ) ;
  }
   
}


extern "C" SEXP DynStateF(SEXP cLndParms, SEXP cDisParms, SEXP eParms, SEXP pParms, SEXP xControl) {

Rprintf("Start of DynStateF\n");

  /*******************************************************************************************************************/
  /* INITIALISE VARIABLES, READ VARIABLES,                                                                           */
  /*******************************************************************************************************************/
  SEXP ReturnObject, Simulations ;
 
  SEXP a       = GET_DIM(cLndParms);
  int kNPatch  = INTEGER(a)[0];
  int kHorizon = INTEGER(a)[1];
  int noInc    = INTEGER(a)[4];
      
  double vSpp1LndQuota     =            REAL(GET_SLOT(xControl,install("spp1LndQuota"    )))[0];
  double vSpp1LndQuotaFine =            REAL(GET_SLOT(xControl,install("spp1LndQuotaFine")))[0];
  double vSpp2LndQuota     =            REAL(GET_SLOT(xControl,install("spp2LndQuota"    )))[0];
  double vSpp2LndQuotaFine =            REAL(GET_SLOT(xControl,install("spp2LndQuotaFine")))[0];
  int    kSimNumber        =    (short) INTEGER(GET_SLOT(xControl,install("simNumber"    )))[0];
  double sigma             =            REAL(GET_SLOT(xControl,install("sigma"           )))[0];
  double kPriceEffort      =            (REAL(GET_SLOT(xControl,install("fuelUse"      )))[0] * REAL(GET_SLOT(xControl,install("fuelPrice")))[0])  + REAL(GET_SLOT(xControl,install("gearMaintenance")))[0] ;
  double sizeSpp1Inc       =            REAL(GET_SLOT(xControl,install("spp1Incs"     )))[0];
  double sizeSpp2Inc       =            REAL(GET_SLOT(xControl,install("spp2Incs"     )))[0];
  double sizeSpp3Inc       =            REAL(GET_SLOT(xControl,install("spp3Incs"     )))[0];
  double sizeSpp4Inc       =            REAL(GET_SLOT(xControl,install("spp4Incs"     )))[0];
  double sizeSpp5Inc       =            REAL(GET_SLOT(xControl,install("spp5Incs"     )))[0];
  int    verbose           =    (short) INTEGER(GET_SLOT(xControl,install("verbose"   )))[0];
  int    numthreads        =    (short) INTEGER(GET_SLOT(xControl,install("numThreads")))[0];
  
  double   kSpp1LndQuota      = vSpp1LndQuota/sizeSpp1Inc;
  double   kSpp1LndQuotaFine  = vSpp1LndQuotaFine * sizeSpp1Inc;
  double   kSpp2LndQuota      = vSpp2LndQuota/sizeSpp2Inc;
  double   kSpp2LndQuotaFine  = vSpp2LndQuotaFine * sizeSpp2Inc;

  Rprintf("kPriceEffort ");     Rprintf("%f \n",kPriceEffort); 

  int Lndspp1, Lndspp2,i,t,s,inc0,inc1,inc2,inc3,inc4; /*for loops */

  /********************************************************************************************************************/
  /* Check if sigma is zero (because this leads to div by zero in prob calcs). If so set to very small num and warn   */
  /********************************************************************************************************************/
  if (sigma == 0){
    Rprintf("Warning, sigma = 0, giving it a small value (0.000001) \n");
    sigma = 0.000001;
  }
  
  /********************************************************************************************************************/
  /* Check if khorizon is not bigger than HORIZON, if so simulations will not extend arrays                           */
  /********************************************************************************************************************/
  if (kHorizon > MAXHORIZON) Rprintf("Warning, horizon dimensions exceed array \n");
  if (noInc > MAXNOINC) Rprintf("Warning, noinc larger than MAXNOINC \n");
 
  /********************************************************************************************************************/
  /* INITIALISE SRAND FOR RANDOM NUMBER GENERATOR                                                                     */
  /********************************************************************************************************************/
  // srand (1);
  srand(time(0));
  
  /**********************************************************************************************************************/
  /* DEFINE INPUT ARRAYS AND CHECK                                                                                      */
  /**********************************************************************************************************************/

  PTYPE theLndParms         = (PTYPE)  malloc((size_t)sizeof(*theLndParms) *  kNPatch);
  PTYPE theDisParms         = (PTYPE)  malloc((size_t)sizeof(*theDisParms) *  kNPatch); 
  ATYPE theLndParmsAgg      = (ATYPE)  malloc((size_t)sizeof(*theLndParmsAgg) * kNPatch); // matalloc(sizeof(float), (void *)0, 4, kNPatch, MAXHORIZON, NOSPEC,  (NOSIZES*noInc) - 1);
  PITYPE thePriceParms      = (PITYPE) malloc((size_t)sizeof(*thePriceParms) * MAXHORIZON); // matalloc(sizeof(float), (void *)0, 3, MAXHORIZON, NOSPEC, NOSIZES);
    
  int theEffortArray[320000][MAXHORIZON];
  double *theIncrementArray     = (double *) malloc((size_t)NOSPEC * sizeof (double));
  double *theShortTermGains     = (double *) malloc((size_t)kNPatch * sizeof (double));
  double *theShortTermCosts     = (double *) malloc((size_t)kNPatch * sizeof (double));
  double *theShortTermCrewShare = (double *) malloc((size_t)kNPatch * sizeof (double));
  double *theShortTermEcon      = (double *) malloc((size_t)kNPatch * sizeof (double));

  if (theLndParms == 0l) Rprintf("error in memory allocation theLndParms \n");
  if (theDisParms == 0l) Rprintf("error in memory allocation theDisParms \n");

  /*********************************************************************************************************************/
  /* PUT Parmdata from xParms in theLndParms and theDisParms arrays                                                    */
  /*********************************************************************************************************************/  
  for (int i = 0; i < kNPatch; i++){
    for (int t = 0; t < kHorizon; t++){
      for (int s = 0; s < NOSPEC; s++){
        for (int si = 0; si < NOSIZES; si++){
	  for (int inc = 0; inc < noInc; inc++){
	      theLndParms[i][t][s][si][inc] = REAL(cLndParms)[i + t*kNPatch + s*kHorizon*kNPatch + si*NOSPEC*kHorizon*kNPatch +
					                      inc*NOSIZES*NOSPEC*kHorizon*kNPatch]; 
	      theDisParms[i][t][s][si][inc] = REAL(cDisParms)[i + t*kNPatch + s*kHorizon*kNPatch + si*NOSPEC*kHorizon*kNPatch +
								inc*NOSIZES*NOSPEC*kHorizon*kNPatch]; 
	  }
        }
      }
    }
  }

//  printPatchDetails ( 0, 0, noInc,  theLndParms);
  
  /*************************************************************************************************************************************/
  /* PUT theLndParms & theDisParms into aggregate arrays DOES NOT SCALE AUTOMATICALLY WITH NOSIZES BUT IT SHOULD IF THOSE ARE ALTERED  */
  /*************************************************************************************************************************************/  
  for (i = 0; i < kNPatch; i++){
    for (t = 0; t < kHorizon; t++){
      for (s = 0; s < NOSPEC; s++){
	for (int inc = 0; inc < ((NOSIZES * noInc)-1); inc++){
	  theLndParmsAgg[i][t][s][inc] =  0; 
	}
      }
    }
  }

  /* turned around inc0 and patch in loo pfor more parallel speedup */
#pragma omp parallel private(i,t,s,inc0,inc1,inc2,inc3,inc4)
{
#pragma omp for schedule(static)    
  for ( inc0 = 0; inc0 < noInc; inc0++){
    for ( i = 0; i < kNPatch; i++){
    for (t = 0; t < kHorizon; t++){
      for (s = 0; s < NOSPEC; s++){
	
	  for (inc1 = 0; inc1 < noInc; inc1++){
	    for (inc2 = 0; inc2 < noInc; inc2++){
	      for (inc3 = 0; inc3 < noInc; inc3++){
	        for (inc4 = 0; inc4 < noInc; inc4++){
	            
		theLndParmsAgg[i][t][s][inc0 + inc1 + inc2 + inc3 + inc4] +=  theLndParms[i][t][s][0][inc0] * theLndParms[i][t][s][1][inc1] *  theLndParms[i][t][s][2][inc2]  * theLndParms[i][t][s][3][inc3] * theLndParms[i][t][s][4][inc4] ;
	      }
	    }
	  }
	}
    }
  }
    }
  }
  }
  
  Rprintf("Generated aggregated distribution functions \n"); R_FlushConsole();
  
  if (verbose == 1){
    Rprintf(" landings probs for size 1-6, choice1 (0), time 0, spec 0\n");
    for (int inc = 0; inc < noInc; inc++){
      Rprintf("%22.22f ",	theLndParms[0][0][0][0][inc]);
      Rprintf("%22.22f ",	theLndParms[0][0][0][1][inc]);
      Rprintf("%22.22f ",	theLndParms[0][0][0][2][inc]);
      Rprintf("%22.22f ",	theLndParms[0][0][0][3][inc]);
      Rprintf("%22.22f ",	theLndParms[0][0][0][4][inc]);
      Rprintf("\n");
    }
    Rprintf(" landings probs for size 1-6, choice2 (1), time 0, spec 0\n");
    for (int inc = 0; inc < noInc; inc++){
      Rprintf("%22.22f ",	theLndParms[1][0][0][0][inc]);
      Rprintf("%22.22f ",	theLndParms[1][0][0][1][inc]);
      Rprintf("%22.22f ",	theLndParms[1][0][0][2][inc]);
      Rprintf("%22.22f ",	theLndParms[1][0][0][3][inc]);
      Rprintf("%22.22f ",	theLndParms[1][0][0][4][inc]);
      Rprintf("\n");
    }
    Rprintf("aggregated landings probs for choice1, choice2, choice3, time 0, spec 0\n");
    for (int inc = 0; inc < ((NOSIZES * noInc)-1); inc++){
      Rprintf("%32.32f ",	theLndParmsAgg[0][0][0][inc]);
      Rprintf("%32.32f ",	theLndParmsAgg[1][0][0][inc]);
      Rprintf("%32.32f ",	theLndParmsAgg[2][0][0][inc]);Rprintf("\n");
    }
  }
  
  /*************************************************************************************************************************************/
  /*  make sure that the aggregated landing probs sum to zero (they sometimes do not because they are result of multiple products)
  /*************************************************************************************************************************************/  
  
  for ( i = 0; i < kNPatch; i++){
    for (t = 0; t < kHorizon; t++){
      for (s = 0; s < NOSPEC; s++){
        double tmp = 0;
        for (int inc = 0; inc <  ((NOSIZES * noInc)-1); inc++){
          tmp +=  theLndParmsAgg[i][t][s][inc];
        }
        for (int inc = 0; inc <  ((NOSIZES * noInc)-1); inc++){
          theLndParmsAgg[i][t][s][inc] = (1./tmp) * theLndParmsAgg[i][t][s][inc]; 
        }
      }
    }
  }
  
    
  Rprintf("Corrected aggregated distribution functions \n"); R_FlushConsole();
  if (verbose == 1){
    Rprintf("corrected aggregated landings probs for choice1, choice2, choice3, time 0, spec 0\n");
    for (int inc = 0; inc < ((NOSIZES * noInc)-1); inc++){
      Rprintf("%32.32f ",	theLndParmsAgg[0][0][0][inc]);
      Rprintf("%32.32f ",	theLndParmsAgg[1][0][0][inc]);
      Rprintf("%32.32f ",	theLndParmsAgg[2][0][0][inc]);Rprintf("\n");
    }
  }
  
  /*************************************************************************************************************************************/
  /*  estimate ranges for which we have nonzeros
  /*************************************************************************************************************************************/  
  
  int whatRangeLT[MAXHORIZON] = {(NOSIZES * noInc) -1};

  nonZeroRanges(kHorizon,(NOSIZES * noInc)-1,kNPatch, theLndParmsAgg, whatRangeLT);
  for ( t = 0; t < kHorizon; t++){
    Rprintf("%d \n", whatRangeLT[t]);
  }
  Rprintf("Defined ranges for distribution functions per timestep\n"); R_FlushConsole();
  
  int kSpp1Capacity = accumulate(whatRangeLT,whatRangeLT+kHorizon,0) +2;                         // calc max number of increments to loop over for this timest
  //int kSpp1Capacity = SPP1CAPACITY;
  int kSpp2Capacity = kSpp1Capacity;
    
  /*********************************************************************************************************************/
  /* DEFINE OUTPUT ARRAYS AND CHECK                                                                                    */
  /*********************************************************************************************************************/
 
  FCTYPE theFF0             = (FCTYPE) malloc((size_t)sizeof(*theFF0) *  kSpp1Capacity); //   (FCTYPE) matalloc(sizeof(float), (void*) 0,3, kSpp1Capacity, kSpp2Capacity,kNPatch) ;
  FCTYPE numerator          = (FCTYPE) malloc((size_t)sizeof(*numerator) * kSpp1Capacity); //   (FCTYPE) matalloc(sizeof(float), (void*) 0,3, kSpp1Capacity, kSpp2Capacity,kNPatch) ;
  FTYPE theFF0Star          = (FTYPE)  malloc((size_t)sizeof(*theFF0Star) * kSpp1Capacity); //matalloc(sizeof(float), (void *)0, 2, kSpp1Capacity, kSpp2Capacity);
  FTYPE theFF1              = (FTYPE)  malloc((size_t)sizeof(*theFF1) * kSpp1Capacity);    //  (FTYPE) matalloc(sizeof(float), (void *)0, 2, kSpp1Capacity, kSpp2Capacity);
  FTYPE deltaSum            = (FTYPE)  malloc((size_t)sizeof(*deltaSum) * kSpp1Capacity); //(FTYPE) matalloc(sizeof(float), (void *)0, 2, kSpp1Capacity, kSpp2Capacity);
  ITYPE theProbChoice       = (ITYPE)  malloc((size_t)sizeof(*theProbChoice) * MAXHORIZON); //matalloc(sizeof(float), (void *)0, 4, MAXHORIZON, kSpp1Capacity,kSpp2Capacity,kNPatch);

  if (theFF0 == 0l) Rprintf("error in memory allocation FF0 \n");
  if (theFF0Star == 0l) Rprintf("error in memory allocation FF0 \n");
  if (theFF1 == 0l) Rprintf("error in memory allocation FF1 \n");
  if (theProbChoice == 0l) Rprintf("error in memory allocation ProbChoice \n");
  
  /*************************************************************************************************************************************/
  /* PUT incement size information per species in a a vector so we can use it when putting price data in thePriceParns                 */
  /*************************************************************************************************************************************/  
  vector <float> sizeSppInc (NOSPEC,0); 

  sizeSppInc[0] = sizeSpp1Inc;
  sizeSppInc[1] = sizeSpp2Inc;
  sizeSppInc[2] = sizeSpp3Inc;
  sizeSppInc[3] = sizeSpp4Inc;
  sizeSppInc[4] = sizeSpp5Inc;
  
  /*************************************************************************************************************************************/
  /* PUT pricedata from pParms in epriceParms array                                                                                    */
  /*************************************************************************************************************************************/  
  for ( t = 0; t < kHorizon; t++){
    for ( s = 0; s < NOSPEC; s++){
      for (int si = 0; si < NOSIZES; si++){
          thePriceParms[t][s][si] = REAL(pParms)[ si  + t*NOSIZES + s*kHorizon*NOSIZES] * sizeSppInc[s]; //check how r part is structured
      }
    }
  } 

//  Rprintf("Price parameters for all species\n");
//  for (int t = 0; t < kHorizon; t++){
//    for (int s = 0; s < NOSPEC; s++){
//      for (int si = 0; si < NOSIZES; si++){
//	Rprintf ("%f ", thePriceParms[t][s][si]);
//      }
//      Rprintf("\n");
//    }
//  }
  
  /*************************************************************************************************************************************/
  /*INITIALISE THE EFFORTCOST FOR THE DIFFERENT PATCHES                                                                                */
  /*************************************************************************************************************************************/
  Rprintf("Effort array initialised\n");
  for ( i = 0; i < kNPatch; i++){ 
    for ( t = 0; t < kHorizon; t++){
      theEffortArray[i][t] = REAL(eParms)[i + t* kNPatch];
    }
  }

  /*************************************************************************************************************************************/
  /*DEBUGCHECK                                                                                                                         */
  /*************************************************************************************************************************************/
  Rprintf("numPatches "); Rprintf("%d\n", kNPatch);
  Rprintf("vSpp1LndQuota ");  Rprintf("%f", vSpp1LndQuota);   Rprintf("; vSpp2LndQuota ");   Rprintf("%f\n", vSpp2LndQuota);
  Rprintf("kSpp1LndQuota ");  Rprintf("%f", kSpp1LndQuota);   Rprintf("; kSpp2LndQuota ");   Rprintf("%f\n", kSpp2LndQuota);
  Rprintf("noInc ");          Rprintf("%d\n", noInc);
  Rprintf("Capacity ");       Rprintf("%d\n", kSpp1Capacity);
  Rprintf("sizeSpp1Inc   ");                  Rprintf("%f\n", sizeSpp1Inc);

  for ( s = 0; s < NOSPEC; s++){ 
    Rprintf("spp");  Rprintf("%d", s+1); Rprintf(" Increment size "); Rprintf("%f\n",  sizeSppInc[s]);
  }

  /*************************************************************************************************************************************/
  /* INITIALISE F0 AND F1 ARRAY TO 0 IF IN OPENMP ENVIRONMENT                                                                          */
  /*************************************************************************************************************************************/
  Rprintf("get_max threads ");
  Rprintf("%d\n",omp_get_max_threads());

  omp_set_num_threads(numthreads);
  Rprintf("OPEN_MP environment: number of threads is "); Rprintf("%d\n",numthreads);
  //Rprintf("%d\n",omp_get_num_threads());

  
  #pragma omp parallel private(Lndspp1,Lndspp2)
  {
#pragma omp for schedule(static)    
    for (Lndspp1 = 0; Lndspp1 < kSpp1Capacity; Lndspp1++) {
      for (Lndspp2 = 0; Lndspp2 < kSpp2Capacity; Lndspp2++) {
	theFF0Star[Lndspp1][Lndspp2] = 0.0;
	theFF1[Lndspp1][Lndspp2] = 0.0; 
      }
    }
  }
  Rprintf("Init memory FF0Star and FF1 \n");   R_FlushConsole();

  #pragma omp parallel private(Lndspp1,Lndspp2)
  {
#pragma omp for schedule(static)    
    for (Lndspp1 = 0; Lndspp1 < kSpp1Capacity; Lndspp1++) {
      for (Lndspp2 = 0; Lndspp2 < kSpp2Capacity; Lndspp2++) {
	for (int ppp = 0; ppp < kNPatch; ppp++) {
	  theFF0[Lndspp1][Lndspp2][ppp] = 0.0; 
	}
      }
    }
  }
  Rprintf("Init memory FF0 \n");   R_FlushConsole();

  
  /*************************************************************************************************************************************/
  /* INITIALISE THE F1 ARRAY AT FINAL TIMESTEP (FINAL FITNESS FUNCTION)                                                                */
  /*************************************************************************************************************************************/
  #pragma omp parallel private(Lndspp1,Lndspp2)
  {
  #pragma omp for schedule(static)
    for (Lndspp1 = 0; Lndspp1 < kSpp1Capacity; Lndspp1++) {
      for (Lndspp2 = 0; Lndspp2 < kSpp2Capacity; Lndspp2++) {
	theFF1[Lndspp1][Lndspp2] = utility(Lndspp1, kSpp1LndQuota, kSpp1LndQuotaFine, Lndspp2, kSpp2LndQuota, kSpp2LndQuotaFine); 
      }
    }
  }
  if (verbose == 1){
     Rprintf("\n",Lndspp1);
    for (Lndspp1 = 0; Lndspp1 < kSpp1Capacity; Lndspp1 +=15) {
      Rprintf("%d ",Lndspp1);
      for (Lndspp2 = 0; Lndspp2 < kSpp2Capacity; Lndspp2 +=50) {
	Rprintf("%4.f ",theFF1[Lndspp1][Lndspp2]);
      }
      Rprintf("\n"); R_FlushConsole();      
    }
  }
  
  Rprintf("Utility function initialised \n");  R_FlushConsole();
  
  /*************************************************************************************************************************************/
  /* MAIN PART OF MODEL :  do backward calculations (RLE encoding, and swapping F arrays)                      */
  /*************************************************************************************************************************************/
  vector <int> maxspp (MAXHORIZON,0);
  vector < vector<UFINT> > val(MAXHORIZON, vector<UFINT> (1,0));                     // for RLE
  vector < vector<unsigned short> > rep(MAXHORIZON, vector<unsigned short> (1,1));   // for RLE
  
  for ( t = kHorizon - 1; t >= 0; t--) {                                          // start the backward calculation

    Rprintf("Timestep %d Calc short term economics ", t+1);  R_FlushConsole();       // for each choice calc mean short term economic gains and costs  
    for ( i = 0; i < kNPatch; i++){     
      theShortTermGains[i] = shortTermGains(t, noInc, theLndParms, i, thePriceParms);
      theShortTermCosts[i] = shortTermCosts(t, i, kPriceEffort, theEffortArray);     //*, noInc, theLndParms, thePriceParms);
      theShortTermCrewShare[i] = theShortTermGains[i]*0.0;     //*, noInc, theLndParms, thePriceParms);
      theShortTermEcon[i]  =  theShortTermGains[i] - theShortTermCosts[i] - theShortTermCrewShare[i];
    }
    
    if (verbose == 1){
      Rprintf("\n"); R_FlushConsole();        
      for ( i = 0; i < kNPatch; i++){     
        Rprintf("%2.2f ",theShortTermGains[i]);
        Rprintf("%2.2f ",theShortTermCosts[i]);
        Rprintf("%2.2f ",theShortTermEcon[i]);
        Rprintf("\n "); R_FlushConsole();
      }
    }
    
    maxspp[t]  = accumulate(whatRangeLT,whatRangeLT+t,0) + 2;                         // calc max number of increments to loop over for this timest
    //maxspp[t]  = (t +1 )* noInc *NOSIZES;                         // calc max number of increments to loop over for this timest 
    Rprintf(" maxspp %d ", maxspp[t]);  R_FlushConsole();
    
    // do backward calcs 
    #pragma omp parallel private(Lndspp1,Lndspp2)
    {
    #pragma omp for schedule(static)    
      for (Lndspp1 = 0; Lndspp1 < maxspp[t]; Lndspp1++) {
	for (Lndspp2 = 0; Lndspp2 < maxspp[t]; Lndspp2++) {
	  FFF(Lndspp1, Lndspp2, whatRangeLT[t],theLndParmsAgg, t, kNPatch,  theShortTermEcon, theFF0, theFF0Star, theFF1, theProbChoice, verbose);
	  //FFF(Lndspp1, Lndspp2, (NOSIZES*noInc)-1,theLndParmsAgg, t, kNPatch,  theShortTermEcon, theFF0, theFF0Star, theFF1, theProbChoice);
	}
      }
    }
  
    if (verbose == 1){
      Rprintf("\n FF0 \n"); R_FlushConsole();          
      for (Lndspp1 = 0; Lndspp1 < maxspp[t]; Lndspp1 +=15) {
	Rprintf("%d ",Lndspp1); 
	for (Lndspp2 = 0; Lndspp2 < maxspp[t]; Lndspp2 +=50) {
	  for (int ppp = 0; ppp < kNPatch; ppp++) {
	    Rprintf("%2.2f,", theFF0[Lndspp1][Lndspp2][ppp]);
	  }
         Rprintf(" "); 
	}
	Rprintf("\n"); R_FlushConsole();      
      }
    } 
    
    Rprintf("Finished FF0 calculations ");  R_FlushConsole();

    /***********************************************************************************************************************************/
    /* CALC delta (note that exp(-deltas/sigmas) are stored in FF0 array first                                                              */
    /***********************************************************************************************************************************/
 
    for (Lndspp1 = 0; Lndspp1 < maxspp[t]; Lndspp1++) {
      for (Lndspp2 = 0; Lndspp2 < maxspp[t]; Lndspp2++) {
	for (int ppp = 0; ppp < kNPatch; ppp++) {
	  numerator[Lndspp1][Lndspp2][ppp] = exp(-(theFF0Star[Lndspp1][Lndspp2] - theFF0[Lndspp1][Lndspp2][ppp])/sigma);
	}}}


    //  Rprintf(" Finished numerator "); R_FlushConsole();
    if (verbose == 1){
      Rprintf("\n FF0star \n"); R_FlushConsole();   
      for (Lndspp1 = 0; Lndspp1 < maxspp[t]; Lndspp1 += 15) {
	Rprintf("%d ",Lndspp1); 
	for (Lndspp2 = 0; Lndspp2 < maxspp[t]; Lndspp2 +=50) {
	  Rprintf("%2.2f ",theFF0Star[Lndspp1][Lndspp2]); 
	}
	Rprintf("\n "); R_FlushConsole();
      }
      Rprintf("\n ");R_FlushConsole(); 
    }
 
    /* first set deltaSum to zero before accumulating numerators into it (otherwise values from other time steps will be in there) */
    for (Lndspp1 = 0; Lndspp1 < maxspp[t]; Lndspp1++) {
      for (Lndspp2 = 0; Lndspp2 < maxspp[t]; Lndspp2++) {
	  deltaSum[Lndspp1][Lndspp2] = 0.0;
	}}

    /* accumulate numerators */
    for (Lndspp1 = 0; Lndspp1 < maxspp[t]; Lndspp1++) {
      for (Lndspp2 = 0; Lndspp2 < maxspp[t]; Lndspp2++) {
	for (int ppp = 0; ppp < kNPatch; ppp++) {
	  deltaSum[Lndspp1][Lndspp2] += numerator[Lndspp1][Lndspp2][ppp];
	}}}

    /* calc probabilities by dividing numerator by deltasum */
    for (Lndspp1 = 0; Lndspp1 < maxspp[t]; Lndspp1++) {
      for (Lndspp2 = 0; Lndspp2 < maxspp[t]; Lndspp2++) {
	for (int ppp = 0; ppp < kNPatch; ppp++) {
	  // Rprintf("\n %d ",Lndspp1); Rprintf(" %d ",Lndspp2); Rprintf(" %d ",ppp);  R_FlushConsole();    
	  theProbChoice[t][Lndspp1][Lndspp2][ppp] = numerator[Lndspp1][Lndspp2][ppp]/  deltaSum[Lndspp1][Lndspp2] ;
          // Rprintf("%.3f ",theProbChoice[t][Lndspp1][Lndspp2][ppp]);  R_FlushConsole();    
	}}}

    if (verbose == 1){
      Rprintf("\n probabilities \n"); R_FlushConsole();    
      for (Lndspp1 = 0; Lndspp1 < maxspp[t]; Lndspp1 += 15) {
	Rprintf("%d ",Lndspp1); 
	for (Lndspp2 = 0; Lndspp2 < maxspp[t]; Lndspp2 +=50) {
	  for (int ppp = 0; ppp < kNPatch; ppp++) {
	    Rprintf("%.2f,",theProbChoice[t][Lndspp1][Lndspp2][ppp]); 
	  }
	  Rprintf(" "); R_FlushConsole();
	}
	Rprintf("\n ");R_FlushConsole(); 
      }
    }
    
    Rprintf("Probability calculations done \n"); R_FlushConsole();
 
    /***********************************************************************************************************************************/
    /* PUT FO ARRAY IN F1 ARRAY AFTER TIMESTEP IS CALCULATED                                                                           */
    /***********************************************************************************************************************************/
#pragma omp parallel private(Lndspp1,Lndspp2)
    {
    #pragma omp for schedule(static)
    for (Lndspp1 = 0; Lndspp1 < maxspp[t]; Lndspp1++) {
      for (Lndspp2 = 0; Lndspp2 < maxspp[t]; Lndspp2++) {
               theFF1[Lndspp1][Lndspp2] =  theFF0Star[Lndspp1][Lndspp2];
      }
    }
    }
  }
  Rprintf("Backward calculations done \n"); R_FlushConsole();
 
  matfree(theFF0);
  matfree(theFF1);
  matfree(theFF0Star);
  matfree(deltaSum);

  /***************************************************************************************************************/
  /* DO MONTE CARLO SIMULATIONS                                                                                  */
  /***************************************************************************************************************/
  Simulations = SimulateF(kSimNumber, kHorizon, theProbChoice, kNPatch, noInc, sizeSppInc, theLndParms, theDisParms, theEffortArray, verbose);

  Rprintf("Simulations done \n"); R_FlushConsole();
  PROTECT(ReturnObject = NEW_OBJECT(MAKE_CLASS("DynState")));
  SET_SLOT(ReturnObject, install("sim"),  Simulations);
  UNPROTECT(31);
  matfree(theProbChoice);
  matfree(theLndParms);
  matfree(theDisParms);
  matfree(theLndParmsAgg);
  return ReturnObject;
}
