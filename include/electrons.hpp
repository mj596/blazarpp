 /**
    @file electrons.hpp
    @author Mateusz Janiak
*/

#ifndef _electrons_hpp
#define _electrons_hpp 1

#include "baseClass.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_linalg.h>

class energyDissProc;
class synchrotron;

/** @class electrons
    Class to hold all detaild for electrons in the jet (all calculated in jet co-moving frame);
    it is responsible for calculating energy losses, solving propagation equation, electron energy spectrum etc.
*/
class electrons : public baseClass
{
  /* requested parameters */
  std::string eleModel, injModel, lumModel;
  double p1, p2, gammaMin, gammaBreak, gammaMax, injRm, injSigma, eleK;
  double eEle, eDiss, Gamma, NeNp;
  double mBH, eDisk, mDot, eJet, sigmaB;
  int evol, saveElectrons, saveElectronsAvg;
  
  /* other */
  double avgGamma;
  double dLogGamma; /* gamma integration steps */  
  double Ljet, Ledd, Ldisk, LjetDiss, Lb, Lprot;

  void allocateGamma( );
  void allocateTridiag( );
  
  void freeGamma( );
  void freeTridiag( );

public:
  /** constructor 
      @param scfgp
      @param jetGeometry
      @param id */
  electrons( scfgp* _cfg, jetGeometry* _r, std::string _id );
  
  /** destructor */
  ~electrons( );

  /** vector to hold active energy dissipation processes */
  std::vector<energyDissProc*> EnDissProc;
  
  /** vectors to store electron spectrum data */
  gsl_vector *gamma, *Ngamma, *NgammaAvg;  

  /** tridag vectors */  
  gsl_vector *A, *B, *C, *S, *U;

  /** electron injection function
      @param g - gamma Lorentz factor
      @param _radius - radius at which injection occurs
      @returns N_gamma */
  double tQ( double g, double _radius );

  /** gaussian radial injection profile
      @param radius */
  double gauss( double _radius );

  /** modified gaussian radial injection profile
      @param radius */
  double gaussmod( double _radius );

  /** triangle radial injection profile
      @param radius */
  double tri( double _radius );

  /** electron injection function */ 
  double inject( );

  /** main electron evolution function; used gsl_tridiag to solve evolution equation */
  void evolve( );

  /** get electron cooling details 
     @param g - electron Lorentz factor 
     @returns d gamma \ d t */
  double dgdr( double g );

  /** get electron gamma factor 
   @param vector index
   @returns electron gamma factor */
  double getGamma( int i );

  /** get electron Ngamma
   @param vector index
   @returns electron Ngamma */
  double getNgamma( int i );

  /** add process to energy dissipation 
      @param energyDissProc */
  void addEnDissProc( energyDissProc* _obj );

  /** list available and active energy dissipation processes */
  void listEnDissProc( );
//  int read( );

  void printInfo( );
  /** print information about cooling details for active processes
      @params electron Lorentz factor */
  void printCoolingInfo( double g );
  
  /** calculate averaged electron distribution */
  void avgNgamma( );

  /** set electron injection parameters */
  void setInjectionParameters( );

  /** set jet energetics */
  void setEnergetics( );

  double getGammaMin( ) { return gammaMin; }
  double getGammaMax( ) { return gammaMax; }

  /** probably doubled from baseClass?? */
  double beta( double x ); 

  /** get d gamma in log scale
      @returns dLogGamma */
  double getdLogGamma( ) { return dLogGamma; }
  
  /** check if we should do electron evolution */
  int ifEvol( ) { return evol; }
  
  /** save injected electron spectrum */
  void saveInjection( );

  /** save current electron spectrum */
  void saveNgamma( );

  /** save averaged electron spectrum */
  void saveNgammaAvg( );
};

/** @class struct gamma_break_params is used to calculate break in electron energy */
struct gamma_break_params
{
  double p1, p2, gamma_min, gamma_max, avg_gamma;
};

/** defines function to be solved to get break in electron energy spectrum
    @param x - electron Lorentz gamma
    @params* - other parameters defined by struct gamma_break_params
    @returns function value */
double functionGammaBreak( double x, void *params );

/** method to calculate gamma break
    @param p1 - spectral index  p1 
    @param p2 - spectral index  p2 
    @gamma_min - minimal electron Lorentz factor
    @gamma_max - max electron Lorentz factor
    @avg_gamma - average electron injecton Lorentz factor
    @return gamma_break
*/
double solveGammaBreak( double p1, double p2, double gamma_min, double gamma_max, double avg_gamma );

/** auxillary for calculating gamma_break */
double f1( double p1, double p2, double gamma_min, double gamma_max, double gamma_break );

/** auxillary for calculating gamma_break */
double f2( double p1, double p2, double gamma_min, double gamma_max, double gamma_break );

#endif /* _electrons_hpp */
