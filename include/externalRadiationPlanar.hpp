/**
   @file externalRadiationPlanar.hpp
   @author Mateusz Janiak
*/

#ifndef _externalRadiationPlanar_hpp
#define _externalRadiationPlanar_hpp 1

#include "energyDissProc.hpp"
#include "logGeometry.hpp"
#include "inverseCompton.hpp"

class electrons;

/** @class externalRadiationPlanar
    Class provides electron energy losses and luminosity calculation for inverse Compton process (ERC) 
    for seed photons comming from accretion disk, BLR and HDR modelled as planar source. Incoming photons are Doppler shifted depending on the geometry. Energy density of seed photons is calculated consistently based on the source characteristics in the accretion disk plane
    @param alpha - radiation vs R slope index
    @param R1 - inner edge
    @param R2 - outer edge
    @param s - radiation vs R slope index
    @param approx - approximation calculations switch
    @param fixedUpe - fix u_ext vs r
    @param saveRm - switch to save information about R_m
    @param R - hirozontal disk plane radius 
    @param rm - R_m vector info
    @param dLdR - vector to store information about dLdR
*/

class externalRadiationPlanar : public energyDissProc {
protected:
  /* requested parameters */
  double alpha, R1, R2, s;
  int approx;
  double fixedUpe;
  double Gamma, thetaObs;
  int KN;
  double mBH, eDisk, mDot;
  int saveRm;
  
  /* own parameters */
  logGeometry* R; /* horizontal disk plane radius R */
  gsl_vector* rm; /* rm vector for storing radius which has maximal Upe(R) */
  gsl_vector* dLdR; /* dLdR vector to store information about dL/dR */
  double Rm, Upe;
  int Rm_index;

 public:
  /** constructor
      @param _cfg - scfgp class object
      @param _r - jetGeometry class object
      @param _ele - electrons class object
      @param _id */
  externalRadiationPlanar( scfgp* _cfg, jetGeometry* _r, electrons* _ele, std::string _id );

  /** destructor */
  ~externalRadiationPlanar( );

  /* virtual methods are always defined in sub-classes */
  virtual void printInfo( ) { }

  /** sets rext if not provided by config file */
  virtual void setRadius( ) { }

 /** get v_ext value in Hz
      @param _R - radius in accretion disk plane
      @returns v_ext */
  virtual double getvext( double _R ) { return 0; }

  /** set dL/dR for particulat source */
  virtual void setdLdR( ) { }

  /* common methods for sub-classes */
  void update( );

  /** caculate u'_ext in the jet co-moving frame using full equations 
      @return u'_ext */
  double calculateUpeFull( );

  /** caculate u'_ext in the jet co-moving frame using approximate equations 
      @param i - index in radius R at which u'_ext is calculated
      @param flag - obsolete ?
      @return u'_ext */
  double calculateUpeApprox( int index, bool flag );

  /** auxillary function to calculate u'_ext 
      @param index - index in radius R at which u'_ext is calculated */
  double function_upe( int index );
  
  /** auxillary function to alculate u'_ext - normalization factor */
  double ctau( );

  /** get fdLdR
      @param index in radius R at which dLdR is calculated */
  double getdLdR( int index );

  /* commmon functions */
  double dotg( double g );
  void setLpv( );

  /** calculate intrinsic luminosity
      @param v - frequency (jet co-moving frame)
      @returns L'_v */
  double calculateLpv( double v, double theta );

  /** function that calculates intrinsic ERC luminosity
      @param v - frequency
      @param theta - observer angle
      @param R - radius in the accretion disk plane
      @returns L' */
  double fsc_integral( double v, double theta, double R );  

  /** auxillary function to calculate L'
      @param index - index in radius R at which L' is calculated
      @param v - frequency
      @param theta - observer angle  */
    double function_lum( int index, double v, double theta );

  /* memory allocation functions */
  void allocateRm( );
  void allocatedLdR( );
  void freeRm( );

  /** get R_m - radius at which contribution of u'_ext is maximal
      @returns R_m */
  double getRm( ) { return gsl_vector_get( rm, r->getIndex( ) ); }

  
  /** set R_m - radius at which contribution of u'_ext is maximal
      @param val - current R_m */
  void setRm( double val ) { gsl_vector_set( rm, r->getIndex( ), val ); }

  /** get vector with R_m values stored */
  gsl_vector* getRmVector( ) { return rm; }

  /** save vector with R_m values stored */
  void saveRmVsR( );
};

#endif /* _externalRadiationPlanar_H */
