/**
   @file externalRadiationGu.hpp
   @author Mateusz Janiak
*/

#ifndef _externalradiationGu_hpp
#define _externalradiationGu_hpp 1

#include "energyDissProc.hpp"
#include "inverseCompton.hpp"

class electrons;

/** @class externalRadiationGu
    Class provides electron energy losses and luminosity calculation for inverse Compton process (ERC) 
    for seed photons comming from BLR and HDR modelled as quasi-spherical sources with u_ext modified by g_u parameter
*/
class externalRadiationGu : public energyDissProc {
 public:
  /* requested parameters */
  double Gamma, thetaObs;
  double mBH, eDisk, mDot;
  int KN;

  /* own parameters */
  //  gsl_vector* upe; // upe vector vs radius r

  /** constructor
      @param _cfg - scfgp class object
      @param _r - jetGeometry class object
      @param _ele - electrons class object
      @param _id */
  externalRadiationGu( scfgp* _cfg, jetGeometry* _r, electrons* _ele, std::string _id );
  
  /** destructor */
  ~externalRadiationGu( );
  
  virtual void printInfo( ) { }

  /** calculate and set Upe every time with new radius r */
  virtual void update( ) { }
  
  /** get v_ext value in Hz
      @param _r - radius
      @returns v_ext */
  virtual double getvext( double _r ) { }

  /** set radial boundaries for processes */
  virtual void setRadius( ) { }

  /* common methods for sub-classes */
  double dotg( double g );

  /** calculate intrinsic luminosity
      @param v - frequency (jet co-moving frame)
      @returns L'_v */
  double calculateLpv( double v, double theta );
  void setLpv( );
};

#endif /* _externalradiationGu_hpp */
