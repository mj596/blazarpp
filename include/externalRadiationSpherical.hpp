/**
   @file externalRadiationSpherical.hpp
   @author Mateusz Janiak
*/

#ifndef _externalradiationspherical_hpp
#define _externalradiationspherical_hpp 1

#include "energyDissProc.hpp"
#include "inverseCompton.hpp"

class electrons;

/** @class externalRadiationSpherical
    Class provides electron energy losses and luminosity calculation for inverse Compton process (ERC) 
    for seed photons comming from BLR and HDR modelled as spherical sources 
*/

class externalRadiationSpherical : public energyDissProc {
 public:
  /* requested parameters */
  double cf, e;
  double Gamma, thetaObs;
  int KN;

  /* own parameters */
  gsl_vector* upe; // upe vector vs radius r  

  /** constructor
      @param _cfg - scfgp class object
      @param _r - jetGeometry class object
      @param _ele - electrons class object
      @param _id */
  externalRadiationSpherical( scfgp* _cfg, jetGeometry* _r, electrons* _ele, std::string _id );
  /** destructor */
  ~externalRadiationSpherical( );

  /* virtual methods are always defined in sub-classes */
  /** sets rext if not provided by config file */
  virtual void setRadius( ) { }
  virtual void printInfo( ) { }

  /** calculate and set Upe every time with new radius r */
  virtual void update( ) { }

  /** now obsolete
      virtual double getdLdlnr( double _r ) { return 0; } */
  
  /** get v_ext value in Hz
      @param _r - radius
      @returns v_ext */
  virtual double getvext( double _r ) { return 0; }

  /* common methods for sub-classes */
  double dotg( double g );

  /** calculate intrinsic luminosity
      @param v - frequency (jet co-moving frame)
      @returns L'_v */
  double calculateLpv( double v, double theta );
  void setLpv( );
};

#endif /* _externalradiationspherical_hpp */
