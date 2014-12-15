/**
   @file DtPlanar.hpp
   @author Mateusz Janiak
*/

#ifndef _DtPlanar_hpp
#define _DtPlanar_hpp 1

#include "externalRadiationPlanar.hpp"

class electrons;

/** @class BlrPlanar 
    Parameters used by class (defaults)
    @param alpha (externalRadiationPlanar) ( alpha = 0.0 )
    @param R1 - (externalRadiationPlanar) dt lower boundary radius (R1 = R_sub)
    @param R2 - (externalRadiationPlanar) dt upper boundary radius (R2 = 10.0*R_sub)
    @param s - (externalRadiationPlanar) stratification index (s = 1.0)
    @param e - energy (ext) in eV (e = 0.6203 @ Rext ~ 1.5e14 Hz)
    @param CF - covering factor (CF = 0.1) */

class DtPlanar : public externalRadiationPlanar {
  /* parameters read from config file */
  double e, cf; 
  
  /* own parameters */
  double vext;

 public:
  /** constructor
      @param _cfg - scfgp class object
      @param _r - jetGeometry class object
      @param _ele - electrons class object
      @param _id */
  DtPlanar( scfgp* _cfg, jetGeometry* _r, electrons* _ele, std::string _id );

  /** destructor */
  ~DtPlanar( );

  void printInfo( ); // prints info on class parameters
  void setRadius( ); // sets R if not provided by config file
  double getvext( double _R );
  void setdLdR( ); 
};

#endif /* _DtPlanar_hpp */
