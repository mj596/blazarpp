/**
   @file AccdPlanar.hpp
   @author Mateusz Janiak
*/

#ifndef _AccdPlanar_hpp
#define _AccdPlanar_hpp 1

#include "externalRadiationPlanar.hpp"

class electrons;

/** @class AccdPlanar
    Parameters used by class (defaults)
    @param alpha (externalRadiationPlanar) ( alpha = -0.3333 )
    @param R1 - (externalRadiationPlanar) accd lower boundary radius ('NEW' stratification) (R1 = 0.1*R_sub)
    @param R2 - (externalRadiationPlanar) accd upper boundary radius ('NEW' stratification) (R2 = R_sub)
    @param s - (externalRadiationPlanar) stratification index ('NEW' stratification) (s = 2.0) */

class AccdPlanar : public externalRadiationPlanar {
  /* parameters read from config file */
  /* -- */

  /* own parameters */
  /* -- */

 public:
  /** constructor
      @param _cfg - scfgp class object
      @param _r - jetGeometry class object
      @param _ele - electrons class object
      @param _id */
  AccdPlanar( scfgp* _cfg, jetGeometry* _r, electrons* _ele, std::string _id );

  /** destructor */
  ~AccdPlanar( );

  void printInfo( ); // prints info on class parameters
  void setRadius( ); // sets R if not provided by config file
  double getvext( double _R );
  void setdLdR( );
};

#endif /* _AccdPlanar_H */
