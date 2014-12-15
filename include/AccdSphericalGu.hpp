/**
   @file AccdSphericalGu.hpp
   @author Mateusz Janiak
*/

#ifndef _AccdSphericalGu_hpp
#define _AccdSphericalGu_hpp 1

#include "externalRadiationGu.hpp"
#include "inverseCompton.hpp"

class electrons;

/** @class AccdSphericalGu
    Parameters used by class (defaults)
    @param rext - accretion disk inner radius */

class AccdSphericalGu : public externalRadiationGu {
  /* Parameters used by class (defaults) */
     double rext;

  /* requested parameters */
  
 public:
  /** constructor
      @param _cfg - scfgp class object
      @param _r - jetGeometry class object
      @param _ele - electrons class object
      @param _id */
  AccdSphericalGu( scfgp* _cfg, jetGeometry* _r, electrons* _ele, std::string _id );

  /** destructor */
  ~AccdSphericalGu( ) { }

  void printInfo( ); // prints info on class parameters
  void update( ); // calculate and set Upe every time with new radius r

  /** get v_ext value in Hz
      @param _r - radius
      @returns v_ext */
  double getvext( double _r );

  /** sets rext if not provided by config file */
  void setRadius( ); // sets rext if not provided by config file
};

#endif /* _AccdSphericalGu_hpp */
