/**
   @file BlrSphericalGu.hpp
   @author Mateusz Janiak
*/

#ifndef _BlrSphericalGu_hpp
#define _BlrSphericalGu_hpp 1

#include "externalRadiationGu.hpp"
#include "inverseCompton.hpp"

class electrons;

/** @class BlrSphericalGu
    Parameters used by class (defaults)
    @param e (externalRadiationGu) - energy (ext) in eV (e = 10.0)
    @param cf (externalRadiationSphericalGu) - covering factor (cf = 0.1)
    @param rext - blr radius (rext = 0.1*r_sub)
    @param k - u_ext vs r index = 3 */

class BlrSphericalGu : public externalRadiationGu {
  /* requested parameters */
  double cf, e;
  double rext, k;

 public:
  /** constructor
      @param _cfg - scfgp class object
      @param _r - jetGeometry class object
      @param _ele - electrons class object
      @param _id */
  BlrSphericalGu( scfgp* _cfg, jetGeometry* _r, electrons* _ele, std::string _id );
  
  /** destructor */
  ~BlrSphericalGu( ) { }

  void printInfo( );

  /** calculate and set Upe every time with new radius r */
  void update( );

  /** get v_ext value in Hz
      @param _r - radius
      @returns v_ext */
  double getvext( double _r );

  /** sets rext if not provided by config file */
  void setRadius( ); 
};

#endif /* _BlrSphericalGu_hpp */



