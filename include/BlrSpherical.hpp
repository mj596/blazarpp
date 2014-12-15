/**
    @file BlrSpherical.hpp
    @author Mateusz Janiak
*/

#ifndef _BlrSpherical_hpp
#define _BlrSpherical_hpp 1

#include "externalRadiationSpherical.hpp"
#include "inverseCompton.hpp"

class electrons;

/** @class BlrSpherical
    Parameters used by class (defaults)
    @param e (externalRadiationSpherical) - energy (ext) in eV (e = 10.0)
    @param cf (externalRadiationSpherical) - covering factor (cf = 0.1)
    @param q1 - lower index in stratification (q1 = 1.0)
    @param q2 - upper index in stratification (q2 = 3.0)
    @param rext - blr radius (rext = 0.1*r_sub) */

class BlrSpherical : public externalRadiationSpherical {
  /* requested parameters */
  double q1, q2, rext; 
  double mBH, eDisk, mDot;
  double Gamma;

  /* own parameters */
  double vext;

 public:
  /** constructor
      @param _cfg - scfgp class object
      @param _r - jetGeometry class object
      @param _ele - electrons class object
      @param _id */
  BlrSpherical( scfgp* _cfg, jetGeometry* _r, electrons* _ele, std::string _id );
  /** destructor */
  ~BlrSpherical( ) { }

  /** calculate and set Upe every time with new radius r */
  void update( );
  void printInfo( );

  /** sets rext if not provided by config file */
  void setRadius( );

  /* now obsolete */
  //  double getdLdlnr( double _r );
  
  /** get v_ext value in Hz
      @param _r - radius
      @returns v_ext */
  double getvext( double _r );
};

#endif /* _BlrSpherical_hpp */





