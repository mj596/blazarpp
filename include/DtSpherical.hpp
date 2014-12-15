/**
    @file DtSpherical.hpp
    @author Mateusz Janiak
*/

#ifndef _DtSpherical_hpp
#define _DtSpherical_hpp 1

#include "externalRadiationSpherical.hpp"
#include "inverseCompton.hpp"

class electrons;

/** @class DtSpherical
    Parameters used by class (defaults)
    @param e (externalRadiationSpherical) - energy (ext) in eV  (e = 0.4136 which 
     corresponds to 1.5e14Hz @ rext)
    @param cf (externalRadiationSpherical) - covering factor (cf = 0.1)
    @param q1 - lower index in stratification (q1 = 0.0)
    @param q2 - upper index in stratification (q2 = -3.0)
    @param rext - dt radius (rext = 1.0*r_sub) */

class DtSpherical : public externalRadiationSpherical {
     /* OBSOLETE */
     /*     - r2 - dt upper boundary radius (r2 = 10.0*r_sub) */

  /* requested parameters */
  /* obsolete  double rext, r2; */
  double q1, q2, rext; 
  double mBH, eDisk, mDot;
  double Gamma;

  /* own parameters */
  double vext; // this is only at rext
  
 public:
  /** constructor
      @param _cfg - scfgp class object
      @param _r - jetGeometry class object
      @param _ele - electrons class object
      @param _id */
  DtSpherical( scfgp* _cfg, jetGeometry* _r, electrons* _ele, std::string _id );

  /** destructor */
  ~DtSpherical( ) { }

  /** calculate and set Upe every time with new radius r */
  void update( );
  void printInfo( );

  /* sets rext and r2 if not provided by config file */
  void setRadius( );
  /* now obsolete 
     double getdLdlnr( double _r ); */

  /** get v_ext value in Hz
      @param _r - radius
      @returns v_ext */
  double getvext( double _r );
};

#endif /* _DtSpherical_hpp */
