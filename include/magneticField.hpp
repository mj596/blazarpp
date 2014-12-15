/**
    @file magneticField.hpp
    @author Mateusz Janiak
*/

#ifndef _magfield_hpp
#define _magfield_hpp 1

#include "baseClass.hpp"

/**
   @class magneticField
   Class defining magnetic field across the jet
*/

class magneticField : public baseClass {
  /* requested parameters */
  double B0, B1, sigmaB, kB, Gamma, thetaJ, injRm;
  double eDiss, mBH, eJet, mDot, eEle;
  std::string magModel;

  /* other */
  double Lb;
  double B0steady;

 public:
  /** consructor
      @param scfgp
      @param jetGeometry
      @param id */
  magneticField( scfgp* _cfg, jetGeometry* _r, std::string _id );

  /** destructor */
  ~magneticField( );

  /** get current value of magnetic field */
  double getB( );

  /** get value of magnetic field at specific radius
      @param radius */
  double getB( double _r );
  
  /** get maximum value of magnetic field (closest radius */
  double getMaxB( );
  
  /** get current value of magnetic energy density */
  double get_uB( );

  /** get value of magnetic energy density at aspecific radius 
      @param radius */
  double get_uB( double _r );

  void printInfo( ); };

#endif /* _magfield_hpp */

