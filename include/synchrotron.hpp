/**
    @file synchrotron.hpp
    @author Mateusz Janiak
*/

#ifndef _synchrotron_hpp
#define _synchrotron_hpp 1

#include "energyDissProc.hpp"
#include "magneticField.hpp"
#include <gsl/gsl_sf_bessel.h>

class electrons;

/** @class synchrotron
    Class provides synchrotron electron energy losses
    and synchrotron luminosity calculation */

class synchrotron : public energyDissProc {
  /* requested parameters */
  std::string lumModel;
  double thetaJ, Gamma;
  int SABS;
  int setSSCBlob;

  /* magentic field */
  magneticField* B;

  /* function to calculate synchrotorn luminosity */
  double FS(double t, double* jS, double* sigmaS);

 public:
  /** constructor
      @param _cfg - scfgp class object
      @param _r - jetGeometry class object
      @param _ele - electrons class object
      @param _B - magnetic field class object
      @param _id */
  synchrotron( scfgp* _cfg, jetGeometry* _r, electrons* _ele, magneticField* _B, std::string _id );

  /** destructor */
  ~synchrotron( );

  void printInfo();
  
  /** iterate synchrotron and SSC calculation to achieve steady state and balance between synchrotron and SSC luminosities */
  double iterate( );

  double dotg( double g );
  void update( );
  void setLpv( );
  double calculateLpv( double v ); };

#endif /* _synchrotron_hpp */


