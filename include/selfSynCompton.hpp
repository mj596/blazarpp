/**
   @file selfSynCompton.hpp
   @author Mateusz Janiak
*/

#ifndef _selfsyncompton_hpp
#define _selfsyncompton_hpp 1

#include "energyDissProc.hpp"
#include "magneticField.hpp"
#include "inverseCompton.hpp"

class electrons;

/** @class selfSynCompton
    Class provides SSC electron energy losses
    and SSC luminosity calculation */

class selfSynCompton : public energyDissProc {
  int KN;
  
 public:
  /** magenetic field B' */
  magneticField* B;

  /** synchrotron process pointer */
  energyDissProc* syn;

  /** constructor
      @param _cfg - scfgp class object
      @param _r - jetGeometry class object
      @param _ele - electrons class object
      @param _B - magneticField class object
      @param _id */
  selfSynCompton( scfgp* _cfg, jetGeometry* _r, electrons* _ele, magneticField* _B, std::string _id );

  /** destructor */
  ~selfSynCompton();
  
  void printInfo();
  
  double dotg( double g );
  void update( );
  double calculateLpv( double v );
  void setLpv( );
};

#endif /* _selfsyncompton_hpp */
