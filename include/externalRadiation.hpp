/**
   @file externalRadiation.hpp
   @author Mateusz Janiak
*/

#ifndef _externalradiation_hpp
#define _externalradiation_hpp 1

#include "energyDissProc.hpp"
#include "inverseCompton.hpp"

class electrons;

/** @class externalRadiationSpherical
    Class provides electron energy losses and luminosity calculation for inverse Compton process (ERC) 
    for seed photons comming from BLR and HDR modelled as spherical sources. Implementation provided by the class is identical to this in original BLAZAR code from Moderski et al. (2003).
    it is not modified OR checked for consistency. Use with caution!
*/

class externalRadiation : public energyDissProc {
  double exte, extu, extr, extk;
  double normA,TempX,ksi;  
  double Gamma, thetaObs;
  int KN;

 public:
  externalRadiation( scfgp* _cfg, jetGeometry* _r, electrons* _ele, std::string _id );
  ~externalRadiation();

  double dotg( double g );
  void update(  );
  void setLpv( );
  double calculateLpv( double v, double theta );
  void printInfo( );
};

#endif /* _externalradiation_hpp */
