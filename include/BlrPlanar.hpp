/**
   @file BlrPlanar.hpp
   @author Mateusz Janiak
*/

#ifndef _BlrPlanar_hpp
#define _BlrPlanar_hpp 1

#include "externalRadiationPlanar.hpp"

class electrons;

/** @class BlrPlanar 
    Parameters used by class (defaults)
    @param alpha (externalRadiationPlanar) ( alpha = 0.0 )
    @param R1 - (externalRadiationPlanar) blr lower boundary radius ('NEW' stratification) (R1 = 0.1*R_sub)
    @param R2 - (externalRadiationPlanar) blr upper boundary radius ('NEW' stratification) (R2 = R_sub)
    @param s - (externalRadiationPlanar) stratification index ('NEW' stratification) (s = 2.0)
    @param e - energy (ext) in eV (e = 10.0)
    @param CF - covering factor (CF = 0.1)
    @param q1 - lower index in stratification ('M03' stratification) (q1 = 1.0)
    @param q2 - upper index in stratification ('M03' stratification) (q2 = 1.0)
    @param Rext - blr upper radius ('M03' stratification) (R2 = 0.1*R_sub)
    @param strat - sets stratification ('NEW' or 'M03') (strat = 'NEW') */

class BlrPlanar : public externalRadiationPlanar {
  /*requested parameters */
  double e, q1, q2, cf, Rext; 
  
  /* own parameters */
  double vext;
  std::string strat;

 public:
  /** constructor
      @param _cfg - scfgp class object
      @param _r - jetGeometry class object
      @param _ele - electrons class object
      @param _id */
  BlrPlanar( scfgp* _cfg, jetGeometry* _r, electrons* _ele, std::string _id );
  
  /** destructor */
  ~BlrPlanar( );

  void printInfo( ); // prints info on class parameters
  void setRadius( ); // sets R if not provided by config file
  double getvext( double _R );
  void setdLdR( );
};

#endif /* _BlrPlanar_hpp */
