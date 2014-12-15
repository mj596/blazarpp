 /**
    @file QuasarAccDisk.hpp
    @author Mateusz Janiak
*/

#ifndef _QuasarAccDisk_hpp
#define _QuasarAccDisk_hpp 1

#include "baseClass.hpp"
#include "logGeometry.hpp"

/** @class QuasarAccDisk
    Class used to process quasar radiation template and add it to calculated jet radiation;
    at the moment quasar radiation is not stricty calculated but a radius-loud quasar radiation template from Elvis et al. (1994) is used.
    Class has some basic methods to calculate standard multitemperature Shakura-Sunyayev disk radiation but it is not used currently as it need 
    much, much detailed study and modelling of dusty torus, coronae etc */

class QuasarAccDisk : public baseClass {
  /* requested parameters */
  double mBH, mDot, eDisk;
  int N;
  double R1, R2;
  int save_Lr;
  std::string templateFilename;
  
  /* other parameters */
  double Rg;
  gsl_vector* v; // frequency
  gsl_vector* Lv; // monochromatic apparent luminosity
  double vMin, vMax;
  int numLines;
  double Ledd, Ldisk;

  double loglogIntegrate( gsl_vector* x, gsl_vector* y );

public:
  /** constructor
     @param scfgp
     @param id */
  QuasarAccDisk( scfgp* _cfg, std::string _id );
  
  /** destructor */
  ~QuasarAccDisk( );

  /** vectors to store data */
  gsl_vector* vTemplate;
  gsl_vector* LvTemplate;
  
  /** logGeometry to store log-spaced data for loglog-interpolation */
  logGeometry* R;

  void printInfo( );
  
  /** get BB temperature at distance r (SS disk)
      @param radius r */
  double getTemperature( double _r );
  
  /** calculate disk luminosity
      @param frequency v
      @param radius_index */
  double calculateLv( double v, int radius_index );
  
  /** calculate a whole disk luminosity */
  void calculateLuminosity( );

  /** set accretion disk radial boundaries */
  void setRadius( );

  /** set energetics for further calculations */
  void setEnergetics( );

  /** read quasar radiation template
      @returns 1 if success; 0 otherwise */
  int readTemplate( );

  /** scale quasar radiation template to match set accretion disk efficiency */
  void scaleTemplate( );

  /** save scaled quasar radiation template */
  int saveTemplate( );

  /** methods to allocate and free memory */
  void allocateLv( );
  void freeLv( );

  void allocateLvTemplate( );
  void freeLvTemplate( );
  
  void set_v( int i, double val ) { gsl_vector_set( v, i, val ); }
  void set_Lv( int i, double val ) { gsl_vector_set( Lv, i, val ); }
  int getN( ) { return N; }

  gsl_vector* get_v( ) { return v; }
  gsl_vector* get_Lv( ) { return Lv; }

  int get_save_Lr( ) { return save_Lr; }
};

#endif /* _QuasarAccDisk_H */
