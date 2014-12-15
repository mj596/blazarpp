#ifndef _jetGeometry_hpp
#define _jetGeometry_hpp 1

#include <bazinga.hpp>
#include <scfgp.hpp>
#include <cmath>

class jetGeometry {
  scfgp* cfg;
  
  /* requested model parameters */
  double r0, rmax, rinjmax;
  std::string Rscale;
  int Ninj, N;
  int nSave;

  /* other */
  gsl_vector* rad;
  int max_index;
  double current_r;
  int current_index;
  double dr;
  double currentPosition;
  
 public:
  jetGeometry( scfgp* _cfg );
  //  jetGeometry( double _r1, double _r2, double _N );
  ~jetGeometry();

  void update( int index );
  double show( );
  double get( );
  double get( int index );
  int getIndex( );
  int getMaxIndex( );
  double getDr( );
  double getDr( int index );
  double getPosition( );
  double getPosition( int index );
  void printInfo( );
  int getNinj( );
  double getR0( );
  double getRInjMax( );
  double getRMax( );
  gsl_vector* getRadius_GSLVector( ) { return rad; }
  int ifSaveRadius( );
};

#endif /* _jetGeometry_hpp */
