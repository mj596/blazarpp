#ifndef _logGeometry_hpp
#define _logGeometry_hpp 1

#include <bazinga.hpp>
#include <cmath>

class logGeometry {
  double r0, rmax;
  int N;
  gsl_vector* rad;
  int max_index;
  double current_r;
  int current_index;
  double dr;
  double currentPosition;
  
 public:
  logGeometry( double _r1, double _r2, double _N );
  ~logGeometry();

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
  double getR0( );
  double getRMax( );
  int getN( ) { return N; }
  gsl_vector* getRadius_GSLVector( ) { return rad; }
};

#endif /* _logGeometry_hpp */
