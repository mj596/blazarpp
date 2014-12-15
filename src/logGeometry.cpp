#include "logGeometry.hpp"

logGeometry::logGeometry( double _r1, double _r2, double N ): r0(_r1), rmax(_r2), max_index(N+1) {
  rad = gsl_vector_alloc(max_index);
  for( int k=0;k<max_index;k++ ) { gsl_vector_set(rad,k,r0*pow(rmax/r0,(double)k/N)); }
  dr = gsl_vector_get(rad,1)/gsl_vector_get(rad,0); }

logGeometry::~logGeometry() {
  gsl_vector_free(rad); }

void logGeometry::update( int index ) {
  current_index = index;
  current_r = gsl_vector_get(rad,index);
  currentPosition = current_r/r0; }

double logGeometry::show( ) {
  std::ostringstream s;
  s << " (" << this->get() << " cm)";
  bazinga::print_info("+ Log radius",this->get()/(double)r0,s.str()); }

double logGeometry::get( ) { return current_r; }
double logGeometry::get( int index ) { return gsl_vector_get(rad,index); }
int logGeometry::getIndex( ) { return current_index; }
int logGeometry::getMaxIndex( ) { return max_index; }
double logGeometry::getPosition( ) { return currentPosition; }
double logGeometry::getPosition( int index ) { return gsl_vector_get(rad,index)/r0; }
double logGeometry::getDr( ) { return gsl_vector_get(rad,current_index)*(dr-1.0); }
double logGeometry::getDr( int index ) { return gsl_vector_get(rad,index)*(dr-1.0); }
double logGeometry::getR0( ) { return r0; }
double logGeometry::getRMax( ) { return rmax; }
