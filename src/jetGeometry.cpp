#include "jetGeometry.hpp"

jetGeometry::jetGeometry( scfgp* _cfg ): cfg(_cfg) {

  cfg->request<double>("R0",1.0e17,&r0);
  cfg->request<double>("RInjMax",2.0e17,&rinjmax);
  cfg->request<double>("RMax",4.0e17,&rmax);
  cfg->request<std::string>("Rscale","linear",&Rscale);
  cfg->request<int>("Ninj",10,&Ninj);
  cfg->request<int>("nSave",100,&nSave);
    
  cfg->updateRequests( );

  if( Ninj < 2 ) {
    bazinga::warning("jet","Can't use Ninj<2; changing to Ninj=2");
    Ninj = 2; }

  currentPosition = 1.0;

  if( Rscale == "linear" ) {
    N = (rmax-r0)/(rinjmax-r0)*Ninj;
    dr = (rmax-r0)/((double)N);
    max_index = N+1;
    rad = gsl_vector_alloc( max_index );
    for( int k=0;k<max_index;k++ ) gsl_vector_set(rad,k,r0+k*dr); }
  
  if( Rscale == "log" ) {
    N = Ninj*log(rmax/r0)/log(rinjmax/r0);
    dr = -1;
    max_index = N+1;
    rad = gsl_vector_alloc( max_index );
    for (int k=0;k<max_index;k++) { gsl_vector_set(rad,k,r0*pow(rinjmax/r0,(double)k/Ninj)); }
    dr = gsl_vector_get(rad,1)/gsl_vector_get(rad,0); }
}

jetGeometry::~jetGeometry() {
  gsl_vector_free(rad); }

void jetGeometry::printInfo( ) {
  bazinga::info("jet","Setting jet geometry");
  bazinga::print_info("jet","Scale",Rscale);
  bazinga::print_info("jet","R0",r0,"cm");
  bazinga::print_info("jet","RInjMax",rinjmax,"cm");
  bazinga::print_info("jet","RMax",rmax,"cm");
  if( Rscale == "linear" ) { bazinga::print_info("jet","dr",dr,"cm"); }
  bazinga::print_info("jet","Number of cells in diss zone",Ninj);
  bazinga::print_info("jet","Number of cells in total",N);
  bazinga::print_info("jet","nSave",nSave); }

void jetGeometry::update( int index ) {
  current_index = index;
  current_r = gsl_vector_get(rad,index);
  currentPosition = current_r/r0; }

double jetGeometry::show( ) {
  std::ostringstream s;
  s << " (" << this->get() << " cm)";
  bazinga::print_header( );
  bazinga::print_info("jet","Jet radius",this->get()/(double)r0,s.str());
  bazinga::print_header( ); }

double jetGeometry::get( ) { return current_r; }
double jetGeometry::get( int index ) { return gsl_vector_get(rad,index); }
int jetGeometry::getIndex( ) { return current_index; }
int jetGeometry::getMaxIndex( ) { return max_index; }
double jetGeometry::getPosition( ) { return currentPosition; }
double jetGeometry::getPosition( int index ) { return gsl_vector_get(rad,index)/r0; }
double jetGeometry::getDr( ) { 
  if( Rscale == "linear" ) { return dr; }
  if( Rscale == "log" ) { return gsl_vector_get(rad,current_index)*(dr-1.0); }
}
double jetGeometry::getDr( int index ) { return gsl_vector_get(rad,index)*(dr-1.0); }
int jetGeometry::getNinj( ) { return Ninj; }
double jetGeometry::getR0( ) { return r0; }
double jetGeometry::getRInjMax( ) { return rinjmax; }
double jetGeometry::getRMax( ) { return rmax; }

int jetGeometry::ifSaveRadius( ) {
  if( nSave <= 0 ) { return 1; }
  int plot = 0;
  double drSave = (this->getMaxIndex( )-2)/((double)nSave-2);

  if( this->getMaxIndex( ) < nSave ) { plot = 1; }
  else {
    if( this->getIndex( ) == 0 || this->getIndex( ) == this->getMaxIndex( )-1 ) { plot = 1; }
    for( int i=1; i<nSave-1; i++ )
      if( this->getIndex( ) == floor(i*drSave) ) { plot = 1; }
  }
  return plot; }
