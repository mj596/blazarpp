/**
   @file DtSphericalGu.cpp
   @author Mateusz Janiak
*/

#include "DtSphericalGu.hpp"
#include "electrons.hpp"

DtSphericalGu::DtSphericalGu( scfgp* cfg, jetGeometry* r, electrons* ele, std::string id ) : externalRadiationGu( cfg, r, ele, id ) {
  /* requested parameters */
  cfg -> request<double>(id + "e", 0.6203, &e );
  cfg -> request<double>(id + "cf", 0.1, &cf );
  cfg -> request<double>(id + "rext", 0.0, &rext );
  cfg -> request<double>(id + "k", 2.0, &k );

  cfg -> updateRequests( );

  if( rext == 0.0 ) { setRadius( ); }
  else { bazinga::warning(id,"Overwritting rext value from config file!"); }
  
  vpMin = 0.01*DSQR( ele -> getGammaMin( ) )*getvext( r->getRMax() );
  vpMax = 100.0*DSQR( ele -> getGammaMax( ) )*getvext( r->getR0() );

  allocateUpeR( ); // this is for storing vector information about Upe(r) 
  allocateLpv( );
  allocateLvPoint( );
  allocateLvPointAvg( );

  for( int i=0;i<N;i++ ) { set_vp( i, vpMin*pow( vpMax/vpMin,(double)i/((double)N-1)) ); }
}

double DtSphericalGu::getvext( double _r ) {
  double vext;
  vext = e*eV2erg/PLANCK_H;
  if( luminosityConstNu ) { return vext; }
  else { return vext*rext/_r; }
}

void DtSphericalGu::update(  ) {
  double Ledd = 1.3e47*mBH;
  double Ldisk = eDisk*mDot*Ledd;
  double temp_upe = 0.0;
  double temp_cf = 0.0;
  temp_cf = cf;
  temp_cf *= pow(r->get()/rext,2)/(1.0+pow(r->get()/rext,k));
  temp_upe = (temp_cf*Ldisk*DSQR(Gamma))/(3.0*M_PI*DSQR(r->get())*LIGHT_SPEED);
  gsl_vector_set( upe_r, r->getIndex( ), temp_upe );
  bazinga::print_info(id,"Upe",gsl_vector_get( upe_r, r->getIndex( ) ) );
  set_upe_r( gsl_vector_get( upe_r, r->getIndex( ) ) );
  flag_upe_r = false; }

void DtSphericalGu::printInfo( ) {
  bazinga::info(id,"Info");
  bazinga::print_info(id,"N",N);
  bazinga::print_info(id,"Avg energy (in external frame)",e,"eV");
  bazinga::print_info(id,"Avg frequency (in external frame)",getvext( r->getR0()),"Hz");
  bazinga::print_info(id,"Clouds covering factor cf",cf);
  bazinga::print_info(id,"radius rext",rext,"cm");
  bazinga::print_info(id,"strat index k",k);
  bazinga::print_info(id,"vpMin",vpMin,"Hz");
  bazinga::print_info(id,"vpMax",vpMax,"Hz");
}

void DtSphericalGu ::setRadius( ) {
  double Ledd = 1.3e47*mBH;
  double Ldisk = eDisk*mDot*Ledd;
  double R_sub = 1.6e-5*sqrt( Ldisk );
  rext = R_sub; }

