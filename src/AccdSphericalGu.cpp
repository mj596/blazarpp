/**
   @file AccdSphericalGu.hpp
   @author Mateusz Janiak
*/

#include "AccdSphericalGu.hpp"
#include "electrons.hpp"

AccdSphericalGu::AccdSphericalGu( scfgp* cfg, jetGeometry* r, electrons* ele, std::string id ) : externalRadiationGu( cfg, r, ele, id ) {
  /* requested parameters */
//  cfg -> request<double>(id + "e", 10.0, &e );
//  cfg -> request<double>(id + "cf", 0.1, &cf );
  cfg -> request<double>(id + "rext", 0.0, &rext );
//  cfg -> request<double>(id + "k", 3.0, &k );
//
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

double AccdSphericalGu::getvext( double _r ) {
  double Ledd = 1.3e47*mBH;
  double Ldisk = eDisk*mDot*Ledd;
  double vext = 0.0;
  double Fdisk;
  if( luminosityConstNu ) { Fdisk = (3.0*G_CONST*mBH*1e09*MSUN)*mDot*Ledd/( DSQR(LIGHT_SPEED)*8.0*M_PI*pow(_r,3.0) ); }
  else { Fdisk = (3.0*G_CONST*mBH*1e09*MSUN)*mDot*Ledd/( DSQR(LIGHT_SPEED)*8.0*M_PI*pow(_r,3.0) ); }
  double Teff = pow( Fdisk/SIGMA_SB, 0.25 );
  vext = 3.92*Teff*K_BOLTZMAN/PLANCK_H;
  return vext; }

void AccdSphericalGu::update(  ) {
  double Ledd = 1.3e47*mBH;
  double Ldisk = eDisk*mDot*Ledd;
  double temp_upe = 0.0;
  double temp_cf = 0.0;
  temp_cf = 0.75*(((0.28*rext)/r->get())+(1.0/(4.0*DSQR(Gamma)*DSQR(Gamma))));
  temp_upe = (temp_cf*Ldisk*DSQR(Gamma))/(3.0*M_PI*DSQR(r->get())*LIGHT_SPEED);
  gsl_vector_set( upe_r, r->getIndex( ), temp_upe );
  bazinga::print_info(id,"Upe",gsl_vector_get( upe_r, r->getIndex( ) ) );
  set_upe_r( gsl_vector_get( upe_r, r->getIndex( ) ) );
  flag_upe_r = false; }

void AccdSphericalGu::printInfo( ) {
  bazinga::info(id,"Info");
  bazinga::print_info(id,"N",N);
  bazinga::print_info(id,"Avg energy (in external frame)",getvext( r->getR0())*PLANCK_H,"eV");
  bazinga::print_info(id,"Avg frequency (in external frame)",getvext( r->getR0()),"Hz");
  bazinga::print_info(id,"radius rext",rext,"cm");
  bazinga::print_info(id,"vpMin",vpMin,"Hz");
  bazinga::print_info(id,"vpMax",vpMax,"Hz");
}

void AccdSphericalGu ::setRadius( ) {
  double Ledd = 1.3e47*mBH;
  double Ldisk = eDisk*mDot*Ledd;
  double R_sub = 1.6e-5*sqrt( Ldisk );
  double Rg = mBH*1e09*MSUN*G_CONST/DSQR(LIGHT_SPEED);
  rext = Rg; }

