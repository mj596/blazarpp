/**
    @file DtSpherical.cpp
    @author Mateusz Janiak
*/

#include "DtSpherical.hpp"
#include "electrons.hpp"

DtSpherical::DtSpherical( scfgp* cfg, jetGeometry* r, electrons* ele, std::string id ) : externalRadiationSpherical( cfg, r, ele, id ) {
  /* request parameters */
  cfg -> request<double>( id+"q1", 1.0, &q1 );
  cfg -> request<double>( id+"q2", 1.0, &q2 );
  cfg -> request<double>( id+"r", 0.0, &rext );
  cfg -> request<double>("mBH",1.0,&mBH);
  cfg -> request<double>("eDisk",0.1,&eDisk);
  cfg -> request<double>("mDot",1.0,&mDot);
  cfg -> request<double>("Gamma",10.0,&Gamma);

  cfg -> updateRequests( );

  if( rext == 0.0 ) { setRadius( ); }
  else { bazinga::warning(id,"Overwritting rext value from config file!"); } 

  vpMin = 0.01*DSQR( ele -> getGammaMin() )*getvext( 10.*rext );
  vpMax = 100.0*DSQR( ele -> getGammaMax() )*getvext( rext );

  allocateUpe( );
  allocateLpv( );
  allocateLvPoint( );
  allocateLvPointAvg( );

  for( int i=0;i<N;i++ ) { set_vp( i, vpMin*pow( vpMax/vpMin,(double)i/((double)N-1)) ); }
}

void DtSpherical::printInfo( ) {
  bazinga::info(id,"Info");
  bazinga::print_info(id,"id",id);
  bazinga::print_info(id,"N",N);
  bazinga::print_info(id,"Avg energy (in external frame) @ rext",e,"eV");
  bazinga::print_info(id,"Avg frequency (in external frame) @ rext",getvext( rext ),"Hz");
  bazinga::print_info(id,"Clouds covering factor cf",cf);
  bazinga::print_info(id,"radius r",rext,"cm");
  //  bazinga::print_info(id,"radius r2",r2,"cm");
  bazinga::print_info(id,"q1",q1);
  bazinga::print_info(id,"q2",q2);
  bazinga::print_info(id,"vpMin",vpMin,"Hz");
  bazinga::print_info(id,"vpMax",vpMax,"Hz");
  if( luminosityConstU ) { bazinga::warning(id,"Using constant u' to calculate luminosity."); }
  if( luminosityConstNu ) { bazinga::warning(id,"Using constant v' to calculate luminosity."); }
 }

void DtSpherical::setRadius( ) {
  double Ledd = 1.3e47*mBH;
  double Ldisk = eDisk*mDot*Ledd;
  double R_sub = 1.6e-5*sqrt(Ldisk);
  rext = R_sub; }
/* obsolete
   r2 = 10.0*R_sub; } */

/* now obsolete
double DtSpherical::getdLdlnr( double _r ) {
  double val = 0.0;
  double Ledd = 1.3e47*mBH;
  double Ldisk = eDisk*mDot*Ledd;

  val = cf*Ldisk/((q1+q2)/(q1*q2));
  if( _r <= rext ) { val *= pow( _r/rext,q1); }
  if( _r > rext ) { val *= pow( _r/rext,-q2); }
  return val; }

  //  val = cf*Ldisk/log(r2/rext);
  //  return val; } */

double DtSpherical::getvext( double _r ) {
  vext = e*eV2erg/PLANCK_H;
  if( luminosityConstNu ) { return vext; }
  else { return vext*rext/_r; }
}

void DtSpherical::update(  ) {
  double Ledd = 1.3e47*mBH;
  double Ldisk = eDisk*mDot*Ledd;
  double u_ext = Gamma*Gamma*cf*Ldisk/(3.*M_PI*LIGHT_SPEED*rext*rext);
  if( r->get( ) <= rext ) { u_ext *= pow( r->get( )/rext, -q1 ); }
  if( r->get( ) > rext ) { u_ext *= pow( r->get( )/rext, -q2 ); }
  gsl_vector_set( upe, r->getIndex( ), u_ext );
  bazinga::print_info(id,"Upe",gsl_vector_get( upe, r->getIndex( ) ) );
  set_upe_r( gsl_vector_get( upe, r->getIndex( ) ) );
  flag_upe_r = false; }

/* obsolete
void DtSpherical::update(  ) {
  if( r->get( ) <= rext ) { gsl_vector_set( upe, r->getIndex( ), (4.0*DSQR(Gamma)*getdLdlnr( rext ))/(3.0*4.0*M_PI*DSQR(rext)*LIGHT_SPEED) ); }
  else { gsl_vector_set( upe, r->getIndex( ), (4.0*DSQR(Gamma)*getdLdlnr( r->get() ))/(3.0*4.0*M_PI*DSQR(r->get())*LIGHT_SPEED) ); }
  bazinga::print_info(id,"Upe",gsl_vector_get( upe, r->getIndex( ) ) );
  set_upe_r( gsl_vector_get( upe, r->getIndex( ) ) );
  flag_upe_r = false; } */
