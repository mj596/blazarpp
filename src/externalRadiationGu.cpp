/**
   @file externalRadiationGu.cpp
   @author Mateusz Janiak
*/

#include "externalRadiationGu.hpp"
#include "electrons.hpp"

externalRadiationGu::externalRadiationGu( scfgp* _cfg, jetGeometry* _r, electrons* _ele, std::string _id ) : energyDissProc( _cfg, _r, _ele, _id ) {
  /* requested parameters */
  cfg -> request<int>("N" + id, 200, &N );
  cfg -> request<double>("Gamma", 10.0, &Gamma );
  cfg -> request<double>("thetaObs", 0.1, &thetaObs );
  cfg -> request<int>("KN", 0, &KN );
  cfg -> request<double>("mBH",1.0,&mBH);
  cfg -> request<double>("eDisk",0.1,&eDisk);
  cfg -> request<double>("mDot",1.0,&mDot);
  cfg -> request<int>(id+"LuminosityConstU",0,&luminosityConstU);  
  cfg -> request<int>(id+"LuminosityConstNu",0,&luminosityConstNu);  

  cfg -> updateRequests( ); }

externalRadiationGu::~externalRadiationGu( ) {
  freeLpv( );
  freeLvPoint( );
  freeLvPointAvg( );
  freeUpeR( ); }

double externalRadiationGu::dotg( double g ) {
  double b, val = 0.0;
  b = 4.0*g*getvext( r->get() )*PLANCK_H*Gamma/mec2;
  if( b > 1.0 ) { set_KN_info( g ); }
  val = gsl_vector_get( upe_r, r->getIndex( ) )*inverseCompton::fKN(b,KN); 
  return val; }

double externalRadiationGu::calculateLpv( double v, double theta ) {
  double Int = 0.0;
  double e, miu, ep;

  e = PLANCK_H*v/mec2;
  miu = -( cos(theta)-beta(Gamma) )/( 1.0-beta(Gamma)*cos(theta) );
  ep = getvext( r->get() )*PLANCK_H*Gamma/mec2;

  for( int k=0;k<ele->getN( );k++ )
    Int  += bazinga::IntCor( k, ele->getN( ) )*inverseCompton::f( ele->getGamma(k), ep, e, miu )*ele->getNgamma(k)/ele->getGamma(k);

  if( luminosityConstU )
    { Int *= gsl_vector_get( upe_r, 0 )*ele->getdLogGamma( )/pow( getvext( r->get() ), 2.0 ); }
  else 
    { Int *= gsl_vector_get( upe_r, r->getIndex( ) )*ele->getdLogGamma( )/pow( getvext( r->get( ) ), 2.0 ); }

  Int *= (3.0*SIGMA_T*LIGHT_SPEED)/(4.0*pow(Gamma,2.0));
  Int *= v;
  return Int; }

void externalRadiationGu::setLpv( ) {
  for( int i=0;i<N;i++ ) { set_Lpv( i, calculateLpv( get_vp(i),thetaObs ) ); }
}
