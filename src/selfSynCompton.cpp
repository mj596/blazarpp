/**
   @file selfSynCompton.cpp
   @author Mateusz Janiak
*/

#include "selfSynCompton.hpp"
#include "electrons.hpp"

selfSynCompton::selfSynCompton( scfgp* _cfg, jetGeometry* _r, electrons* _ele, magneticField* _B, std::string _id ) : energyDissProc( _cfg, _r, _ele, _id  ), B(_B) {
  /* requested parameters */
  cfg -> request<int>("N"+id,200,&N);
  cfg -> request<int>("KN",0,&KN);
  
  cfg -> updateRequests( );

  vpMin = 1.0;
  /* here we set vpMAX to a maximum value specified by maximal value of magnetic field at R0 */
  vpMax = 100.0*A43*DSQR(ele->getGammaMax( ))*DSQR(ele->getGammaMax( ))*( B->getMaxB( )/B_CR )*mec2h;

  allocateUpe( );
  allocateUpeR( );
  allocateLpv( );
  allocateLvPoint( );
  allocateLvPointAvg( );

  for( int i=0;i<N;i++ ) { set_vp( i, vpMin*pow( vpMax/vpMin,(double)i/((double)N-1)) ); }
 }

selfSynCompton::~selfSynCompton( ) {
  freeLpv( );
  freeLvPoint( );
  freeLvPointAvg( );
  freeUpe( ); }

void selfSynCompton::printInfo( ) {
  bazinga::info(id,"Info");
  bazinga::print_info(id,"N",N);
  bazinga::print_info(id,"KN",KN); }

void selfSynCompton::update( ) { flag_upe_r = false; }

double selfSynCompton::dotg( double g ) {
  double b,sum, val = 0.0;

  for( int i=0;i<syn->getN( );i++ ) {
    b = 4.0*syn->get_ep(i)*g;
    if( b > 1.0 ) { set_KN_info( g ); }
    sum += bazinga::IntCor( i, syn->getN( ) )*syn->get_ep(i)*syn->get_upe(i)*inverseCompton::fKN( b, KN ); }
  
  val = syn->getdLogE( )*sum;
  
  set_upe_r( val );  
  return( val ); }

void selfSynCompton::setLpv( ) {
  for( int i=0;i<N;i++ ) { set_Lpv( i, calculateLpv( get_vp(i) ) ); } 
}

double selfSynCompton::calculateLpv( double v ) {
  double dg,de,j;
  double Int1;
  double e;
  Int1 = 0.0;
  j = 0.0;

  e = PLANCK_H*v/(ELECTRON_MASS*LIGHT_SPEED*LIGHT_SPEED);
  
  dg = ele->getdLogGamma( );
  de = syn->getdLogE( );

  for( int i=0; i<syn->getN( );i++ ) {
    Int1 = 0.0;
    /* integration over electron energy (gamma) */
    for( int k=0;k<ele->getN( );k++ ) { 
      Int1  += bazinga::IntCor( k, ele->getN( ) )*inverseCompton::fiso( ele->getGamma(k), syn->get_ep(i),e )*ele->getNgamma(k)/ele->getGamma(k); }
    
      /* integration over seed photon energy */
    j += bazinga::IntCor( i, syn->getN( ) )*syn->get_upe(i)/syn->get_ep(i)*Int1*dg; }

  j *= constC1*e*de;
  return j; }
