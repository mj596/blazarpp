/**
   @file externalRadiation.cpp
   @author Mateusz Janiak
*/

#include "externalRadiation.hpp"
#include "electrons.hpp"

externalRadiation::externalRadiation( scfgp* _cfg, jetGeometry* _r, electrons* _ele, std::string _id ) : energyDissProc( _cfg, _r, _ele, _id ) {
  /* requested parameters */
  // std::string parameter;
  //  parameter = "N" + id;
  cfg -> request<int>("N" + id,200,&N);
  cfg -> request<double>(id+"E",1.0,&exte);
  cfg -> request<double>(id+"U",1.0e-2,&extu);
  cfg -> request<double>(id+"R",0.0,&extr);
  cfg -> request<double>(id+"k",2.0,&extk);
  cfg -> request<double>("Gamma",10.0,&Gamma);
  cfg -> request<int>("KN",0,&KN);
  cfg -> request<double>("thetaObs",0.1,&thetaObs);  

  cfg -> updateRequests( );  

  exte *= Gamma;
  extu *= 4.0/3.0*Gamma*Gamma;

  allocateUpe( );
  allocateUpeR( );
  allocateLpv( );
  allocateLvPoint( );
  allocateLvPointAvg( );

  /* set parameters for calculating BB radiation */
  TempX = constA1*exte;
  normA = constA2*extu/TempX;
  epMin = (1.0e-6*TempX>PLANCK_H/mec2) ? 1.0e-6*TempX : PLANCK_H/mec2;
  epMax = 100.0*TempX;

  /* set frequency for luminosities */
  vpMin = 0.01*DSQR( ele->getGammaMin( ) )*mec2h*epMin;
  vpMax = 100.1*DSQR( ele->getGammaMax( ) )*mec2h*epMax;
  
  for( int i=0;i<N;i++ ) {
    set_ep( i, epMin*pow(epMax/epMin,(double)i/(double)(N-1)) );
    set_vp( i, vpMin*pow( vpMax/vpMin,(double)i/((double)N-1)) ); }    

  /* set integration step over ep */
  dLogE = log( get_ep(1)/get_ep(0) ); }

externalRadiation::~externalRadiation( ) {
  freeLpv( );
  freeLvPoint( );
  freeLvPointAvg( );
  freeUpe( );
  freeUpeR( ); }

void externalRadiation::printInfo( ) {
  bazinga::info(id,"Info");
  bazinga::print_info(id,"N",N);
  bazinga::print_info(id,"radius R",extr,"cm");
  bazinga::print_info(id,"Index kERC",extk);
  bazinga::info(id,"(in electrons co-moving frame)");
  bazinga::print_info(id,"Avg energy",exte,"eV");
  bazinga::print_info(id,"Energy density",extu,"erg cm-3"); }

void externalRadiation::update(  ) {
  for (int i=0;i<N;i++ ) {
    ksi = get_ep(i)/TempX;
    set_upe( i, normA*pow(ksi,3)/(exp(ksi)-1.0)/(1.0e0+pow(r->get()/extr,extk) ) ); }
  flag_upe_r = false; }

double externalRadiation::dotg( double g ) {
  double b, sum, val = 0.0;
  for( int i=0;i<N;i++ ) {
      b = 4.0*get_ep(i)*g;
      if( b > 1.0 ) { set_KN_info( g ); }
      sum += bazinga::IntCor(i,N)*get_ep(i)*get_upe(i)*inverseCompton::fKN(b,KN); }

  val = dLogE*sum;
  set_upe_r( val );
  return val; }

void externalRadiation::setLpv( ) {
  for( int i=0;i<N;i++ ) { set_Lpv( i, calculateLpv( get_vp(i), thetaObs ) ); }
}

double externalRadiation::calculateLpv( double v, double theta ) {
  double j;
  double Int1;
  int i,k;
  double e,miu;
  
  j = 0.0;
  e = PLANCK_H*v/mec2;
  miu = -( cos(theta)-beta(Gamma) )/( 1.0-beta(Gamma)*cos(theta) );
    
  for( int i=0;i<N;i++ ) {
    Int1 = 0.0;
    for( int k=0;k<ele->N;k++ ) { Int1  += bazinga::IntCor( k, ele->getN() )*inverseCompton::f( ele->getGamma(k), get_ep(i), e, miu )*ele->getNgamma(k)/ele->getGamma(k); }
    j += bazinga::IntCor(i,N)*get_upe(i)/get_ep(i)*Int1*ele->getdLogGamma( );	
  }
  j *= constC1*e*dLogE;
  return j; }
