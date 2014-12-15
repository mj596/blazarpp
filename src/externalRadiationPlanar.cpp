/**
   @file externalRadiationPlanar.cpp
   @author Mateusz Janiak
*/

#include "externalRadiationPlanar.hpp"
#include "electrons.hpp"

externalRadiationPlanar::externalRadiationPlanar( scfgp* _cfg, jetGeometry* _r, electrons* _ele, std::string _id ) : energyDissProc( _cfg, _r, _ele, _id ) {
  /* requested parameters */
  /* R1, R2, s, alpha - needed by Accd, BlrPl, DtPl */
  cfg -> request<int>("N" + id, 200, &N );
  cfg -> request<double>(id +  "R1", 0.0, &R1 );
  cfg -> request<double>(id +  "R2", 0.0, &R2 );
  cfg -> request<double>(id +  "s", 0.0, &s );
  cfg -> request<double>(id +  "alpha", 0.0, &alpha );
  cfg -> request<double>("Gamma", 10.0, &Gamma );
  cfg -> request<double>("thetaObs", 0.1, &thetaObs );
  cfg -> request<int>("KN", 0, &KN );
  cfg -> request<double>("mBH",1.0,&mBH);
  cfg -> request<double>("eDisk",0.1,&eDisk);
  cfg -> request<double>("mDot",1.0,&mDot);
  cfg -> request<int>(id+"LuminosityConstU",0,&luminosityConstU);  
  cfg -> request<int>(id+"LuminosityConstNu",0,&luminosityConstNu);  

  /* works for both stratifications and for Accd */
  cfg -> request<int>("extPlApp",0,&approx);
  cfg -> request<int>("saveRmVsR",0,&saveRm);
  
  cfg -> updateRequests( );
  
  allocatedLdR( );
  allocateRm( ); }

externalRadiationPlanar::~externalRadiationPlanar( ) {
  freeLpv( );
  freeLvPoint( );
  freeLvPointAvg( );
  freeUpe( );
  freeRm( );
  
  delete R;
  R = NULL; }

/* this function basically calculates Upe before going further */
void externalRadiationPlanar::update(  ) {
  /* even though we use full version we need to calculate the Rm for gamma_dot */
  double oldUpe, currentUpe = 0.0;
  for( int i=0;i<R->getMaxIndex();i++ ) {
    R -> update( i );	  
    oldUpe = currentUpe;
    currentUpe = calculateUpeApprox( i, false );
    set_upe( i, currentUpe );
    if( currentUpe > oldUpe ) {
      Rm = R->get( );
      Rm_index = i; }
    
    if( Rm < R1 ) { Rm = R1; Rm_index = 0; }
    if( Rm > R2 ) { Rm = R2; Rm_index = R -> getMaxIndex( ); }
  }
  
  setRm( Rm );
  bazinga::print_info(id,"Radius Rm", getRm( ) );
  
  if( approx ) { Upe = calculateUpeApprox( Rm_index, true ); }  /* approximate version */
  else { Upe = calculateUpeFull( ); }  /* full integral version */
  bazinga::print_info(id,"Upe",Upe);
  
  set_upe_r( Upe ); /* this is to make a plot gamma_dot vs r */
  flag_upe_r = false; }

double externalRadiationPlanar::calculateUpeFull( ) {
  double _Int = 0.0;
  gsl_vector *_temp_function_upe_vector = gsl_vector_alloc( N );
  double _function_upe_value = 0;

  for( int i=0;i<R->getMaxIndex();i++ ) {
    gsl_vector_set( _temp_function_upe_vector, i, function_upe( i ) );
    _Int += bazinga::IntCor( i, R->getN( ) )*gsl_vector_get( _temp_function_upe_vector, i)*R->getDr(i);
  }
  
  /* before quit save _temp_function_upe_vector */
  bazinga::save_GSLVector( "dUpedR_"+this->whoAmI( ), R->getRadius_GSLVector( ), _temp_function_upe_vector, R->getPosition( ), cfg->get<std::string>("output") );
  gsl_vector_free( _temp_function_upe_vector );
  _temp_function_upe_vector = NULL;
  
  return _Int; }

double externalRadiationPlanar::calculateUpeApprox( int index, bool flag ) { return function_upe( index )*R->get(index); }

double externalRadiationPlanar::function_upe( int index ) {
  double val;
  double r2R2 = pow(r->get(),2.0)+pow(R->get(index),2);
  val = getdLdR(index)*pow(1.0-((beta(Gamma)*r->get())/sqrt(r2R2)),2)/r2R2;
  val *= pow(Gamma,2)/(4.0*M_PI*LIGHT_SPEED);
  return val; }

double externalRadiationPlanar::ctau( ) {
  double val = 0.0;
  if( s == 1 ) { val = 1.0/log(R2/R1); }
  else { val = (1.0-s)/(pow(R2,1.0-s)-pow(R1,1.0-s)); }
  return val; }

double externalRadiationPlanar::getdLdR( int index ) { return gsl_vector_get( dLdR, index ); }

double externalRadiationPlanar::dotg( double g ) {
  double val = 0.0;
  double b = 0.0;
  double r2R2 = pow(r->get(),2.0)+pow(Rm,2.0);
  
  double doppler_ext = 1.0/( Gamma*( 1.0-beta(Gamma)*(r->get()/sqrt(r2R2)) ) );  
  /* as e_ext max we take e_ext @ Rm where it is maximum; in case of BLRPL it makes no difference at all */
  b = 4.0*g*getvext( Rm )*PLANCK_H/(mec2*doppler_ext);

  if( b > 1.0 ) { set_KN_info( g ); }
  val = Upe*inverseCompton::FKN( b, alpha, KN ); /* fKN not used anymore; using FKN instead */

  return val; }

void externalRadiationPlanar::setLpv( ) { for( int i=0;i<N;i++ ) set_Lpv( i, calculateLpv( get_vp(i), thetaObs )); }

double externalRadiationPlanar::fsc_integral( double v, double theta, double R ) {
  double Int = 0.0;
  double ee = PLANCK_H*v/(ELECTRON_MASS*LIGHT_SPEED*LIGHT_SPEED);
  double miu = -(cos(theta)-beta(Gamma))/(1.0-beta(Gamma)*cos(theta));

  double r2R2, doppler_ext;
  if( luminosityConstU ) {
    r2R2 = pow(r->getRInjMax( ),2.0)+pow(R,2.0);
    doppler_ext = 1.0/( Gamma*( 1.0-beta(Gamma)*(r->getRInjMax( )/sqrt(r2R2)) ) ); }
  else {    
    r2R2 = pow(r->get(),2.0)+pow(R,2.0);
    doppler_ext = 1.0/( Gamma*( 1.0-beta(Gamma)*(r->get()/sqrt(r2R2)) ) ); }

  /* integration loop over electrons */
  for( int k=0;k<ele->getN( );k++ ) { Int += bazinga::IntCor( k, ele->getN( ) )*inverseCompton::f(ele->getGamma(k), PLANCK_H*getvext( R )/(mec2*doppler_ext), ee, miu )*ele->getNgamma(k)/ele->getGamma(k); }
  
  Int *= ele->getdLogGamma( );
  return Int; }

double externalRadiationPlanar::function_lum( int index, double v, double theta ) {
  double r2R2, doppler_ext;
  if( luminosityConstU ) {
    r2R2 = pow(r->getRInjMax( ),2.0)+pow(R->get(index),2.0);
    doppler_ext = 1.0/( Gamma*( 1.0-beta(Gamma)*(r->getRInjMax( )/sqrt(r2R2)) ) ); }
  else {    
     r2R2 = pow(r->get(),2.0)+pow(R->get(index),2.0);
    doppler_ext = 1.0/( Gamma*( 1.0-beta(Gamma)*(r->get()/sqrt(r2R2)) ) ); }

  return getdLdR(index)*fsc_integral( v, theta, R->get(index) )/(r2R2*pow( getvext( R->get(index) ), 2 )); }

double externalRadiationPlanar::calculateLpv( double v, double theta ) {
  double val = 0.0;
  if( approx ) { val = function_lum( Rm_index, v, theta )*Rm; }
  else {
    for( int i=0;i<R->getMaxIndex();i++ ) {
      R -> update( i );
      val += bazinga::IntCor( i, R->getN( ) )*function_lum( i, v, theta )*R->getDr(i); }
  }

  val *= 3.0*SIGMA_T/(16.0*M_PI);
  val *= v;
  return val; }

void externalRadiationPlanar::allocateRm( ) {
  rm = gsl_vector_alloc(r->getMaxIndex());
  gsl_vector_set_zero( rm ); }
 
void externalRadiationPlanar::allocatedLdR( ) {
  dLdR = gsl_vector_alloc( N );
  gsl_vector_set_zero( dLdR ); }

void externalRadiationPlanar::freeRm( ) { gsl_vector_free( rm ); }

void externalRadiationPlanar::saveRmVsR( ) {
 if( saveRm ) { bazinga::save_GSLVector( "Rm_"+ this->whoAmI( ), r->getRadius_GSLVector( ), getRmVector( ), cfg->get<std::string>("output") ); }
}

