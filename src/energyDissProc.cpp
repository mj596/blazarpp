 /**
    @file energyDissProc.cpp
    @author Mateusz Janiak
*/

#include "energyDissProc.hpp"
#include "electrons.hpp"

energyDissProc::energyDissProc( scfgp* _cfg, jetGeometry* _r, electrons* _ele, std::string _id ) : baseClass(_cfg, _r, _id ), ele( _ele ), flag_upe_r( false ) {
  
  /* request parameters */
  cfg -> request<int>("saveLum",0,&saveLum);
  cfg -> request<int>("saveUpeVsR",0,&saveUpeVsR);
  cfg -> request<double>("injRm",2.0e17,&injRm);
  
  cfg -> updateRequests( );

  /* allocate common memory */
  allocateGammaKN(); }

energyDissProc::~energyDissProc( ) {
  freeGammaKN( );
  ele = NULL; }

void energyDissProc::allocateGammaKN( ) {
  gammaKN = gsl_vector_alloc(r->getMaxIndex());
  bazinga::print_GSLVector_allocated_memory( id, gammaKN );
  gsl_vector_set_zero( gammaKN ); }

void energyDissProc::allocateLpv( ) {
  vp = gsl_vector_alloc( N );
  Lpv = gsl_vector_alloc( N );
  bazinga::print_GSLVector_allocated_memory( id, vp );
  bazinga::print_GSLVector_allocated_memory( id, Lpv );
  gsl_vector_set_zero( vp );
  gsl_vector_set_zero( Lpv ); }

void energyDissProc::allocateLvPoint( ) {
  vPoint = gsl_vector_alloc( N );
  LvPoint = gsl_vector_alloc( N );
  bazinga::print_GSLVector_allocated_memory( id, vPoint );
  bazinga::print_GSLVector_allocated_memory( id, LvPoint );
  gsl_vector_set_zero( vPoint );
  gsl_vector_set_zero( LvPoint ); }

void energyDissProc::allocateLvPointAvg( ) {
  LvPointAvg = gsl_vector_alloc( N );
  bazinga::print_GSLVector_allocated_memory( id, LvPointAvg );
  gsl_vector_set_zero( LvPointAvg ); }

void energyDissProc::allocateUpe( ) {
  ep = gsl_vector_alloc( N );
  upe = gsl_matrix_alloc( r->getMaxIndex(), N );
  bazinga::print_GSLVector_allocated_memory( id, ep );
  bazinga::print_GSLMatrix_allocated_memory( id, upe );  
  gsl_vector_set_zero( ep );
  gsl_matrix_set_zero( upe ); }

void energyDissProc::allocateUpeR( ) {
  upe_r = gsl_vector_alloc( r->getMaxIndex() );
  bazinga::print_GSLVector_allocated_memory( id, upe_r );
  gsl_vector_set_zero( upe_r ); }

void energyDissProc::freeLpv( ) {
  gsl_vector_free( vp );
  gsl_vector_free( Lpv ); }

void energyDissProc::freeLvPoint( ) {
  gsl_vector_free( vPoint );
  gsl_vector_free( LvPoint ); }

void energyDissProc::freeLvPointAvg( ) {
  gsl_vector_free( LvPointAvg ); }

void energyDissProc::freeUpe( ) {
  gsl_vector_free( ep );
  gsl_matrix_free( upe ); }

void energyDissProc::freeUpeR( ) {
  gsl_vector_free( upe_r ); }

void energyDissProc::freeGammaKN( ) {
  gsl_vector_free( gammaKN ); }

void energyDissProc::set_upe_r( double val ) {
  if( !flag_upe_r ) {
      flag_upe_r = true;
      gsl_vector_set( upe_r, r->getIndex( ), val ); }
}

void energyDissProc::set_KN_info( double g ) {
  if( get_gammaKN( ) == 0 ) { set_gammaKN( g ); }
  else if( get_gammaKN( ) > g ) { set_gammaKN( g ); }
}

void energyDissProc::print_KN_info( ) { 
  std::cout << " KN info "+this->whoAmI() << "\tradius: " << r->get() << "\tgammaKN: ";
  if( get_gammaKN() > ele->getGammaMax( ) || get_gammaKN( ) == 0.0 ) {
    std::cout << "> gammaMax = " << std::scientific << std::setprecision(2) << ele->getGammaMax( ) << std::endl; }
  else { std::cout << std::scientific << std::setprecision(2) << get_gammaKN() << std::endl; }
}

void energyDissProc::saveLuminosity( ) {
  std::string type = "Lpv_"+this->whoAmI( );
  bazinga::info(this->whoAmI( ),"Saving luminosity.");
  if( saveLum && r->ifSaveRadius( ) ) { bazinga::save_GSLVector( type, this->vp, ele->gamma, this->Lpv, r->getPosition( ), cfg->get<std::string>("output") ); }
}

void energyDissProc::saveUpeR( ) {
  std::string type = "UpeR_"+this->whoAmI( );
  bazinga::info(this->whoAmI( ),"Saving upe vs radius info.");
  if( saveUpeVsR ) { bazinga::save_GSLVector( type, r->getRadius_GSLVector( ), this->upe_r, cfg->get<std::string>("output") ); }
}
