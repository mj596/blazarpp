/**
   @file AccdPlanar.cpp
   @author Mateusz Janiak
*/

#include "AccdPlanar.hpp"
#include "electrons.hpp"

AccdPlanar::AccdPlanar( scfgp* cfg, jetGeometry* r, electrons* ele, std::string id ) : externalRadiationPlanar( cfg, r, ele, id ) {
  /* request parameters */
  /* -- */

  if( R1 == 0.0 || R2 == 0.0 ) { setRadius( ); }
  else { bazinga::warning(id,"Overwritting R1 and R2 values from config file!"); }
  
  vpMin = 0.01*DSQR( ele -> getGammaMin( ) )*getvext( R2 );
  vpMax = 100.0*DSQR( ele -> getGammaMax( ) )*getvext( R1 );

  /* initialize disk radius R;
     here we use N-1 to use the same Upe as in other processes */
  R = new logGeometry( R1, R2, N-1 );

  allocateUpe( ); // this is for storing matrix information about dUpe(r)/dR (R) 
  allocateUpeR( ); // this is for storing vector information about Upe(r) 
  allocateLpv( );
  allocateLvPoint( );
  allocateLvPointAvg( );

  for( int i=0;i<N;i++ ) { set_vp( i, vpMin*pow( vpMax/vpMin,(double)i/((double)N-1)) ); }

  /* set the values of dLdR */
  setdLdR( ); }

AccdPlanar::~AccdPlanar( ) {
  freeUpe( );
  freeLpv( );
  freeLvPoint( );
  freeLvPointAvg( ); }

void AccdPlanar::printInfo( ) {
  bazinga::info(id,"Info");
  bazinga::print_info(id,"N",N);
  bazinga::info(id,"Using stratification version NEW (only option)");
  bazinga::print_info(id,"radius R1",R1,"cm");
  bazinga::print_info(id,"radius R2",R2,"cm");
  bazinga::print_info(id,"stratification index",s);
  bazinga::print_info(id,"alpha",alpha);
  bazinga::print_info(id,"vpMin",vpMin,"Hz");
  bazinga::print_info(id,"vpMax",vpMax,"Hz");
  if( approx ) { bazinga::info(id,"Using approximate dependence on R"); } 
  if( luminosityConstU ) { bazinga::warning(id,"Using constant u' to calculate luminosity."); }
  if( luminosityConstNu ) { bazinga::warning(id,"Using constant v_ext to calculate luminosity."); }
}

void AccdPlanar::setRadius( ) {
  double Ledd = 1.3e47*mBH;
  double Ldisk = eDisk*mDot*Ledd;
  double R_sub = 1.6e-5*sqrt( Ldisk );
  double Rg = mBH*1e09*MSUN*G_CONST/DSQR(LIGHT_SPEED);
  R1 = Rg/(2.0*eDisk);
  R2 = R_sub; }

void AccdPlanar::setdLdR( ) {
  double Ledd = 1.3e47*mBH;
  double Ldisk = eDisk*mDot*Ledd;
  /* fill in dLdR vector that will be used for the rest of the calculations */
  bazinga::info(id,"Seeting dLdR.");
  for( int i=0;i<R->getMaxIndex();i++ ) {
      R->update( i );
      double val = 0.0;
      val = (3.0*G_CONST*mBH*1e09*MSUN)*mDot*Ledd*pow(R->get(),-s)/( 2.0*DSQR(LIGHT_SPEED) );
      gsl_vector_set( dLdR, R->getIndex( ), val ); }
  /* Now when we have all set let us save what we have just calculated - dLdR */
  bazinga::info(id,"Saving dLdR.");
  bazinga::save_GSLVector( "dLdR_"+this->whoAmI( ), R->getRadius_GSLVector( ), dLdR, cfg->get<std::string>("output") ); }

double AccdPlanar::getvext( double _R ) {
  double Ledd = 1.3e47*mBH;
  double Ldisk = eDisk*mDot*Ledd;
  double vext = 0.0;
  double Fdisk;
  if( luminosityConstNu ) { Fdisk = (3.0*G_CONST*mBH*1e09*MSUN)*mDot*Ledd/( DSQR(LIGHT_SPEED)*8.0*M_PI*pow(R2,3.0) ); }
  else { Fdisk = (3.0*G_CONST*mBH*1e09*MSUN)*mDot*Ledd/( DSQR(LIGHT_SPEED)*8.0*M_PI*pow(_R,3.0) ); }
  double Teff = pow( Fdisk/SIGMA_SB, 0.25 );
  vext = 3.92*Teff*K_BOLTZMAN/PLANCK_H;
  return vext; }
