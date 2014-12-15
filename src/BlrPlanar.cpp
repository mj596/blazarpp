/**
   @file BlrPlanar.cpp
   @author Mateusz Janiak
*/

#include "BlrPlanar.hpp"
#include "electrons.hpp"

BlrPlanar::BlrPlanar( scfgp* cfg, jetGeometry* r, electrons* ele, std::string id ) : externalRadiationPlanar( cfg, r, ele, id ) {
  /* request parameters */
  cfg -> request<double>( id+"e", 10.0, &e );
  cfg -> request<double>( id+"cf", 0.1, &cf );
  cfg -> request<double>( id+"q1", 1.0, &q1 );
  cfg -> request<double>( id+"q2", 1.0, &q2 );
  cfg -> request<double>( id+"R", 0.0, &Rext );
  cfg -> request<std::string>( id+"Strat", "NEW", &strat );

  cfg -> updateRequests( );

  if( strat == "NEW" ) {
    if( R1 == 0.0 || R2 == 0.0 ) { setRadius( ); }
    else { bazinga::warning(id,"Overwritting R1 and R2 values from config file!"); }
  }
  else if( strat == "M03" ) {
    if( Rext == 0.0 ) { setRadius( ); }
    else { bazinga::warning(id,"Overwritting Rext value from config file!"); }
  }

  vpMin = 0.01*DSQR( ele -> getGammaMin( ) )*getvext( Rext );
  vpMax = 100.0*DSQR( ele -> getGammaMax( ) )*getvext( Rext );

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

BlrPlanar::~BlrPlanar( ) {
  freeUpe( );
  freeLpv( );
  freeLvPoint( );
  freeLvPointAvg( ); }

void BlrPlanar::printInfo( ) {
  bazinga::info(id,"Info");
  bazinga::print_info(id,"N",N);
  bazinga::print_info(id,"Using stratification version",strat);
  if( strat == "NEW" ) {
    bazinga::print_info(id,"radius R1",R1,"cm");
    bazinga::print_info(id,"radius R2",R2,"cm");
    bazinga::print_info(id,"stratification index",s); }
  if( strat == "M03" ) {
    bazinga::print_info(id,"radius Rext",Rext,"cm");
    bazinga::print_info(id,"radius R1",R1,"cm");
    bazinga::print_info(id,"radius R2",R2,"cm");
    bazinga::print_info(id,"stratfication index q1",q1);
    bazinga::print_info(id,"stratfication index q2",q2); }
  
  bazinga::print_info(id,"Avg energy (in external frame)",e,"eV");
  bazinga::print_info(id,"Avg frequency (in external frame)",getvext( Rext ),"Hz");
  bazinga::print_info(id,"Clouds covering factor cf",cf);
  bazinga::print_info(id,"alpha",alpha);
  bazinga::print_info(id,"vpMin",vpMin,"Hz");
  bazinga::print_info(id,"vpMax",vpMax,"Hz");
  if( approx ) { bazinga::print_info(id,"Using approximate dependence on R"); }
  if( luminosityConstU ) { bazinga::warning(id,"Using constant u' to calculate luminosity."); }
  if( luminosityConstNu ) { bazinga::warning(id,"Using constant v_ext to calculate luminosity."); }
}

void BlrPlanar::setRadius( ) {
  double Ledd = 1.3e47*mBH;
  double Ldisk = eDisk*mDot*Ledd;
  double R_sub = 1.6e-5*sqrt( Ldisk );
  if( strat == "NEW" ) {
    R1 = 0.1*R_sub;
    R2 = R_sub;
      /* Rext is set only to avoid 'if's in constructor when setting vpMIN and vpMAX; Rext is just set to R1 */
    Rext = R1; }
  else if( strat == "M03" ) {
    Rext = 0.1*R_sub;
    /* values R1 and r2 in M03 are only bounds for numerical reasons */
    R1 = 1.0e-5*R_sub; // ???
    R2 = 10.0*R_sub; }
}

void BlrPlanar::setdLdR( ) {
  double Ledd = 1.3e47*mBH;
  double Ldisk = eDisk*mDot*Ledd;
  /* fill in dLdR vector that will be used for the rest of the calculations */
  bazinga::info(id,"Seeting dLdR.");
  for( int i=0;i<R->getMaxIndex();i++ ) {
    R -> update( i );
    double val = 0.0;
    if( strat == "NEW" ) { val = cf*Ldisk*ctau( )*pow(R->get( ),-s); }
    if( strat == "M03" ) { 
      val = cf*Ldisk/(R->get()*(q1+q2)/(q1*q2));
      if( R->get() <= Rext ) { val *= pow(R->get()/Rext,q1); }
      if( R->get() > Rext ) { val *= pow(R->get()/Rext,-q2); }
    }
    gsl_vector_set( dLdR, R->getIndex( ), val );
  }
  /* Now when we have all set let us save what we have just calculated - dLdR */
  bazinga::info(id,"Saving dLdR.");
  bazinga::save_GSLVector( "dLdR_"+this->whoAmI( ), R->getRadius_GSLVector( ), dLdR, cfg->get<std::string>("output") ); }

double BlrPlanar::getvext( double _R ) {
  vext = e*eV2erg/PLANCK_H;
  return vext; }
