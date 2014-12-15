/**
    @file synchrotron.cpp
    @author Mateusz Janiak
*/

#include "synchrotron.hpp"
#include "electrons.hpp"

synchrotron::synchrotron( scfgp* _cfg, jetGeometry* _r, electrons* _ele, magneticField* _B, std::string _id ) : energyDissProc( _cfg, _r, _ele, _id ), B(_B) {
  /* requested parameters */
  cfg -> request<int>("N"+id,200,&N);
  cfg -> request<std::string>("lumModel","blob",&lumModel);
  cfg -> request<double>("thetaJ",0.1,&thetaJ);
  cfg -> request<double>("Gamma",10.0,&Gamma);
  cfg -> request<int>("SABS",1,&SABS);
  cfg -> request<int>(id+"LuminosityConstU",0,&luminosityConstU);  
  cfg -> request<int>(id+"LuminosityConstNu",0,&luminosityConstNu);  
  cfg -> request<int>("setSSCBlob",0,&setSSCBlob);  

  cfg -> updateRequests( );

  vpMin = 1.0;
  /* here we set vpMAX to a maximum value specified by maximal value of magnetic field at R0 */
  vpMax = 100.0*A43*DSQR(ele->getGammaMax( ))*( B->getMaxB( )/B_CR )*mec2h;
  epMin = PLANCK_H/mec2;
  epMax = vpMax/mec2h;
  
  allocateUpe( );
  allocateUpeR( );
  allocateLpv( );
  allocateLvPoint( );
  allocateLvPointAvg( );

  for( int i=0;i<N;i++ ) {
    set_ep( i, epMax*pow( epMax/epMin,(double)i/(double)(N-1)) );
    set_vp( i, vpMin*pow( vpMax/vpMin,(double)i/((double)N-1)) ); }

  dLogE = log(get_ep(1)/get_ep(0)); }

synchrotron::~synchrotron( ) {
  freeLpv( );
  freeLvPoint( );
  freeLvPointAvg( );
  freeUpe( ); }

void synchrotron::printInfo( ) {
  bazinga::info(id,"Info");
  bazinga::print_info(id,"N",N);
  bazinga::print_info(id,"SABS",SABS);
  if( setSSCBlob ) { 
    bazinga::print_info(id,"setSSCBlob",setSSCBlob); }

  if( luminosityConstU || luminosityConstNu ) { 
    bazinga::warning(id,"Using constant u' to calculate luminosity.");
    bazinga::warning(id,"Using constant v' to calculate luminosity."); }
}

double synchrotron::dotg( double g ) { return B->get_uB( ); }

void synchrotron::update( ) {
  setLpv( );
  double lum_to_upe = 0.0;
  if( lumModel == "blob" || setSSCBlob ) { lum_to_upe = 2.0*M_PI*thetaJ*thetaJ*r->get( )*r->get( )*LIGHT_SPEED; }
  if( lumModel == "steady" && !setSSCBlob ) { lum_to_upe = 2.0*M_PI*r->get( )*r->getDr( )*thetaJ*Gamma*LIGHT_SPEED; }

  for( int i=0;i<N;i++ ) {
      set_ep( i, epMin*pow(epMax/epMin,(double)i/((double)N-1)) );
      set_upe( i, mec2h*get_Lpv(i)/lum_to_upe ); }

  set_upe_r( B->get_uB( ) );
  //  std::cout << "BFIELD " << r->get() << " " << B->getB();
  flag_upe_r = false; }

void synchrotron::setLpv( ) {
  for( int i=0;i<N;i++ ) {
    set_vp( i, vpMin*pow( vpMax/vpMin,(double)i/((double)N-1)) );
    set_Lpv( i, calculateLpv( get_vp(i) ) ); }
}

double synchrotron::FS( double t, double* jS, double* sigmaS ) {
  double K13,K43;  
  if(t<50.0){
    K13 = gsl_sf_bessel_Knu(1.0e0/3.0e0,t);
    K43 = gsl_sf_bessel_Knu(4.0e0/3.0e0,t);
    (*sigmaS) = t*(K43*K43-K13*K13);
    (*jS) = t*t*(K43*K13-3.0e0/5.0e0*(*sigmaS));
  } else
    (*sigmaS) = (*jS) = 0.0;
  return 0; }

double synchrotron::calculateLpv( double v ) {
  double a,jS,sigmaS,corr;
  double vB,h,t,L,tau;
  
  vB = ELECTRON_CHARGE*B->getB( )/(2.0*M_PI*ELECTRON_MASS*LIGHT_SPEED);
  h = ele->getdLogGamma( );
  L = tau = 0.0e0;

  for( int i=0;i<ele->getN( );i++) {
      t = v/(3.0*DSQR( ele->getGamma(i))*vB );
      FS(t,&jS,&sigmaS);
      tau += bazinga::IntCor( i, ele->getN( ) )*sigmaS*ele->getNgamma(i)/pow(ele->getGamma(i),4);
      L += bazinga::IntCor( i, ele->getN( ) )*ele->getGamma(i)*ele->getNgamma(i)*jS; }

  
  if( SABS && (r->get( )>=r->getR0( )) ) {
      /* since in a there is r-r0 then first luminosity wont have absorption; original blazar code had different naming: 
	 all quantities calculated at r were tagged as r+rd - we do not follow this approach thus r>*(model->R0) and not 
	 r>=R0 (it gives errors) */
      /* I cheat; I set DR to be dr for r=R0; is it a bad approximation? I dunno */
      double DR = 0;
      if( r->get() == r->getR0( ) ) { DR = r->getDr(); }
      else { DR = r->get( )-(r->getR0( )); }

      a = 2.0*DR/(beta(Gamma)*Gamma);
      tau *= h*3.0/(4.0*M_PI*a*a)*2.0*sqrt(3.0)*M_PI/15.0*ELECTRON_CHARGE/B->getB( );
      corr = (tau>1.0e-5) ? (1.0-exp(-tau))/tau : (1.0-0.5*tau+1.0/6.0*tau*tau-1.0/24.0*pow(tau,3)+1.0/120.0*pow(tau,4)-1.0/720.0*pow(tau,5));
  } else { corr = 1.0e0; }
  
  if( luminosityConstU || luminosityConstNu ) {   L *= h*3.0*sqrt(3.0)*SIGMA_T*LIGHT_SPEED*B->get_uB( injRm )/(M_PI*vB)*corr; }
  else {   L *= h*3.0*sqrt(3.0)*SIGMA_T*LIGHT_SPEED*B->get_uB( )/(M_PI*vB)*corr; }
  
  return L; }

double synchrotron::iterate( ) {
  double sumN = 0.0;
  double err = 0.0;
  gsl_vector* temp_upe = gsl_vector_alloc( N );
  
  /* copy current upe to temp_upe */
  for( int i=0;i<N;i++ ) { gsl_vector_set( temp_upe, i, get_upe( i ) ); }

  /* calculate synchrotron energy densities and luminosities with new Ngamma*/
  update( );
	  
  for( int i=0;i<N;i++ ) {
      sumN += get_upe( i );
      err += fabs( get_upe( i ) - gsl_vector_get( temp_upe, i ) ); }
  err /= sumN;
  gsl_vector_free( temp_upe );
  temp_upe = NULL;

  return err; }

//void synchrotron::setJetFrequency( )
//{
//  double theta, Doppler;
//  
//  /* part of the jet that is closer to observer */
//  theta = (*(model->ThetaObs)<*(model->ThetaJ)) ? 0.0 : (*(model->ThetaObs)-*(model->ThetaJ));
//  
//  Doppler = 1.0/(*(model->Gamma)*(1.0-cos(theta)*model->beta(*(model->Gamma))));
//  print_info(" Close Doppler",Doppler);
//  if( vjetMAX < Doppler*vpMAX ) vjetMAX = Doppler*vpMAX;
//
//  /* part of the jet that is further to observer */
//  theta = *(model->ThetaObs)+*(model->ThetaJ);
//  Doppler = 1.0/(*(model->Gamma)*(1.0-cos(theta)*model->beta(*(model->Gamma))));
//  print_info(" Far Doppler",Doppler);
//
//  if( vjetMIN > Doppler*vpMIN ) vjetMIN = Doppler*vpMIN;
//
//  print_info(" vjetMIN", vjetMIN);
//  print_info(" vjetMAX", vjetMAX);
//
//  for( int i=0;i<N;i++ )
//    set_vJet( i, vjetMIN*pow( vjetMAX/vjetMIN,(double)i/((double)N-1)) );
//}
