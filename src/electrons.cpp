 /**
    @file electrons.cpp
    @author Mateusz Janiak
*/

#include "electrons.hpp"
#include "energyDissProc.hpp"
#include "synchrotron.hpp"

electrons::electrons( scfgp* _cfg, jetGeometry* _r, std::string _id ) : baseClass(_cfg, _r, _id ) {
  /* requested parameters */
  cfg -> request<int>("N"+id,200,&N);
  cfg -> request<std::string>("eleModel","blob",&eleModel);
  cfg -> request<std::string>("injModel","REC",&injModel);
  cfg -> request<std::string>("lumModel","blob",&lumModel);
  cfg -> request<double>("p1",2.0,&p1);
  cfg -> request<double>("p2",2.0,&p2);
  cfg -> request<double>("gammaMin",1.0,&gammaMin);
  cfg -> request<double>("gammaBreak",0.0,&gammaBreak);
  cfg -> request<double>("gammaMax",1.0e2,&gammaMax);
  cfg -> request<double>("injRm",2.0e17,&injRm);
  cfg -> request<double>("injSigma",0.6,&injSigma);
  cfg -> request<double>("eleK",0.0,&eleK);
  cfg -> request<double>("eEle",0.1,&eEle);
  cfg -> request<double>("eDiss",0.1,&eDiss);
  cfg -> request<double>("Gamma",10.0,&Gamma);
  cfg -> request<double>("NeNp",0.1,&NeNp);
  cfg -> request<double>("mBH",1.0,&mBH);
  cfg -> request<double>("eDisk",0.1,&eDisk);
  cfg -> request<double>("mDot",1.0,&mDot);
  cfg -> request<double>("eJet",0.3,&eJet);
  cfg -> request<double>("sigmaB",0.1,&sigmaB);
  cfg -> request<int>("evol",0,&evol);
  cfg -> request<int>("saveElectrons",1,&saveElectrons);
  cfg -> request<int>("saveElectronsAvg",1,&saveElectronsAvg);

  cfg -> updateRequests( );

  /* check if spectral indices are OK; for steady model p1 must be lower than 1 and p2 must be greater than 2; otherwise calculation of break electron energy will fail;
     for blob model it is all good to set any indices you like */
  if( (p1 > 1.0 || p2 < 2.0) && eleModel == "steady" && gammaBreak == 0.0 ) {
    bazinga::warning(id,"Can't use steady model with p1>1.0 or p2<2.0!");
    bazinga::warning(id,"To overwrite that you can set gammaBreak explicity.");
    bazinga::warning(id,"Quit.");
    exit(0); }

  allocateGamma( );
  allocateTridiag( );
  
  for( int i=0;i<=N+1;i++ ) { gsl_vector_set( gamma, i, pow( gammaMax,(double)i/(double)N )); }

  setEnergetics( );
  setInjectionParameters( );
  dLogGamma = log( gsl_vector_get(gamma,1)/gsl_vector_get(gamma,0) ); }


//int electrons::read( ){
//  FILE *input;
//  stringstream s;
//
//  radius* r = new radius( model );
//  double radiusToSave;
//  int currentRadiusIndex;
//  double _gamma, _Ngamma;
//
//  for( int i=0;i<r->get_max_index();i++ )
//    {
//      r->update( i );
//      radiusToSave = r->get()/(*(model->R0));
//      currentRadiusIndex = r->get_index();
//      s << *(model->Output) << "/" << "Ngamma" << "_" << radiusToSave;
//      input = fopen(s.str().c_str(),"r");
//      print_info("Reading electron data from file ",s.str().c_str());
//      s.str("");
//
//      if( input == NULL )
//	{
//	  print_info("[***] Error opening file.");
//	  return 1;
//	}
//      else
//	{
//	  for( int i=0;i<N;i++ )
//	    {
//	      if( !fscanf( input,"%le %le\n",&_gamma,&_Ngamma ) ) print_info("[***] Error reading file.");
//	      gsl_vector_set( gamma, i+1, _gamma );
//	      gsl_vector_set( Ngamma, i, _Ngamma );
//	    }
//	  
//	  print_info("Done.");
//	  fclose(input);
//	}
//    }
//  
//  return 0;
//}
//

void electrons::printInfo( ) {
  bazinga::info(id,"Info");
  bazinga::print_info(id,"evol",evol);
  bazinga::print_info(id,"Electron model",eleModel);
  bazinga::print_info(id,"Injection model",injModel);
  bazinga::print_info(id,"N",N);
  bazinga::print_info(id,"Index p 1",p1);
  bazinga::print_info(id,"Index p 2",p2);
  bazinga::print_info(id,"Gamma MIN",gammaMin);
  bazinga::print_info(id,"Gamma MAX",gammaMax);
  
  if( injModel == "GAUSS" || injModel == "GAUSSMOD" ) {
    bazinga::print_info(id,"injection maximium R_m", injRm);
    bazinga::print_info(id,"injection disspersion sigma",injSigma); }
  else if( injModel == "TRI") {
    bazinga::print_info(id,"injection maximium R_m",injRm); }
  
  if( eleModel == "steady" && gammaBreak ) { bazinga::print_info(id,"Injection average gamma",avgGamma); }
  
  bazinga::print_info(id,"Gamma break",gammaBreak);
  bazinga::print_info(id,"Electron normalization",eleK);

  if( eleModel == "steady") {
    bazinga::print_info(id,"Luminosity model",lumModel);
    bazinga::print_info(id,"eta electrons",eEle);
    bazinga::print_info(id,"eta dissipation",eDiss);
    bazinga::print_info(id,"Jet Lorentz factor",Gamma);
    bazinga::print_info(id,"n_e/n_p",NeNp);
    bazinga::print_info(id,"BH mass",mBH);
    bazinga::print_info(id,"eta disk",eDisk);
    bazinga::print_info(id,"m dot",mDot);
    bazinga::print_info(id,"eta jet",eJet);
    bazinga::print_info(id,"magnetic sigma",sigmaB);
    bazinga::print_info(id,"Eddington luminosity",Ledd,"erg/s");
    bazinga::print_info(id,"Accretion disk luminosity",Ldisk,"erg/s");
    bazinga::print_info(id,"Jet power before dissipation",Ljet,"erg/s");
    bazinga::print_info(id,"Jet power after dissipation",LjetDiss,"erg/s");
    bazinga::print_info(id,"Magnetic energy flux",Lb,"erg/s");
    bazinga::print_info(id,"Kinetic energy of cold protons",Lprot,"erg/s"); }
}

electrons::~electrons( ) {
  freeGamma( );
  freeTridiag( ); }

/* evolve electrons given energy dissipation processes */
void electrons::evolve( ) {
  double err = 0.0;
  bazinga::info(id,"Cooling details:");
  /* print info about cooling */
  bazinga::info(id,"@ gammaMIN:");
  printCoolingInfo( gsl_vector_get(gamma,0) );
  bazinga::info(id,"@ gammaMAX:");
  printCoolingInfo( gsl_vector_get(gamma,N-1) );

  if( cfg->get<int>("Ssc") ) { bazinga::info(id,"Ssc loop:"); }
  do {
    /* solve equations for electrons */
    for( int i=0;i<N;i++ )
      {	gsl_vector_set( B, i, 1.0+r->getDr()*dgdr( 0.5*(gsl_vector_get(gamma,i)+gsl_vector_get(gamma,i+1)) )/(0.5*(gsl_vector_get(gamma,i+2)-gsl_vector_get(gamma,i))) ); }
    for( int i=0;i<N-1;i++ ) 
      {	gsl_vector_set( C, i, -r->getDr()*dgdr(0.5*(gsl_vector_get(gamma,i+1)+gsl_vector_get(gamma,i+2)) )/(0.5*(gsl_vector_get(gamma,i+2)-gsl_vector_get(gamma,i))));	}
    
    gsl_linalg_solve_tridiag(B,C,A,S,U);

    gsl_vector_memcpy( Ngamma, U );      
    
    if ( cfg->get<int>("Ssc") ) { /* in order to calculate corrected model for ssc after electron evolution we need to first identify which pointer in array corresponds to synchrotron */
      for( int j=0; j<EnDissProc.size(); j++ )
	if( EnDissProc[j]->whoAmI( ) == "syn" ) { 
	  err = ( static_cast<synchrotron*>( EnDissProc[j] ) ) -> iterate( );
	  bazinga::print_info(id,"  iteration error",err); }
    }
  } while ( err > ESP ); }

/* inject accelerated electrons */
double electrons::inject( ) {
  for( int i=1;i<=N;i++ )
    if( evol ) { gsl_vector_set( S, i-1, gsl_vector_get(Ngamma,i-1)+r->getDr()*tQ( gsl_vector_get(gamma, i), r->get( )) ); }
    else { gsl_vector_set( Ngamma, i-1, gsl_vector_get(Ngamma,i-1)+r->getDr()*tQ( gsl_vector_get(gamma, i), r->get()) ); }
}

double electrons::tQ( double g, double _radius ) {
  double x;

  if( eleModel == "steady" ) {
    if( g >= gammaMin && g <= gammaMax && _radius >= r->getR0( ) && _radius <= r->getRInjMax( ) ) {
      if( g <= gammaBreak ) { x = eleK*pow( gammaBreak, p1-p2 )*pow( g, -p1 ); }
      else { x = eleK*pow( g, -p2 ); }
    }
    else { x = 0.0; }
    return x/(LIGHT_SPEED*beta(Gamma)*Gamma); }

  if( eleModel == "blob" ) {
    double corr = 1.0;
    if( injModel == "GAUSS" ) { corr = gauss( _radius ); }
    if( injModel == "GAUSSMOD" ) { corr = gaussmod( _radius ); }
    if( injModel == "TRI" ) { corr = tri( _radius ); }

    if( g >= gammaMin && g <= gammaMax && _radius >= r->getR0( ) && _radius <= r->getRInjMax( ) ) {
      x = eleK*pow(g,-p1)*pow(1.0 + pow(g/gammaBreak,4.0),(p1-p2)/4.0)/(LIGHT_SPEED*beta(Gamma)*Gamma); }
    else { x = 0.0; }
    return corr*x; }
}

double electrons::gauss( double _radius ) {
  double x = _radius/injRm;
  return exp( -pow(x-1.0,2)/(pow( injSigma,2)) ); }

double electrons::gaussmod( double _radius ) {
  double x = _radius/injRm;
  return x*exp( (1-x)*(x+pow( injSigma,2)-1.0 )/pow( injSigma,2) ) ; }

double electrons::tri( double _radius ) {
  if( _radius < injRm ) { return 0.5+0.5*(_radius-r->getR0( ))/(injRm-r->getR0( )); }
  else if( _radius >= injRm ) { return 1.0+0.5*(_radius-injRm)/(injRm-r->getRInjMax( )); }
}

void electrons::printCoolingInfo( double g ) {
  for( int i=0; i<EnDissProc.size(); i++ ) { bazinga::print_info(id,EnDissProc[i]->whoAmI(),EnDissProc[i]->dotg( g )*A43sigTmec*g*g/(beta(Gamma)*LIGHT_SPEED*Gamma)); } 
  if( cfg->get<int>("Adiabatic") ) { bazinga::print_info(id,"adiabatic",cfg->get<double>("AdiabaticABG")*g/r->get( )); }
}

double electrons::dgdr( double g ) {
  double dotg = 0.0;
  for( int i=0; i<EnDissProc.size(); i++ ) { dotg += EnDissProc[i]->dotg( g ); }  
  dotg *= A43sigTmec; /* multiply by common factor */
  dotg *= g*g;
  dotg *= 1.0/(beta(Gamma)*LIGHT_SPEED*Gamma); /* convert from time to radius */  
  if( cfg->get<int>("Adiabatic") ) { dotg += 0.666667*g/r->get(); } /* adiabatic cooling */
  return dotg; }

double electrons::getGamma( int i ) { return gsl_vector_get( gamma,i+1 ); }
double electrons::getNgamma( int i ) { return gsl_vector_get( Ngamma,i ); }

void electrons::addEnDissProc( energyDissProc* _obj ) { EnDissProc.push_back( _obj ); }

void electrons::listEnDissProc( ) {
  std::stringstream s;
  for( int i=0; i<EnDissProc.size(); i++ ) { s << EnDissProc[i]->whoAmI() << " "; }
  bazinga::print_info(id,"Processes used in electron cooling",s.str( )); }

void electrons::allocateGamma( ) {
  gamma = gsl_vector_alloc( N+2 );
  Ngamma = gsl_vector_alloc( N );
  NgammaAvg = gsl_vector_alloc( N );
  bazinga::print_GSLVector_allocated_memory( id, gamma );
  bazinga::print_GSLVector_allocated_memory( id, Ngamma );
  bazinga::print_GSLVector_allocated_memory( id, NgammaAvg );
  gsl_vector_set_zero( gamma );
  gsl_vector_set_zero( Ngamma );
  gsl_vector_set_zero( NgammaAvg ); }

void electrons::allocateTridiag( ) {
  A = gsl_vector_alloc(N-1);
  B = gsl_vector_alloc(N);
  C = gsl_vector_alloc(N-1);
  S = gsl_vector_alloc(N);
  U = gsl_vector_alloc(N);
  bazinga::print_GSLVector_allocated_memory( id, A );
  bazinga::print_GSLVector_allocated_memory( id, B );
  bazinga::print_GSLVector_allocated_memory( id, C );
  bazinga::print_GSLVector_allocated_memory( id, S );
  bazinga::print_GSLVector_allocated_memory( id, U );
  gsl_vector_set_zero( B );
  gsl_vector_set_zero( S );
  gsl_vector_set_zero( U );
  gsl_vector_set_zero( A );
  gsl_vector_set_zero( C ); }

void electrons::freeGamma( ) {
  gsl_vector_free( gamma );
  gsl_vector_free( Ngamma );
  gsl_vector_free( NgammaAvg ); }

void electrons::freeTridiag( ) {
  gsl_vector_free( A );
  gsl_vector_free( B );
  gsl_vector_free( C );
  gsl_vector_free( S );
  gsl_vector_free( U ); }

void electrons::setInjectionParameters( ) {
  if( eleModel == "steady" ) {
    if( !gammaBreak || !eleK ) {
      avgGamma = mpme*eEle*eDiss*(Gamma-1.0)*(1.0+sigmaB)/( NeNp*(1.0-eDiss)*Gamma );
      avgGamma += 1.0;
      
      gammaBreak = solveGammaBreak( p1, p2, gammaMin, gammaMax, avgGamma );
      
      eleK = eEle*eDiss*Ljet/( (avgGamma-1.0)*mec2 );
      eleK /= f2( p1, p2, gammaMin, gammaMax, gammaBreak ); }
    else { bazinga::error(id,"Set both gammaBreak and eleK or do not set any of them."); }

    if( lumModel == "steady" ) { eleK /= (r->getNinj( )); }
  }
}

double f1( double p1, double p2, double gamma_min, double gamma_max, double gamma_break ) {
  return ( (1.0-pow(gamma_min/gamma_break,2.0-p1))/(2.0-p1) + (1.0-pow(gamma_break/gamma_max,p2-2.0))/(p2-2.0) )*pow(gamma_break,1.0-p2); }

double f2( double p1, double p2, double gamma_min, double gamma_max, double gamma_break ) {
  return ( (1.0-pow(gamma_min/gamma_break,1.0-p1))/(1.0-p1) + (1.0-pow(gamma_break/gamma_max,p2-1.0))/(p2-1.0) )*pow(gamma_break,1.0-p2); }

double functionGammaBreak( double x, void *params ) {
  /* x <- gamma_break */
  struct gamma_break_params *p = (struct gamma_break_params *) params;

  double p1 = p->p1;
  double p2 = p->p2;
  double gamma_min = p->gamma_min;
  double gamma_max = p->gamma_max;
  double avg_gamma = p->avg_gamma;
  
  return x*f1(p1, p2, gamma_min, gamma_max, x)/f2(p1, p2, gamma_min, gamma_max, x) - avg_gamma; }

double solveGammaBreak( double p1, double p2, double gamma_min, double gamma_max, double avg_gamma ) {
  bazinga::info("","Entering solve_gamma_break");
  //  std::cout << "Using p1: " << p1 << " p2: " << p2 << " gmin: " << gamma_min << " gmax: " << gamma_max << " avg gamma: " << avg_gamma << std::endl;
  int status;
  int iter = 0;
  int max_iter = 100;

  double r = 0;
  double x_low = 0;
  double x_high = 0;

  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc( T );

  gsl_function F;
  struct gamma_break_params params = { p1, p2, gamma_min, gamma_max, avg_gamma };

  F.function = &functionGammaBreak;
  F.params = &params;
  
  gsl_root_fsolver_set( s, &F, gamma_min, gamma_max );
  
  bazinga::print_header();
  bazinga::info("","Trying to find gamma_break ...");

  printf("%5s [%8s, %7s] %8s\n", "iter", "lower", "upper", "root");

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate( s );
      r = gsl_root_fsolver_root( s );
      x_low = gsl_root_fsolver_x_lower( s );
      x_high = gsl_root_fsolver_x_upper( s );
      status = gsl_root_test_interval( x_low, x_high, 0, 0.001 );
      
      if( status == GSL_SUCCESS ) { bazinga::print_info("Converged:"); }
      
      printf("%5d [%.2e %.2e] %.2e\n", iter, x_low, x_high, r);
    }
  
  while( status == GSL_CONTINUE && iter <= max_iter );

  gsl_root_fsolver_free( s );
  return r; }

void electrons::avgNgamma( ) {
  for( int i=0;i<N;i++ ) { gsl_vector_set( NgammaAvg, i, gsl_vector_get( NgammaAvg, i ) + getNgamma( i ) ); } 
}

void electrons::setEnergetics( ) {
   Ledd = 1.3e47*mBH;
   Ldisk = eDisk*mDot*Ledd;
   Ljet = 0.5*eJet*mDot*Ledd;
   LjetDiss = (1.0-eDiss)*Ljet;
   Lprot = (1.0-eDiss)*Ljet/(1.0+sigmaB);
   Lb = sigmaB*(1.0-eDiss)*Ljet/(1.0+sigmaB); }

double electrons::beta( double x ){ return sqrt(1.0-1.0/(x*x)); }

void electrons::saveInjection( ) {
  if( saveElectrons && r->ifSaveRadius( ) ) { 
    bazinga::info(id,"Saving injection");
    bazinga::save_GSLVectorEle( "Injection", gamma, S, r->getPosition( ), cfg->get<std::string>("output") ); }
}

void electrons::saveNgamma( ) {
  if( saveElectrons && r->ifSaveRadius( ) ) {
    bazinga::info(id,"Saving Ngamma");
    bazinga::save_GSLVectorEle( "Ngamma", gamma, Ngamma, r->getPosition( ), cfg->get<std::string>("output") ); }
}

void electrons::saveNgammaAvg( ) {
  if( saveElectronsAvg ) { 
    bazinga::info(id,"Saving NgammaAvg");
    bazinga::save_GSLVectorEle( "Ngamma_Avg", gamma, NgammaAvg, cfg->get<std::string>("output") ); }
}
