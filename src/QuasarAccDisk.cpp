 /**
    @file QuasarAccDisk.cpp
    @author Mateusz Janiak
*/

#include "QuasarAccDisk.hpp"

QuasarAccDisk::QuasarAccDisk( scfgp* _cfg, std::string _id ) : baseClass( _cfg, NULL, _id ) {
  /* requested parameters */
  cfg -> request<double>("mBH",1.0,&mBH);
  cfg -> request<double>("mDot",1.0,&mDot);
  cfg -> request<double>("eDisk",0.1,&eDisk);
  cfg -> request<int>("N" + id, 200, &N );
  cfg -> request<double>(id + "R1", 0.0, &R1 );
  cfg -> request<double>(id + "R2", 0.0, &R2 );
  cfg -> request<std::string>(id + "TemplateFilename", "dat/elvisRL.dat", &templateFilename );

  cfg -> request<int>(id + "SaveLr", 0, &save_Lr);

  cfg -> updateRequests( );

  setEnergetics( );

  if( R1 == 0.0 || R2 == 0.0 ) { setRadius( ); }

  vMin = 1.0e12;
  vMax = 1.0e18;

  numLines = 0;

  allocateLv( );

  for( int i=0;i<N;i++ ) { set_v( i, vMin*pow( vMax/vMin,(double)i/((double)N-1)) ); }
  R = new logGeometry( R1, R2, 1000 );
}

QuasarAccDisk::~QuasarAccDisk( ) {
  freeLv( );
  freeLvTemplate( );
}

double QuasarAccDisk::loglogIntegrate( gsl_vector* x, gsl_vector* y ){
  double integral = 0.0;
  for( int i = 0; i < x->size-1; ++i ) {
    integral += 0.5*( gsl_vector_get(y, i+1)/gsl_vector_get(x, i+1)+gsl_vector_get(y, i)/gsl_vector_get(x, i) )*( gsl_vector_get(x,i+1)-gsl_vector_get(x,i) ); }
  return integral;
}

void QuasarAccDisk::scaleTemplate( ) {
  double integral = loglogIntegrate( vTemplate, LvTemplate );
  double correctionFactor = Ldisk/integral;
  bazinga::print_info( id, "Template correction factor", correctionFactor );
  gsl_vector_scale( LvTemplate, correctionFactor );
  integral = loglogIntegrate( vTemplate, LvTemplate );
  bazinga::print_info( id, "Template Lbol after correction", integral );
  
  // Find peak in IR, Optical and X-ray
  double IRmin = 1.0e10;
  double IRmax = 2.0e14;
  double OPTmin = IRmax;
  double OPTmax = 7.0e16;
  double Xmin = OPTmax;
  double Xmax = 1.0e20;

  double vIRpeak = 0.0;
  double vOPTpeak = 0.0;
  double vXpeak = 0.0;

  double vLvIRpeak = 0.0;
  double vLvOPTpeak = 0.0;
  double vLvXpeak = 0.0;
  
  for( int i=0; i<vTemplate->size-1; ++i ) {
    if( gsl_vector_get( vTemplate, i ) >= IRmin && gsl_vector_get( vTemplate, i ) <= IRmax )
      if( gsl_vector_get( LvTemplate, i ) > vLvIRpeak ) {
	vLvIRpeak = gsl_vector_get( LvTemplate, i );
	vIRpeak = gsl_vector_get( vTemplate, i ); }
  }
  bazinga::print_info( id, "IR template v peak", vIRpeak );
  bazinga::print_info( id, "IR template vLv peak", vLvIRpeak );

  for( int i=0; i<vTemplate->size-1; ++i ) {
    if( gsl_vector_get( vTemplate, i ) >= OPTmin && gsl_vector_get( vTemplate, i ) <= OPTmax )
      if( gsl_vector_get( LvTemplate, i ) > vLvOPTpeak ) {
	vLvOPTpeak = gsl_vector_get( LvTemplate, i );
	vOPTpeak = gsl_vector_get( vTemplate, i ); }
  }
  bazinga::print_info( id, "OPT template v peak", vOPTpeak );
  bazinga::print_info( id, "OPT template vLv peak", vLvOPTpeak );

  for( int i=0; i<vTemplate->size-1; ++i ) {
    if( gsl_vector_get( vTemplate, i ) >= Xmin && gsl_vector_get( vTemplate, i ) <= Xmax )
      if( gsl_vector_get( LvTemplate, i ) > vLvXpeak ) {
	vLvXpeak = gsl_vector_get( LvTemplate, i );
	vXpeak = gsl_vector_get( vTemplate, i ); }
  }
  bazinga::print_info( id, "X template v peak", vXpeak );
  bazinga::print_info( id, "X template vLv peak", vLvXpeak );
}
  
int QuasarAccDisk::saveTemplate( ) {
  std::string type = "LvTemplate_"+this->whoAmI( );
  bazinga::info(this->whoAmI( ),"Saving luminosity.");

  gsl_vector* LvTemplateTemp = gsl_vector_alloc( numLines );
  gsl_vector_set_zero( LvTemplateTemp );
  gsl_vector_memcpy( LvTemplateTemp, LvTemplate );
  gsl_vector_div( LvTemplateTemp, vTemplate );
    
  bazinga::save_GSLVector( type, this->vTemplate, LvTemplateTemp, cfg->get<std::string>("output") );
  gsl_vector_free( LvTemplateTemp );
  LvTemplateTemp = NULL;
 }

int QuasarAccDisk::readTemplate( ){
  std::ifstream input( templateFilename.c_str() );
  bazinga::print_info( id, "Reading quasar template data from file ",templateFilename.c_str());
     if( input == NULL )
	{
	  bazinga::error(id,"Opening file:",templateFilename);
	  return 1;
	}
      else
	{
	  // check how many lines there are
	  std::string unused;
	  while ( std::getline(input, unused) )
	    ++numLines;
	  
	  bazinga::print_info(id,"Found",numLines,"lines");
	  
	  // allocate space for vectors
	  allocateLvTemplate( );

	  // read data from template
	  input.clear();
	  input.seekg(0, std::ios::beg);
	  double x, y;
	  int i = 0; // gsl counter
	  while( std::getline(input, unused) ) {
	  std::istringstream ss(unused);
	  std::istream_iterator<std::string> begin(ss), end;
	  
	  //putting all the tokens in the vector
	  std::vector<std::string> arrayTokens(begin, end); 
	  
	  gsl_vector_set( vTemplate, i, atof( arrayTokens[0].c_str() ) );
	  gsl_vector_set( LvTemplate, i, atof( arrayTokens[1].c_str() ) );			 
	  ++i; }
	  return 0;
	}
}

void QuasarAccDisk::printInfo( ) {
  bazinga::info( id, "Info" );
  bazinga::print_info( id, "N", N );
  bazinga::print_info( id, "Gravitational radius", Rg, "cm" );
  bazinga::print_info( id, "R1", R1, "cm" );
  bazinga::print_info( id, "R2", R2, "cm" );
  bazinga::print_info( id, "Save luminosity", save_Lr );
  bazinga::print_info(id,"BH mass",mBH);
  bazinga::print_info(id,"eta disk",eDisk);
  bazinga::print_info(id,"m dot",mDot);
  bazinga::print_info(id,"Eddington luminosity",Ledd,"erg/s");
  bazinga::print_info( id, "Accretion disk luminosity", Ldisk, "erg/s");
}

double QuasarAccDisk::getTemperature( double _r ) {
  double Risco = 6.0*Rg;
  return pow( ( 3.0*G_CONST*1.3e56*pow(mBH,2)*MSUN*mDot )*( 1.0-sqrt(Risco/_r) )/( 8.0*M_PI*pow(_r, 3)*SIGMA_SB*pow(LIGHT_SPEED,2) ), 0.25 );
}

void QuasarAccDisk::allocateLv( ) {
  v = gsl_vector_alloc( N );
  Lv = gsl_vector_alloc( N );
  bazinga::print_GSLVector_allocated_memory( id, v );
  bazinga::print_GSLVector_allocated_memory( id, Lv );
  gsl_vector_set_zero( v );
  gsl_vector_set_zero( Lv ); }

void QuasarAccDisk::freeLv( ) {
  gsl_vector_free( v );
  gsl_vector_free( Lv ); }

void QuasarAccDisk::allocateLvTemplate( ) {
  vTemplate = gsl_vector_alloc( numLines );
  LvTemplate = gsl_vector_alloc( numLines );
  bazinga::print_GSLVector_allocated_memory( id, vTemplate );
  bazinga::print_GSLVector_allocated_memory( id, LvTemplate );
  gsl_vector_set_zero( vTemplate );
  gsl_vector_set_zero( LvTemplate ); }

void QuasarAccDisk::freeLvTemplate( ) {
  gsl_vector_free( vTemplate );
  gsl_vector_free( LvTemplate ); }

void QuasarAccDisk::calculateLuminosity( ) {
  bazinga::info(id,"Calculating luminosity.");
  gsl_vector* tempLv = gsl_vector_alloc( N );

  for( int i=0;i<R->getMaxIndex();i++ ) {
    R -> update( i );
    gsl_vector_set_zero( tempLv );

    for( int j=0;j<N;j++ ) {
      if( save_Lr ) {
	gsl_vector_set( tempLv, j, calculateLv( gsl_vector_get(v,j), i ) ); }
      gsl_vector_set( Lv, j, gsl_vector_get( Lv, j ) + calculateLv( gsl_vector_get(v,j), i ) );
    }

    if( save_Lr ) {
      std::string type = "Lv_";
      type += this->whoAmI( );
      bazinga::print_info(id,"Saving QLuminosity.",R->get( )/R->getR0( ));
      bazinga::save_GSLVector( type, v, tempLv, R->get( )/R->getR0( ), cfg->get<std::string>("output") ); }
  }

  gsl_vector_free( tempLv );
  tempLv = NULL;

  std::string type = "Lv_";
  type += this->whoAmI( );
  bazinga::info(id,"Saving QLuminosity.");
  bazinga::save_GSLVector( type, v, Lv, cfg->get<std::string>("output") );

  double integral = 0.0;
  tempLv = gsl_vector_alloc( N );
  gsl_vector_set_zero( tempLv );
  gsl_vector_memcpy( tempLv, Lv );
  gsl_vector_mul( tempLv, v );
  integral = loglogIntegrate( v, tempLv );
  bazinga::print_info( id, "Acc Disk Lbol", integral );    
}

double QuasarAccDisk::calculateLv( double v, int i ) {
  double val = ( 16.0*pow(M_PI,2)*PLANCK_H*cos(0.0)*pow(v,3) )/( pow(LIGHT_SPEED,2) );
  return val*( R->get(i)*R->getDr(i) )/( exp( (PLANCK_H*v)/(K_BOLTZMAN*getTemperature(R->get(i)) ) ) - 1.0 );
}

void QuasarAccDisk::setRadius( ) {
  /* gravitational radius */
  Rg = (2.0*G_CONST*mBH*1.0e9*MSUN)/pow(LIGHT_SPEED,2);
  R1 = 6.0*Rg;
  double R_sub = 1.6e-5*sqrt( Ldisk );
  R2 = R_sub; }

void QuasarAccDisk::setEnergetics( ) {
   Ledd = 1.3e47*mBH;
   Ldisk = eDisk*mDot*Ledd; }
