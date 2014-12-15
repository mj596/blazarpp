/**
    @file observer.cpp
    @author Mateusz Janiak
*/

#include "observer.hpp"
#include "electrons.hpp"
#include "energyDissProc.hpp"

observer::observer( scfgp* _cfg, jetGeometry* _r, std::string _id ) : baseClass(_cfg, _r, _id ) { 
  /* request parameters */
  cfg -> request<int>("N"+id,1000,&N);
  cfg -> request<double>("thetaObs",0.1,&thetaObs);
  cfg -> request<double>("thetaJ",0.1,&thetaJ);
  cfg -> request<double>("nu1",0,&nu1);
  cfg -> request<double>("nu2",0,&nu2);
  cfg -> request<double>("nu3",0,&nu3);
  cfg -> request<double>("nu4",0,&nu4);
  cfg -> request<double>("nu5",0,&nu5);
  cfg -> request<double>("nu6",0,&nu6);
  cfg -> request<double>("Gamma",10.0,&Gamma);
  cfg -> request<std::string>("lumModel","blob",&lumModel);
  cfg -> request<int>("saveLumPoint",0,&saveLumPoint);
  cfg -> request<int>("saveLumPointAvg",0,&saveLumPointAvg);
  cfg -> request<int>("saveExtPlVsRm",0,&saveExtPlVsRm);

  cfg -> updateRequests( );
  
  allocateLvPointSum( );
  allocateLvPointAvgSum( ); }

observer::~observer( ) {
  freeLvPointSum( );
  freeLvPointAvgSum( ); }

int observer::getFreqListSize( ){ return freqList.size(); }

double observer::getFreqList( int i ){ return freqList[i]; }

void observer::calculateAllFlares( ) {
  if( nu1 || nu2 || nu3 || nu4 || nu5 || nu6 ) {
    bazinga::info(id,"Calculating flares");
    for( int j=0; j<sizeEnDissProcPoint(); j++ ) {
      bazinga::info(id,"Calculating point source flare",EnDissProcPoint[j]->whoAmI( ));
      for( int k=0; k<getFreqListSize(); k++ ) {
	calculateFlare( EnDissProc[j], getFreqList(k) );
	bazinga::print_info(id,"frequency",getFreqList(k),"Hz"); }
    }
  }
}

void observer::printInfo( ) {
  bazinga::info(id,"Info");
  bazinga::print_info(id,"N",N);
  bazinga::print_info(id,"Observer @",thetaObs,"rad");
  if( nu1 || nu2 || nu3 || nu4 || nu5 || nu6 ) {
      bazinga::info(id,"Calculating flares");
      if( nu1 ) {
	freqList.push_back( nu1 );
	bazinga::print_info(id,"@ frequency",nu1); }
      if( nu2 ) {
	freqList.push_back( nu2 );
	bazinga::print_info(id,"@ frequency",nu2); }
      if( nu3 ) {
	freqList.push_back( nu3 );
	bazinga::print_info(id,"@ frequency",nu3); }
      if( nu4 ) {
	freqList.push_back( nu4 );
	bazinga::print_info(id,"@ frequency",nu4); }
      if( nu5 ) {
	freqList.push_back( nu5 );
	bazinga::print_info(id,"@ frequency",nu5); }
      if( nu6 ) {
	freqList.push_back( nu6 );
	bazinga::print_info(id,"@ frequency",nu6); }
  }
}

void observer::addEnDissProc( energyDissProc* _obj ) { EnDissProc.push_back( _obj ); }
void observer::addEnDissProcPoint( energyDissProc* _obj ) { EnDissProcPoint.push_back( _obj ); }
void observer::addExtPl( energyDissProc* _obj ) { ExtPl.push_back( _obj ); }
void observer::addExtSp( energyDissProc* _obj ) { ExtSp.push_back( _obj ); }
void observer::addExtGu( energyDissProc* _obj ) { ExtGu.push_back( _obj ); }
int observer::sizeEnDissProc( ) { return EnDissProc.size(); }
int observer::sizeEnDissProcPoint( ) { return EnDissProcPoint.size(); }
int observer::sizeExtPl( ) { return ExtPl.size(); }
int observer::sizeExtSp( ) { return ExtSp.size(); }
int observer::sizeExtGu( ) { return ExtGu.size(); }

void observer::listExtPl( ) {
  std::stringstream s;
  for( int i=0; i<ExtPl.size(); i++ ) s << ExtPl[i]->whoAmI() << " ";
  bazinga::print_info(id,"External Radiation Planar Sources",s.str()); }

void observer::listExtSp( ) {
  std::stringstream s;
  for( int i=0; i<ExtSp.size(); i++ ) s << ExtSp[i]->whoAmI() << " ";
  bazinga::print_info(id,"External Radiation Spherical Sources",s.str()); }

void observer::listExtGu( ) {
  std::stringstream s;
  for( int i=0; i<ExtGu.size(); i++ ) s << ExtGu[i]->whoAmI() << " ";
  bazinga::print_info(id,"External Radiation Gu Sources",s.str()); }

void observer::listEnDissProc( ) {
  std::stringstream s;
  for( int i=0; i<EnDissProc.size(); i++ ) s << EnDissProc[i]->whoAmI() << " ";
  bazinga::print_info(id,"Processes used to calculate intrinsic lumonisities",s.str()); }

void observer::listEnDissProcPoint( ) {
  std::stringstream s;
  for( int i=0; i<EnDissProcPoint.size(); i++ ) s << EnDissProcPoint[i]->whoAmI() << " ";
  bazinga::print_info(id,"Processes used to calculate point lumonisities",s.str()); }

void observer::calculatePointLuminosity( energyDissProc* x ) {
  bazinga::info(id,"Calculating point source luminosity", x -> whoAmI( ) );
  double Doppler = 1.0/(Gamma*(1.0-cos( thetaObs )*beta(Gamma)));
  for( int i=0;i<x->N;i++ ) {
    x->set_vPoint( i, x->get_vp( i )*Doppler );
    if( lumModel == "blob" ) { x->set_LvPoint( i, x->get_Lpv( i )*pow(Doppler, 3) ); }
    if( lumModel == "steady" ) { x->set_LvPoint( i, x->get_Lpv( i )*pow(Doppler, 2)/Gamma ); }
  }
}

void observer::avgPointLuminosity( energyDissProc* x ) {
  bazinga::info(id,"Averaging spectra", x -> whoAmI( ) );
  for( int i=0;i<x->N;i++ ) {
    if( lumModel == "blob" ) { x->set_LvPointAvg( i, x->get_LvPointAvg( i ) + x->get_LvPoint( i )/r->getNinj( ) ); }
    if( lumModel == "steady" ) { x->set_LvPointAvg( i, x->get_LvPointAvg( i ) + x->get_LvPoint( i ) ); }
  }
}

void observer::sumPointLuminosity( ) {
  /* preparing vector for interpolated data */
  gsl_vector* xx = gsl_vector_alloc( N );
  gsl_vector* yy = gsl_vector_alloc( N );
  
  gsl_vector_set_zero( xx );
  gsl_vector_set_zero( yy );
  
  for( int i=0;i<N;i++ ) { gsl_vector_set( xx, i, vmin*pow(vmax/vmin,(double)i/(double)N) ); }
  
  for( int i=0; i<sizeEnDissProcPoint(); i++ ) {
    bazinga::print_info(id,"Add",EnDissProcPoint[i]->whoAmI());
    
    /* we will copy actual values stores in matrices to gsl_vectors; y vector is nu L nu!! */
    gsl_vector* temp_x = gsl_vector_alloc( EnDissProcPoint[i]->N );
    gsl_vector* temp_y = gsl_vector_alloc( EnDissProcPoint[i]->N );
    
    gsl_vector_memcpy( temp_x, EnDissProcPoint[i]->vPoint );
    gsl_vector_memcpy( temp_y, EnDissProcPoint[i]->LvPoint );
      
    /* nu F nu! */
    gsl_vector_mul( temp_y, temp_x );
      
    /* fill in the final vector and add only those values which are within boundaries */
    for( int k=0; k<N; k++ ) { gsl_vector_set( yy, k, gsl_vector_get( yy, k ) + interpolate( temp_x, temp_y, gsl_vector_get( xx, k ) ) ); }
    
    gsl_vector_free( temp_x );
    gsl_vector_free( temp_y ); }
  
  /* copy vector to matrix */
  for( int i=0;i<N;i++ ) {
    gsl_vector_set( vPointSum, i, gsl_vector_get( xx, i ) );    
    gsl_vector_set( LvPointSum, i, gsl_vector_get( yy, i ) );
 }
  
  gsl_vector_free( xx );
  gsl_vector_free( yy ); }

void observer::addQAccdLuminosity( QuasarAccDisk* QAccd ) {
  bazinga::info(id,"Adding QAccd spectra to averaged point source luminosities ...");

  double QAccd_vmin = gsl_vector_min( QAccd -> vTemplate );
  double QAccd_vmax = gsl_vector_max( QAccd -> vTemplate );
  bazinga::print_info(id, "QAccd vmin", QAccd_vmin );
  bazinga::print_info(id, "QAccd vmax", QAccd_vmax );

  double interpolatedValue = 0.0;
  for( int i=0;i<N;i++ ) {
    interpolatedValue = gsl_vector_get( LvPointAvgSum, i );
    if( gsl_vector_get( vPointSum, i ) > QAccd_vmin && gsl_vector_get( vPointSum, i ) < QAccd_vmax )
      interpolatedValue += interpolate( QAccd -> vTemplate, QAccd -> LvTemplate, gsl_vector_get( vPointSum, i ) );
    
    gsl_vector_set( LvPointAvgSum, i, interpolatedValue );
  }
}

void observer::sumPointAvgLuminosity( )
{
  bazinga::info(id,"Summing calculated averaged point source luminosities ...");

  /* preparing vector for interpolated data */
  gsl_vector* xx = gsl_vector_alloc( N );
  gsl_vector* yy = gsl_vector_alloc( N );
  
  gsl_vector_set_zero( xx );
  gsl_vector_set_zero( yy );
  
  for( int i=0;i<N;i++ ) { gsl_vector_set( xx, i, vmin*pow(vmax/vmin,(double)i/(double)N) ); }
    
  for( int i=0; i<sizeEnDissProcPoint(); i++ ) {
    bazinga::print_info(id,"Add",EnDissProcPoint[i]->whoAmI());
    
    /* we will copy actual values stores in matrices to gsl_vectors; y vector is nu L nu!! */
    gsl_vector* temp_x = gsl_vector_alloc( EnDissProcPoint[i]->N );
    gsl_vector* temp_y = gsl_vector_alloc( EnDissProcPoint[i]->N );
      
    gsl_vector_memcpy( temp_x, EnDissProcPoint[i]->vPoint );
    gsl_vector_memcpy( temp_y, EnDissProcPoint[i]->LvPointAvg );
    
    /* nu F nu! */
    gsl_vector_mul( temp_y, temp_x );
    
    /* fill in the final vector and add only those values which are within boundaries */
    for( int k=0; k<N; k++ ) { gsl_vector_set( yy, k, gsl_vector_get( yy, k ) + interpolate( temp_x, temp_y, gsl_vector_get( xx, k ) ) ); }

    gsl_vector_free( temp_x );
    gsl_vector_free( temp_y ); }
  
  /* copy vector to matrix */
  for( int i=0;i<N;i++ ) {
    gsl_vector_set( vPointSum, i, gsl_vector_get( xx, i ) );    
    gsl_vector_set( LvPointAvgSum, i, gsl_vector_get( yy, i ) ); }
  
  gsl_vector_free( xx );
  gsl_vector_free( yy ); }

double observer::calculateFlare( energyDissProc* obj, double nu ) {
  gsl_vector* temp_x = gsl_vector_alloc( obj->N );
  gsl_vector* temp_y = gsl_vector_alloc( obj->N );
  
  gsl_vector* rad = gsl_vector_alloc( r->getMaxIndex() );
  gsl_vector* flare = gsl_vector_alloc( r->getMaxIndex() );
  
  for( int i=0;i<r->getMaxIndex();i++ ) {
    gsl_vector_set_zero( temp_x );
    gsl_vector_set_zero( temp_y );
    
    r->update( i );
      
    gsl_vector_memcpy( temp_x, obj->vPoint );
    gsl_vector_memcpy( temp_y, obj->LvPoint );
    gsl_vector_mul( temp_y, temp_x );
    
    gsl_vector_set( rad, i, r->get() );
    gsl_vector_set( flare, i, interpolate( temp_x, temp_y, nu ) ); }
  
  std::string type = "FlarePoint_";
  type += obj->whoAmI( );
  bazinga::info(id,"Saving flare.");
  bazinga::save_GSLVectorFlare( type, r->getRadius_GSLVector( ), flare, nu, cfg->get<std::string>("output") );

  gsl_vector_free( temp_x );
  gsl_vector_free( temp_y );
  
  gsl_vector_free( rad );
  gsl_vector_free( flare ); }

double observer::interpolate( gsl_vector* x, gsl_vector* y, double _x ) {
  double interpolatedValue;
  double xmin, xmax = 0.0;
  
  if( _x < gsl_vector_min( x ) || _x > gsl_vector_max( x ) ) { interpolatedValue = 1.0; }
  
  for( int i=0; i<x->size-1; i++ ) {
    if( gsl_vector_get( x, i ) > _x ) {
      if( gsl_vector_get( y, i ) == 0 || gsl_vector_get( y, i+1 ) == 0 ) { interpolatedValue = 0.0; }
      else {
	interpolatedValue = log( gsl_vector_get( y, i ) ) + log(gsl_vector_get( y, i+1 )/gsl_vector_get( y, i ))/log(gsl_vector_get( x, i+1 )/gsl_vector_get( x, i ))*log(_x/gsl_vector_get( x, i ) );
      }
      interpolatedValue = exp(interpolatedValue);
      //      interpolatedValue = 0.5*( gsl_vector_get( y, i ) + gsl_vector_get( y, i+1 ) );
      break; }
  }	  

  return interpolatedValue; }

void observer::allocateLvPointSum( ) {
  vPointSum = gsl_vector_alloc(N);
  bazinga::print_GSLVector_allocated_memory( id,vPointSum );
  LvPointSum = gsl_vector_alloc(N);
  bazinga::print_GSLVector_allocated_memory( id,LvPointSum );  
  gsl_vector_set_zero( vPointSum );
  gsl_vector_set_zero( LvPointSum ); }

void observer::allocateLvPointAvgSum( ) {
  LvPointAvgSum = gsl_vector_alloc(N);
  bazinga::print_GSLVector_allocated_memory( id,LvPointAvgSum );
  gsl_vector_set_zero( LvPointAvgSum ); }

 void observer::freeLvPointSum( ) {
  gsl_vector_free( vPointSum );
  gsl_vector_free( LvPointSum ); }

void observer::freeLvPointAvgSum( ) {
  gsl_vector_free( LvPointAvgSum ); }

double observer::findr( double r, double theta ) {
  double theta_min = (thetaObs<thetaJ) ? 0.0 : (thetaObs-thetaJ);
  return r*(1.0-beta(Gamma)*cos(theta_min))/(1.0-beta(Gamma)*cos(theta)); }

 int observer::locate( gsl_vector* v, double val ) {
  int im;
  int il = 0;
  int iu = v->size;
  
  while(iu-il>1) {
      im=(iu+il) >> 1;
      if(val>=gsl_vector_get(v,im)) { il=im; }
      else { iu=im; }
  }
  return il; }

void observer::savePointLuminosity( energyDissProc* x ) {
  if( saveLumPoint && r->ifSaveRadius( ) ) {
    bazinga::info(id,"Saving point luminosities.");
    std::string type = "LvPoint_" + x -> whoAmI( );
    bazinga::save_GSLVector( type, x->vPoint, x->ele->gamma, x->LvPoint, r->getPosition( ), cfg->get<std::string>("output") ); }
}

void observer::savePointLuminositySum( ) {
  if( saveLumPoint && r->ifSaveRadius( ) ) {
    bazinga::info(id,"Saving summed point luminosities.");
    std::string type = "LvPoint_" + this -> whoAmI( );
    bazinga::save_GSLVector( type, vPointSum, LvPointSum, r->getPosition( ), cfg->get<std::string>("output") ); }
}

void observer::savePointLuminositySum( gsl_vector* gamma ) {
  if( saveLumPoint && r->ifSaveRadius( ) ) {
    bazinga::info(id,"Saving summed point luminosities.");
    std::string type = "LvPoint_" + this -> whoAmI( );
    bazinga::save_GSLVector( type, vPointSum, gamma, LvPointSum, r->getPosition( ), cfg->get<std::string>("output") ); }
}

void observer::saveAveragedPointLuminosity( energyDissProc* x ) {
  if( saveLumPointAvg ) { 
    bazinga::info(id,"Saving averaged point luminosities.");
    std::string type = "LvPointAvg_" + x->whoAmI( );
    bazinga::save_GSLVector( type, x->vPoint, x->ele->gamma, x->LvPointAvg, cfg->get<std::string>("output") ); }
}

void observer::saveAveragedPointLuminositySum( ) {
  if( saveLumPointAvg  ) {
    bazinga::info(id,"Saving summed averaged point luminosities.");
    std::string type = "LvPointAvg_" + this->whoAmI( );
    bazinga::save_GSLVector( type, vPointSum, LvPointAvgSum, cfg->get<std::string>("output") ); }
}

void observer::saveAveragedPointLuminositySum( gsl_vector* gamma ) {
  if( saveLumPointAvg  ) {
    bazinga::info(id,"Saving summed averaged point luminosities.");
    std::string type = "LvPointAvg_" + this->whoAmI( );
    bazinga::save_GSLVector( type, vPointSum, gamma, LvPointAvgSum, cfg->get<std::string>("output") ); }
}

int observer::ifCalculateFlares( ) {
  if( nu1 || nu2 || nu3 || nu4 || nu5 || nu6 ) { return 1; }
  else { return 0; }
}

//void observer::calculateJetLuminosity( energyDissProc* obj )
//{
//  print_info("+ Starting Jet ... ");
//
//  int NMI = 20;
//  double theta, Doppler;
//  
//  double dtheta = 2.0*(*(model->ThetaJ)/NMI);
//  
//  /* real stuff starts here */
//  double Lvp, robs, vpobs;
//  int p, q;
//  double u,t,dfi;
//  
//  for( int i=0;i<obj->N;i++ ){
//    Lvp = 0.0;
//    theta = *(model->ThetaObs)+*(model->ThetaJ);
//    do{
//      //      print_info(" theta",theta);
//      robs = findr( r->get(), theta );
//      //print_info("r",r->get());
//      //      print_info("robs",robs);
//      //      print_info("R0",*(model->R0));
//      //      print_info("RMax",*(model->RMax));
//      
//      if( robs > *(model->R0) && robs < *(model->RMax) ){
//	//	print_info("Promien OK");
//	Doppler = 1.0/(*(model->Gamma)*(1.0-cos(theta)*model->beta(*(model->Gamma))));
//	vpobs = obj->get_vJet(i)/Doppler;
//	p = locate( r->rad, robs );
//	if( vpobs > obj->get_vp(0) && vpobs < obj->get_vp(obj->N-1) )
//	  {
//	    q = locate( obj->vp, vpobs );
//	    t = (vpobs-obj->get_vp(q))/(obj->get_vp(q+1)-obj->get_vp(q));
//	    u = (robs-gsl_vector_get(r->rad,p))/(gsl_vector_get(r->rad,p+1)-gsl_vector_get(r->rad,p));
//	    Lvp = (1.0-t)*(1.0-u)*gsl_matrix_get(obj->Lpv,p,q)+t*(1.0-u)*gsl_matrix_get(obj->Lpv,p,q+1)+t*u*gsl_matrix_get(obj->Lpv,p+1,q+1)+(1.0-t)*u*gsl_matrix_get(obj->Lpv,p+1,q);
//
//	    //  print_info("Lvp",Lvp);
//
//	    if( theta<(*(model->ThetaObs)-(*(model->ThetaJ))) || *(model->ThetaObs)==0.0 )
//	      dfi=M_PI;
//	    else{
//	      if( fabs(((cos(*(model->ThetaJ))-cos(theta)*cos(*(model->ThetaObs)))/(sin(theta)*sin(*(model->ThetaObs)))))>1.0e0 )
//		dfi = 2.0*M_PI;
//	      else
//		dfi = 2.0e0*acos((cos(*(model->ThetaJ))-cos(theta)*cos(*(model->ThetaObs)))/(sin(theta)*sin(*(model->ThetaObs))));
//	    }
//	    obj->set_LvJet( i, obj->get_LvJet(i)+dfi*Lvp*pow(Doppler,3)*fabs(sin(theta)) );
//	  }
//      }
//      theta -= dtheta;
//    }
//    while ( theta>( *(model->ThetaObs)-*(model->ThetaJ)) );
//    obj->set_LvJet( i, obj->get_LvJet( i )*dtheta/(2.0*M_PI*(1.0-cos(*(model->ThetaJ)))) );  }
//}
//

