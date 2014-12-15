/**
    @file magneticField.cpp
    @author Mateusz Janiak
*/

#include "magneticField.hpp"

magneticField::magneticField( scfgp* _cfg, jetGeometry* _r, std::string _id ) : baseClass(_cfg, _r, _id ) { 

  /* request parameters */
  cfg -> request<double>("B0",0.0,&B0);
  cfg -> request<double>("B1",1.0,&B1);
  cfg -> request<double>("sigmaB",0.01,&sigmaB);
  cfg -> request<double>("kB",2.0,&kB);
  cfg -> request<std::string>("magModel","blob",&magModel);
  cfg -> request<double>("Gamma",10.0,&Gamma);
  cfg -> request<double>("thetaJ",0.1,&thetaJ);
  cfg -> request<double>("injRm",2.0e17,&injRm);
  cfg -> request<double>("eDiss",0.1,&eDiss);
  cfg -> request<double>("eEle",0.1,&eEle);
  cfg -> request<double>("mBH",1.0,&mBH);
  cfg -> request<double>("eJet",0.3,&eJet);
  cfg -> request<double>("mDot",1.0,&mDot);

  cfg -> updateRequests( );

  /* set energetics and magnetic flux value if in 'steady' model */
  if( magModel == "steady" ) {
    double Ledd, Ljet;
    Ledd = 1.3e47*mBH;
    Ljet = 0.5*eJet*mDot*Ledd;
    Lb = sigmaB*(1.0-eDiss)*Ljet/(1.0+sigmaB);
    B0steady = sqrt(6.0*Lb/(LIGHT_SPEED*beta(Gamma)))/(thetaJ*Gamma); }
}

void magneticField::printInfo( ) {
  bazinga::info(id,"Info");
  bazinga::print_info(id,"Magnetic field model",magModel);
  if( magModel == "blob" ) {
      bazinga::print_info(id,"B0",B0);
      bazinga::print_info(id,"B1 @ InjRm",B1);
      bazinga::print_info(id,"index kB",kB); }
  
  if( magModel == "steady" ) {
      bazinga::print_info(id,"sigmaB",sigmaB);
      bazinga::print_info(id,"B @ injRm",B0steady/injRm);
      bazinga::print_info(id,"Magnetic field flux", Lb);
      bazinga::print_info(id,"thetaJ", thetaJ); } }
 
double magneticField::getB( ) { 
  if( magModel == "steady" ) { return B0steady/r->get( ); };
  if( magModel == "blob" ) { return B0+(B1*pow(injRm/r->get( ),0.5*kB)); }
}

double magneticField::getB( double _r ) { 
  if( magModel == "steady" ) { return B0steady/_r; };
  if( magModel == "blob" ) { return B0+(B1*pow(injRm/_r,0.5*kB)); }
}

double magneticField::getMaxB( ) {
  if( magModel == "steady" ) { return B0steady/r->getR0( ); }
  if( magModel == "blob" ) { return B0+(B1*pow( injRm/r->getR0( ),0.5*kB)); }
}

double magneticField::get_uB( ) { return( DSQR( getB( ) )/(8.0*M_PI) ); }
double magneticField::get_uB( double _r ) { return( DSQR( getB( _r ) )/(8.0*M_PI) ); }
