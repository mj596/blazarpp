/**
    @file    blazar.cpp
    Main blazar++ file

    @author Mateusz Janiak
*/

#include <scfgp.hpp>
#include <bazinga.hpp>
#include "blazar.hpp"
#include "baseClass.hpp"
#include "jetGeometry.hpp"
#include "electrons.hpp"
#include "magneticField.hpp"
#include "observer.hpp"
#include "synchrotron.hpp"
#include "plots.hpp"
#include "selfSynCompton.hpp"
#include "externalRadiation.hpp"
#include "BlrPlanar.hpp"
#include "DtPlanar.hpp"
#include "AccdPlanar.hpp"
#include "BlrSpherical.hpp"
#include "DtSpherical.hpp"
#include "BlrSphericalGu.hpp"
#include "DtSphericalGu.hpp"
#include "AccdSphericalGu.hpp"
#include "QuasarAccDisk.hpp"

/**   Main programme starts here */
int main( int argc, char *argv[] )
{
  bazinga::print_section( );
  bazinga::print_info("\t\tBlazar++");
  bazinga::print_section( );

  scfgp *cfg = new scfgp( ); /* Initialize config file reader */
  
  if( argc != 1 ) { cfg->addConfigFile( argv[1] ); }
  else {
    bazinga::error("main","No config file given. You forgot about it?");;
    exit( 0 ); }

  /* Gather information about processes used */
  getConfigSwitches( cfg );

  /* initialize jet geometry */
  jetGeometry* rJet = new jetGeometry( cfg );
  rJet -> printInfo( );

  /* initialize electrons*/
  electrons* ele = new electrons( cfg, rJet, "ele" );
  ele -> printInfo( );

  /* initialize magnetic field */
  magneticField* B = new magneticField( cfg, rJet, "mag" );
  B -> printInfo( );

  /* initialize observer to caluclate luminosities  */
  observer* obs = new observer( cfg, rJet, "obs" );
  obs -> printInfo( );

  /* gather information about processes used */
  safetyConfigSwitches( cfg, ele );

  /* initialize Synchroton */
  energyDissProc* Syn = cfg->get<int>("Syn") ? new synchrotron( cfg, rJet, ele, B, "syn" ) : NULL;
  
  if( cfg->get<int>("Syn") ) { Syn -> printInfo( ); }
  if( cfg->get<int>("Syn") ) { ele -> addEnDissProc( Syn ); }
  if( cfg->get<int>("lumSyn") ) { obs -> addEnDissProc( Syn ); }
  if( cfg->get<int>("lumSynPoint") ) { obs -> addEnDissProcPoint( Syn ); }

  /* initialize synchrotron self Compton */
  energyDissProc* Ssc = cfg->get<int>("Ssc") ? new selfSynCompton( cfg, rJet, ele, B, "ssc" ) : NULL;

  /* I probably should not be doing that dynamic_cast but I really have no better idea ... */
  if( cfg->get<int>("Ssc") ) { dynamic_cast<selfSynCompton*>(Ssc)->syn=Syn; }

  if( cfg->get<int>("Ssc") ) { Ssc -> printInfo( ); }
  if( cfg->get<int>("Ssc") ) { ele -> addEnDissProc( Ssc ); }
  if( cfg->get<int>("lumSsc") ) { obs -> addEnDissProc( Ssc ); }
  if( cfg->get<int>("lumSscPoint") ) { obs -> addEnDissProcPoint( Ssc ); }

  /* initialize External Radiation - Blr */
  energyDissProc* Ext1 = cfg->get<int>("Ext1") ? new externalRadiation( cfg, rJet, ele, "ext1" ) : NULL;  
  if( cfg->get<int>("Ext1") ) { Ext1 -> printInfo( ); }
  if( cfg->get<int>("Ext1") ) { ele->addEnDissProc( Ext1 ); }
  if( cfg->get<int>("lumExt1") ) { obs->addEnDissProc( Ext1 ); }
  if( cfg->get<int>("lumExt1Point") ) { obs->addEnDissProcPoint( Ext1 ); }

  /* initialize External Radiation - Dt */
  energyDissProc* Ext2 = cfg->get<int>("Ext2") ? new externalRadiation( cfg, rJet, ele, "ext2" ) : NULL;  
  if( cfg->get<int>("Ext2") ) { Ext2 -> printInfo( ); }
  if( cfg->get<int>("Ext2") ) { ele->addEnDissProc( Ext2 ); }
  if( cfg->get<int>("lumExt2") ) { obs->addEnDissProc( Ext2 ); }
  if( cfg->get<int>("lumExt2Point") ) { obs->addEnDissProcPoint( Ext2 ); }

  /* initialize External Radiation Planar Geometry  - BlrPl */
  energyDissProc* BlrPl = cfg->get<int>("BlrPl") ? new BlrPlanar( cfg, rJet, ele, "blrPl" ) : NULL;
  if( cfg->get<int>("BlrPl") ) { BlrPl -> printInfo( ); }
  if( cfg->get<int>("BlrPl") ) { ele -> addEnDissProc( BlrPl ); }
  if( cfg->get<int>("BlrPl") ) { obs -> addExtPl( BlrPl ); }
  if( cfg->get<int>("lumBlrPl") ) { obs -> addEnDissProc( BlrPl ); }
  if( cfg->get<int>("lumBlrPlPoint") ) { obs -> addEnDissProcPoint( BlrPl ); }

  /* initialize External Radiation Planar Geometry - DtPl */
  energyDissProc* DtPl = cfg->get<int>("DtPl") ? new DtPlanar( cfg, rJet, ele, "dtPl" ) : NULL;
  if( cfg->get<int>("DtPl") ) { DtPl -> printInfo( ); }
  if( cfg->get<int>("DtPl") ) { ele -> addEnDissProc( DtPl ); }
  if( cfg->get<int>("DtPl") ) { obs -> addExtPl( DtPl ); }
  if( cfg->get<int>("lumDtPl") ) { obs -> addEnDissProc( DtPl ); }
  if( cfg->get<int>("lumDtPlPoint") ) { obs -> addEnDissProcPoint( DtPl ); }

  /* initialize External Radiation Planar Geometry - Accd */
  energyDissProc* Accd = cfg->get<int>("Accd") ? new AccdPlanar( cfg, rJet, ele, "accd" ) : NULL;
  if( cfg->get<int>("Accd") ) { Accd -> printInfo( ); }
  if( cfg->get<int>("Accd") ) { ele -> addEnDissProc( Accd ); }
  if( cfg->get<int>("Accd") ) { obs -> addExtPl( Accd ); }
  if( cfg->get<int>("lumAccd") ) { obs -> addEnDissProc( Accd ); }
  if( cfg->get<int>("lumAccdPoint") ) { obs -> addEnDissProcPoint( Accd ); }

  /* initialize External Radiation Spherical Geometry  - BlrSp */
  energyDissProc* BlrSp = cfg->get<int>("BlrSp") ? new BlrSpherical( cfg, rJet, ele, "blrSp" ) : NULL;
  if( cfg->get<int>("BlrSp") ) { BlrSp -> printInfo( ); }
  if( cfg->get<int>("BlrSp") ) { ele -> addEnDissProc( BlrSp ); }
  if( cfg->get<int>("BlrSp") ) { obs -> addExtSp( BlrSp ); }
  if( cfg->get<int>("lumBlrSp") ) { obs -> addEnDissProc( BlrSp ); }
  if( cfg->get<int>("lumBlrSpPoint") ) { obs -> addEnDissProcPoint( BlrSp ); }

  /* initialize External Radiation Spherical Geometry - DtSp */
  energyDissProc* DtSp = cfg->get<int>("DtSp") ? new DtSpherical( cfg, rJet, ele, "dtSp" ) : NULL;
  if( cfg->get<int>("DtSp") ) { DtSp -> printInfo( ); }
  if( cfg->get<int>("DtSp") ) { ele -> addEnDissProc( DtSp ); }
  if( cfg->get<int>("DtSp") ) { obs -> addExtSp( DtSp ); }
  if( cfg->get<int>("lumDtSp") ) { obs -> addEnDissProc( DtSp ); }
  if( cfg->get<int>("lumDtSpPoint") ) { obs -> addEnDissProcPoint( DtSp ); }

  /* initialize External Radiation Gu - Blr */
  energyDissProc* BlrGu = cfg->get<int>("BlrGu") ? new BlrSphericalGu( cfg, rJet, ele, "blrGu" ) : NULL;  
  if( cfg->get<int>("BlrGu") ) { BlrGu -> printInfo( ); }
  if( cfg->get<int>("BlrGu") ) { ele->addEnDissProc( BlrGu ); }
  if( cfg->get<int>("BlrGu") ) { obs -> addExtGu( BlrGu ); }
  if( cfg->get<int>("lumBlrGu") ) { obs->addEnDissProc( BlrGu ); }
  if( cfg->get<int>("lumBlrGuPoint") ) { obs->addEnDissProcPoint( BlrGu ); }

  /* initialize External Radiation Gu - Dt */
  energyDissProc* DtGu = cfg->get<int>("DtGu") ? new DtSphericalGu( cfg, rJet, ele, "dtGu" ) : NULL;  
  if( cfg->get<int>("DtGu") ) { DtGu -> printInfo( ); }
  if( cfg->get<int>("DtGu") ) { ele->addEnDissProc( DtGu ); }
  if( cfg->get<int>("DtGu") ) { obs -> addExtGu( DtGu ); }
  if( cfg->get<int>("lumDtGu") ) { obs->addEnDissProc( DtGu ); }
  if( cfg->get<int>("lumDtGuPoint") ) { obs->addEnDissProcPoint( DtGu ); }

  /* initialize External Radiation Gu - Accd */
  energyDissProc* AccdGu = cfg->get<int>("AccdGu") ? new AccdSphericalGu( cfg, rJet, ele, "accdGu" ) : NULL;  
  if( cfg->get<int>("AccdGu") ) { AccdGu -> printInfo( ); }
  if( cfg->get<int>("AccdGu") ) { ele->addEnDissProc( AccdGu ); }
  if( cfg->get<int>("AccdGu") ) { obs -> addExtGu( AccdGu ); }
  if( cfg->get<int>("lumAccdGu") ) { obs->addEnDissProc( AccdGu ); }
  if( cfg->get<int>("lumAccdGuPoint") ) { obs->addEnDissProcPoint( AccdGu ); }

  /* initialize Quasar Accretion Disk */
  QuasarAccDisk* QAccd = cfg->get<int>("QAccd") ? new QuasarAccDisk( cfg, "QAccd" ) : NULL;  
  if( cfg->get<int>("QAccd") ) { QAccd -> printInfo( ); }

  /* read QAccd template data */
  if( cfg->get<int>("QAccd") ) { QAccd -> readTemplate( ); }
  if( cfg->get<int>("QAccd") ) { QAccd -> scaleTemplate( ); }
  if( cfg->get<int>("QAccd") ) { QAccd -> saveTemplate( ); }

  bazinga::print_section( );
  /* Before proceeding further show us what you're going to do */
  ele -> listEnDissProc( );
  
  /* Before proceeding further show us what luminosities are you going to calculate */
  obs -> listEnDissProc( );
  obs -> listEnDissProcPoint( );
  obs -> listExtPl( );
  obs -> listExtSp( );
  obs -> listExtGu( );
  
  /* ---------------------------------------------------------------------------------------------- */
  /* BEGIN root Canvases and Graphs */
  /* ---------------------------------------------------------------------------------------------- */

#ifdef USE_ROOT
/* Create TApplication */
TApplication* theApp = cfg->get<int>("root") ? new TApplication("App",0,0) : NULL;

  /* Create all the necessary TCanvas */
  /* First row - electrons*/
TCanvas *eleInjCanvas = cfg->get<int>("root") && ele->ifEvol() ? new TCanvas("EleInj","Electron injection",0,0,400,400) : NULL;   
  TCanvas *eleCanvas = cfg->get<int>("root") ? new TCanvas("EleEvo","Electron evolution",400,0,400,400) : NULL;
  /* Second row - luminosities */
  TCanvas *intLumCanvas = cfg->get<int>("root") && obs->sizeEnDissProc( ) ? new TCanvas("IntLum","Intrinsic luminosities",0,430,400,400) : NULL;
  TCanvas *lumPointCanvas = cfg->get<int>("root") && obs->sizeEnDissProcPoint( ) ? new TCanvas("LumPoint","Point luminosities",400,430,400,400) : NULL;
  TCanvas *lumPointAvgCanvas = cfg->get<int>("root") && obs->sizeEnDissProcPoint( ) ? new TCanvas("LumPointAvg","Avg Point luminosities",800,430,400,400) : NULL;

  gSystem -> ProcessEvents();

  /* Create all the necessary TMultiGraph */
  TMultiGraph* mgInjNgamma = cfg->get<int>("root") && ele->ifEvol( ) ? new TMultiGraph("EleInj","Electron injection;log #gamma;log #gamma^{2} N_{#gamma}") : NULL;
  TMultiGraph* mgNgamma = cfg->get<int>("root") ? new TMultiGraph("EleEvo","Electron evolution;log #gamma;log #gamma^{2} N_{#gamma}") : NULL;
  TMultiGraph* mgIntLum = cfg->get<int>("root") && obs->sizeEnDissProc( ) ? new TMultiGraph("IntLum","Intrinsic luminosities;log #nu;log #nu L_{#nu}") : NULL;
  TMultiGraph* mgLumPoint = cfg->get<int>("root") && obs->sizeEnDissProcPoint( ) ? new TMultiGraph("LumPoint","Point luminosities;log #nu;log #nu L_{#nu}") : NULL;
  TMultiGraph* mgLumAvgPoint = cfg->get<int>("root") && obs->sizeEnDissProcPoint( ) ? new TMultiGraph("LumPointAvg","Avg Point luminosities;log #nu;log #nu L_{#nu}") : NULL;
#endif

  /* ---------------------------------------------------------------------------------------------- */
  /* END root Canvases and Graphs */
  /* ---------------------------------------------------------------------------------------------- */
  
  bazinga::print_section( );

  if( ele->ifEvol( ) ) {
    bazinga::info( "main", "Starting electron evolution loop ...");
    bazinga::print_section( ); }
  else {
    bazinga::info("main", "No electron evolution. Injection only.");
    bazinga::print_section( ); }

  /* ---------------------------------------------------------------------------------------------- */
  /* MAIN LOOP STARTS HERE */
  /* ---------------------------------------------------------------------------------------------- */

  for( int i=0;i<rJet->getMaxIndex( );i++ ) {
    rJet->update( i );
    rJet->show( );
    /* inject electrons */
    ele -> inject( );
    
    /* update the data about processes that will cool electrons */
    for( int i=0; i<ele->EnDissProc.size(); i++ ) {
      bazinga::info( ele->EnDissProc[i]->whoAmI( ), "Updating" );
      ele -> EnDissProc[i] -> update( ); }
    
    /* show some info about magnetic field */
    bazinga::print_info( B->whoAmI(), "Magnetic field", B->getB(), "Gauss");

    /* electrons evolution taking part here*/
    if( ele -> ifEvol( ) ) {
      /* save injected electrons */
      ele -> saveInjection( );
      
#ifdef USE_ROOT
      /* make ROOT plot */
      if( cfg->get<int>("root") && ele -> ifEvol( ) ) { makeNgammaPlotROOT( ele->gamma, ele->S, eleInjCanvas, mgInjNgamma, ele->whoAmI(), rJet ); }
#endif
      
      /* solve evolution equations for electrons */
      ele -> evolve( ); }
	  
    /* save electrons */
    ele -> saveNgamma( );
    /* calculate AvgNgamma; this is just an update with the newest Ngamma so is has to be run everytime in a loop */
    ele -> avgNgamma();
    
#ifdef USE_ROOT
    /* make ROOT plot */
    if( cfg->get<int>("root") ) { makeNgammaPlotROOT( ele->gamma, ele->Ngamma, eleCanvas, mgNgamma, ele->whoAmI(), rJet ); }
#endif
    
    /* ---------------------------------------------------------------------------------------------- */
    /* BEGIN Calculate intrinsic luminosities */
    /* ---------------------------------------------------------------------------------------------- */
    /* take electrons from the previous step - either evolved, only injected or read from a file */
    for( int j=0; j<obs->sizeEnDissProc(); j++ ) {
      bazinga::info( obs -> EnDissProc[j ]-> whoAmI( ), "Calculating intrinsic luminosity." );
      /* set luminosities */
      obs -> EnDissProc[j] -> setLpv( );
      obs -> EnDissProc[j] -> saveLuminosity( );

#ifdef USE_ROOT
      if( cfg->get<int>("root") ) { makeLPlotROOT( obs->EnDissProc[j]->vp, obs->EnDissProc[j]->Lpv, intLumCanvas, mgIntLum, obs->EnDissProc[j]->whoAmI( ), rJet ); }
#endif
    }
    /* ---------------------------------------------------------------------------------------------- */
    /* END Calculate intrinsic luminosities */
    /* ---------------------------------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------------------------------- */
    /* BEGIN Calculate point luminosities */
    /* ---------------------------------------------------------------------------------------------- */
    for( int j=0; j<obs->sizeEnDissProcPoint(); j++ ) {
      obs -> calculatePointLuminosity( obs->EnDissProcPoint[j] );
      obs -> savePointLuminosity( obs->EnDissProcPoint[j] );

#ifdef USE_ROOT
      if( cfg->get<int>("root") ) { makeLPlotROOT( obs->EnDissProcPoint[j]->vPoint, obs->EnDissProcPoint[j]->LvPoint, lumPointCanvas, mgLumPoint, obs->EnDissProcPoint[j]->whoAmI( ), rJet ); }
#endif
      
      /* Averaging spectra */
      obs -> avgPointLuminosity( obs->EnDissProcPoint[j] ); }
    
    /* Calculate accretion disk luminosity - QAccd */
    if( cfg->get<int>("QAccd") ) { QAccd -> calculateLuminosity( ); }
    
    /* ---------------------------------------------------------------------------------------------- */
    /* END Calculate point luminosities */
    /* ---------------------------------------------------------------------------------------------- */
    
    /* ---------------------------------------------------------------------------------------------- */
    /* BEGIN Sum point luminosities */
    /* ---------------------------------------------------------------------------------------------- */
    if( obs->sizeEnDissProcPoint( ) ) {
      obs -> sumPointLuminosity( );
      obs -> savePointLuminositySum( ele->gamma );

#ifdef USE_ROOT
      if( cfg->get<int>("root") ) { makeLSumPlotROOT( obs->vPointSum, obs->LvPointSum, lumPointCanvas, mgLumPoint, obs->whoAmI( ), rJet ); }
#endif
    }
    /* ---------------------------------------------------------------------------------------------- */
    /* END Sum point luminosities */
    /* ---------------------------------------------------------------------------------------------- */  

    /* save information on upe vs radius */
    for( int j=0; j < obs->sizeEnDissProc( ); j++ ) { obs ->  EnDissProc[j] -> saveUpeR( ); }

    /* save information on Rm vs radius for ExtPl */
    for( int j=0; j < obs->sizeExtPl( ); j++ ) { dynamic_cast<externalRadiationPlanar*>( obs->ExtPl[j] ) -> saveRmVsR( ); }
  }
  /* ---------------------------------------------------------------------------------------------- */
  /* MAIN LOOP ENDS HERE */
  /* ---------------------------------------------------------------------------------------------- */
  
  if( ele->ifEvol( ) ) {
    bazinga::info( ele->whoAmI( ), "Ending electron evolution loop.");
    bazinga::print_section( ); }
  else {
    bazinga::info( ele->whoAmI( ), "Ending electron loop.");
    bazinga::print_section( ); }

  /* ---------------------------------------------------------------------------------------------- */
  
  /* save information about NgammaAvg */
  ele -> saveNgammaAvg( );
  
  /* ---------------------------------------------------------------------------------------------- */
  /* BEGIN Save average spectra */
  /* ---------------------------------------------------------------------------------------------- */
  for( int j=0; j < obs -> sizeEnDissProcPoint(); j++ ) {
    obs -> saveAveragedPointLuminosity( obs -> EnDissProcPoint[j] );
    
#ifdef USE_ROOT
  if( cfg->get<int>("root") ) { makeLPlotROOT( obs->EnDissProc[j]->vPoint, obs->EnDissProc[j]->LvPointAvg, lumPointAvgCanvas, mgLumAvgPoint, obs->EnDissProcPoint[j]->whoAmI( ), rJet ); }
#endif
  }

  /* ---------------------------------------------------------------------------------------------- */
  /* END Save average spectra */
  /* ---------------------------------------------------------------------------------------------- */

  /* ---------------------------------------------------------------------------------------------- */
  /* BEGIN Sum and save averaged point sum spectra */
  /* ---------------------------------------------------------------------------------------------- */
  
    if( obs -> sizeEnDissProcPoint( ) ) {
      obs -> sumPointAvgLuminosity( );
      
      if( cfg->get<int>("QAccd") ) { obs -> addQAccdLuminosity( QAccd ); }

      obs -> saveAveragedPointLuminositySum( ele->gamma );
      
#ifdef USE_ROOT
    if( cfg->get<int>("root") ) { makeLSumPlotROOT( obs->vPointSum, obs->LvPointAvgSum, lumPointAvgCanvas, mgLumAvgPoint, obs->whoAmI( ), rJet ); }
#endif
    }
    
  /* ---------------------------------------------------------------------------------------------- */
  /* END Sum and save averaged sum spectra */
  /* ---------------------------------------------------------------------------------------------- */

  /* ---------------------------------------------------------------------------------------------- */
  /* BEGIN Calculating flares */
  /* ---------------------------------------------------------------------------------------------- */
    obs -> calculateAllFlares( );
  /* ---------------------------------------------------------------------------------------------- */
  /* END Calculating flares */
  /* ---------------------------------------------------------------------------------------------- */

  /* ---------------------------------------------------------------------------------------------- */
  /* BEGIN make gnuplot plots */
  /* ---------------------------------------------------------------------------------------------- */
    bazinga::print_section( );
    bazinga::info("main","Preparing gnuplot files.");
    if( cfg->get<int>("saveElectrons") ) { makeNgammaPlots( ele, cfg, rJet ); } /* electron plot */
    if( cfg->get<int>("saveLumPoint") ) { makeLPointPlots( obs, cfg, rJet ); } /* point source luminosities */
    if( cfg->get<int>("saveLumPoint") ) { makeLPointVsElePlots( obs, cfg, rJet ); } /* point source luminosities */
    if( cfg->get<int>("saveLumPoint") ) { makeLPointSumPlots( obs, cfg, rJet ); } /* point source luminosity sum */
    if( obs->sizeEnDissProc( ) ) { makeUpeRPlots( obs, cfg, rJet ); } /* make Upe vs radius plots */
    if( cfg->get<int>("QAccd") ) { makeLQAccdPlots( QAccd, cfg ); } /* make QAccd plots -  template and accretion disk multicolor blackbody */

    if( obs->sizeExtPl( ) ) {
      makeRmPlots( obs, cfg, rJet ); /* make Rm vs radius plots */
      makeExtPldLdRPlots( obs, cfg, rJet ); } /* make dLdR vs R plots for ExtPl */
 
    /* flare plots */
    if( obs -> ifCalculateFlares( ) ) { makePointFlare( obs, cfg, rJet ); }
    
    bazinga::info("main","Gnuplot files ready.");
  /* ---------------------------------------------------------------------------------------------- */
  /* END make gnuplot plots */
  /* ---------------------------------------------------------------------------------------------- */

#ifdef USE_ROOT
    if( cfg->get<int>("root") && cfg->get<int>("printRoot") ) {
    bazinga::print_section( );
    bazinga::info("main","Priniting ROOT canvases");
    eleCanvas -> Print();
    if( ele->ifEvol( ) ) { eleInjCanvas -> Print(); }
    if( obs->sizeEnDissProc( ) ) { intLumCanvas -> Print(); }
    if( obs->sizeEnDissProcPoint( ) ) { lumPointCanvas -> Print(); }
    if( obs->sizeEnDissProcPoint( ) ) { lumPointAvgCanvas -> Print(); }
  }
#endif

  bazinga::print_section( );
    /* print KN info */
    /* go along the jet once again */
  for( int i=0;i<rJet->getMaxIndex();i++ ) {
    rJet -> update( i );
    for( int j=0; j < ele->EnDissProc.size(); j++ ) { ele->EnDissProc[j]->print_KN_info( ); }
  }

bazinga::print_section( );
bazinga::info("main","Merci beaucoup. I'm done here.");

#ifdef USE_ROOT
 if( cfg->get<int>("root") ) { bazinga::info("main","PRESS 'ENTER' TO QUIT.");
  std::cin.ignore(); }
#endif

  /* free memory */
 bazinga::info("BLAZAR","Let's free some memory:");
  for( int j=0; j<ele->EnDissProc.size(); j++ ) { delete ele->EnDissProc[j]; }
  delete obs;
  delete ele;
  delete rJet;
  delete cfg;
  delete QAccd;
#ifdef USE_ROOT
  delete theApp;
#endif
return 0;
}

