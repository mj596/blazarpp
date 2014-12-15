#include "plots.hpp"

#ifdef USE_ROOT
double MIN_ROOT=1.0e60;
double MAX_ROOT=0;
#endif

double MultiplyPeak = 5;
double DividePeak = 10;

#ifdef USE_ROOT
int checkJetGeometryToRootPlot( jetGeometry *r )
{
  int plot = 0;

  int rMaxId = r->getMaxIndex( );
  int rCurrentId = r->getIndex( );
  int maxRootDisplay = 5;
  double drRoot = (rMaxId-2)/((double)maxRootDisplay-2);

  if( rMaxId < maxRootDisplay )
      plot = 1;
  else
    {
      if( rCurrentId == 0 || rCurrentId == rMaxId-1 )
	plot = 1;

      for( int i=1; i<maxRootDisplay-1; i++ )
	if( rCurrentId == floor(i*drRoot) )
	  plot = 1;
    }

  return plot;
}

void makeLSumPlotROOT( gsl_vector* x, gsl_vector* y, TCanvas* c1, TMultiGraph* mg, std::string type, jetGeometry* r )
{
  if( checkJetGeometryToRootPlot( r ) )
    {
      c1->cd();
 
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      gPad->SetLogy(1);
      gPad->SetLogx(1);
      gPad->Modified();
      gPad->Update();
      
      int n = y->size-1;
      
      double max = 0;
      double min = 0;

      double xx[n], yy[n];
      for (int i=0;i<n;i++)
	{
	  xx[i] = gsl_vector_get(x,i);
	  yy[i] = gsl_vector_get(y,i);
	  if( yy[i] > max )
	      max = yy[i];
	}

      min = max/DividePeak;
      max *= MultiplyPeak;

      TGraph *gr = new TGraph(n,xx,yy);
      gr->SetMarkerStyle(6);
      gr->SetLineColor(chooseColorRoot(type));

      if( max > MAX_ROOT ) MAX_ROOT = max;
      if( min < MIN_ROOT ) MIN_ROOT = min;

      mg->SetMinimum(MIN_ROOT);
      mg->SetMaximum(MAX_ROOT);

      mg->Add( gr );
      c1->Clear();
      mg->Draw("al");
      mg->GetXaxis()->SetLimits(1.0e05,1.0e28);

      c1->Modified();
      c1->Update();
      gSystem->ProcessEvents();
    }
}

void makeNgammaPlotROOT( gsl_vector* x, gsl_vector* y, TCanvas* c1, TMultiGraph* mg, std::string type, jetGeometry* r )
{
  if( checkJetGeometryToRootPlot( r ) )
    {
      c1->cd();

      gPad->SetGridx(1);
      gPad->SetGridy(1);
      gPad->SetLogy(1);
      gPad->SetLogx(1);
      gPad->Modified();
      gPad->Update();
      
      int n = y->size-1;
      
      double max = 0;

      double xx[n], yy[n];
      for (int i=0;i<n;i++)
	{
	  xx[i] = gsl_vector_get(x,i+1);
	  yy[i] = gsl_vector_get(y,i)*pow(xx[i],2.0);
	  if( yy[i] > max )
	    max = yy[i];
	}
      
      TGraph *gr = new TGraph(n,xx,yy);
      gr->SetMarkerStyle(6);
      gr->SetLineColor(chooseColorRoot(type)); 
      mg->SetMaximum(max*MultiplyPeak);
      mg->Add( gr );

      c1->Clear();
      mg->Draw("al");
      
      c1->Modified();
      c1->Update();
      gSystem->ProcessEvents();
    }
}

void makeLPlotROOT( gsl_vector* x, gsl_vector* y, TCanvas* c1, TMultiGraph* mg, std::string type, jetGeometry* r )
{
  if( checkJetGeometryToRootPlot( r ) )
    {
      c1->cd();
 
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      gPad->SetLogy(1);
      gPad->SetLogx(1);
      gPad->Modified();
      gPad->Update();
      
      int n = y->size-1;
      
      double max = 0;
      double min = 0;

      double xx[n], yy[n];
      for (int i=0;i<n;i++)
	{
	  xx[i] = gsl_vector_get(x,i);
	  yy[i] = gsl_vector_get(y,i)*xx[i];
	  if( yy[i] > max )
	    max = yy[i];
	}

      min = max/DividePeak;
      max *= MultiplyPeak;

      TGraph *gr = new TGraph(n,xx,yy);
      gr->SetMarkerStyle(6);
      gr->SetLineColor(chooseColorRoot(type));

      if( max > MAX_ROOT ) MAX_ROOT = max;
      if( min < MIN_ROOT ) MIN_ROOT = min;

      mg->SetMinimum(MIN_ROOT);
      mg->SetMaximum(MAX_ROOT);

      mg->Add( gr );
      c1->Clear();
      mg->Draw("al");
      mg->GetXaxis()->SetLimits(1.0e05,1.0e28);

      c1->Modified();
      c1->Update();
      gSystem->ProcessEvents();
    }
}


Color_t chooseColorRoot( std::string type ) {
  /* see http://root.cern.ch/root/html/TColor.html */
  Color_t color;
  color = kBlack;
  if( type == "ele" ) color = kBlack;
  if( type == "syn" ) color = kBlue+1;
  if( type == "ssc" ) color = kCyan+1;
  if( type == "ext1" ) color = kYellow-6;
  if( type == "ext2" ) color = kOrange+8;
  if( type == "blrPl" ) color = kGreen-8;
  if( type == "dtPl" ) color = kRed-8;
  if( type == "accd" ) color = kMagenta;
  if( type == "blrSp" ) color = kGreen+2;
  if( type == "dtSp" ) color = kRed+2;
  return color; }
#endif

std::string chooseLineGnuplot( std::string type ){
  std::string color;
  color = "lc rgb\"#000000\"";
  if( type == "QAccd" ) color = "lt 1 lw 1 lc rgb\"#26A32A\"";
  if( type == "ele" ) color = "lc rgb\"#000000\"";
  if( type == "syn" ) color = "lt 1 lw 1 lc rgb\"#26A32A\"";
  if( type == "ssc" ) color = "lt 1 lw 1 lc rgb\"#49CDF5\""; /* Light-Blue */
  if( type == "ext1" ) color = "lt 1 lw 1 lc rgb\"#FFB145\""; /* Orange */
  if( type == "ext2" ) color = "lt 3 lw 1 lc rgb\"#FFB145\"";
  if( type == "blrPl" ) color = "lt 1 lw 1 lc rgb\"#BA2727\""; /* Red */
  if( type == "dtPl" ) color = "lt 3 lw 1 lc rgb\"#BA2727\"";
  if( type == "accd" ) color = "lt 1 lw 1 lc rgb\"#BD35B4\""; /* Purple */
  if( type == "blrSp" ) color = "lt 1 lw 1 lc rgb\"#1F2BCF\""; /* Blue */
  if( type == "dtSp" ) color = "lt 3 lw 1 lc rgb\"#1F2BCF\"";
  return color; }

void makeLQAccdPlots( QuasarAccDisk* obj, scfgp* _cfg )
{
  double nuLnuMaxValue;
  std::ofstream gnuplotFile;
  std::stringstream filename;
  std::string dir = _cfg->get<std::string>("output");
  
  bazinga::print_info("plt","Plot",obj->whoAmI());
  filename.str("");
  filename << dir << "/" << obj->whoAmI() << ".gp";
  gnuplotFile.open( filename.str().c_str() );

  /* determine roughly maximum nuLnu value to make nice looking plots */
  gsl_vector* temp_x = gsl_vector_alloc( obj->getN( ) );
  gsl_vector* temp_y = gsl_vector_alloc( obj->getN( ) );
  gsl_vector_memcpy( temp_x, obj->get_v( ) );
  gsl_vector_memcpy( temp_y, obj->get_Lv( ) );
  gsl_vector_mul( temp_y, temp_x );
  nuLnuMaxValue = gsl_vector_max( temp_y );

  /* set up gnuplot settings */
  gnuplotFile << "reset" << std::endl;
  gnuplotFile << "set terminal postscript enhanced color" << std::endl;
  gnuplotFile << "set output \"" << obj->whoAmI() << ".eps\"" << std::endl;
  gnuplotFile << "set size square" << std::endl;
  gnuplotFile << "set border lw 2" << std::endl;
  gnuplotFile << "set grid" << std::endl;
  gnuplotFile << "set title \"" << obj->whoAmI() << " spectra\" " << std::endl;
  gnuplotFile << "set xlabel \"{/Symbol n}\" " << std::endl;
  gnuplotFile << "set ylabel \"{/Symbol n} L_{/Symbol n}\"" << std::endl;
  gnuplotFile << "set logscale xy" << std::endl;
  gnuplotFile << "set yrange [" << nuLnuMaxValue/yLowerScale << ":" << nuLnuMaxValue*yUpperScale << "]" << std::endl;
  gnuplotFile << "plot\\" << std::endl;
  gnuplotFile << "\"Lv_" << obj->whoAmI() << "\" using ($1):($1*$2) w l lw 3 notitle,\\" << std::endl;
  gnuplotFile << "\"LvTemplate_" << obj->whoAmI() << "\" using ($1):($1*$2) w l lw 3 notitle,\\" << std::endl;

  for( int k=1;k<obj->R->getMaxIndex();k++ ) 
    {
      obj->R->update( k );
      if( obj->get_save_Lr( ) ) {
	gnuplotFile << "\"Lv_" << obj->whoAmI()<< "_" << obj->R->getPosition( ) << "\" using ($1):($1*$2) lc rgb \"grey\" w l notitle";
	if( k != obj->R->getMaxIndex()-1 ) 
	  gnuplotFile << ",\\" << std::endl; }
    }
      
  gnuplotFile.close();
}

void makeNgammaPlots( electrons* obj, scfgp* _cfg, jetGeometry* _r )
{
  std::ofstream gnuplotFile;
  std::stringstream filename;
  std::string dir = _cfg->get<std::string>("output");

  bazinga::print_info("plt","Plot",obj->whoAmI());

  filename << dir << "/" << "Ngamma.gp";
  gnuplotFile.open( filename.str().c_str() );
  
  /* set up gnuplot settings */
  gnuplotFile << "reset" << std::endl;
  gnuplotFile << "set terminal postscript enhanced color" << std::endl;
  gnuplotFile << "set output \"Ngamma.eps\"" << std::endl;
  gnuplotFile << "set size square" << std::endl;
  gnuplotFile << "set border lw 2" << std::endl;
  gnuplotFile << "set grid" << std::endl;
  gnuplotFile << "set title \"Electron evolution\" " << std::endl;
  gnuplotFile << "set xlabel \"{/Symbol g} \" " << std::endl;
  gnuplotFile << "set ylabel \"{/Symbol g}^2 N_{/Symbol g} \"" << std::endl;
  gnuplotFile << "set logscale xy" << std::endl;
  gnuplotFile << "plot\\" << std::endl;

  gnuplotFile << "\"Ngamma_Avg\" using ($1):($1**2*$2) w l lw 2 title \"Averaged\",\\" << std::endl;

  for( int i=1;i<_r->getMaxIndex();i++ )
    {
      _r->update( i );
      if( _r->ifSaveRadius( ) ) { 
	gnuplotFile << "\"Ngamma_" << _r->getPosition( ) << "\" using ($1):($1**2*$2) w l notitle";
	if( i != _r->getMaxIndex()-1 ) 
	  gnuplotFile << ",\\" << std::endl; }
    }
  
  gnuplotFile.close();
}

void makeLPointPlots( observer* obs, scfgp* _cfg, jetGeometry* _r )
{
  double nuLnuMaxValue;
  std::ofstream gnuplotFile;
  std::stringstream filename;
  std::string dir = _cfg->get<std::string>("output");

  for( int i=0; i<obs->EnDissProcPoint.size(); i++ )
    {
      bazinga::print_info("plt","Plot (point source)",obs->EnDissProcPoint[i]->whoAmI());
      filename.str("");
      filename << dir << "/" << obs->EnDissProcPoint[i]->whoAmI() << "Point.gp";
      gnuplotFile.open( filename.str().c_str() );

       /* determine roughly maximum nuLnu value to make nice looking plots */
      gsl_vector* temp_x = gsl_vector_alloc( obs->EnDissProcPoint[i]->N );
      gsl_vector* temp_y = gsl_vector_alloc( obs->EnDissProcPoint[i]->N );
      gsl_vector_memcpy( temp_x, obs->EnDissProcPoint[i]->vPoint );
      gsl_vector_memcpy( temp_y, obs->EnDissProcPoint[i]->LvPointAvg );
      gsl_vector_mul( temp_y, temp_x );
      nuLnuMaxValue = gsl_vector_max( temp_y );

      /* set up gnuplot settings */
      gnuplotFile << "reset" << std::endl;
      gnuplotFile << "set terminal postscript enhanced color" << std::endl;
      gnuplotFile << "set output \"" << obs->EnDissProcPoint[i]->whoAmI() << "Point.eps\"" << std::endl;
      gnuplotFile << "set size square" << std::endl;
      gnuplotFile << "set border lw 2" << std::endl;
      gnuplotFile << "set grid" << std::endl;
      gnuplotFile << "set title \"" << obs->EnDissProcPoint[i]->whoAmI() << " spectra evolution (POINT SOURCE)\" " << std::endl;
      gnuplotFile << "set xlabel \"{/Symbol n}\" " << std::endl;
      gnuplotFile << "set ylabel \"{/Symbol n} L_{/Symbol n}\"" << std::endl;
      gnuplotFile << "set logscale xy" << std::endl;
      gnuplotFile << "set yrange [" << nuLnuMaxValue/yLowerScale << ":" << nuLnuMaxValue*yUpperScale << "]" << std::endl;
      gnuplotFile << "plot\\" << std::endl;
      gnuplotFile << "\"LvPointAvg_" << obs->EnDissProcPoint[i]->whoAmI() << "\" using ($1):($1*$3) w l lw 3 notitle,\\" << std::endl;

      for( int k=1;k<_r->getMaxIndex();k++ ) 
	{
	  _r->update( k );
	  if( _r->ifSaveRadius( ) ) { 
	    gnuplotFile << "\"LvPoint_" << obs->EnDissProcPoint[i]->whoAmI()<< "_" << _r->getPosition( ) << "\" using ($1):($1*$3) " << chooseLineGnuplot( obs->EnDissProcPoint[i]->whoAmI() ) << " w l notitle";
	    if( k != _r->getMaxIndex()-1 ) 
	      gnuplotFile << ",\\" << std::endl; }
	}
      
      gnuplotFile.close();
    }
}

void makeLPointSumPlots( observer* obj, scfgp* _cfg, jetGeometry* _r )
{
  double nuLnuMaxValue;
  std::ofstream gnuplotFile;
  std::stringstream filename;
  std::string dir = _cfg->get<std::string>("output");

  bazinga::print_info("plt","Plot (point source sum)",obj->whoAmI());

  if( obj->sizeEnDissProcPoint() )
    {
      filename << dir << "/" << obj->whoAmI() << "Point.gp";
      gnuplotFile.open( filename.str().c_str() );
    }
  
  /* determine roughly maximum nuLnu value to make nice looking plots */
  gsl_vector* temp_x = gsl_vector_alloc( obj->N );
  gsl_vector* temp_y = gsl_vector_alloc( obj->N );
  gsl_vector_memcpy( temp_x, obj->vPointSum );
  gsl_vector_memcpy( temp_y, obj->LvPointAvgSum );
  nuLnuMaxValue = gsl_vector_max( temp_y );
  
  /* set up gnuplot settings */
  gnuplotFile << "reset" << std::endl;
  gnuplotFile << "set terminal postscript enhanced color" << std::endl;
  gnuplotFile << "set output \"" << obj->whoAmI() << "Point.eps\"" << std::endl;
  gnuplotFile << "set size square" << std::endl;
  gnuplotFile << "set border lw 2" << std::endl;
  gnuplotFile << "set grid" << std::endl;
  gnuplotFile << "set title \"" << obj->whoAmI() << " spectra evolution\" " << std::endl;
  gnuplotFile << "set xlabel \"{/Symbol n}\" " << std::endl;
  gnuplotFile << "set ylabel \"{/Symbol n} L_{/Symbol n}\"" << std::endl;
  gnuplotFile << "set logscale xy" << std::endl;
  gnuplotFile << "set yrange [" << nuLnuMaxValue/yLowerScale << ":" << nuLnuMaxValue*yUpperScale << "]" << std::endl;
  gnuplotFile << "plot\\" << std::endl;
  gnuplotFile << "\"LvPointAvg_" << obj->whoAmI() << "\" using ($1):($3) w l lw 3 notitle,\\" << std::endl;

  for( int i=1;i<_r->getMaxIndex();i++ ) 
    {
      _r->update( i );
      if( _r->ifSaveRadius( ) ) { 
	gnuplotFile << "\"LvPoint_" << obj->whoAmI()<< "_" << _r->getPosition( ) << "\" using ($1):($3) " << chooseLineGnuplot( obj->whoAmI() ) << " w l notitle";
	if( i != _r->getMaxIndex()-1 ) 
	  gnuplotFile << ",\\" << std::endl; }
    }
  
  gnuplotFile.close();
}

void makeLPointVsElePlots( observer* obs, scfgp* _cfg, jetGeometry* _r )
{
  double nuLnuMaxValue;
  std::ofstream gnuplotFile;
  std::stringstream filename;
  std::string dir = _cfg->get<std::string>("output");

  for( int i=0; i<obs->EnDissProcPoint.size(); i++ )
    {
      bazinga::print_info("plt","Plot (point source)",obs->EnDissProcPoint[i]->whoAmI());
      filename.str("");
      filename << dir << "/" << obs->EnDissProcPoint[i]->whoAmI() << "PointVsEle.gp";
      gnuplotFile.open( filename.str().c_str() );

       /* determine roughly maximum nuLnu value to make nice looking plots */
      gsl_vector* temp_x = gsl_vector_alloc( obs->EnDissProcPoint[i]->N );
      gsl_vector* temp_y = gsl_vector_alloc( obs->EnDissProcPoint[i]->N );
      gsl_vector_memcpy( temp_x, obs->EnDissProcPoint[i]->vPoint );
      gsl_vector_memcpy( temp_y, obs->EnDissProcPoint[i]->LvPointAvg );
      gsl_vector_mul( temp_y, temp_x );
      nuLnuMaxValue = gsl_vector_max( temp_y );

      /* set up gnuplot settings */
      gnuplotFile << "reset" << std::endl;
      gnuplotFile << "set terminal postscript enhanced color" << std::endl;
      gnuplotFile << "set output \"" << obs->EnDissProcPoint[i]->whoAmI() << "PointVsEle.eps\"" << std::endl;
      gnuplotFile << "set size square" << std::endl;
      gnuplotFile << "set border lw 2" << std::endl;
      gnuplotFile << "set grid" << std::endl;
      gnuplotFile << "set title \"" << obs->EnDissProcPoint[i]->whoAmI() << " spectra evolution (POINT SOURCE)\" " << std::endl;
      gnuplotFile << "set xlabel \"{/Symbol g}\" " << std::endl;
      gnuplotFile << "set ylabel \"{/Symbol n} L_{/Symbol n}\"" << std::endl;
      gnuplotFile << "set logscale xy" << std::endl;
      gnuplotFile << "set yrange [" << nuLnuMaxValue/yLowerScale << ":" << nuLnuMaxValue*yUpperScale << "]" << std::endl;
      gnuplotFile << "plot\\" << std::endl;
      gnuplotFile << "\"LvPointAvg_" << obs->EnDissProcPoint[i]->whoAmI() << "\" using ($2):($1*$3) w l lw 3 notitle,\\" << std::endl;

      for( int k=1;k<_r->getMaxIndex();k++ ) 
	{
	  _r->update( k );
	  if( _r->ifSaveRadius( ) ) { 
	    gnuplotFile << "\"LvPoint_" << obs->EnDissProcPoint[i]->whoAmI()<< "_" << _r->getPosition( ) << "\" using ($2):($1*$3) " << chooseLineGnuplot( obs->EnDissProcPoint[i]->whoAmI() ) << " w l notitle";
	    if( k != _r->getMaxIndex()-1 ) 
	      gnuplotFile << ",\\" << std::endl; }
	}
      
      gnuplotFile.close();
    }
}

void makePointFlare( observer* obs, scfgp* _cfg, jetGeometry* _r  )
{
  std::ofstream gnuplotFile;
  std::stringstream filename;
  std::string dir = _cfg->get<std::string>("output");

  for( int i=0; i<obs->EnDissProcPoint.size(); i++ )
    {
      filename.str("");
      bazinga::print_info("plt","Plot (point source flare)",obs->EnDissProcPoint[i]->whoAmI());
      filename << dir << "/" << obs->EnDissProcPoint[i]->whoAmI() << "Flare.gp";
      gnuplotFile.open( filename.str().c_str() );
      
      /* set up gnuplot settings */
      gnuplotFile << "reset" << std::endl;
      gnuplotFile << "set terminal postscript enhanced color" << std::endl;
      gnuplotFile << "set output \"" << obs->EnDissProcPoint[i]->whoAmI() << "Flare.eps\"" << std::endl;
      gnuplotFile << "set size square" << std::endl;
      gnuplotFile << "set border lw 2" << std::endl;
      gnuplotFile << "set grid" << std::endl;
      gnuplotFile << "set title \"" << obs->EnDissProcPoint[i]->whoAmI() << " flares\" " << std::endl;
      gnuplotFile << "set xlabel \"r/r_{0}\" " << std::endl;
      gnuplotFile << "set ylabel \"L_{/Symbol n}\"" << std::endl;
      gnuplotFile << "plot\\" << std::endl;
      for( int k=0; k<obs->getFreqListSize(); k++ )
	{
	  gnuplotFile << "\"FlarePoint_" << obs->EnDissProcPoint[i]->whoAmI()<< "_" << obs->getFreqList( k ) << "Hz" << "\" using ($1):($2) " << chooseLineGnuplot( obs->EnDissProcPoint[i]->whoAmI() ) << " w l title\"{/Symbol n} = " << obs->getFreqList( k ) << "Hz \"";
	  if( k != obs->getFreqListSize()-1 ) 
	    gnuplotFile << ",\\" << std::endl;
	}
      gnuplotFile.close();
    }
}

void makeUpeRPlots( observer* obs, scfgp* _cfg, jetGeometry* _r )
{
  std::ofstream gnuplotFile;
  std::stringstream filename;
  std::string dir = _cfg->get<std::string>("output");
  
  bazinga::info("plt","Plot upe vs jetGeometry");
  filename.str("");
  filename << dir << "/" << "UpeR.gp";
  gnuplotFile.open( filename.str().c_str() );
  
  /* set up gnuplot settings */
  gnuplotFile << "reset" << std::endl;
  gnuplotFile << "set terminal postscript enhanced color" << std::endl;
  gnuplotFile << "set output \"UpeR.eps\"" << std::endl;
  gnuplotFile << "set size square" << std::endl;
  gnuplotFile << "set border lw 2" << std::endl;
  gnuplotFile << "set grid" << std::endl;
  gnuplotFile << "set title \"Upe vs jetGeometry\"" << std::endl;
  gnuplotFile << "set xlabel \"jetGeometry\" " << std::endl;
  gnuplotFile << "set ylabel \"Upe\"" << std::endl;
  gnuplotFile << "set logscale xy" << std::endl;
  
  gnuplotFile << "plot\\" << std::endl;
  
  for( int i=0; i<obs->EnDissProc.size(); i++ )
    {
      bazinga::print_info("plt","Add",obs->EnDissProc[i]->whoAmI());
      gnuplotFile << "\"UpeR_" << obs->EnDissProc[i]->whoAmI() << "\" using ($1):($2) " << chooseLineGnuplot( obs->EnDissProc[i]->whoAmI() ) << " w l title \"" << obs->EnDissProc[i]->whoAmI() << "\"";
      if( i != obs->EnDissProc.size()-1 )
	gnuplotFile << ",\\" << std::endl;
    }

  gnuplotFile.close();
}

void makeRmPlots( observer* obs, scfgp* _cfg, jetGeometry* _r ) {
  std::ofstream gnuplotFile;
  std::stringstream filename;
  std::string dir = _cfg->get<std::string>("output");
  
  bazinga::info("plt","Plot Rm vs jetGeometry");
  filename.str("");
  filename << dir << "/" << "Rm.gp";
  gnuplotFile.open( filename.str().c_str() );
  
  /* set up gnuplot settings */
  gnuplotFile << "reset" << std::endl;
  gnuplotFile << "set terminal postscript enhanced color" << std::endl;
  gnuplotFile << "set output \"Rm.eps\"" << std::endl;
  gnuplotFile << "set size square" << std::endl;
  gnuplotFile << "set border lw 2" << std::endl;
  gnuplotFile << "set grid" << std::endl;
  gnuplotFile << "set title \"Rm vs jetGeometry\"" << std::endl;
  gnuplotFile << "set xlabel \"jetGeometry\" " << std::endl;
  gnuplotFile << "set ylabel \"Rm\"" << std::endl;
  gnuplotFile << "set logscale xy" << std::endl;
  
  gnuplotFile << "plot\\" << std::endl;
  
  for( int i=0; i<obs->ExtPl.size(); i++ )
    {
      bazinga::print_info("plt","Add",obs->ExtPl[i]->whoAmI());
      gnuplotFile << "\"Rm_" << obs->ExtPl[i]->whoAmI() << "\" using ($1):($2) " << chooseLineGnuplot( obs->EnDissProcPoint[i]->whoAmI() ) << " w l title \"" << obs->ExtPl[i]->whoAmI() << "\"";
      if( i != obs->ExtPl.size()-1 )
	gnuplotFile << ",\\" << std::endl;
    }

  gnuplotFile.close();
}

void makeExtPldLdRPlots( observer* obs, scfgp* _cfg, jetGeometry* _r ) {
  std::ofstream gnuplotFile;
  std::stringstream filename;
  std::string dir = _cfg->get<std::string>("output");
  
  bazinga::info("plt","Plot dLdR vs R");
  filename.str("");
  filename << dir << "/" << "dLdR.gp";
  gnuplotFile.open( filename.str().c_str() );
  
  /* set up gnuplot settings */
  gnuplotFile << "reset" << std::endl;
  gnuplotFile << "set terminal postscript enhanced color" << std::endl;
  gnuplotFile << "set output \"dLdR.eps\"" << std::endl;
  gnuplotFile << "set size square" << std::endl;
  gnuplotFile << "set border lw 2" << std::endl;
  gnuplotFile << "set title \"dLdR vs R\"" << std::endl;
  gnuplotFile << "set grid" << std::endl;
  gnuplotFile << "set xlabel \"R\" " << std::endl;
  gnuplotFile << "set ylabel \"dLdR\"" << std::endl;
  gnuplotFile << "set logscale xy" << std::endl;
  
  gnuplotFile << "plot\\" << std::endl;
  
  for( int i=0; i<obs->ExtPl.size(); i++ )
    {
      bazinga::print_info("plt","Add",obs->ExtPl[i]->whoAmI());
      gnuplotFile << "\"dLdR_" << obs->ExtPl[i]->whoAmI() << "\" using ($1):($2) " << chooseLineGnuplot( obs->EnDissProcPoint[i]->whoAmI() ) << " w l title \"" << obs->ExtPl[i]->whoAmI() << "\"";
      if( i != obs->ExtPl.size()-1 )
	gnuplotFile << ",\\" << std::endl;
    }

  gnuplotFile.close();
}





