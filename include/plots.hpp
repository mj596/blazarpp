#ifndef _plots_hpp
#define _plots_hpp 1

#include <scfgp.hpp>
#include "jetGeometry.hpp"
#include "energyDissProc.hpp"
#include "observer.hpp"
#include "electrons.hpp"
#include "QuasarAccDisk.hpp"

#ifdef USE_ROOT
#include "root.hpp"
#endif

const static double xLowerScale = 5;
const static double xUpperScale = 5;
const static double yLowerScale = 1.0e3;
const static double yUpperScale = 1.3;

void makeNgammaPlots( electrons* obj, scfgp* _model, jetGeometry* _r );
void makeLPointPlots( observer* obj, scfgp* _model, jetGeometry* _r );
void makeLPointVsElePlots( observer* obj, scfgp* _model, jetGeometry* _r );
void makeLPointSumPlots( observer* obj, scfgp* _model, jetGeometry* _r );
void makePointFlare( observer* obs, scfgp* _model, jetGeometry* r );
void makeUpeRPlots( observer* obj, scfgp* _model, jetGeometry* _r );
void makeRmPlots( observer* obs, scfgp* _model, jetGeometry* _r );
void makeExtPldLdRPlots( observer* obj, scfgp* _model, jetGeometry* _r );
void makeLQAccdPlots( QuasarAccDisk* obj, scfgp* model );

#ifdef USE_ROOT
/* ROOT plots */
void makeNgammaPlotROOT( gsl_vector* x, gsl_vector* y, TCanvas* c1, TMultiGraph* mg, std::string type, jetGeometry* r );
void makeLPlotROOT( gsl_vector* x, gsl_vector* y, TCanvas* c1, TMultiGraph* mg, std::string type, jetGeometry* r );
void makeLSumPlotROOT( gsl_vector* x, gsl_vector* y, TCanvas* c1, TMultiGraph* mg, std::string type, jetGeometry* r );
Color_t chooseColorRoot( std::string type );
int checkJetGeometryToRootPlot( jetGeometry *r );
#endif

/* auxillary */
std::string chooseLineGnuplot( std::string type );

#endif /* _plots_H */
