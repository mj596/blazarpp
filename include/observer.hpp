/**
    @file observer.hpp
    @author Mateusz Janiak
*/

#ifndef _observer_hpp
#define _observer_hpp 1

#include "baseClass.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "QuasarAccDisk.hpp"

class energyDissProc;

/** @class observer
    Class is mainly used to calculate observed values from intrinsic jet properties such as observed luminosities instead of luminosities in jet co-moving frame or flare */

class observer : public baseClass {
  /* requested parameters */
  double thetaObs, thetaJ;
  double Gamma;
  int saveLumPoint, saveLumPointAvg, saveExtPlVsRm;
  double nu1, nu2, nu3, nu4, nu5, nu6;
  std::string lumModel;

  /* min and max nu hard coded to avoid searching for max and min nu in data; to be done later */
  const static double vmin = 1.0e05;
  const static double vmax = 1.0e30;
  
  /* interpolation */
  const static double nuLnuMax = 1.0e06;
  
  /* method for allocating and freeing space */
  void allocateLvPointSum( );
  void allocateLvPointAvgSum( );
  
  void freeLvPointSum( );
  void freeLvPointAvgSum( );

  /* not used! methods for caculating extended jet luminosity (to be done later) */
  double findr( double r, double theta );
  int locate(gsl_vector* v, double val);
  
 public:
  /** constructor 
      @param scfgp
      @param jetgeometry
      @param id */
  observer( scfgp* _cfg, jetGeometry* _r, std::string _id );
  
  /** destructor */
  ~observer( );

  /** vectors to hold processes being used in calculations */
  std::vector<energyDissProc*> EnDissProc;
  std::vector<energyDissProc*> EnDissProcPoint;
  std::vector<energyDissProc*> ExtPl;
  std::vector<energyDissProc*> ExtSp;
  std::vector<energyDissProc*> ExtGu;

  /** vectors to store summed luminosities */
  gsl_vector* vPointSum;
  gsl_vector* LvPointSum;
  gsl_vector* LvPointAvgSum;

  /** frequency list */
  std::vector<double> freqList;

  void printInfo( );

  /** check if you should calculate flares
      @returns 1 if yes, 0 if no **/
  int ifCalculateFlares( );
  
  /** check how many frequencies to look at wen calculating flares 
      @returns number of frequencies */
  int getFreqListSize( );

  /** get frequency
      @param index
      @returns i-th frequency */
  double getFreqList( int i );
  
  /** add energy dissipation process to list of active ones
      @param energyDissProc */
  void addEnDissProc( energyDissProc* _obj );

  /** show active energy dissipation processes */
  void listEnDissProc( );

  /** show how many active energy dissipation processes there are 
      @returns int */
  int sizeEnDissProc( );

  /** add energy dissipation process to list of those for which we have to calculate 'point source' luminosity
      @param energyDissProc */
  void addEnDissProcPoint( energyDissProc* _obj );

  /** show active energy dissipation processes for which we have to calculate 'point source' luminosity */
  void listEnDissProcPoint( );

  /** show how many active energy dissipation processes for which we have to calculate 'point source' luminosity there are 
      @returns int */
  int sizeEnDissProcPoint( );

  /** add ExtPl process to list of active ones 
      @param energyDissProc */
  void addExtPl( energyDissProc* _obj );

  /** show active ExtPl energy dissipation processes */
  void listExtPl( );

  /** show how many active ExtPl energy dissipation processes there are 
      @returns int */
  int sizeExtPl( );

  /** add ExtSp process to list of active ones 
      @param energyDissProc */
  void addExtSp( energyDissProc* _obj );

  /** show active ExtSp energy dissipation processes */
  void listExtSp( );

  /** show how many active ExtSp energy dissipation processes there are 
      @returns int */
  int sizeExtSp( );

  /** add ExtGu process to list of active ones 
      @param energyDissProc */
  void addExtGu( energyDissProc* _obj );

  /** show active ExtGu energy dissipation processes */
  void listExtGu( );

  /** show how many active ExtGu energy dissipation processes there are 
      @returns int */
  int sizeExtGu( );

  /** caluclate 'point-source' luminosity
      @param energyDissProc */
  void calculatePointLuminosity( energyDissProc* _obj );
  //  void calculateJetLuminosity( energyDissProc* _obj );

  /** caculate a sum of all 'point-source' luminosities from every process involved */
  void sumPointLuminosity( );

  /** calculate averaged 'point-source' luminosity; avergaed over all s-cells in the jet avtive region
      @param energyDissProc */
  void avgPointLuminosity( energyDissProc* _obj );

  /** caculate a sum of all averaged 'point-source' luminosities from every process involved */
  void sumPointAvgLuminosity( );

  /** add Quasar Template to calculated luminosities
      @param QuasarAccDisk */
  void addQAccdLuminosity( QuasarAccDisk* QAccd );

  /** calculate flares
      @param energyDissProc
      @param frequency \nu */
  double calculateFlare( energyDissProc* _obj, double nu );

  /** calculate all flares */
  void calculateAllFlares( );

  /** technical: interpolate data y( x )
     @param gsl_vector x
     @param gsl_vector y
     @param x_value to calculate y( x_value ) */
  double interpolate( gsl_vector* x, gsl_vector* y, double xToInterpolate );

  /** save 'point-source' luminosity
      @param energyDissProc */
  void savePointLuminosity( energyDissProc* x );

  /** save summed 'point-source' luminosity */
  void savePointLuminositySum( );

  /** save summed 'point-source' luminosity along with electron gamma values (obsolete)
      @param gsl_vector gamma */
  void savePointLuminositySum( gsl_vector* gamma );

  /** save averaged 'point-source' luminosity
      @param energyDissProc */
  void saveAveragedPointLuminosity( energyDissProc* x );

  /** save summed averaged 'point-source' luminosity */ 
  void saveAveragedPointLuminositySum( );

  /** save summed averaged 'point-source' luminosity along with electron gamma values (obsolete)
      @param gsl_vector gamma */
  void saveAveragedPointLuminositySum( gsl_vector* gamma );
};

#endif /* _observer_hpp */
