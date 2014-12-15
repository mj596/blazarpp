 /**
    @file energyDissProc.hpp
    @author Mateusz Janiak
*/

#ifndef _energyDissProc_hpp
#define _energyDissProc_hpp 1

#include "baseClass.hpp"

class electrons;

/** @class energyDissProc
    Interface class for every energy dissipation process; provides memory management 
    and other common methods shared by every process */

class energyDissProc : public baseClass {
  int saveLum, saveUpeVsR;

public:
  /** class needs to have access to electrons to calculate luminosities etc */
  electrons* ele;
  
  /** data matrices and vectors;
      first index is radius */

  /** photon field energy - primed */
  gsl_vector* ep;
  /** photon field monochromatic energy density - primed */
  gsl_matrix* upe;
  /** frequency - primed */
  gsl_vector* vp; 
  /** monochromatic apparent luminosity - primed */
  gsl_vector* Lpv;
  /** freqency - point source */
  gsl_vector* vPoint;
  /** monochromatic apparent luminosity for point source */
  gsl_vector* LvPoint;
  /** monochromatic apparent luminosity for point source averaged over radius */
  gsl_vector* LvPointAvg;

  /** data vector for storing energy density vs radius */
  gsl_vector* upe_r;

  /** used to specify radial maximum of injected electrons - used only for non-uniform injections */
  double injRm;

  /** set this to 1 if luminosity is to be calculated with constant u_ext;
      in such case u' will be set to a constant value given by u'(injrm) */
  int luminosityConstU;

  /** set this to 1 if luminosity is to be calculated with constant v';
      in such case v' will be set to a constant value given by v'(injrm) */
  int luminosityConstNu;

  /** data vector to store information on electrons energies that cool in the Klein-Nishina regime
      for each radius a data on electrons in KN in being saved */
  gsl_vector* gammaKN;

  /** boundary energies */
  double epMin, epMax;
  /** boundary frequencies */
  double vpMin, vpMax;
  /** ep integration step */
  double dLogE;

  /** set this to true after setting upe_r */
  bool flag_upe_r;
  
  /** constructor
      @param scfgp
      @param jetGeometry
      @param electrons
      @param id */
  energyDissProc( scfgp* cfg, jetGeometry* r, electrons* ele, std::string _id );
  
  /** destructor */
  ~energyDissProc( );

  /** set intrinsic luminosities */
  virtual void setLpv( ) { };

  /** calculate intrinsic luminosity
      @param v - frequency (jet co-moving frame)
      @returns L'_v */
  virtual double calculateLpv( double v ) { };

  /** update all process internal parameters to current radius */
  virtual void update( ) { };

  /** get d gamma \ d t for particular process
      @param g - electron Lorentz factor
      @returns d gamma \ d t */
  virtual double dotg( double g ) { };
  virtual void printInfo( ) { };

  /** save calulated luminosity */
  void saveLuminosity( );
  
  /** save energy density vs radius info */
  void saveUpeR( );

  /* methods for allocating memeory */
  void allocateLpv( );
  void allocateLvPoint( );
  void allocateLvPointAvg( );
  void allocateUpeR( );
  void allocateUpe( );
  void allocateGammaKN( );

  /* freeing memory methods */
  void freeLpv( );
  void freeLvPoint( );
  void freeLvPointAvg( );
  void freeUpe( );
  void freeUpeR( );
  void freeGammaKN( );

  /* common getters */
  double get_ep( int i ) { return gsl_vector_get( ep, i ); }
  double get_upe( int i ) { return gsl_matrix_get( upe, r->getIndex( ), i ); }
  double get_upe_r( ) { return gsl_vector_get( upe_r, r->getIndex( ) ); }
  double get_vPoint( int i ) { return gsl_vector_get( vPoint, i ); }
  double get_LvPoint( int i ) {  return gsl_vector_get( LvPoint, i ); }
  double get_LvPointAvg( int i ) {  return gsl_vector_get( LvPointAvg, i ); }
  double get_vp( int i ) { return gsl_vector_get( vp, i ); }
  double get_Lpv( int i ) { return gsl_vector_get( Lpv, i ); }

  double get_gammaKN( ) { return gsl_vector_get( gammaKN, r->getIndex( ) ); }  

  /* common setters */
  void set_ep( int i, double val ) { gsl_vector_set( ep, i, val ); }
  void set_upe( int i, double val ) { gsl_matrix_set( upe, r->getIndex( ), i, val ); }
  void set_upe_r( double val );
  void set_vPoint( int i, double val ) { gsl_vector_set( vPoint, i, val ); }
  void set_LvPoint( int i, double val ) { gsl_vector_set( LvPoint, i, val ); }
  void set_LvPointAvg( int i, double val ) { gsl_vector_set( LvPointAvg, i, val ); }
  void set_vp( int i, double val ) { gsl_vector_set( vp, i, val ); }
  void set_Lpv( int i, double val ) { gsl_vector_set( Lpv, i, val ); }

  double set_gammaKN( double val ) { gsl_vector_set( gammaKN, r->getIndex( ), val ); }  

  /** set information about Klein-Nishina regime
      @param - electron Lorentz factor */
  void set_KN_info( double g );

  /** print info whether particular process for particular electron energy is withing Klein-Nishina region already */
  void print_KN_info( );
  
  /** get log-energy step
      @returns log-energy step */
  double getdLogE( ) { return dLogE; }
};

#endif


