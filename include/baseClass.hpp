/**
    @file baseClass.hpp
    @author Mateusz Janiak
*/

#ifndef _baseClass_hpp
#define _baseClass_hpp 1

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <scfgp.hpp>
#include <bazinga.hpp>
#include "jetGeometry.hpp"

/**
   @class baseClass
   Base class for all other classes used in blazar++. Provides unified way of sharing scfgp, id and jetGeometry
*/

class baseClass {
public:
  /** id characterizing object */
  std::string id;

  /** every class needs to have an information about parameters available - this is realized by attaching a pointer to scfgp class */
  scfgp* cfg;

  /** every class needs to have an information of current radius and overal geomtery - this is realized by attaching a pointer to radius class */
  jetGeometry* r;

  /** vectors size - since this is often used will be copied from model at the beginning */
  int N;

  baseClass( );
  baseClass( scfgp* cfg, jetGeometry* r, std::string id );
  ~baseClass( );

  /** tell me what is my id */
  std::string whoAmI( );

  /** print basic information about myself (virtual) */
  virtual void printInfo( ) { };

  /** get beta from Lorentz factor \Gamma; it shouldn't probably be here */
  double beta( double x ); 
  
  /** get N */
  int getN( ) { return N; }
};

#endif /* _baseClass_hpp */

