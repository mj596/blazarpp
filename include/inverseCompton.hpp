#ifndef _inverseCompton_hpp
#define _inverseCompton_hpp 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_dilog.h>

#include <cmath>

namespace inverseCompton {
  
  /* unnormalized FKN in total electron cooling rate due to inverse
     compton scattering */
  double fKN( double b, int KN );
  double FKN( double b, double alpha, int KN  );
  
  /* fiso and fiso*/
  double f( double g, double e0, double e, double miu );
  double fiso( double g, double e0, double e );
  
}

#endif /* _inverseCompton_hpp */
