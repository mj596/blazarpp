#include "inverseCompton.hpp"

namespace inverseCompton {
  
#define B1  (63.0/40.0)
#define B2  (441.0/200.0)
#define B3  (101.0/35.0)
#define B4  (705.0/196.0)
#define B5  (969.0/224.0)
#define B6  (3647.0/720.0)
#define B7  (1598.0/275.0)
#define B8  (3969.0/605.0)
#define B9  (8365.0/1144.0)
#define B10 (76329.0/9464.0)
#define B11 (28089.0/3185.0)
#define B12 (13403.0/1400.0)
#define B13 (112371.0/10880.0)
#define B14 (102495.0/9248.0)
#define B15 (34412.0/2907.0)
  
  double fKN( double b, int KN ) {
    double fg;
    if( KN ) {
      if (b<0.1) { fg = 1.0 - b*(B1 + b*(B2 - b*(B3 + b*(B4 - B5*b)))); }
      else { fg = ((0.5*b+6.0+6.0/b)*log(1.0+b) - ((11.0/12.0)*b*b*b + 6.0*b*b + 9.0*b + 4)/((1.0+b)*(1.0+b)) - 2.0 + 2.0e0*gsl_sf_dilog(-b));
	fg *= 9.0/(b*b*b); }
    }
    else { fg = 1.0; }
    return fg; }
  
#undef B15
#undef B14
#undef B13
#undef B12
#undef B11
#undef B10
#undef B9
#undef B8
#undef B7
#undef B6
#undef B5
#undef B4
#undef B3
#undef B2
#undef B1

  double FKN( double b, double alpha, int KN ) {
    if( KN ) { return pow( 1.0+b, alpha-1.0 ); }
    else { return 1.0; }
  }

  double f( double g, double e0, double e, double miu ) {
    double w,wp,b,t,fx;
    
    w = e/g;
    b = 2.0*(1.0-miu)*e0*g;
  
    if (e>e0 && e<b*g/(1+b)) {
      wp = 1.0-w;
      fx = (1.0+0.5*w*w/wp-2.0*w/(b*wp)+2.0*w*w/(b*b*wp*wp));
      // fx = 1.0;
      if (fx>0.0) { return fx; }
      else { return 0.0; }
    }
    else { return 0.0; }
  }

  double fiso( double g, double e0, double e ) {
    double w,b,t,fx;
    if ( e<g ) {
      w = e/g;
      b = 4.0*e0*g;
      t = w/((1.0-w)*b);
      if((t>(1.0/(4.0*g*g))) && (t<1.0)) {
	fx = 2.0*t*log(t) + t + 1.0 - 2.0*t*t + 0.5*b*t*b*t/(1.0+b*t)*(1.0-t);
	return fx; }
    }
    return 0.0; }
}
