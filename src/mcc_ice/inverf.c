#include <math.h>
//#include "nr.h"

#include "inverf.h"

Doub inverfc(Doub p) {
  Doub x, err, t, pp;

  if (p >= 2.) return -100.;
  if (p <= 0.0) return 100.;

  pp=(p < 1.0)? p:2.-p;
  t=sqrt(-2.*log(pp/2));

  x= -0.70711*((2.30753+t*0.27061)/(1+t*(0.99229+t*0.04481))-t);

  for (int j=0;j<2;j++) {
    err=erfc(x)-pp;
    x += err/(1.12837916709551257*exp(-x*x)-x*err);
  }
  return (p<1.0 ? x:-x);
}


