#include <math.h>
//#include "nr.h"

typedef double Doub;

Doub inverfc(Doub p);

inline Doub inverf(Doub p) {return inverfc(1.-p);}

