#include "ice_model_util.h"

class ice_model {
  protected:
    real_t freq;

    real_t tsky;		//sky temperature
    real_t tice;		//ice temperature
    real_t tw;			//water temperature
    complex_t nice;		//refractive index of ice
    complex_t nsea;		//refractive index of water

    real_t * s1;		//abcissa for ice surface
    real_t * h;			//coordinates ice surface
    long n1;

    real_t * s2;		//abcissa for ice-water interface
    real_t * z;			//coordinates ice-water interface
    long n2;

    real_t * hp;		//derivatives for cubic spline
    real_t * zp;

    long seed;			//seed for random number generator

  public:
    ice_model();
    ~ice_model();

    long set_surface(real_t *x1a, real_t *ha, long n1a, 
		    real_t *x2a, real_t *za, long n2a);

    //one interseption:
    real_t intercept(real_t &x0, real_t &y0, real_t &vx, real_t &vy, int &stype, int pol);
    
    //one iteration:
    real_t do_iter(real_t x0, real_t y0, real_t vx, real_t vy, int pol, long &nb, int &stype);

    //multiple iterations:
    real_t tb(real_t x0, real_t y0, real_t zangle, real_t fov, int pol, long niter, real_t &err);

};
