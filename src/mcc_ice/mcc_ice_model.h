#include <stdio.h>

#include "gsl/gsl_spline.h"
#include "gsl/gsl_rng.h"

#include "mcc_ice_util.h"

class mcc_ice_model {
  protected:
    real_t freq;

    real_t tsky;		//sky temperature
    real_t tice;		//ice temperature
    real_t tw;			//water temperature
    //complex_t nice;		//refractive index of ice
    complex_t nsea;		//refractive index of water

    //spline describing ice surface:
    gsl_spline *ice_air;
    gsl_interp_accel *haccel;

    //spline describing ice-water interface:
    gsl_spline *ice_water;
    gsl_interp_accel *zaccel;

    //absorption coeffs. and refractive indices:
    real_t *abs1, *abs2;
    real_t *rindex1, *rindex2;

    //minimum and maximum y-values for each surface:
    real_t hmin, hmax;
    real_t zmin, zmax;

    gsl_rng *rng;		//seed for random number generator

    FILE *logfs;

    real_t fov;

    real_t eps;

    long maxnb;			//maximum number of bounces

  public:
    mcc_ice_model();
    ~mcc_ice_model();

    int set_surface(double *x1a, double *ha, complex_t *eps1, long n1a, 
		    double *x2a, double *za, complex_t *eps2, long n2a);

    void setlogfile(char *fname);
    void set_fov(real_t f1);
    void set_parm(real_t nu, real_t ti, real_t tw1, real_t ts, real_t tsky=5.);

    //one intersection:
    real_t intercept(real_t &x0, real_t &y0, real_t &vx, real_t &vy, int &stype, 
		    int &interface, double &int_ind, int pol);
    
    //one iteration:
    real_t do_iter(real_t x0, real_t y0, real_t vx, real_t vy, int pol, long &nb, int &stype);

    //multiple iterations:
    real_t tb(real_t x0, real_t y0, real_t zangle, real_t fov, int pol, long niter, real_t &err);

    void print_fit(char *fname, real_t res);

};
