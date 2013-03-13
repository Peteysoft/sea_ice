#include <math.h>

#include "mcc_ice_model.h"

int mcc_ice_idl_wrapper (int argc, void **argv) {

  mcc_ice_model *ice;

  //parameters:
  real_t ti;		//ice temp.
  real_t tw;		//water temp.
  real_t sw;		//sea salinity
  real_t nu;		//frequency
  real_t ts;		//sky temp.

  real_t *x0, *y0;	//initial conditions
  real_t *zangle;
  real_t *offset;

  long ninit;		//number of initial conditions
  long niter;		//number of iterations

  real_t *s, *h, *z;	//ice surface
  long np;

  real_t fov;

  real_t *tbv;
  real_t *tbh;

  real_t *verr;
  real_t *herr;

  real_t *epsp, *epspp;
  complex_t *eps;

  char * logfile=NULL;

  int err=0;

  //argument list:
  //0  x grid
  //1  ice surface
  //2  ice-water interface
  //3  real permittivity
  //4  imaginary permittivity
  //5  number of samples
  //6  initial x coords
  //7  initial y coords
  //8  viewing angle
  //9  offset
  //10 number of initial conditions
  //11 number of iterations
  //12 frequency
  //13 instrument fov
  //14 ice temperature
  //15 water temp.
  //16 water salinity
  //17 sky temperature
  //18 return tbv
  //19 return tbh
  //20 return e_v
  //21 return e_h

  //make sure the number of arguments is right:
  if (argc != 22) return -1;

  //start setting all the values from the argument list:
  //ice topography:
  s=(double *) argv[0];
  h=(double *) argv[1];
  z=(double *) argv[2];
  epsp=(double *) argv[3];
  epspp=(double *) argv[4];
  np=* (long*) argv[5];

  //initial conditions:
  x0=(double *) argv[6];
  y0=(double *) argv[7];
  zangle=(double *) argv[8];
  offset=(double *) argv[9];
  ninit=* (long *) argv[10];

  //parameters:
  niter=* (long *) argv[11];
  nu=* (double *) argv[12];
  fov=* (double *) argv[13];
  ti=* (double *) argv[14];
  tw=* (double *) argv[15];
  sw=* (double *) argv[16];
  ts=* (double *) argv[17];

  //return values:
  tbv=(double *) argv[18];
  tbh=(double *) argv[19];
  verr=(double *) argv[20];
  herr=(double *) argv[21];

  eps=new complex_t[np];

  for (long i=0; i<np; i++) {
    eps[i]=complex_t(epsp[i], epspp[i]);
  }

  ice=new mcc_ice_model();

  if (logfile!=NULL) ice->setlogfile(logfile);
  ice->set_fov(fov);

  ice->set_parm(nu, ti, tw, sw, ts);

  err=ice->set_surface(s, h, eps, np, s, z, eps, np);
  if (err != 0) return err;

  //ice->print_fit("pissfit.txt", 0.01);

  for (long i=0; i<ninit; i++) {
    tbv[i]=ice->tb(x0[i], y0[i], zangle[i], offset[i], 0, niter, verr[i]);
    tbh[i]=ice->tb(x0[i], y0[i], zangle[i], offset[i], 1, niter, herr[i]);

    printf("%f +/- %f\n", tbv[i], verr[i]);
    printf("%f +/- %f\n", tbh[i], herr[i]);

  }

  delete ice;

  delete [] eps;

  return err;

}
    

