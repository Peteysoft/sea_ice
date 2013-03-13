#include <limits>

#include <sys/timeb.h>

#include "gsl/gsl_spline.h"
#include "gsl/gsl_rng.h"

#include "mcc_ice_model.h"

mcc_ice_model::mcc_ice_model() {
  long iseed;
  timeb now;

  freq=1.4;
	
  tsky=5.;
  tice=269.;
  tw=271.;
  //nice=complex_t(2, 0.15);
  nsea=complex_t(9, 2.5);

  haccel=gsl_interp_accel_alloc();
  zaccel=gsl_interp_accel_alloc();

  ice_air=NULL;
  ice_water=NULL;

  //initialize random number generator:
  rng=gsl_rng_alloc(gsl_rng_mt19937);

  //seed random number generator:
  ftime(&now);
  iseed=(long) now.time + ((long) now.millitm) << 7;
  gsl_rng_set(rng, iseed);

  //logfs=fopen("piss_traj.txt", "w");
  logfs=NULL;

  fov=13.6*M_PI/180;

  eps=1e-7;

  maxnb=100;

  abs1=NULL;
  abs2=NULL;
  rindex1=NULL;
  rindex2=NULL;

}

void mcc_ice_model::set_fov(real_t fov) {
  fov=fov*M_PI/180;
}

extern "C" void ckw_(double *, double *, double *, double *, double *);

void mcc_ice_model::set_parm(real_t nu, real_t ti, real_t tw1, real_t sw, real_t ts) {
  complex_t eps;
  real_t eps1, eps2;

  tice=ti;
  tw=tw1;
  tsky=ts;
  freq=nu;

  //determine the complex permittivity of the sea water:
  ckw_(&freq, &tw, &sw, &eps1, &eps2);
  eps=complex_t(eps1, eps2);
  //printf("(%f, %f)\n", eps1, eps2);
  nsea=sqrt(eps);
  //printf("(%f, %f)\n", nsea.real(), nsea.imag());
}

mcc_ice_model::~mcc_ice_model() {
  if (ice_air != NULL) gsl_spline_free(ice_air);
  if (ice_water != NULL) gsl_spline_free(ice_water);

  gsl_interp_accel_free(haccel);
  gsl_interp_accel_free(zaccel);

  if (abs1 != NULL) delete [] abs1;
  if (abs2 != NULL) delete [] abs2;
  if (rindex1 != NULL) delete [] rindex1;
  if (rindex2 != NULL) delete [] rindex2;

  gsl_rng_free(rng);

  if (logfs!=NULL) fclose(logfs);
}

//print out the spline fit:
void mcc_ice_model::print_fit(char *fname, real_t res) {
  FILE *fs;
  long n;

  real_t x;
  real_t y;

  fs=fopen(fname, "w");

  n=(ice_air->x[ice_air->size-1]-ice_air->x[0])/res+1;

  fprintf(fs, "%ld\n", n);
  for (long i=0; i<n; i++) {
    x=ice_air->x[0]+res*i;
    y=gsl_spline_eval(ice_air, x, haccel);
    fprintf(fs, "%20.12lg %lg\n", x, y);
  }

  n=(ice_water->x[ice_air->size-1]-ice_water->x[0])/res+1;

  fprintf(fs, "%d\n", n);
  for (long i=0; i<n; i++) {
    x=ice_water->x[0]+res*i;
    y=gsl_spline_eval(ice_water, x, zaccel);
    fprintf(fs, "%20.12lg %lg\n", x, y);
  }

  fclose(fs);

}

void mcc_ice_model::setlogfile(char *fname) {
  logfs=fopen(fname, "w");
}

int mcc_ice_model::set_surface(double *x1a, double *ha, complex_t *eps1, long n1a, 
		double *x2a, double *za, complex_t *eps2, long n2a) {

  complex_t rind;
  int err=0;

  if (ice_air != NULL) gsl_spline_free(ice_air);
  if (ice_water != NULL) gsl_spline_free(ice_water);
	
  if (abs1 != NULL) delete [] abs1;
  if (abs2 != NULL) delete [] abs2;
  if (rindex1 != NULL) delete [] rindex1;
  if (rindex2 != NULL) delete [] rindex2;

  ice_air=gsl_spline_alloc(gsl_interp_cspline, n1a);
  ice_water=gsl_spline_alloc(gsl_interp_cspline, n2a);

  err=gsl_spline_init(ice_water, x2a, za, n2a);
  if (err != 0) {
    fprintf(stderr, "Cubic spline fitting failed\n");
    return err;
  }
  err=gsl_spline_init(ice_air, x1a, ha, n1a);
  if (err != 0) {
    fprintf(stderr, "Cubic spline fitting failed\n");
    return err;
  }

  abs1=new real_t[n1a];
  rindex1=new real_t[n1a];
  abs2=new real_t[n2a];
  rindex2=new real_t[n2a];

  //find minimum and maximum y values for each surface:
  hmin=ha[0];
  hmax=ha[0];

  rind=sqrt(eps1[0]);
  rindex1[0]=rind.real();
  abs1[0]=4*M_PI*rind.imag()*freq/0.3;
  for (long i=1; i<n1a; i++) {
    if (ha[i] < hmin) hmin=ha[i];
    if (ha[i] > hmax) hmax=ha[i];

    rind=sqrt(eps1[i]);
    rindex1[i]=rind.real();
    abs1[i]=4*M_PI*rind.imag()*freq/0.3;
    //printf("eps'=%lf, eps''=%lf\n", eps1[i].real(), eps1[i].imag());
    //printf("n=%lf, a=%lf\n", rindex1[i], abs1[i]);
  }

  zmin=za[0];
  zmax=za[0];

  rind=sqrt(eps2[0]);
  rindex2[0]=rind.real();
  abs2[0]=4*M_PI*rind.imag()*freq/0.3;
  for (long i=1; i<n2a; i++) {
    if (za[i] < zmin) zmin=za[i];
    if (za[i] > zmax) zmax=za[i];

    rind=sqrt(eps2[i]);
    rindex2[i]=rind.real();
    abs2[i]=4*M_PI*rind.imag()*freq/0.3;
  }

  return err;

}

real_t mcc_ice_model::intercept(real_t &x, real_t &y, 	//current coords
		real_t &vx, real_t &vy, 		//direction vector
		int &stype,				//substance type (0=air, 1=ice, 2=water)
		int &interface,			//for interpolation of aux. data...
		double &int_ind,
		int pol)				//polarisation

{
  real_t x2, y2;	//resultant coords
  real_t vx2, vy2;	//resultant dir. vector
  complex_t vx2a, vy2a;
  int st2;		//new substance
	
  real_t t, t1, t2;
  real_t r;		//reflection coefficient
  real_t m;		//slope
  real_t x2a, y2a, x2b, y2b;
  //int interface;

  real_t rind;		//index of refractivity

  long indd;
  double frac;

  double cy1, cy2;

  if (stype == 0) {
    t=spl_intercept(x, y, vx, vy, ice_air, haccel, hmin, hmax, x2, y2);
    if (t < 0) {
      x2=x+vx;
      y2=y+vy;
      vx2=vx;
      vy2=vy;
      interface=-1;
      st2=0;
    } else {
      int_ind=interpolate(ice_air->x, ice_air->size, x2, -1);
      indd=(long) int_ind;
      frac=int_ind-(double) indd;
      rind=(1-frac)*rindex1[indd]+frac*rindex1[indd+1];

      //printf("n=%f\n", rind);

      m=gsl_spline_eval_deriv(ice_air, x2, haccel);
      cy1=gsl_spline_eval(ice_air, x2, haccel);
      cy2=gsl_spline_eval(ice_air, x2+eps, haccel);
      //printf("m (anal)=%g, m (num)=%g\n", m, (cy2-cy1)/eps);
      r=m_refract(vx, vy, m, real_t(1.), rind, vx2, vy2, pol);
      //printf("%f\n", r);
      if (gsl_rng_uniform(rng) > r) {
        //ray is refracted--now in ice
        st2=1;
      } else {
        //ray is reflected--still in air
        st2=0;
        m_reflect(vx, vy, m, vx2, vy2);
      }
      //interpolate refractive index and absorption coeffs:
      interface=0;
    }
  } else if (stype == 1) {
    t1=spl_intercept(x, y, vx, vy, ice_air, haccel, hmin, hmax, x2a, y2a);
    t2=spl_intercept(x, y, vx, vy, ice_water, zaccel, zmin, zmax, x2b, y2b);
    if (t1 < eps && t2 < eps) {
      x2=x+vx;
      y2=y+vy;
      vx2=vx;
      vy2=vy;
      st2=1;
      t=-1;
      interface=-1;
    } else if (t1 < eps) {
      t=t2;
      interface=1;
    } else if (t2 < eps) {
      t=t1;
      interface=0;
    } else if (t1 < t2) {
      t=t1;
      interface=0;
    } else if (t2 < t1) {
      t=t2;
      interface=1;
    }

    if (interface == 0) {
      x2=x2a;
      y2=y2a;

      int_ind=interpolate(ice_air->x, ice_air->size, x2, -1);
      indd=(long) int_ind;
      frac=int_ind-(double) indd;
      rind=(1-frac)*rindex1[indd]+frac*rindex1[indd+1];

      //printf("n=%f\n", rind);

      m=gsl_spline_eval_deriv(ice_air, x2, haccel);
      r=m_refract(vx, vy, m, rind, real_t(1.), vx2, vy2, pol);
      if (gsl_rng_uniform(rng) > r) {
        //ray is refracted--now in air
        st2=0;
      } else {
        //ray is reflected--still in ice
        st2=1;
	m_reflect(vx, vy, m, vx2, vy2);
      }
    } else {
      x2=x2b;
      y2=y2b;

      int_ind=interpolate(ice_water->x, ice_water->size, x2, -1);
      indd=(long) int_ind;
      frac=int_ind-(double) indd;
      rind=(1-frac)*rindex2[indd]+frac*rindex2[indd+1];

      //printf("n=%f\n", rind);

      m=gsl_spline_eval_deriv(ice_water, x2, zaccel);
      r=abs(m_refract(complex_t(vx, 0), complex_t(vy, 0), m, 
		      complex_t(rind, 0), nsea, vx2a, vy2a, pol));
      //printf("(%f, %f)\n", nsea.real(), nsea.imag());
      //r=m_refract(vx, vy, m, rind, abs(nsea), vx2, vy2, pol);
      //printf("%f\n", r);
      if (gsl_rng_uniform(rng) > r) {
        //ray is refracted--now in water
	vx2=vx2a.real();
	vy2=vy2a.real();
	st2=2;
      } else {
        //ray is reflected--still in ice
	st2=1;
	m_reflect(vx, vy, m, vx2, vy2);
      }
    }
  }

  x=x2;
  y=y2;
  vx=vx2;
  vy=vy2;
  stype=st2;

  return t;

}

real_t mcc_ice_model::do_iter(real_t x, real_t y, 
		real_t vx, real_t vy, 
		int pol,
		long &nb,
		int &stf)
{
  long i;
  real_t tb;
  real_t t[maxnb-1];
  int st[maxnb];
  real_t absarr[maxnb];

  //for interpolation of auxiliary data:
  double int_ind, frac;
  long ind;
  int interface;

  //for de-bugging:
  real_t xall[maxnb];
  real_t yall[maxnb];

  if (logfs!=NULL) fprintf(logfs, "%lg %lg\n", x, y);

  st[0]=0;
  //absarr[0]=0;

//#define debug 1
  
#ifdef debug
  xall[0]=x;
  yall[0]=y;
#endif

  for (i=1; i<maxnb; i++) {
    //printf("Intercept #%d\n", i);
    st[i]=st[i-1];
    t[i-1]=intercept(x, y, vx, vy, st[i], interface, int_ind, pol);

#ifdef debug
    xall[i]=x;
    yall[i]=y;
#endif

    //printf("t=%f\n", t[i-1]);
    //interpolate the attenuation coefficient:
    ind=(long) int_ind;
    frac=int_ind-(double) ind;
    if (interface == 0) {
      absarr[i-1]=(1-frac)*abs1[ind]+frac*abs1[ind+1];
    } else if (interface == 1) {
      absarr[i-1]=(1-frac)*abs2[ind]+frac*abs2[ind+1];
    }
    if (logfs!=NULL) fprintf(logfs, "%lg %lg\n", x, y);
    //printf("%d ", st[i]);
    if (t[i-1] < 0 || st[i] == 2) break;
  }
  //printf("\n");
  if (logfs!=NULL) {
    fprintf(logfs, "%d\n", st[i]);	  
    if (st[i] == 2) {
      fprintf(logfs, "* %lg %lg\n", x+vx, y+vy);
    }
    fprintf(logfs, "\n");
  }

  //that's it for the ray-tracing, now for the radiative transfer:
  nb=i;
  stf=st[nb];
  if (nb == 1 && t[0] < 0) {
    fprintf(stderr, "Premature ray-trace exit\n");
    //printf("Premature ray-trace exit\n\n");
    return tsky;
  }

  if (st[nb] == 1) {
    fprintf(stderr, "Ray-trace exit in ice\n");
    //printf("Ray-trace exit in ice\n");
  }
  //printf("\n");
    
  if (st[nb]==0) tb=tsky; 
  	else if (st[nb]==1) tb=tice;
  	else if (st[nb]==2) tb=tw;
  //printf("tb=%f\n", tb);
  //printf("%d\n", nb);
  //for (long j=0; j<=nb; j++) printf("%d ", st[j]);
  //printf("\n");

#ifdef debug
  printf("       t         x         y       alpha      s\n");
  for (long k=0; k<nb; k++) {
    printf("%10.2f %10.2f %10.2f %10.2f %5d\n", t[k], xall[k], yall[k], absarr[k], st[k]);
  }
  printf("%10.2f %10.2f %10.2f %10.2f %5d\n", 0., xall[nb], yall[nb], 0, st[nb]);
#endif

  for (long j=nb-1; j>0; j-=1) {
    //if (absarr[j] < 0 || absarr[j-1] < 0) {
    if (st[j]==1) {
      //printf("alpha=%f, dt=%f\n", (absarr[j-1]+absarr[j])/2, t[j]);
      tb=tice-(tice-tb)*exp(-t[j]*(absarr[j-1]+absarr[j])/2);
      //printf("tb=%f\n", tb);
    }
  }

  return tb;

}

#include "inverf.h"

real_t mcc_ice_model::tb(real_t x0, real_t y0,
		real_t zangle, real_t offset,
		int pol,
		long niter,
		real_t &err)
{
  real_t *tb;
  real_t vx, vy;
  real_t x, y;
  real_t t0;
  long nb;
  int sf;
  real_t tbave;
  real_t zangle1;
  real_t minrnum;
  real_t dza;		//perturbation to zenith angle
  real_t rnum;
  real_t tmp;
  int sgn;

  tb=new real_t[niter];

  minrnum=erf(asin(offset)/sqrt(2)/fov);

  for (long i=0; i<niter; i++) {
    //set direction vector:
    //random perturbation to zenith angle:
    rnum=minrnum+gsl_rng_uniform(rng)*(1-minrnum);
    tmp=sin(inverf(rnum)*sqrt(2)*fov);
    sgn=((int) (gsl_rng_uniform(rng)*2))*2-1;
    dza=sgn*asin(sqrt(tmp*tmp-offset*offset));

    //printf("dtheta=%f\n", dza*180/M_PI);

    //printf("zangle=%f, ", zangle);
    //zangle1=zangle+2*gsl_rng_uniform(rng)*fov-fov;
    zangle1=zangle+dza;
    //printf("%f\n", zangle1);
    vx=sin(zangle1);
    vy=-cos(zangle1);

    //printf("x0: (%g, %g), v=(%g, %g)\n", x0, y0, vx, vy);

    //to avoid unnecessary computation, we advance the ray
    //until it is within spitting distance of the ice surface:
    //t0=(maxh-y0)/vy;
    //y=maxh;
    //x=x0+vx*t0;

    tb[i]=do_iter(x0, y0, vx, vy, pol, nb, sf);
    //printf("%f\n", tb[i]);
    //if (nb == 1) {
      //i--;
      //continue;
    //}
  }

  //calculate statistics :
  tbave=0;
  for (long i=0; i<niter; i++) tbave+=tb[i];
  tbave/=niter;
  err=0;
  for (long i=0; i<niter; i++) {
    real_t diff=tb[i]-tbave;
    err+=diff*diff;
  }

  err/=(niter-1);
  err=sqrt(err);

  delete [] tb;

  return tbave;

}

