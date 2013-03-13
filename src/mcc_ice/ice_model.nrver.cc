#include <limits>

#include "nr.h"
#include "ice_model.h"

ice_model::ice_model() {
  freq=1.4;
	
  tsky=5.;
  tice=269.;
  tw=271.;
  nice=complex_t(2, 0.15);
  nsea=complex_t(9, 2.5);

  n1=0;
  n2=0;

}

ice_model::~ice_model() {
  if (n1 != 0) {
    delete [] s1;
    delete [] h;
    delete [] hp;
  }
  if (n2 != 0) {
    delete [] s2;
    delete [] z;
    delete [] zp;
  }
}

long ice_model::set_surface(real_t *x1a, real_t *ha, long n1a, 
		real_t *x2a, real_t *za, long n2a) {

//  numeric_limits<float> nl;
	
  if (n1 != 0) {
    delete [] s1;
    delete [] h;
    delete [] hp;
  }
  if (n2 != 0) {
    delete [] s2;
    delete [] z;
    delete [] zp;
  }
  
  n1=n1a;
  n2=n2a;
  
  s1=new real_t[n1];
  h=new real_t[n1];
  hp=new real_t[n1];
  
  for (long i=0; i<n1; i++) {
    //printf("%f %f\n", x1a[i], ha[i]);
    s1[i]=x1a[i];
    h[i]=ha[i];
  }
  spline(s1-1, h-1, n1, 0, 0, hp-1);
  for (long i=0; i<n1; i++) {
    if (isnan(hp[i])) {
      fprintf(stderr, "Cubic spline fitting failed\n");
      return 1;
    }
    //printf("%f\n", hp[i]);
  }

  s2=new real_t[n1];
  z=new real_t[n1];
  zp=new real_t[n1];
  
  for (long i=0; i<n1; i++) {
    //printf("%f %f\n", x2a[i], za[i]);
    s2[i]=x2a[i];
    z[i]=za[i];
  }
  spline(s2-1, z-1, n2, 0, 0, zp-1);
  for (long i=0; i<n1; i++) {
    if (isnan(hp[i])) {
      fprintf(stderr, "Cubic spline fitting failed\n");
      return 1;
    }
    //printf("%f\n", zp[i]);
  }

  return 0;

}

real_t ice_model::intercept(real_t &x, real_t &y, 	//current coords
		real_t &vx, real_t &vy, 		//direction vector
		int &stype,				//substance type
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
  int interface;
  real_t eps;

  eps=1e-7;		//should be set at the object-level

  if (stype == 0) {
    t=spl_intercept(x, y, vx, vy, s1, h, hp, n1, x2, y2);
    if (t < 0) {
      x2=x+vx;
      y2=y+vy;
      vx2=vx;
      vy2=vy;
    } else {
      m=spl_slope(x2, s1, h, hp, n1);
      r=m_refract(vx, vy, m, real_t(1.), nice.real(), vx2, vy2, pol);
      if (ran2(&seed) > r) {
        //ray is refracted--now in ice
        st2=1;
      } else {
        //ray is reflected--still in air
        st2=0;
        m_reflect(vx2, vy2, m, vx2, vy2);
      }
    }
  } else if (stype == 1) {
    t1=spl_intercept(x, y, vx, vy, s1, h, hp, n1, x2a, y2a);
    t2=spl_intercept(x, y, vx, vy, s2, z, zp, n2, x2b, y2b);
    if (t1 < eps && t2 < eps) {
      x2=x+vx;
      y2=y+vy;
      vx2=vx;
      vy2=vy;
      st2=1;
      t=-1;
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
      m=spl_slope(x2, s1, h, hp, n1);
      r=m_refract(vx, vy, m, nice.real(), real_t(1.), vx2, vy2, pol);
      if (ran2(&seed) > r) {
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
      m=spl_slope(x2, s2, z, zp, n2);
      r=(m_refract(complex_t(vx), complex_t(vy), m, 
		      nice, nsea, vx2a, vy2a, pol)).real();
      if (ran2(&seed) > r) {
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

real_t ice_model::do_iter(real_t x, real_t y, 
		real_t vx, real_t vy, 
		int pol,
		long &nb,
		int &stf)
{
  long maxnb;		//maximum number of bounces
  real_t *t;
  int *st;
  long i;
  real_t tb;
  real_t abs_coef;

  maxnb=100;
  
  t=new real_t[maxnb-1];
  st=new int[maxnb];

  abs_coef=2*3.1415926537*nice.imag()*freq/0.3;

  st[0]=0;
  for (i=1; i<maxnb; i++) {
    st[i]=st[i-1];
    t[i-1]=intercept(x, y, vx, vy, st[i], pol);
    printf("%f %f\n", x, y);
    if (t[i-1] < 0 || st[i] == 2) break;
  }
  printf("\n");

  //that's it for the ray-tracing, now for the radiative transfer:
  nb=i;
  if (nb == 1) {
    fprintf(stderr, "Premature ray-trace exit\n");
    return tsky;
  }

  if (st[nb] == 1) {
    fprintf(stderr, "Ray-trace exit in ice\n");
  }
    
  if (st[nb]==0) tb=tsky; 
  	else if (st[nb]==1) tb=tice;
  	else if (st[nb]==2) tb=tw;

  for (long j=nb-1; j>=0; j-=1) {
    if (st[j]==1) tb=tice-(tice-tb)*exp(-abs_coef*t[j]);
  }

  return tb;

}

real_t ice_model::tb(real_t x0, real_t y0,
		real_t zangle, real_t fov,
		int pol,
		long niter,
		real_t &err)
{
  real_t maxh;
  real_t *tb;
  real_t vx, vy;
  real_t x, y;
  real_t t0;
  long nb;
  int sf;
  real_t tbave;

  //bit of a hack:
  maxh=h[0];
  for (long i=1; i<n1; i++) if (h[i] > maxh) maxh=h[i];
  maxh+=0.1;

  tb=new real_t[niter];

  for (long i=0; i<niter; i++) {
    //set direction vector:
    //random perturbation to zenith angle:
    zangle=zangle+2*ran2(&seed)*fov-fov;
    vx=sin(zangle);
    vy=-cos(zangle);

    printf("%g %g %g %g\n", x0, y0, vx, vy);

    //to avoid unnecessary computation, we advance the ray
    //until it is within spitting distance of the ice surface:
    t0=(maxh-y0)/vy;
    y=maxh;
    x=x0+vx*t0;

    tb[i]=do_iter(x, y, vx, vy, pol, nb, sf);
    if (nb == 1) {
      continue;
      i--;
    }
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

}

