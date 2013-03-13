#include <assert.h>

#include <stdio.h>

#include "gsl/gsl_spline.h"
#include "gsl/gsl_poly.h"

#include "mcc_ice_util.h"

inline real_t sabs(complex_t &x) {return abs(x);}
inline real_t sabs(real_t &x) {return fabs(x);}

//find the nearest intercept of a line with a cubic spline:
real_t spl_intercept(real_t x0_1, 		//initial x coord
		real_t y0_1, 			//initial y coord
		real_t vx, 			//unit direction vector, x coord
		real_t vy, 			//unit direction vector, y coord
		gsl_spline *spline,		//defines the surface of intersection
		gsl_interp_accel *accel,	//why the fuck does this thing have to be separate??
		real_t miny,
		real_t maxy,
		real_t &xint,
		real_t &yint)
{
  long i;
  //real_t miny, maxy;
  long ind1, ind2;
  int dir;
  long iend;
  real_t eps;
  real_t eps2;
  real_t h;
  real_t t0, t1, dt;
  real_t x0, y0;
  real_t y1;
  real_t yp1;
  real_t ypp1, ypp2;
  real_t dypp;
  real_t h2, x02;
  real_t a, b, c, d;
  double desc;
  double root[3];
  int nroots;

  double sqrtdesc;

  real_t t;
  real_t vx2, vx3;
  real_t xt;

  double *x;
  double *y;
  size_t n;

  int fflag;

  //to reduce typing:
  x=spline->x;
  y=spline->y;
  n=spline->size;

  //first we want to find the range of intervals to check:
  if (vy > 0) {
    if (y0_1 < miny) {
      x0=x0_1+vx*(miny-y0_1)/vy;
    } else {
      x0=x0_1;
    }
    x02=x0_1+vx*(maxy-y0_1)/vy;
  } else {
    if (y0_1 > maxy) {
      x0=x0_1+vx*(maxy-y0_1)/vy;
    } else {
      x0=x0_1;
    }
    x02=x0_1+vx*(miny-y0_1)/vy;
  }
  ind1=bin_search(x, n, x0);
  ind2=bin_search(x, n, x02);

  if (ind1 < 0 || ind2 >= n) return -1;

  dir=vx/fabs(vx);

  eps=1e-6;
  eps2=1e-10;
  fflag=0;

  for (i=ind1; dir*(ind2-i) >= 0; i+=dir) {
    h=x[i+1]-x[i];
    dt=h/vx;

    //do this generally, so we're always solving from x[i] which we set at 0:
    t0=(x[i]-x0_1)/vx;
    //this is y for the ray on the left-hand side of the interval:
    y0=y0_1+vy*t0;
    //x_0 is not needed because it's always zero...

    //generate polynomial coefficients:
    //expand about x[i]:
    y1=y[i];
    yp1=gsl_spline_eval_deriv(spline, x[i], accel);
    ypp1=gsl_spline_eval_deriv2(spline, x[i], accel);
    ypp2=gsl_spline_eval_deriv2(spline, x[i+1], accel);

    h2=h*h;
    vx2=vx*vx;
    vx3=vx2*vx;
    dypp=ypp2-ypp1;

    a=dypp*vx3/h/6;
    b=vx2*ypp1/2;
    c=(vx*yp1-vy);
    d=y1-y0;

    if (fabs(a) < eps2) {
      if (fabs(b) == 0) {
        root[0]=-d/c;
	nroots=1;
      } else {
        //calculate the discriminant:
        desc=c*c-4*b*d;
	if (desc < 0) {
	  nroots=0;
	} else {
	  sqrtdesc=sqrt(desc);
	  root[0]=(-c-sqrtdesc)/2/b;
	  root[1]=(-c+sqrtdesc)/2/b;
	  nroots=2;
	}
      }
    } else {
      nroots=gsl_poly_solve_cubic(b/a, c/a, d/a, root, root+1, root+2);
    }

    for (int j=0; j<nroots; j++) {
      if (root[j]/dt <= 1 && root[j]/dt >= 0) {
        t1=root[j];
        t=t1+t0;
	if (t > eps) {
          fflag=1;
          break;
	}
      }
    }

    if (fflag) break;

  }

  if (fflag) {
    real_t yintc, xintc;
    xint=x[i]+t1*vx;
    yint=gsl_spline_eval(spline, xint, accel);
    yintc=y0+t*vy;
  } else {
    t=-1;
  }

  return t;

}

//calculate the direction vector after refraction:
template <class numeric>
numeric m_refract(numeric vx, 
		numeric vy, 
		real_t m, 
		numeric n1, 
		numeric n2, 
		numeric &vxnew, 
		numeric &vynew, 
		int pol)
{
  numeric r;		//reflection coefficient
	
  real_t num;
  real_t tx;
  real_t ty;
  numeric proj;
  numeric tanlen;
  numeric normlen1, normlen;
  numeric cost0, cost1;

  //tangent vector:
  num=sqrt(1+m*m);
  tx=1./num;
  ty=m/num;
  
  //project the slope onto the tangent vector:
  proj=vx*tx+vy*ty;

  //length of component along tangent:
  tanlen=proj*n1/n2;

  //if this is greater than one, we have total internal reflection:
  if (sabs(tanlen) > 1) {
    m_reflect(vx, vy, m, vxnew, vynew);
    return 1;
  }

  //length of component along normal:
  normlen1=-vx*ty+vy*tx;
  normlen=-sqrt(numeric(1.)-tanlen*tanlen)*normlen1/sabs(normlen1);

  vxnew=tanlen*tx+normlen*ty;
  vynew=-normlen*tx+tanlen*ty;

  //if finite(vxnew) ne 1 then stop

  r=sqrt(vxnew*vxnew+vynew*vynew);

  vxnew=vxnew/r;
  vynew=vynew/r;

  cost0=sabs(normlen1);
  cost1=sabs(normlen);

  if (pol == 1) {
    r=(n1*cost0-n2*cost1)/(n1*cost0+n2*cost1);
  } else {
    r=(n1*cost1-n2*cost0)/(n1*cost1+n2*cost0);
  }

  return r*r;

}

//calculate the direction vector after reflection:
template <class numeric>
void m_reflect(numeric vx, numeric vy, 
		real_t m, 
		numeric &vx2, numeric &vy2)
{
  numeric proj1;
  numeric proj2;
  real_t tx, ty;
  real_t num;

  //tangent vector:
  num=sqrt(1+m*m);
  tx=1./num;
  ty=m/num;

  proj1=vx*tx+vy*ty;
  proj2=-vx*ty+vy*tx;

  vx2=proj1*tx+proj2*ty;
  vy2=proj1*ty-proj2*tx;

}

template complex_t m_refract(complex_t vx, complex_t vy, 
		real_t m,
		complex_t n1, complex_t n2, 
		complex_t &vxnew, complex_t &vynew,
		int pol);

template real_t m_refract(real_t vx, real_t vy, 
		real_t m,
		real_t n1, real_t n2, 
		real_t &vxnew, real_t &vynew,
		int pol);

template void m_reflect(real_t vx, real_t vy, real_t m, real_t &vx2, real_t &vy2);

template void m_reflect(complex_t vx, complex_t vy, real_t m, complex_t &vx2, complex_t &vy2);

