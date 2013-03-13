#include <assert.h>

#include "peteys_tmpl_lib.h"
#include "nrcomplex.h"
#include "nr.h"

#include "ice_model_util.h"

inline real_t sabs(complex_t &x) {return abs(x);}
inline real_t sabs(real_t &x) {return fabs(x);}

//find the slope of a cubic spline:
float spl_slope(real_t x0, 		//interpolation point
		real_t *x, 		//abscissa
		real_t *y, 		//ordinates
		real_t *yp,		//derivatives
		long n) 		//number of points
{
  double xind;
  long lind;
  double frac1, frac2;
  double h;
  real_t result;

  xind=interpolate(x, n, x0);

  lind=long(xind);
  frac1=xind-lind;
  frac2=1-frac1;
  h=x[lind+1]-x[lind];
  result=h*((frac1*frac1-1./3)*yp[lind+1]-(frac2*frac2-1./3)*yp[lind])/2 + 
		(y[lind+1]-y[lind])/h;
  return result;
}

//find the nearest intercept of a line with a cubic spline:
real_t spl_intercept(real_t x0_1, 
		real_t y0_1, 
		real_t vx, 
		real_t vy, 
		real_t *x, 
		real_t *y, 
		real_t *yp, 
		long n,
		real_t &xint,
		real_t &yint)
{
  long i;
  real_t miny;
  real_t maxy;
  long ind;
  int dir;
  long iend;
  real_t eps;
  real_t eps2;
  real_t h;
  real_t t0;
  real_t x0, y0;
  real_t y1, y2;
  real_t yp1, yp2;
  real_t h2, x02, x03;
  real_t a, b, c, d;
  fcomplex coeff[4];
  real_t desc;
  fcomplex root[3];
  fcomplex sqrtdesc;
  fcomplex tmp;
  real_t t;
  int ngood;		//number of good roots
  real_t vx2, vx3;
  real_t xt;

  miny=y[0];
  maxy=y[0];
  for (long i=1; i<n; i++) {
    if (y[i] < miny) miny=y[i];
    if (y[i] > maxy) maxy=y[i];
  }

  ind=bin_search(x, n, x0_1);

  dir=vx/fabs(vx);
  if (dir < 0) iend=0; else iend=n-2;

  eps=1e-6;
  eps2=1e-10;

  for (i=ind; i<iend;i+=dir) {
    h=x[i+1]-x[i];
    if (x0_1 < x[i] || x0_1 > x[i+1]) {
      if (dir < 0) {
        t0=(x[i+1]-x0_1)/vx;
        x0=h;
        y0=y0_1+vy*t0;
      } else {
        t0=(x[i]-x0_1)/vx;
        x0=0.;
        y0=y0_1+vy*t0;
      }
    } else {
      x0=x0_1-x[i];
      t0=0.;
      y0=y0_1;
    }
    //check to see if intersections are even possible:
    if ((y0 > maxy && vy > 0) || (y0 < miny && vy < 0)) {
      t=-1;
      break;
    }

    //generate polynomial coefficients:
    y1=y[i];
    y2=y[i+1];
    yp1=yp[i];
    yp2=yp[i+1];

    x02=x0*x0;
    x03=x02*x0;
    h2=h*h;
    vx2=vx*vx;
    vx3=vx2*vx;

    d=((x03-h2*x0)*yp2+(-x03+3*h*x02-2*h2*x0)*yp1+6*x0*y2+(6*h-6*x0)*y1-6*h*y0)/h/6;
    c=((3*vx*x02-h2*vx)*yp2+(-3*vx*x02+6*h*vx*x0-2*h2*vx)*yp1+6*vx*y2-6*vx*y1-6*h*vy)/h/6;
    b=(3*vx2*x0*yp2+(3*h*vx2-3*vx2*x0)*yp1)/h/6;
    a=(vx3*yp2-vx3*yp1)/h/6;

    root[1].r=-1; root[1].i=-1; root[2].r=-1; root[2].i=-1;
    coeff[0].r=a; coeff[1].r=b; coeff[2].r=c; coeff[3].r=d;
    coeff[0].i=0; coeff[1].i=0; coeff[2].i=0; coeff[3].i=d;

    if (fabs(a) < eps2) {
      if (fabs(b) == 0) {
        root[0].r=-d/c;
      } else {
        desc=c*c-4*b*d;
	sqrtdesc.r=desc;
	sqrtdesc=Csqrt(sqrtdesc);
	tmp.r=-coeff[2].r;
	tmp.i=-coeff[2].i;
        root[0]=Cdiv(Cadd(tmp, sqrtdesc), coeff[1]);
	root[0].r=root[0].r/2.;
	root[0].i=root[0].i/2;
	sqrtdesc.r=-sqrtdesc.r;
	sqrtdesc.i=-sqrtdesc.i;
        root[1]=Cdiv(Cadd(tmp, sqrtdesc), coeff[1]);
	root[1].r=root[1].r/2.;
	root[1].i=root[1].i/2;
      }
    } else {
      zroots(coeff, 4, root, 1);
      //t=solve_cubic(a, b, c, d)
    }

    if (t0 == 0.) {
      t=100000;		//yet another hack...
      ngood=0;
      for (long j=0; j<2; j++) {
        if (fabs(root[j].i) < eps && root[j].r > eps) {
	  ngood++;
	  if (root[j].r < t) t=root[j].r;
        }
      }
    } else {
      t=100000;
      ngood=0;
      for (long j=0; j<2; j++) {
        if (fabs(root[j].i) < eps && root[j].r > 0) {
	  ngood++;
	  if (root[j].r < t) t=root[j].r;
        }
      }
    }

    if (ngood == 0) {
      t=-1;		//continue may knock us out of the loop...
      continue;
    }

    xt=vx*t+x0;
    if ((xt > h) || (xt < 0)) {
      t=-1;
      continue;
    }

    break;
  }

  if (t > 0) {
    xint=x[i]+xt;
    splint(x, y, yp, n, xint, &yint);
    t=t+t0;
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

