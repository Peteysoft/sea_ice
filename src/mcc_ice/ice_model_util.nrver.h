#include <complex.h>

//#define complex_t complex
typedef complex<float> complex_t;
typedef double real_t;

//find the slope of a cubic spline:
float spl_slope(real_t x0, 		//interpolation point
		real_t *x, 		//abscissa
		real_t *y, 		//ordinates
		real_t *yp,		//derivatives
		long n); 		//number of points

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
		real_t &yint);

//calculate the direction vector after refraction:
template<class numeric>
numeric m_refract(numeric vx, 
		numeric vy, 
		real_t m, 
		numeric n1, 
		numeric n2, 
		numeric &vxnew, 
		numeric &vynew, 
		int pol);

//calculate the direction vector after reflection:
template<class numeric>
void m_reflect(numeric vx, numeric vy, 
		real_t m, 
		numeric &vx2, numeric &vy2);


