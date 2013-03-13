;fit the Fresnel equations to a given pair of emissivity profiles:
function cost_fit_fresnel, n, params=params
  zangle=*params[0]
  ev=*params[1]
  eh=*params[2]
  
  theta=asin(sin(zangle)/n)
  eh_fr=1-((cos(zangle)-n*cos(theta))/(cos(zangle)+n*cos(theta)))^2
  ev_fr=1-((n*cos(zangle)-cos(theta))/(n*cos(zangle)+cos(theta)))^2

  return, total((ev_fr-ev)^2+(eh_fr-eh)^2)

end

;fit linearly scaled Fresnel equations to a given pair of emissivity profiles:
function cost_fit_fresnel1, n, params=params
  zangle=*params[0]
  ev=*params[1]
  eh=*params[2]
  
  theta=asin(sin(zangle)/n)
  eh_fr=1-((cos(zangle)-n*cos(theta))/(cos(zangle)+n*cos(theta)))^2
  ev_fr=1-((n*cos(zangle)-cos(theta))/(n*cos(zangle)+cos(theta)))^2

  a=(regress([ev_fr, eh_fr], [ev, eh], const=b))[0]

  return, total((a*ev_fr+b-ev)^2+(a*eh_fr+b-eh)^2)
end

;fit the Fresnel equations to a given pair of emissivity profiles:
function cost_fit_fresnel2, n1
  common frefit, zangle, ev, eh

;  zangle=*params[0]
;  ev=*params[1]
;  eh=*params[2]

  n=complex(n1[0], n1[1])
  
  theta=asin(sin(zangle)/n)
  eh_fr=1-((cos(zangle)-n*cos(theta))/(cos(zangle)+n*cos(theta)))^2
  ev_fr=1-((n*cos(zangle)-cos(theta))/(n*cos(zangle)+cos(theta)))^2

  return, total((ev_fr-ev)^2+(eh_fr-eh)^2)

end

function fit_fresnel, zangle, ev, eh, minn=minn, maxn=maxn, tol=tol

;+*****************************************************************************
;               FIT_FRESNEL
;******************************************************************************
;
; usage:        n=fit_fresnel(zangle, ev, eh [,minn=minn] [,maxn=maxn] [,tol=tol])
;
; purpose:      Given vertically and horizontally polarized emissivities as a
;		function of zenith angle, determines the real refractive index
;		producing the best fit to the Fresnel equations.
;
; parameters:
;	zangle	Array of zenith angles.
;
;	ev	Emissivity in vertical polarization.
;
;	eh	Emissivity in horizontal polarization.
;
;	tol	Tolerance of fit.
;
;	minn	Lowest permissible value of refractive index with which 
;		to bracket the minimum.
;
;	maxn	Highest permissible value of refractive index.
;
; output:	A real refractive index.
;
; author:	Peter Mills (pmills@iup.physik.uni-bremen.de)
;
; history:	2009-6-02, documented
;
;-*************************************************************************************

  if n_elements(minn) ne 1 then minn=1.
  if n_elements(maxn) ne 1 then maxn=10.
  if n_elements(tol) ne 1 then tol=0.001

  params=ptrarr(3)
  params[0]=ptr_new(zangle)
  params[1]=ptr_new(minn)
  params[2]=ptr_new(maxn)

  n=golden("cost_fit_fresnel", minn, maxn, tol, params=params)

  return, n

end

function fit_fresnel1, zangle, ev, eh, a, b, minn=minn, maxn=maxn, tol=tol

;+*****************************************************************************
;               FIT_FRESNEL1
;******************************************************************************
;
; usage:        n=fit_fresnel1(zangle, ev, eh, a, b [,minn=minn] [,maxn=maxn] [,tol=tol])
;
; purpose:      Given vertically and horizontally polarized emissivities as a
;		function of zenith angle, determines the real refractive index
;		plus linear scaling coefficients producing the best fit to 
;		the Fresnel equations.
;
; parameters:
;	zangle	Array of zenith angles.
;
;	ev	Emissivity in vertical polarization.
;
;	eh	Emissivity in horizontal polarization.
;
;	tol	Tolerance of fit.
;
;	minn	Lowest permissible value of refractive index with which 
;		to bracket the minimum.
;
;	maxn	Highest permissible value of refractive index.
;
; output:	Returns a real refractive index, plus:
;
;	a	Multiplication coefficient.
;
;	b	Constant offset.
;
; author:	Peter Mills (pmills@iup.physik.uni-bremen.de)
;
; history:	2009-6-02, documented
;
; *notes:
; This routine will be useful if we have several radiance measurements from SMOS or
; some other instrument from the same ground location, but several different viewing
; angles.  By first fitting the re-scaled Fresnel equations, we can extrapolate the
; brightness temperature to a constant and arbitrary viewing angle (e.g., 50 degrees)
; thus the retrieval algorithm need only be designed for single viewing angle.
;
;-*************************************************************************************
  if n_elements(minn) ne 1 then minn=1.
  if n_elements(maxn) ne 1 then maxn=10.
  if n_elements(tol) ne 1 then tol=0.001

  params=ptrarr(3)
  params[0]=ptr_new(zangle)
  params[1]=ptr_new(ev)
  params[2]=ptr_new(eh)

  n=golden("cost_fit_fresnel1", minn, (minn+maxn)/2, maxn, tol, params=params)

  theta=asin(sin(zangle)/n)
  eh_fr=1-((cos(zangle)-n*cos(theta))/(cos(zangle)+n*cos(theta)))^2
  ev_fr=1-((n*cos(zangle)-cos(theta))/(n*cos(zangle)+cos(theta)))^2

  a=(regress([ev_fr, eh_fr], [ev, eh], const=b))[0]

  return, n

end

