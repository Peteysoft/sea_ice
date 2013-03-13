
pro min2d, fname, tol, x, y, a, b, rrat, ntrial, $
		params=params, plot=plt, maxiter=maxiter

;+*****************************************************************************
;		MIN2D
;******************************************************************************
;
; min2d, fname, tol, x, y, a, b, rrat, ntrial [,params=params] [,maxiter=maxiter]
;
; purpose:
;	minimize a function in two dimensions;
;	just picks a bunch of random points, takes the lowest,
;	picks a bunch more random points in a smaller diameter
;	and so on...
;
; parameters:
;	fname	Name of function to minimize.  Should take two inputs and 
;		a set of parameters (param=param).
;
;	x	Starting value of first (x) variable.
;
;	y	Starting value of second (y) variable.
;
;	a	Search radius in x direction.
;
;	b	Search radius in y direction.
;
;	rrat	Reduces search radius by this much each iteration.
;
;	ntrial	Number of trials each iteration.
;
; keywords:
;	params	Set of parameters to pass to cost function.
;
;	maxiter	Maximum number of iterations (default is 100).
;
; author:	Peter Mills (pmills@iup.physik.uni-bremen.de)
;
; history:	2009-6-02 documented.
;
;-*****************************************************************************

  if n_elements(maxiter) ne 1 then maxiter=100

  if keyword_set(plt) then begin
    r=indgen(maxiter)*255/(maxiter-1)
    colind=rgbtoindex(r, fltarr(maxiter), 255-r)
  endif
  
  z=fltarr(ntrial)

  if keyword_set(plt) then begin
    plot, findgen(10), xrange=[x-a/2, x+a/2], yrange=[y-b/2, y+b/2], /nodata
  endif 

  iter=0L
  
  repeat begin
    theta=randomu(seed, ntrial)*2*!pi
    r=sqrt(randomu(seed, ntrial))

    x=a*r*cos(theta)+x
    y=b*r*sin(theta)+y

    ;for i=0, ntrial-1 do begin
    ;  z[i]=call_function(fname, x[i], y[i], params=params)
    ;endfor
    z=call_function(fname, x, y, params=params)

    zmin=min(z, j)

    if keyword_set(plt) then oplot, x, y, color=colind[iter], psym=3

    x=x[j]
    y=y[j]

    a=a/rrat
    b=b/rrat

    iter=iter+1

    ;stop

  endrep until zmin lt tol or iter ge maxiter

  tol=zmin

end
    
