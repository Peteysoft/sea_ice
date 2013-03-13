@three_layer

function cost3layer, eps
  common c3l, tbv, tbh, nu, zangle, h, tice, tw, sal, fname, rparm, next_gen

  eps1=complex(exp(eps[0]), exp(eps[1]))
  call_procedure, fname, nu, zangle, h, tice, eps1, tbv1, tbh1, $
		tw=tw, sal=sal, rparm=rparm, next_gen=next_gen
  ;stop

  return, (tbv-tbv1)^2+(tbh-tbh1)^2
end

function cost3layer2, x, y, params=params
  common c3l, tbv, tbh, nu, zangle, h, tice, tw, sal, fname, rparm, next_gen

  eps1=complex(exp(x), exp(y))
  call_procedure, fname, nu, zangle, h, tice, eps1, tbv1, tbh1, $
		tw=tw, sal=sal, rparm=rparm, next_gen=next_gen

  return, (tbv-tbv1)^2+(tbh-tbh1)^2
end

function invert3layer, nu1, zangle1, h1, tice1, tbv1, tbh1, $
		tw=tw1, sal=sal1, tol=tol, rparm=rparm1, $
		menashi=menashi, next_gen=next_gen1, $
		stochastic=stochastic, eps0=eps0, eps_range=eps_range, $
		iter=iter, reduction_factor=facter

;+*****************************************************************************
;		INVERT3LAYER
;******************************************************************************
;
; usage:	eps=invert3layer(nu, zangle, h, tice, tbv, tbh)
;
; purpose:	Inverts the three-layer, plane-parallel, radiative-transfer
;		ice emissivity model to retrieve effective permittivity from
;		brightness temperature, ice thickness, ice temperature and
;		zenith angle.
;
; input:	nu	Frequency in GHz
;
;		zangle	Zenith angle in radians.
;
;		h	Ice thickness in m.
;
;		tice	Ice temperature in K.
;
;		tbv	Vertically polarized brightness temperature in K.
;
;		tbh	Horizontally polarized brightness temperature in K.
;
; output:	Effective relative complex permittivity.
;
; keywords:	tol	Tolerance of inversion scheme.  Returns error on exit.
;
;		/stochastic	Use stochastic inversion algorithm.
;			Default is downhill simplex.
;
;		iter	Number of trials per step in stochastic algorithm.
;
;		reduction_factor	Factor by which to reduce the search
;			radius between each step in stochastic algorithm.
;
;		eps0	Starting complex permittivity.
;
;		eps_range	Range of permittivities (search space).
;
; dependencies: three_layer
;
; author:	Peter Mills (pmills@iup.physik.uni-bremen.de)
;
; history:	2009-6-02, documented.
;
; *notes:
; Suppose we have an inversion algorithm that is designed for a single viewing
; angle and uses constant tie-points.  We can extrapolate the algorithm to
; other viewing angles by first deriving the effective permittivities corresponding
; to the brightness temperatures of the tie point and then run the model forward
; again for an array of viewing angles using the "three_layer" routine.
;
;-*****************************************************************************

  common c3l, tbv, tbh, nu, zangle, h, tice, tw, sal, fname, rparm, next_gen
  tbv=tbv1
  tbh=tbh1
  nu=nu1
  zangle=zangle1
  h=h1
  tice=tice1
  if n_elements(rparm1) eq 1 then rparm=rparm1
  if keyword_set(next_gen1) then next_gen=1
  if n_elements(sal1) ne 0 then sal=sal1 else sal=35.
  if n_elements(tw1) ne 0 then tw=tw1 else tw=271.

  if n_elements(tol) eq 0 then tol=1e-4

  eps_im_min=0.0001
  eps_re_min=2.


  p0=[alog(eps_re_min), alog(eps_im_min)]
  
  if keyword_set(menashi) then begin 
    fname="menashi"
  endif else begin
    fname="three_layer"
  endelse

  if keyword_set(stochastic) then begin
    if n_elements(iter) ne 1 then iter=200
    if n_elements(factor) ne 1 then factor=1.5
    if n_elements(eps0) ne 2 then p0=[1.84, -4.6] else p0=alog(eps0)
    if n_elements(eps_range) ne 2 then scale=[3., 9.21] else scale=alog(eps_range)
    x=p0[0]
    y=p0[1]
    ;stop
    min2d, "cost3layer2", tol, x, y, scale[0], scale[1], factor, iter
    eps=complex(exp(x), exp(y))
  endif else begin
    if n_elements(eps0) ne 2 then p0=alog([2, 0.0001]) else p0=alog(eps0)
    if n_elements(eps_range) ne 2 then begin
      eps_ow=eps_oc2(nu, tw, sal)
      scale=[alog(float(eps_ow))-p0[0], alog(imaginary(eps_ow))-p0[1]]
    endif else begin
      scale=alog(eps_range)
    endelse
    eps=amoeba(tol, function_name="cost3layer", p0=p0, scale=scale)
    if eps[0] eq -1 then return, eps
    tol=sqrt(cost3layer(eps))
    eps=complex(exp(eps[0]), exp(eps[1]))
  endelse

  return, eps

end

pro test_invert3lay, h, tice

  nu=1.4
  zangle=50*!dtor

  ntest=100
  eps1=complexarr(ntest)
  eps2=complexarr(ntest)

  for i=0L, ntest-1 do begin
    eps1[i]=complex(randomu(seed, 1)*9+1, randomu(seed, 1))
    three_layer, nu, zangle, h, tice, eps1[i], tbv, tbh

    eps2[i]=invert3layer(nu, zangle, h, tice, tbv, tbh)

  endfor

  stop

end


