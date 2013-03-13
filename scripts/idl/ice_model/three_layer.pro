@ice_model_util

;directly calculates the brightness temperatures for a three-layer
;(water, ice and air) plane-parallel model:

pro three_layer, nu, zangle, h, tice, eps, tbv, tbh, $
		tw=tw, eps_oc=epsw, sal=sal, tsky=tsky, $
		rparm=rparm, $
		next_gen=next_gen

;+*****************************************************************************
;		THREE_LAYER
;******************************************************************************
;
; usage:	three_layer, nu, zangle, h, tice, eps, tbv, tbh
;
; purpose:	Performs a radiative transfer simulation for a sheet of sea ice
;		having uniform properties throughout.
;
; input:	nu	frequence in GHz
;
;		zangle 	viewing angle in radians
;
;		h	ice thickness [m]
;
;		tice	ice temperature [K]
;
;		eps	complex relative permittivity of sea ice
;
; output:	tbv	emitted brightness temperature, vertical polarization [K]
;
;		tbh	emitted brightness temperature, horizontal polarization [K]
;
; optional:	tw	water temperature (default is 271 K)
;
;		eps_oc	for directly specifying the permittivity of sea water
;
;		sal	salinity of the sea water (default is 35 psu)
;
;		tsky	sky brightness temperature 
;
;		rparm	roughness parameter (used in combination with /next_gen)
;
;		/next_gen	combines radiative transfer solution with model from
;			Menashi (JGR 20 (C12): 22569 (1993))
;
;		dependencies:	ice_model_util.pro
;
; author: Peter Mills (pmills@iup.physik.uni-bremen.de)
;
; history: 2009-6-02, formally documented
;          2010-8-14 PM got rid of /coherent keyword
;
; *notes: see "invert3layer"
;______________________________________________________________________________

  if n_elements(tw) ne 1 then tw=271.		;water temp.
  if n_elements(tsky) ne 1 then tsky=3.		;sky temp.
  if n_elements(sal) ne 1 then sal=35.		;water salinity

  if n_elements(rparm) ne 1 then rparm=h*0.1

  if n_elements(epsw) ne 1 then epsw=eps_oc2(nu, tw, sal)

  nice=sqrt(eps)
  nw=sqrt(epsw)
  c=0.3

  ;this is consistent with MEMLS:
  nice1=sqrt(float(eps))
  nw1=sqrt(float(epsw))

  theta1=asin(sin(zangle)/nice1)
  theta2=asin(nice1*sin(theta1)/nw1)

  rh1=(cos(zangle)-nice1*cos(theta1))/(cos(zangle)+nice1*cos(theta1))
  rv1=(nice1*cos(zangle)-cos(theta1))/(nice1*cos(zangle)+cos(theta1))

  rh2=(nice1*cos(theta1)-nw1*cos(theta2))/(nice1*cos(theta1)+nw1*cos(theta2))
  rv2=(nw1*cos(theta1)-nice1*cos(theta2))/(nw1*cos(theta1)+nice1*cos(theta2))

  ;print, "Rwi:", rh2, rv2

    abs_coef=4*!pi*imaginary(nice)*nu/c
    ;transmission coefficient:
    sigma=exp(-abs_coef*h/cos(abs(theta1)))
    rh1=rh1*rh1
    rv1=rv1*rv1
    rh2=rh2*rh2
    rv2=rv2*rv2

  ;solution to the matrix equation:
  sigma2=sigma*sigma
  th2=((rh2-1)*sigma*tw+(rh1-1)*rh2*sigma2*tsky+ $
		(rh2*sigma2+(1-rh2)*sigma-1)*tice)/(rh1*rh2*sigma2-1)

  tv2=((rv2-1)*sigma*tw+(rv1-1)*rh2*sigma2*tsky+ $
		(rv2*sigma2+(1-rv2)*sigma-1)*tice)/(rv1*rv2*sigma2-1)

  if keyword_set(next_gen) then begin
    ;Menashi "dielectric slab" model:
    beta=2*!pi*nice*nu/c
    int1=sqrt(1-sin(zangle)^2/eps)
    l=h/int1
    sigma_l=rparm/int1
    int2=exp(-beta*sigma_l)
    tbh=(1-rh1)*th2*(1-sigma*sqrt(rh1*rh2)*int2)/(1+sigma*sqrt(rh1*rh2)*int2)
    tbv=(1-rv1)*tv2*(1-sigma*sqrt(rv1*rv2)*int2)/(1+sigma*sqrt(rv1*rv2)*int2)
  endif else begin
    tbh=(1-rh1)*th2
    tbv=(1-rv1)*tv2
  endelse

end


