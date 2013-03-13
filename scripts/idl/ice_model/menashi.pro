@ice_model_util

;directly calculates the brightness temperatures for a three-layer
;(water, ice and air) plane-parallel model:

pro menashi, nu, zangle, h, tice, eps, tbv, tbh, rparm=rparm, $
		tw=tw, eps_oc=epsw, sal=sal

  if n_elements(tw) ne 1 then tw=271.		;water temp.
  if n_elements(sal) ne 1 then sal=35.		;water salinity

  if n_elements(rparm) ne 1 then rparm=h*0.1

  if n_elements(epsw) ne 1 then epsw=eps_oc2(nu, tw, sal)

  nice=sqrt(eps)
  nw=sqrt(epsw)
  c=0.3
  abs_coef=2*!pi*imaginary(nice)*nu/c

  theta1=asin(sin(zangle)/nice)
  theta2=asin(nice*sin(theta1)/nw)

  rh1=(cos(zangle)-nice*cos(theta1))/(cos(zangle)+nice*cos(theta1))
  rv1=(nice*cos(zangle)-cos(theta1))/(nice*cos(zangle)+cos(theta1))

  rh1=abs(rh1*rh1)
  rv1=abs(rv1*rv1)

  rh2=(nice*cos(theta1)-nw*cos(theta2))/(nice*cos(theta1)+nw*cos(theta2))
  rv2=(nw*cos(theta1)-nice*cos(theta2))/(nw*cos(theta1)+nice*cos(theta2))

  rh2=abs(rh2*rh2)
  rv2=abs(rv2*rv2)

  ;Menashi "dielectric slab" model:
  beta=2*!pi*nice*nu/c
  int1=sqrt(1-sin(zangle)^2/eps)
  l=h/int1
  sigma_l=rparm/int1
  int2=exp(-beta*sigma_l)
  a=exp(-4*abs_coef*l)

  emh=(1-rh1)*(1-a*rh2)*(1-sqrt(a*rh1*rh2)*int2) / $
		(1 - a*rh1*rh2)/(1+sqrt(a*rh1*rh2)*int2)
  emv=(1-rv1)*(1-a*rv2)*(1-sqrt(a*rv1*rv2)*int2) / $
		(1 - a*rv1*rv2)/(1+sqrt(a*rv1*rv2)*int2)

  tbh=emh*tice
  tbv=emv*tice

end


