@ice_model_util

;directly calculates the brightness temperatures for a three-layer
;(water, ice and air) plane-parallel model:

pro three_layer, nu, zangle, h, tice, eps, tbv, tbh, $
		tw=tw, eps_oc=epsw, sal=sal, tsky=tsky, $
		rparm=rparm, coherent=coherent, $
		next_gen=next_gen

  if n_elements(tw) ne 1 then tw=271.		;water temp.
  if n_elements(tsky) ne 1 then tsky=3.		;sky temp.
  if n_elements(sal) ne 1 then sal=35.		;water salinity

  if n_elements(rparm) ne 1 then rparm=h*0.1

  if n_elements(epsw) ne 1 then epsw=eps_oc2(nu, tw, sal)

  nice=sqrt(eps)
  nw=sqrt(epsw)
  c=0.3

  theta1=asin(sin(zangle)/nice)
  theta2=asin(nice*sin(theta1)/nw)

  rh1=(cos(zangle)-nice*cos(theta1))/(cos(zangle)+nice*cos(theta1))
  rv1=(nice*cos(zangle)-cos(theta1))/(nice*cos(zangle)+cos(theta1))

  rh2=(nice*cos(theta1)-nw*cos(theta2))/(nice*cos(theta1)+nw*cos(theta2))
  rv2=(nw*cos(theta1)-nice*cos(theta2))/(nw*cos(theta1)+nice*cos(theta2))

  if keyword_set(coherent) then begin
    abs_coef=4*!pi*nice*nu/c
    ;transmission coefficient:
    sigma=exp(-complex(0, 1)*abs_coef*h/cos(theta1))
    ;stop
    rh1=rh1*rh1
    rv1=rv1*rv1
    rh2=rh2*rh2
    rv2=rv2*rv2
  endif else begin
    abs_coef=4*!pi*imaginary(nice)*nu/c
    ;transmission coefficient:
    sigma=exp(-abs_coef*h/cos(float(theta1)))
    rh1=abs(rh1*rh1)
    rv1=abs(rv1*rv1)
    rh2=abs(rh2*rh2)
    rv2=abs(rv2*rv2)
  endelse

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


