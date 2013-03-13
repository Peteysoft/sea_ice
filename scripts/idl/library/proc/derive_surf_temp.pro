@net_heat_flux

;simple thermodynamic model for determining surface temperature
;as a function of thickness:

;computes Qnet*d/k + Tw - Ts
;where:
;Qnet is surface heat flux
;d is ice thickness
;k is thermal conductivity
;Tw is water temperature
;Ts is surface temperature
function tserr, Ts, params=params
  d=params[0]		;ice thickness
  k=params[1]		;thermal conductivity
  v=params[2]		;wind speed
  Ta=params[3]		;air temperature
  rh=params[4]		;relative humidity
  rho=params[5]		;air density
  qsw=params[6]		;short wave heat flux
  tw=params[7]		;water temperature
  cc=params[8]		;cloud-cover
  ccc=params[9]		;cloud-cover coefficient

  ;looks a bit weird, now that the code doesn't really do anything...
  q=net_heat_flux(ts, ta, v, rh, cc, tw=tw, rho=rho)+qsw

  params[10]=q

  ;print, "Net flux:", q

  ;stop

  return, q*d/k+tw-ts

end

;solves for surface temperature based on thermodynamic equilibrium
;parameters:
;d   = ice thickness
;ta  = air temperature
;v   = wind speed
;rh  = relative humidity
;cc  = cloud cover
;sw  = shortwave flux

function derive_surf_temp, d, ta, v, rh, cc, sw, qstar=q, $
			tw=tw, rho=rho, alpha=alpha

  params=fltarr(11)

  if n_elements(tw) ne 1 then tw=271.
  if n_elements(rho) ne 1 then rho=1.2

  ;cloud-cover coefficient:
  ccc=0.8

  ;solar "constant:"
  ;s=1366.

  ;calculate short-wave heat flux:
  ;sw=s*sw_flux(lat, t)*(1-0.62*cc)*(1-alpha)

  ;print, "shortwave flux:", sw

  params[1]=2.		;thermal conductivity
  params[2]=v		;wind speed
  params[3]=ta		;air temperature
  params[4]=rh		;relative humidity
  params[5]=rho		;air density
  params[6]=sw		;short wave heat flux
  params[7]=tw		;water temperature
  params[8]=cc		;cloud-cover
  params[9]=ccc		;cloud-cover coefficient

  params[0]=d		;ice thickness

  ;use bisection to determine surface temperature:
  ts=bisection( "tserr", 300., 200., 0.001, params=params)

  ;q=params[10]
  if arg_present(q) then q=params[1]*(tserr(ts, params=params)-tw+ts)/d

  return, ts

end


