@wv_formulas

;simple thermodynamic model for determining surface temperature
;as a function of thickness:

;computes Qnet*d/k + Tw - Ts
;where:
;Qnet is surface heat flux
;d is ice thickness
;k is thermal conductivity
;Tw is water temperature
;Ts is surface temperature
function net_heat_flux, ts, ta, v, rh, cc, tw=tw, rho=rho

  if n_elements(tw) ne 1 then tw=271.
  if n_elements(rho) ne 1 then rho=1.2

  ;cloud-cover coefficient:
  ccc=0.8

  k=2.			;thermal conductivity

  cs=0.003		;transfer coeff. sensible
  ce=0.003		;transfer coeff. latent

  l=2833000.		;latent heat of sublimation [J/kg]
  cpa=1025.		;heat capacity of air

  stefbolt=5.670e-8	;stefan boltzmann const.
  em_ice=0.9		;emissivity of ice

  ;sensible heat flux:
  qh=rho*cpa*cs*v*(ta-ts)

  ;latent heat flux:
  eqa=eq_vp(Ta)/101000./1.6	;saturation vapour pressure in air
  eq0=eq_vp_ice(Ts)/101000./1.6	;vapour pressure at surface
  qe=rho*L*ce*v*(rh*eqa-eq0)

  ;longwave heat flux:
  qlw=-em_ice*stefbolt*(0.39*(1-ccc*cc*cc)*ts^4+ts^3*(ts-ta))

  ;net heat flux (excluding short-wave):
  q=qh+qe+qlw

  ;print, qh, qe, qlw
  ;stop
  ;print, "Net flux:", q

  ;stop

  return, q

end



