
pro tsmod_const_weather, ta, v, rh, cc, qs, d, savg, s, nd,  $
		sw=sw, flux=q, latent=l, rho=rho, growth=g, ts=ts

  ;implements a dirt-simple growth model based on the
  ;heat flux for different thicknesses but constant weather conditions...

  ;ta=270.		;air temperature
  ;v=2.		;wind speed
  ;rh=0.5		;relative humidity
  ;cc=0.2		;cloud-cover
  ;qs=0.		;short wave flux

  if n_elements(l) ne 1 then l=333.55	;latent heat of freezing (J/g)
  if n_elements(rho) ne 1 then rho=0.92	;specific mass of ice

  if n_elements(sw) ne 1 then sw=35.		;salinity of parent water

  if n_elements(nd) eq 0 then nd=100
  d=3*(findgen(nd)+1)/nd
  q=fltarr(nd)

  ts=fltarr(nd)

  for i=0L, nd-1 do begin
    ts[i]=derive_surf_temp(d[i], ta, v, rh, cc, qs, qstar=q1)
    q[i]=q1
  endfor

  ;use the heat flux to calculate saline content:

  ;first we convert to growth rates (cm/s):
  g=-q/l/10000./rho

  ;salinity is a function of growth rate (Cox & Weeks, 1975, 1988)
  k=fltarr(nd)
  ind=where(g gt 3.6e-5, cnt)
  if cnt gt 0 then k[ind]=0.26/(0.26+0.74*exp(-7243*g[ind]))
  ind=where(g lt 3.6e-5 and g gt 2.e-6, cnt)
  if cnt gt 0 then k[ind]=0.8925+0.0568*alog(g[ind])
  ind=where(g lt 2.e-6, cnt)
  if cnt gt 0 then k[ind]=0.12

  s=k*sw

  ;finally, we calculate the bulk salinity across different thicknesses:
  dz=d-[0, d[0:nd-2]]
  savg=total(dz*s, /cum)/d

end


