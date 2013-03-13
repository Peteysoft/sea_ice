@wv_formulas
@ice_model_util
@plane_parallel

;calculates brightness temperature as function of ice thickness:
pro tbvsd_func, lon, lat, date, d, tbv, tbh, t0, nu=nu, theta=theta, $
		tbsky=tbsky, sw=sw, tw=tw, scoef=scoef, tsmax=tsmax, $
		tsurf=tsurf, allbound=allbound, next_gen=next_gen

  if n_elements(nu) ne 1 then nu=1.4
  if n_elements(theta) ne 1 then theta=35.*!dtor

  alpha=0.7		;shortwave surface albedo
  solar=1366.

  if n_elements(d) eq 0 then begin
    nn=81L
    d=findgen(nn)/20
  endif

  ;get relevant parameters:
  ;surface temperature:
  get_ecmwf_profile, "TT", lon, lat, date, tprof
  ;humidity:
  get_ecmwf_profile, "SQ", lon, lat, date, qprof
  ;to convert to relative humidity:
  get_ecmwf_profile, "PP", lon, lat, date, pprof
  ;cloud cover:
  get_ecmwf_profile, "CC", lon, lat, date, ccprof
  ;wind speed:
  get_ecmwf_profile, "UU", lon, lat, date, uprof
  get_ecmwf_profile, "VV", lon, lat, date, vprof

  nlev=n_elements(tprof)
  ta=tprof[nlev-1]
  rh=qprof[nlev-1]*1.6*pprof[nlev-1]/eq_vp(ta)*100

  ;stop

  ;should really incorporate layer thickness:
  cc=1.
  for i=0, nlev-1 do begin
    cc=cc*(1-ccprof[i])
  endfor
  cc=1.-cc

  v=sqrt(uprof[nlev-1]^2+vprof[nlev-1]^2)

  tbv=fltarr(nn)
  tbh=fltarr(nn)

  ;calculate short-wave heat flux:
  swf=solar*sw_flux(lat, date)*(1-0.62*cc)*(1-alpha)

  print, "Parameters:"
  print, "air temperature:", ta
  print, "relative humidity:", rh
  print, "cloud cover:", cc
  print, "wind speed:", v
  print, "shortwave flux:", swf

  nz=10						;number of ice layers
  if n_elements(tw) ne 1 then tw=271.		;water temperature
  if n_elements(tbsky) ne 1 then tbsky=3.	;sky temperature
  if n_elements(sw) ne 1 then sw=35.		;water salinity

  ;coeffs. for salinity vs. thickness curve:
  if n_elements(scoef) ne 3 then begin
    ;scoef=[-0.274, 2.08, 0]
    scoef=[-2.66879, 1.80280, 0.404441]
  endif

  ;maximum value for surface temp.:
  if n_elements(tsmax) ne 1 then tsmax=272.

  t0=fltarr(nn)

  for i=0, nn-1 do begin
    dz=d[i]/(nz-1)
    z=findgen(nz)*dz

    if n_elements(tsurf) ne 1 then begin
      ;derive surface temperature:
      ts=derive_surf_temp(d[i], ta, v, rh, cc, swf, tw=tw)
      t0[i]=ts
      ;print, "surface temperature:", ts
      if ts gt tsmax then ts=tsmax
    endif else begin
      ts=tsurf
    endelse

    ;linear temperature profile:
    tprof=findgen(nz)*(tw-ts)/(nz-1)+ts

    ;create salinity profile:
    s=exp(scoef[0]*d[i]+scoef[1])+scoef[2]
    ;print, "bulk salinity:", s
    sprof=calc_sprof(s, nz)

    ;print, tprof
    ;print, sprof

    dz1=[1, fltarr(nz)+dz, 1]
    tprof1=[tbsky, tprof, tw]
    sprof1=[0, sprof, sw]
    rho=fltarr(nz+2)		;density irrelevant
    type=[0, fltarr(nz)+3, 5]
    if keyword_set(allbound) then begin
      bound=lindgen(nz+1)
    endif else begin
      bound=[0, nz]
    endelse

    plane_parallel1, nu, theta, dz1, tprof1, type, rho, sprof1, bound, tbv1, tbh1
    tbv[i]=tbv1[0]
    tbh[i]=tbh1[0]
    ;stop
  endfor

end


