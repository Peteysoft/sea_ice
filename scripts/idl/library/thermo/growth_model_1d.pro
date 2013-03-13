@derive_surf_temp

;more sophisticated growth model that takes an actual time series:

function growth_rate, t, d
  common growth_model_common, data, tw, l, rho, eps

  nt=(size(data, /dim))[1]
  n=n_elements(d)

  g=fltarr(n)

  for i=0L, n-1 do begin
    ta=interpol(data[i, *].t2m, findgen(nt), t)
    v=interpol(data[i, *].v2m, findgen(nt), t)
    rh=interpol(data[i, *].rh, findgen(nt), t)/100.
    cc=interpol(data[i, *].cc, findgen(nt), t)/100.
    qs=interpol(data[i, *].sw, findgen(nt), t)

    if d[i] lt eps then begin
      q1=net_heat_flux(tw, ta, v, rh, cc)+qs
    endif else begin
      ts=derive_surf_temp(d[i], ta, v, rh, cc, qs, qstar=q1)
      print, q1
    endelse

    ;convert to a growth rate in cm/s:
    g[i]=-q1/l/10000./rho
  endfor

  ;stop

  ;convert to m/day:
  return, g*864

end

function growth_model_1d, data1, tw=tw1, dt=dt

  common growth_model_common, data, tw, l, rho, eps

  ;some constants:
  l=333.55        ;latent heat of freezing (J/g)
  rho=0.92        ;specific mass of ice
  eps=1e-7

  data=data1
  dim=size(data, /dim)
  n=dim[0]
  ntgrid=dim[1]

  d=fltarr(n)

  if n_elements(tw1) eq 0 then tw=271. else tw=tw1
  if n_elements(dt) eq 0 then dt=0.25
  nt=ntgrid

  dall=fltarr(n, nt)

  for i=0, nt-1 do begin
    print, i
    t=i*dt
    dddt=growth_rate(t, d)
    dnew=rk4(d, dddt, t, dt, "growth_rate")
    ;stop
    ind=where(dnew lt 0, cnt)
    if cnt gt 0 then dnew[ind]=0.
    ;print, dnew
    dall[*, i]=dnew
    d=dnew
  endfor

  return, dall

end

