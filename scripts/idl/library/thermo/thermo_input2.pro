@interpolate_ncep_flux
@wv_formulas

;get the input to the thermodynamic model:

function thermo_input2, t0, lon, lat, nt, path=path
  if n_elements(path) ne 1 then path="/vapor/pmills/ncep_flux/"

  ;use 2 m height:
  h=2.
  alpha=0.7
  solar=1366.

  date0={time_class, t0, 0}
  datef={time_class, t0+nt, 0}
  sixhour=time_class_convert('6:0')

  ntot=time_class_divide(time_class_subtract(datef, date0), sixhour)

  np=n_elements(lon)

  data=replicate({date:{time_class}, rh:0., t2m:0., v2m:0., cc:0., p2m:0., sw:0.}, $
	  	np, ntot)

  k=0
  for i=0L, nt-1 do begin
    print, "For year, ", i, ":"
    tstr=string(t0+i, format="(i4.4)")
    print, "Getting surface temperature..."
    t2m=interpolate_flux_tseries(path+"air.2m.gauss."+tstr+".nc", lon, lat, date)
    nper=n_elements(date)

    for j=0, np-1 do begin
      data[j, k:k+nper-1].date=reform(date, 1, nper)
    endfor

    data[*, k:k+nper-1].t2m=t2m

    print, "Getting surface pressure..."
    data[*, k:k+nper-1].p2m=interpolate_flux_tseries(path+"pres.sfc."+tstr+".nc", lon, lat)
    print, "Getting cloud cover..."
    data[*, k:k+nper-1].cc=interpolate_flux_tseries(path+"tcdc.eatm.gauss."+tstr+".nc", lon, lat)/100

    print, "Getting surface winds..."
    u2m=interpolate_flux_tseries(path+"uwnd.10m.gauss."+tstr+".nc", lon, lat)
    v2m=interpolate_flux_tseries(path+"vwnd.10m.gauss."+tstr+".nc", lon, lat)

    data[*, k:k+nper-1].v2m=sqrt(u2m*u2m+v2m*v2m)

    print, "Getting surface humidity..."
    q2m=interpolate_flux_tseries(path+"shum.2m.gauss."+tstr+".nc", lon, lat)

    data[*, k:k+nper-1].rh=q2m*1.606*data[*, k:k+nper-1].p2m/eq_vp(data[*, k:k+nper-1].t2m)

    swfrac=fltarr(np, nper)
    for j=0, np-1 do begin
      swfrac[j, *]=sw_flux(lat[j], date)
    endfor

    print, "Generating shortwave flux.."
    data[*, k:k+nper-1].sw=solar*swfrac*(1-0.62*data[*, k:k+nper-1].cc)*(1-alpha)

    ;stop

    k=k+nper

  endfor

return, data

end

