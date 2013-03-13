@time_class

;interpolates ncep surface fluxes for specified lon-lat
;returns the whole time series

function interpolate_flux_tseries, filename, lon, lat, t

  nc=ncdf_open(filename)
          
  ncdf_varget, nc, "lon", longrid
  ncdf_varget, nc, "lat", latgrid
  ncdf_varget, nc, "time", t1

  ncdf_attget, nc, 3, "add_offset", add_offset
  ncdf_attget, nc, 3, "scale_factor", scale_factor
  ncdf_varget, nc, 3, data1

  ncdf_attget, nc, 3, "add_offset", add_offset
  ncdf_attget, nc, 3, "scale_factor", scale_factor

  data=data1*scale_factor+add_offset

  ;stop
  ncdf_close, nc

  nlon=n_elements(longrid)
  nlat=n_elements(latgrid)
  nt=n_elements(t1)

  xind=interpol(findgen(nlon), longrid, lon)
  yind=interpol(findgen(nlat), latgrid, lat)

  np=n_elements(lon)

  xind1=rebin(xind, np, nt)
  yind1=rebin(yind, np, nt)
  tind1=rebin(transpose(findgen(nt)), np, nt)

  result=interpolate(data, xind1, yind1, tind1)

  if arg_present(t) then begin
    t=time_class_convert(t1/24, ref_date={time_class, 1, 3})
  endif

  ;stop

  return, result

end
		  
