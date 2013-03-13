
pro get_landmask, xrange, yrange, dst, lon, lat, fname=fname
  if n_elements(fname) ne 1 then fname="../data/NAVO-lsmask-world8-var.dist5.5.nc"

  nc=ncdf_open("NAVO-lsmask-world8-var.dist5.5.nc")
  ncdf_varget, nc, "lon", lon
  ncdf_varget, nc, "lat", lat

  nlon=n_elements(lon)
  nlat=n_elements(lat)

  lonind_min=long(interpol(findgen(nlon), lon, lonmin))
  lonind_max=long(interpol(findgen(nlon), lon, lonmax))+1
  latind_min=long(interpol(findgen(nlat), lat, latmax))
  latind_max=long(interpol(findgen(nlat), lat, latmin))+1

  ncdf_varget, nc, "dst", dst, count=[lonind_max-lonind_min+1, latind_max-latind_min+1], $
                offset=[lonind_min, latind_min]
  ncdf_close, nc

  lon=lon[lonind_min:lonind_max]
  lat=lat[latind_min:latind_max]

end

pro plot_landmask, color=color, thick=thick, fname=fname
  
  get_landmask, !x.crange, !y.crange, dst, lon, lat

  contour, dst, lon[lonind_min:lonind_max], lat[latind_min:latind_max], levels=0.5, $
                c_thick=thick, /overplot, color=color
end

