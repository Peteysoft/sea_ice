;plot all of the laser altimetry data:
spawn, "ls laser/*.fb", flist

latmin=64.5
latmax=65.5
lonmin=22.5
lonmax=25.

ps_start, "map_laser.ps", /color

plot, findgen(10), xrange=[lonmin, lonmax], yrange=[latmin, latmax], $
		title="Pol-Ice 2007 Laser altimetry data", $
		xtitle="lon", ytitle="lat"

tvlct, [250], [150], [50], 1

brown=1

loadct, 13, ncolor=250, bottom=2

round_sym, 0.5

;high-resolution coastline:
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

contour, dst, lon[lonind_min:lonind_max], lat[latind_min:latind_max], levels=0.5, $
                c_thick=5, /overplot, color=brown

resample=20
nf=n_elements(flist)

xleg0=22.6
yleg0=64.8
dyleg=-0.35/nf

for i=0, nf-1 do begin
  print, i
  read_laser, flist[i], z, lon, lat, t
  n=n_elements(z)
  sub=lindgen(n/resample)*resample

  lon1=lon[sub]
  lat1=lat[sub]
  z1=(smooth(z, resample))[sub]
  plots, lon1, lat1, color=fix(z1*100)+20, psym=3
  plots, lon[0], lat[0], psym=2
  xyouts, lon[0], lat[0], string(i, format="(i3)")

  xyouts, xleg0, yleg0+dyleg*i, flist[i], charsize=0.4

  ;stop
endfor

ps_end

end

