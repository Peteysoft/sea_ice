@time_lib

restore, "smos_v_asi.idlsav"

;land:
;first we need the "high-resolution" landmask:
nc=ncdf_open("/rinax/pmills/landmask/NAVO-lsmask-world8-var.dist5.5.nc")
ncdf_varget, nc, "lon", lon
ncdf_varget, nc, "lat", lat

nlon=n_elements(lon)
nlat=n_elements(lat)

lonind=interpol(findgen(nlon), lon, [20.0, 26])
latind=interpol(findgen(nlat), lat, [61, 66.])

lonind_min=long(lonind[0])
lonind_max=long(lonind[1]+1)
latind_min=long(latind[1])
latind_max=long(latind[0]+1)

;lonind_min=long(min(lonind))
;lonind_max=long(max(lonind[where(lon40 lt 360)])+1)

;latind_min=long(min(latind[where(lat40 lt 90)]))
;latind_max=long(max(latind[where(lat40 gt 0)])+1)

lonind=interpol(findgen(nlon), lon, lon40)
latind=interpol(findgen(nlat), lat, lat40)

ncdf_varget, nc, "dst", dst, count=[lonind_max-lonind_min+1, latind_max-latind_min+1], $
		offset=[lonind_min, latind_min]
ncdf_close, nc

mask=interpolate(float(dst), lonind-lonind_min, latind-latind_min)

;make a nice little plot showing the overland flights:
ps_start, "flight_tracks5a.lf.bw.ps", /portrait;, /color

;set color table:
tvlct, [250, 230, 200, 200, 0, 246], [0, 200, 220, 0, 200, 214], [0, 150, 250, 250, 150, 160], 1

red=1
red=0
light_brown=2
light_brown=240
light_blue=3
light_blue=255
ice_blue=4
ice_blue=0
water_blue=5
water_blue=0

mask_brown=6
mask_brown=250

contour, dst, lon[lonind_min:lonind_max], lat[latind_min:latind_max], levels=[0, 1], $
		xstyle=1, ystyle=1, title="Pol-Ice EMIRAD data", $
		xtitle="longitude [deg.]", ytitle="latitude [deg.]", $
		c_color=[light_brown, light_blue], /fill, charsize=1.5, $
		xmargin=[6, 1]
contour, dst, lon[lonind_min:lonind_max], lat[latind_min:latind_max], levels=1, $
		c_thick=2, /overplot
ind=where(mask gt 0.5)
oplot, lon40[ind], lat40[ind], psym=3, color=ice_blue
ind=where(asi40 lt 1)
oplot, lon40[ind], lat40[ind], psym=3, color=water_blue

ind=where(mask le 0.5)
oplot, lon40[ind], lat40[ind], psym=3, color=red

xleg_min=23.5
xleg_max=24.7
yleg_min=61.3
yleg_max=61.8
xmarg=0.2
legdist=0.1
ymarg=0.1
;polyfill, [xleg_min, xleg_max, xleg_max, xleg_min], $
;		[yleg_min, yleg_min, yleg_max, yleg_max], color=255
;oplot, [xleg_min, xleg_max, xleg_max, xleg_min], $
;		[yleg_min, yleg_min, yleg_max, yleg_max]
;draw_legend, [xleg_min+xmarg, xleg_max+xmarg+legdist], [yleg_max-ymarg, yleg_min-ymarg], $
;		["  ice", "  land", "  open water"], $
;		psym=3, color=[ice_blue, red, water_blue]

;do this the easy way:
x0leg=20.3
ycorr=0.4
;xyouts, x0leg, 65.25+ycorr, "ice", color=ice_blue, charsize=1.4, charthick=4
;xyouts, x0leg, 65.1+ycorr, "land", color=red, charsize=1.4, charthick=4
;xyouts, x0leg, 64.95+ycorr, "open water", color=water_blue, charsize=1.4, charthick=4

lonleg=22.56
lonleg1=22.85
latleg0=63.1
dlatleg=-1.95
annotsize=0.65
;annotsize=0.5

;we also want to label the flight tracks:
uind=uniq(fno40)
nf=n_elements(uind)

;mask out the bit containing the legend:
polyfill, [22.55, 25.8, 25.8, 22.55], [61.1, 61.1, 63.3, 63.3], color=mask_brown

xyouts, lonleg1, 63.2, "   date               number    name", charsize=annotsize, charthick=2

for i=0, nf-1 do begin
  pathfields=strsplit(flist[fno40[uind[i]]], '/', /extract)

  flightname=pathfields[2]
  flightname=byte(flightname)
  ind=where(flightname eq (byte("-"))[0] or $
		flightname eq (byte("_"))[0], cnt)
  if cnt ne 0 then flightname[ind]=byte(' ')
  flightname=string(flightname)
  flightdate=date40[uind[i]]
  flightdate.second=0
  flightno=(strsplit(pathfields[3], '.', /extract))[0]

  annot=time_string(flightdate)+", "+flightno+", "+flightname

  xyouts, lon40[uind[i]], lat40[uind[i]], string(i, format="(i3)"), charsize=annotsize
  xyouts, lonleg, latleg0+dlatleg*i/(nf-1), /data, charsize=annotsize, $
		string(i+1, format="(i3)")+"."
  xyouts, lonleg1, latleg0+dlatleg*i/(nf-1), annot, charsize=annotsize
endfor

ps_end



end

