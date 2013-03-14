@read_PolIce
@read_thickness

;goto, legend

rlist1=strarr(8)
rlist1[0]='20070313_TrackAndCircles/Kruunupyy_to_P3X/07216200.r31'
rlist1[1]='20070312_TrackAndCircles/Marjaniemi_to_P2X/07117100.r31'
rlist1[2]='20070312_TrackAndCircles/Marjaniemi_to_P7X/07117410.r31'
rlist1[3]='20070312_TrackAndCircles/P2X_to_Marjaniemi/07117250.r31'
rlist1[4]='20070312_TrackAndCircles/P7X_to_Marjaniemi/07116520.r31'
rlist1[5]='20070313_TrackAndCircles/P2A_to_P4X/07217260.r31'
rlist1[6]='20070313_TrackAndCircles/P3X_to_P2A/07217010.r31'
rlist1[7]='20070313_TrackAndCircles/P4X_to_P2A/07217460.r31'

;base_dir="/smiles_local/pmills/Police2007/"

lonmin=23.3
lonmax=24.8
latmin=64.7
latmax=65.4

;spawn, "ls "+base_dir+"*/*/*.r31", rlist1
;spawn, "ls "+base_dir+"*/*.dat", tlist
spawn, "ls ice_thickness/*.dat", tlist

set_plot, "ps"
device, filename="EM-bird_flights.lf.ps", /color, bits_per_pixel=8
;device, /encapsulated
device, /landscape, xsize=27.7, ysize=19.0, xoffset=1, yoffset=28.7
device, /times
;device, /landscape, xsize=27.7, ysize=19.0, xoffset=1, yoffset=28.7

tvlct, [160, 150, 200, 100], [130, 170, 0, 100], [70, 250, 0, 100], 1

;grey=4

plot, findgen(10), xrange=[lonmin, lonmax], yrange=[latmin, latmax], /nodata, $
		xstyle=1, ystyle=1, xtitle="longitude [deg.]", ytitle="latitude [deg.]", $
		title="Pol-Ice 2007 campaign E-M bird flights", $
		charsize=1.7, charthick=2, font=0

;high resolution coast-line:
nc=ncdf_open("/rinax/pmills/landmask/NAVO-lsmask-world8-var.dist5.5.nc")
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

;contour, dst, lon[lonind_min:lonind_max], lat[latind_min:latind_max], levels=[0.1, 1], $
;		c_thick=5, /overplot, c_color=[1, 2]
contour, dst, lon[lonind_min:lonind_max], lat[latind_min:latind_max], levels=0.5, $
		c_thick=5, /overplot, color=1


red=rgbtoindex(200, 0, 0)
blue=rgbtoindex(0, 0, 200)
grey=rgbtoindex(100, 100, 100)

n1=100
;yellow to red:
cind1=rgbtoindex(255, 240-findgen(n1)*240/(n1-1), 0)

n2=100
;green to blue:
cind2=rgbtoindex(0, 255-findgen(n2)*255/(n2-1), findgen(n2)*255/(n2-1))

;create the color table:
r=[intarr(n1)+255, intarr(n2), 100]
g=[240-findgen(n1)*240/(n1-1), 255-findgen(n2)*255/(n2-1), 100]
b=[intarr(n1), findgen(n2)*255/(n2-1), 100]

;tvlct, fix(r), fix(g), fix(b), 1

cind1=indgen(n1)+1
cind2=indgen(n2)+n1+1
grey=n1+n2+1

;number of days in campaign:
ndays=3.

;start of campaign:
campaign_start=convert_date('2007/3/11')

thick1=10
thick2=3

for i=0, n_elements(rlist1)-1 do begin
  print, rlist1[i]
  ;date=long(strmid(rlist1[i], strlen(base_dir), 8))
  date=long(strmid(rlist1[i], 0, 8))
  t0={time_str}
  t0.year=date/10000L
  t0.day=date mod 100
  t0.month=(date-t0.year*10000L)/100

  if time_compare(t0, campaign_start) lt 0 then continue
 
  read_PolIce, "corrected/"+rlist1[i], t0, tbs, lon, lat, t, err=err
  if err ne 0 then continue

  if time_compare(t[0], campaign_start) lt 0 then begin
    oplot, lon, lat, color=grey
  endif else begin
    ;convert dates to number of days from start of campaign:
    day=day_value(t[0], campaign_start)
    if day gt 3 then begin
      oplot, lon, lat, color=grey
    endif else begin
      cind=cind1[day*n1/ndays]
      ;oplot, lon, lat, color=cind
      oplot, lon, lat, color=2, thick=thick1
      n=n_elements(t)
      arrow, lon[n/4], lat[n/4], lon[n/4+10], lat[n/4+10], /data, color=2, thick=2
      arrow, lon[n/2], lat[n/2], lon[n/2+10], lat[n/2+10], /data, color=2, thick=2
      arrow, lon[3*n/4], lat[3*n/4], lon[3*n/4+10], lat[3*n/4+10], /data, color=2, thick=2
    endelse
  endelse

endfor

for i=0, n_elements(tlist)-1 do begin
  print, tlist[i]
  read_thickness, tlist[i], thickness, lon, lat, t

  ;convert dates to number of days from start of campaign:
  day=day_value(t[0], campaign_start)

  cind=cind2[day*n2/ndays]
  ;oplot, lon, lat, color=cind, thick=2
  oplot, lon, lat, thick=thick2, color=3
  n=n_elements(t)
  arrow, lon[n/4], lat[n/4], lon[n/4+10], lat[n/4+10], /data, color=3, thick=2
  arrow, lon[n/2], lat[n/2], lon[n/2+10], lat[n/2+10], /data, color=3, thick=2
  arrow, lon[3*n/4], lat[3*n/4], lon[3*n/4+10], lat[3*n/4+10], /data, color=3, thick=2

endfor

;plot the locations of the stations:
round_sym, 7

annot_thick=5
annot_size=1.5
annot_col=1

;Marjaniemi:
x=24.6
y=65.045
plots, x, y, /data, psym=8, thick=2
xyouts, x, y, "Marjaniemi", /data, color=annot_col, $
		charthick=annot_thick, charsize=annot_size

;P2X:
x=23.77
y=65.15
plots, x, y, /data, psym=8, thick=2
xyouts, x, y, "P2X", /data, color=annot_col, $
		charthick=annot_thick, charsize=annot_size

;P7X:
x=23.45
y=64.935
plots, x, y, /data, psym=8, thick=2
xyouts, x, y, "P7X", /data, color=annot_col, $
		charthick=annot_thick, charsize=annot_size

;P2A
x=24.35
y=64.8
plots, x, y, /data, psym=8, thick=2
xyouts, x, y, "P2A", /data, color=annot_col, $
		charthick=annot_thick, charsize=annot_size

;P4X:
x=24.02
y=65.34
plots, x, y, /data, psym=8, thick=2
xyouts, x, y, "P4X", /data, color=annot_col, $
		charthick=annot_thick, charsize=annot_size

;P3X:
x=23.92
y=65.28
plots, x, y, /data, psym=8, thick=2
xyouts, x, y, "P3X", /data, color=annot_col, $
		charthick=annot_thick, charsize=annot_size

;P7X:
x=23.45
y=64.935
plots, x, y, /data, psym=8, thick=2
xyouts, x, y, "P7X", /data, color=annot_col, $
		charthick=annot_thick, charsize=annot_size

oplot, [24.4, 24.475], [65.32, 65.32], thick=10, color=2
oplot, [24.4, 24.475], [65.28, 65.28], thick=3, color=3

round_sym, 3
plots, 24.435, 65.23, psym=8, thick=2, /data

annot_thick2=2
annot_size2=1.6
xyouts, 24.5, 65.32, "radiometer", $
		charthick=annot_thick2, charsize=annot_size2
xyouts, 24.5, 65.28, "E-M bird", $
		charthick=annot_thick2, charsize=annot_size2
xyouts, 24.5, 65.225, "waypoint", $
		charthick=annot_thick2, charsize=annot_size2

device, /close_file

legend:

;plot the legend:
device, filename="lband_leg.eps"
device, /portrait, xsize=4.5, ysize=6, /inches

ytickname=time_string(time_array(campaign_start, convert_date('0/1/1'), ndays+1))

contour, rebin(transpose(findgen(n1)), 2, n1), [0, 1], findgen(n1), $
		title="radiometer", xstyle=4, ytitle=date, $
		ytickv=findgen(4)*33, ytickname=ytickname, ystyle=1, yticks=4, $
		levels=findgen(n1), c_color=indgen(n1)+1, /fill, $
		charsize=2, charthick=2, xmargin=[12, 10]

device, /close_file

device, filename="embird_leg.eps"
contour, rebin(transpose(findgen(n2)), 2, n2), [0, 1], findgen(n1), $
		title="E-M bird", xstyle=4, ytitle=date, $
		ytickv=findgen(4)*33, ytickname=ytickname, ystyle=1, yticks=4, $
		levels=findgen(n2), c_color=indgen(n2)+n1+1, /fill, $
		charsize=2, charthick=2, xmargin=[12, 10]

device, /close_file

stop

device, /close
set_plot, "x"

end

