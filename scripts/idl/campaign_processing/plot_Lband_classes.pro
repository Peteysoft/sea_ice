restore, "smos_v_asi.idlsav"

goto, plt

emh=tbs40[0, *]/skt40
pol=(tbs40[1, *]-tbs40[0, *])/(tbs40[1, *]+tbs40[0, *])

ncls=13
nh=100
np=200

n=n_elements(pol)

hist=lonarr(ncls, np, nh)

flag=bytarr(ncls, n)

;0% ice concentration:
lflag=zangle40 gt 36 and zangle40 le 40
flag[0, *]=asi40 le 1 and lflag
;1-10% ice concentration:
flag[1, *]=asi40 gt 1 and asi40 le 10 and lflag

;all fractional ice concentrations from 10-20 to 80-90%:
for i=2, 9 do begin
  flag[i, *]=asi40 gt (i-1)*10 and asi40 le i*10 and lflag
endfor

;80-99% ice concentration:
flag[10, *]=asi40 gt 80 and asi40 le 99 and lflag
;1-10% ice concentration:
flag[11, *]=asi40 gt 99 and lflag

;land:
;first we need the "high-resolution" landmask:
nc=ncdf_open("NAVO-lsmask-world8-var.dist5.5.nc")
ncdf_varget, nc, "lon", lon
ncdf_varget, nc, "lat", lat

nlon=n_elements(lon)
nlat=n_elements(lat)

lonind=interpol(findgen(nlon), lon, [20.0, 25])
latind=interpol(findgen(nlat), lat, [61, 65.5])

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

flag[12, *]=mask lt 0.5 and lflag

for k=0, ncls-1 do begin
  ind1=where(flag[k, *], cnt)
  for i=0L, cnt-1 do begin
    j=ind1[i]
    indy=emh[j]*nh
    ;indx=tbs40[0, j]-50
    indx=(pol[j]+1)*np/2
    if indy ge nh then continue
    if indx ge np then continue
    if indx lt 0 then continue
    hist[k, indx, indy]=hist[k, indx, indy]+1
  endfor

endfor

plt:

;set color table:
nscale=20
r=[240-indgen(nscale)*155/(nscale-1), intarr(ncls), 50]
g=[240-indgen(nscale)*155/(nscale-1), 255-indgen(ncls)*255/(ncls-1), 50]
b=[240-indgen(nscale)*155/(nscale-1), indgen(ncls)*255/(ncls-2), 0]

ps_start, "conc_ret.ps", /color
tvlct, r, g, b, 1
loadct, 13, bottom=nscale+1, ncol=ncls
;tvlct, 50, 50, 0, nscale

x=2*findgen(np)/(np-1)-1
y=findgen(nh)/(nh-1)
histt=reform(total(hist, 1))
ntt=total(histt)
maxh=max(histt)
levels=[1, indgen(nscale-1)*maxh/(nscale-1)+maxh/(nscale-1)/10]
contour, histt, x, y, levels=levels, xrange=[0, 0.4], $
		c_color=indgen(nscale)+2, /cell_fill

for i=0, ncls-1 do begin
  hist1=hist[i, *, *]
  sind=sort(hist1)
  nhist=total(hist1[sind], /cum)
  per=fltarr(np, nh)
  nt=nhist[n_elements(nhist)-1]
  per[sind]=1.*nhist/nt
  members=hist[i, *, *] gt 1.*histt/ncls
  contour, members, x, y, levels=0.5, /overplot, c_color=i+nscale+1, /fill
endfor

ps_end

ps_start, "conc_ret2.ps", /color

mind=intarr(np, nh)
norm=total(hist1[0, *, *])
hist1=hist
;hist1[0, *, *]=hist[0, *, *]*total(hist[11, *, *])/total(hist[0, *, *])
nzind=where(hist1[0, *, *] gt 0)
;mind[nzind]=0
xind=rebin(indgen(np), np, nh)
yind=rebin(transpose(indgen(nh)), np, nh)
for i=1, ncls-1 do begin
  norm1=total(hist1[i, *, *])
  gtind=where(1.*hist1[i, *, *]/norm1 gt 1.*hist1[mind, xind, yind]/norm, ngt)
;  gtind=where(hist1[i, *, *] gt hist1[mind, xind, yind], ngt)
  if ngt gt 0 then mind[gtind]=i
  norm=norm1
endfor

zind=where(histt eq 0)
mind[zind]=-1

;contour, mind, x, y, levels=indgen(ncls), xrange=[0, 0.4], $
;		c_color=indgen(ncls)+nscale+2, /cell_fill, $
;		xtitle="polarization", ytitle="emissivity", $
;		title="Ice concentration from AMSR-E vs. L-band emission from PolIce campaign"

memb=mind eq 0
contour, memb, x, y, levels=0.5, xrange=[0, 0.4], $
		c_color=nscale+2, /cell_fill, /nodata, $
		xtitle="polarization", ytitle="emissivity", $
		title="Ice concentration from AMSR-E vs. L-band emission from PolIce campaign"

annot=["open water", '1-10%', '10-20%', '20-30%', '30-40%', '40-50%', '50-60%', $
		'60-70%', '70-80%', '80-90%', '90-99%', '100%', 'land']

for i=0, nh-2 do begin
  for j=0, np-2 do begin
    if mind[j, i] ne -1 then begin
      x1=[x[j], x[j], x[j+1], x[j+1]]
      y1=[y[i], y[i+1], y[i+1], y[i]]
      polyfill, x1, y1, color=mind[j, i]+nscale+2, /data
    endif
  endfor
endfor

for i=0, ncls-1 do begin
;  memb=mind eq i
;  contour, memb, x, y, levels=0.5, c_color=i+nscale+2, /cell_fill, /overplot
  x1=[0.3, 0.32, 0.32, 0.3]
  y1=0.55+i*0.03+[0., 0., 0.03, 0.03]
  polyfill, x1, y1, /data, color=i+nscale+2
  xyouts, 0.33, 0.56+i*0.03, annot[i], /data
endfor

ps_end

;make a nice little plot showing the overland flights:
ps_start, "overland.ps", /portrait, /color

;set color table:
tvlct, [250, 230, 200, 200, 0], [0, 200, 220, 0, 200], [0, 150, 250, 250, 150], 1

red=1
light_brown=2
light_blue=3
ice_blue=4
water_blue=5

contour, dst, lon[lonind_min:lonind_max], lat[latind_min:latind_max], levels=[0, 1], $
		xstyle=1, ystyle=1, title="Pol-Ice data by surface type", $
		xtitle="longitude [deg.]", ytitle="latitude [deg.]", $
		c_color=[light_brown, light_blue], /fill, charsize=1.4
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
xyouts, x0leg, 65.25, "ice", color=ice_blue, charsize=1.2, charthick=2
xyouts, x0leg, 65.1, "land", color=red, charsize=1.2, charthick=2
xyouts, x0leg, 64.95, "open water", color=water_blue, charsize=1.2, charthick=2

ps_end



end

