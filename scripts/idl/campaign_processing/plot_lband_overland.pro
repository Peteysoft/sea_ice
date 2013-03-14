;restore, "smos_v_asi.idlsav"

goto, plt

emh=tbs40[0, *]/skt40
pol=(tbs40[1, *]-tbs40[0, *])/(tbs40[1, *]+tbs40[0, *])

nh=230
np=200

n=n_elements(pol)

hist=lonarr(np, nh)
histt=lonarr(np, nh)

;land:
;first we need the "high-resolution" landmask:
nc=ncdf_open("NAVO-lsmask-world8-var.dist5.5.nc")
ncdf_varget, nc, "lon", lon
ncdf_varget, nc, "lat", lat

nlon=n_elements(lon)
nlat=n_elements(lat)

lonind=interpol(findgen(nlon), lon, lon40)
latind=interpol(findgen(nlat), lat, lat40)

lonind_min=long(min(lonind))
lonind_max=long(max(lonind)+1)

latind_min=long(min(latind))
latind_max=long(max(latind)+1)

ncdf_varget, nc, "dst", dst, count=[lonind_max-lonind_min+1, latind_max-latind_min+1], $
		offset=[lonind_min, latind_min]
ncdf_close, nc

mask=interpolate(float(dst), lonind-lonind_min, latind-latind_min)

ind1=where(mask lt 0.5 and zangle40 gt 36 and zangle40 le 40, cnt)

  for i=0L, cnt-1 do begin
    j=ind1[i]
    ;indy=emh[j]*nh
    indy=tbs40[0, j]-50
    indx=(pol[j]+1)*np/2
    if indy ge nh then continue
    if indx ge np then continue
    if indx lt 0 then continue
    hist[indx, indy]=hist[indx, indy]+1
  endfor

ind2=where(zangle40 gt 36 and zangle40 le 40, cnt)

  for i=0L, cnt-1 do begin
    j=ind2[i]
    ;indy=emh[j]*nh
    indy=tbs40[0, j]-50
    indx=(pol[j]+1)*np/2
    if indy ge nh then continue
    if indx ge np then continue
    if indx lt 0 then continue
    histt[indx, indy]=histt[indx, indy]+1
  endfor

plt:

ps_start, "land_hist.ps", /color, /square

nlev=10
loadct, 13, bottom=1, ncolors=nlev
tvlct, 255, 255, 255, 1

xleg_min=0.3
xleg_max=0.325
yleg_min=170
yleg_max=280
;yleg_min=170
;yleg_max=270
dyleg=(yleg_max-yleg_min)/(nlev)

pol1=2*findgen(np)/(np-1)-1
tb=findgen(nh)+50

cind=indgen(nlev)+1
  maxn=max(hist)
  l2=maxn/10/(nlev-2)
  levels=[0, 1, l2, (lindgen(nlev-3)+1)*maxn/(nlev-2)]
  contour, hist, pol1, tb, xrange=[0, 0.4], $
                levels=levels, title="L-band radiometric signatures over land", $
                c_color=cind, /fill, font=0, $
                xtitle="Polarisation", ytitle="Tbh [K]"
  contour, histt, pol1, tb, levels=[1, 18, 250, 1000, 2500], c_label=indgen(5)+1, /overplot
  ;plot the legend:
  for j=0, nlev-1 do begin
    x=[xleg_min, xleg_max, xleg_max, xleg_min]
    y=yleg_min+[dyleg*j, dyleg*j, dyleg*(j+1), dyleg*(j+1)]
    polyfill, x, y, color=cind[j], /data
    xyouts, xleg_max+0.014, yleg_min+dyleg*j, /data, $
                string(levels[j], format="(i4.1)"), font=0, charsize=0.8
  endfor
  oplot, [0.29, 0.380, 0.380, 0.29, 0.29], $
                [yleg_min-5, yleg_min-5, yleg_max+5, yleg_max+5, $
                yleg_min-5]

ps_end


end


