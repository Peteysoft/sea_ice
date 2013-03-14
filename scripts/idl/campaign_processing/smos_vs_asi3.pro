@read_PolIce
@read_thickness

goto, plt

;restore, "smos_v_asi.idlsav"
n=n_elements(skt40)
;em40=tbs40/rebin(transpose(skt40), 4, n)

;make a histogram of emissivities:
nh=200
np=200
;emh=findgen(nh)/(nh-1)
tb=findgen(nh)+50
pol=2*findgen(np)/(np-1)-1
hist=lonarr(4, nh, np)

flag=bytarr(4, n)
;all data:
flag[0, *]=zangle40 gt 36 and zangle40 lt 40 and finite(asi40)
;100% ice concentration:
flag[1, *]=zangle40 gt 36 and zangle40 lt 40 and asi40 gt 99 $
		and finite(asi40)
;open water:
flag[2, *]=zangle40 gt 36 and zangle40 lt 40 and asi40 lt 1 $
		and finite(asi40)
;fractional ice concentration:
flag[3, *]=zangle40 gt 36 and zangle40 lt 40 and asi40 gt 1 and asi40 lt 99 $
		and finite(asi40)

polf=(tbs40[1, *]-tbs40[0, *])/(tbs40[1, *]+tbs40[0, *])

for k=0, 3 do begin
  ind1=where(flag[k, *], cnt)
  for i=0L, cnt-1 do begin
    j=ind1[i]
    ;indx=em40[0, j]*nh
    indx=tbs40[0, j]-50
    indy=tbs40[1, j]-100
    ;indy=(polf[j]+1)*np/2
    if indx ge nh then continue
    if indy ge np then continue
    if indy lt 0 then continue
    hist[k, indx, indy]=hist[k, indx, indy]+1
  endfor
  
endfor

plt:

title=strarr(4)
title[0]="Pol-Ice 2007 L-Band Measurements"
title[1]='100% Ice Concentration'
title[2]='0% Ice Concentration'
title[3]="Fractional Ice Concentration"

set_plot, "ps"
device, filename="tb_hist.ps"
device, /landscape, xsize=9, ysize=6.5, xoffset=1, yoffset=10, /inches
device, /color, bits=8
;device, set_font="Times"
device, /times, /isolatin1

;nlev=255
nlev=10
loadct, 13, bottom=1, ncolors=nlev
;tvlct, [255, indgen(nlev-1)*255/(nlev-2)], [255, intarr(nlev-1)], $
;		[255, 255-indgen(nlev-1)*255/(nlev-1)], 1
tvlct, 255, 255, 255, 1

!p.multi=[0, 2, 2]

;xleg_min=0.3
;xleg_max=0.325
xleg_min=270
xleg_max=275
yleg_min=60
yleg_max=170
;yleg_min=170
;yleg_max=270
dyleg=(yleg_max-yleg_min)/(nlev)

cind=indgen(nlev)+1

for i=0, 3 do begin
  maxn=max(hist[i, *, *])
  l2=maxn/10/(nlev-2)
  levels=[0, 1, l2, (lindgen(nlev-3)+1)*max(hist[i, *, *])/(nlev-2)]
  contour, transpose(hist[i, *, *]), tb+50, tb, $ xrange=[0, 0.4], $
		levels=levels, title=title[i], $
		c_color=cind, /fill, font=0, $
		xtitle="Tbv [K]", ytitle="Tbh [K]"
  ;plot the legend:
  for j=0, nlev-1 do begin
    x=[xleg_min, xleg_max, xleg_max, xleg_min]
    y=yleg_min+[dyleg*j, dyleg*j, dyleg*(j+1), dyleg*(j+1)]
    polyfill, x, y, color=cind[j], /data
    xyouts, xleg_max+3, yleg_min+dyleg*j, /data, $
		string(levels[j], format="(i4.1)"), font=0, charsize=0.8
  endfor
;  oplot, [0.29, 0.380, 0.380, 0.29, 0.29], $
  oplot, [xleg_min-5, xleg_max+20, xleg_max+20, xleg_min-5, xleg_min-5], $
		[yleg_min-5, yleg_min-5, yleg_max+5, yleg_max+5, $
		yleg_min-5]
endfor

!p.multi=[0, 1, 1]
device, /close_file
device, /close
set_plot, "x"

end


