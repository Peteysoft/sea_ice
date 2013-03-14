lf_flag=1

;for flagging shallow water:
max_sw_flag=0.01
;maximum separation between measurements:
;max_sep=50.
max_sep=1000.

minw=50.

eps_t="eps=4.0+0.1i"
eps_t="p5"

;fbase="model_results/model_run."+eps_t+".0101.c"
fbase="model_results/paper_rev."+eps_t

corr_vic=500.

savefile=fbase+".idlsav"
tseries_outfile=fbase+".tseries.ps"
scatter_outfile=fbase+".scatter.ps"
;plotfile=fbase+".ps"

restore, savefile

print, rfile
print, tfile
print, lfile

;rname="../campaign_data/Police2007/corrected/20070313_TrackAndCircles/P4X_to_P2A/07217460.r31"
t0='2007/3/13'
read_police, rfile, t0, tbs, lonr, latr, /cor, angle=-40

;tname="../campaign_data/Police2007/ice_thickness/20070313_P2A_P4X.dat"
read_thickness, tfile, thick, lont, latt
;lname="../campaign_data/Police2007/laser/200703131539.fb"
read_laser, lfile, z, lonl, latl

rotate_flight_track, lonr, latr, xr, yr, /new
rotate_flight_track, lont, latt, xt, yt
rotate_flight_track, lonl, latl, xl, yl

dy=yt-interpol(yr, xr, xt)
dy2=interpol(dy, xt, s)

ind_sw2=where(swflag gt max_sw_flag, cnt_sw2)
;ind_sw=where(swflag lt max_sw_flag and abs(dy2) lt max_sep and finite(tbv) and finite(tbh), cnt_sw)
ind_sw=where(swflag lt max_sw_flag and wfactor gt minw and finite(tbv) and finite(tbh), cnt_sw)

ps_start, tseries_outfile, /portrait, /color
;ps_start, plotfile, /portrait, /color

tvlct, [250, 0, 0, 190], [0, 250, 0, 190], [0, 0, 250, 190], 1
red=1
green=2
blue=3
grey=4

!p.multi=[0, 1, 3]

symsize=0.1
round_sym, 1

xmarg=[10, 8]

sw_width=10

charsize=1.8

annotsize=1.2

;rfile1=strjoin((strsplit(rfile, "/", /extract))[6:7], " ")
;tfile1=strjoin((strsplit(tfile, "/", /extract))[5], " ")

;ttlbase=rfile1+", "+tfile1+", "+eps_t

ttlbase="P2A to P3X"

plot, s, [tbv_meas, tbh_meas], title="a. "+ttlbase+": plane-parallel", $
		xtitle="distance [m]", ytitle="Tb [K]", charsize=charsize, $
		/nodata, /ynozero, yrange=[160, 300], xmargin=xmarg, xstyle=1

for i=0L, cnt_sw2-1 do begin
  s1=s[ind_sw2[i]]
  polyfill, [s1-sw_width, s1-sw_width, s1+sw_width, s1+sw_width], $
	  	[!y.crange, reverse(!y.crange)], color=grey
endfor

oplot, s, tbv_meas, psym=8, symsize=symsize
oplot, s, tbh_meas, psym=8, symsize=symsize, color=blue
oplot, s, tbv_pp, psym=1, symsize=symsize, color=green
oplot, s, tbh_pp, psym=1, symsize=symsize, color=red

x0leg=4000
y0leg=290
dyleg=10
xyouts, x0leg, y0leg, "Tbv", charsize=annotsize, charthick=annot_thick
xyouts, x0leg, y0leg-dyleg, "Tbh", color=blue, charsize=annotsize, charthick=annot_thick
xyouts, x0leg, y0leg-2*dyleg, "Tbv modelled", color=green, charsize=annotsize, charthick=annot_thick
xyouts, x0leg, y0leg-3*dyleg, "Tbh modelled", color=red, charsize=annotsize, charthick=annot_thick
axis, /yaxis, ytitle="Tb [K]", yrange=!y.crange, charsize=charsize

!y.crange=[160, 280]

plot, s, [tbv_meas, tbh_meas], title="b. "+ttlbase+": plane-parallel ensemble", $
		xtitle="distance [m]", ytitle="Tb [K]", charsize=charsize, $
		/nodata, /ynozero, yrange=!y.crange, xmargin=xmarg, xstyle=1

for i=0L, cnt_sw2-1 do begin
  s1=s[ind_sw2[i]]
  polyfill, [s1-sw_width, s1-sw_width, s1+sw_width, s1+sw_width], $
	  	[!y.crange, reverse(!y.crange)], color=grey
endfor

oplot, s, tbv_meas, psym=8, symsize=symsize
oplot, s, tbh_meas, psym=8, symsize=symsize, color=blue
oplot, s, tbv_ppe, psym=1, symsize=symsize, color=green
oplot, s, tbh_ppe, psym=1, symsize=symsize, color=red
axis, /yaxis, ytitle="Tb [K]", yrange=!y.crange, charsize=charsize

plot, s, [tbv_meas, tbh_meas], title="c. "+ttlbase+": ridged Monte Carlo", $
		xtitle="distance [m]", ytitle="Tb [K]", charsize=charsize, $
		/nodata, /ynozero, yrange=!y.crange, xmargin=xmarg, xstyle=1

for i=0L, cnt_sw2-1 do begin
  s1=s[ind_sw2[i]]
  polyfill, [s1-sw_width, s1-sw_width, s1+sw_width, s1+sw_width], $
	  	[!y.crange, reverse(!y.crange)], color=grey
endfor

oplot, s, tbv_meas, psym=8, symsize=symsize
oplot, s, tbh_meas, psym=8, symsize=symsize, color=blue
oplot, s, tbv, psym=1, symsize=symsize, color=green
oplot, s, tbh, psym=1, symsize=symsize, color=red
axis, /yaxis, ytitle="Tb [K]", yrange=!y.crange, charsize=charsize

goto, dontplotthis

plot, xt, thick, title="d. "+ttlbase+": auxiliary data", xtitle="distance [m]", $
			ytitle="ice thickness [m]", xrange=!x.crange, $
			ystyle=8, xmargin=xmarg, charsize=charsize, /nodata

sw_width=10
for i=0L, cnt_sw2-1 do begin
  s1=s[ind_sw2[i]]
  polyfill, [s1-sw_width, s1-sw_width, s1+sw_width, s1+sw_width], $
	  	[!y.crange, reverse(!y.crange)], color=grey
endfor
oplot, xt, thick

;z=z-min(z)
oplot, xl, z, color=red
axis, /yaxis, ytitle="measurement separation [m]", yrange=[min(dy), max(dy)], $
		/save, charsize=charsize
oplot, xt, dy, color=blue
oplot, !x.crange, [0, 0], linestyle=1
xyouts, 35000, 30, "ice thickness", $
		charsize=annotsize2, charthick=annot_thick
xyouts, 35000, 5, "surface height", color=red, $
		charsize=annotsize2, charthick=annot_thick
xyouts, 35000, -20, "measurement separation", color=blue, $
		charsize=annotsize2, charthick=annot_thick

dontplotthis:

;ps_end

;device, /portrait, xsize=6.5, ysize=6.5, xoffset=1, yoffset=2.5, /inches

charsize2=1.4
ps_start, scatter_outfile, /square

!p.multi=[0, 3, 3]

round_sym, 0.5

tbvall=[tbv_pp, tbv_ppe, tbv]
yrange=[min(tbv_ppe), max(tbv_ppe)]

r=correlate(tbv_meas[ind_sw], tbv_pp[ind_sw])
subtitle="r="+string(r, format="(f8.3)")
plot, tbv_meas[ind_sw], tbv_pp[ind_sw], psym=8, title="a. Plane-parallel model", $
		xtitle="Tbv (measured) [K]", ytitle="Tbv (modelled) [K]", $
		subtitle=subtitle, /ynozero, charsize=charsize2, $
		yrange=yrange
oplot, findgen(300), findgen(300)

r=correlate(tbv_meas[ind_sw], tbv_ppe[ind_sw])
subtitle="r="+string(r, format="(f8.3)")
plot, tbv_meas[ind_sw], tbv_ppe[ind_sw], psym=8, title="b. Plane-parallel ensemble", $
		xtitle="Tbv (measured) [K]", ytitle="Tbv (modelled) [K]", $
		subtitle=subtitle, /ynozero, charsize=charsize2, $
		yrange=yrange
oplot, findgen(300), findgen(300)

r=correlate(tbv_meas[ind_sw], tbv[ind_sw])
subtitle="r="+string(r, format="(f8.3)")
plot, tbv_meas[ind_sw], tbv[ind_sw], psym=8, title="c. Ridged Monte Carlo model", $
		xtitle="Tbv (measured) [K]", ytitle="Tbv (modelled) [K]", $
		subtitle=subtitle, /ynozero, charsize=charsize2, $
		yrange=yrange
oplot, findgen(300), findgen(300)

tbhall=[tbh_pp, tbh_ppe, tbh]
yrange=[min(tbh_ppe), max(tbh_pp)]

r=correlate(tbh_meas[ind_sw], tbh_pp[ind_sw])
subtitle="r="+string(r, format="(f8.3)")
plot, tbh_meas[ind_sw], tbh_pp[ind_sw], psym=8, title="d. Plane-parallel model", $
		xtitle="Tbh (measured) [K]", ytitle="Tbh (modelled) [K]", $
		subtitle=subtitle, /ynozero, charsize=charsize2, $
		yrange=yrange
oplot, findgen(300), findgen(300)

r=correlate(tbh_meas[ind_sw], tbh_ppe[ind_sw])
subtitle="r="+string(r, format="(f8.3)")
plot, tbh_meas[ind_sw], tbh_ppe[ind_sw], psym=8, title="e. Plane-parallel ensemble", $
		xtitle="Tbh (measured) [K]", ytitle="Tbh (modelled) [K]", $
		subtitle=subtitle, /ynozero, charsize=charsize2, $
		yrange=yrange
oplot, findgen(300), findgen(300)

r=correlate(tbh_meas[ind_sw], tbh[ind_sw])
subtitle="r="+string(r, format="(f8.3)")
plot, tbh_meas[ind_sw], tbh[ind_sw], psym=8, title="f. Ridged Monte Carlo model", $
		xtitle="Tbh (measured) [K]", ytitle="Tbh (modelled) [K]", $
		subtitle=subtitle, /ynozero, charsize=charsize2, $
		yrange=yrange
oplot, findgen(300), findgen(300)

r=correlate(tbv_meas[ind_sw]-tbh_meas[ind_sw], tbv_pp[ind_sw]-tbh_pp[ind_sw])
subtitle="r="+string(r, format="(f8.3)")
plot, tbv_meas[ind_sw]-tbh_meas[ind_sw], tbv_pp[ind_sw]-tbh_pp[ind_sw], $
		psym=8, title="g. Plane-parallel model", $
		xtitle="Q (measured) [K]", ytitle="Q (modelled) [K]", $
		subtitle=subtitle, charsize=charsize2, yrange=[0, 60]
oplot, findgen(300), findgen(300)

r=correlate(tbv_meas[ind_sw]-tbh_meas[ind_sw], tbv_ppe[ind_sw]-tbh_ppe[ind_sw])
subtitle="r="+string(r, format="(f8.3)")
plot, tbv_meas[ind_sw]-tbh_meas[ind_sw], tbv_ppe[ind_sw]-tbh_ppe[ind_sw], $
		psym=8, title="h. Plane-parallel ensemble", $
		xtitle="Q (measured) [K]", ytitle="Q (modelled) [K]", $
		subtitle=subtitle, charsize=charsize2, yrange=!y.crange
oplot, findgen(300), findgen(300)

r=correlate(tbv_meas[ind_sw]-tbh_meas[ind_sw], tbv[ind_sw]-tbh[ind_sw])
subtitle="r="+string(r, format="(f8.3)")
plot, tbv_meas[ind_sw]-tbh_meas[ind_sw], tbv[ind_sw]-tbh[ind_sw], $
		psym=8, title="i. Ridged Monte Carlo model", $
		xtitle="Q (measured) [K]", ytitle="Q (modelled) [K]", $
		subtitle=subtitle, charsize=charsize2, yrange=!y.crange
oplot, findgen(300), findgen(300)

ps_end

!p.multi=[0, 1, 1]

;I'm interested in the cross-correlation matrices:
a=transpose([[tbv], [tbh], [tbv_meas], [tbh_meas]])
print, "Monte Carlo:"
print, correlate(a)

a=transpose([[tbv_ppe], [tbh_ppe], [tbv_meas], [tbh_meas]])
print, "Ensemble:"
print, correlate(a)

;try it in transformed coordinates:
n=n_elements(tbv_meas)
tbv_ave=total(tbv_meas)/n
tbh_ave=total(tbh_meas)/n
tbv_std=stddev(tbv_meas)
tbh_std=stddev(tbh_meas)
b=[[(tbv_meas-tbv_ave)/tbv_std], [(tbh_meas-tbh_ave)/tbh_std]]
bb=transpose(b)#b
lam=eigenql(bb, eigenvectors=v)

;can probably do it by transforming above, 
;but I'm too lazy to do the math:
a=[transpose([[(tbv-tbv_ave)/tbv_std], [(tbh-tbh_ave)/tbh_std]]#v), transpose(b#v)]
print, "Monte Carlo:"
print, correlate(a)
a=[transpose([[(tbv_ppe-tbv_ave)/tbv_std], [(tbh_ppe-tbh_ave)/tbh_std]]#v), transpose(b#v)]
print, "Ensemble:"
print, correlate(a)

end

