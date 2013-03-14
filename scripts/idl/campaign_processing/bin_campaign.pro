@time_lib
@boundary_geometry

;average campaign data into bins using the natural groupings of the
;measurements:

;goto, plt

restore, "collall2.idlsav"

dave=double(dave)
dstd=double(dstd)

;ds=point_spacing(lon, lat, /nowrap)
t0=convert_date('2007/3/12')
tind=day_value(t, t0)*86400.

nall=n_elements(dave)
;ind=[-1, where(ds gt 2., cnt), nall-1]
dt=tind[1:nall-1]-tind[0:nall-2]
ind=[-1, where(abs(dt) gt 0.1, cnt), n_elements(dave)-1]

dbin=dblarr(cnt+1)
stdbin=fltarr(cnt+1)
lon_bin=fltarr(cnt+1)
lat_bin=fltarr(cnt+1)
tbs_bin=fltarr(4, cnt+1)
tbs_std=fltarr(4, cnt+1)

ts_bin=fltarr(cnt+1)		;total number of thickness samples
nr_bin=lonarr(cnt+1)		;total number of radiometer samples

za_bin=fltarr(cnt+1)
za_std=fltarr(cnt+1)

fno_bin=intarr(cnt+1)

badbin=fltarr(cnt+1)

pindex1=lonarr(cnt+1)
pindex2=lonarr(cnt+1)

tbin=replicate({time_str}, cnt+1)

for i=0L, cnt do begin
  l1=ind[i]+1
  l2=ind[i+1]

  ;we want to weight the averages by the number of samples in each:
  ns=nsamp[l1:l2]
  ts=total(ns)
  print, ts

  dave1=dave[l1:l2]
  dstd1=dstd[l1:l2]
  dbin[i]=total(dave1*ns)/ts
  ind1=where(dstd1 lt 0, cnt)
  if cnt ne 0 then dstd1[ind1]=0
  stdbin[i]=(total((dstd1^2+dave1^2)*(ns-1))-dbin[i]^2)/(ts-1)

  ;this algorithm should really be general;
  ;I shouldn't be messing around with special cases...
  ;one more reason to stop using IDL...
  if l2-l1 eq 0 then begin
    tbs_bin[*, i]=tbs[*, l1]
    tbs_std[*, i]=0
  endif else begin
    tb1=tbs[*, l1:l2]
    ns1=rebin(transpose(ns), 4, l2-l1+1)
    ;weighted average:
    if ts gt 0 then begin
      tbs_bin[*, i]=total(tb1*ns1, 2)/ts
      ;weighted standard-deviation:
      tbs_std[*, i]=sqrt(total(((tb1-rebin(tbs_bin[*, i], 4, l2-l1+1))*ns1)^2, 2)/(ts-1))
    endif else begin
      tbs_bin[*, i]=total(tb1)/(l2-l1+1)
      tbs_std[*, i]=sqrt(total((tb1-rebin(tbs_bin[*, i], 4, l2-l1+1))^2)/(l2-l1))
    endelse
  endelse

  lon_bin[i]=total(lon[l1:l2]*ns)/ts
  lat_bin[i]=total(lat[l1:l2]*ns)/ts

  badbin[i]=total(badness[l1:l2]*ns)/ts

  za_bin[i]=total(zangle[l1:l2]*ns)/ts
  za_std[i]=sqrt(total(((zangle[l1:l2]-za_bin[i])*ns)^2)/(ts-1))

  ts_bin[i]=ts
  nr_bin[i]=l2-l1+1

  ;they should all be from the same flight:
  check=where(flight_no[l1] ne flight_no[l1:l2], cnt)
  if cnt ne 0 then stop

  fno_bin[i]=flight_no[l1]

  pindex1[i]=pindex[l1]
  pindex2[i]=pindex[l2]

  tbin[i]=t[(l1+l2)/2]

endfor

save, filename="collave2.idlsav", dbin, stdbin, tbs_bin, tbs_std, lon_bin, lat_bin, $
		za_bin, za_std, ts_bin, nr_bin, fno_bin, $
		rname, pindex1, pindex2, badbin

plt:

;plot some of the data:
round_sym, 0.5
ps_start, "averaged.ps", /color
loadct, 13, ncol=9, bottom=1
tvlct, r, g, b, /get

;quality check:
ind=where(za_std lt 1 and za_bin lt 38 and za_bin gt 34 and tbs_std[0, *] lt 10 $
		and tbs_std[1, *] lt 10)

plot, tbs_bin[0, ind], tbs_bin[1, ind], psym=8, xtitle="Tbh [K]", ytitle="Tbv [K]", $
		/ynozero, /nodata
plots, tbs_bin[0, ind], tbs_bin[1, ind], psym=8, color=fno_bin[ind]+2

draw_legend, [108, 111], [270, 245], rname, psym=8, color=indgen(8)+2, charsize=0.5

plot, tbs_bin[1, ind]-tbs_bin[0, ind], tbs_bin[1, ind], psym=8, $
		xtitle="Q [K]", ytitle="Tbv [K]", /ynozero, /nodata
draw_legend, [4, 5], [190, 165], rname, psym=8, color=indgen(8)+2, charsize=0.5
plots, tbs_bin[1, ind]-tbs_bin[0, ind], tbs_bin[1, ind], psym=8, color=fno_bin[ind]+2

plot, dbin[ind], tbs_bin[1, ind], psym=8, xtitle="ice thickness [m]", ytitle="Tbv [K]", $
		/ynozero, /nodata
plots, dbin[ind], tbs_bin[1, ind], psym=8, color=fno_bin[ind]+2

draw_legend, [1.7, 1.75], [195, 165], rname, psym=8, color=indgen(8)+2, charsize=0.5

plot, dbin[ind], tbs_bin[1, ind]-tbs_bin[0, ind], psym=8, xtitle="ice thickness [m]", $
		ytitle="Q [K]", /ynozero, /nodata
plots, dbin[ind], tbs_bin[1, ind]-tbs_bin[0, ind], psym=8, color=fno_bin[ind]+2

draw_legend, [1.8, 1.85], [57, 46], rname, psym=8, color=indgen(8)+2, charsize=0.5

ps_end

end

