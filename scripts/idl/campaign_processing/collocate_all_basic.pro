@time_lib
@vector_util

;simplest collocation method, used for illustrative purposes only...

pro coll, rname, tname1, tname2=tname2, tbs, d, nsamp, tr, date=date, $
		angle0=angle0, zangle=zangle, $
		lon=lonr, lat=latr, flag=flag
  if n_elements(sigma) ne 1 then sigma=0.23
  if n_elements(date) ne 1 then date='2007/03/13'
  if n_elements(angle0) ne 1 then angle0=-40

  print, "Reading radiance measurements"
  read_police, rname, date, tbs, lonr, latr, tr, /cor, $
                altitude=altitude, angle0=angle0, zangle=zangle
  rotate_flight_track, lonr, latr, xr, yr, /new
  nr=n_elements(lonr)

  print, "Reading ice thickness measurements"
  read_thickness, tname1, thick, lont, latt, tt1, flag
  rotate_flight_track, lont, latt, xt, yt
  nt=n_elements(xt)
  if xt[0] gt xt[nt-1] then begin
    xt=reverse(xt)
    yt=reverse(yt)
    thick=reverse(thick)
  endif

  if n_elements(tname2) eq 1 then begin
    read_thickness, tname2, thick2, lont2, latt2, tt2, flag
    rotate_flight_track, lont2, latt2, xt2, yt2
    round_sym, 0.1
    nt2=n_elements(xt2)
    ;plot, xr, yr, psym=3
    ;oplot, xt, yt, psym=8
    ;oplot, xt2, yt2, psym=8
    ;stop
    if xt2[0] gt xt2[nt2-1] then begin
      xt2=reverse(xt2)
      yt2=reverse(yt2)
      thick2=reverse(thick2)
    endif
  endif

  d=fltarr(nr)-1
  if arg_present(alpha) then alpha=fltarr(nr)

  nsamp=lonarr(nr)

  for i=0L, nr-1 do begin
    if i mod 100 eq 0 then print, 100.*i/(nr-1)

    if xr[i] ge xt[0] or xr[i] le xt[nt-1] then begin
      d[i]=interpol(thick, xt, xr[i])
      nsamp[i]=2
    endif

    if n_elements(tname2) eq 1 then begin
      if xr[i] ge xt2[0] or xr[i] le xt[nt-1] then begin
        if d[i] eq -1 then begin
          d[i]=interpol(thick2, xt2, xr[i])
          nsamp[i]=2
        endif else begin
          d2=interpol(thick2, xt2, xr[i])
          dy1=abs(interpol(yt, xt, xr[i])-yr[i])
          dy2=abs(interpol(yt2, xt2, xr[i])-yr[i])
          d[i]=d2+dy1*(d[i]-d2)/(dy1+dy2)
          nsamp[i]=4
        endelse
      endif
    endif
  endfor

end

;faster and just as good as the previous, "definitive" versions

;goto, save1

date='2007/03/13'
fov=7.5

rname="corrected/"+[$
	'20070313_TrackAndCircles/P4X_to_P2A/07217460.r31', $
	'20070313_TrackAndCircles/Kruunupyy_to_P3X/07216200.part.r31', $
	'20070313_TrackAndCircles/P2A_to_P4X/07217260.r31', $
	'20070313_TrackAndCircles/P3X_to_P2A/07217010.r31', $
	'20070312_TrackAndCircles/Marjaniemi_to_P7X/07117410.r31', $
	'20070312_TrackAndCircles/P7X_to_Marjaniemi/07116520.r31', $
	'20070312_TrackAndCircles/P2X_to_Marjaniemi/07117250.r31', $
	'20070312_TrackAndCircles/Marjaniemi_to_P2X/07117100.r31']

tname="ice_thickness/20070313_P2A_P4X.dat"
coll, rname[0], tname, tbs1, dave1, n1, t1, $
		date=date, $
		zangle=zangle1, $
		lon=lon1, lat=lat1

tname="ice_thickness/20070313_P3X_P2A.dat"
coll, rname[1], tname, tbs2, dave2, n2, t2, $
		date=date, $
		zangle=zangle2, $
		lon=lon2, lat=lat2

tname="ice_thickness/20070313_P2A_P4X.dat"
coll, rname[2], tname, tbs3, dave3, n3, t3, $
		date=date, $
		zangle=zangle3, $
		lon=lon3, lat=lat3

tname="ice_thickness/20070313_P3X_P2A.dat"
coll, rname[3], tname, tbs4, dave4, n4, t4, $
		date=date, $
		zangle=zangle4, $
		lon=lon4, lat=lat4

date='2007/03/12'
tname1="ice_thickness/20070312_P7X_MRJN_afternoon.dat"
tname2="ice_thickness/20070312_P7X_MRJN_morning.dat"
coll, rname[4], tname1, tname2=tname2, tbs5, dave5, n5, t5, $
		date=date, $
		zangle=zangle5, $
		lon=lon5, lat=lat5

tname1="ice_thickness/20070312_P7X_MRJN_afternoon.dat"
tname2="ice_thickness/20070312_P7X_MRJN_morning.dat"
coll, rname[5], tname1, tname2=tname2, tbs6, dave6, n6, t6, $
		date=date, $
		zangle=zangle6, $
		lon=lon6, lat=lat6

tname1="ice_thickness/20070312_MRJN_P2X_afternoon.dat"
tname2="ice_thickness/20070312_MRJN_P2X_morning.dat"
coll, rname[6], tname1, tname2=tname2, tbs7, dave7, n7, t7, $
		date=date, $
		zangle=zangle7, $
		lon=lon7, lat=lat7

tname1="ice_thickness/20070312_MRJN_P2X_afternoon.dat"
tname2="ice_thickness/20070312_MRJN_P2X_morning.dat"
coll, rname[7], tname1, tname2=tname2, tbs8, dave8, n8, t8, $
		date=date, $
		zangle=zangle8, $
		lon=lon8, lat=lat8

save1:

flight_no=[intarr(n_elements(n1)), $
		intarr(n_elements(n2))+1, $
		intarr(n_elements(n3))+2, $
		intarr(n_elements(n4))+3, $
		intarr(n_elements(n5))+4, $
		intarr(n_elements(n6))+5, $
		intarr(n_elements(n7))+6, $
		intarr(n_elements(n8))+7]

pindex=[lindgen(n_elements(n1)), $
		lindgen(n_elements(n2)), $
		lindgen(n_elements(n3)), $
		lindgen(n_elements(n4)), $
		lindgen(n_elements(n5)), $
		lindgen(n_elements(n6)), $
		lindgen(n_elements(n7)), $
		lindgen(n_elements(n8))]

tbs=[[tbs1], [tbs2], [tbs3], [tbs4], [tbs5], [tbs6], [tbs7], [tbs8]]
dave=[dave1, dave2, dave3, dave4, dave5, dave6, dave7, dave8]
nsamp=[n1, n2, n3, n4, n5, n6, n7, n8]
zangle=[zangle1, zangle2, zangle3, zangle4, zangle5, zangle6, zangle7, zangle8]
lon=[lon1, lon2, lon3, lon4, lon5, lon6, lon7, lon8]
lat=[lat1, lat2, lat3, lat4, lat5, lat6, lat7, lat8]
t=[t1, t2, t3, t4, t5, t6, t7, t8]

goto, skip
ind1=where(nsamp gt 0)
nsamp=nsamp[ind1]
tbs=tbs[*, ind1]
dave=dave[ind1]
zangle=zangle[ind1]
flight_no=flight_no[ind1]
pindex=pindex[ind1]
lon=lon[ind1]
lat=lat[ind1]
skip:

ind2=where(nsamp eq 1, cnt)
if cnt ne 0 then dstd[ind2]=-1

save, filename="collall_basic.idlsav", rname, nsamp, $
		tbs, dave, zangle, lon, lat, flight_no, pindex, t

end


