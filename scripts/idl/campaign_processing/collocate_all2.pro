@time_lib
@boundary_geometry
@vector_util

pro coll, rname, tname1, tname2=tname2, tbs, dave, dstd, n, tr, date=date, $
		angle0=angle0, sigma=sigma, zangle=zangle, pitch=pitch, roll=roll, $
		alpha=alpha, lon=lonr, lat=latr, footprint=footprint, $
		badness=badness
  if n_elements(sigma) ne 1 then sigma=0.23
  if n_elements(footprint) ne 1 then footprint=400.
  if n_elements(date) ne 1 then date='2007/03/13'
  if n_elements(angle0) ne 1 then angle0=-40

  print, "Reading radiance measurements"
  read_police, rname, date, tbs, lonr, latr, tr, /cor, $
                altitude=altitude, angle0=angle0, zangle=zangle, $
		pitch=pitch, roll=roll, lon0=lon0, lat0=lat0
  rotate_flight_track, lon0, lat0, x0, y0, /new
  rotate_flight_track, lonr, latr, xr, yr

  print, "Reading ice thickness measurements"
  read_thickness, tname1, thick, lont, latt, tt1, flag
  rotate_flight_track, lont, latt, xt, yt
  if n_elements(tname2) eq 1 then begin
    read_thickness, tname2, thick2, lont2, latt2, tt2, flag2
    rotate_flight_track, lont2, latt2, xt2, yt2
    round_sym, 0.1
    ;plot, xr, yr, psym=3
    ;oplot, xt, yt, psym=8
    ;oplot, xt2, yt2, psym=8
    ;stop

    xt=[xt, xt2]
    yt=[yt, yt2]
    thick=[thick, thick2]
    flag=[flag, flag2]

  endif

  nr=n_elements(lonr)

  dave=fltarr(nr)
  dstd=fltarr(nr)
  n=fltarr(nr)
  badness=fltarr(nr)
  if arg_present(alpha) then alpha=fltarr(nr)

  ;sort the thickness data by x coord:
  sind=sort(xt)
  xt=xt[sind]
  yt=yt[sind]
  thick=thick[sind]
  flag=flag[sind]

  nt=n_elements(xt)

  for i=0L, nr-1 do begin
    if i mod 100 eq 0 then print, 100.*i/(nr-1)
    ;ind=where(abs(xt-xr[i]) lt footprint, cnt)
    ;if cnt eq 0 then begin
    ;  n[i]=0
    ;  continue
    ;endif
    ;thick1=thick[ind]

    dx1=xr[i]-x0[i]
    dy1=yr[i]-y0[i]
    dx2=xt-x0[i]
    dy2=yt-y0[i]
    z2=altitude[i]*altitude[i]

    theta=acos((dx1*dx2+dy1*dy2+z2)/sqrt(dx1*dx1+dy1*dy1+z2)/sqrt(dx2*dx2+dy2*dy2+z2))
    w=exp(-theta^2/sigma^2/2)
    ;print, 100.*i/(nr-1), min(w), min(abs(yr[i]-yt[ind]))
    tw=total(w)
    n[i]=tw

    dave[i]=total(thick*w)/tw
    dstd[i]=sqrt(total(((thick-dave[i])*w)^2)/tw)

    badness[i]=total(flag*w)/tw

    ;there are a ton of redundant calculations here;
    ;one of these days I'm going to unravel this tangled mess...
    ;one of these days...
    if arg_present(alpha) then begin
      ;angle between polarisation planes defined by the ground, and
      ;those defined by the instrument (assumes a flat ground):
      roll1=-roll[i]*!dtor
      pitch1=pitch[i]*!dtor
      vangle=angle0*!dtor

      ;normal vector:
      cosroll=cos(roll1)
      norm=[sin(pitch1)*cosroll, sin(roll1), cos(pitch1)*cosroll]
      v=[sin(vangle), 0, cos(vangle)]
      ncrossv=cross(norm, v)
      alpha1=acos(ncrossv[1]/sqrt(dot(ncrossv, ncrossv)))
      alpha[i]=alpha1*abs(norm[1])/norm[1]
      ;stop
    endif
  endfor

end  

;faster and just as good as the previous, "definitive" versions

;goto, save1

vtype=""

if vtype eq "nadir" then suffix=".r31" else suffix=".r32"

rname="corrected/"+[$
	'20070313_TrackAndCircles/P4X_to_P2A/07217460', $
	'20070313_TrackAndCircles/Kruunupyy_to_P3X/07216200.part', $
	'20070313_TrackAndCircles/P2A_to_P4X/07217260', $
	'20070313_TrackAndCircles/P3X_to_P2A/07217010', $
	'20070312_TrackAndCircles/Marjaniemi_to_P7X/07117410', $
	'20070312_TrackAndCircles/P7X_to_Marjaniemi/07116520', $
	'20070312_TrackAndCircles/P2X_to_Marjaniemi/07117250', $
	'20070312_TrackAndCircles/Marjaniemi_to_P2X/07117100']
rname=rname+suffix

date='2007/03/13'
tname="ice_thickness/20070313_P2A_P4X.dat"
coll, rname[0], tname, tbs1, dave1, dstd1, n1, t1, $
		date=date, badness=bad1, $
		zangle=zangle1, alpha=alpha1, $
		lon=lon1, lat=lat1

tname="ice_thickness/20070313_P3X_P2A.dat"
coll, rname[1], tname, tbs2, dave2, dstd2, n2, t2, $
		date=date, badness=bad2, $
		zangle=zangle2, alpha=alpha2, $
		lon=lon2, lat=lat2

tname="ice_thickness/20070313_P2A_P4X.dat"
coll, rname[2], tname, tbs3, dave3, dstd3, n3, t3, $
		date=date, badness=bad3, $
		zangle=zangle3, alpha=alpha3, $
		lon=lon3, lat=lat3

tname="ice_thickness/20070313_P3X_P2A.dat"
coll, rname[3], tname, tbs4, dave4, dstd4, n4, t4, $
		date=date, badness=bad4, $
		zangle=zangle4, alpha=alpha4, $
		lon=lon4, lat=lat4

date='2007/03/12'
tname1="ice_thickness/20070312_P7X_MRJN_afternoon.dat"
tname2="ice_thickness/20070312_P7X_MRJN_morning.dat"
coll, rname[4], tname1, tname2=tname2, tbs5, dave5, dstd5, n5, t5, $
		date=date, badness=bad5, $
		zangle=zangle5, alpha=alpha5, $
		lon=lon5, lat=lat5

tname1="ice_thickness/20070312_P7X_MRJN_afternoon.dat"
tname2="ice_thickness/20070312_P7X_MRJN_morning.dat"
coll, rname[5], tname1, tname2=tname2, tbs6, dave6, dstd6, n6, t6, $
		date=date, badness=bad6, $
		zangle=zangle6, alpha=alpha6, $
		lon=lon6, lat=lat6

tname1="ice_thickness/20070312_MRJN_P2X_afternoon.dat"
tname2="ice_thickness/20070312_MRJN_P2X_morning.dat"
coll, rname[6], tname1, tname2=tname2, tbs7, dave7, dstd7, n7, t7, $
		date=date, badness=bad7, $
		zangle=zangle7, alpha=alpha7, $
		lon=lon7, lat=lat7

tname1="ice_thickness/20070312_MRJN_P2X_afternoon.dat"
tname2="ice_thickness/20070312_MRJN_P2X_morning.dat"
coll, rname[7], tname1, tname2=tname2, tbs8, dave8, dstd8, n8, t8, $
		date=date, badness=bad8, $
		zangle=zangle8, alpha=alpha8, $
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
dstd=[dstd1, dstd2, dstd3, dstd4, dstd5, dstd6, dstd7, dstd8]
nsamp=[n1, n2, n3, n4, n5, n6, n7, n8]
zangle=[zangle1, zangle2, zangle3, zangle4, zangle5, zangle6, zangle7, zangle8]
alpha=[alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8]
lon=[lon1, lon2, lon3, lon4, lon5, lon6, lon7, lon8]
lat=[lat1, lat2, lat3, lat4, lat5, lat6, lat7, lat8]
t=[t1, t2, t3, t4, t5, t6, t7, t8]
badness=[bad1, bad2, bad3, bad4, bad5, bad6, bad7, bad8]

goto, skip
ind1=where(nsamp gt 0)
nsamp=nsamp[ind1]
tbs=tbs[*, ind1]
dave=dave[ind1]
dstd=dstd[ind1]
zangle=zangle[ind1]
alpha=alpha[ind1]
flight_no=flight_no[ind1]
pindex=pindex[ind1]
lon=lon[ind1]
lat=lat[ind1]
skip:

ind2=where(nsamp eq 1, cnt)
if cnt ne 0 then dstd[ind2]=-1

save, filename="collall2."+vtype+".idlsav", nsamp, rname, $
		tbs, dave, dstd, badness, $
		zangle, alpha, zangle, alpha, $
		lon, lat, flight_no, pindex, t

end


