@read_PolIce
@read_thickness

flist=strarr(100)
openr, 1, "flist2a.txt"
on_ioerror, finish
readf, 1, flist


finish:
  close, 1
  flist=flist[where(flist ne "")]

on_ioerror, null

flist="corrected/"+flist

;plot, findgen(10), xrange=[0, 50], yrange=[-0.1, 0.5], /nodata

nmax=2000000

tbs0=fltarr(4, nmax)
tbs40=fltarr(4, nmax)

asi0=fltarr(nmax)
asi40=fltarr(nmax)

skt0=fltarr(nmax)
skt40=fltarr(nmax)

lon0=fltarr(nmax)
lat0=fltarr(nmax)

lon40=fltarr(nmax)
lat40=fltarr(nmax)

zangle0=fltarr(nmax)
zangle40=fltarr(nmax)

track0=fltarr(nmax)
track40=fltarr(nmax)

head0=fltarr(nmax)
head40=fltarr(nmax)

speed0=fltarr(nmax)
speed40=fltarr(nmax)

date0=replicate({time_str}, nmax)
date40=replicate({time_str}, nmax)

fno0=intarr(nmax)
fno40=intarr(nmax)

roll0=fltarr(nmax)
roll40=fltarr(nmax)

pitch0=fltarr(nmax)
pitch40=fltarr(nmax)

j0=0L
j40=0L
for i=0, n_elements(flist)-1 do begin
  print, i

  date=long(strmid(flist[i], 10, 8))
  t0={time_str}
  t0.year=date/10000L
  t0.day=date mod 100
  t0.month=(date-t0.year*10000L)/100

  fname1=flist[i]
  print, fname1
  fname2=flist[i]
  strput, fname1, ".r31", strlen(fname2)-4
 
  read_PolIce, fname1, t0, tbs, lon, lat, t, err=err, $
		angle0=-40, /cor, zangle=zangle, $
		track=track, head=head, speed=speed, $
		roll=roll, pitch=pitch
  if err ne 0 then continue
  nel=n_elements(lon)

  tbs40[*, j40:j40+nel-1]=tbs

  asi40[j40]=get_ice_conc(t0, lon, lat, 1)
  skt40[j40]=get_skin_temp(t[nel/2], lon, lat)

  lon40[j40]=lon
  lat40[j40]=lat
  zangle40[j40]=zangle

  track40[j40]=track
  head40[j40]=head
  speed40[j40]=speed
  roll40[j40]=roll
  pitch40[j40]=pitch

  date40[j40]=t

  fno40[j40:j40+nel-1]=i

  j40=j40+nel
  print, nel, " elements found in ", fname1

  read_PolIce, fname2, t0, tbs, lon, lat, t, err=err, $
		zangle=zangle, angle0=0, /cor, $
		track=track, head=head, speed=speed, $
		roll=roll, pitch=pitch
  if err ne 0 then continue
  nel=n_elements(lon)
  tbs0[*, j0:j0+nel-1]=tbs
  asi0[j0]=get_ice_conc(t0, lon, lat, 1)
  skt0[j0]=get_skin_temp(t[nel/2], lon, lat)

  lon0[j0]=lon
  lat0[j0]=lat
  zangle0[j0]=zangle

  track0[j0]=track
  head0[j0]=head
  speed0[j0]=speed
  roll0[j0]=roll
  pitch0[j0]=pitch

  date0[j0]=t

  fno0[j0:j0+nel-1]=i


  j0=j0+nel



  ;oplot, zangle40, (tbs40[1, *]-tbs40[0, *])/skt, psym=3
;  plot, zangle40, (tbs40[1, *]-tbs40[0, *])/skt, psym=3, $
;		xrange=[0, 50], yrange=[0, 0.5]
;  oplot, zangle0, (tbs0[0, *]-tbs0[1, *])/skt, psym=3

;  stop

endfor

tbs0=tbs0[*, 0:j0-1]
tbs40=tbs40[*, 0:j40-1]
asi0=asi0[0:j0-1]
asi40=asi40[0:j40-1]
skt0=skt0[0:j0-1]
skt40=skt40[0:j40-1]

lon0=lon0[0:j0-1]
lat0=lat0[0:j0-1]
lon40=lon40[0:j40-1]
lat40=lat40[0:j40-1]
zangle0=zangle0[0:j0-1]
zangle40=zangle40[0:j40-1]

speed0=speed0[0:j0-1]
track0=track0[0:j0-1]
head0=head0[0:j0-1]
roll0=roll0[0:j0-1]
pitch0=pitch0[0:j0-1]

speed40=speed40[0:j40-1]
track40=track40[0:j40-1]
head40=head40[0:j40-1]
roll40=roll40[0:j40-1]
pitch40=pitch40[0:j40-1]

date0=date0[0:j0-1]
date40=date40[0:j40-1]

fno0=fno0[0:j0-1]
fno40=fno40[0:j40-1]

em0=tbs0/rebin(transpose(skt0), 4, j0)
em40=tbs40/rebin(transpose(skt40), 4, j40)

save, tbs0, tbs40, asi0, asi40, skt0, skt40, date0, date40, $
		lon0, lon40, lat0, lat40, zangle0, zangle40, $
		speed0, speed40, track0, track40, head0, head40, $
		filename="smos_v_asi.idlsav", fno0, fno40, flist, $
		roll0, roll40, pitch0, pitch40

end


