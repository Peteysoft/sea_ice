@time_lib

pro read_eicken_prof, filename, type, data, $
		date=date, lon=lon, lat=lat, $
		shipno=shipno, tair=tair, tw=tw, $
		snow_depth=snow_depth, tsnow=tsnow

  maxn=10
  data=ptrarr(maxn)
  type=lonarr(maxn)

  on_ioerror, finish
  openr, lun, filename, /get_lun
  
  ;read the headers:
  header=""
  readf, lun, header
  ;print, header
  ;just lop off the date field so it doesn't get in our way
  ;add a dummy field so we don't have to updated the rest of the code...
  field=[strmid(header, 0, 16), strsplit(strmid(header, 16), /extract)]
  nf=n_elements(field)
  datestr=strtrim(field[0], 2)

  date={time_str, strmid(datestr, 6, 4), strmid(datestr, 3, 2), strmid(datestr, 0, 2), $
		strmid(datestr, 10, 2), strmid(datestr, 13, 2), 0}
  if date.hour eq 99 then date.hour=0
  if date.minute eq 99 then date.minute=0

  lon=float(field[1])
  lonlon=long(lon)
  lon=lonlon+(lon-lonlon)*100./60.

  lat=float(field[2])
  lonlat=long(lat)
  lat=lonlat+(lat-lonlat)*100./60.

  shipno=field[3]
  if nf gt 4 then begin
    tair=float(field[4])
    if nf gt 6 then begin
      snow_depth=float(field[6])
      if nf gt 7 then begin
        tsnow=float(field[7])
        if nf gt 8 then begin
          tw=float(field[8])
        endif
      endif else begin
        tsnow=99.
      endelse
    endif
  endif

  ;print, time_string(date), " (", lon, ", ", lat, ")"

  ;read the comments:
  ncom=0L
  readf, lun, ncom
  for i=0, ncom-1 do begin
    readf, lun, header
    ;print, header
  endfor

  ltype=0L
  dtype=0L
  n=0L
  for i=0L, maxn-1 do begin
    readf, lun, ltype, dtype, n
    ;print, ltype, dtype, n
    if n eq 0 then continue
    type[i]=dtype
    if ltype eq 2 then begin
      data1=fltarr(3, n)
      readf, lun, data1
    endif else begin
      data1=fltarr(2, n)
      readf, lun, data1
    endelse
    data[i]=ptr_new(data1)
  endfor

  finish: free_lun, lun

end

pro prepare_eicken_prof, filename, dz, t, type, s, err=err, $
		lon=lon, lat=lat, date=date
  read_eicken_prof, filename, type, data, $
		date=date, lon=lon, lat=lat, $
		shipno=shipno, tair=tair, tw=tw, $
		snow_depth=snow_depth, tsnow=tsnow

  err=1
  indg=where(ptr_valid(data), cnt)
  if cnt eq 0 then return

  ;use salinity measurements as base:
  err=2
  ind=where(type eq 2, cnt)
  if cnt eq 0 then return
  z1=reform((*data[ind[0]])[0, *])
  z2=reform((*data[ind[0]])[1, *])
  s=reform((*data[ind[0]])[2, *])

  err=3
  ind=where(type eq 1, cnt)
  if cnt eq 0 then return
  zt=reform((*data[ind[0]])[0, *])
  t1=reform((*data[ind[0]])[1, *])

  ;interpolate temperatures to salinity profiles:
  zmid=(z1+z2)/2
  t=interpol(t1, zt, zmid)

  dz=z2-z1
  type=intarr(n_elements(dz))+3

  if snow_depth ge 0 and snow_depth lt 99 then begin
    if tsnow ge 99 then begin
      if tair ge 99 then begin
        tsnow=(tair+t[0])/2
      endif else begin
        tsnow=t[0]
      endelse
    endif
    t=[tsnow, t]
    dz=[snow_depth, dz]
    s=[0, s]
    type=[1, type]
  endif

  t=273.2+t

  err=0

end
    
