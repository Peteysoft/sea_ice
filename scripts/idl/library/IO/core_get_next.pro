;reads the next core sample from an ASCII file of compiled
;samples from H. Eicken

;z = vertical levels for temperature
;z1, z2 = vertical bounds for salinity

pro core_get_next, lun, date, lon, lat, z, t, z1, z2, s, err=err

  on_ioerror, done

  id=""
  readf, lun, id
  print, id
  date=""
  readf, lun, date
  print, date
  comment=""
  readf, lun, comment
  lon=0. 
  lat=0.
  readf, lun, lon, lat
  print, "(", lon, ",", lat, ")"
  header=strarr(5)
  readf, lun, header
  ;print, header
  n1=0L
  readf, lun, n1
  print, n1, " temperature levels found"
  if n1 gt 0 then begin
    data1=fltarr(2, n1)
    readf, lun, data1
    z=reform(data1[0, *])
    t=reform(data1[1, *])
  endif
  n2=0L
  readf, lun, n2
  print, n2, " salinity levels found"
  data2=fltarr(4, n2)
  readf, lun, data2
  z1=reform(data2[0, *])
  z2=reform(data2[1, *])
  ind=where(z2-z1 lt 0, cnt)
  if cnt gt 0 then stop
  s=reform(data2[3, *])

  n3=0L
  readf, lun, n3
  print, n3, " brine volume levels found"
  if n3 gt 0 then begin
    stuff=strarr(n3)
    readf, lun, stuff
  endif
  junk=""
  readf, lun, junk

  on_ioerror, null
  err=0

  return

  done:
    print, !error_state.msg
    ;stop
    err=1
end


