;rotates a roughly straight flight track so that it now 
;consists of a distance, x, in metres along the track and
;another distance, y, in metres perpendicular

;uses the angle of the first converted track as the basis
;to convert subsequent tracks

pro rotate_flight_track, lon, lat, x, y, new=new
  common rft_base, lon0, lat0, sinangle, cosangle

  mperdeg=1000000./9.

  if n_elements(lon0) eq 0 then new=1

  n=n_elements(lon)

  coslat=cos(!dtor*lat)

  if keyword_set(new) then begin
    lon0=lon[0]
    lat0=lat[0]
    dlat=lat[n-1]-lat[0]
    coslat2=cos(!dtor*(lat[0]+lat[n-1])/2)
    dx=(lon[n-1]-lon[0])*coslat2
    r=sqrt(dlat*dlat+dx*dx)
    sinangle=dlat/r
    cosangle=dx/r
    xoffset=0.
    yoffset=0.
  endif else begin
    ;xoffset0=(lon[0]-lon0)*coslat
    ;yoffset0=(lat[0]-lat0)
    ;xoffset=(xoffset0*cosangle+yoffset0*sinangle)*mperdeg
    ;yoffset=(-xoffset0*sinangle+yoffset0*cosangle)*mperdeg
  endelse

  x0=(lon-lon0)*coslat
  y0=lat-lat0

  x=(x0*cosangle+y0*sinangle)*mperdeg;+xoffset
  y=(-x0*sinangle+y0*cosangle)*mperdeg;+yoffset

  ;stop

end

;the inverse transformation:
pro unrotate_flight_track, x, y, lon, lat
  common rft_base, lon0, lat0, sinangle, cosangle

  mperdeg=1000000./9.

  x1=x*cosangle-y*sinangle
  y1=x*sinangle+y*cosangle

  lat=y1/mperdeg+lat0
  coslat=cos(lat*!dtor)
  lon=x1/mperdeg/coslat+lon0

end

