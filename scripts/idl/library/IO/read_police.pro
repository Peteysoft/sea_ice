@time_lib

;file format:
;column	quantity
;1	Measurement time stamp, UTC, [sec. since midnight]
;2	Channel1 (horizontal) data output [Kelvin]
;3	Channel2 (vertical) data output [Kelvin]
;4	Correlation coefficient, between channel1 and channel2, Real part.
;5	Correlation coefficient, between channel1 and channel2, Imaginary part.
;6	Number of basic (8 msec.) integration periods integrated.
;7	Time stamp, when latest sampling of internal temperatures was performed, UTC.
;8	Internal temperature 1 [ºC]: Vertical load
;9	Internal temperature 2 [ºC]: Vertical LNA
;10	Internal temperature 3 [ºC]: Vertical AMP2
;11	Internal temperature 4 [ºC]: Vertical AMP3
;12	Internal temperature 5 [ºC]: Noise diode
;13	Internal temperature 6 [ºC]: Horizontal load
;14	Internal temperature 7 [ºC]: Horizontal LNA
;15	Internal temperature 8 [ºC]: Horizontal AMP2
;16	Internal temperature 9 [ºC]: Horizontal AMP3
;17	Internal temperature 10 [ºC]: Baseplate
;18	Internal temperature 11 [ºC]: Vertical DFE ADC
;19	Internal temperature 12 [ºC]: Horizontal DFE ADC
;20	Internal temperature 13 [ºC]  Aft-looking antenna horizontal channel antenna cable
;21	Internal temperature 14 [ºC]: Aft-looking antenna vertical channel antenna cable
;22	Internal temperature 15 [ºC]: Nadir antenna horizontal channel antenna cable
;23	Internal temperature 16 [ºC]: Nadir antenna vertical channel antenna cable
;24	Time stamp, when latest sampling of internal temperatures was performed, UTC.
;25	Aircraft position latitude [degrees]
;26	Aircraft position longitude [degrees]
;27	Aircraft altitude [m]
;28	Aircraft track, relative to Earth North [degrees]
;29	Aircraft ground speed [m/sec.]
;30	Aircraft true heading, relative to Earth North [degrees]
;31	Platform azimuth, relative to Earth North [degrees]
;32	Aircraft roll [degrees, positive numbers correspond to right turn]
;33	Aircraft pitch [degrees, positive numbers correspond to nose up]

pro read_PolIce, rfile, t00, tbs, lon, lat, t, err=err, speed=speed, $
		head=head, azimuth=az, roll=roll, pitch=pitch, $
		angle0=angle0, zangle=zangle, altitude=alt, $
		cor=cor, offset=offset, track=track, $
		instr_plane=instrangle, lon0=lon0, lat0=lat0

  maxn=100000
  mperdeg=1000000./9.

  t0=convert_date(t00)

  on_ioerror, finish2
  data2=dblarr(33, maxn)
  openr, lun, rfile, /get_lun
  readf, lun, data2

  finish2: ind=where(data2[0, *] ne 0, n2)
  err=1
  if n2 eq 0 then begin
    if n_elements(lun) ne 0 then free_lun, lun
    return
  endif
  err=0
  data2=data2[*, ind]
  free_lun, lun

  lon=reform(data2[25, *])
  lat=reform(data2[24, *])

  tbs=reform(data2[1:4, *])

  t=replicate(t0, n2)
  ts=reform(data2[0, *])
  t.second=ts mod 60
  t.hour=ts/3600
  t.minute=ts/60-t.hour*60

  head=reform(data2[29, *])
  az=reform(data2[30, *])

  speed=reform(data2[28, *])
  roll=reform(data2[31, *])
  pitch=reform(data2[32, *])
  alt=reform(data2[26, *])

  track=reform(data2[27, *])

  if keyword_set(cor) then begin
    pitch1=(angle0+pitch)*!dtor
    roll1=roll*!dtor
    ;correction along the direction of flight:
    dx=alt*tan(pitch1)
    ;correction perpendicular to flight line:
    dy=alt*sin(roll1)
    ;correction to position:
    head1=head*!dtor
    lon0=lon
    lat0=lat
    lon=lon+(dx*sin(head1)+dy*cos(head1))/mperdeg/cos(lat*!dtor)
    lat=lat+(dx*cos(head1)-dy*sin(head1))/mperdeg
    offset=sqrt(dx*dx+dy*dy)
    ;stop
  endif

  if arg_present(zangle) then begin
    if n_elements(angle0) eq 0 then angle0=0.
    pitch1=pitch*!dtor
    cosroll=cos(roll*!dtor)
    za1=-angle0*!dtor
    zangle=acos(sin(za1)*sin(pitch1)*cosroll + cos(za1)*cos(pitch1)*cosroll)
    zangle=zangle/!dtor
    ind=where(zangle gt 90, cnt)
    if cnt gt 0 then zangle[ind]=180-zangle[ind]
  endif

  on_ioerror, null

  ;stop

end

