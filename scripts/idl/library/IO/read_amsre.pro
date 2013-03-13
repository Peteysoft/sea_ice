FUNCTION _split_coefficients, coeff
RETURN, FLOAT(STRSPLIT(coeff,'([ABG]+\-|,)',/EXTRACT,/REGEX))
END ;_split_coefficients

FUNCTION _split_nl_coefficients, nl_coeff
RETURN, float(STRSPLIT(nl_coeff,'[0-9]{1,2}G[AB]?[VH]-',/EXTRACT,/regex))
END ;_split_coefficients

;simple read routine for AMSR-E data...

;channel numbers:
;1.	6 GHz v
;2.	6 GHz h
;3.	10 GHz v
;4.	10 GHz h
;5.	18 GHz v
;6.	18 GHz h
;7.	23 GHz v
;8.	23 GHz h
;9.	36 GHz v
;10.	36 GHz h
;11.	50 GHz v
;12.	53 GHz v
;13.	89 (a) GHz v
;14.	89 (a) GHz h
;(currently the lons and lats don't match for these...)
;15.	89 (b) GHz v
;16.	89 (b) GHz h
pro read_amsre, filename, chan, tbv, tbh, lon=lon, lat=lat, $
		date=date, time=time, err=err

err=0

catch, err
if err ne 0 then begin
  message, !error_state.msg, /continue
  catch, /cancel
  return
end

id=hdf_sd_start(filename)

hdf_sd_attrinfo, id, 7, data=time
hdf_sd_attrinfo, id, 8, data=date
print, date, "   ", time

;get coefficients:
if arg_present(tbv) then begin
  sid=hdf_sd_select(id, 21)
  hdf_sd_getdata, sid, coeff

  sid=hdf_sd_select(id, chan)
  hdf_sd_getinfo, sid, name=name
  print, name
  hdf_sd_attrinfo, sid, 0, data=sf
  hdf_sd_attrinfo, sid, 1, data=unit
  print, unit

  hdf_sd_getdata, sid, counts
  counts=counts*sf[0]
  hdf_sd_endaccess, sid

  ;use coefficients to return brightness temperatures:
  sz=size(counts, /dim)
  ch1=(chan-1)*2
  ch2=chan*2-1
  print, ch1, ch2
  tb=rebin(coeff[ch2, *], sz[0], sz[1])*counts+$
		rebin(coeff[ch1, *], sz[0], sz[1])

  ;stop

  ;get non-linear correction coefficients:
  c0_id=HDF_SD_ATTRFIND(id,'CalibrationCurveCoefficient#1')
  if c0_id ne -1 then begin
        c1_id=HDF_SD_ATTRFIND(id,'CalibrationCurveCoefficient#2')
        c2_id=HDF_SD_ATTRFIND(id,'CalibrationCurveCoefficient#3')
        c3_id=HDF_SD_ATTRFIND(id,'CalibrationCurveCoefficient#4')
        c4_id=HDF_SD_ATTRFIND(id,'CalibrationCurveCoefficient#5')
        HDF_SD_ATTRINFO,id,c0_id, DATA=c0
        HDF_SD_ATTRINFO,id,c1_id, DATA=c1
        HDF_SD_ATTRINFO,id,c2_id, DATA=c2
        HDF_SD_ATTRINFO,id,c3_id, DATA=c3
        HDF_SD_ATTRINFO,id,c4_id, DATA=c4
        c0=_split_nl_coefficients(c0)
        c1=_split_nl_coefficients(c1)
        c2=_split_nl_coefficients(c2)
        c3=_split_nl_coefficients(c3)
        c4=_split_nl_coefficients(c4)

    ;apply all of these fucking things...
    tb=c0[chan-1]+tb*(c1[chan-1]+tb*(c2[chan-1]+tb*(c3[chan-1]+c4[chan-1]*tb)))
  endif 

  ;get the other polarisation (and apply calibration correction):
  if arg_present(tbh) then begin
    hind=2*((fix(chan)+1)/2)
    vind=hind-1

    if (chan eq vind) then begin
      tbv=tb
      chan=hind
    endif else begin
      tbh=tb
      chan=vind
    endelse

    sid=hdf_sd_select(id, chan)
    hdf_sd_getinfo, sid, name=name
    print, name
    hdf_sd_attrinfo, sid, 0, data=sf
    hdf_sd_attrinfo, sid, 1, data=unit
    print, unit

    hdf_sd_getdata, sid, counts
    counts=counts*sf[0]
    hdf_sd_endaccess, sid

    ;use coefficients to return brightness temperatures:
    sz=size(counts, /dim)
    ch1=(chan-1)*2
    ch2=chan*2-1
    print, ch1, ch2
    tb=rebin(coeff[ch2, *], sz[0], sz[1])*counts+$
		rebin(coeff[ch1, *], sz[0], sz[1])
    if n_elements(c0) ne 0 then begin
      tb=c0[chan-1]+tb*(c1[chan-1]+tb*(c2[chan-1]+tb*(c3[chan-1]+c4[chan-1]*tb)))
    endif

    if chan eq vind then tbv=tb else tbh=tb

    ;more calibration coefficients:
    ;;Coefficients
    Avvind=HDF_SD_ATTRFIND(id, 'CoefficientAvv')
    HDF_SD_ATTRINFO,id,Avvind, DATA=Avv
    Avv=_split_coefficients(Avv)

    Ahvind=HDF_SD_ATTRFIND(id, 'CoefficientAhv')
    HDF_SD_ATTRINFO,id,Ahvind, DATA=Ahv
    Ahv=_split_coefficients(Ahv)

    Aovind=HDF_SD_ATTRFIND(id, 'CoefficientAov')
    HDF_SD_ATTRINFO,id,Aovind, DATA=Aov
    Aov=_split_coefficients(Aov)

    Ahhind=HDF_SD_ATTRFIND(id, 'CoefficientAhh')
    HDF_SD_ATTRINFO,id,Ahhind, DATA=Ahh
    Ahh=_split_coefficients(Ahh)

    Avhind=HDF_SD_ATTRFIND(id, 'CoefficientAvh')
    HDF_SD_ATTRINFO,id,Avhind, DATA=Avh
    Avh=_split_coefficients(Avh)

    Aohind=HDF_SD_ATTRFIND(id, 'CoefficientAoh')
    HDF_SD_ATTRINFO,id,Aohind, DATA=Aoh
    Aoh=_split_coefficients(Aoh)

    tbv=avv[vind]*tbv+ahv[vind]*tbh+2.7*aov[vind]
    tbh=ahh[vind]*tbh+avh[vind]*tbv+2.7*aoh[vind]
  endif

endif

;get longitudes and latitudes:
if arg_present(lon) then begin
  if chan eq 15 or chan eq 16 then begin
    lonid=28
    latid=27
  endif else begin
    lonid=26
    latid=25
  endelse

  sid=hdf_sd_select(id, latid)
  hdf_sd_attrinfo, sid, 0, data=sf
  hdf_sd_attrinfo, sid, 1, data=unit
  print, unit

  hdf_sd_getdata, sid, lat
  lat=lat*sf[0]
  hdf_sd_endaccess, sid

  sid=hdf_sd_select(id, lonid)
  hdf_sd_attrinfo, sid, 0, data=sf
  hdf_sd_attrinfo, sid, 1, data=unit
  print, unit

  hdf_sd_getdata, sid, lon
  lon=lon*sf[0]
  hdf_sd_endaccess, sid

  if chan lt 13 then begin
    sz2=size(lon, /dim)
    ind=lindgen(sz2[0]/2)*2
    lat=lat[ind, *]
    lon=lon[ind, *]
  endif
endif

hdf_sd_end, id

end

