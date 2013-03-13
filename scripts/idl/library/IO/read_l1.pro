FUNCTION _split_coefficients, coeff
RETURN, FLOAT(STRSPLIT(coeff,'([ABG]+\-|,)',/EXTRACT,/REGEX))
END ;_split_coefficients

FUNCTION _split_nl_coefficients, nl_coeff
RETURN, float(STRSPLIT(nl_coeff,'[0-9]{1,2}G[AB]?[VH]-',/EXTRACT,/regex))
END ;_split_coefficients


;;-----------------------------------------------------------------------------
;; ----------------------------------------------------------------------------

FUNCTION read_l1, datafile, counts=counts, SMOOTH=smooth, LIN=lin, nogeoloc=nogeoloc, crop=crop, compare=compare

;+
; NAME:
;   read_l1a
;
; Reads in Level 1A AMSR-E data and returns a structure with
; brightness temperatures calculated from the data (see end of this
; program for definition of this structure)
;
; INPUT: datafile   :  Name with path of the l1a file
;
; OUTPUT: ae        : Structure with TBs, lat, lon etc.
;-
; Gunnar Spreen, September 2003
; LML, 2005-04-19 introduced correct _split_coefficients to split
;      coefficients correctly: - is minus, not separator.
; SEE, 2006-01-13 read SensorCounts additional.
; GH,  2006-02-21 additional keyword:
; Jens Borgmann, August 2007
;      combine read_l1a and read_l1b routines wil all comfort of
;      the read_l1a routine (swathname, qualityflag,
;      equatorcrossingtime) and automatic detecting using the
;      processlevelid
;
; smooth
;      This keyword has TWO purposes:
;      1. control smoothing of calibration constants in
;         field coeff over smooth scan lines. Typical
;         value: smooth = 80.
;         Returned Tbs are caculated with smoothed
;         calibrations constants. Smoothing is done with IDL 
;         routine SMOOTH
;         smooth = 1: No smoothing.
;         
;      2. If smooth is set to a nonpositive value, a plot of
;         the calibration coefficients is written to file
;         caliplot.ps. ATTENTION: Although the plot shows
;         also the smoothed calibtation constants, the
;         returned TBs are caculated using the original,
;         unsmoothed calibration constants.
;         smooth = -1: Make plot of unsmoothed values only
;                         

if not keyword_set(counts)   then counts  =0 else counts  =1
if not keyword_set(smooth  ) then smooth  =0
if not keyword_set(swathtv ) then swathtv =0
if not keyword_set(lin  ) then lin  =0
if not keyword_set(nogeoloc) then nogeoloc=0;  1=old 0=new
if not keyword_set(crop) then crop=0;  1=old 0=new
if not keyword_set(compare) then compare=0
PRINT, 'read_l1:'
PRINT, 'smooth=', smooth



sw_name=strmid(datafile,strlen(datafile)-20,5)
date=strmid(datafile,strlen(datafile)-26,6)
mon=fix(strmid(datafile,strlen(datafile)-24,2))
day=fix(strmid(datafile,strlen(datafile)-22,2))

;-----------------------------------------
;;use new geolocation for the ivorygull maps
;year=strmid(datafile,strlen(datafile)-26,2)
;jd=julday(mon,day,'20'+year)
;syst=systime(/JULIAN,/UTC)-5
;if jd LT syst then nogeoloc=0
;-----------------------------------------

; Read in data
PRINT, 'Reading data'
PRINT, 'Testing if datafile is readable:'
result = {name:'', exist:0L, read:0L, size:-1L}
openr, lun, datafile, /get_lun, error=error
if (error eq 0) then begin
  finfo = fstat(lun)
  result.name = finfo.name
  result.size = finfo.size
  result.read = 1
  free_lun, lun
endif
help,result,/struct

IF result.read EQ 1 THEN BEGIN
    PRINT, '# datafile is detected as readable #'
ENDIF ELSE BEGIN
    PRINT, '# datafile is detected as not readable #'
ENDELSE

; Establish error handler. When errors occur, the index of the
; error is returned in the variable Error_status:
CATCH, Error_status
;This statement begins the error handler:
IF Error_status NE 0 THEN BEGIN
    PRINT, '##################################################################################'
    PRINT, 'Error index: ', Error_status
    PRINT, 'Error message: ', !ERROR_STATE.MSG
    PRINT, '##################################################################################'
    RETURN,-1
ENDIF




sd_id=HDF_SD_START(datafile,/READ)

;detecting processing level
Plid=HDF_SD_ATTRFIND(sd_id, 'ProcessingLevelID')
HDF_SD_ATTRINFO,sd_id,Plid, DATA=plevel
processinglevel=STRMID(plevel,0,3)
print,'Processing Level:',processinglevel
;

IF processinglevel EQ 'L1A' THEN BEGIN


    ;;Coefficients
    Avvind=HDF_SD_ATTRFIND(sd_id, 'CoefficientAvv')
    HDF_SD_ATTRINFO,sd_id,Avvind, DATA=Avv
    Avv=_split_coefficients(Avv)

    Ahvind=HDF_SD_ATTRFIND(sd_id, 'CoefficientAhv')
    HDF_SD_ATTRINFO,sd_id,Ahvind, DATA=Ahv
    Ahv=_split_coefficients(Ahv)

    Aovind=HDF_SD_ATTRFIND(sd_id, 'CoefficientAov')
    HDF_SD_ATTRINFO,sd_id,Aovind, DATA=Aov
    Aov=_split_coefficients(Aov)

    Ahhind=HDF_SD_ATTRFIND(sd_id, 'CoefficientAhh')
    HDF_SD_ATTRINFO,sd_id,Ahhind, DATA=Ahh
    Ahh=_split_coefficients(Ahh)

    Avhind=HDF_SD_ATTRFIND(sd_id, 'CoefficientAvh')
    HDF_SD_ATTRINFO,sd_id,Avhind, DATA=Avh
    Avh=_split_coefficients(Avh)

    Aohind=HDF_SD_ATTRFIND(sd_id, 'CoefficientAoh')
    HDF_SD_ATTRINFO,sd_id,Aohind, DATA=Aoh
    Aoh=_split_coefficients(Aoh)

    IF lin EQ 0 THEN BEGIN
        c0_id=HDF_SD_ATTRFIND(sd_id,'CalibrationCurveCoefficient#1')
        c1_id=HDF_SD_ATTRFIND(sd_id,'CalibrationCurveCoefficient#2')
        c2_id=HDF_SD_ATTRFIND(sd_id,'CalibrationCurveCoefficient#3')
        c3_id=HDF_SD_ATTRFIND(sd_id,'CalibrationCurveCoefficient#4')
        c4_id=HDF_SD_ATTRFIND(sd_id,'CalibrationCurveCoefficient#5')
        HDF_SD_ATTRINFO,sd_id,c0_id, DATA=c0
        HDF_SD_ATTRINFO,sd_id,c1_id, DATA=c1
        HDF_SD_ATTRINFO,sd_id,c2_id, DATA=c2
        HDF_SD_ATTRINFO,sd_id,c3_id, DATA=c3
        HDF_SD_ATTRINFO,sd_id,c4_id, DATA=c4
        c0=_split_nl_coefficients(c0)
        c1=_split_nl_coefficients(c1)
        c2=_split_nl_coefficients(c2)
        c3=_split_nl_coefficients(c3)
        c4=_split_nl_coefficients(c4)
    ENDIF

    ;;read AutomaticQAFlag
    qaid=HDF_SD_ATTRFIND(sd_id, 'AutomaticQAFlag')
    HDF_SD_ATTRINFO,sd_id,qaid, DATA=qaflag
    qaflag=strmid(QAFlag,0,4);;FAIL or PASS
    
    ;;read Equator Crossing Time
    ectid=HDF_SD_ATTRFIND(sd_id, 'EquatorCrossingTime')
    HDF_SD_ATTRINFO,sd_id,ectid, DATA=ect
    ect=STRMID(ect,0,11)        ; because 

    ;;read Scan_Time
    ;;only for reading of the scan_time
    fid=hdf_open(datafile,/READ)
    refur=hdf_vd_find(fid,'Scan_Time')
    vd=hdf_vd_attach(fid,refur)
    uread=hdf_vd_read(vd,x)
    hdf_close,fid
    scan_time=DOUBLE(x)
    
    IF NOT nogeoloc THEN BEGIN

        print,'NEW GEOLOCATION..........'
        ;;navigation data
        sds_id=HDF_SD_SELECT(sd_id,23) & HDF_SD_GETDATA,sds_id,navigation_data
        ;;attitude data
        sds_id=HDF_SD_SELECT(sd_id,24) & HDF_SD_GETDATA,sds_id,attitude_data
        ;;data quality
        sds_id=HDF_SD_SELECT(sd_id,37) & HDF_SD_GETDATA,sds_id,data_quality
        ;;
        SECONDS_PER_DAY=24.*3600
        j1993t0=julday(1,1,1993,0,0,0) ; julian day of the date 1.1.1993
        j2006t0=julday(1,1,2006,0,0,0) ; julian day
        j2006_time = (j2006t0 - j1993t0)*SECONDS_PER_DAY
        jd_utc=dblarr(n_elements(scan_time))
        idt0=where(scan_time LT j2006_time)
        IF idt0[0] NE -1 THEN jd_utc[idt0]=j1993t0 + (scan_time[idt0] - 5.0)/SECONDS_PER_DAY
        idt1=where(scan_time GE j2006_time)
        IF idt1[0] NE -1 THEN jd_utc[idt1]=j1993t0 + (scan_time[idt1] - 6.0)/SECONDS_PER_DAY

        navigation=navigation_data
        attitude=attitude_data*0.017453292519943295
        tacho_offset=data_quality[3,*]-45.0 ; python data_quality[:,3]
        tacho_offset=tacho_offset[*]

        ;;write L1TMP
        tmpdir='/tmp'
        l1tmp=tmpdir+'/L1TMP'

        OPENW,outunit,l1tmp,/GET_LUN
        sjd=(size(jd_utc))[1]

        IF CROP THEN BEGIN
            ;crop done by lothar
            WRITEU,outunit,[326],[FIX(sjd)-21+1]
            WRITEU,outunit,[jd_utc[10:sjd-11]]
            WRITEU,outunit,[tacho_offset[10:sjd-11]]
            WRITEU,outunit,[navigation[*,10:sjd-11]]
            WRITEU,outunit,[attitude[*,10:sjd-11]]
        ENDIF ELSE BEGIN
            WRITEU,outunit,[326],[FIX(sjd)]
            WRITEU,outunit,[jd_utc]
            WRITEU,outunit,[tacho_offset]
            WRITEU,outunit,[navigation]
            WRITEU,outunit,[attitude]
        ENDELSE

        WRITEU,outunit,[326]
        FREE_LUN,outunit
    
        ;; call external c-routine pte2
        nadir='47.082'
        scan='75.410'
        roll='0.09'
        fname_lon=tmpdir+'/lon'
        fname_lat=tmpdir+'/lat'
        int_time='1.3e-3'
        offset='47'             ; cuts off l1a data
        lines='392'
        c_dir='/home/borgmann/python/geolocation/c'

        SPAWN,'echo 1 '+nadir+' '+int_time+' 240.0 '+offset+' '+lines+' '+scan+' '+roll+' | '+c_dir+'/pte2 '+l1tmp+' '+fname_lon+' '+fname_lat
    
        ;;read file produce by the pte2 c-code

        IF CROP THEN BEGIN
            ;crop done by lothar
            lat=intarr(392,FIX(sjd)-21)
            lon=intarr(392,FIX(sjd)-21)
        ENDIF ELSE BEGIN
            ;crop the last scanline
            ;last crop of lothars routine
            lat=intarr(392,FIX(sjd)-1)
            lon=intarr(392,FIX(sjd)-1)
        ENDELSE
        openr,inunit,fname_lat,/get_lun
        readu,inunit,lat
        free_lun,inunit

        lat=lat/128.

        openr,inunit,fname_lon,/get_lun
        readu,inunit,lon
        free_lun,inunit
        
        lon=lon/128.
    ENDIF  ; not nogeoloc 

    ;;Antenna Temperature Coefficient (Offset + Slope)
    sds_id=HDF_SD_SELECT(sd_id,21)
    HDF_SD_GETDATA,sds_id,coeff
    ;;6 GHz channels
    sds_id=HDF_SD_SELECT(sd_id,1)
    HDF_SD_GETDATA,sds_id,tbv06
    sds_id=HDF_SD_SELECT(sd_id,2)
    HDF_SD_GETDATA,sds_id,tbh06
    ;;10
    sds_id=HDF_SD_SELECT(sd_id,3)
    HDF_SD_GETDATA,sds_id,tbv10
    sds_id=HDF_SD_SELECT(sd_id,4)
    HDF_SD_GETDATA,sds_id,tbh10
    ;;18
    sds_id=HDF_SD_SELECT(sd_id,5)
    HDF_SD_GETDATA,sds_id,tbv18
    sds_id=HDF_SD_SELECT(sd_id,6)
    HDF_SD_GETDATA,sds_id,tbh18
    ;;23
    sds_id=HDF_SD_SELECT(sd_id,7)
    HDF_SD_GETDATA,sds_id,tbv23
    sds_id=HDF_SD_SELECT(sd_id,8)
    HDF_SD_GETDATA,sds_id,tbh23
    ;;36
    sds_id=HDF_SD_SELECT(sd_id,9)
    HDF_SD_GETDATA,sds_id,tbv36
    sds_id=HDF_SD_SELECT(sd_id,10)
    HDF_SD_GETDATA,sds_id,tbh36
    ;;89
    sds_id=HDF_SD_SELECT(sd_id,13)
    HDF_SD_GETDATA,sds_id,tbv89a
    sds_id=HDF_SD_SELECT(sd_id,15)
    HDF_SD_GETDATA,sds_id,tbv89b
    sds_id=HDF_SD_SELECT(sd_id,14)
    HDF_SD_GETDATA,sds_id,tbh89a
    sds_id=HDF_SD_SELECT(sd_id,16)
    HDF_SD_GETDATA,sds_id,tbh89b
    ;;lat/lon
    sds_id=HDF_SD_SELECT(sd_id,25)
    HDF_SD_GETDATA,sds_id,lat_a
    sds_id=HDF_SD_SELECT(sd_id,26)
    HDF_SD_GETDATA,sds_id,lon_a
    sds_id=HDF_SD_SELECT(sd_id,27)
    HDF_SD_GETDATA,sds_id,lat_b
    sds_id=HDF_SD_SELECT(sd_id,28)
    HDF_SD_GETDATA,sds_id,lon_b
    ;;land
    sds_id=HDF_SD_SELECT(sd_id,33)
    HDF_SD_GETDATA,sds_id,land
    HDF_SD_ENDACCESS,sds_id

    HDF_SD_END, sd_id

    ;; Get Brightness Temperatures in correct format
    ;; Ta  = offset  + slope*count
    ;; INSERT HERE NONLINEAR CORRECTION for data revision 2 onwards!!
    ;; Tbv = Avv*Tav + Ahv*Tah + 2.7*Aov
    ;; Tbh = Ahh*Tah + Avh*Tav + 2.7*Aoh

    ;; Meaning of first index of coeff:
    ;;     offset count
    ;; 06v   0      1    
    ;; 06h   2      3
    ;; 10v   4      5
    ;; 10h   6      7
    ;; 18v   8      9
    ;; 18h  10     11
    ;; 23v  12     13
    ;; 23h  14     15
    ;; 36v  16     17
    ;; 36h  18     19
    ;; 89vA 24     25
    ;; 89hA 26     27
    ;; 89vB 28     29
    ;; 89hB 30     31

    ;; smooth calibration coefficients over several scanlines
    IF smooth GT 1 THEN BEGIN
        FOR i = 0, 19 DO BEGIN  ; low frequency channels
            coeff[i,*] = SMOOTH(coeff[i,*], smooth,/EDGE_TRUNCATE)
        ENDFOR
        FOR i = 24, 31 DO BEGIN ; 89 GHz channels
            coeff[i,*] = SMOOTH(coeff[i,*], smooth,/EDGE_TRUNCATE)
        ENDFOR
    ENDIF

    ;;06 GHz channels
    sz=SIZE(tbv06)
    if counts eq 0 then begin $
      tav06= REBIN(coeff[0,*],sz[1],sz[2],/SAMPLE)+REBIN(coeff[1,*],sz[1],sz[2],/SAMPLE)*tbv06
        tah06= REBIN(coeff[2,*],sz[1],sz[2],/SAMPLE)+REBIN(coeff[3,*],sz[1],sz[2],/SAMPLE)*tbh06
        if lin eq 0 then begin
            tav06=c0[0]+c1[0]*tav06+c2[0]*tav06^2+c3[0]*tav06^3+c4[0]*tav06^4
            tah06=c0[1]+c1[1]*tah06+c2[1]*tah06^2+c3[1]*tah06^3+c4[1]*tah06^4
        endif
        
        tbv06= Avv[1]*tav06+Ahv[1]*tah06+2.7*Aov[1]
        tbh06= Ahh[1]*tah06+Avh[1]*tav06+2.7*Aoh[1]
    endif else begin $
      tbv06=float(tbv06)  ;; Ausgabe der SensorCounts für v06.
        tbh06=float(tbh06)  ;; Ausgabe der SensorCounts für h06.
    endelse
    missdata=WHERE(tbv06 LE -999)
    IF missdata[0] NE -1 THEN tbv06[missdata] = !VALUES.F_NAN
    missdata=WHERE(tbh06 LE -999)
    IF missdata[0] NE -1 THEN tbh06[missdata] = !VALUES.F_NAN

    ;;10
    sz=SIZE(tbv10)
    if counts eq 0 then begin $
      tav10= REBIN(coeff[4,*],sz[1],sz[2],/SAMPLE)+REBIN(coeff[5,*],sz[1],sz[2],/SAMPLE)*tbv10
        tah10= REBIN(coeff[6,*],sz[1],sz[2],/SAMPLE)+REBIN(coeff[7,*],sz[1],sz[2],/SAMPLE)*tbh10
        if lin eq 0 then begin
            tav10=c0[2]+c1[2]*tav10+c2[2]*tav10^2+c3[2]*tav10^3+c4[2]*tav10^4
            tah10=c0[3]+c1[3]*tah10+c2[3]*tah10^2+c3[3]*tah10^3+c4[3]*tah10^4
        endif
        
        tbv10= Avv[3]*tav10+Ahv[3]*tah10+2.7*Aov[3]
        tbh10= Ahh[3]*tah10+Avh[3]*tav10+2.7*Aoh[3]
    endif else begin $
      tbv10=float(tbv10)  ;; Ausgabe der SensorCounts für v10.
        tbh10=float(tbh10)  ;; Ausgabe der SensorCounts für h10.
    endelse
    missdata=WHERE(tbv10 LE -999)
    IF missdata[0] NE -1 THEN tbv10[missdata] = !VALUES.F_NAN
    missdata=WHERE(tbh10 LE -999)
    IF missdata[0] NE -1 THEN tbh10[missdata] = !VALUES.F_NAN
    
    ;;18
    sz=SIZE(tbv18)
    if counts eq 0 then begin $
      tav18= REBIN(coeff[8,*],sz[1],sz[2],/SAMPLE)+REBIN(coeff[9,*],sz[1],sz[2],/SAMPLE)*tbv18
        tah18= REBIN(coeff[10,*],sz[1],sz[2],/SAMPLE)+REBIN(coeff[11,*],sz[1],sz[2],/SAMPLE)*tbh18
        if lin eq 0 then begin
            tav18=c0[4]+c1[4]*tav18+c2[4]*tav18^2+c3[4]*tav18^3+c4[4]*tav18^4
            tah18=c0[5]+c1[5]*tah18+c2[5]*tah18^2+c3[5]*tah18^3+c4[5]*tah18^4
        endif

        tbv18= Avv[5]*tav18+Ahv[5]*tah18+2.7*Aov[5]
        tbh18= Ahh[5]*tah18+Avh[5]*tav18+2.7*Aoh[5]
    endif else begin $
      tbv18=float(tbv18)  ;; Ausgabe der SensorCounts für v18.
        tbh18=float(tbh18)  ;; Ausgabe der SensorCounts für h18.
    endelse
    missdata=WHERE(tbv18 LE -999)
    IF missdata[0] NE -1 THEN tbv18[missdata] = !VALUES.F_NAN
    missdata=WHERE(tbh18 LE -999)
    IF missdata[0] NE -1 THEN tbh18[missdata] = !VALUES.F_NAN

    ;;23
    sz=SIZE(tbv23)
    if counts eq 0 then begin $
      tav23= REBIN(coeff[12,*],sz[1],sz[2],/SAMPLE)+REBIN(coeff[13,*],sz[1],sz[2],/SAMPLE)*tbv23
        tah23= REBIN(coeff[14,*],sz[1],sz[2],/SAMPLE)+REBIN(coeff[15,*],sz[1],sz[2],/SAMPLE)*tbh23
        if lin eq 0 then begin
            tav23=c0[6]+c1[6]*tav23+c2[6]*tav23^2+c3[6]*tav23^3+c4[6]*tav23^4
            tah23=c0[7]+c1[7]*tah23+c2[7]*tah23^2+c3[7]*tah23^3+c4[7]*tah23^4
    endif
    
    tbv23= Avv[7]*tav23+Ahv[7]*tah23+2.7*Aov[7]
    tbh23= Ahh[7]*tah23+Avh[7]*tav23+2.7*Aoh[7]
    endif else begin $
      tbv23=float(tbv23)  ;; Ausgabe der SensorCounts für v23.
        tbh23=float(tbh23)  ;; Ausgabe der SensorCounts für h23.
    endelse
    missdata=WHERE(tbv23 LE -999)
    IF missdata[0] NE -1 THEN tbv23[missdata] = !VALUES.F_NAN
    missdata=WHERE(tbh23 LE -999)
    IF missdata[0] NE -1 THEN tbh23[missdata] = !VALUES.F_NAN
          
    ;;36
    sz=SIZE(tbv36)
    if counts eq 0 then begin $
      tav36= REBIN(coeff[16,*],sz[1],sz[2],/SAMPLE)+REBIN(coeff[17,*],sz[1],sz[2],/SAMPLE)*tbv36
        tah36= REBIN(coeff[18,*],sz[1],sz[2],/SAMPLE)+REBIN(coeff[19,*],sz[1],sz[2],/SAMPLE)*tbh36
        if lin eq 0 then begin
            tav36=c0[8]+c1[8]*tav36+c2[8]*tav36^2+c3[8]*tav36^3+c4[8]*tav36^4
            tah36=c0[9]+c1[9]*tah36+c2[9]*tah36^2+c3[9]*tah36^3+c4[9]*tah36^4
        endif
        
        tbv36= Avv[9]*tav36+Ahv[9]*tah36+2.7*Aov[9]
        tbh36= Ahh[9]*tah36+Avh[9]*tav36+2.7*Aoh[9]
    endif else begin $
      tbv36=float(tbv36)  ;; Ausgabe der SensorCounts für v36.
        tbh36=float(tbh36)  ;; Ausgabe der SensorCounts für h36.
    endelse
    missdata=WHERE(tbv36 LE -999)
    IF missdata[0] NE -1 THEN tbv36[missdata] = !VALUES.F_NAN
    missdata=WHERE(tbh36 LE -999)
    IF missdata[0] NE -1 THEN tbh36[missdata] = !VALUES.F_NAN
    
    ;;89
    sz=SIZE(tbv89a)
    if counts eq 0 then begin $
      tav89a= REBIN(coeff[24,*],sz[1],sz[2],/SAMPLE)+REBIN(coeff[25,*],sz[1],sz[2],/SAMPLE)*tbv89a
        tah89a= REBIN(coeff[26,*],sz[1],sz[2],/SAMPLE)+REBIN(coeff[27,*],sz[1],sz[2],/SAMPLE)*tbh89a
        if lin eq 0 then begin
            tav89a=c0[12]+c1[12]*tav89a+c2[12]*tav89a^2+c3[12]*tav89a^3+c4[12]*tav89a^4
            tah89a=c0[13]+c1[13]*tah89a+c2[13]*tah89a^2+c3[13]*tah89a^3+c4[13]*tah89a^4
        endif
        
        tbv89a= Avv[13]*tav89a+Ahv[13]*tah89a+2.7*Aov[13]
        tbh89a= Ahh[13]*tah89a+Avh[13]*tav89a+2.7*Aoh[13]
        tav89b= REBIN(coeff[28,*],sz[1],sz[2],/SAMPLE)+REBIN(coeff[29,*],sz[1],sz[2],/SAMPLE)*tbv89b
        tah89b= REBIN(coeff[30,*],sz[1],sz[2],/SAMPLE)+REBIN(coeff[31,*],sz[1],sz[2],/SAMPLE)*tbh89b
        if lin eq 0 then begin
            tav89b=c0[14]+c1[14]*tav89b+c2[14]*tav89b^2+c3[14]*tav89b^3+c4[14]*tav89b^4
            tah89b=c0[15]+c1[15]*tah89b+c2[15]*tah89b^2+c3[15]*tah89b^3+c4[15]*tah89b^4
        endif
        
        
        tbv89b= Avv[15]*tav89b+Ahv[15]*tah89b+2.7*Aov[15]
        tbh89b= Ahh[15]*tah89b+Avh[15]*tav89b+2.7*Aoh[15]
                                ;PRINT, 'tbv89a[200:201,1000:1003]:'
                                ;print,tbv89a[200:202,1000:1002]
    endif else begin $
      tbv89a=float(tbv89a)  ;; Ausgabe der SensorCounts für v89a.
        tbh89a=float(tbh89a)  ;; Ausgabe der SensorCounts für h89a.
        tbv89b=float(tbv89b)  ;; Ausgabe der SensorCounts für v89b.
        tbh89b=float(tbh89b)  ;; Ausgabe der SensorCounts für h89b.
    endelse
    missdata=WHERE(tbv89a LE -999)
    IF missdata[0] NE -1 THEN tbv89a[missdata] = !VALUES.F_NAN
    missdata=WHERE(tbv89b LE -999)
    IF missdata[0] NE -1 THEN tbv89b[missdata] = !VALUES.F_NAN 
    missdata=WHERE(tbh89a LE -999)
    IF missdata[0] NE -1 THEN tbh89a[missdata] = !VALUES.F_NAN
    missdata=WHERE(tbh89b LE -999)
    IF missdata[0] NE -1 THEN tbh89b[missdata] = !VALUES.F_NAN
    


    ;; Here do analysis and display of calibration constants
    
    IF smooth LT 0 THEN BEGIN
        PRINT, 'Plot coeffs, return unsmoothed TBs, smooth window=', -smooth    
        ;;    smooth=ABS(smooth)
        LOADCT, 13
        HELP, tbh06, tbh89b, coeff
        PRINT, coeff[*, 0]
        
        SET_PLOT, "PS"
        !P.FONT=0
        DEVICE, FILENAME='caliplot.ps'
        DEVICE, /PORTRAIT,/COLOR,BITS=8
        DEVICE, /TIMES,/ISOLATIN1
        DEVICE, XSIZE=20,YSIZE=30
        DEVICE, /ENCAPSULATED
        !P.THICK=2.0
        !X.THICK=2.0
        !Y.THICK=2.0
        !P.CHARTHICK=1.2
        !P.CHARSIZE=1.2
        !P.SYMSIZE=1
        !P.MULTI=[0,1,2]
        
        cf = 255/7              ; color factor. 6 = # colors needed here

        ;;   Plot slope coefficients
        PLOT, [0,0], XRANGE=[0,2000], YRANGE=[0.395,0.415], $
              TITLE=STRMID(datafile, STRLEN(datafile)-31,22) + ':  Calibration Slopes'
        legend=['6GHz', '  10', ' 18', '23', '36', '89A', ' 89B']
        i_coef=[  0,    4,    8,   12,   16,   24,    28  ] 
        i_col =[ 20,   40,   90,  110,  160,  220,  240  ] 
        coeff[1,*] = coeff[1,*] +0.31 ; additional offset
        coeff[3,*] = coeff[3,*] +0.31 ; for 6v, 6h
        FOR i=0,6 DO BEGIN
            OPLOT, coeff[i_coef[i]+1,*], COLOR=i_col[i]
            OPLOT, coeff[i_coef[i]+3,*], COLOR=i_col[i], LINESTYLE=1 ; h
            IF smooth LT -1 THEN BEGIN
                OPLOT, SMOOTH(coeff[i_coef[i]+1,*],-smooth,/EDGE_TRUNCATE), COLOR=i_col[i] ; v smoothed
                OPLOT, SMOOTH(coeff[i_coef[i]+3,*],-smooth,/EDGE_TRUNCATE), COLOR=i_col[i] ; h smoothed
            ENDIF
            XYOUTS, 200+150*i, 0.3955, legend[i], COLOR=i_col[i] ; legend
        ENDFOR
        
        ;;   Plot offset coefficients
        PLOT, [0,0], XRANGE=[0,2000], YRANGE=[145,185.], YSTYLE=1, TITLE='Calibration Offsets'
        FOR i=0,6 DO BEGIN
            OPLOT, coeff[i_coef[i]  ,*], COLOR=i_col[i] ; v
            OPLOT, coeff[i_coef[i]+2,*], COLOR=i_col[i], LINESTYLE=1 ; h
            IF smooth LT -1 THEN BEGIN
                OPLOT, SMOOTH(coeff[i_coef[i]  ,*],-smooth,/EDGE_TRUNCATE), COLOR=i_col[i] ; v smoothed
                OPLOT, SMOOTH(coeff[i_coef[i]+2,*],-smooth,/EDGE_TRUNCATE), COLOR=i_col[i] ; h smoothed
            ENDIF
            XYOUTS, 200+150*i, 149,   legend[i], COLOR=i_col[i] ; legend
        ENDFOR
        
        
        
        DEVICE, /CLOSE
        SET_PLOT, 'X'
    ENDIF
        
ENDIF


IF processinglevel EQ 'L1B' THEN BEGIN
    sw_name=strmid(datafile,strlen(datafile)-20,5)
    date=strmid(datafile,strlen(datafile)-26,6)
    mon=fix(strmid(datafile,strlen(datafile)-24,2))
    day=fix(strmid(datafile,strlen(datafile)-22,2))
    
    
    ;sd_id=HDF_SD_START(datafile,/READ)
    sds_id=HDF_SD_SELECT(sd_id,1) & HDF_SD_GETDATA,sds_id,tbv06
    sds_id=HDF_SD_SELECT(sd_id,2) & HDF_SD_GETDATA,sds_id,tbh06
    sds_id=HDF_SD_SELECT(sd_id,3) & HDF_SD_GETDATA,sds_id,tbv10
    sds_id=HDF_SD_SELECT(sd_id,4) & HDF_SD_GETDATA,sds_id,tbh10
    sds_id=HDF_SD_SELECT(sd_id,5) & HDF_SD_GETDATA,sds_id,tbv18
    sds_id=HDF_SD_SELECT(sd_id,6) & HDF_SD_GETDATA,sds_id,tbh18
    sds_id=HDF_SD_SELECT(sd_id,7) & HDF_SD_GETDATA,sds_id,tbv23
    sds_id=HDF_SD_SELECT(sd_id,8) & HDF_SD_GETDATA,sds_id,tbh23
    sds_id=HDF_SD_SELECT(sd_id,9) & HDF_SD_GETDATA,sds_id,tbv36
    sds_id=HDF_SD_SELECT(sd_id,10) & HDF_SD_GETDATA,sds_id,tbh36
    sds_id=HDF_SD_SELECT(sd_id,13) & HDF_SD_GETDATA,sds_id,tbv89a
    sds_id=HDF_SD_SELECT(sd_id,14) & HDF_SD_GETDATA,sds_id,tbh89a
    sds_id=HDF_SD_SELECT(sd_id,15) & HDF_SD_GETDATA,sds_id,tbv89b
    sds_id=HDF_SD_SELECT(sd_id,16) & HDF_SD_GETDATA,sds_id,tbh89b
    
    sds_id=HDF_SD_SELECT(sd_id,25) & HDF_SD_GETDATA,sds_id,lat_a
    sds_id=HDF_SD_SELECT(sd_id,26) & HDF_SD_GETDATA,sds_id,lon_a
    sds_id=HDF_SD_SELECT(sd_id,27) & HDF_SD_GETDATA,sds_id,lat_b
    sds_id=HDF_SD_SELECT(sd_id,28) & HDF_SD_GETDATA,sds_id,lon_b
    
    sds_id=HDF_SD_SELECT(sd_id,33) & HDF_SD_GETDATA,sds_id,land

    ;;read AutomaticQAFlag
    qaid=HDF_SD_ATTRFIND(sd_id, 'AutomaticQAFlag')
    HDF_SD_ATTRINFO,sd_id,qaid, DATA=qaflag
    qaflag=strmid(QAFlag,0,4);;FAIL or PASS
    
    ;;read Equator Crossing Time
    ectid=HDF_SD_ATTRFIND(sd_id, 'EquatorCrossingTime')
    HDF_SD_ATTRINFO,sd_id,ectid, DATA=ect
    ect=STRMID(ect,0,11)        ; because 
    
    ;;read Scan_Time
    ;;only for reading of the scan_time
    fid=hdf_open(datafile,/READ)
    refur=hdf_vd_find(fid,'Scan_Time')
    vd=hdf_vd_attach(fid,refur)
    uread=hdf_vd_read(vd,x)
    hdf_close,fid
    scan_time=DOUBLE(x)

    IF NOT nogeoloc THEN BEGIN
        print,'NEW GEOLOCATION..........'
        ;;navigation data
        sds_id=HDF_SD_SELECT(sd_id,23) & HDF_SD_GETDATA,sds_id,navigation_data
        ;;attitude data
        sds_id=HDF_SD_SELECT(sd_id,24) & HDF_SD_GETDATA,sds_id,attitude_data
        ;;data quality
        sds_id=HDF_SD_SELECT(sd_id,37) & HDF_SD_GETDATA,sds_id,data_quality
        
        SECONDS_PER_DAY=24.*3600
        j1993t0=julday(1,1,1993,0,0,0); julian day of the date 1.1.1993
        j2006t0=julday(1,1,2006,0,0,0) ; julian day
        j2006_time = (j2006t0 - j1993t0)*SECONDS_PER_DAY
        jd_utc=dblarr(n_elements(scan_time))
        idt0=where(scan_time LT j2006_time)
        IF idt0[0] NE -1 THEN jd_utc[idt0]=j1993t0 + (scan_time[idt0] - 5.0)/SECONDS_PER_DAY
        idt1=where(scan_time GE j2006_time)
        IF idt1[0] NE -1 THEN jd_utc[idt1]=j1993t0 + (scan_time[idt1] - 6.0)/SECONDS_PER_DAY
        
        navigation=navigation_data
        attitude=attitude_data*0.017453292519943295
        tacho_offset=data_quality[3,*]-45.0 ; python data_quality[:,3]
        tacho_offset=tacho_offset[*]
        
        ;;write L1TMP
        tmpdir='/tmp'
        l1tmp=tmpdir+'/L1TMP'
        
        OPENW,outunit,l1tmp,/GET_LUN
        sjd=(size(jd_utc))[1]
        
        

        IF CROP THEN BEGIN
            ;crop done by lothar
            WRITEU,outunit,[326],[FIX(sjd)-21+1]
            WRITEU,outunit,[jd_utc[10:sjd-11]]
            WRITEU,outunit,[tacho_offset[10:sjd-11]]
            WRITEU,outunit,[navigation[*,10:sjd-11]]
            WRITEU,outunit,[attitude[*,10:sjd-11]]
        ENDIF ELSE BEGIN
            WRITEU,outunit,[326],[FIX(sjd)]
            WRITEU,outunit,[jd_utc]
            WRITEU,outunit,[tacho_offset]
            WRITEU,outunit,[navigation]
            WRITEU,outunit,[attitude]
        ENDELSE
        
        WRITEU,outunit,[326]
        FREE_LUN,outunit

    
        ;; call external c-routine pte2
        nadir='47.082'
        scan='75.410'
        roll='0.09'
        fname_lon=tmpdir+'/lon'
        fname_lat=tmpdir+'/lat'
        int_time='1.3e-3'
        offset='47'             
        lines='392'
        c_dir='/home/borgmann/python/geolocation/c'

        SPAWN,'echo 1 '+nadir+' '+int_time+' 240.0 '+offset+' '+lines+' '+scan+' '+roll+' | '+c_dir+'/pte2 '+l1tmp+' '+fname_lon+' '+fname_lat
    
        ;;read file produce by the pte2 c-code
        IF CROP THEN BEGIN
            ;crop done by lothar
            lat=intarr(392,FIX(sjd)-21)
            lon=intarr(392,FIX(sjd)-21)
        ENDIF ELSE BEGIN
            ;crop the last scanline
            ;last crop of lothars routine
            lat=intarr(392,FIX(sjd)-1)
            lon=intarr(392,FIX(sjd)-1)
        ENDELSE

        openr,inunit,fname_lat,/get_lun
        readu,inunit,lat
        free_lun,inunit
        
        lat=lat/128.

        openr,inunit,fname_lon,/get_lun
        readu,inunit,lon
        free_lun,inunit
        
        lon=lon/128.
    ENDIF
    
    HDF_SD_ENDACCESS,sds_id
    
    HDF_SD_END, sd_id

    
    tbv06=tbv06/10.
    ind=WHERE(tbv06 LE -999.9)
    IF ind[0] NE -1 THEN tbv06[ind]=!VALUES.F_NAN
    tbh06=tbh06/10.
    ind=WHERE(tbh06 LE -999.9)
    IF ind[0] NE -1 THEN tbh06[ind]=!VALUES.F_NAN       
    tbv10=tbv10/10.
    ind=WHERE(tbv10 LE -999.9)
    IF ind[0] NE -1 THEN tbv10[ind]=!VALUES.F_NAN
    tbh10=tbh10/10.
    ind=WHERE(tbh10 LE -999.9)
    IF ind[0] NE -1 THEN tbh10[ind]=!VALUES.F_NAN
    tbv18=tbv18/10.
    ind=WHERE(tbv18 LE -999.9)
    IF ind[0] NE -1 THEN tbv18[ind]=!VALUES.F_NAN
    tbh18=tbh18/10.
    ind=WHERE(tbh18 LE -999.9)
    IF ind[0] NE -1 THEN tbh18[ind]=!VALUES.F_NAN
    tbv23=tbv23/10.
    ind=WHERE(tbv23 LE -999.9)
    IF ind[0] NE -1 THEN tbv23[ind]=!VALUES.F_NAN
    tbh23=tbh23/10.
    ind=WHERE(tbh23 LE -999.9)
    IF ind[0] NE -1 THEN tbh23[ind]=!VALUES.F_NAN
    tbv36=tbv36/10.
    ind=WHERE(tbv36 LE -999.9)
    IF ind[0] NE -1 THEN tbv36[ind]=!VALUES.F_NAN
    tbh36=tbh36/10.
    ind=WHERE(tbh36 LE -999.9)
    IF ind[0] NE -1 THEN tbh36[ind]=!VALUES.F_NAN
    tbv89a=tbv89a/10.
    ind=WHERE(tbv89a LE -999.9)
    IF ind[0] NE -1 THEN tbv89a[ind]=!VALUES.F_NAN
    tbv89b=tbv89b/10.
    ind=WHERE(tbv89b LE -999.9)
    IF ind[0] NE -1 THEN tbv89b[ind]=!VALUES.F_NAN 
    tbh89a=tbh89a/10.
    ind=WHERE(tbh89a LE -999.9)
    IF ind[0] NE -1 THEN tbh89a[ind]=!VALUES.F_NAN
    tbh89b=tbh89b/10.
    ind=WHERE(tbh89b LE -999.9)
    IF ind[0] NE -1 THEN tbh89b[ind]=!VALUES.F_NAN
    
ENDIF

IF nogeoloc THEN BEGIN
    ;no geolocation
    ;;Define structures
    tb = {v6:tbv06,                $
          h6:tbh06,                $
          v10:tbv10,               $
          h10:tbh10,               $
          v18:tbv18,               $
          h18:tbh18,               $
          v23:tbv23,               $
          h23:tbh23,               $
          v36:tbv36,               $
          h36:tbh36,               $
          v89a:tbv89a,             $
          h89a:tbh89a,             $
          v89b:tbv89b,             $
          h89b:tbh89b}             

    lat=lat_b/100.
    lon=lon_b/100.
    ae = {tb:tb,            $
          lat_a:lat_a/100., $
          lon_a:lon_a/100., $
          lat_b:lat_b/100., $
          lon_b:lon_b/100., $
          lat:lat, $
          lon:lon, $
          land:land,        $
          date:date,        $
          mon:mon,          $
          day:day,          $
          sw_name:sw_name,  $
          scan_time:scan_time, $	
          ect:ect,          $
          qaflag:qaflag}
    
ENDIF ELSE BEGIN
    ;geolocation
    IF CROP THEN BEGIN
        start_line=10
        tb = {v6:tbv06[*,start_line:start_line+sjd-21-1],            $
              h6:tbh06[*,start_line:start_line+sjd-21-1],                $
              v10:tbv10[*,start_line:start_line+sjd-21-1],               $
              h10:tbh10[*,start_line:start_line+sjd-21-1],               $
              v18:tbv18[*,start_line:start_line+sjd-21-1],               $
              h18:tbh18[*,start_line:start_line+sjd-21-1],               $
              v23:tbv23[*,start_line:start_line+sjd-21-1],               $
              h23:tbh23[*,start_line:start_line+sjd-21-1],               $
              v36:tbv36[*,start_line:start_line+sjd-21-1],               $
              h36:tbh36[*,start_line:start_line+sjd-21-1],               $
              v89a:tbv89a[*,start_line:start_line+sjd-21-1],             $
              h89a:tbh89a[*,start_line:start_line+sjd-21-1],             $
              v89b:tbv89b[*,start_line:start_line+sjd-21-1],             $
              h89b:tbh89b[*,start_line:start_line+sjd-21-1]}             
    

        land=land[*,start_line:start_line+sjd-21-1,*]
        scan_time=scan_time[start_line:start_line+sjd-21-1]

    ENDIF ELSE BEGIN
        IF processinglevel EQ 'L1A' THEN BEGIN
            tb = {v6:tbv06[23:218,0:sjd-2],                $
                  h6:tbh06[23:218,0:sjd-2],                $
                  v10:tbv10[23:218,0:sjd-2],               $
                  h10:tbh10[23:218,0:sjd-2],               $
                  v18:tbv18[23:218,0:sjd-2],               $
                  h18:tbh18[23:218,0:sjd-2],               $
                  v23:tbv23[23:218,0:sjd-2],               $
                  h23:tbh23[23:218,0:sjd-2],               $
                  v36:tbv36[23:218,0:sjd-2],               $
                  h36:tbh36[23:218,0:sjd-2],               $
                  v89a:tbv89a[47:438,0:sjd-2],             $
                  h89a:tbh89a[47:438,0:sjd-2],             $
                  v89b:tbv89b[47:438,0:sjd-2],             $
                  h89b:tbh89b[47:438,0:sjd-2]}             

            
            land=land[23:218,0:sjd-2,*]
        ENDIF
        IF processinglevel EQ 'L1B' THEN BEGIN
            tb = {v6:tbv06[*,0:sjd-2],                $
                  h6:tbh06[*,0:sjd-2],                $
                  v10:tbv10[*,0:sjd-2],               $
                  h10:tbh10[*,0:sjd-2],               $
                  v18:tbv18[*,0:sjd-2],               $
                  h18:tbh18[*,0:sjd-2],               $
                  v23:tbv23[*,0:sjd-2],               $
                  h23:tbh23[*,0:sjd-2],               $
                  v36:tbv36[*,0:sjd-2],               $
                  h36:tbh36[*,0:sjd-2],               $
                  v89a:tbv89a[*,0:sjd-2],             $
                  h89a:tbh89a[*,0:sjd-2],             $
                  v89b:tbv89b[*,0:sjd-2],             $
                  h89b:tbh89b[*,0:sjd-2]}             

             land=land[*,0:sjd-2,*]
        ENDIF
        scan_time=scan_time[0:sjd-2]

    ENDELSE
    ;;since 4 NOV 2004 no 89a data therfore no lat_a, lon_a
    IF COMPARE THEN BEGIN
        ae = {tb:tb,            $
              lat_a:lat_a/100., $
              lon_a:lon_a/100., $
              lat_b:lat_b/100., $
              lon_b:lon_b/100., $
              lat:lat, $
              lon:lon, $
              land:land,        $
              date:date,        $
              mon:mon,          $
              day:day,          $
              sw_name:sw_name,  $
              scan_time:scan_time, $	
              ect:ect,          $
              qaflag:qaflag}
    ENDIF ELSE BEGIN
        ae = {tb:tb,            $
              lat_a:lat, $
              lon_a:lon, $
              lat_b:lat, $
              lon_b:lon, $
              lat:lat, $
              lon:lon, $
              land:land,        $
              date:date,        $
              mon:mon,          $
              day:day,          $
              sw_name:sw_name,  $
              scan_time:scan_time, $	
              ect:ect,          $
              qaflag:qaflag}
    ENDELSE

ENDELSE
RETURN,ae
END
