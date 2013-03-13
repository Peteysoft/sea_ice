PRO MAPXY, X, Y, ALAT, ALONG, SGN, SLat=slat, delta=delta
;+    
;NAME:            MapXY
;    
;PURPOSE:         This subroutine converts from Polar Stereographic (X, Y) coordinates
;                 to geodetic latitude and longitude for the polar regions.
;    
;CATEGORY:        Misc
;    
;USAGE:           MapXY, X, Y, Latitude, Longitude, SGN
;
;INPUTS:          X          (Float) = Polar Stereographic X Coordinate (km)
;                 Y          (Float) = Polar Stereographic Y Coordinate (km)
;                 Latitude   (Float) = Geodetic Latitude  (will be calculated for output)  
;                 Longitude  (Float) = Geodetic Longitude (will be calculated for output)
;                 SGN                = Sign of latitude
;                                      +1 : north latitude
;                                      -1 : south latitude
;                 SLat       (Float) = Latitude for the grid with no
;                                      distortion (default: 70.); never set zero.
;                 delta      (Float) = Angle of rotation of the grid
;                                      to Greenwich 0 degree longitude
;                                      [default: NSIDC grid]
;
;
;OUTPUTS:         Latitude   (Float) = Geodetic Latitude
;                 Longitude  (Float) = Geodetic Longitude
;    
;KEYWORDS:        NONE    
;    
;SIDE EFFECTS:    NONE  
;
;RESTRICTIONS:    NONE    
;    
;EXAMPLE:         
;
;AUTHOR:          The original equations are from Snyder, J. P.,
;                 1982,  Map Projections Used by the U.S. Geological Survey,
;                 Geological Survey Bulletin 1532, U.S. Government Printing 
;                 Office.  See JPL Technical Memorandum 3349-85-101 for 
;                 further details.
;
;                 Original FORTRAN program written by C. S. Morris,
;                 April 1985, Jet Propulsion Laboratory, California
;                 Institute of Technology
;
;                 IDL conversion by Helmut Schottmueller, September 1995,
;                 Institute for environmental physics, University of
;                 Bremen
;
;                 Added some options: 2003-10-06 Gunnar Spreen
;-

; Conversion constant from degrees to radians
  CDR  = 57.29577951
; Standard latitude for the SSM/I grid with no distorsion
  IF NOT KEYWORD_SET(slat) THEN SLAT = 70.
; Radius of the earth in kilometers
  RE   = 6378.273
; Eccentricity of the Hughes ellipsoid squared
  E2   = .006693883
; Eccentricity of the Hughes ellipsoid
  E    =  sqrt(E2)
 
  IF (NOT VAR_EXISTS(SGN)) THEN SGN = -1

  IF NOT KEYWORD_SET(delta) THEN BEGIN
      IF (SGN EQ 1) THEN BEGIN
          delta = -45.
      ENDIF ELSE BEGIN
          delta = 0.0
      ENDELSE
  ENDIF

; slat must be positiv, correct this
  slat= ABS(slat)
  
  SL = SLAT * !PI/180.
  RHO = SQRT(X^2 + Y^2)
  CM = COS(SL) / SQRT(1.0 - E2 * (SIN(SL)^2))
  T = TAN((!PI / 4.0) - (SL / 2.0)) / ((1.0 - E * SIN(SL)) / (1.0 + E * SIN(SL)))^(E / 2.0)
  IF (ABS(SLAT-90.) LT 1.E-5) THEN BEGIN
    T = RHO * SQRT((1. + E)^(1. + E) * (1. - E)^(1. - E)) / 2. / RE
  ENDIF ELSE BEGIN
    T = RHO * T / (RE * CM)
  ENDELSE
  CHI = (!PI / 2.0) - 2.0 * ATAN(T)
  ALAT = CHI + ((E2 / 2.0) + (5.0 * E2^2.0 / 24.0) + (E2^3.0 / 12.0)) * $
         SIN(2 * CHI) + ((7.0 * E2^2.0 / 48.0) + (29.0 * E2^3 / 240.0)) * $
         SIN(4.0 * CHI) + (7.0 * E2^3.0 / 120.0) * SIN(6.0 * CHI)
  ALAT = SGN * ALAT
  ALONG = ATAN(SGN * X, -SGN * Y)
  ALONG = SGN * ALONG
  Result = WHERE(RHO LE 0.1, Count)
  IF (Count GT 0) THEN BEGIN
    ALAT(Result) = 90. * SGN
    ALONG(Result) = 0.0
  ENDIF

  aLong = aLong * 180. / !PI
  aLat  = aLat * 180. / !PI
  aLong = aLong + delta 
END
