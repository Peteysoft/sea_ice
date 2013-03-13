PRO MapLL, X, Y, Lat, Lon, SGN
;+    
;NAME:            MapLL
;    
;PURPOSE:         This subroutine converts from geodetic latitude and
;                 longitude to Polar Stereographic (X,Y) coordinates
;                 for the polar regions.
;    
;CATEGORY:        Misc
;    
;USAGE:           MapLL, X, Y, Lat, Lon, SGN
;
;INPUTS:          X   (Float)        = Polar Stereographic X Coordinate
;                                      (will be calculated for output)
;                 Y   (Float)        = Polar Stereographic Y Coordinate
;                                      (will be calculated for output)
;                 Lat        (Float) = Geodetic Latitude  (degrees, +90 to -90)  
;                 Lon        (Float) = Geodetic Longitude (degrees, 0 to 360)
;                 SGN                = Sign of latitude
;                                      +1 : north latitude
;                                      -1 : south latitude
;
;
;OUTPUTS:         X   (Float) = Polar Stereographic X Coordinate (km)
;                 Y   (Float) = Polar Stereographic Y Coordinate (km)
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
;                 Institute for environmental physics, University of Bremen
;
;                 Error correction for RHO: 1997-11-03 by LML
;-

; Conversion constant from degrees to radians
  CDR  = 57.29577951
; Standard latitude for the SSM/I grid with no distorsion
  SLAT = 90.
; Radius of the earth in kilometers
  RE   = 6378.273
; Eccentricity of the Hughes ellipsoid squared
  E2   = .006693883
; Eccentricity of the Hughes ellipsoid
  E    =  sqrt(E2)

  IF (NOT VAR_EXISTS(SGN)) THEN SGN = -1

  IF (SGN EQ 1) THEN BEGIN
    delta = 0.0; 45
  ENDIF ELSE BEGIN
    delta = 0.0
  ENDELSE

  Latitude  = ABS(Lat) * !PI/180.
  Longitude = (Lon + delta) * !PI/180.

; Compute X and Y in grid coordinates.
  T = TAN(!PI/4-Latitude/2) / $
         ((1-E*SIN(Latitude))/(1+E*SIN(Latitude)))^(E/2)
  IF ((90-SLAT) LT 1.E-5) THEN BEGIN
    RHO = 2.*RE*T/SQRT((1.+E)^(1.+E)*(1.-E)^(1.-E))
  ENDIF ELSE BEGIN
    SL  = SLAT*!PI/180.
    TC  = TAN(!PI/4.-SL/2.)/((1.-E*SIN(SL))/(1.+E*SIN(SL)))^(E/2.)
    MC  = COS(SL)/SQRT(1.0-E2*(SIN(SL)^2))
    RHO = RE*MC*T/TC
  ENDELSE
  Y=-RHO*SGN*COS(SGN*Longitude)
  X= RHO*SGN*SIN(SGN*Longitude)
  Result = WHERE(Latitude GE !PI/2, Count)
  IF (Count GT 0) THEN BEGIN
    X(Result) = 0.0
    Y(Result) = 0.0
  ENDIF
END

