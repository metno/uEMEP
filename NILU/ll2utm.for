      SUBROUTINE LL2UTM(IUTM,ISONE_IN,LAT,LON,UTMN,UTME)

      INTEGER IUTM
      INTEGER ISONE
      INTEGER ISONE_IN
      !DOUBLE PRECISION LAT
      !DOUBLE PRECISION LON
      !DOUBLE PRECISION UTMN
      !DOUBLE PRECISION UTME
      REAL LAT
      REAL LON
      REAL UTMN
      REAL UTME

C I IUTM  - UTM coordinate system indicator
C I ISONE - UTM sone
C I LAT   - Latitude  in decimal degrees
C I LON   - Longitude in decimal degrees
C O UTMN  - UTM north-coordinate (X) (meter from equator)
C O UTME  - UTM  east-coordinate (Y) (meter from west border)

C Local variables

      DOUBLE PRECISION A
      DOUBLE PRECISION B
      DOUBLE PRECISION DEAST
      DOUBLE PRECISION E
      DOUBLE PRECISION E2
      DOUBLE PRECISION F
      DOUBLE PRECISION M
      DOUBLE PRECISION N
      DOUBLE PRECISION PI
      DOUBLE PRECISION SCALE
      DOUBLE PRECISION LATV
      DOUBLE PRECISION LONV

C A     - Big semiaxis
C B     - Intermediate value
C DEAST - East movement UTM
C E     - Intermediate value 
C E2    - Intermediate value
C F     - Flattening
C LON0  - Central meridian of UTM sone
C M     - Intermediate value
C N     - Intermediate value
C PI    - The mathematical constant Pi
C SCALE - Scale UTM
C LATV  - Scaled latitude
C LONV  - Scaled longitude

C Define constants

      IF (IUTM .EQ. 1) THEN

C UTM WGS84 EUREF 89 (AirQUIS)

        A = 6378137.0
        F = 1./298.257222101        

      ENDIF

      IF (IUTM .EQ. 2) THEN

C UTM WGS84 OLD

        A = 6378137.0
        F = 1./298.257223563

      ENDIF

      IF (IUTM .EQ. 3) THEN

C UTM ED50

        A = 6378388.0
        F = 1./297.

      ENDIF

      DEAST = 500000.
      PI    = 3.141592653589793
      SCALE = 0.9996

C Scale coordinates
      ISONE=ABS(ISONE_IN)
      LATV = LAT*PI/180.0
      LON0 = (ISONE - 30)*6. - 3.
      LONV = (LON - LON0)*PI/180.

C Calculate some intermediate quantities

      E2  = F*(2. - F)

      N   = A/DSQRT(1. - E2*DSIN(LATV)*DSIN(LATV))
      E   = DSQRT(E2*DCOS(LATV)*DCOS(LATV)/(1. - E2))
      M   = N/(1. + E*E)

      B = ((1.0 - F/2.0 + F*F/16.0 + F*F*F/32.0)*LATV -
     .     (3.0*F/4.0 - 3.0*F*F*F/128.0)*DSIN(2.0*LATV) +
     .     (15.0*F*F/64.0 + 15.0*F*F*F/128.0)*DSIN(4.0*LATV) -
     .     (35.0*F*F*F/384.0)*DSIN(6.0*LATV))*A

C Calculate UTM North coordinate

      UTMN = (B + LONV*LONV*N*DSIN(LATV)*DCOS(LATV)/2.0 +
     .        LONV*LONV*LONV*LONV*N*DSIN(LATV)*
     .        DCOS(LATV)*DCOS(LATV)*DCOS(LATV)*
     .        (5.0 - DTAN(LATV)*DTAN(LATV) + 9.0*E*E + 4.0*E*E*E*E)
     .        /24.0 + LONV*LONV*LONV*LONV*LONV*LONV*N*DSIN(LATV)*
     .        DCOS(LATV)*DCOS(LATV)*DCOS(LATV)*DCOS(LATV)*DCOS(LATV)*
     .        (61.0 - 58.0*DTAN(LATV)*DTAN(LATV) +
     .        DTAN(LATV)*DTAN(LATV)*DTAN(LATV)*DTAN(LATV))/
     .        720.0)*SCALE

      IF (LAT.LT.0) then
        UTMN=UTMN+10000000.
      ENDIF

C Calculate UTM East coordinate
       
      UTME = (LONV*N*DCOS(LATV) +
     .        LONV*LONV*LONV*N*DCOS(LATV)*DCOS(LATV)*DCOS(LATV)*
     .        (1.0 - DTAN(LATV)*DTAN(LATV) + E*E)/6.0 +
     .        LONV*LONV*LONV*LONV*LONV*N*
     .        DCOS(LATV)*DCOS(LATV)*DCOS(LATV)*DCOS(LATV)*DCOS(LATV)*
     .        (5.0 - 18.0*DTAN(LATV)*DTAN(LATV) +
     .        DTAN(LATV)*DTAN(LATV)*DTAN(LATV)*DTAN(LATV))/120.0)*
     .        SCALE + DEAST

      RETURN

C End of subroutine LL2UTM

      END
