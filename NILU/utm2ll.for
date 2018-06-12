      SUBROUTINE UTM2LL(ISONE_IN,UTMN_IN,UTME,LAT,LON)

C The subroutine converts UTM north- and east-coordinates to latitude
C and longitude

      INTEGER ISONE
      REAL UTMN,UTMN_IN,UTME,LAT,LON
      INTEGER ISONE_IN
C ISONE - UTM sone
C UTMN  - UTM north-coordinate (X) (meter from equator)
C UTME  - UTM  east-coordinate (Y) (meter from west border)
C LAT   - Latitude  in decimal degrees
C LON   - Longitude in decimal degrees

C Local variables

      DOUBLE PRECISION A,BB0,DEAST,E,E2,F,FI,LA0,M,N,PI,SCALE,X,Y

C A     - Store halvakse
C BB0   - Intermediate value
C DEAST - \stforskyvning UTM
C E     - Intermediate value 
C E2    - Intermediate value
C F     - Flattrykning
C FI    - Intermediate value
C LA0   - Tangeringsmeridian
C M     - Intermediate value
C N     - Intermediate value
C PI    - The mathematical constant Pi
C SCALE - Scale UTM
C X     - Scaled north-coordinate
C Y     - Scaled  east-coordinate

C Define constants

      A     = 6378388.
      DEAST = 500000.
      F     = 1./297.
      PI    = 3.1415927
      SCALE = 0.9996

C Adjust for Southern Hemisphere, specified by negative ISONE
      IF (ISONE_IN.LT.0) then
        UTMN=UTMN_IN-10000000.
      ELSE
        UTMN=UTMN_IN
      ENDIF
      ISONE=ABS(ISONE_IN)

C Scale coordinates

      X   = UTMN/SCALE
      Y   = (UTME - DEAST)/SCALE
      LA0 = (ISONE - 30)*6. - 3.

C Calculate some intermediate quantities

      E2  = F*(2. - F)
      BB0 = (1. - F/2. + F*F/16. + F*F*F/32.)*A

      FI  = X/BB0 +
     .      (3.*F/4. + 3.*F*F/8. + 21.*F*F*F/256.)*DSIN(2.*X/BB0) +
     .      (21.*F*F/64. + 21.*F*F*F/64.)*DSIN(4.*X/BB0) +
     .      (151.*F*F*F/768.)*DSIN(6.*X/BB0)

      N   = A/DSQRT(1. - E2*DSIN(FI)*DSIN(FI))
      E   = DSQRT(E2*DCOS(FI)*DCOS(FI)/(1. - E2))
      M   = N/(1. + E*E)

C Calculate latitude and longitude in radians

      LAT = FI - (Y*Y*DTAN(FI))/(2.*M*N) +
     .      (Y*Y*Y*Y*DTAN(FI))/(24.*M*N*N*N)*
     .      (5. + 3.*DTAN(FI)*DTAN(FI) + E*E -
     .       9.*E*E*DTAN(FI)*DTAN(FI) - 4.*E*E*E*E) -
     .      (Y*Y*Y*Y*Y*Y*DTAN(FI))/(720.*M*N*N*N*N*N)*
     .      (61 + 90*DTAN(FI)*DTAN(FI) +
     .       45*DTAN(FI)*DTAN(FI)*DTAN(FI)*DTAN(FI))

      LON = Y/(N*DCOS(FI)) -
     .      (Y*Y*Y*(1. + 2.*DTAN(FI)*DTAN(FI) + E*E))/
     .      (6.*N*N*N*DCOS(FI)) +
     .      (Y*Y*Y*Y*Y*(5. + 28.*DTAN(FI)*DTAN(FI) + 
     .      24.*DTAN(FI)*DTAN(FI)*DTAN(FI)*DTAN(FI)))/
     .      (120.*N*N*N*N*N*DCOS(FI)) + LA0*PI/180.

C Convert from radians to degrees

      LAT = LAT*180./PI
      LON = LON*180./PI

      RETURN

C End of subroutine UTM2LL

      END
