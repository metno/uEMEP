      SUBROUTINE DISTRL(X0,Y0,X1,Y1,X2,Y2,XM,YM,DM,WM)

C The subroutine calculates the minimum distance from a given receptor
C point to a given line source.
      implicit none

      REAL X0
      REAL Y0
      REAL X1
      REAL Y1
      REAL X2
      REAL Y2
      REAL XM
      REAL YM
      REAL DM
      REAL WM
      
C X0 - Receptor point   x-coordinate
C Y0 - Receptor point   y-coordinate
C X1 - Line source      x-coordinate 1
C Y1 - Line source      y-coordinate 1
C X2 - Line source      x-coordinate 2
C Y2 - Line source      y-coordinate 2
C XM - Minimum distance x-coordinate
C YM - Minimum distance y-coordinate
C DM - Minimum distance value

C Local variables

      !REAL WM

C Calculate minimum distance (m) from the current receptor point to
C the current line source

      IF (X1 .EQ. X2 .AND. Y1 .EQ. Y2) THEN
          WM = 0.5
      ELSE
          WM = ((X0 - X1)*(X2 - X1) + (Y0 - Y1)*(Y2 - Y1))/
     .         ((X2 - X1)*(X2 - X1) + (Y2 - Y1)*(Y2 - Y1))
      ENDIF

      WM = MIN(WM,1.)
      WM = MAX(WM,0.)

      XM = (1. - WM)*X1 + WM*X2
      YM = (1. - WM)*Y1 + WM*Y2
	
      DM = SQRT((X0 - XM)*(X0 - XM) + (Y0 - YM)*(Y0 - YM))

C End of subroutine DISTRL

      END

      SUBROUTINE DISTRL_SQR(X0,Y0,X1,Y1,X2,Y2,XM,YM,DM_SQR,WM)

C The subroutine calculates the minimum distance from a given receptor
C point to a given line source.
      implicit none

      REAL X0
      REAL Y0
      REAL X1
      REAL Y1
      REAL X2
      REAL Y2
      REAL XM
      REAL YM
      REAL DM_SQR
      REAL WM
      
C X0 - Receptor point   x-coordinate
C Y0 - Receptor point   y-coordinate
C X1 - Line source      x-coordinate 1
C Y1 - Line source      y-coordinate 1
C X2 - Line source      x-coordinate 2
C Y2 - Line source      y-coordinate 2
C XM - Minimum distance x-coordinate
C YM - Minimum distance y-coordinate
C DM - Minimum distance value

C Local variables

      !REAL WM

C Calculate minimum distance (m) from the current receptor point to
C the current line source

      IF (X1 .EQ. X2 .AND. Y1 .EQ. Y2) THEN
          WM = 0.5
      ELSE
          WM = ((X0 - X1)*(X2 - X1) + (Y0 - Y1)*(Y2 - Y1))/
     .         ((X2 - X1)*(X2 - X1) + (Y2 - Y1)*(Y2 - Y1))
      ENDIF

      WM = MIN(WM,1.)
      WM = MAX(WM,0.)

      XM = (1. - WM)*X1 + WM*X2
      YM = (1. - WM)*Y1 + WM*Y2
	
      DM_SQR = (X0 - XM)*(X0 - XM) + (Y0 - YM)*(Y0 - YM)

C End of subroutine DISTRL_SQR

      END