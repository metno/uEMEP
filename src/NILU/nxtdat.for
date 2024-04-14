      SUBROUTINE NXTDAT(UN,LEOF)

C The subroutine prepares for reading the next uncommented 
C line of data from file.
C Adapted by Bruce Rolstad Denby to include '{' and '#' to skip

C Scalar arguments

      INTEGER UN
      CHARACTER*256 TXTSTR
      LOGICAL LEOF

C UN     - Fileunit
C TXTSTR - Textstring
C LEOF   - If end of file then true else false 

C If fileunit is nonpositive then just return

      IF (UN .LE. 0) RETURN

C Fileunit is positive

      LEOF = .FALSE.

C Read lines from file

  100 CONTINUE
      READ (UN,1000,END=999) TXTSTR
      IF (TXTSTR(1:1).EQ.'*'.OR.TXTSTR(1:1).EQ.'{') GOTO 100
      IF (TXTSTR(1:1).EQ.'#') GOTO 100
      BACKSPACE(UN)

      RETURN

  999 CONTINUE

      LEOF = .TRUE.
      RETURN

 1000 FORMAT (A256)

C End of subroutine NXTDAT

      END
