      PROGRAM CALLBHMIE
      IMPLICIT NONE
C Parameters:
      INTEGER MXNANG
      PARAMETER(MXNANG=1000)
C Variables:
      INTEGER IREADEP,IREADX,J,NAN,NANG,NANG0
      REAL AJ,ANG,DANG,GSCA,PI,POL,
     &     QABS,QBACK,QEXT,QSCA,RAD,REFMED,
     &     S11,S12,S33,S34,WAVEL,X
      COMPLEX REFREL,CXEPS,S1(2*MXNANG-1),S2(2*MXNANG-1)
C***********************************************************************
C Program to interactively call Bohren-Huffman Mie theory program
C
C CALLBHMIE will interactively prompt for:
C 1. refractive index of surrounding medium
C 2. either refractive index or dielectric constant of sphere
C 3. radius of sphere
C 4. wavelength (in vacuo)
C 5. number of angles at which to calculate scattering intensities
C
C CALLBHMIE will return:
C 1. Q_ext, Q_abs, Q_sca, g, Q_back
C 2. If NANG>0, then will also return scattering matrix elements
C    S_11, S_12, S_21, S_22, S_33, S_34, S_43, S_44 and POL
C    by symmetry, S_3 = S_4 = 0, so that
c           S_13 = S_14 = S_23 = S_24 = S_31 = S_32 = S_41 = S_42 = 0
C           S_21 = S_12
C           S_22 = S_11
C           S_44 = S_33
C           S_43 = -S_34
C Adapted by B.T.Draine, Princeton Univ. Obs.
C History:
C 97.08.08 (BTD): changed sign of POL
C 01.02.16 (BTD): added IMPLICIT NONE
C                 modified input to allow direct input of x
C                 added new output file callbhmie.outS_j to compare with 
C                 Wiscombe code results
C end history
C***********************************************************************
      PI=REAL(4.D0*ATAN(1.D0))
      OPEN(UNIT=7,FILE='callbhmie.out',STATUS='UNKNOWN')
      OPEN(UNIT=8,FILE='callbhmie.outS_j',STATUS='UNKNOWN')
      WRITE(*,*)' Enter (real) refractive index of surrounding medium'
      READ(*,*)REFMED
      WRITE(*,*)' Wish to enter refr.index or epsilon? (0 or 1)'
      READ(*,*)IREADEP
 1000 IF(IREADEP.LE.0)THEN
         WRITE(*,*)' Enter complex refractive index of sphere in form',
     &      ' (a,b)   [ (0,0) to stop]'
         READ(*,*)REFREL
      ELSE
         WRITE(*,*)' Enter complex epsilon of sphere in form (a,b) ',
     &            '[(0,0) to stop]'
         READ(*,*)CXEPS
         REFREL=SQRT(CXEPS)
      ENDIF
      IF(REAL(REFREL).LE.0.)STOP
      REFREL=REFREL/REFMED
      WRITE(*,6012)REFREL
      WRITE(7,6012)REFREL
 6012 FORMAT(' Complex refractive index=',1PE10.3,' +',E10.3,'*i')
 2000 WRITE(0,*)' Enter 0 to change refractive index'
      WRITE(0,*)'       1 to input x = 2*pi*a/lambda'
      WRITE(0,*)'       2 to input a and lambda separately'
      READ(*,*)IREADX
 3000 IF(IREADX.EQ.0)THEN
         GOTO 1000
      ELSEIF(IREADX.EQ.1)THEN
         WRITE(*,*)' Enter x = 2*pi*a/lambda (0 to change refr. index)'
         READ(*,*)X
         IF(X.LE.0.)GOTO 1000
         RAD=X
         WAVEL=2.*PI*RAD/X
      ELSEIF(IREADX.EQ.2)THEN
         WRITE(*,*)' Enter radius (0 to change refr. index)'
         READ(*,*)RAD
         IF(RAD.LE.0.)GOTO 1000
         WRITE(*,*)' Enter wavelength (in vacuo)'
         READ(*,*)WAVEL
         X=2.E0*PI*RAD*REFMED/WAVEL
      ELSE
         WRITE(0,*)IREADX,' is invalid entry'
         GOTO 2000
      ENDIF
      WRITE(*,*)' Enter NANG = number of angles between 0 and 90'
      READ(*,*)NANG0
      IF(NANG0.GT.MXNANG)STOP'***Error: NANG > MXNANG'
      NANG=NANG0
      IF(NANG0.LT.2)NANG=2
      WRITE(7,6013)RAD,WAVEL,X
      WRITE(*,6013)RAD,WAVEL,X
      WRITE(8,6013)RAD,WAVEL,X
C**********
C NANG=number of angles between 0 and 90 degrees (incl. 0 and 90)
C Scattering matrix elements are calculated for 2*NANG-1 angles
C including 0, 90, and 180 degrees.
C**********
      IF(NANG.GT.1)DANG=0.5E0*PI/FLOAT(NANG-1)
      CALL BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA)
      QABS=QEXT-QSCA
      WRITE(*,6065)QEXT,QSCA,QABS,GSCA,QBACK
      WRITE(7,6065)QEXT,QSCA,QABS,GSCA,QBACK
C
C POL=degree of polarization (for incident upolarized light)
C POL > 0 when polarization is perpendicular to scattering plane
C
      IF(NANG0.GT.1)THEN
         NAN=2*NANG-1
         WRITE(*,6017)
         WRITE(7,6017)
         WRITE(8,8017)
         DO 355 J=1,NAN
            AJ=J
            S11=0.5*REAL(S1(J)*CONJG(S1(J))+S2(J)*CONJG(S2(J)))
            S12=0.5*REAL(S2(J)*CONJG(S2(J))-S1(J)*CONJG(S1(J)))
            S33=REAL(S1(J)*CONJG(S2(J)))
            S34=AIMAG(S2(J)*CONJG(S1(J)))
            POL=-S12/S11
            ANG=DANG*(AJ-1.E0)*180.E0/PI
            WRITE(7,6075)ANG,S11,S12,S33,S34,POL
            WRITE(*,6075)ANG,S11,S12,S33,S34,POL
            WRITE(8,8075)ANG,S1(J),S2(J),S11,-POL
  355    CONTINUE
      ENDIF
      GOTO 3000
 6013 FORMAT(' radius=',1PE11.4,' lambda=',E11.4,' x=',E11.4)
 6017 FORMAT(2X,'theta',7X,'S11',9X,'S12=S21',7X,'S33=S44',6X,
     &   'S34=-S43',7X,'P')
 8017 FORMAT(2X,'theta','  --------- S_1 ----------',
     &                  '  --------- S_2 ----------',
     &                  ' --  S_11 --',
     &                  '    -pol.')
 6065 FORMAT(/,'Qext=',1PE11.4,' Qsca=',E11.4,' Qabs=',E11.4,
     &' <cos>=',E11.4,/,17X,'Qbk =',E11.4)
 6075 FORMAT(1X,F6.2,2X,1PE12.5,2X,E12.5,2X,E12.5,2X,E12.5,0PF9.5)
 8075 FORMAT(1X,F6.2,1P4E13.5,1PE12.5,0PF9.4)
      END

