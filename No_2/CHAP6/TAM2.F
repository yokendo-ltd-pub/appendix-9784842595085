C  /*
C   (1/8) (1/4) (1/2) (1/1) SPHERE SURFACE ONLY
C   FILE NAME = LGCYYF0671/KAYO/TAM2
C  */
C MAIN()
      implicit double precision (A-H,O-Z)
      WRITE(*,600)
  600 FORMAT(1H0,'FLAT1  OR  FLAT3 ?')
      READ(*,*)ID
      IF(ID.NE.1.AND.ID.NE.3)STOP' AHO !'
      IF(ID.EQ.1) CALL FLAT1
      IF(ID.EQ.3) CALL FLAT3
      CALL JOINT
      CALL OUTPUT
      STOP
      END
C
C FLAT1()
      SUBROUTINE FLAT1
      implicit double precision (A-H,O-Z)
      COMMON /TWO/NP2,NE2,Z2(1000,3),NOE2(2000,3)
      DIMENSION NUE(40),NST(40)
      PI=3.141592654
      PRINT*,'SUB FLAT1()'
      PRINT*,'NSD ?'
      NSD=5
      READ(5,*,END=510)NSD
  510 CONTINUE
      N=0
      DO 230 I=1,NSD+1
      PHI=(0.5*PI)*FLOAT(I-1)/FLOAT(NSD)
      ND=(I-1)
      IF(ND.EQ.0) ND=1
      DO 220 J=1,I
      N=N+1
      THE=(0.5*PI)*FLOAT(J-1)/FLOAT(ND)
      Z2(N,1)=COS(THE)*SIN(PHI)
      IF(ABS(Z2(N,1)).LT.1E-7) Z2(N,1)=0.
      Z2(N,2)=SIN(THE)*SIN(PHI)
      IF(ABS(Z2(N,2)).LT.1E-7) Z2(N,2)=0.
      Z2(N,3)=COS(PHI)
      IF(ABS(Z2(N,3)).LT.1E-7) Z2(N,3)=0.
  220 CONTINUE
  230 CONTINUE
      NP2=N
      N=0
      DO 270 I=1,NSD
C
      NUE0=(I*(I-1)/2)
      DO 240 K=1,I
      NUE(K)=NUE0+K
  240 CONTINUE
      IF(I.EQ.1)NUE(1)=1
C
      NST0=(I*(I+1)/2)
      DO 250 K=1,I+1
      NST(K)=NST0+K
  250 CONTINUE
C
      N=N+1
      NOE2(N,1)=NST(1)
      NOE2(N,2)=NST(2)
      NOE2(N,3)=NUE(1)
      IF(I.EQ.1) GO TO 270
      DO 260 L=1,I-1
      N=N+1
      NOE2(N,1)=NST(L+1)
      NOE2(N,2)=NUE(L+1)
      NOE2(N,3)=NUE(L)
      N=N+1
      NOE2(N,1)=NST(L+1)
      NOE2(N,2)=NST(L+2)
      NOE2(N,3)=NUE(L+1)
  260 CONTINUE
  270 CONTINUE
      NE2=N
      RETURN
      END
C
C FLAT3()
      SUBROUTINE FLAT3
      implicit double precision (A-H,O-Z)
      COMMON /TWO/NP2,NE2,Z2(1000,3),NOE2(2000,3)
      DIMENSION NUE(40),NST(40),TH(40),PH(40)
      PI=3.141592654
      PRINT*,'SUB FLAT3()'
      PRINT*,'ND TH(AROUND X/Y-AXIS) ?'
      READ(*,*)NTH
      DO 110 I=1,NTH+1
      TH(I)=(90.)*FLOAT(I-1)/FLOAT(NTH)
  110 CONTINUE
      PRINT*,'ND PHI(AROUND Z-AXIS) ?'
      READ(*,*)NPH
      DO 120 I=1,NPH+1
      PH(I)=(90.)*FLOAT(I-1)/FLOAT(NPH)
  120 CONTINUE
      N=0
      DO 230 I=1,NTH+1
      THE=(PI/180.)*TH(I)
      ND=NPH+1
      IF(I.EQ.1) ND=1
      DO 220 J=1,ND
      N=N+1
      PHI=(PI/180.)*PH(J)
      Z2(N,1)=COS(PHI)*SIN(THE)
      Z2(N,2)=SIN(PHI)*SIN(THE)
      Z2(N,3)=COS(THE)
  220 CONTINUE
  230 CONTINUE
      NP2=N
      N=0
      DO 240 K=1,NPH+1
  240 NST(K)=1
      DO 270 I=1,NTH
      DO 250 K=1,NPH+1
  250 NUE(K)=NST(K)
      NST0=1+(NPH+1)*(I-1)
      DO 260 K=1,NPH+1
  260 NST(K)=NST0+K
      DO 270 L=1,NPH
      N=N+1
      NOE2(N,1)=NST(L)
      NOE2(N,2)=NST(L+1)
      NOE2(N,3)=NUE(L)
      IF(I.EQ.1)GOTO 270
      N=N+1
      NOE2(N,1)=NST(L+1)
      NOE2(N,2)=NUE(L+1)
      NOE2(N,3)=NUE(L)
  270 CONTINUE
      NE2=N
      RETURN
      END
C
C OUTPUT()
      SUBROUTINE OUTPUT
      implicit double precision (A-H,O-Z)
      COMMON /TWO/NP2,NE2,Z2(1000,3),NOE2(2000,3)
      COMMON /THR/NP3,NE3,Z3(2000,3),NOE3(3000,3)
      CHARACTER*20 FNAM
      PRINT*,'RADIUS ?'
      READ(*,*)RMAX
      WRITE(*,*)'FILE17='
      FNAM='f17'
      READ(*,'(A20)',END=111)
  111 OPEN(07,FILE=FNAM,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     &                                      FORM='FORMATTED')
      WRITE(07,1)NP3,NE3
      WRITE(07,2)((RMAX*Z3(I,J),J=1,3),I=1,NP3)
      WRITE(07,1)((NOE3(I,J),J=1,3),I=1,NE3)
      CLOSE(07)
      WRITE(*,610)7,NP3,NE3
    1 FORMAT(I6)
    2 FORMAT(E15.8)
  610 FORMAT(1H0,A4,'NP =',I5,' , NE =',I5)
      RETURN
      END
C
C JOINT()
      SUBROUTINE JOINT
      implicit double precision (A-H,O-Z)
      COMMON /TWO/NP2,NE2,Z2(1000,3),NOE2(2000,3)
      COMMON /THR/NP3,NE3,Z3(2000,3),NOE3(3000,3)
      DIMENSION ICHK(2000)
      PRINT*,'X=0 SYM PLANE ? (1/ANY)'
      ISX=0
      READ(5,*,END=110)ISX
  110 PRINT*,'Y=0 SYM PLANE ? (1/ANY)'
      ISY=0
      READ(5,*,END=120)ISY
  120 PRINT*,'Z=0 SYM PLANE ? (1/ANY)'
      ISZ=0
      READ(5,*,END=130)ISZ
  130 CONTINUE
      NP3=0
      NE3=0
      LAD=0
      DO 190 IX=1,ISX+1
      DO 190 IY=1,ISY+1
      DO 190 IZ=1,ISZ+1
      LAD=LAD+1
      DO 140 II=1,NP2
      NP3=NP3+1
      Z3(NP3,1)=Z2(II,1)*(-1.)**(IX-1)
      Z3(NP3,2)=Z2(II,2)*(-1.)**(IY-1)
      Z3(NP3,3)=Z2(II,3)*(-1.)**(IZ-1)
  140 CONTINUE
      NREV=IX-1+IY-1+IZ-1
      NREV=MOD(NREV,2)
      IF(NREV.EQ.1)GO TO 160
      DO 150 MM=1,NE2
      NE3=NE3+1
      DO 150 J=1,3
      NOE3(NE3,J)=NOE2(MM,J)+(LAD-1)*NP2
  150 CONTINUE
      GO TO 190
  160 CONTINUE
      DO 170 MM=1,NE2
      NE3=NE3+1
      NOE3(NE3,1)=NOE2(MM,1)+(LAD-1)*NP2
      NOE3(NE3,2)=NOE2(MM,3)+(LAD-1)*NP2
      NOE3(NE3,3)=NOE2(MM,2)+(LAD-1)*NP2
  170 CONTINUE
  190 CONTINUE
C
      DO 210 I=1,NP3
  210 ICHK(I)=I
      DO 250 I1=1,NP3-1
      IF(ICHK(I1).NE.I1)GO TO 250
      DO 240 I2=I1+1,NP3
      R=SQRT((Z3(I1,1)-Z3(I2,1))**2
     &      +(Z3(I1,2)-Z3(I2,2))**2
     &      +(Z3(I1,3)-Z3(I2,3))**2)
      IF(R.GT.1.E-06)GO TO 240
      ICHK(I2)=I1
  240 CONTINUE
  250 CONTINUE
C
      N=0
      DO 270 II=1,NP3
      IF(ICHK(II).EQ.II)GO TO 255
      ICHK(II)=ICHK(ICHK(II))
      GO TO 270
  255 N=N+1
      DO 260 J=1,3
  260 Z3(N,J)=Z3(II,J)
      ICHK(II)=N
  270 CONTINUE
      NP3=N
C
      DO 280 MM=1,NE3
      DO 280 J=1,3
      NOE3(MM,J)=ICHK(NOE3(MM,J))
  280 CONTINUE
      RETURN
      END

