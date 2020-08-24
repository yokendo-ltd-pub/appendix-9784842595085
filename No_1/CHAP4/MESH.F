***************************************************************
      PROGRAM MESH
***************************************************************
      INCLUDE'MESH.H'     
      DIMENSION XZM(10,2,3),XZC(4,2,3),XZD(5,2,3)

      CALL INPUT (FREQ,OMG,TURN,XZM,NRMTL,NRDEF,XZC,NRCOIL,XZD,
     *            MODEL,NK,NL,NM,NKL,NEKL,YC,ZC,SMAX,ITUBE)
      IF(MODEL.EQ.1) CALL MESH1 (NK,NL,NM,NKL,NEKL)  
      IF(MODEL.EQ.2) CALL MESH2 (NK,NL,NM,NKL,NEKL)     

      WRITE(6,610)NK,NL,NM
      WRITE(6,620)NPOIN
      WRITE(6,630)NELEM
      WRITE(6,*)' '
  610 FORMAT(7X,'NK=',I2,'  NL=',I2,'  NM=',I2)
  620 FORMAT(7X,'NUMBER OF NODE=',I5)
  630 FORMAT(7X,'NUMBER OF ELEMENT=',I5)
  640 FORMAT(7X,'NUMBER OF UNKNOWN VALUE=',I5)

      CALL MAKEDAT(XZM,NRMTL,XZC,NRCOIL,XZD,NRDEF,
     *              NK,NL,NM,NKL,NEKL,FREQ,OMG,RL,ZL,A0,Z0)
     
      IF(ITUBE.NE.0) CALL TUBE (YC,ZC,SMAX,-RL-A0,RL+A0,Z0,Z0+ZL)

      CALL BOUND  (NK,NL,NM,NKL,MODEL)
      CALL OUTPUT (NK,NL,NM,NKL,NEKL,FREQ,OMG,TURN,
     *             MODEL,NEQ,RL,ZL,A0,Z0)

      STOP
      END

***************************************************************
      SUBROUTINE INPUT (FREQ,OMG,TURN,XZM,NRMTL,NRDEF,XZC,NRCOIL,XZD,
     *                  MODEL,NK,NL,NM,NKL,NEKL,YC,ZC,SMAX,ITUBE)
*************************************************          
      INCLUDE'MESH.H'
      DIMENSION XZM(10,2,3),XZC(4,2,3),XZD(5,2,3)
      CHARACTER*20 DUMY,FDAT1,FDAT2
            
      OPEN (1,FILE='MESH.FIL',STATUS='UNKNOWN')
         READ(1,999) DUMY
         READ(1,998) FDAT1
         READ(1,999) DUMY
         READ(1,998) FDAT2
         WRITE(6,997) FDAT1
         WRITE(6,997) FDAT2
      CLOSE(1,STATUS='KEEP')
  997 FORMAT(7X,'DATA FILE NAME = ',A20)
  998 FORMAT(A20)

      OPEN(50,FILE=FDAT1,STATUS='UNKNOWN')

         READ(50,999) DUMY
         READ(50,*)NTT,YY
         DO 10 I=1,NTT
            READ(50,*) NNT(I),YT(I)
  10     CONTINUE

         READ(50,999) DUMY
         READ(50,*)NUU,ZZ
         DO 20 I=1,NUU
            READ(50,*) NNU(I),ZU(I)
  20    CONTINUE
         READ(50,999) DUMY
         READ(50,*)XYL,ZL1,ZL2
  999    FORMAT(A50)
      CLOSE(50)

      OPEN( 3,FILE=FDAT2,STATUS='UNKNOWN')
      READ(3,99) DUMY
         READ(3,560) FREQ,CUR1(1),TURN
      READ(3,99) DUMY
         READ(3,*) NRMTL,NRCOIL,NRDEF
            IF(NRMTL.GE.10) THEN
               WRITE(*,*)'  ERR 7'
               STOP
            END IF
            IF(NRCOIL.GE.4) THEN
               WRITE(*,*)'  ERR 8'
               STOP
            END IF
            IF(NRDEF.GE.5) THEN
               WRITE(*,*)'  ERR 9'
              STOP
            END IF
      DO 30 I=1,NRCOIL
         READ(3,99) DUMY
         READ(3,*)    XZC(I,1,1),XZC(I,2,1),
     *                  XZC(I,1,3),XZC(I,2,3)
   30 CONTINUE
      DO 40 I=1,NRMTL
         READ(3,99) DUMY
         READ(3,*)    XZM(I,1,1),XZM(I,2,1),
     *                  XZM(I,1,2),XZM(I,2,2),
     *                  XZM(I,1,3),XZM(I,2,3)
         READ(3,*)    ZM(I),SI(I)
   40 CONTINUE
      DO 50 I=1,NRDEF
      READ(3,99) DUMY
         READ(3,*)  XZD(I,1,1),XZD(I,2,1),
     *                XZD(I,1,2),XZD(I,2,2),
     *                XZD(I,1,3),XZD(I,2,3)
   50 CONTINUE
         READ(3,99) DUMY
         READ(3,*) MODEL,ITUBE,YC,ZC,SMAX

      CLOSE(3)

   99 FORMAT(A70)
  500 FORMAT(7I7)
  510 FORMAT(6E12.4)
  520 FORMAT(9I7)
  530 FORMAT(5F10.3)
  540 FORMAT(13I5)
  550 FORMAT(12I5)
  560 FORMAT(6E10.3)
  561 FORMAT(4E10.3)
  562 FORMAT(2E10.3)
  570 FORMAT(2I5,6E10.3)

      RETURN
      END

***************************************************************
      SUBROUTINE MESH2 (NK,NL,NM,NKL,NEKL)
***************************************************************
      INCLUDE'MESH.H'   
      DIMENSION NAD(1000)

      J=1
      T(J)=YY
      DO 11 K=1,NTT
      DO 11 I=1,NNT(K)
      J=J+1
   11 T(J)=T(J-1)+YT(K)
      NL=J

      DO 12 I=1,NL
   12 S(I)=T(I)
      NK=J

      J=1
      U(J)=ZZ
      DO 13 K=1,NUU
      DO 13 I=1,NNU(K)
      J=J+1
   13 U(J)=U(J-1)+ZU(K)
      NM=J

      NP=0
      NU=0
      DO 10 I=1,NM
         NS=0
         NU=NU+1
         DO 10 K=1,NK
            NT=0
            NS=NS+1
            DO 10 J=1,NL
               NT=NT+1
               NP=NP+1
               XYZ(NP,1)=S(NS)
               XYZ(NP,2)=T(NT)
               XYZ(NP,3)=U(NU)
   10 CONTINUE

      NELEM=(NK-1)*(NL-1)*(NM-1)
      NPOIN=NK*NL*NM
      NKL=NK*NL
      NEKL=(NK-1)*(NL-1)
      IF(NP.NE.NPOIN)WRITE(6,600)    
  600 FORMAT(7X,'FAILED MESH GENERATION')

      M=0
      DO 20 K=1,NM-1
      DO 20 I=1,NK-1
      DO 20 J=1,NL-1
         M=M+1
         NOD(M,1)=J+(I-1)*NL+(K-1)*(NKL)
         NOD(M,2)=NOD(M,1)+NL
         NOD(M,3)=NOD(M,2)+1
         NOD(M,4)=NOD(M,1)+1
         NOD(M,5)=NOD(M,1)+NKL
         NOD(M,6)=NOD(M,2)+NKL
         NOD(M,7)=NOD(M,3)+NKL
         NOD(M,8)=NOD(M,4)+NKL
   20 CONTINUE

      PAI= 3.1415926535898
      C1=COS(PAI/4)
      C2=COS(PAI*3/4)
      DO 1 I=1,NM
	 M=(I-1)*NKL+1
         IF(XYZ(M,3).GE.ZL1.AND.XYZ(M,3).LE.ZL2) THEN
         DO 2 J=1,NL
            IF(XYZ(J,2).GT.XYL) GOTO 3
    2    CONTINUE
    3    CONTINUE
      N11=J-1

         DO 4 J=2,N11
            M=0
         DO 5 K=1,J-1
            M=M+1
            NAD(M)=J+(K-1)*NL+(I-1)*NKL
            M=M+1
            NAD(M)=1+(J-1)*NL+K-1+(I-1)*NKL
    5    CONTINUE
            M=M+1
            NAD(M)=NL*(J-1)+J+(I-1)*NKL
         NR=J
         NN=M
         DO 6 L=1,NN
            R0=100.
            DO 7 NT=1,91
               DEG=(NT-1)*PAI/180.        
               R=SQRT((XYZ(NAD(L),1)-T(NR)*COS(DEG))**2+
     *                (XYZ(NAD(L),2)-T(NR)*SIN(DEG))**2)
                  IF(R.LE.R0) THEN
                     R0=R
                     XR=T(NR)*COS(DEG)
                     YR=T(NR)*SIN(DEG)
                  END IF
    7       CONTINUE
            DO 102 JJ=1,NM
            DO 102 II=1,NK  
            DO 102 KK=1,NL  
               IF( NAD(L) .EQ. (JJ-1)*NKL+2+(II-1)*NL) GOTO 112
               IF( NAD(L) .EQ. (JJ-1)*NKL+KK+NL) GOTO 112
  102       CONTINUE
               XYZ(NAD(L),1)=XR
               XYZ(NAD(L),2)=YR
  112       CONTINUE
    6    CONTINUE
    4 CONTINUE
      END IF
    1 CONTINUE
   88 CONTINUE

      RETURN
      END


***************************************************************
      SUBROUTINE MESH1 (NK,NL,NM,NKL,NEKL)
***************************************************************
      INCLUDE'MESH.H'   
      DIMENSION NAD(1000)

      J=1
      T(J)=YY
      DO 11 K=1,NTT
      DO 11 I=1,NNT(K)
      J=J+1
   11 T(J)=T(J-1)+YT(K)
      NL=J

      J=0
      DO 13 I=NL,1,-1
      J=J+1
   13 S(J)=-T(I)
      DO 14 I=2,NL
      J=J+1
   14 S(J)=T(I)
      NK=J

      J=1
      U(J)=ZZ
      DO 12 K=1,NUU
      DO 12 I=1,NNU(K)
      J=J+1
   12 U(J)=U(J-1)+ZU(K)
      NM=J

      NP=0
      NU=0
      DO 10 I=1,NM
         NS=0
         NU=NU+1
         DO 10 K=1,NK
         NT=0
         NS=NS+1
         DO 10 J=1,NL
            NT=NT+1
            NP=NP+1
            XYZ(NP,1)=S(NS)
            XYZ(NP,2)=T(NT)
            XYZ(NP,3)=U(NU)
   10 CONTINUE

      NELEM=(NK-1)*(NL-1)*(NM-1)
      NPOIN=NK*NL*NM
      NKL=NK*NL
      NEKL=(NK-1)*(NL-1)
      IF(NP.EQ.NPOIN)WRITE(6,600)
      WRITE(6,610)NK,NL,NM
      WRITE(6,620)NPOIN,NELEM
  600 FORMAT('   << SUCCESS GENERATING MESHS >>')
  610 FORMAT('  NK=',I2,'  NL=',I2,'  NM=',I2)
  620 FORMAT('  NPOIN=',I5,'  NELEM=',I5)

      M=0
      DO 20 K=1,NM-1
      DO 20 I=1,NK-1
      DO 20 J=1,NL-1
         M=M+1
         NOD(M,1)=J+(I-1)*NL+(K-1)*(NKL)
         NOD(M,2)=NOD(M,1)+NL
         NOD(M,3)=NOD(M,2)+1
         NOD(M,4)=NOD(M,1)+1
         NOD(M,5)=NOD(M,1)+NKL
         NOD(M,6)=NOD(M,2)+NKL
         NOD(M,7)=NOD(M,3)+NKL
         NOD(M,8)=NOD(M,4)+NKL
   20 CONTINUE

      PAI= 3.1415926535898
      DO 1 I=1,NM
         M=(I-1)*NKL+1
         IF(XYZ(M,3).GE.ZL1.AND.XYZ(M,3).LE.ZL2) THEN
         DO 2 J=1,NL
            IF(XYZ(J+M ,2).GT.XYL) GOTO 3
    2    CONTINUE
    3    CONTINUE
         N11=J-1
         N=((NK-1)/2)*NL
         DO 4 J=2,N11
            M=0
            DO 5 K=1,J-1
               M=M+1
               NAD(M)=N+J+(K-1)*NL+(I-1)*NKL
               M=M+1
               NAD(M)=N+1+(J-1)*NL+K-1+(I-1)*NKL
    5       CONTINUE
            M=M+1
            NAD(M)=N+NL*(J-1)+J+(I-1)*NKL
               IF( M .GE. I1 ) STOP
            NR=J
            NN=M
            DO 6 L=1,NN
               R0=1.E5
               DO 7 NT=1,91
                  DEG=(NT-1)*1*PAI/180.        
                  R=SQRT((XYZ(NAD(L),1)-T(NR)*COS(DEG))**2+
     *                (XYZ(NAD(L),2)-T(NR)*SIN(DEG))**2)
                  ABSR=ABS(R)
                  IF(ABSR.LT.R0) THEN
                     R0=ABSR
                     XR=T(NR)*COS(DEG)
                     YR=T(NR)*SIN(DEG)
                  END IF
    7       CONTINUE
          
            DO 8 K1=1,NM
               N11=(K1-1)*NKL
               IF(XYZ(NAD(L),3).EQ.XYZ(N11+1,3)) THEN
                  DO 9 K2=1,NKL
                  N22=N11+K2
                  IF( XYZ(NAD(L),1) .EQ.-XYZ(N22,1).AND.
     *                XYZ(NAD(L),2) .EQ. XYZ(N22,2) ) THEN
                   DO 101 JJ=1,NM
                   DO 101 II=1,NK 
                   DO 101 KK=1,NL
                     IF( N22 .EQ. (JJ-1)*NKL+2+(II-1)*NL) GOTO 111
                     IF( N22 .EQ. (JJ-1)*NKL+KK+NL*(NK-3)/2) GOTO 111
  101             CONTINUE
                     XYZ(N22,1)=-XR
                     XYZ(N22,2)=YR
  111             CONTINUE
                  END IF
    9            CONTINUE
               END IF
    8    CONTINUE
             DO 102 JJ=1,NM
             DO 102 II=1,NK
             DO 102 KK=1,NL
               IF( NAD(L) .EQ. (JJ-1)*NKL+2+(II-1)*NL) GOTO 112
               IF( NAD(L) .EQ. (JJ-1)*NKL+KK+NL*(NK+1)/2) GOTO 112
  102        CONTINUE
                        XYZ(NAD(L),1)=XR
                        XYZ(NAD(L),2)=YR
  112                CONTINUE
    6    CONTINUE
    4 CONTINUE
      END IF
    1 CONTINUE
      RETURN
      END

***************************************************************
      SUBROUTINE  MAKEDAT(XZM,NRMTL,XZC,NRCOIL,XZD,NRDEF,
     *                    NK,NL,NM,NKL,NEKL,FREQ,OMG,RL,ZL,A0,Z0)
***************************************************************
      INCLUDE'MESH.H'
      DIMENSION XZM(10,2,3),XZC(4,2,3),XZD(5,2,3)
      DIMENSION XCL(4,3),ZCL(4,3)
      DIMENSION XML(10,3),YML(10,3),ZML(10,3)
      DIMENSION XDF(5,3),YDF(5,3),ZDF(5,3)

      PAI=3.141592654
      OMG=2.*PAI*FREQ
      
      IF(1.EQ.0.) GOTO 105
      M=0
      XCL(1,1)= 10000.
      ZCL(1,1)= 10000.
      XCL(1,2)=-10000.
      ZCL(1,2)=-10000.
      DO 140 I=1,NELEM
         II=0
         DO 120 K=1,8
            N=NOD(I,K)
            Z0=XYZ(N,3)
            Z00=XYZ(NOD(I,1),3)
            R0=SQRT(XYZ(N,1)**2+XYZ(N,2)**2)
            IF(R0.GE.XZC(1,1,1).AND.R0.LE.XZC(1,2,1)) II=II+1
            IF(Z0.GE.XZC(1,1,3).AND.Z0.LE.XZC(1,2,3)) II=II+1
  120    CONTINUE
      IF(II.EQ.16) THEN
        DO 110 J=1,8
            N=NOD(I,J) 
            IF( R0 .LE.XCL(1,1) ) XCL(1,1)=R0
            IF( Z00.LE.ZCL(1,1) ) ZCL(1,1)=Z00
            IF( R0 .GE.XCL(1,2) ) XCL(1,2)=R0
  110       IF( Z0 .GE.ZCL(1,2) ) ZCL(1,2)=Z0
         M=M+1
         NNCL(M)=I
         CENT1=0.
         CENT2=0.
         DO 130 L=1,8
            N=NOD(NNCL(M),L)
            CENT1=CENT1+XYZ(N,1)
  130       CENT2=CENT2+XYZ(N,2)
         CENT1=CENT1/8
         CENT2=CENT2/8
         CRT=SQRT(CENT1**2+CENT2**2)
         CUR(M,1)=CUR1(1)*CENT2/CRT
         CUR(M,2)=-CUR1(1)*CENT1/CRT
         CUR(M,3)=0.0
      END IF
  140 CONTINUE
      NCL=M
  105 CONTINUE
         
      DO 35 I=1,NELEM
         SIGMA(I)=0.
         ZMU(I)= 4.0E-7*PAI
   35 CONTINUE
      M=0
      DO 30 K=1,NRMTL
      XML(K,1)= 10000.
      YML(K,1)= 10000.
      ZML(K,1)= 10000.
      XML(K,2)=-10000.
      YML(K,2)=-10000.
      ZML(K,2)=-10000.
      DO 30 I=1,NELEM
      II=0
      DO 40 J=1,8
         N=NOD(I,J) 
         IF(XYZ(N,1).GE.XZM(K,1,1).AND.XYZ(N,1).LE.XZM(K,2,1)) II=II+1
         IF(XYZ(N,2).GE.XZM(K,1,2).AND.XYZ(N,2).LE.XZM(K,2,2)) II=II+1
         IF(XYZ(N,3).GE.XZM(K,1,3).AND.XYZ(N,3).LE.XZM(K,2,3)) II=II+1
   40 CONTINUE

         IF(II.EQ.24) THEN
            DO 45 J=1,8
                N=NOD(I,J) 
                IF( XYZ(N,1).LE.XML(K,1) ) XML(K,1)=XYZ(N,1)
                IF( XYZ(N,2).LE.YML(K,1) ) YML(K,1)=XYZ(N,2)
                IF( XYZ(N,3).LE.ZML(K,1) ) ZML(K,1)=XYZ(N,3)
                IF( XYZ(N,1).GE.XML(K,2) ) XML(K,2)=XYZ(N,1)
                IF( XYZ(N,2).GE.YML(K,2) ) YML(K,2)=XYZ(N,2)
   45           IF( XYZ(N,3).GE.ZML(K,2) ) ZML(K,2)=XYZ(N,3)
            M=M+1
            NNMTL(M)=I
            SIGMA(I)=SI(K)
            ZMU(I)=ZM(K)*4.E-7*PAI
         END IF
   30 CONTINUE
      NMTL=M
      M1=0    
      DO 1 I=1,NMTL
            DO 3 K=1,8
               MM=NOD(NNMTL(I),K)
               DO 4 L=1,M1+1
                  IF(MM.EQ.NMTLP(L)) GOTO 6
    4          CONTINUE
               M1=M1+1
               NMTLP(M1)=MM
    6          CONTINUE
    3      CONTINUE
    1    CONTINUE 
      MTLP=M1

      M=0
      DO 50 K=1,NRDEF
      XDF(K,1)= 10000.
      YDF(K,1)= 10000.
      ZDF(K,1)= 10000.
      XDF(K,2)=-10000.
      YDF(K,2)=-10000.
      ZDF(K,2)=-10000.
      DO 50 I=1,NELEM
      II=0
      DO 60 J=1,8
         N=NOD(I,J) 
         IF(XYZ(N,1).GE.XZD(K,1,1).AND.XYZ(N,1).LE.XZD(K,2,1)) II=II+1
         IF(XYZ(N,2).GE.XZD(K,1,2).AND.XYZ(N,2).LE.XZD(K,2,2)) II=II+1
         IF(XYZ(N,3).GE.XZD(K,1,3).AND.XYZ(N,3).LE.XZD(K,2,3)) II=II+1
   60 CONTINUE
         IF(II.EQ.24) THEN
         DO 23 J=1,8
                N=NOD(I,J) 
                IF( XYZ(N,1).LE.XDF(K,1) ) XDF(K,1)=XYZ(N,1)
                IF( XYZ(N,2).LE.YDF(K,1) ) YDF(K,1)=XYZ(N,2)
                IF( XYZ(N,3).LE.ZDF(K,1) ) ZDF(K,1)=XYZ(N,3)
                IF( XYZ(N,1).GE.XDF(K,2) ) XDF(K,2)=XYZ(N,1)
                IF( XYZ(N,2).GE.YDF(K,2) ) YDF(K,2)=XYZ(N,2)
   23           IF( XYZ(N,3).GE.ZDF(K,2) ) ZDF(K,2)=XYZ(N,3)
            M=M+1
            NNDEF(M)=I
            SIGMA(I)=0.
            ZMU(I)= 4.0E-7*PAI
         END IF
   50 CONTINUE
      NDEF=M

      WRITE(*,*)'  '
      DO 80 I=1,NRCOIL
         WRITE(*,*)' ******* COIL REGION ******* '
         WRITE(*,604) XCL(I,1),XCL(I,2)
   80    WRITE(*,603) ZCL(I,1),ZCL(I,2)
      ZL=(ZCL(1,2)-ZCL(1,1))*1.E-3
      RL=(XCL(1,2)-XCL(1,1))*1.E-3
      Z0=ZCL(1,1)*1.E-3
      A0=XCL(1,1)*1.E-3
      WRITE(*,*)' FREQUENCY=',FREQ

      DO 90 I=1,NRMTL

         WRITE(*,*)' ******* METAL REGION ******* '
         WRITE(*,601) XML(I,1),XML(I,2)
         WRITE(*,602) YML(I,1),YML(I,2)
         WRITE(*,603) ZML(I,1),ZML(I,2)
   90    WRITE(*,605) ZM(I),SI(I)

      DO 100 I=1,NRDEF
         WRITE(*,*)' ******* DEFECT REGION ******* '
         WRITE(*,601) XDF(I,1),XDF(I,2)
         WRITE(*,602) YDF(I,1),YDF(I,2)
  100    WRITE(*,603) ZDF(I,1),ZDF(I,2)

  601 FORMAT('  ( X1 , X2 ) = ',' ( ',E12.5,' , ',E12.5,' )')
  602 FORMAT('  ( Y1 , Y2 ) = ',' ( ',E12.5,' , ',E12.5,' )')
  603 FORMAT('  ( Z1 , Z2 ) = ',' ( ',E12.5,' , ',E12.5,' )')
  604 FORMAT('  ( R1 , R2 ) = ',' ( ',E12.5,' , ',E12.5,' )')
  605 FORMAT('   ZMU = ',E12.5,' , ','   SIGMA = ',E12.5)
         
      CUR1(2)= 1000.
      CUR1(3)=-1000. 
      CUR1(4)= 1000.
      CUR1(5)=-1000.     
      DO 21 I=1,NCL
      DO 21 J=1,8
      M=NOD(NNCL(I),J)
      IF(XYZ(M,2).EQ.0.) THEN
         R1=XYZ(M,1)
         Z1=XYZ(M,3)
         IF(R1.LT.CUR1(2)) CUR1(2)=R1
         IF(R1.GT.CUR1(3)) CUR1(3)=R1
         IF(Z1.LT.CUR1(4)) CUR1(4)=Z1
         IF(Z1.GT.CUR1(5)) CUR1(5)=Z1 
      END IF    
   21 CONTINUE

      CUR1(6)=(CUR1(3)-CUR1(2))*(CUR1(5)-CUR1(4))
      CUR1(7)=CUR1(6)*CUR1(1)
      CUR1(8)=REAL(M)
      
      RETURN
      END

***************************************************************
      SUBROUTINE TUBE (YC,ZC,SMAX,XMIN,XMAX,ZMIN,ZMAX)
***************************************************************
      INCLUDE'MESH.H'
   
      XMIN0=XMIN*1.E3
      XMAX0=XMAX*1.E3
      ZMIN0=ZMIN*1.E3
      ZMAX0=ZMAX*1.E3

      PAI=3.141592
      STA0=SMAX*PAI/180.
      ERR=(ZMAX0-ZMIN0)/20
      DO 10 I=1,NPOIN
          IF(XYZ(I,3).GE.ZC) GOTO 15
            Y1=XYZ(I,2)-YC
            Z1=ZC-XYZ(I,3)
            RD=3.1415*Z1/2
          II=0
          IF(XYZ(I,1).GE.XMIN0-ERR.AND.XYZ(I,1).LE.XMAX0+ERR) II=II+1
          IF(XYZ(I,2).GE.XMIN0-ERR.AND.XYZ(I,2).LE.XMAX0+ERR) II=II+1
          IF(XYZ(I,3).GE.ZMIN0-ERR.AND.XYZ(I,3).LE.ZMAX0+ERR) II=II+1
          IF(II.EQ.3) GOTO 15


          IF(Y1.NE.0) THEN
             STA=Y1/Z1
          ELSE
             STA=0.
          END IF
         
          IF(STA.LE.STA0) THEN
             STA1=(STA)/2
             AA=2*Z1*SIN(STA1)
             Y2=AA*COS(STA1)
             Z2=AA*SIN(STA1)          
             XYZ(I,2)=Y2+YC
             XYZ(I,3)=-Z1+Z2+ZC
          ELSE
             STA1=(STA0)/2
             AA=2*Z1*SIN(STA1)
             Y2=AA*COS(STA1)
             Z2=AA*SIN(STA1)          
             XYZ(I,3)=-Z1+Z2+ZC
          END IF

   15    CONTINUE
   10 CONTINUE

      RETURN
      END

***************************************************************
      SUBROUTINE BOUND (NK,NL,NM,NKL,MODEL)
***************************************************************
      INCLUDE'MESH.H'

      FLX=0.0
      NBD=0
      DO 10 I=1,NM

         DO 20 J=1,NL
         DO 20 L=MODEL,4
            NBD=NBD+1
            M=J+(I-1)*NKL
            NNBD(NBD)=M+(L-1)*NPOIN
            CALL VA(M,(L-1)*NPOIN,VAP,FLX)
            VNBD(NBD)=VAP
   20    CONTINUE

         DO 30 J=1,NK
         DO 30 L=1,4
         IF(L.EQ.2) GOTO 31
            NBD=NBD+1
            M=1+(J-1)*NL+(I-1)*NKL
            NNBD(NBD)=M+(L-1)*NPOIN
            CALL VA(M,(L-1)*NPOIN,VAP,FLX)
            VNBD(NBD)=VAP
   31    CONTINUE
   30    CONTINUE   
   10 CONTINUE

      DO 2 K=1,NM
         DO 3 J=1,NK
         DO 3 L=1,4
            NBD=NBD+1
            M=J*NL+(K-1)*NKL
            NNBD(NBD)=M+(L-1)*NPOIN
            CALL VA(M,(L-1)*NPOIN,VAP,FLX)
            VNBD(NBD)=VAP
    3    CONTINUE
         DO 4 J=1,NL
         DO 4 L=1,4
            NBD=NBD+1
            M=J+(NK-1)*NL+(K-1)*NKL
            NNBD(NBD)=M+(L-1)*NPOIN
            CALL VA(M,(L-1)*NPOIN,VAP,FLX)
            VNBD(NBD)=VAP
    4    CONTINUE
    2 CONTINUE

      DO 5 I=1,NKL
      DO 5 J=1,2
      DO 5 L=1,4
            NBD=NBD+1
            M=I+(J-1)*(NM-1)*NKL
            NNBD(NBD)=M+(L-1)*NPOIN
            CALL VA(M,(L-1)*NPOIN,VAP,FLX)
            VNBD(NBD)=VAP
    5 CONTINUE


      M=0
      DO 50 I=1,NBD-1       
      DO 50 J=I+1,NBD
         IF( NNBD(I).EQ.NNBD(J) ) NNBD(J)=0
   50 CONTINUE
      DO 55 I=1,NBD
         IF( NNBD(I).EQ.0 ) M=M+1
   55 CONTINUE
      DO 60 K=1,M
      DO 60 I=1,NBD
      IF( NNBD(I).EQ.0 ) THEN
         DO 65 J=I+1,NBD
            NNBD(J-1)=NNBD(J)
            VNBD(J-1)=VNBD(J)
   65    CONTINUE
      END IF
   60 CONTINUE
      NBD=NBD-M

      RETURN
      END


***************************************************************
      SUBROUTINE VA (M,NP,VAP,FLX)
***************************************************************
      INCLUDE'MESH.H'
      
      NP1=NPOIN
      NP2=NPOIN*2
      NP3=NPOIN*3
      NP4=NPOIN*4
      VAP=0.
      IF(FLX.EQ.0.) RETURN
      IF(NP.EQ.0)     VAP= FLX * XYZ(M,2)*1.E-3/2.
      IF(NP.EQ.NPOIN) VAP=-FLX * XYZ(M,1)*1.E-3/2.

      RETURN
      END 
    
***************************************************************
      SUBROUTINE ERR (N,M,MESSG)
***************************************************************
      CHARACTER*20 MESSG

      IF(M.GT.N) THEN
         WRITE(*,900)MESSG,M
         STOP
      END IF
  900 FORMAT(A20,I5)

      RETURN
      END


***************************************************************
      SUBROUTINE OUTPUT (NK,NL,NM,NKL,NEKL,FREQ,OMG,TURN,
     *                   MODEL,NEQ,RL,ZL,A0,Z0)
***************************************************************
      INCLUDE'MESH.H'

      OPEN (60,FILE='MESH.DAT',STATUS='UNKNOWN')
      OPEN (65,FILE='DOMAIN.DAT',STATUS='UNKNOWN')

         DO 200 I=1,NPOIN
         DO 200 J=1,3
            XYZ(I,J)=XYZ(I,J)*1.E-3
  200    CONTINUE
         WRITE(60,710)NPOIN,NELEM
         IF(NELEM.GT.N1) THEN
            WRITE(*,*)' NOD ',NELEM
            STOP
         END IF
         IF(NPOIN.GT.N2) THEN
            WRITE(*,*)' XYZ ',NPOIN
            STOP
         END IF
         
         DO 110 I=1,NPOIN
            WRITE(60,720) XYZ(I,1),XYZ(I,2),XYZ(I,3)
  110    CONTINUE
         DO 120 I=1,NELEM
            WRITE(60,710) (NOD(I,J),J=1,8)
  120    CONTINUE
      CLOSE (60,STATUS='KEEP')

  710 FORMAT (10I5)
  720 FORMAT (6E10.3)
     
      WRITE(65,610) NCL,NMTL,MTLP,NDEF,NBD,NEQ,MODEL      
      CALL ERR(N4,NCL,' ERR NNCL')
      CALL ERR(N4,NDEF,' ERR NNDEF ')
      CALL ERR(N6,NMTL,' ERR NNMTL ')
      CALL ERR(N7,MTLP,' ERR NMTLP ')
      CALL ERR(N3,NBD,' ERR NNBD ')
      CALL ERR(N3,NPOIN*3+MTLP,' ERR NNPT,UNKV')
      
      WRITE(65,600) FREQ,CUR1(1),TURN,RL,ZL,A0,Z0

      DO 10 I=1,NCL,12
         WRITE(65,620) (NNCL(I+J),J=0,11)
   10 CONTINUE
      DO 15 I=1,NCL,2
         WRITE(65,600) (CUR(I+J,1),CUR(I+J,2),CUR(I+J,3),J=0,1)
   15 CONTINUE

      DO 20 I=1,NMTL,12
         WRITE(65,620) (NNMTL(I+J),J=0,11)
   20 CONTINUE      
      DO 21 I=1,MTLP,12
         WRITE(65,620) (NMTLP(I+J),J=0,11)
   21 CONTINUE
      DO 22 I=1,NELEM,6
         WRITE(65,600) (ZMU(I+J),J=0, 5)
   22 CONTINUE
      DO 23 I=1,NELEM,6
         WRITE(65,600) (SIGMA(I+J),J=0, 5)
   23 CONTINUE

      DO 30 I=1,NDEF,12
         WRITE(65,620) (NNDEF(I+J),J=0,11)
   30 CONTINUE

      DO 40 I=1,NBD,12
         WRITE(65,620) (NNBD(I+J),J=0,11)
   40 CONTINUE
      DO 45 I=1,NBD,6
         WRITE(65,600) (VNBD(I+J),J=0, 5)
   45 CONTINUE

      CLOSE (65,STATUS='KEEP')

  600 FORMAT(6E12.5)
  610 FORMAT(10I7)
  620 FORMAT(12I5)

      RETURN
      END

