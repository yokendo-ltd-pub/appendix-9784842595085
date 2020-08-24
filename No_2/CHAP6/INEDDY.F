      PROGRAM MAIN
      PARAMETER (NPPN=1000,NEEN=2000,NVVN=3000,IFMX=2)
      implicit double precision ( A-H, O-Z)
      COMPLEX*16 FEPS1,FEPS2,SA,SB,SC
      COMPLEX*16 ALP,JO,JOME1,JOME2,JK1,JK2
      COMPLEX*16 BI,EI,CCE
      COMMON /TR1/NP,NE3,ZZ(NPPN,3),NOE3(NEEN,6),ZNP(NEEN,3)
      COMMON /VCT/ZLN(NVVN,6),NOV(NVVN,2),NV
      COMMON /TR4/BN(NEEN,3),T1(NEEN,3),T2(NEEN,3),ZI(NEEN,3)
      COMMON /KNOWN/BI(NVVN+NVVN),EI(NVVN+NVVN)
      COMMON /UNKNOWN1/NPI,IO(1000)
      COMMON /UNKNOWN2/ZZZ(1000,3)
      COMMON /E0B0/B0(1000,3),A0(1000,3),CCE(1000,6)
      COMMON /TR6/SA(6,3*(NVVN+NVVN)),SB(6,3*(NVVN+NVVN)),
     &                                SC(6,3*(NVVN+NVVN))
      COMMON /MTE0/AMU0,AMU1,AMU2,BETA,EPS0,EPS1,EPS2,OMG
      COMMON /MTE1/SGM1,SGM2,FEPS1,FEPS2,ALP
      COMMON /MTE2/JO,JOME1,JOME2,JK1,JK2
C
      CALL FILEIN
C
      FEPS1=CMPLX(EPS1,-SGM1/OMG)
      FEPS2=CMPLX(EPS2,-SGM2/OMG)
      JO=CMPLX(0.,OMG)
      ALP=FEPS2/FEPS1
      BETA=AMU2/AMU1
      JOME1=JO*AMU1*FEPS1
      JOME2=JO*AMU2*FEPS2
      JK1=JO*ZSQRT(AMU1*FEPS1)
      JK2=JO*ZSQRT(AMU2*FEPS2)
      WRITE(*,*)'SGM1,SGM2,AMU1,AMU2=',SGM1,SGM2,AMU1,AMU2
      WRITE(*,*)'ALP,BETA=',ALP,BETA
      WRITE(*,*)'OMG=',OMG,' JOME1=',JOME1
      CALL FLAT
C Calculate Integral Point
C               +---> R<0.0001(m) Triangle Division
C               | +---> R=0(m) Triangle Division
      CALL SGAL(2,3)
      CALL CEINT(1)
      CALL VCTPREP
      CALL MAT
      CALL OUTDATA
C
      END
C
      SUBROUTINE FLAT
      PARAMETER (NPPN=1000,NEEN=2000,NVVN=3000,IFMX=2)
      implicit double precision ( A-H, O-Z)
      COMMON /TR1/NP,NE3,ZZ(NPPN,3),NOE3(NEEN,6),ZNP(NEEN,3)
      COMMON /VCT/ZLN(NVVN,6),NOV(NVVN,2),NV
      COMMON /TR4/BN(NEEN,3),T1(NEEN,3),T2(NEEN,3),ZI(NEEN,3)
      COMMON /HSI/HOSEI
C
      SS=0.
      DO 10 I=1,NE3
         I1=NOE3(I,1)
         I2=NOE3(I,2)
         I3=NOE3(I,3)
         X1=ZZ(I1,1)
         X2=ZZ(I2,1)
         X3=ZZ(I3,1)
         Y1=ZZ(I1,2)
         Y2=ZZ(I2,2)
         Y3=ZZ(I3,2)
         Z1=ZZ(I1,3)
         Z2=ZZ(I2,3)
         Z3=ZZ(I3,3)
         RA=SQRT((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2)
         RB=SQRT((X3-X1)**2+(Y3-Y1)**2+(Z3-Z1)**2)
         RC=SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
C Inner Center Point
         ZI(I,1)=(RA*X1+RB*X2+RC*X3)/(RA+RB+RC)
         ZI(I,2)=(RA*Y1+RB*Y2+RC*Y3)/(RA+RB+RC)
         ZI(I,3)=(RA*Z1+RB*Z2+RC*Z3)/(RA+RB+RC)
         AX=X2-X1
         AY=Y2-Y1
         AZ=Z2-Z1
         BX=X3-X2
         BY=Y3-Y2
         BZ=Z3-Z2
         BNX=AY*BZ-AZ*BY
         BNY=AZ*BX-AX*BZ
         BNZ=AX*BY-AY*BX
         BB=SQRT(BNX*BNX+BNY*BNY+BNZ*BNZ)
         SS=SS+BB
         BNX=BNX/BB
         BNY=BNY/BB
         BNZ=BNZ/BB
         IF(ABS(BNX).GT.ABS(BNY).AND.ABS(BNX).GT.ABS(BNZ))THEN
            T1X=-(BNY+BNZ)/BNX
            T1Y=1.
            T1Z=1.
         ELSEIF(ABS(BNY).GT.ABS(BNX).AND.ABS(BNY).GT.ABS(BNZ))THEN
            T1X=1.
            T1Y=-(BNX+BNZ)/BNY
            T1Z=1.
         ELSE
            T1X=1.
            T1Y=1.
            T1Z=-(BNX+BNY)/BNZ
         ENDIF
         TT1=SQRT(T1X*T1X+T1Y*T1Y+T1Z*T1Z)
         T1X=T1X/TT1
         T1Y=T1Y/TT1
         T1Z=T1Z/TT1
         T2X=BNY*T1Z-BNZ*T1Y
         T2Y=BNZ*T1X-BNX*T1Z
         T2Z=BNX*T1Y-BNY*T1X
         BN(I,1)=BNX
         BN(I,2)=BNY
         BN(I,3)=BNZ
         T1(I,1)=T1X
         T1(I,2)=T1Y
         T1(I,3)=T1Z
         T2(I,1)=T2X
         T2(I,2)=T2Y
         T2(I,3)=T2Z
   10 CONTINUE
C
C for the vector element
C
      DO 20 I=1,NV
         I1=NOV(I,1)
         I2=NOV(I,2)
         DO 30 J=1,3
            ZLN(I,J+3)=(ZZ(I2,J)+ZZ(I1,J))*.5
            ZLN(I,J  )=(ZZ(I2,J)-ZZ(I1,J))*.5
   30    CONTINUE
         AL=SQRT(ZLN(I,1)**2+ZLN(I,2)**2+ZLN(I,3)**2)
         DO 40 J=1,3
            ZLN(I,J)=ZLN(I,J)/AL
   40    CONTINUE
   20 CONTINUE
      HOSEI=(3.14159265E-2)/(SS*.5)
      RETURN
      END
C
C Read File
      SUBROUTINE FILEIN
      PARAMETER (NPPN=1000,NEEN=2000,NVVN=3000,IFMX=2)
      implicit double precision ( A-H, O-Z)
      COMPLEX*16 FEPS1,FEPS2
      COMPLEX*16 ALP,JO,JOME1,JOME2,JK1,JK2
      COMPLEX*16 BI,EI,CCE
      COMMON /TR1/NP,NE3,ZZ(NPPN,3),NOE3(NEEN,6),ZNP(NEEN,3)
      COMMON /KNOWN/BI(NVVN+NVVN),EI(NVVN+NVVN)
      COMMON /UNKNOWN1/NPI,IO(1000)
      COMMON /UNKNOWN2/ZZZ(1000,3)
      COMMON /MTE0/AMU0,AMU1,AMU2,BETA,EPS0,EPS1,EPS2,OMG
      COMMON /MTE1/SGM1,SGM2,FEPS1,FEPS2,ALP
      COMMON /MTE2/JO,JOME1,JOME2,JK1,JK2
      COMMON /E0B0/B0(1000,3),A0(1000,3),CCE(1000,6)
      COMMON /VCT/ZLN(NVVN,6),NOV(NVVN,2),NV
      COMMON /LIN/NOL(NVVN),NUML,NUMOLD,NUMXYZ
      COMMON /FNM/FN1,FN2,FN3,FN4,FN5,FN6,FN7,FN8,FN9,FN10,FN11,FN12
      CHARACTER*20 FN,FN1,FN2,FN3,FN4,FN5,FN6,
     &                FN7,FN8,FN9,FN10,FN11,FN12,FN13
      CHARACTER*4 FN0
C
      WRITE(*,*)'FIL='
      FN='FFF'
      READ(*,'(A20)',END=111)FN
  111 OPEN(08,FILE=FN,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      READ(08,'(A20)')FN1
      READ(08,'(A20)')FN2
      READ(08,'(A20)')FN3
      READ(08,'(A4)' )FN0
      READ(08,'(A20)')FN5
      READ(08,'(A20)')FN6
      READ(08,'(A20)')FN7
      READ(08,'(A20)')FN8
      READ(08,'(A20)')FN9
      READ(08,'(A20)')FN10
      READ(08,'(A20)')FN11
      READ(08,'(A20)')FN12
      READ(08,'(A20)')FN13
      CLOSE(08)
      WRITE(*,*)'Triangular & Vector Element       =',FN1
      WRITE(*,*)'B,E on the Boundary Surface       =  NTTBE'
      WRITE(*,*)'( x, y, z ) at the arbitary point =',FN7
      WRITE(*,*)'B(impressed)                      =',FN8
      WRITE(*,*)'E(impressed)                      =',FN9
      WRITE(*,*)'E(Calculate)                      =',FN10
      WRITE(*,*)'B(Calculate)                      =',FN11
C Form File
      OPEN(17,FILE=FN1,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      READ(17,3)NP,NE3
      READ(17,4)((ZZ(I,J),J=1,3),I=1,NP)
      READ(17,3)((NOE3(I,J),J=1,3),I=1,NE3)
C for the vector element
      READ(17,1)((NOE3(I,J),J=4,6),I=1,NE3)
      READ(17,1)NV
      READ(17,1)((NOV(I,J),J=1,2),I=1,NV)
      CLOSE(17)
C
      OPEN(37,FILE=FN7,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      READ(37,1)NPI
      READ(37,2)((ZZZ(I,J),J=1,3),I=1,NPI)
      READ(37,1)(IO(I),I=1,NPI)
      CLOSE(37)
      OPEN(24,FILE='NTTBE',ACCESS='SEQUENTIAL',FORM='FORMATTED')
      READ(24,2)(EI(I),I=1,NV)
      READ(24,2)(BI(I),I=1,NV)
      CLOSE(24)
      OPEN(49,FILE=FN8,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      READ(49,1) NN
      READ(49,2)(X,Y,Z,I=1,NN)
      READ(49,2)((B0(I,J),J=1,3),I=1,NN)
      CLOSE(49)
      OPEN(50,FILE=FN9,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      READ(50,1) NN
      READ(50,2)(X,Y,Z,I=1,NN)
      READ(50,2)((A0(I,J),J=1,3),I=1,NN)
      CLOSE(50)
C Frequency, SGM1, SGM2
C
      PI=3.141592654
      EPS0=1./36.E+9/PI
      AMU0=4.E-7*PI
C
      OPEN(41,FILE=FN13,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      READ(41,2)F
      READ(41,2)AMU1
      READ(41,2)EPS1
      READ(41,2)SGM1
      READ(41,2)AMU2
      READ(41,2)EPS2
      READ(41,2)SGM2
      CLOSE(41)
      OMG=2.*PI*F
      AMU1=AMU1*AMU0
      AMU2=AMU2*AMU0
      EPS1=EPS1*EPS0
      EPS2=EPS2*EPS0
C
    1 FORMAT(I6)
    2 FORMAT(E15.8)
    3 FORMAT(3I5)
    4 FORMAT(3F10.3)
      RETURN
      END
C
C Triangle Division
      SUBROUTINE SGAL(NN,NHK)
      PARAMETER (NPPN=1000,NEEN=2000,NVVN=3000,IFMX=2)
      implicit double precision ( A-H, O-Z)
      COMMON /TR1/NP,NE3,ZZ(NPPN,3),NOE3(NEEN,6),ZNP(NEEN,3)
      COMMON /TRS1/SAL1(500),SAL2(500),SWWW(500)
      COMMON /TRS2/NGS
      COMMON /TRSB1/SBL1(1,500),SBL2(1,500),SBWW(1,500)
      COMMON /TRSB2/NGBS
      COMMON /TRG/ALS(500,4),ALB(1,3,2)
      COMMON /TAL/ALZN(NEEN,3)
      COMMON /INTEG/NIG
      DIMENSION EL1(6),EL2(6),Z(20,2),NOE(20,3),NNA(20,20)
      DIMENSION EK1(6),EK2(6)
      DIMENSION NMG(4,3)
C
      NMG(1,1)=1 
      NMG(1,2)=4 
      NMG(1,3)=6
      NMG(2,1)=4 
      NMG(2,2)=2 
      NMG(2,3)=5
      NMG(3,1)=6 
      NMG(3,2)=5 
      NMG(3,3)=3
      NMG(4,1)=4 
      NMG(4,2)=5 
      NMG(4,3)=6
      NIG=7
      G1=0.79742699
      G2=0.10128651
      G3=0.05971587
      G4=0.47014206
      G5=0.33333333
      W1=0.12593918
      W2=0.13239415
      W3=0.22500000
      ALS(1,1)=G1 
      ALS(1,2)=G2 
      ALS(1,3)=G2 
      ALS(1,4)=W1
      ALS(2,1)=G2 
      ALS(2,2)=G1 
      ALS(2,3)=G2 
      ALS(2,4)=W1
      ALS(3,1)=G2 
      ALS(3,2)=G2 
      ALS(3,3)=G1 
      ALS(3,4)=W1
      ALS(4,1)=G3 
      ALS(4,2)=G4 
      ALS(4,3)=G4 
      ALS(4,4)=W2
      ALS(5,1)=G4 
      ALS(5,2)=G3 
      ALS(5,3)=G4 
      ALS(5,4)=W2
      ALS(6,1)=G4 
      ALS(6,2)=G4 
      ALS(6,3)=G3 
      ALS(6,4)=W2
      ALS(7,1)=G5 
      ALS(7,2)=G5 
      ALS(7,3)=G5 
      ALS(7,4)=W3
      DO 10 M=1,NN+1
         DO 10 N=1,NN-M+2
   10       NNA(M,N)=(N+M-1)*(N+M-2)/2+M
      DO 30 M=1,NN
         DO 30 N=1,NN-M+1
            IA=(N+M-2)**2+2*M-1
            NOE(IA,1)=NNA(M,N)
            NOE(IA,2)=NNA(M,N+1)
   30       NOE(IA,3)=NNA(M,N+1)+1
      DO 50 M=1,NN-1
         DO 50 N=1,NN-M
            IB=(N+M-1)**2+2*M
            NOE(IB,1)=NNA(M+1,N)
            NOE(IB,2)=NNA(M+1,N)-1
   50       NOE(IB,3)=NNA(M+1,N+1)
      NE=NN*NN
      DO 70 N=2,NN+1
         NA=N*(N+1)/2
         NB=N*(N-1)/2+1
         Z(NB,1)=FLOAT(N-1)/FLOAT(NN)
         Z(NB,2)=0.
         Z(NA,1)=-Z(NB,2)
         Z(NA,2)= Z(NB,1)
         DO 70 I=NB+1,NA-1
            Z(I,1)=((NA-I)*Z(NB,1)+(I-NB)*Z(NA,1))/(NA-NB)
   70       Z(I,2)=((NA-I)*Z(NB,2)+(I-NB)*Z(NA,2))/(NA-NB)
      NGS=0
      DO 100 I=1,NE
         DO 110 J=1,3
            EL1(J)=Z(NOE(I,J),1)
  110       EL2(J)=Z(NOE(I,J),2)
         DO 100 IG=1,NIG
            NGS=NGS+1
            SAL1(NGS)=ALS(IG,1)*EL1(1)+ALS(IG,2)*EL1(2)+ALS(IG,3)*EL1(3)
            SAL2(NGS)=ALS(IG,1)*EL2(1)+ALS(IG,2)*EL2(2)+ALS(IG,3)*EL2(3)
  100       SWWW(NGS)=ALS(IG,4)/FLOAT(NE)
      WRITE(*,*)'......R<0.0001(m) Rediv. Tri. ',NGS/NIG
C Triangle in case of Source Point near Field Point
C
C      DO 956 JK=1,NE3
      JK=1
         N=0
         SSS=0
         EL1(1)=1.
         EL2(1)=0.
         EL1(2)=0.
         EL2(2)=1.
         EL1(3)=0.
         EL2(3)=0.
         DO 905 IQ=1,NHK
            IF(IQ.NE.1)THEN
               DO 225 I=1,3
                  EL1(I)=EL1(NMG(IAX,I))
  225             EL2(I)=EL2(NMG(IAX,I))
            ENDIF
            EL1(4)=(EL1(1)+EL1(2))/2.
            EL2(4)=(EL2(1)+EL2(2))/2.
            EL1(5)=(EL1(2)+EL1(3))/2.
            EL2(5)=(EL2(2)+EL2(3))/2.
            EL1(6)=(EL1(3)+EL1(1))/2.
            EL2(6)=(EL2(3)+EL2(1))/2.
            DLS=0.25**(IQ+1)
            DO 305 IA=1,4
               EK1(1)=EL1(NMG(IA,1))
               EK2(1)=EL2(NMG(IA,1))
               EK1(2)=EL1(NMG(IA,2))
               EK2(2)=EL2(NMG(IA,2))
               EK1(3)=EL1(NMG(IA,3))
               EK2(3)=EL2(NMG(IA,3))
               EK1(4)=0.5*(EK1(1)+EK1(2))
               EK2(4)=0.5*(EK2(1)+EK2(2))
               EK1(5)=0.5*(EK1(2)+EK1(3))
               EK2(5)=0.5*(EK2(2)+EK2(3))
               EK1(6)=0.5*(EK1(3)+EK1(1))
               EK2(6)=0.5*(EK2(3)+EK2(1))
               AX=1./3.
               AY=1./3.
               X1=EK1(1)
               X2=EK1(2)
               X3=EK1(3)
               Y1=EK2(1)
               Y2=EK2(2)
               Y3=EK2(3)
               SS=ABS((X1-X2)*(Y2-Y3)+(Y1-Y2)*(X3-X2))
               S1=ABS((X1-X2)*(Y2-AY)+(Y1-Y2)*(AX-X2))
               S2=ABS((X1-X3)*(Y3-AY)+(Y1-Y3)*(AX-X3))
               S3=ABS((X2-X3)*(Y3-AY)+(Y2-Y3)*(AX-X3))
               SA=ABS(S1+S2+S3-SS)
               IF(SA.LT.1.E-7) THEN
                  IAX=IA
               ELSE
                  DO 307 IB=1,4
                     DO 307 IG=1,NIG
                        N=N+1
                        SBL1(JK,N)=ALS(IG,1)*EK1(NMG(IB,1))+ALS(IG,2)
     &                  *EK1(NMG(IB,2))+ALS(IG,3)*EK1(NMG(IB,3))
                        SBL2(JK,N)=ALS(IG,1)*EK2(NMG(IB,1))+ALS(IG,2)
     &                  *EK2(NMG(IB,2))+ALS(IG,3)*EK2(NMG(IB,3))
  307                   SBWW(JK,N)=ALS(IG,4)*DLS
               ENDIF
  305       CONTINUE
  905    CONTINUE
         NGBS=N
C Small Triangle around Null Ponit in Triangle in case of which
C    Field Point equal to Source Point
         ALB(JK,1,1)=EL1(NMG(IAX,1))
         ALB(JK,1,2)=EL2(NMG(IAX,1))
         ALB(JK,2,1)=EL1(NMG(IAX,2))
         ALB(JK,2,2)=EL2(NMG(IAX,2))
         ALB(JK,3,1)=EL1(NMG(IAX,3))
         ALB(JK,3,2)=EL2(NMG(IAX,3))
C Null Point in Triangle in case of which Field Point equal Source Point
C 956 CONTINUE
      WRITE(*,*)'......R=0     (m) Rediv. Tri. ',NGBS/NIG
      RETURN
      END
C
      SUBROUTINE CEINT(NDS)
      PARAMETER (NPPN=1000,NEEN=2000,NVVN=3000,IFMX=2)
      implicit double precision ( A-H, O-Z)
      COMMON /ITGCE/SL1(3,500),SL2(3,500),SL3(3,500),SWW(3,500),NCE
C
      COMMON /TT0/ALGS(7,4),NOE(500,3),ZZ(800,2)
C
      ALGS(1,1)=0.333333E0
      ALGS(1,2)=0.333333E0
      ALGS(1,3)=0.333333E0
      ALGS(1,4)=0.225000E0
      ALGS(2,1)=0.797426E0
      ALGS(2,2)=0.101286E0
      ALGS(2,3)=0.101286E0
      ALGS(2,4)=0.125929E0
      ALGS(3,1)=ALGS(2,3)
      ALGS(3,2)=ALGS(2,1)
      ALGS(3,3)=ALGS(2,2)
      ALGS(3,4)=ALGS(2,4)
      ALGS(4,1)=ALGS(2,2)
      ALGS(4,2)=ALGS(2,3)
      ALGS(4,3)=ALGS(2,1)
      ALGS(4,4)=ALGS(2,4)
      ALGS(5,1)=0.059715E0
      ALGS(5,2)=0.470142E0
      ALGS(5,3)=0.470142E0
      ALGS(5,4)=0.132394E0
      ALGS(6,1)=ALGS(5,3)
      ALGS(6,2)=ALGS(5,1)
      ALGS(6,3)=ALGS(5,2)
      ALGS(6,4)=ALGS(5,4)
      ALGS(7,1)=ALGS(5,2)
      ALGS(7,2)=ALGS(5,3)
      ALGS(7,3)=ALGS(5,1)
      ALGS(7,4)=ALGS(5,4)
C
C MENSEKI SEKIBUN TOKUITEN SHORI NO OMOMI TO ZAHYOU
C
      IP=1
      WRITE(*,*) 'MENSEKI BUNKATSU=',NDS
      ZP=0.E0
      ZL=0.E0
      DO 30 I=1,NDS
         ZP=ZP+ZL
         ZL=0.5E0**(I+1)
         ZZ(IP,1)=ZP
         ZZ(IP,2)=ZP
         IP=IP+1
         DO 40 J=1,3
            ZZ(IP,1)=ZP+ZL*FLOAT(J+1)
            ZZ(IP,2)=ZP
            IP=IP+1
            ZZ(IP,1)=ZZ(IP-1,2)
            ZZ(IP,2)=ZZ(IP-1,1)
            IP=IP+1
   40    CONTINUE
         ZZ(IP,1)=ZP+ZL
         ZZ(IP,2)=ZP+ZL
         IP=IP+1
         DO 50 J=1,2
            ZZ(IP,1)=ZP+ZL*FLOAT(J+1)
            ZZ(IP,2)=ZP+ZL
            IP=IP+1
            ZZ(IP,1)=ZZ(IP-1,2)
            ZZ(IP,2)=ZZ(IP-1,1)
            IP=IP+1
   50    CONTINUE
   30 CONTINUE
      ZZ(IP,1)=ZP+2.E0*ZL
      ZZ(IP,2)=ZP+2.E0*ZL
      NP=IP
      IP=0
      DO 60 I=1,NDS
         INIP=12*(I-1)+1
C 1
         IP=IP+1
         NOE(IP,1)=INIP
         NOE(IP,2)=INIP+1
         NOE(IP,3)=INIP+2
C 2
         IP=IP+1
         NOE(IP,1)=INIP+1
         NOE(IP,2)=INIP+8
         NOE(IP,3)=INIP+7
C 3      
         IP=IP+1
         NOE(IP,1)=INIP+1
         NOE(IP,2)=INIP+3
         NOE(IP,3)=INIP+8
C 4
         IP=IP+1
         NOE(IP,1)=INIP+3
         NOE(IP,2)=INIP+10
         NOE(IP,3)=INIP+8
C 5
         IP=IP+1
         NOE(IP,1)=INIP+3
         NOE(IP,2)=INIP+5
         NOE(IP,3)=INIP+10
C 6
         IP=IP+1
         NOE(IP,1)=INIP+2
         NOE(IP,2)=INIP+7
         NOE(IP,3)=INIP+9
C 7
         IP=IP+1
         NOE(IP,1)=INIP+2
         NOE(IP,2)=INIP+9
         NOE(IP,3)=INIP+4
C 8
         IP=IP+1
         NOE(IP,1)=INIP+4
         NOE(IP,2)=INIP+9
         NOE(IP,3)=INIP+11
C 9      
         IP=IP+1
         NOE(IP,1)=INIP+4
         NOE(IP,2)=INIP+11
         NOE(IP,3)=INIP+6
   60 CONTINUE
      IP=IP+1
      NOE(IP,1)=INIP+7
      NOE(IP,2)=INIP+8
      NOE(IP,3)=INIP+9
      IP=IP+1
      NOE(IP,1)=INIP+8
      NOE(IP,2)=INIP+10
      NOE(IP,3)=INIP+12
      IP=IP+1
      NOE(IP,1)=INIP+8
      NOE(IP,2)=INIP+12
      NOE(IP,3)=INIP+9
      IP=IP+1
      NOE(IP,1)=INIP+9
      NOE(IP,2)=INIP+12
      NOE(IP,3)=INIP+11
      NE=IP
C
      IP=1
      DO 70 I=1,NE
         DO 80 J=1,7
             SL1(1,IP)=ZZ(NOE(I,1),1)*ALGS(J,1)
     &          +ZZ(NOE(I,2),1)*ALGS(J,2)
     &          +ZZ(NOE(I,3),1)*ALGS(J,3)
             SL2(1,IP)=ZZ(NOE(I,1),2)*ALGS(J,1)
     &          +ZZ(NOE(I,2),2)*ALGS(J,2)
     &          +ZZ(NOE(I,3),2)*ALGS(J,3)
             SL3(1,IP)=1.E0-SL1(1,IP)-SL2(1,IP)
             X1=ZZ(NOE(I,1),1)
             Y1=ZZ(NOE(I,1),2)
             X2=ZZ(NOE(I,2),1)
             Y2=ZZ(NOE(I,2),2)
             X3=ZZ(NOE(I,3),1)
             Y3=ZZ(NOE(I,3),2)
             SWW(1,IP)=((X2-X1)*(Y3-Y1)-(Y2-Y1)*(X3-X1))
     &            *ALGS(J,4)*1.E0
             IP=IP+1
   80    CONTINUE
   70 CONTINUE
      NCE=IP-1
C
      DO 90 I=1,NCE
         SL1(2,I)=SL3(1,I)
         SL2(2,I)=SL1(1,I)
         SL3(2,I)=SL2(1,I)
         SWW(2,I)=SWW(1,I)
         SL1(3,I)=SL1(1,I)
         SL2(3,I)=SL3(1,I)
         SL3(3,I)=SL2(1,I)
         SWW(3,I)=SWW(1,I)
   90 CONTINUE
      RETURN
      END
C
      SUBROUTINE OUTDATA
C
      PARAMETER (NPPN=1000,NEEN=2000,NVVN=3000,IFMX=2)
      implicit double precision ( A-H, O-Z)
      COMPLEX*16 CCE
      COMMON /UNKNOWN1/NPI,IO(1000)
      COMMON /UNKNOWN2/ZZZ(1000,3)
C      COMMON /UNKNOWN/NPI,ZZZ(1000,3),IO(1000)
      COMMON /E0B0/B0(1000,3),A0(1000,3),CCE(1000,6)
      COMMON /FNM/FN1,FN2,FN3,FN4,FN5,FN6,FN7,FN8,FN9,FN10,FN11,FN12
      CHARACTER*20 FN1,FN2,FN3,FN4,FN5,
     &             FN6,FN7,FN8,FN9,FN10,FN11,FN12
C
      OPEN( 21,FILE=FN10,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      WRITE(21,'(I4)')NPI
      WRITE(21,'(E15.8)')((ZZZ(I,J),J=1,3),I=1,NPI)
      WRITE(21,'(E15.8)')((CCE(I,J),J=1,3),I=1,NPI)
      CLOSE(21)
      OPEN( 21,FILE=FN11,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      WRITE(21,'(I4)')NPI
      WRITE(21,'(E15.8)')((ZZZ(I,J),J=1,3),I=1,NPI)
      WRITE(21,'(E15.8)')((CCE(I,J),J=4,6),I=1,NPI)
      CLOSE(21)
C
      RETURN
      END
C
      SUBROUTINE VCTPREP
      PARAMETER (NPPN=1000,NEEN=2000,NVVN=3000,IFMX=2)
      implicit double precision ( A-H, O-Z)
      COMMON /VCT/ZLN(NVVN,6),NOV(NVVN,2),NV
      COMMON /TR1/NP,NE3,ZZ(NPPN,3),NOE3(NEEN,6),ZNP(NEEN,3)
      COMMON /TR4/BN(NEEN,3),T1(NEEN,3),T2(NEEN,3),ZI(NEEN,3)
      COMMON /VPR0/E1XE2(NEEN),E2XE3(NEEN),E1XE3(NEEN),E2XE1(NEEN)
      COMMON /VPR1/SXYZ(NEEN),UNXE2(NEEN,3), UNXE1(NEEN,3)
C
      DO 100 II=1,NE3
C
      I1=NOE3(II,4)
      I2=NOE3(II,5)
      I3=NOE3(II,6)
      A1=I1/IABS(I1)
      A2=I2/IABS(I2)
      A3=I3/IABS(I3)
      I1=IABS(I1)
      I2=IABS(I2)
      I3=IABS(I3)
C
      E1X=ZLN(I1,1)*A1
      E1Y=ZLN(I1,2)*A1
      E1Z=ZLN(I1,3)*A1
      E2X=ZLN(I2,1)*A2
      E2Y=ZLN(I2,2)*A2
      E2Z=ZLN(I2,3)*A2
      E3X=ZLN(I3,1)*A3
      E3Y=ZLN(I3,2)*A3
      E3Z=ZLN(I3,3)*A3
      X1=ZZ(NOE3(II,1),1)
      Y1=ZZ(NOE3(II,1),2)
      Z1=ZZ(NOE3(II,1),3)
      X2=ZZ(NOE3(II,2),1)
      Y2=ZZ(NOE3(II,2),2)
      Z2=ZZ(NOE3(II,2),3)
      X3=ZZ(NOE3(II,3),1)
      Y3=ZZ(NOE3(II,3),2)
      Z3=ZZ(NOE3(II,3),3)
      AX=X2-X1
      AY=Y2-Y1
      AZ=Z2-Z1
      BX=X3-X2
      BY=Y3-Y2
      BZ=Z3-Z2
      BNX=AY*BZ-AZ*BY
      BNY=AZ*BX-AX*BZ
      BNZ=AX*BY-AY*BX
      SXYZ(II)=SQRT(BNX*BNX+BNY*BNY+BNZ*BNZ)*.5
      BNX=BN(II,1)
      BNY=BN(II,2)
      BNZ=BN(II,3)
C
      CALL CROS(V1,V2,V3,E1X,E1Y,E1Z,E2X,E2Y,E2Z)
      CALL DOT( E1XE2(II), V1, V2, V3, BNX,BNY,BNZ)
      CALL CROS(V1,V2,V3,E2X,E2Y,E2Z,E3X,E3Y,E3Z)
      CALL DOT( E2XE3(II), V1, V2, V3, BNX,BNY,BNZ)
      CALL CROS(V1,V2,V3,E1X,E1Y,E1Z,E3X,E3Y,E3Z)
      CALL DOT( E1XE3(II), V1, V2, V3, BNX,BNY,BNZ)
      CALL CROS(V1,V2,V3,E2X,E2Y,E2Z,E1X,E1Y,E1Z)
      CALL DOT( E2XE1(II), V1, V2, V3, BNX,BNY,BNZ)
C
      CALL CROS( UNXE2(II,1),UNXE2(II,2),UNXE2(II,3),
     &           BNX,      BNY,   BNZ, E2X, E2Y, E2Z)
C
      CALL CROS( UNXE1(II,1),UNXE1(II,2),UNXE1(II,3),
     &           BNX,      BNY,   BNZ, E1X, E1Y, E1Z)
C
  100 CONTINUE
      RETURN
      END

