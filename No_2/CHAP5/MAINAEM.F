C  *********************************************************************
C  *****  MAINAEMX                                          V4.10  *****
C  *********************************************************************
C                                               ( Apr 26 16:51:42 1994 )
C  ==(......)***************************************************(....)==
      INCLUDE"AEMX0"
C  =====================================================================
      PROGRAM MAINAEMX
      INCLUDE"AEMX1"
      CALL DATSUB(59)
      WRITE(6,*) 'FREQ=',FRQ,'Hz'
      CALL KPTSUB
      CALL NODSUB
      CALL MEDSUB
      WRITE(6,*) 'UNKNOWN =',NQM
      CALL ISMSUB
      CALL QMXA1C
      CALL SRCSUB
      CALL QMXS1C(1.D-10,2*NQM,NQM/20)
      CALL FLDELM(31,5)
      STOP
      END
C  =====================================================================
      SUBROUTINE DATSUB(IO)
      INCLUDE"AEMX1"
      READ (IO  ,*   ) NVX,NEM,MDV
      DO 100 N =1,NVX
        READ  (IO  ,*   ) (XVX(N,M),M=1,3)
  100 CONTINUE
      DO 110 L =1,NEM
        READ  (IO  ,*   ) (KVX(L,M),M=1,MDV)
     &                   ,(NBC(L,M),M=1,MDV),NMD(L)
  110 CONTINUE
      READ(IO  ,*   ) IMD,FRQ
      DO 200 M =1,IMD
        READ(IO  ,*   ) PEL(M),PMG(M),PCD(M)
  200 CONTINUE
      READ(IO  ,*   ,END=900) KSC
      BACKSPACE IO
      READ(IO  ,*   ) KSC
      GO TO (31,32),KSC
   31 READ(IO  ,*   ) ICS
      DO 310 I =1,ICS
        READ(IO  ,*   ) NCS(I),ACS(I)
        READ(IO  ,*   ) (LVX(I,J),J=1,NCS(I))
  310 CONTINUE
      RETURN
   32 READ(IO  ,*   ) IPS
      DO 320 I =1,IPS
        READ(IO  ,*   ) (VPS(I,J),J=1,3)
  320 CONTINUE
      RETURN
  900 RETURN
      END
C  =====================================================================
      SUBROUTINE KPTSUB
      INCLUDE"AEMX1"
      NPT=0
      DO 100 L =1,NEM
        DO 110 I =1,12
          NN=NPT+1
          DO 120 M =1,3
            XPT(NN,M)=0.D0
            DO 130 N =1,4
              XPT(NN,M)=XPT(NN,M)+XVX(KVX(L,N),M)*RPT(I,N)
  130       CONTINUE
  120     CONTINUE
          DO 140 N =1,NPT
            DO 150 M =1,3
              IF (ABS(XPT(NN,M)-XPT(N,M)).GT.1.D-10)
     &                                   GO TO 140
  150       CONTINUE
            KPT(L,I)=N
            GO TO 110
  140     CONTINUE
          NPT=NPT+1
          KPT(L,I)=NPT
  110   CONTINUE
  100 CONTINUE
      RETURN
      END
C  =====================================================================
      SUBROUTINE NODSUB
      INCLUDE"AEMX1"
      DO 100 N =1,NPT
        NOD(N)=1
  100 CONTINUE
      DO 200 L =1,NEM
        DO 210 I =1,4
          IF (NBC(L,I)) 21,210,22
   21     DO 220 M =1,12
            IF (NOD(KPT(L,M)).EQ.0) GO TO 220
            IF (RPT(M,I).GT.0.D0) GO TO 220
            NOD(KPT(L,M))=NBC(L,I)
  220     CONTINUE
          GO TO 210
   22     DO 230 M =1,12
            IF (RPT(M,I).GT.0.D0) GO TO 230
            NOD(KPT(L,M))=0
  230     CONTINUE
  210   CONTINUE
  200 CONTINUE
      NQM=0
      DO 300 N =1,NPT
        IF (NOD(N).LE.0) GO TO 300
        NQM=NQM+1
        NOD(N)=NQM
  300 CONTINUE
      NQA=NQM
      DO 310 N =1,NPT
        IF (NOD(N).GE.0) GO TO 310
        NQA=NQA+1
        NOD(N)=-10*NQA+NOD(N)
  310 CONTINUE
      RETURN
      END
C  =====================================================================
      SUBROUTINE MEDSUB(IO)
      INCLUDE"AEMX1"
      FOG=(2.D0*PAI*FRQ)**2
      DO 100 M =1,IMD
        RKT(M)=1.D0/PHM/PMG(M)
        RMS(M)=PHE*PEL(M)
        IF (FRQ.GT.0.D0) RMS(M)=RMS(M)
     &                   -(0.D0,1.D0)*PCD(M)/2.D0/PAI/FRQ
  100 CONTINUE
      RETURN
      END
C  =====================================================================
      SUBROUTINE ISMSUB
      INCLUDE"AEMX1"
      DIMENSION IW(OQM)
      NQL=0
      DO 100 N =1,NQM
        ISM(N)=NQL+1
        DO 110 I =1,N
          IW(I)=0
  110   CONTINUE
        DO 120 L =1,NEM
          DO 130 I =1,12
            NA=NOD(KPT(L,I))
            IF (NA.LE.0) GO TO 130
            IF (NA.NE.N) GO TO 130
            DO 140 J =1,12
              NB=NOD(KPT(L,J))
              IF (NB.LE.0) GO TO 140
              IW(NB)=1
  140       CONTINUE
  130     CONTINUE
  120   CONTINUE
        DO 150 I =1,N
          IF (IW(I).EQ.0) GO TO 150
          NQL=NQL+1
          JSM(NQL)=I
  150   CONTINUE
  100 CONTINUE
      ISM(NQM+1)=NQL+1
      RETURN
      END
C  =====================================================================
      SUBROUTINE QMXA1C
      INCLUDE"AEMX1"
      DO 100 N =1,NQL
        QKT(N)=(0.D0,0.D0)
  100 CONTINUE
      DO 200 L =1,NEM
        MD=NMD(L)
        CALL DMXSUB(L)
        CALL PMXSUB
        DO 210 M =1,12
          MM=NOD(KPT(L,M))
          IF (MM.LE.0) GO TO 210
          DO 220 N =1,12
            NN=NOD(KPT(L,N))
            IF (NN.LE.0) GO TO 220
            IF (NN.GT.MM) GO TO 220
            DO 250 II=ISM(MM),ISM(MM+1)-1
              MN=JSM(II)
              IF(MN.EQ.NN) GO TO 260
  250       CONTINUE
            GO TO 220
  260       CONTINUE
            QKT(II)=QKT(II)+RKT(MD)*PKT(M,N)/DAR/6.D0
     &             -FOG*RMS(MD)*PMS(M,N)*DAR/120.D0
  220     CONTINUE
  210   CONTINUE
  200 CONTINUE
      END
C  =====================================================================
      SUBROUTINE DMXSUB(L)
      INCLUDE"AEMX1"
      DIMENSION IA(12),IB(12)
      DO 100 M =1,4
        NS=KVX(L,M)
        DMX(1,M)=1.D0
        DMX(2,M)=XVX(NS,1)
        DMX(3,M)=XVX(NS,2)
        DMX(4,M)=XVX(NS,3)
  100 CONTINUE
      CALL DBM2DT(DMX,12,4,DAR)
      DO 200 I =1,12
        XX=0.D0
        DO 210 M =1,3
          QHS(I,M)=XVX(KVX(L,NDR(I,2)),M)
     &            -XVX(KVX(L,NDR(I,1)),M)
          XX=XX+QHS(I,M)*QHS(I,M)
  210   CONTINUE
        XX=SQRT(XX)
        DO 220 M =1,3
          QHS(I,M)=QHS(I,M)/XX
  220   CONTINUE
        DO 230 M =1,3
          DO 240 N =1,4
            DMX(4*(M-1)+N,I)=RPT(I,N)*QHS(I,M)
  240     CONTINUE
  230   CONTINUE
  200 CONTINUE
      CALL DBM2IV(DMX,12,12,IA,IB)
      DO 300 N =1,12
        DO 310 M =1,3
          BMX(N,M)=0.D0
  310   CONTINUE
  300 CONTINUE
      DO 320 N =1,4
        N1=MOD(N  ,4)+1
        N2=MOD(N+1,4)+1
        N3=MOD(N+2,4)+1
        DO 330 M =1,3
          M1=MOD(M  ,3)+1
          M2=MOD(M+1,3)+1
          BB= XVX(KVX(L,N2),M1)*XVX(KVX(L,N3),M2)
     &       -XVX(KVX(L,N3),M1)*XVX(KVX(L,N2),M2)
     &       +XVX(KVX(L,N3),M1)*XVX(KVX(L,N1),M2)
     &       -XVX(KVX(L,N1),M1)*XVX(KVX(L,N3),M2)
     &       +XVX(KVX(L,N1),M1)*XVX(KVX(L,N2),M2)
     &       -XVX(KVX(L,N2),M1)*XVX(KVX(L,N1),M2)
          BB=(-1)**N*BB
          BMX(4*(M2-1)+N,M1)= BB
          BMX(4*(M1-1)+N,M2)=-BB
  330   CONTINUE
  320 CONTINUE
      RETURN
      END
C  =====================================================================
      SUBROUTINE PMXSUB
      INCLUDE"AEMX1"
      DIMENSION W1(12,12),W2(12,12)
      CALL DAM2MU(DMX,12,12,12,SMX,12,12,W1 ,12)
      CALL DAMT2M( W1,12,12,12,DMX,12,12,PMS,12)
      CALL DAMT2M(BMX,12,12,3,BMX,12,12,W1 ,12)
      CALL DAM2MU(DMX,12,12,12, W1,12,12,W2 ,12)
      CALL DAMT2M( W2,12,12,12,DMX,12,12,PKT,12)
      RETURN
      END
C  =====================================================================
      SUBROUTINE SRCSUB
      INCLUDE"AEMX1"
      DO 100 N =1,NQA
        QM1(N)=(0.D0,0.D0)
  100 CONTINUE
      IF (KSC.EQ.1) CALL SRCLCR
      IF (KSC.EQ.2) CALL SRCCVP
      RETURN
      END
C  =====================================================================
      SUBROUTINE SRCCVP
      INCLUDE"AEMX1"
      DO 100 L =1,NEM
        MD=NMD(L)
        CALL DMXSUB(L)
        CALL PMXSUB
        DO 110 M =1,12
          MM=NOD(KPT(L,M))
          IF (MM.LE.0) GO TO 110
          DO 120 N =1,12
            NN=NOD(KPT(L,N))
            IF (NN.GE.0) GO TO 120
            NA=MOD(NN,10)
            CF=0.D0
            DO 130 I =1,3
              CF=CF+QHS(N,I)*VPS(-NA,I)
  130       CONTINUE
            QM2(-NN/10)=CF
            QM1(MM)=QM1(MM)-RKT(MD)*PKT(M,N)/DAR/6.D0*CF
  120     CONTINUE
  110   CONTINUE
  100 CONTINUE
      RETURN
      END
C  =====================================================================
      SUBROUTINE SRCLCR
      INCLUDE"AEMX1"
      DIMENSION XP(3)
      DIMENSION IW(OQM)
      DO 100 I =1,ICS
        DO 200 J =1,NCS(I)-1
          N1=LVX(I,J)
          N2=LVX(I,J+1)
          DO 210 L =1,NEM
            DO 212 M =1,4
              IF (KVX(L,M).EQ.N1) GO TO 20
  212       CONTINUE
            GO TO 210
   20       M1=M
            DO 214 M =1,4
              IF (KVX(L,M).EQ.N2) GO TO 22
  214       CONTINUE
  210     CONTINUE
          GO TO 200
   22     M2=M
          XL=0.D0
          DO 220 M =1,3
            XL=XL+(XVX(N1,M)-XVX(N2,M))**2
  220     CONTINUE
          XL=SQRT(XL)
          DO 300 M =1,12
            IF (NDR(M,1).NE.M1) GO TO 300
            IF (NDR(M,2).EQ.M2) GO TO 30
  300     CONTINUE
   30     CONTINUE
          MM=M
          DO 310 M =1,3
            XP(M)=0.D0
            DO 312 K =1,4
              XP(M)=XP(M)+XVX(KVX(L,K),M)*RPT(MM,K)
  312       CONTINUE
  310     CONTINUE
          DO 320 N =1,NPT
            DO 322 M =1,3
              IF (ABS(XPT(N,M)-XP(M)).GT.1.D-10) GO TO 320
  322       CONTINUE
            GO TO 32
  320     CONTINUE
   32     ND=NOD(N)
          IF (ND.LE.0) GO TO 200
          QM1(ND)=QM1(ND)-ACS(I)*XL/2.D0
          DO 400 M =1,12
            IF (NDR(M,1).NE.M2) GO TO 400
            IF (NDR(M,2).EQ.M1) GO TO 40
  400     CONTINUE
   40     CONTINUE
          MM=M
          DO 410 M =1,3
            XP(M)=0.D0
            DO 412 K =1,4
              XP(M)=XP(M)+XVX(KVX(L,K),M)*RPT(MM,K)
  412       CONTINUE
  410     CONTINUE
          DO 420 N =1,NPT
            DO 422 M =1,3
              IF (XPT(N,M).NE.XP(M)) GO TO 420
  422       CONTINUE
            GO TO 42
  420     CONTINUE
   42     ND=NOD(N)
          IF (ND.LE.0) GO TO 200
          QM1(ND)=QM1(ND)+ACS(I)*XL/2.D0
  200   CONTINUE
  100 CONTINUE
      RETURN
      END
C  =====================================================================
      SUBROUTINE QMXS1C(R1,I1,I2)
      INCLUDE"AEMX1"
      COMPLEX*16 R(OQM),P(OQM),AP(OQM),S(OQM),DU,BT,AL,WW,WA,WB,WC
      I3=0
      DO 100 I =1,NQM
        II=ISM(I+1)-1
        CC=(1.D0,0.D0)/CDSQRT(QKT(II))
        S(I)=CC
        QM1(I)=QM1(I)*CC
        DO 110 IJ=ISM(I),ISM(I+1)-1
          J=JSM(IJ)
          QKT(IJ)=QKT(IJ)*CC*S(J)
  110   CONTINUE
  100 CONTINUE
      R3=0.D0
      DO 120 I =1,NQM
        QM2(I)=(0.D0,0.D0)
        P(I)=(0.D0,0.D0)
        R(I)=QM1(I)
        CC=QM1(I)
        R3=R3+CONJG(CC)*CC
  120 CONTINUE
      C0=(1.D0,0.D0)
   60 CONTINUE
      I3=I3+1
      DU=(0.D0,0.D0)
      DO 210 I =1,NQM
        DU=DU+R(I)*R(I)
  210 CONTINUE
      BT=DU/C0
      C0=DU
      DO 220 I =1,NQM
        P(I)=R(I)+BT*P(I)
  220 CONTINUE
      DO 310 I =1,NQM
        WA=P(I)
        WW=(0.D0,0.D0)
        DO 312 IJ=ISM(I),ISM(I+1)-2
          J=JSM(IJ)
          WC=QKT(IJ)
          WW=WW+WC*P(J)
          AP(J)=AP(J)+WC*WA
  312   CONTINUE
        II=ISM(I+1)-1
        AP(I)=QKT(II)*WA+WW
  310 CONTINUE
      WA=(0.D0,0.D0)
      WB=(0.D0,0.D0)
      DO 320 I =1,NQM
        WA=WA+P(I)*R(I)
        WB=WB+P(I)*AP(I)
  320 CONTINUE
      AL=WA/WB
      DO 330 I =1,NQM
        QM2(I)=QM2(I)+AL*P(I)
        R(I)=R(I)-AL*AP(I)
  330 CONTINUE
      R2=0.D0
      DO 400 I =1,NQM
        CC=R(I)
        R2=R2+CONJG(CC)*CC
  400 CONTINUE
      R2=DSQRT(R2/R3)
      CALL LOGIND(R1,R2,I2,I3)
      IF (R2.LT.R1) GO TO 90
      IF (I3.GE.I1) GO TO 90
      GO TO 60
   90 CONTINUE
      DO 200 I =1,NQM
        QM2(I)=QM2(I)*S(I)
  200 CONTINUE
      RETURN
      END
C  =====================================================================
      SUBROUTINE LOGIND(R1,R2,I2,I3)
      INCLUDE"AEMX1"
      CHARACTER*1 ICH(100)
      IF (MOD(I3,I2).NE.0) RETURN
      IM=66
      DO 100 I =1,IM
        ICH(I)='.'
  100 CONTINUE
      DO 110 I =1,IM/5+1
        ICH(5*I+1)='+'
  110 CONTINUE
      II=-5.D0*LOG10(R1)+1.5D0
      II=MIN(MAX(II,1),IM)
      ICH(II)='>'
      JJ=-5.D0*LOG10(R2)+1.5D0
      JJ=MIN(MAX(JJ,1),IM)
      ICH(JJ)='*'
      WRITE (6   ,2000) I3,(ICH(I),I=II,JJ,-1)
 2000 FORMAT(3X,I7,1X,100A1)
      RETURN
      END
C  =====================================================================
      SUBROUTINE FLDELM(IO,IS)
      INCLUDE"AEMX1"
      DIMENSION QX(2,12),WW(2,12)
      DIMENSION FD(2,3)
      DIMENSION PS(OEM,3)
      DIMENSION CF(OEM,3)
      CIM=(0.D0,1.D0)
      DO 100 L =1,NEM
        MD=NMD(L)
        DO 102 M =1,3
          PS(L,M)=0.D0
          DO 104 N =1,4
            PS(L,M)=PS(L,M)+XVX(KVX(L,N),M)
  104     CONTINUE
          PS(L,M)=PS(L,M)/4
  102   CONTINUE
        CALL DMXSUB(L)
        DO 110 M =1,12
          QX(1,M)=0.D0
          QX(2,M)=0.D0
          MM=NOD(KPT(L,M))
          IF (MM) 10,110,12
   10     QX(1,M)=REAL(QM2(-MM/10))
          QX(2,M)=IMAG(QM2(-MM/10))
          GO TO 110
   12     QX(1,M)=REAL(QM2(MM))
          QX(2,M)=IMAG(QM2(MM))
  110   CONTINUE
        CALL DAM2MU(QX,2,2,12,DMX,12,12,WW,2)
        GO TO (21,22,23,24,25,26),IS
   21   CALL DAM2MU(WW,2,2,12,RMX,12,3,FD,2)
        DO 210 M =1,3
          CC=-CIM*2.D0*PAI*FRQ*(FD(1,M)+CIM*FD(2,M))
          FD(1,M)=REAL(CC)
          FD(2,M)=IMAG(CC)
  210   CONTINUE
        GO TO 29
   22   CALL DAM2MU(WW,2,2,12,RMX,12,3,FD,2)
        DO 220 M =1,3
          CC=-CIM*2.D0*PAI*FRQ*RMS(MD)*(FD(1,M)+CIM*FD(2,M))
          FD(1,M)=REAL(CC)
          FD(2,M)=IMAG(CC)
  220   CONTINUE
        GO TO 29
   23   CALL DAM2MU(WW,2,2,12,BMX,12,3,FD,2)
        DO 230 M =1,3
          FD(1,M)=RKT(MD)*FD(1,M)/DAR
          FD(2,M)=RKT(MD)*FD(2,M)/DAR
  230   CONTINUE
        GO TO 29
   24   CALL DAM2MU(WW,2,2,12,BMX,12,3,FD,2)
        DO 240 M =1,3
          FD(1,M)=FD(1,M)/DAR
          FD(2,M)=FD(2,M)/DAR
  240   CONTINUE
        GO TO 29
   25   CALL DAM2MU(WW,2,2,12,RMX,12,3,FD,2)
        DO 250 M =1,3
          CC=-CIM*2.D0*PAI*FRQ*PCD(MD)*(FD(1,M)+CIM*FD(2,M))
          FD(1,M)=REAL(CC)
          FD(2,M)=IMAG(CC)
  250   CONTINUE
        GO TO 29
   26   CALL DAM2MU(WW,2,2,12,RMX,12,3,FD,2)
   29   CONTINUE
        DO 290 M =1,3
          CF(L,M)=CMPLX(FD(1,M),FD(2,M))
  290   CONTINUE
  100 CONTINUE
      DO 400 L =1,NEM
        WRITE (IO  ,2000) (PS(L,M),M=1,3),(REAL(CF(L,M)),M=1,3)
     &                                   ,(IMAG(CF(L,M)),M=1,3)
  400 CONTINUE
      RETURN
 2000 FORMAT(1X,9D10.3)
      END
C  =====================================================================
      SUBROUTINE DAM2MU(A,LMA,M,N,B,LMB,L,C,LMC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(LMA,N),B(LMB,L),C(LMC,L)
      DO 100 I =1,M
        DO 110 J =1,L
          W=0.D0
          DO 120 K =1,N
            W=W+A(I,K)*B(K,J)
  120     CONTINUE
          C(I,J)=W
  110   CONTINUE
  100 CONTINUE
      RETURN
      END
C  =====================================================================
      SUBROUTINE DAMT2M(A,LMA,M,N,B,LMB,L,C,LMC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(LMA,N),B(LMB,N),C(LMC,L)
      DO 100 I =1,M
        DO 110 J =1,L
          W=0.D0
          DO 120 K =1,N
            W=W+A(I,K)*B(J,K)
  120     CONTINUE
          C(I,J)=W
  110   CONTINUE
  100 CONTINUE
      RETURN
      END
C  =====================================================================
      SUBROUTINE DBM2DT(A,LMA,N,DET)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(LMA,N)
      PV=0.D0
      DO 100 I =1,N
        DO 110 J =1,N
          PW=ABS(A(I,J))
          IF (PW.LE.PV) GO TO 110
          PV=PW
          L=I
          M=J
  110   CONTINUE
  100 CONTINUE
      EPS=PV*DEXP(DBLE(-14)*DLOG(10.D0))
      DET=1.D0
      DO 200 K =1,N-1
        IF (K.NE.L) THEN
          DET=-DET
          DO 210 J =1,N
            FF=A(K,J)
            A(K,J)=A(L,J)
            A(L,J)=FF
  210     CONTINUE
        ENDIF
        IF (K.NE.M) THEN
          DET=-DET
          DO 220 I =1,N
            FF=A(I,K)
            A(I,K)=A(I,M)
            A(I,M)=FF
  220     CONTINUE
        ENDIF
        PP=1.D0/A(K,K)
        DET=DET*A(K,K)
        DO 230 I =K+1,N
          FF=A(I,K)*PP
          DO 240 J =K+1,N
            A(I,J)=A(I,J)-A(K,J)*FF
  240     CONTINUE
  230   CONTINUE
        PV=0.D0
        DO 250 I =K+1,N
          DO 260 J =K+1,N
            PW=ABS(A(I,J))
            IF (PW.LE.PV) GO TO 260
            PV=PW
            L=I
            M=J
  260     CONTINUE
  250   CONTINUE
  200 CONTINUE
      DET=DET*A(N,N)
      RETURN
      END
C  =====================================================================
      SUBROUTINE DBM2IV(A,LMA,N,IW1,IW2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(LMA,N),IW1(N),IW2(N)
      PV=0.D0
      DO 100 I =1,N
        DO 110 J =1,N
          PW=ABS(A(I,J))
          IF (PW.LE.PV) GO TO 110
          PV=PW
          L=I
          M=J
  110   CONTINUE
  100 CONTINUE
      EPS=PV*DEXP(DBLE(-14)*DLOG(10.D0))
      DO 300 K =1,N
        IW1(K)=L
        IF (K.NE.L) THEN
          DO 210 J =1,N
            FF=A(K,J)
            A(K,J)=A(L,J)
            A(L,J)=FF
  210     CONTINUE
        ENDIF
        IW2(K)=M
        IF (K.NE.M) THEN
          DO 220 I =1,LMA
            FF=A(I,K)
            A(I,K)=A(I,M)
            A(I,M)=FF
  220     CONTINUE
        ENDIF
        PP=1.D0/A(K,K)
        A(K,K)=1.D0
        DO 310 J =1,N
          A(K,J)=A(K,J)*PP
  310   CONTINUE
        IF (K.EQ.N) GO TO 300
        DO 320 I =K+1,N
          FF=A(I,K)
          A(I,K)=0.D0
          DO 330 J =1,N
            A(I,J)=A(I,J)-A(K,J)*FF
  330     CONTINUE
  320   CONTINUE
        PV=0.D0
        DO 340 I =K+1,N
          DO 350 J =K+1,N
            PW=ABS(A(I,J))
            IF (PW.LE.PV) GO TO 350
            PV=PW
            L=I
            M=J
  350     CONTINUE
  340   CONTINUE
  300 CONTINUE
      DO 400 K =2,N
        DO 410 I =1,K-1
          FF=A(I,K)
          A(I,K)=0.D0
          DO 420 J =1,N
            A(I,J)=A(I,J)-A(K,J)*FF
  420     CONTINUE
  410   CONTINUE
  400 CONTINUE
      DO 500 K =N-1,1,-1
        L=IW2(K)
        IF (K.NE.L) THEN
          DO 510 J =1,N
            FF=A(K,J)
            A(K,J)=A(L,J)
            A(L,J)=FF
  510     CONTINUE
        ENDIF
        M=IW1(K)
        IF (K.NE.M) THEN
          DO 520 I =1,LMA
            FF=A(I,K)
            A(I,K)=A(I,M)
            A(I,M)=FF
  520     CONTINUE
        ENDIF
  500 CONTINUE
      RETURN
      END

