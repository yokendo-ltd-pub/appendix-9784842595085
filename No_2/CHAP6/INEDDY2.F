
C
      SUBROUTINE MAT
      PARAMETER (NPPN=1000,NEEN=2000,NVVN=3000,IFMX=2)
      implicit double precision (A-H, O-Z)
      COMPLEX*16 SA,SB,SC
      COMPLEX*16 BI,EI,CCE
      COMMON /TR1/NP,NE3,ZZ(NPPN,3),NOE3(NEEN,6),ZNP(NEEN,3)
      COMMON /TR4/BN(NEEN,3),T1(NEEN,3),T2(NEEN,3),ZI(NEEN,3)
      COMMON /VCT/ZLN(NVVN,6),NOV(NVVN,2),NV
      COMMON /KNOWN/BI(NVVN+NVVN),EI(NVVN+NVVN)
      COMMON /UNKNOWN1/NPI,IO(1000)
      COMMON /UNKNOWN2/ZZZ(1000,3)
      COMMON /E0B0/B0(1000,3),A0(1000,3),CCE(1000,6)
      COMMON /TR6/SA(6,3*(NVVN+NVVN)),SB(6,3*(NVVN+NVVN)),
     &                                SC(6,3*(NVVN+NVVN))
C
      DO 10 I=1,NPI
         DO 10 J=1,6
            CCE(I,J)=0.
   10 CONTINUE
C
      DO 200 I=1,NPI
         II=I
         FX=ZZZ(II,1)
         FY=ZZZ(II,2)
         FZ=ZZZ(II,3)
         CALL ZERO
         WRITE(*,'(1H ,50X,I4)')I
         DO 450 J=1,NE3
            JJ=J
            CALL SFT(FX,FY,FZ,JJ)
  450    CONTINUE
         CALL SETT(II,FX,FY,FZ)
  200 CONTINUE
      RETURN
      END
C
      SUBROUTINE SETT(I,FLX,FLY,FLZ)
      PARAMETER (NPPN=1000,NEEN=2000,NVVN=3000,IFMX=2)
      implicit double precision (A-H, O-Z)
      COMPLEX*16 FEPS1,FEPS2,SA,SB,SC
      COMPLEX*16 ALP,JO,JOME1,JOME2,JK1,JK2
      COMPLEX*16 BI,EI,CCE
C      COMPLEX*16 GFX,GFY,GFZ
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
      COMMON /LIN/NOL(NVVN),NUML,NUMOLD,NUMXYZ
C
      IF( IO(I) .EQ. 1)THEN
         DO 200 J=1,NV
            CCE(I,1)=CCE(I,1)-SB(1,J)*EI(J)-SC(1,J)*BI(J)/JOME1
            CCE(I,2)=CCE(I,2)-SB(2,J)*EI(J)-SC(2,J)*BI(J)/JOME1
            CCE(I,3)=CCE(I,3)-SB(3,J)*EI(J)-SC(3,J)*BI(J)/JOME1
            CCE(I,4)=CCE(I,4)-SB(1,J)*BI(J)-SC(1,J)*EI(J)/(-JO)
            CCE(I,5)=CCE(I,5)-SB(2,J)*BI(J)-SC(2,J)*EI(J)/(-JO)
            CCE(I,6)=CCE(I,6)-SB(3,J)*BI(J)-SC(3,J)*EI(J)/(-JO)
  200    CONTINUE
         DO 300 J=1,NV
            CCE(I,1)=CCE(I,1)+SA(1,J)*BI(J)*JO
            CCE(I,2)=CCE(I,2)+SA(2,J)*BI(J)*JO
            CCE(I,3)=CCE(I,3)+SA(3,J)*BI(J)*JO
            CCE(I,4)=CCE(I,4)+SA(1,J)*EI(J)*(-JOME1)
            CCE(I,5)=CCE(I,5)+SA(2,J)*EI(J)*(-JOME1)
            CCE(I,6)=CCE(I,6)+SA(3,J)*EI(J)*(-JOME1)
  300    CONTINUE
      ELSE
         DO 400 J=1,NV
            CCE(I,1)=CCE(I,1)+SB(4,J)*EI(J)
     &                       +SC(4,J)*BI(J)/(JOME1*ALP)*0.
            CCE(I,2)=CCE(I,2)+SB(5,J)*EI(J)
     &                       +SC(5,J)*BI(J)/(JOME1*ALP)*0.
            CCE(I,3)=CCE(I,3)+SB(6,J)*EI(J)
     &                       +SC(6,J)*BI(J)/(JOME1*ALP)*0.
            CCE(I,4)=CCE(I,4)+SB(4,J)*BI(J)*BETA
     &                       +SC(4,J)*EI(J)/(-JO)
            CCE(I,5)=CCE(I,5)+SB(5,J)*BI(J)*BETA
     &                       +SC(5,J)*EI(J)/(-JO)
            CCE(I,6)=CCE(I,6)+SB(6,J)*BI(J)*BETA
     &                       +SC(6,J)*EI(J)/(-JO)
  400    CONTINUE
         DO 500 J=1,NV
            CCE(I,1)=CCE(I,1)-SA(4,J)*BI(J)*JO*BETA
            CCE(I,2)=CCE(I,2)-SA(5,J)*BI(J)*JO*BETA
            CCE(I,3)=CCE(I,3)-SA(6,J)*BI(J)*JO*BETA
            CCE(I,4)=CCE(I,4)+SA(4,J)*BI(J)*(JOME2)
            CCE(I,5)=CCE(I,5)+SA(5,J)*BI(J)*(JOME2)
            CCE(I,6)=CCE(I,6)+SA(6,J)*BI(J)*(JOME2)
  500    CONTINUE
         CCE(I,1)=CCE(I,1)+A0(I,1)*(-JO)
         CCE(I,2)=CCE(I,2)+A0(I,2)*(-JO)
         CCE(I,3)=CCE(I,3)+A0(I,3)*(-JO)
         CCE(I,4)=CCE(I,4)+B0(I,1)
         CCE(I,5)=CCE(I,5)+B0(I,2)
         CCE(I,6)=CCE(I,6)+B0(I,3)
      ENDIF
      RETURN
      END
C
      SUBROUTINE SFT(FLX,FLY,FLZ,II)
C
      PARAMETER (NPPN=1000,NEEN=2000,NVVN=3000,IFMX=2)
      implicit double precision (A-H, O-Z)
      COMPLEX*16 SA,SB,SC
      COMPLEX*16 FEPS1,FEPS2
      COMPLEX*16 ALP,JO,JOME1,JOME2,JK1,JK2
      COMPLEX*16 FA1,GF1,GF1X,GF1Y,GF1Z
      COMPLEX*16 FA2,GF2,GF2X,GF2Y,GF2Z
      COMPLEX*16 SC11,SC12,SC13,SC21,SC22,SC23,SC31,SC32,SC33,D1,D2,D3
      COMPLEX*16 SD11,SD12,SD13,SD21,SD22,SD23,SD31,SD32,SD33
      COMMON /ITGCE/SL1(3,500),SL2(3,500),SL3(3,500),SWW(3,500),NCE
      COMMON /TRG/ALS(500,4),ALB(1,3,2)
      COMMON /TRS1/SAL1(500),SAL2(500),SWWW(500)
      COMMON /TRS2/NGS
      COMMON /TRSB1/SBL1(1,500),SBL2(1,500),SBWW(1,500)
      COMMON /TRSB2/NGBS
      COMMON /INTEG/NIGN
      COMMON /TR6/SA(6,3*(NVVN+NVVN)),SB(6,3*(NVVN+NVVN)),
     &                                SC(6,3*(NVVN+NVVN))
      COMMON /TR4/BN(NEEN,3),T1(NEEN,3),T2(NEEN,3),ZI(NEEN,3)
      COMMON /TR1/NP,NE3,ZZ(NPPN,3),NOE3(NEEN,6),ZNP(NEEN,3)
      COMMON /VCT/ZLN(NVVN,6),NOV(NVVN,2),NV
      COMMON /MTE0/AMU0,AMU1,AMU2,BETA,EPS0,EPS1,EPS2,OMG
      COMMON /MTE1/SGM1,SGM2,FEPS1,FEPS2,ALP
      COMMON /MTE2/JO,JOME1,JOME2,JK1,JK2
      COMMON /EEE/EE(3,3),DD(3,3)
C
      II3=3*II
      II2=II3-1
      II1=II2-1
C
      IV1=NOE3(II,4)
      IV2=NOE3(II,5)
      IV3=NOE3(II,6)
      IVA1=IABS(IV1)
      IVA2=IABS(IV2)
      IVA3=IABS(IV3)
      A1=IV1/IVA1
      A2=IV2/IVA2
      A3=IV3/IVA3
      X1=ZZ(NOE3(II,1),1)
      X2=ZZ(NOE3(II,2),1)
      X3=ZZ(NOE3(II,3),1)
      Y1=ZZ(NOE3(II,1),2)
      Y2=ZZ(NOE3(II,2),2)
      Y3=ZZ(NOE3(II,3),2)
      Z1=ZZ(NOE3(II,1),3)
      Z2=ZZ(NOE3(II,2),3)
      Z3=ZZ(NOE3(II,3),3)
      CTX=ZI(II,1)
      CTY=ZI(II,2)
      CTZ=ZI(II,3)
      BNX=BN(II,1)
      BNY=BN(II,2)
      BNZ=BN(II,3)
      ED1X=ZLN(IVA1,4)
      ED1Y=ZLN(IVA1,5)
      ED1Z=ZLN(IVA1,6)
      ED2X=ZLN(IVA2,4)
      ED2Y=ZLN(IVA2,5)
      ED2Z=ZLN(IVA2,6)
      ED3X=ZLN(IVA3,4)
      ED3Y=ZLN(IVA3,5)
      ED3Z=ZLN(IVA3,6)
C
      AX=X1-X3
      AY=Y1-Y3
      AZ=Z1-Z3
      BX=X2-X3
      BY=Y2-Y3
      BZ=Z2-Z3
      CX=AY*BZ-AZ*BY
      CY=AZ*BX-AX*BZ
      CZ=AX*BY-AY*BX
      CC=SQRT(CX*CX+CY*CY+CZ*CZ)
      R1=SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
      R2=SQRT((X3-X2)**2+(Y3-Y2)**2+(Z3-Z2)**2)
      R3=SQRT((X1-X3)**2+(Y1-Y3)**2+(Z1-Z3)**2)
      RL=SQRT((CTX-FLX)**2+(CTY-FLY)**2+(CTZ-FLZ)**2)
      RL1=SQRT((ED1X-FLX)**2+(ED1Y-FLY)**2+(ED1Z-FLZ)**2)
      RL2=SQRT((ED2X-FLX)**2+(ED2Y-FLY)**2+(ED2Z-FLZ)**2)
      RL3=SQRT((ED3X-FLX)**2+(ED3Y-FLY)**2+(ED3Z-FLZ)**2)
      NIGS=NIGN
      IF(RL.LT.R1)     NIGS = NGS
      IF(RL.LT.R2)     NIGS = NGS
      IF(RL.LT.R3)     NIGS = NGS
      IS1=0
      IF(RL1.LT.1.E-20) IS1=3
      IF(RL2.LT.1.E-20) IS1=2
      IF(RL3.LT.1.E-20) IS1=1
      IF(IS1.NE.0) NIGS=NCE
      DO 10 IG=1,NIGS
         AL1=ALS(IG,1)
         AL2=ALS(IG,2)
         W=ALS(IG,4)
         IF(NIGS.EQ.NGS)THEN
            AL1=SAL1(IG)
            AL2=SAL2(IG)
            W=SWWW(IG)
         ENDIF
         IF(NIGS.EQ.NGBS)THEN
            AL1= SBL1(1,IG)
            AL2= SBL2(1,IG)
            W  = SBWW(1,IG)
         ENDIF
         IF(IS1.NE.0)THEN
            AL1= SL1(IS1,IG)
            AL2= SL2(IS1,IG)
            W  = SWW(IS1,IG)
         ENDIF
         AL3=1.-AL1-AL2
         RX=FLX-(X3*AL3+X1*AL1+X2*AL2)
         RY=FLY-(Y3*AL3+Y1*AL1+Y2*AL2)
         RZ=FLZ-(Z3*AL3+Z1*AL1+Z2*AL2)
         WP=CC*W*.5
         R=SQRT(RX*RX+RY*RY+RZ*RZ)
         FA1=ZEXP(-JK1*R)/(4.*3.141592654*R)*WP
         GF1=(1./R+JK1)*FA1/R
         FA2=ZEXP(-JK2*R)/(4.*3.141592654*R)*WP
         GF2=(1./R+JK2)*FA2/R
         GF1X=GF1*RX
         GF1Y=GF1*RY
         GF1Z=GF1*RZ
         GF2X=GF2*RX
         GF2Y=GF2*RY
         GF2Z=GF2*RZ
         CALL CALJC(II,EE,DD,AL1,AL2)
C
C        ( n x F ) x grad fai
C
         SC11=BNZ*GF1Z+BNY*GF1Y
         SC12=-BNX*GF1Y
         SC13=-BNX*GF1Z
         SC21=-BNY*GF1X
         SC22=BNZ*GF1Z+BNX*GF1X
         SC23=-BNY*GF1Z
         SC31=-BNZ*GF1X
         SC32=-BNZ*GF1Y
         SC33=BNY*GF1Y+BNX*GF1X
         SD11=BNZ*GF2Z+BNY*GF2Y
         SD12=-BNX*GF2Y
         SD13=-BNX*GF2Z
         SD21=-BNY*GF2X
         SD22=BNZ*GF2Z+BNX*GF2X
         SD23=-BNY*GF2Z
         SD31=-BNZ*GF2X
         SD32=-BNZ*GF2Y
         SD33=BNY*GF2Y+BNX*GF2X
C
         SB(1,IVA1)=SB(1,IVA1)+
     &  (SC11*EE(1,1)+SC12*EE(2,1)+SC13*EE(3,1))*A1
         SB(1,IVA2)=SB(1,IVA2)+
     &  (SC11*EE(1,2)+SC12*EE(2,2)+SC13*EE(3,2))*A2
         SB(1,IVA3)=SB(1,IVA3)+
     &  (SC11*EE(1,3)+SC12*EE(2,3)+SC13*EE(3,3))*A3
         SB(2,IVA1)=SB(2,IVA1)+
     &  (SC21*EE(1,1)+SC22*EE(2,1)+SC23*EE(3,1))*A1
         SB(2,IVA2)=SB(2,IVA2)+
     &  (SC21*EE(1,2)+SC22*EE(2,2)+SC23*EE(3,2))*A2
         SB(2,IVA3)=SB(2,IVA3)+
     &  (SC21*EE(1,3)+SC22*EE(2,3)+SC23*EE(3,3))*A3
         SB(3,IVA1)=SB(3,IVA1)+
     &  (SC31*EE(1,1)+SC32*EE(2,1)+SC33*EE(3,1))*A1
         SB(3,IVA2)=SB(3,IVA2)+
     &  (SC31*EE(1,2)+SC32*EE(2,2)+SC33*EE(3,2))*A2
         SB(3,IVA3)=SB(3,IVA3)+
     &  (SC31*EE(1,3)+SC32*EE(2,3)+SC33*EE(3,3))*A3
C
         SB(4,IVA1)=SB(4,IVA1)+
     &  (SD11*EE(1,1)+SD12*EE(2,1)+SD13*EE(3,1))*A1
         SB(4,IVA2)=SB(4,IVA2)+
     &  (SD11*EE(1,2)+SD12*EE(2,2)+SD13*EE(3,2))*A2
         SB(4,IVA3)=SB(4,IVA3)+
     &  (SD11*EE(1,3)+SD12*EE(2,3)+SD13*EE(3,3))*A3
         SB(5,IVA1)=SB(5,IVA1)+
     &  (SD21*EE(1,1)+SD22*EE(2,1)+SD23*EE(3,1))*A1
         SB(5,IVA2)=SB(5,IVA2)+
     &  (SD21*EE(1,2)+SD22*EE(2,2)+SD23*EE(3,2))*A2
         SB(5,IVA3)=SB(5,IVA3)+
     &  (SD21*EE(1,3)+SD22*EE(2,3)+SD23*EE(3,3))*A3
         SB(6,IVA1)=SB(6,IVA1)+
     &  (SD31*EE(1,1)+SD32*EE(2,1)+SD33*EE(3,1))*A1
         SB(6,IVA2)=SB(6,IVA2)+
     &  (SD31*EE(1,2)+SD32*EE(2,2)+SD33*EE(3,2))*A2
         SB(6,IVA3)=SB(6,IVA3)+
     &  (SD31*EE(1,3)+SD32*EE(2,3)+SD33*EE(3,3))*A3
C
C ( n x F )fai
         SA(1,IVA1)=SA(1,IVA1)+
     &              (-BNZ*EE(2,1)+BNY*EE(3,1) )*FA1*A1
         SA(1,IVA2)=SA(1,IVA2)+
     &              (-BNZ*EE(2,2)+BNY*EE(3,2) )*FA1*A2
         SA(1,IVA3)=SA(1,IVA3)+
     &              (-BNZ*EE(2,3)+BNY*EE(3,3) )*FA1*A3
         SA(2,IVA1)=SA(2,IVA1)+
     &              ( BNZ*EE(1,1)-BNX*EE(3,1) )*FA1*A1
         SA(2,IVA2)=SA(2,IVA2)+
     &              ( BNZ*EE(1,2)-BNX*EE(3,2) )*FA1*A2
         SA(2,IVA3)=SA(2,IVA3)+
     &              ( BNZ*EE(1,3)-BNX*EE(3,3) )*FA1*A3
         SA(3,IVA1)=SA(3,IVA1)+
     &              (-BNY*EE(1,1)+BNX*EE(2,1) )*FA1*A1
         SA(3,IVA2)=SA(3,IVA2)+
     &              (-BNY*EE(1,2)+BNX*EE(2,2) )*FA1*A2
         SA(3,IVA3)=SA(3,IVA3)+
     &              (-BNY*EE(1,3)+BNX*EE(2,3) )*FA1*A3
C
         SA(4,IVA1)=SA(4,IVA1)+
     &              (-BNZ*EE(2,1)+BNY*EE(3,1) )*FA2*A1
         SA(4,IVA2)=SA(4,IVA2)+
     &              (-BNZ*EE(2,2)+BNY*EE(3,2) )*FA2*A2
         SA(4,IVA3)=SA(4,IVA3)+
     &              (-BNZ*EE(2,3)+BNY*EE(3,3) )*FA2*A3
         SA(5,IVA1)=SA(5,IVA1)+
     &              ( BNZ*EE(1,1)-BNX*EE(3,1) )*FA2*A1
         SA(5,IVA2)=SA(5,IVA2)+
     &              ( BNZ*EE(1,2)-BNX*EE(3,2) )*FA2*A2
         SA(5,IVA3)=SA(5,IVA3)+
     &              ( BNZ*EE(1,3)-BNX*EE(3,3) )*FA2*A3
         SA(6,IVA1)=SA(6,IVA1)+
     &              (-BNY*EE(1,1)+BNX*EE(2,1) )*FA2*A1
         SA(6,IVA2)=SA(6,IVA2)+
     &              (-BNY*EE(1,2)+BNX*EE(2,2) )*FA2*A2
         SA(6,IVA3)=SA(6,IVA3)+
     &              (-BNY*EE(1,3)+BNX*EE(2,3) )*FA2*A3
C (F * n) GRAD fai
         D1=BNX*DD(1,1)+BNY*DD(2,1)+BNZ*DD(3,1)
         D2=BNX*DD(1,2)+BNY*DD(2,2)+BNZ*DD(3,2)
         D3=BNX*DD(1,3)+BNY*DD(2,3)+BNZ*DD(3,3)
         SC(1,IVA1)=SC(1,IVA1)+GF1X*D1*A1
         SC(2,IVA1)=SC(2,IVA1)+GF1Y*D1*A1
         SC(3,IVA1)=SC(3,IVA1)+GF1Z*D1*A1
         SC(1,IVA2)=SC(1,IVA2)+GF1X*D2*A2
         SC(2,IVA2)=SC(2,IVA2)+GF1Y*D2*A2
         SC(3,IVA2)=SC(3,IVA2)+GF1Z*D2*A2
         SC(1,IVA3)=SC(1,IVA3)+GF1X*D3*A3
         SC(2,IVA3)=SC(2,IVA3)+GF1Y*D3*A3
         SC(3,IVA3)=SC(3,IVA3)+GF1Z*D3*A3
C
         SC(4,IVA1)=SC(4,IVA1)+GF2X*D1*A1
         SC(5,IVA1)=SC(5,IVA1)+GF2Y*D1*A1
         SC(6,IVA1)=SC(6,IVA1)+GF2Z*D1*A1
         SC(4,IVA2)=SC(4,IVA2)+GF2X*D2*A2
         SC(5,IVA2)=SC(5,IVA2)+GF2Y*D2*A2
         SC(6,IVA2)=SC(6,IVA2)+GF2Z*D2*A2
         SC(4,IVA3)=SC(4,IVA3)+GF2X*D3*A3
         SC(5,IVA3)=SC(5,IVA3)+GF2Y*D3*A3
         SC(6,IVA3)=SC(6,IVA3)+GF2Z*D3*A3
C
   10 CONTINUE
   11 CONTINUE
      RETURN
      END
C
      SUBROUTINE CALJC(II,EE,DD,AL1,AL2)
      PARAMETER (NPPN=1000,NEEN=2000,NVVN=3000,IFMX=2)
      implicit double precision (A-H, O-Z)
      COMMON /TR1/NP,NE3,ZZ(NPPN,3),NOE3(NEEN,6),ZNP(NEEN,3)
      COMMON /VCT/ZLN(NVVN,6),NOV(NVVN,2),NV
      COMMON /TR4/BN(NEEN,3),T1(NEEN,3),T2(NEEN,3),ZI(NEEN,3)
      COMMON /VPR0/E1XE2(NEEN),E2XE3(NEEN),E1XE3(NEEN),E2XE1(NEEN)
      COMMON /VPR1/SXYZ(NEEN),UNXE2(NEEN,3), UNXE1(NEEN,3)
      DIMENSION EE(3,3),DD(3,3)
C
      S=SXYZ(II)
      E1CE2DN=E1XE2(II)
      E2CE3DN=E2XE3(II)
      E1CE3DN=E1XE3(II)
      E2CE1DN=E2XE1(II)
      BNCE2X =UNXE2(II,1)
      BNCE2Y =UNXE2(II,2)
      BNCE2Z =UNXE2(II,3)
      BNCE1X =UNXE1(II,1)
      BNCE1Y =UNXE1(II,2)
      BNCE1Z =UNXE1(II,3)
      X1=ZZ(NOE3(II,1),1)
      Y1=ZZ(NOE3(II,1),2)
      Z1=ZZ(NOE3(II,1),3)
      X2=ZZ(NOE3(II,2),1)
      Y2=ZZ(NOE3(II,2),2)
      Z2=ZZ(NOE3(II,2),3)
      X3=ZZ(NOE3(II,3),1)
      Y3=ZZ(NOE3(II,3),2)
      Z3=ZZ(NOE3(II,3),3)
      X=AL1*X1+AL2*X2+(1.-AL1-AL2)*X3
      Y=AL1*Y1+AL2*Y2+(1.-AL1-AL2)*Y3
      Z=AL1*Z1+AL2*Z2+(1.-AL1-AL2)*Z3
      S1X=(Z2-Z3)*Y+(Y3-Y2)*Z+Z3*Y2-Y3*Z2
      S1Y=(X2-X3)*Z+(Z3-Z2)*X+X3*Z2-Z3*X2
      S1Z=(Y2-Y3)*X+(X3-X2)*Y+Y3*X2-X3*Y2
C
      SSAL1=1./(S*S*AL1*4.)
      RL1RX=(S1Y*(Z3-Z2)+S1Z*(Y2-Y3))*SSAL1
      RL1RY=(S1X*(Z2-Z3)+S1Z*(X3-X2))*SSAL1
      RL1RZ=(S1X*(Y3-Y2)+S1Y*(X2-X3))*SSAL1
C
      S2X=(Z3-Z1)*Y-(Y3-Y1)*Z-Z3*Y1+Y3*Z1
      S2Y=(X3-X1)*Z-(Z3-Z1)*X-X3*Z1+Z3*X1
      S2Z=(Y3-Y1)*X-(X3-X1)*Y-Y3*X1+Y1*X3
C
      SSAL2=1./(S*S*AL2*4.)
      RL2RX=(S2Z*(Y3-Y1)-S2Y*(Z3-Z1))*SSAL2
      RL2RY=(S2X*(Z3-Z1)-S2Z*(X3-X1))*SSAL2
      RL2RZ=(S2Y*(X3-X1)-S2X*(Y3-Y1))*SSAL2
C
      EE(1,1) = (1.-AL2)*BNCE2X/E2CE1DN
     &          -AL1*E2CE3DN*BNCE1X/E1CE3DN/E2CE1DN
      EE(1,2) = -AL2*E1CE3DN*BNCE2X/E1CE2DN/E2CE3DN+
     &          (1.-AL1)*BNCE1X/E1CE2DN
      EE(1,3) = AL2*BNCE2X/E2CE3DN+AL1*BNCE1X/E1CE3DN
      EE(2,1) = (1.-AL2)*BNCE2Y/E2CE1DN
     &          -AL1*E2CE3DN*BNCE1Y/E1CE3DN/E2CE1DN
      EE(2,2) = -AL2*E1CE3DN*BNCE2Y/E1CE2DN/E2CE3DN+
     &          (1.-AL1)*BNCE1Y/E1CE2DN
      EE(2,3) = AL2*BNCE2Y/E2CE3DN+AL1*BNCE1Y/E1CE3DN
      EE(3,1) = (1.-AL2)*BNCE2Z/E2CE1DN
     &          -AL1*E2CE3DN*BNCE1Z/E1CE3DN/E2CE1DN
      EE(3,2) = -AL2*E1CE3DN*BNCE2Z/E1CE2DN/E2CE3DN+
     &          (1.-AL1)*BNCE1Z/E1CE2DN
      EE(3,3) = AL2*BNCE2Z/E2CE3DN+AL1*BNCE1Z/E1CE3DN
C
      RU1YRX = -RL2RX*BNCE2Y/E2CE1DN
     &         -RL1RX*E2CE3DN*BNCE1Y/E1CE3DN/E2CE1DN
      RU2YRX = -RL2RX*E1CE3DN*BNCE2Y/E1CE2DN/E2CE3DN
     &         -RL1RX*BNCE1Y/E1CE2DN
      RU3YRX =  RL2RX*BNCE2Y/E2CE3DN+RL1RX*BNCE1Y/E1CE3DN
      RU1ZRX = -RL2RX*BNCE2Z/E2CE1DN
     &         -RL1RX*E2CE3DN*BNCE1Z/E1CE3DN/E2CE1DN
      RU2ZRX = -RL2RX*E1CE3DN*BNCE2Z/E1CE2DN/E2CE3DN
     &         -RL1RX*BNCE1Z/E1CE2DN
      RU3ZRX =  RL2RX*BNCE2Z/E2CE3DN+RL1RX*BNCE1Z/E1CE3DN
      RU1XRY = -RL2RY*BNCE2X/E2CE1DN
     &         -RL1RY*E2CE3DN*BNCE1X/E1CE3DN/E2CE1DN
      RU2XRY = -RL2RY*E1CE3DN*BNCE2X/E1CE2DN/E2CE3DN
     &         -RL1RY*BNCE1X/E1CE2DN
      RU3XRY =  RL2RY*BNCE2X/E2CE3DN+RL1RY*BNCE1X/E1CE3DN
      RU1ZRY = -RL2RY*BNCE2Z/E2CE1DN
     &         -RL1RY*E2CE3DN*BNCE1Z/E1CE3DN/E2CE1DN
      RU2ZRY = -RL2RY*E1CE3DN*BNCE2Z/E1CE2DN/E2CE3DN
     &         -RL1RY*BNCE1Z/E1CE2DN
      RU3ZRY =  RL2RY*BNCE2Z/E2CE3DN+RL1RY*BNCE1Z/E1CE3DN
      RU1XRZ = -RL2RZ*BNCE2X/E2CE1DN
     &         -RL1RZ*E2CE3DN*BNCE1X/E1CE3DN/E2CE1DN
      RU2XRZ = -RL2RZ*E1CE3DN*BNCE2X/E1CE2DN/E2CE3DN
     &         -RL1RZ*BNCE1X/E1CE2DN
      RU3XRZ =  RL2RZ*BNCE2X/E2CE3DN+RL1RZ*BNCE1X/E1CE3DN
      RU1YRZ = -RL2RZ*BNCE2Y/E2CE1DN
     &         -RL1RZ*E2CE3DN*BNCE1Y/E1CE3DN/E2CE1DN
      RU2YRZ = -RL2RZ*E1CE3DN*BNCE2Y/E1CE2DN/E2CE3DN
     &         -RL1RZ*BNCE1Y/E1CE2DN
      RU3YRZ =  RL2RZ*BNCE2Y/E2CE3DN+RL1RZ*BNCE1Y/E1CE3DN
      DD(1,1)=-(RU1YRZ-RU1ZRY)
      DD(1,2)=-(RU2YRZ-RU2ZRY)
      DD(1,3)=-(RU3YRZ-RU3ZRY)
      DD(2,1)=-(RU1ZRX-RU1XRZ)
      DD(2,2)=-(RU2ZRX-RU2XRZ)
      DD(2,3)=-(RU3ZRX-RU3XRZ)
      DD(3,1)=-(RU1XRY-RU1YRX)
      DD(3,2)=-(RU2XRY-RU2YRX)
      DD(3,3)=-(RU3XRY-RU3YRX)
   20 CONTINUE
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE CROS( RX,RY,RZ,V1X,V1Y,V1Z,V2X,V2Y,V2Z)
C
      implicit double precision (A-H, O-Z)
      RX=V1Y*V2Z-V1Z*V2Y
      RY=V1Z*V2X-V1X*V2Z
      RZ=V1X*V2Y-V1Y*V2X
C
      RETURN
      END
C
      SUBROUTINE DOT(R,V1X,V1Y,V1Z,V2X,V2Y,V2Z)
      implicit double precision (A-H, O-Z)
      R=V1X*V2X+V1Y*V2Y+V1Z*V2Z
      RETURN
      END
C
      SUBROUTINE ZERO
      PARAMETER (NPPN=1000,NEEN=2000,NVVN=3000,IFMX=2)
      implicit double precision (A-H, O-Z)
      COMPLEX*16 SA,SB,SC
      COMMON /TR1/NP,NE3,ZZ(NPPN,3),NOE3(NEEN,6),ZNP(NEEN,3)
      COMMON /VCT/ZLN(NVVN,6),NOV(NVVN,2),NV
      COMMON /TR6/SA(6,3*(NVVN+NVVN)),SB(6,3*(NVVN+NVVN)),
     &                                SC(6,3*(NVVN+NVVN))
      DO 20 K=1,(NV+NV)
         DO 20 J=1,6
            SA(J,K)=0.
            SB(J,K)=0.
            SC(J,K)=0.
   20 CONTINUE
      RETURN
      END
C
      SUBROUTINE CALDIA( AKAKU, II)
      PARAMETER (NPPN=1000,NEEN=2000,NVVN=3000,IFMX=2)
C
      implicit double precision (A-H, O-Z)
      COMMON /VCT/ZLN(NVVN,6),NOV(NVVN,2),NV
      COMMON /TR1/NP,NE3,ZZ(NPPN,3),NOE3(NEEN,6),ZNP(NEEN,3)
      COMMON /TR4/BN(NEEN,3),T1(NEEN,3),T2(NEEN,3),ZI(NEEN,3)
      COMMON /HSI/HOSEI
C
      II1 = 0
      II2 = 0
      DO 10 I=1,NE3
         I1 = IABS(NOE3(I,4))
         I2 = IABS(NOE3(I,5))
         I3 = IABS(NOE3(I,6))
         IF( (I1.EQ.II .OR. I2.EQ.II .OR. I3.EQ.II )
     & .AND. (II1.NE.0 .AND. II2.EQ.0) ) II2=I
C
         IF( (I1.EQ.II .OR. I2.EQ.II .OR. I3.EQ.II )
     &                              .AND. II1.EQ.0 ) II1=I
C
   10 CONTINUE
      ANAI = BN(II1,1)*BN(II2,1)+BN(II1,2)*BN(II2,2)+
     &       BN(II1,3)*BN(II2,3)
      AKAKU = (3.141592 - ACOS( ANAI ))/(2.*3.141592)
      AKAKU=AKAKU*HOSEI
      RETURN
      END

