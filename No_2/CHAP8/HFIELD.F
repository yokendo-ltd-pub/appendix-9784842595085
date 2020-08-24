C***********************************************************************
C*                << 3-D EDDY CURRENT ANALYSIS >>                      *
C*                       << hfield.F77 >>                              *
C***********************************************************************
C     Eddy currents are formulated by volume integral equations whose 
C     unkowns are the surface electric currents Js defined by Js=nxHs 
C     with the surface magnetic field Hs and the unit normal n. 
C
C     Setting PARAMETER 
      PARAMETER (NEO=120)
      DOUBLE PRECISION DW
      CHARACTER FNAME*12
      COMMON /VNT/ VX,VY,VZ
      COMMON /XYZO/ XXO,YYO,ZZO
      COMMON /XYZ/ X00,Y00,Z00,X10,Y10,Z10,X01,Y01,Z01,X11,Y11,Z11
      DIMENSION B(NEO*8)
      DIMENSION XYZ(NEO*2,3),NOE(NEO,4),NO(NEO)
      DIMENSION SXYZ(40,3),NSE(40),NOS(40,4)
      DIMENSION XO(NEO),YO(NEO),ZO(NEO)
      DIMENSION VN(NEO,3),VT(2,NEO,3)
C
      PI=3.14159
C
C     KSYM=0: 3-D,    KSYM=1: AXIALLY-SYMMETRICAL
                   KSYM=0
C
 6666 FORMAT('   Input file name as file.dat')
      WRITE(6,6666)
 5555 FORMAT(A12)
      READ(5,5555) FNAME
C
      OPEN(7,FILE=FNAME)
C
C     NP: Total Nodes,   NE: Total Elements
      READ(7,*) NP,NE
      WRITE(6,6000) NP,NE
 6000 FORMAT(' Total Nodes=',I5,'   Total Elements=',I5)
C
C     XYZ: Nodal Positions
      DO 1000 I=1,NP
      READ(7,*) (XYZ(I,J),J=1,3)
      DO 1000 K=1,3
C     [m]--->[cm]
      XYZ(I,K)=XYZ(I,K)*100.
 1000 CONTINUE
C
C     NO: Number of Nodes in a Element
C     NOE: Order of Number Counting in Anticlockwise way
      DO 1100 I=1,NE
      READ(7,*) NO(I),(NOE(I,J),J=1,NO(I))
 1100 CONTINUE
C
C     F...Frequency [Hz]
C     PS...relative magnetic permeability
C     D...electric conductivity [S/m]
C
      READ(7,*) F,PS,D
C
C     Identifier to know whether source is current or field
C     KSC=1: Source is field,   KSC=2: Source is current.
      READ(7,*) KSC
C
      IF(KSC.NE.1) THEN
C
C     NS: Total Nodes of Source Current Loops
C     NTURN: Total Source Current Loops
      READ(7,*) NS,NTURN
      WRITE(6,6001) NS
      WRITE(6,6002) NTURN
 6001 FORMAT(1H ,'Total Nodes of Current Loops=',I5)
 6002 FORMAT(1H ,'Total Current Loops=',I5)
C
C     SXYZ: Nodal Positions of Source Current Loop
      DO 1010 I=1,NS
      READ(7,*) (SXYZ(I,J),J=1,3)
C     [m]--->[cm]
      DO 1010 K=1,3
      SXYZ(I,K)=SXYZ(I,K)*100.
 1010 CONTINUE
C
C     NSE: Number of Nodes in a Element
C     NOE: Order of Number Counting in Anticlockwise way
      DO 1020 I=1,NTURN
      READ(7,*) NSE(I),(NOS(I,J),J=1,NSE(I))
 1020 CONTINUE
C
C     Current Flowing in a Source Current Loop
      READ(7,*) AMP
C
      ELSE 
C     HSX,HSY,HSZ: Source Magnetic Fields [A/m]
 6003 FORMAT(1H ,'Source Magnetic Fields')
 6004 FORMAT(1H ,'HSX=',E10.3,'[A/m]','  HSY=',E10.3,'[A/m]',
     &'  HSZ=',E10.3,'[A/m]')
      READ(7,*) HSXO,HSYO,HSZO
      WRITE(6,6003)
      WRITE(6,6004) HSXO,HSYO,HSZO
      ENDIF
      CLOSE(7)
C
      KEND=2-KSYM
      N1=NE*KEND
      N2=N1+N1
C
      OPEN(7,FILE='j.out',FORM='UNFORMATTED')
      READ(7) (B(J),J=1,N2)
      CLOSE (7)
C      OPEN(7,FILE='jro.out',FORM='FORMATTED')
C      WRITE(7,630) (B(J),J=1,N1)
C      CLOSE (7)
C      OPEN(7,FILE='jio.out',FORM='FORMATTED')
C      WRITE(7,630) (B(J+N1),J=1,N1)
C      CLOSE (7)
C
C     [A/m]--->[A/cm]
      DO 2410 J=1,N1
      B(J)=B(J)/100.
      B(J+N1)=B(J+N1)/100.
 2410 CONTINUE
C
C     D...electric conductivity [S/cm]
      D=D/100.
C     P0...magnetic permeability [H/cm]
      P0=4.*PI*1.E-9
C     U...inverse of skin depth [1/cm]
      WO=2.*PI*F
      U=SQRT(WO*P0*PS*D/2.)
      SDP=1./U
      WRITE(6,631) SDP
      TS=PI/U
C     TT....Conductor Thickness
      TT=5./U
C
C     XO,YO,ZO: Center of Element
      DO 1001 J=1,NE
      XO(J)=0.
      YO(J)=0.
      ZO(J)=0.
      DO 1002 I=1,NO(J)
      XO(J)=XO(J)+XYZ(NOE(J,I),1)
      YO(J)=YO(J)+XYZ(NOE(J,I),2)
      ZO(J)=ZO(J)+XYZ(NOE(J,I),3)
 1002 CONTINUE
      XO(J)=XO(J)/4.
      YO(J)=YO(J)/4.
      ZO(J)=ZO(J)/4.
C      WRITE(6,630) XO(J),YO(J),ZO(J)
 1001 CONTINUE
C
C     VT: Unit Vector of Tangential Direction
C     VN: Unit Vector of Normal Direction
      CALL VECTN(NE,NOE,XYZ,VN,VT,NEO)
C
C
   55 WRITE(6,6665)
 6665 FORMAT('   Input calculating point xo, yo, zo')
      READ(5,*) XXO,YYO,ZZO
C     [m]-->[cm]
      XXO=XXO*100.
      YYO=YYO*100.
      ZZO=ZZO*100.
C
      HXRO=0.
      HYRO=0.
      HZRO=0.
      HXIO=0.
      HYIO=0.
      HZIO=0.
C
C     OBSERVATION POINT
C
      II=0
      DO 20 K=1,KEND
      DO 20 I=1,NE
      II=II+1
C      IF(TT.LT.TS) WRITE(6,630) TT
      VX=VN(I,1)
      VY=VN(I,2)
      VZ=VN(I,3)
C
      IF(K.EQ.2) GO TO 23
      U01=XYZ(NOE(I,4),1)
      U00=XYZ(NOE(I,1),1)
      V01=XYZ(NOE(I,4),2)
      V00=XYZ(NOE(I,1),2)
      W01=XYZ(NOE(I,4),3)
      W00=XYZ(NOE(I,1),3)
      U11=XYZ(NOE(I,3),1)
      U10=XYZ(NOE(I,2),1)
      V11=XYZ(NOE(I,3),2)
      V10=XYZ(NOE(I,2),2)
      W11=XYZ(NOE(I,3),3)
      W10=XYZ(NOE(I,2),3)
      GO TO 22
   23 U01=XYZ(NOE(I,1),1)
      U00=XYZ(NOE(I,2),1)
      V01=XYZ(NOE(I,1),2)
      V00=XYZ(NOE(I,2),2)
      W01=XYZ(NOE(I,1),3)
      W00=XYZ(NOE(I,2),3)
      U11=XYZ(NOE(I,4),1)
      U10=XYZ(NOE(I,3),1)
      V11=XYZ(NOE(I,4),2)
      V10=XYZ(NOE(I,3),2)
      W11=XYZ(NOE(I,4),3)
      W10=XYZ(NOE(I,3),3)
C
   22 DL1=SQRT((U11+U10-U01-U00)**2+(V11+V10-V01-V00)**2+
     &(W11+W10-W01-W00)**2)/2.
      DL2=SQRT((U01+U11-U00-U10)**2+(V01+V11-V00-V10)**2+
     &(W01+W11-W00-W10)**2)/2.
      DL12=(DL1+DL2)/2.
      DL1O=SQRT((XO(I)-XXO)**2+(YO(I)-YYO)**2+(ZO(I)-ZZO)**2)
      NDO=1
      IF(DL1O.LE.DL12) NDO=2
C     Number of Divisions of the Element at Observation Point NDO*NDO
      ND=NDO*2-1
C
      DU31=(U10-U00)/ND
      DV31=(V10-V00)/ND
      DW31=(W10-W00)/ND
      DU42=(U11-U01)/ND
      DV42=(V11-V01)/ND
      DW42=(W11-W01)/ND
C
      DO 25 JD=1,ND
C
      U10J=DU31*JD+U00
      U00J=DU31*(JD-1)+U00
      V10J=DV31*JD+V00
      V00J=DV31*(JD-1)+V00
      W10J=DW31*JD+W00
      W00J=DW31*(JD-1)+W00
      U11J=DU42*JD+U01
      U01J=DU42*(JD-1)+U01
      V11J=DV42*JD+V01
      V01J=DV42*(JD-1)+V01
      W11J=DW42*JD+W01
      W01J=DW42*(JD-1)+W01
C
      DU12=(U01J-U00J)/ND
      DV12=(V01J-V00J)/ND
      DW12=(W01J-W00J)/ND
      DU34=(U11J-U10J)/ND
      DV34=(V11J-V10J)/ND
      DW34=(W11J-W10J)/ND
C
      DO 25 ID=1,ND
      X01=DU12*ID+U00J
      X00=DU12*(ID-1)+U00J
      Y01=DV12*ID+V00J
      Y00=DV12*(ID-1)+V00J
      Z01=DW12*ID+W00J
      Z00=DW12*(ID-1)+W00J
      X11=DU34*ID+U10J
      X10=DU34*(ID-1)+U10J
      Y11=DV34*ID+V10J
      Y10=DV34*(ID-1)+V10J
      Z11=DW34*ID+W10J
      Z10=DW34*(ID-1)+W10J
C
C     Making Matrix |A| by Magnetic Charges
C
      IF(PS.LT.1.01) GO TO 28
      CALL ABYM(U,TT,HXR,HYR,HZR,HXI,HYI,HZI)
      HXRO=HXRO+(HXR*B(II)-HXI*B(II+N1))*(PS-1.)
      HYRO=HYRO+(HYR*B(II)-HYI*B(II+N1))*(PS-1.)
      HZRO=HZRO+(HZR*B(II)-HZI*B(II+N1))*(PS-1.)
      HXIO=HXIO+(HXI*B(II)+HXR*B(II+N1))*(PS-1.)
      HYIO=HYIO+(HYI*B(II)+HYR*B(II+N1))*(PS-1.)
      HZIO=HZIO+(HZI*B(II)+HZR*B(II+N1))*(PS-1.)
C
C     Making Matrix |A| by Currents
C
   28 CALL ABYI(U,TT,HXR,HYR,HZR,HXI,HYI,HZI)
      HXRO=HXRO+HXR*B(II)-HXI*B(II+N1)
      HYRO=HYRO+HYR*B(II)-HYI*B(II+N1)
      HZRO=HZRO+HZR*B(II)-HZI*B(II+N1)
      HXIO=HXIO+HXI*B(II)+HXR*B(II+N1)
      HYIO=HYIO+HYI*B(II)+HYR*B(II+N1)
      HZIO=HZIO+HZI*B(II)+HZR*B(II+N1)
   25 CONTINUE
C
   20 CONTINUE
C
      HXRO=-HXRO
      HYRO=-HYRO
      HZRO=-HZRO
      HXIO=-HXIO
      HYIO=-HYIO
      HZIO=-HZIO
C
C     Making Matrix |B| by Currents
C
      IF(KSC.EQ.1) GO TO 310
C
C     Magnetic Field by source current
C
      DO 300 NTN=1,NTURN
      DO 304 J=1,NSE(NTN)
      J1=J+1
      IF(J.EQ.NSE(NTN)) J1=1
      X1=SXYZ(NOS(NTN,J),1)
      Y1=SXYZ(NOS(NTN,J),2)
      Z1=SXYZ(NOS(NTN,J),3)
      X2=SXYZ(NOS(NTN,J1),1)
      Y2=SXYZ(NOS(NTN,J1),2)
      Z2=SXYZ(NOS(NTN,J1),3)
C
      CALL MFBYI(XXO,YYO,ZZO,X1,Y1,Z1,X2,Y2,Z2,HX,HY,HZ)
      HXRO=HXRO+HX*AMP
      HYRO=HYRO+HY*AMP
      HZRO=HZRO+HZ*AMP
C
  304 CONTINUE
  300 CONTINUE
C
      GO TO 320
C
C     Making Matrix |B| by Magnetic Fields
C
C     Mgnetic Field by source magnetic field H [A/m]
C     [A/m]-->[A/cm]
  310 HSX=HSXO/100.
      HSY=HSYO/100.
      HSZ=HSZO/100.
C
      HXRO=HXRO+HSX
      HYRO=HYRO+HSY
      HZRO=HZRO+HSZ
C
C
C     [A/cm]-->[A/m],   [V/cm]-->[V/m]
  320 HXRO=HXRO*100.
      HYRO=HYRO*100.
      HZRO=HZRO*100.
      HXIO=HXIO*100.
      HYIO=HYIO*100.
      HZIO=HZIO*100.
      WRITE(6,632) HXRO,HYRO,HZRO
  632 FORMAT(1H ,'HXR=',E11.3,'   HYR=',E11.3,'   HZR=',E11.3)
      WRITE(6,633) HXIO,HYIO,HZIO
  633 FORMAT(1H ,'HXI=',E11.3,'   HYI=',E11.3,'   HZI=',E11.3)
C
 6664 FORMAT('   Input 0  when END,   Input 1  when CONTINUE')
      WRITE(6,6664)
      READ(5,*) IEND
C
      IF(IEND.EQ.1) GO TO 55
C
      STOP
  600 FORMAT(' MAKING MAT.',I5,'/',I4)
  601 FORMAT(' L.S. method',I5,'/',I4)
  610 FORMAT(2I5,6E11.3)
  630 FORMAT(6E11.3)
  631 FORMAT(' skin depth=',E11.3,'[cm]')
  677 FORMAT(15F5.1)
      END
C
C
      SUBROUTINE VECTN(NE,NOE,XYZ,VN,VT,NEO)
C     Obtaining Tangential and Normal Directions ar Observation Point
      DIMENSION XYZ(NEO*2,3),NOE(NEO,4),VN(NEO,3),VT(2,NEO,3)
      DO 10 J=1,NE
      X1=XYZ(NOE(J,1),1)
      X2=XYZ(NOE(J,2),1)
      X3=XYZ(NOE(J,3),1)
      X4=XYZ(NOE(J,4),1)
      X23=(X2+X3)/2.
      X14=(X1+X4)/2.
      X34=(X3+X4)/2.
      X12=(X1+X2)/2.
      Y1=XYZ(NOE(J,1),2)
      Y2=XYZ(NOE(J,2),2)
      Y3=XYZ(NOE(J,3),2)
      Y4=XYZ(NOE(J,4),2)
      Y23=(Y2+Y3)/2.
      Y14=(Y1+Y4)/2.
      Y34=(Y3+Y4)/2.
      Y12=(Y1+Y2)/2.
      Z1=XYZ(NOE(J,1),3)
      Z2=XYZ(NOE(J,2),3)
      Z3=XYZ(NOE(J,3),3)
      Z4=XYZ(NOE(J,4),3)
      Z23=(Z2+Z3)/2.
      Z14=(Z1+Z4)/2.
      Z34=(Z3+Z4)/2.
      Z12=(Z1+Z2)/2.
C
C     TANGENTIAL DIRECTION
      R1=SQRT((X23-X14)**2+(Y23-Y14)**2+(Z23-Z14)**2)
      VT(1,J,1)=(X23-X14)/R1
      VT(1,J,2)=(Y23-Y14)/R1
      VT(1,J,3)=(Z23-Z14)/R1
      R2=SQRT((X34-X12)**2+(Y34-Y12)**2+(Z34-Z12)**2)
      VT(2,J,1)=(X34-X12)/R2
      VT(2,J,2)=(Y34-Y12)/R2
      VT(2,J,3)=(Z34-Z12)/R2
C     NORMAL DIRECTION
      VN(J,1)=VT(1,J,2)*VT(2,J,3)-VT(1,J,3)*VT(2,J,2)
      VN(J,2)=VT(1,J,3)*VT(2,J,1)-VT(1,J,1)*VT(2,J,3)
      VN(J,3)=VT(1,J,1)*VT(2,J,2)-VT(1,J,2)*VT(2,J,1)
      ABVN=SQRT(VN(J,1)**2+VN(J,2)**2+VN(J,3)**2)
      VN(J,1)=VN(J,1)/ABVN
      VN(J,2)=VN(J,2)/ABVN
      VN(J,3)=VN(J,3)/ABVN
C      WRITE(6,630) J,VN(J,1),VN(J,2),VN(J,3)
C      WRITE(6,630) J,VT(1,J,1),VT(1,J,2),VT(1,J,3)
C      WRITE(6,630) J,VT(2,J,1),VT(2,J,2),VT(2,J,3)
   10 CONTINUE
      RETURN
  630 FORMAT(1H ,I5,6E12.3)
      END
C
      SUBROUTINE ABYI(U,T,HXR,HYR,HZR,HXI,HYI,HZI)
C
C     Making Matrix Element of |A| by Current
C
      COMMON /VNT/ VX,VY,VZ
      COMMON /XYZO/ XXO,YYO,ZZO
      COMMON /XYZ/ X00,Y00,Z00,X10,Y10,Z10,X01,Y01,Z01,X11,Y11,Z11
C     << BIO-SAVART FIELD BY ELECTRIC CURRENT >>
C
      PI=3.14159
      NSD=6
      NWD=4
C
      U0=(X11+X10-X01-X00)/2.
      V0=(Y11+Y10-Y01-Y00)/2.
      W0=(Z11+Z10-Z01-Z00)/2.
      D0=SQRT(U0**2+V0**2+W0**2)
      U32=X11-X10
      V32=Y11-Y10
      W32=Z11-Z10
      U41=X01-X00
      V41=Y01-Y00
      W41=Z01-Z00
C
      SX1=V0*W32-W0*V32
      SY1=W0*U32-U0*W32
      SZ1=U0*V32-V0*U32
      S1=SQRT(SX1**2+SY1**2+SZ1**2)
      SX2=V0*W41-W0*V41
      SY2=W0*U41-U0*W41
      SZ2=U0*V41-V0*U41
      S2=SQRT(SX2**2+SY2**2+SZ2**2)
C
      IF(S1.GE.S2) THEN
        IR=1
        DW=S2/D0/NWD
        X1=X00
        Y1=Y00
        Z1=Z00
        X4=X01
        Y4=Y01
        Z4=Z01
        RATIO=S2/S1
        X2=((1.+RATIO)*X10+(1.-RATIO)*X11)/2.
        Y2=((1.+RATIO)*Y10+(1.-RATIO)*Y11)/2.
        Z2=((1.+RATIO)*Z10+(1.-RATIO)*Z11)/2.
        X3=((1.+RATIO)*X11+(1.-RATIO)*X10)/2.
        Y3=((1.+RATIO)*Y11+(1.-RATIO)*Y10)/2.
        Z3=((1.+RATIO)*Z11+(1.-RATIO)*Z10)/2.
      ELSE
        IR=0
        DW=S1/D0/NWD
        X2=X10
        Y2=Y10
        Z2=Z10
        X3=X11
        Y3=Y11
        Z3=Z11
        RATIO=S1/S2
        X1=((1.+RATIO)*X00+(1.-RATIO)*X01)/2.
        Y1=((1.+RATIO)*Y00+(1.-RATIO)*Y01)/2.
        Z1=((1.+RATIO)*Z00+(1.-RATIO)*Z01)/2.
        X4=((1.+RATIO)*X01+(1.-RATIO)*X00)/2.
        Y4=((1.+RATIO)*Y01+(1.-RATIO)*Y00)/2.
        Z4=((1.+RATIO)*Z01+(1.-RATIO)*Z00)/2.
      ENDIF
      IF(RATIO.GT.1.01) STOP 100
C
      TT=T
      UT=U*TT
      IF(UT.GT.20.) TT=20./U
      TS=PI/U
      IF(T.LT.TS) TS=TT
      UT=U*TT
      DT=TS/NSD
C
      EC1=EXP(UT)*COS(UT)
      EC2=EXP(-UT)*COS(-UT)
      ES1=EXP(UT)*SIN(UT)
      ES2=EXP(-UT)*SIN(-UT)
      CR=EC1-EC2
      CI=ES1-ES2
      C2=CR**2+CI**2
      CR=CR/C2
      CI=-CI/C2
C
      HXR=0.
      HYR=0.
      HZR=0.
      HXI=0.
      HYI=0.
      HZI=0.
C
      DO 30 ID=1,NSD
      SD=DT*(ID-.5)
      SDX=VX*SD
      SDY=VY*SD
      SDZ=VZ*SD
      USD=U*(TT-SD)
      ER=(EXP(USD)*COS(USD)+EXP(-USD)*COS(-USD))
      EI=(EXP(USD)*SIN(USD)+EXP(-USD)*SIN(-USD))
      EXR=ER*CR-EI*CI
      EXI=ER*CI+EI*CR
C
      CALL HBYI(NWD,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     &SDX,SDY,SDZ,DW,HX,HY,HZ)
C
      HXR=HXR+HX*U*(EXR-EXI)*DT
      HYR=HYR+HY*U*(EXR-EXI)*DT
      HZR=HZR+HZ*U*(EXR-EXI)*DT
      HXI=HXI+HX*U*(EXR+EXI)*DT
      HYI=HYI+HY*U*(EXR+EXI)*DT
      HZI=HZI+HZ*U*(EXR+EXI)*DT
   30 CONTINUE
C
      IF(T.LT.TS) GO TO 300
C
      USD=U*(TT-TS)
      ER=(EXP(USD)*COS(USD)-EXP(-USD)*COS(-USD))
      EI=(EXP(USD)*SIN(USD)-EXP(-USD)*SIN(-USD))
      EXR=ER*CR-EI*CI
      EXI=ER*CI+EI*CR
C
      SDX=VX*TS
      SDY=VY*TS
      SDZ=VZ*TS
C
C     Magnetic Field by Current on Rectangular Surface
      CALL HBYI(NWD,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     &SDX,SDY,SDZ,DW,HX,HY,HZ)
C
      HXR=HXR+HX*EXR
      HYR=HYR+HY*EXR
      HZR=HZR+HZ*EXR
      HXI=HXI+HX*EXI
      HYI=HYI+HY*EXI
      HZI=HZI+HZ*EXI
C
  300 IF(RATIO.GT.0.98) RETURN
C
      DO 4000 K=1,2
C
      IF(IR.EQ.1) THEN
        PN=1.
      IF(K.EQ.1) THEN
        U1O=X00
        V1O=Y00
        W1O=Z00
        U2O=X10
        V2O=Y10
        W2O=Z10
        U3O=X2
        V3O=Y2
        W3O=Z2
      ELSE
C
        U1O=X01
        V1O=Y01
        W1O=Z01
        U2O=X11
        V2O=Y11
        W2O=Z11
        U3O=X3
        V3O=Y3
        W3O=Z3
      ENDIF
C
      ELSE
        PN=-1.
      IF(K.EQ.1) THEN
        U1O=X10
        V1O=Y10
        W1O=Z10
        U2O=X00
        V2O=Y00
        W2O=Z00
        U3O=X1
        V3O=Y1
        W3O=Z1
C
      ELSE
        U1O=X11
        V1O=Y11
        W1O=Z11
        U2O=X01
        V2O=Y01
        W2O=Z01
        U3O=X4
        V3O=Y4
        W3O=Z4
      ENDIF
      ENDIF
C
        U31=U3O-U1O
        V31=V3O-V1O
        W31=W3O-W1O
        D0=SQRT(U31**2+V31**2+W31**2)
        U32=U3O-U2O
        V32=V3O-V2O
        W32=W3O-W2O
C
        SX1=V31*W32-W31*V32
        SY1=W31*U32-U31*W32
        SZ1=U31*V32-V31*U32
        S1=SQRT(SX1**2+SY1**2+SZ1**2)
        DW=S1/D0/NWD
C
      DO 40 ID=1,NSD
      SD=DT*(ID-.5)
      SDX=VX*SD
      SDY=VY*SD
      SDZ=VZ*SD
      USD=U*(TT-SD)
      ER=(EXP(USD)*COS(USD)+EXP(-USD)*COS(-USD))
      EI=(EXP(USD)*SIN(USD)+EXP(-USD)*SIN(-USD))
      EXR=ER*CR-EI*CI
      EXI=ER*CI+EI*CR
C
C     Magnetic Field by Current on Triangular Surface
      CALL HBYIT(NWD,U1O,V1O,W1O,U2O,V2O,W2O,U3O,V3O,W3O,
     &SDX,SDY,SDZ,DW,PN,HX,HY,HZ)
C
      HXR=HXR+HX*U*(EXR-EXI)*DT
      HYR=HYR+HY*U*(EXR-EXI)*DT
      HZR=HZR+HZ*U*(EXR-EXI)*DT
      HXI=HXI+HX*U*(EXR+EXI)*DT
      HYI=HYI+HY*U*(EXR+EXI)*DT
      HZI=HZI+HZ*U*(EXR+EXI)*DT
   40 CONTINUE
C
      IF(T.LT.TS) RETURN
C
      USD=U*(TT-TS)
      ER=(EXP(USD)*COS(USD)-EXP(-USD)*COS(-USD))
      EI=(EXP(USD)*SIN(USD)-EXP(-USD)*SIN(-USD))
      EXR=ER*CR-EI*CI
      EXI=ER*CI+EI*CR
C
      SDX=VX*TS
      SDY=VY*TS
      SDZ=VZ*TS
C
C     Magnetic Field by Current on Triangular Surface
      CALL HBYIT(NWD,U1O,V1O,W1O,U2O,V2O,W2O,U3O,V3O,W3O,
     &SDX,SDY,SDZ,DW,PN,HX,HY,HZ)
C
      HXR=HXR+HX*EXR
      HYR=HYR+HY*EXR
      HZR=HZR+HZ*EXR
      HXI=HXI+HX*EXI
      HYI=HYI+HY*EXI
      HZI=HZI+HZ*EXI
C
 4000 CONTINUE
C
      RETURN
  600 FORMAT(1H ,10I5)
  610 FORMAT(1H ,2I5,6E11.3)
      END
C
      SUBROUTINE HBYI(NWD,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     &SDX,SDY,SDZ,DW,HX,HY,HZ)
C
C     Magnetic Field by Current on Rectangular Surface
C
      COMMON /XYZO/ XXO,YYO,ZZO
C
      HX=0.
      HY=0.
      HZ=0.
C
      U1=X1-SDX
      V1=Y1-SDY
      W1=Z1-SDZ
      U2=X2-SDX
      V2=Y2-SDY
      W2=Z2-SDZ
      U3=X3-SDX
      V3=Y3-SDY
      W3=Z3-SDZ
      U4=X4-SDX
      V4=Y4-SDY
      W4=Z4-SDZ
      DU32=(U3-U2)/NWD
      DV32=(V3-V2)/NWD
      DW32=(W3-W2)/NWD
      DU41=(U4-U1)/NWD
      DV41=(V4-V1)/NWD
      DW41=(W4-W1)/NWD
C
      DO 31 JD=1,NWD
C      WRITE(6,600) ID,JD
      XX1=U1+DU41*(JD-.5)
      YY1=V1+DV41*(JD-.5)
      ZZ1=W1+DW41*(JD-.5)
      XX2=U2+DU32*(JD-.5)
      YY2=V2+DV32*(JD-.5)
      ZZ2=W2+DW32*(JD-.5)
C
      CALL MFBYI(XXO,YYO,ZZO,XX1,YY1,ZZ1,XX2,YY2,ZZ2,DHX,DHY,DHZ)
      HX=HX+DHX*DW
      HY=HY+DHY*DW
      HZ=HZ+DHZ*DW
C      WRITE(6,630) XX1,YY1,ZZ1
C      WRITE(6,630) XX2,YY2,ZZ2
   31 CONTINUE
C      WRITE(6,630) XXO,YYO,ZZO
C      WRITE(6,630) HX,HY,HZ
C
  630 FORMAT(1H ,6E11.3)
      RETURN
      END
C
      SUBROUTINE HBYIT(NWD,U1O,V1O,W1O,U2O,V2O,W2O,U3O,V3O,W3O,
     &SDX,SDY,SDZ,DW,PN,HX,HY,HZ)
C
C     Magnetic Field by Current on Triangular Surface
C
      COMMON /XYZO/ XXO,YYO,ZZO
C
      HX=0.
      HY=0.
      HZ=0.
C
      U1=U1O-SDX
      V1=V1O-SDY
      W1=W1O-SDZ
      U2=U2O-SDX
      V2=V2O-SDY
      W2=W2O-SDZ
      U3=U3O-SDX
      V3=V3O-SDY
      W3=W3O-SDZ
      DU32=(U2-U3)/NWD
      DV32=(V2-V3)/NWD
      DW32=(W2-W3)/NWD
      DU12=(U2-U1)/NWD
      DV12=(V2-V1)/NWD
      DW12=(W2-W1)/NWD
C
      DO 41 JD=1,NWD
      XX1=U1+DU12*(JD-.5)
      YY1=V1+DV12*(JD-.5)
      ZZ1=W1+DW12*(JD-.5)
      XX2=U3+DU32*(JD-.5)
      YY2=V3+DV32*(JD-.5)
      ZZ2=W3+DW32*(JD-.5)
C
      CALL MFBYI(XXO,YYO,ZZO,XX1,YY1,ZZ1,XX2,YY2,ZZ2,DHX,DHY,DHZ)
      HX=HX+DHX*DW*PN
      HY=HY+DHY*DW*PN
      HZ=HZ+DHZ*DW*PN
   41 CONTINUE
C
      RETURN
  630 FORMAT(1H ,6E11.3)
      END
C
      SUBROUTINE ABYM(U,T,HXR,HYR,HZR,HXI,HYI,HZI)
C
C     Making Matrix Element of |A| by Magnetic Charge
C     Coulomb Field by Magnetic Charge
C
      COMMON /VNT/ VX,VY,VZ
      COMMON /XYZO/ XXO,YYO,ZZO
      COMMON /XYZ/ X00,Y00,Z00,X10,Y10,Z10,X01,Y01,Z01,X11,Y11,Z11
      DIMENSION SP(4,3)
C
      PI=3.14159
      NSD=6
C
      U0=(X11+X10-X01-X00)/2.
      V0=(Y11+Y10-Y01-Y00)/2.
      W0=(Z11+Z10-Z01-Z00)/2.
      RC=SQRT(U0**2+V0**2+W0**2)
C
      SP(1,1)=X00
      SP(2,1)=X10
      SP(3,1)=X11
      SP(4,1)=X01
      SP(1,2)=Y00
      SP(2,2)=Y10
      SP(3,2)=Y11
      SP(4,2)=Y01
      SP(1,3)=Z00
      SP(2,3)=Z10
      SP(3,3)=Z11
      SP(4,3)=Z01
C
      TT=T
      UT=U*TT
      IF(UT.GT.20.) TT=20./U
      TS=PI/U
      IF(T.LT.TS) TS=TT
      UT=U*TT
      DT=TS/NSD
C
      EC1=EXP(UT)*COS(UT)
      EC2=EXP(-UT)*COS(-UT)
      ES1=EXP(UT)*SIN(UT)
      ES2=EXP(-UT)*SIN(-UT)
      CR=EC1-EC2
      CI=ES1-ES2
      C2=CR**2+CI**2
      CR=CR/C2
      CI=-CI/C2
C
      HXR=0.
      HYR=0.
      HZR=0.
      HXI=0.
      HYI=0.
      HZI=0.
C
      DO 40 J=1,4
C
      J1=J+1
      IF(J1.EQ.5) J1=1
      X1=SP(J,1)
      Y1=SP(J,2)
      Z1=SP(J,3)
      X2=SP(J1,1)
      Y2=SP(J1,2)
      Z2=SP(J1,3)
C
      U21=X2-X1
      V21=Y2-Y1
      W21=Z2-Z1
      RS=SQRT(U21**2+V21**2+W21**2)
      CS210=(U0*U21+V0*V21+W0*W21)/RS/RC
C
      DO 30 ID=1,NSD
      SD=DT*(ID-.5)
      SDX=VX*SD
      SDY=VY*SD
      SDZ=VZ*SD
      USD=U*(TT-SD)
      ER=(EXP(USD)*COS(USD)-EXP(-USD)*COS(-USD))
      EI=(EXP(USD)*SIN(USD)-EXP(-USD)*SIN(-USD))
      EXR=ER*CR-EI*CI
      EXI=ER*CI+EI*CR
      XX1=X1-SDX
      YY1=Y1-SDY
      ZZ1=Z1-SDZ
      XX2=X2-SDX
      YY2=Y2-SDY
      ZZ2=Z2-SDZ
C
C     Magnetic Field by Line Magnetic Charge
      CALL MFBYM(XXO,YYO,ZZO,XX1,YY1,ZZ1,XX2,YY2,ZZ2,DHX,DHY,DHZ)
      HXR=HXR+DHX*EXR*DT*CS210
      HYR=HYR+DHY*EXR*DT*CS210
      HZR=HZR+DHZ*EXR*DT*CS210
      HXI=HXI+DHX*EXI*DT*CS210
      HYI=HYI+DHY*EXI*DT*CS210
      HZI=HZI+DHZ*EXI*DT*CS210
   30 CONTINUE
C
C
      IF(T.LT.TS) GO TO 40
C
      USD=U*(TT-TS)
      ER=(EXP(USD)*COS(USD)+EXP(-USD)*COS(-USD))-2.
      EI=(EXP(USD)*SIN(USD)+EXP(-USD)*SIN(-USD))
      EXR=ER*CR-EI*CI
      EXI=ER*CI+EI*CR
C
      SDX=VX*TS
      SDY=VY*TS
      SDZ=VZ*TS
      XX1=X1-SDX
      YY1=Y1-SDY
      ZZ1=Z1-SDZ
      XX2=X2-SDX
      YY2=Y2-SDY
      ZZ2=Z2-SDZ
C
      CALL MFBYM(XXO,YYO,ZZO,XX1,YY1,ZZ1,XX2,YY2,ZZ2,DHX,DHY,DHZ)
      HXR=HXR+DHX*(EXR+EXI)/2./U*CS210
      HYR=HYR+DHY*(EXR+EXI)/2./U*CS210
      HZR=HZR+DHZ*(EXR+EXI)/2./U*CS210
      HXI=HXI+DHX*(EXI-EXR)/2./U*CS210
      HYI=HYI+DHY*(EXI-EXR)/2./U*CS210
      HZI=HZI+DHZ*(EXI-EXR)/2./U*CS210
   40 CONTINUE
      RETURN
  610 FORMAT(2I5,6E11.3)
      END
C
      SUBROUTINE MFBYI(X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,HX,HY,HZ)
C
C     Obtaining Magnetic Field by Straight Line Current
C     Bio-Savart Field by Current
C
C     X0,Y0,Z0...Observation Point
C     X1,Y1,Z1...One End Position of Line Current
C     X2,Y2,Z2...Another End Position of Line Current
C     HX,HY,HZ...Magnetic Field by Current
C
      PI=3.14159
      DDLO=1.E-10
      X21=X2-X1
      Y21=Y2-Y1
      Z21=Z2-Z1
      DL=SQRT(X21*X21+Y21*Y21+Z21*Z21)
      IF(DL.LE.DDLO) GO TO 1000
      DDL=DL*1.E-6
C
      XXL=X21/DL
      YYL=Y21/DL
      ZZL=Z21/DL
C
      XXO=X0-X1
      YYO=Y0-Y1
      ZZO=Z0-Z1
C
      A=XXL**2+YYL**2+ZZL**2
      B=-2.*(XXO*XXL+YYO*YYL+ZZO*ZZL)
      C=XXO**2+YYO**2+ZZO**2
      ABC=4.*A*C-B*B
      IF(ABC.LE.DDL) GO TO 1000
C
      HH=2.*(2.*A*DL+B)/SQRT(A*DL*DL+B*DL+C)-2.*B/SQRT(C)
      HH=HH/4./PI/DL/(4.*A*C-B*B)
      HX=HH*(Y21*ZZO-Z21*YYO)
      HY=HH*(Z21*XXO-X21*ZZO)
      HZ=HH*(X21*YYO-Y21*XXO)
      RETURN
 1000 WRITE(6,6000)
      WRITE(6,600) X0,Y0,Z0
      WRITE(6,600) X1,Y1,Z1,X2,Y2,Z2
      STOP
  600 FORMAT(1H ,6E11.3)
 6000 FORMAT(1H ,'Position data are something wrong!!!')
      END
C
      SUBROUTINE MFBYM(X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,HX,HY,HZ)
C
C     Obtaining Magnetic Field by Straight Line Current
C     Bio-Savart Field by Current
C
C     X0,Y0,Z0...Observation Point
C     X1,Y1,Z1...One End Position of Line Current
C     X2,Y2,Z2...Another End Position of Line Current
C     HX,HY,HZ...Magnetic Field by Current
C
      PI=3.14159
      DDLO=1.E-10
      X21=X2-X1
      Y21=Y2-Y1
      Z21=Z2-Z1
      DL=SQRT(X21*X21+Y21*Y21+Z21*Z21)
      IF(DL.EQ.DDLO) GO TO 1000
      DDL=DL*1.E-6
C
      XXL=X21/DL
      YYL=Y21/DL
      ZZL=Z21/DL
C
      XXO=X0-X1
      YYO=Y0-Y1
      ZZO=Z0-Z1
C
      A=XXL**2+YYL**2+ZZL**2
      B=-2.*(XXO*XXL+YYO*YYL+ZZO*ZZL)
      C=XXO**2+YYO**2+ZZO**2
      ABC=4.*A*C-B*B
      IF(ABC.LE.DDL) GO TO 1000
C
      HH=2.*(2.*A*DL+B)/SQRT(A*DL*DL+B*DL+C)-2.*B/SQRT(C)
      HH=HH/4./PI/(4.*A*C-B*B)
      HX=HH*XXO
      HY=HH*YYO
      HZ=HH*ZZO
C
      HH=2.*(B*DL+2*C)/SQRT(A*DL*DL+B*DL+C)-4.*C/SQRT(C)
      HH=HH/4./PI/(4.*A*C-B*B)
      HX=HX+HH*XXL
      HY=HY+HH*YYL
      HZ=HZ+HH*ZZL
      RETURN
 1000 WRITE(6,6000)
      WRITE(6,600) X0,Y0,Z0
      WRITE(6,600) X1,Y1,Z1,X2,Y2,Z2
      STOP
  600 FORMAT(1H ,6E11.3)
 6000 FORMAT(1H ,'Position data are something wrong!!!')
      END
