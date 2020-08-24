      DIMENSION SP(122,3),NOE(120,4),NO(120)
C
      PI=3.14159
      R=.05
      NT=12
      NZ=6
      NT1=NT+1
      NZ1=NZ+1
      DT=2.*PI/NT
      DZ=PI/NZ1
C
      NP=0
      DO 10 J=1,NT1
      JJ=J-1
      DO 10 I=1,NZ1
      II=I-1
      NP=NP+1
      SP(NP,1)=R*SIN(DZ*(II+.5))*COS(DT*JJ)
      SP(NP,2)=R*SIN(DZ*(II+.5))*SIN(DT*JJ)
      SP(NP,3)=R*COS(DZ*(II+.5))
   10 CONTINUE
C
      NE=0
      DO 110 J=1,NT
      DO 110 I=1,NZ
      NE=NE+1
      NO(NE)=4
      NOE(NE,1)=I+1+NZ1*(J-1)
      NOE(NE,2)=I+1+NZ1*J
      NOE(NE,3)=I+NZ1*J
      NOE(NE,4)=I+NZ1*(J-1)
  110 CONTINUE
C
C
      F=50.
      PS=200.
      D=6485000.
      P0=4.*PI*1.E-7
      KSC=1
      HSX=0.
      HSY=0.
      HSZ=1./P0
C
C
      OPEN(7,FILE='SPHERE.DAT')
      WRITE(7,*) NP,NE
      DO 100 I=1,NP
      WRITE(7,*) (SP(I,J),J=1,3)
  100 CONTINUE
      DO 200 I=1,NE
      WRITE(7,*) NO(I),(NOE(I,J),J=1,NO(I))
  200 CONTINUE
C
C
      WRITE(7,*) F,PS,D
      WRITE(7,*) KSC
      WRITE(7,*) HSX,HSY,HSZ
      CLOSE (7)
C
      STOP
  600 FORMAT(1H ,10I5)
  610 FORMAT(1H ,I5,3E12.3)
      END
