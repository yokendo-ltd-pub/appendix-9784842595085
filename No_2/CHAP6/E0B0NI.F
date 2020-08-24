      DIMENSION ZZI(200,3),B0(200,3),A0(200,3),IO(200)
C
      R0=0.05
      NTH=1
      NR=20
C      WRITE(*,*)'KAKU HANKEI'
C      READ(*,*,END=10)NTH,NR

  10  DTH= 1. * 3.141592 / FLOAT(NTH)
      DR= 0.10 / FLOAT(NR)
      NP=0
      DO 20 I=0,NTH-1
         TH=DTH*I
         DO 30 J=0,NR
            R=DR*J
            NP=NP+1
            ZZI(NP,1) = R * COS( TH )
            ZZI(NP,3) = R * SIN( TH )
            ZZI(NP,2) = 0.
            RR=SQRT(ZZI(NP,1)**2 + ZZI(NP,2)**2 + ZZI(NP,3)**2)
            IF((R0-RR) .GT. 0.) THEN
              B0(NP,1)=0.
              B0(NP,2)=0.
              B0(NP,3)=0.
              A0(NP,1)=0.
              A0(NP,2)=0.
              A0(NP,3)=0.
              IO(NP)=1
            ELSE
              B0(NP,1)=0.
              B0(NP,2)=0.
              B0(NP,3)=1.
              A0(NP,1)= -ZZI(NP,2)*.5
              A0(NP,2)=  ZZI(NP,1)*.5
              A0(NP,3)=0.
              IO(NP)=0
            ENDIF
            IF( ABS(R0-RR)/R0*100. .LT. 1.)THEN
                NP=NP-1
            ENDIF
   30    CONTINUE
   20 CONTINUE
      DO 40 I=1,NP
         ZZI(I+NP,1)=ZZI(I,1)
         ZZI(I+NP,2)=ZZI(I,3)
         ZZI(I+NP,3)=ZZI(I,2)
   40 CONTINUE
C      NP=NP*2
      OPEN(37,FILE='f37',ACCESS='SEQUENTIAL',FORM='FORMATTED')
      WRITE(37,1)NP
      WRITE(37,2)((ZZI(I,J),J=1,3),I=1,NP)
      WRITE(37,1)(IO(I),I=1,NP)
      CLOSE(37)
      OPEN(49,FILE='f49',ACCESS='SEQUENTIAL',FORM='FORMATTED')
      WRITE(49,1)NP
      WRITE(49,2)((ZZI(I,J),J=1,3),I=1,NP)
      WRITE(49,2)((B0(I,J),J=1,3),I=1,NP)
      CLOSE(49)
      OPEN(50,FILE='f50',ACCESS='SEQUENTIAL',FORM='FORMATTED')
      WRITE(50,1)NP
      WRITE(50,2)((ZZI(I,J),J=1,3),I=1,NP)
      WRITE(50,2)((A0(I,J),J=1,3),I=1,NP)
      CLOSE(50)
    1 FORMAT(I5)
    2 FORMAT(E15.8)
      END

