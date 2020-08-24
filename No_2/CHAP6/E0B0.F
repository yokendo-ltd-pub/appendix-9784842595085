C
      PARAMETER (NPPN=8000,NEEN=6000)
      implicit double precision (A-H,O-Z)
      COMMON /TR1/NP,NE3,ZZ(NPPN,3),NOE3(NEEN,6),ZNP(NEEN,3)
      COMMON /B00/B0(3*NEEN,3),A0(3*NEEN,3)
      COMMON /VCT/ZLN(3*NEEN,6),NOV(NEEN*3,2),NV
      REAL ZI(3*NEEN,3)
      CHARACTER*20 FNAME
C
      WRITE(*,*)'FILE 17='
      FNAME='f17'
      READ(*,'(A20)',END=111)FNAME
  111 OPEN(17,FILE=FNAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     & FORM='FORMATTED')
      READ(17,3)NP,NE3
      DO 100 I=1,NP
 100     READ(17,4)(ZZ(I,J),J=1,3)
      DO 200 I=1,NE3
 200     READ(17,3)(NOE3(I,J),J=1,3)
C
      READ(17,1)((NOE3(I,J),J=4,6),I=1,NE3)
      READ(17,1)NV
      READ(17,1)((NOV(I,J),J=1,2),I=1,NV)
C
      CLOSE(17)
      DO 10 I=1,NE3
         I1=NOE3(I,1)
         I2=NOE3(I,2)
         I3=NOE3(I,3)
         ZI(I,1)=(ZZ(I1,1)+ZZ(I2,1)+ZZ(I3,1))*0.33333333
         ZI(I,2)=(ZZ(I1,2)+ZZ(I2,2)+ZZ(I3,2))*0.33333333
         ZI(I,3)=(ZZ(I1,3)+ZZ(I2,3)+ZZ(I3,3))*0.33333333
         B0(I,1)=0.
         B0(I,2)=0.
         B0(I,3)=1.
         A0(I,1)= -ZI(I,2)*.5
         A0(I,2)=  ZI(I,1)*.5
         A0(I,3)=    .0
   10 CONTINUE
      DO 20 I=1,NV
         I1=NOV(I,1)
         I2=NOV(I,2)
         A0(I+NE3,1)=-( (ZZ(I2,2)+ZZ(I1,2))*.5 )*.5
         A0(I+NE3,2)= ( (ZZ(I2,1)+ZZ(I1,1))*.5 )*.5
         A0(I+NE3,3)=   .0
         B0(I+NE3,1)=0.
         B0(I+NE3,2)=0.
         B0(I+NE3,3)=1.
   20 CONTINUE
C
      WRITE(*,*)'FILE 39='
      FNAME='f39'
      READ(*,'(A20)',END=222)FNAME
  222 OPEN(39,FILE=FNAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     & FORM='FORMATTED')
      WRITE(39,1) NE3
      WRITE(39,2)((ZI(I,J),J=1,3),I=1,NE3)
      WRITE(39,2)((B0(I,J),J=1,3),I=1,NE3)
      WRITE(39,2)((B0(I+NE3,J),J=1,3),I=1,NV)
      WRITE(*,*)'FILE 40='
      FNAME='f40'
      READ(*,'(A20)',END=333)FNAME
  333 OPEN(40,FILE=FNAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     & FORM='FORMATTED')
      WRITE(40,1) NE3
      WRITE(40,2)((ZI(I,J),J=1,3),I=1,NE3)
      WRITE(40,2)((A0(I,J),J=1,3),I=1,NE3)
      WRITE(40,2)((A0(I+NE3,J),J=1,3),I=1,NV)
    1 FORMAT(I6)
    2 FORMAT(E15.8)
    3 FORMAT(3I5)
    4 FORMAT(3F10.3)
      END

