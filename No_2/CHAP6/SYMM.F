C
C MAKE VECTOR ELEMENT for tangential conponent of unknown vector
C
C     SYMMETRY CONDITION VERSION
C
      PARAMETER (NPPN=4000,NEEN=8000,NVVN=12000)
C
      implicit double precision ( A-H, O-Z)
      COMMON /TR1/NP,NE3,ZZ(NPPN,3),NOE3(NEEN,6),ZNP(NEEN,3)
      COMMON /TR2/NOV(NVVN,2),NV
      COMMON /TR3/NPOLD,NE3OLD,NVOLD,IISYM(8,3),ITBL(16,NVVN)
      COMMON /TR4/IESYM(8,3),IBSYM(8,3)
C
      CHARACTER*20 FN7
C
      FN7='f17'
      WRITE(*,*)'INPUT FILE:'
      READ(*,'(A20)',END=110)
C
  110 OPEN(17,FILE=FN7,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     &                 status='unknown')
      rewind 17
      READ(17,100)NP,NE3
      READ(17,200)((ZZ(I,J),J=1,3),I=1,NP)
      write(*,*)'aaaaaa'
      READ(17,100)((NOE3(I,J),J=1,3),I=1,NE3)
      CLOSE(17)
  100 FORMAT(I6)
  200 FORMAT(E15.8)
 300  FORMAT(3I5)
 400  FORMAT(3F10.3)
C
      NXYZ=1
      CALL DUP(NXYZ)
C
      NV=1
      NOV(NV,1)=NOE3(1,1)
      NOV(NV,2)=NOE3(1,2)
      DO 10 I=1,NE3
C
         I1=NOE3(I,1)
         I2=NOE3(I,2)
         I3=NOE3(I,3)
         IH=0
         DO 20 J=1,NV
            IF( (I1.EQ.NOV(J,1)).AND.(I2.EQ.NOV(J,2)) ) IH=IH+1
            IF( (I1.EQ.NOV(J,2)).AND.(I2.EQ.NOV(J,1)) ) IH=IH+1
   20    CONTINUE
         IF( IH .EQ. 0 )THEN
            NV = NV +1
            NOV(NV,1)=I1
            NOV(NV,2)=I2
         ENDIF
         IH=0
         DO 30 J=1,NV
            IF( (I2.EQ.NOV(J,1)).AND.(I3.EQ.NOV(J,2))) IH=IH+1
            IF( (I2.EQ.NOV(J,2)).AND.(I3.EQ.NOV(J,1))) IH=IH+1
   30    CONTINUE
         IF( IH .EQ. 0 )THEN
            NV = NV +1
            NOV(NV,1)=I2
            NOV(NV,2)=I3
         ENDIF
         IH=0
         DO 40 J=1,NV
            IF( (I3.EQ.NOV(J,1)).AND.(I1.EQ.NOV(J,2))) IH=IH+1
            IF( (I3.EQ.NOV(J,2)).AND.(I1.EQ.NOV(J,1))) IH=IH+1
   40    CONTINUE
         IF( IH .EQ. 0 )THEN
            NV = NV +1
            NOV(NV,1)=I1
            NOV(NV,2)=I3
         ENDIF
C
   10 CONTINUE
C
      DO 50 I=1,NE3
         I1=NOE3(I,1)
         I2=NOE3(I,2)
         I3=NOE3(I,3)
         DO 60 J=1,NV
            IF((I3.EQ.NOV(J,1)).AND.(I1.EQ.NOV(J,2)))NOE3(I,4)= J
            IF((I3.EQ.NOV(J,2)).AND.(I1.EQ.NOV(J,1)))NOE3(I,4)=-J
            IF((I3.EQ.NOV(J,1)).AND.(I2.EQ.NOV(J,2)))NOE3(I,5)= J
            IF((I3.EQ.NOV(J,2)).AND.(I2.EQ.NOV(J,1)))NOE3(I,5)=-J
            IF((I1.EQ.NOV(J,1)).AND.(I2.EQ.NOV(J,2)))NOE3(I,6)= J
            IF((I1.EQ.NOV(J,2)).AND.(I2.EQ.NOV(J,1)))NOE3(I,6)=-J
   60    CONTINUE
   50 CONTINUE
      WRITE(*,*)'NP,NE3=',NP,NE3
C
      CALL MKTBL1(NXYZ)
C
      OPEN(17,FILE=FN7,ACCESS='SEQUENTIAL',FORM='FORMATTED')
C IJK 1:OLD Format Else:NEW Format
      IJK=0
      IF(IJK.EQ.1)THEN
         WRITE(17,100)NP,NE3
         WRITE(17,200)((ZZ(I,J),J=1,3),I=1,NP)
         WRITE(17,100)((NOE3(I,J),J=1,3),I=1,NE3)
      ELSE
         WRITE(17,300)NP,NE3
         DO 80 I=1,NP
            WRITE(17,400)(ZZ(I,J),J=1,3)
 80      CONTINUE
         DO 70 I=1,NE3
            WRITE(17,300)(NOE3(I,J),J=1,3)
 70      CONTINUE
      ENDIF
C for the vector element
      WRITE(17,100)((NOE3(I,J),J=4,6),I=1,NE3)
      WRITE(17,100)NV
      WRITE(17,100)((NOV(I,J),J=1,2),I=1,NV)
      WRITE(17,100)NXYZ,NPOLD,NE3OLD,NVOLD
      WRITE(17,100)((IISYM(I,J),J=1,3),I=1,8)
      WRITE(17,100)((ITBL(I,J),J=1,NV),I=1,NXYZ)
      WRITE(17,100)((ITBL(I,J),J=1,NE3),I=9,8+NXYZ)
      WRITE(17,100)((IBSYM(I,J),J=1,3),I=1,8)
      WRITE(17,100)((IESYM(I,J),J=1,3),I=1,8)
      CLOSE(17)
C
      WRITE(10,*)'NP=',NP,'NE3=',NE3,'NV=',NV,'SYM cond=',NXYZ
      WRITE(10,*)'(NP,NE,NV) in the region 0',NPOLD,NE3OLD,NVOLD
      END
C
      SUBROUTINE DUP( NSYM )
C
      PARAMETER (NPPN=4000,NEEN=8000,NVVN=12000)
C
      implicit double precision ( A-H, O-Z)
      COMMON /TR1/NP,NE3,ZZ(NPPN,3),NOE3(NEEN,6),ZNP(NEEN,3)
      COMMON /TR2/NOV(NVVN,2),NV
      COMMON /TR3/NPOLD,NE3OLD,NVOLD,IISYM(8,3),ITBL(16,NVVN)
      COMMON /TR4/IESYM(8,3),IBSYM(8,3)
C
      DATA ((IISYM(I,J),I=1,8),J=1,3)/
     & +1,-1,-1,+1,+1,-1,-1,+1,
     & +1,+1,-1,-1,+1,+1,-1,-1,
     & +1,+1,+1,+1,-1,-1,-1,-1/
C
      DATA ((IBSYM(I,J),I=1,8),J=1,3)/
     & +1,-1,-1,+1,-1,+1,+1,-1,
     & +1,+1,-1,-1,-1,-1,+1,+1,
     & +1,+1,+1,+1,+1,+1,+1,+1/
C
      DATA ((IESYM(I,J),I=1,8),J=1,3)/
     & +1,+1,-1,-1,+1,+1,-1,-1,
     & +1,-1,-1,+1,+1,-1,-1,+1,
     & +1,-1,+1,-1,-1,+1,-1,+1/
C
      NPOLD=NP
      NE3OLD=NE3
C
      NEWNP=NP
      DO 10 NN=2,NSYM
         WRITE(*,*)'NP SYM=',NN
         DO 20 I=1,NP
            NEWNP=NEWNP+1
            ZZ(NEWNP,1)=ZZ(I,1)*IISYM(NN,1)
            ZZ(NEWNP,2)=ZZ(I,2)*IISYM(NN,2)
            ZZ(NEWNP,3)=ZZ(I,3)*IISYM(NN,3)
            IDUP=0
            DO 25 J=1,NN-1
               ALX=ZZ(NEWNP,1)-ZZ(I,1)*IISYM(J,1)
               ALY=ZZ(NEWNP,2)-ZZ(I,2)*IISYM(J,2)
               ALZ=ZZ(NEWNP,3)-ZZ(I,3)*IISYM(J,3)
               AL1=SQRT(ALX**2+ALY**2+ALZ**2)
               IF(AL1.LT.1.E-5)IDUP=1
   25       CONTINUE
            IF( IDUP .EQ. 1 ) NEWNP=NEWNP-1
   20   CONTINUE
   10 CONTINUE
      NP=NEWNP
C
      NEWNE=NE3
      DO 30 NN=2,NSYM
         WRITE(*,*)'NE3 SYM=',NN
         DO 40 I=1,NE3
            I1=NOE3(I,1)
            I2=NOE3(I,2)
            I3=NOE3(I,3)
            O1X=ZZ(I1,1)*IISYM(NN,1)
            O1Y=ZZ(I1,2)*IISYM(NN,2)
            O1Z=ZZ(I1,3)*IISYM(NN,3)
            O2X=ZZ(I2,1)*IISYM(NN,1)
            O2Y=ZZ(I2,2)*IISYM(NN,2)
            O2Z=ZZ(I2,3)*IISYM(NN,3)
            O3X=ZZ(I3,1)*IISYM(NN,1)
            O3Y=ZZ(I3,2)*IISYM(NN,2)
            O3Z=ZZ(I3,3)*IISYM(NN,3)
            CALL FIND(O1X,O1Y,O1Z,NEW1)
            CALL FIND(O2X,O2Y,O2Z,NEW2)
            CALL FIND(O3X,O3Y,O3Z,NEW3)
            NEWNE=NEWNE+1
            IF(IISYM(NN,1)*IISYM(NN,2)*IISYM(NN,3) .NE. -1)THEN
               NOE3(NEWNE,1)=NEW1
               NOE3(NEWNE,2)=NEW2
               NOE3(NEWNE,3)=NEW3
            ELSE
               NOE3(NEWNE,1)=NEW1
               NOE3(NEWNE,2)=NEW3
               NOE3(NEWNE,3)=NEW2
            ENDIF
            ITBL( 8+NN, NEWNE) = I
C           WRITE(*,*)NEWNE,I
            IF( NEW1*NEW2*NEW3 .EQ.0)THEN
               NEWNE=NEWNE-1
            ENDIF
   40    CONTINUE
   30 CONTINUE
      NE3=NEWNE
      RETURN
      END
C
      SUBROUTINE FIND( OX, OY, OZ, IRTN )
      PARAMETER (NPPN=4000,NEEN=8000,NVVN=12000)
C
      implicit double precision ( A-H, O-Z)
      COMMON /TR1/NP,NE3,ZZ(NPPN,3),NOE3(NEEN,6),ZNP(NEEN,3)
C
      IRTN=0
      DO 10 I=1,NP
         ALX=ZZ(I,1)-OX
         ALY=ZZ(I,2)-OY
         ALZ=ZZ(I,3)-OZ
         AL1=SQRT(ALX**2+ALY**2+ALZ**2)
         IF(AL1.LT.1.E-5)THEN
            IRTN=I
            RETURN
         ENDIF
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE MKTBL1( NSYM )
C
      PARAMETER (NPPN=4000,NEEN=8000,NVVN=12000)
C
      implicit double precision ( A-H, O-Z)
      COMMON /TR1/NP,NE3,ZZ(NPPN,3),NOE3(NEEN,6),ZNP(NEEN,3)
      COMMON /TR2/NOV(NVVN,2),NV
      COMMON /TR3/NPOLD,NE3OLD,NVOLD,IISYM(8,3),ITBL(16,NVVN)
C
      NVOLD=0
      DO 1 I=1,NE3OLD
         DO 1 J=1,3
            IF( IABS(NOE3(I,J+3)).GT.NVOLD)THEN
                NVOLD=IABS(NOE3(I,J+3))
            ENDIF
    1 CONTINUE
      WRITE(*,*)'NVOLD=',NVOLD,'  NV=',NV
      WRITE(*,*)
C
      DO 35 I=1,NVOLD
   35    ITBL(1,I)=I
C
      DO 10 I=1,NV
         I1=NOV(I,1)
         I2=NOV(I,2)
         WRITE(*,'(1H+,A15,I4,A1,I4,I4,A1)'
     &                      )'NV search=',I,'(',I1,I2,')'
         DO 20 J=1,NVOLD
            NFLG=0
            J1=NOV(J,1)
            J2=NOV(J,2)
            DO 30 NN=1,NSYM
               B1X=ZZ(J1,1)*IISYM(NN,1)-ZZ(I1,1)
               B1Y=ZZ(J1,2)*IISYM(NN,2)-ZZ(I1,2)
               B1Z=ZZ(J1,3)*IISYM(NN,3)-ZZ(I1,3)
               BB1=SQRT( B1X**2 + B1Y**2 + B1Z**2 )
               B2X=ZZ(J2,1)*IISYM(NN,1)-ZZ(I2,1)
               B2Y=ZZ(J2,2)*IISYM(NN,2)-ZZ(I2,2)
               B2Z=ZZ(J2,3)*IISYM(NN,3)-ZZ(I2,3)
               BB2=SQRT( B2X**2 + B2Y**2 + B2Z**2 )
               IF( BB1.LT.1.E-5 .AND. BB2.LT.1.E-5) THEN
                 JJ1=J
               ELSE
                 JJ1=0
               ENDIF
               B1X=ZZ(J1,1)*IISYM(NN,1)-ZZ(I2,1)
               B1Y=ZZ(J1,2)*IISYM(NN,2)-ZZ(I2,2)
               B1Z=ZZ(J1,3)*IISYM(NN,3)-ZZ(I2,3)
               BB1=SQRT( B1X**2 + B1Y**2 + B1Z**2 )
               B2X=ZZ(J2,1)*IISYM(NN,1)-ZZ(I1,1)
               B2Y=ZZ(J2,2)*IISYM(NN,2)-ZZ(I1,2)
               B2Z=ZZ(J2,3)*IISYM(NN,3)-ZZ(I1,3)
               BB2=SQRT( B2X**2 + B2Y**2 + B2Z**2 )
               IF( BB1.LT.1.E-5 .AND. BB2.LT.1.E-5) THEN
                 JJ2=J
               ELSE
                 JJ2=0
               ENDIF
               IF( ((JJ1.NE.0) .OR. (JJ2.NE.0)).AND.
     &                                      (NFLG.EQ.0))THEN
                  ITBL(NN,I)=J
                  NFLG=1
               ENDIF
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
      RETURN
      END

