      CHARACTER*20 FN
      DIMENSION NOE(10000,3)
      DIMENSION Z(10000,3),VR(10000,3),VI(10000,3)
C
      WRITE(*,*)'File Format Converter'
      WRITE(*,*)'1:Triangular File & Vector FIle'
      WRITE(*,*)'2:Vector File'
      Write(*,*)'?'
      READ(*,*)IYN
      if( IYN .EQ. 1 )THEN
         CALL GETFNAME( 'File Name(Triangular)?                 ', FN)
         OPEN(10,FILE=FN,ACCESS='SEQUENTIAL',FORM='FORMATTED')
         READ(10,6) NP,NE
         DO 100 I=1,NP
 100        READ(10,7) (Z(I,J),J=1,3)
         DO 200 I=1,NE
 200     READ(10,6) (NOE(I,J),J=1,3)
         CLOSE(10)
         CALL GETFNAME( 'OutputName(Triangular)?                ', FN)
         OPEN(30,FILE=FN,ACCESS='SEQUENTIAL',FORM='FORMATTED')
         WRITE(30,1)NP,NE
         WRITE(30,2)((Z(I,J),J=1,3),I=1,NP)
         WRITE(30,1)((NOE(I,J),J=1,3),I=1,NE)
         CLOSE(30)
      ENDIF
      CALL GETFNAME( 'File Name(Vector)?                   ', FN)
      OPEN(40,FILE=FN,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      READ(40,1)NP
      READ(40,2)((Z(I,J),J=1,3),I=1,NP)
      READ(40,'(E15.8)')((VR(I,J),VI(I,J),J=1,3),I=1,NP)
      CLOSE(40)
      CALL GETFNAME( 'OutputName?                         ', FN)
      OPEN(20,FILE=FN,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      WRITE(20,3)NP
      DO 10 I=1,NP
         WRITE(20,4)(Z(I,J),J=1,3)
 10   CONTINUE
      DO 20 I=1,NP
         DO 30 J=1,3
            WRITE(20,5)VR(I,J),VI(I,J)
 30      CONTINUE
 20   CONTINUE
      CLOSE(20)
 1    FORMAT(I6)
 2    FORMAT(E15.8)
 3    FORMAT(I5)
 4    FORMAT(E15.8)
 5    FORMAT(E15.8)
 6    FORMAT(3I5)
 7    FORMAT(3F10.3)
      END
C
      SUBROUTINE GETFNAME( MESG , FN)
      CHARACTER*40 MESG
      CHARACTER*20 FN
      WRITE(*,'(1H ,A20)') MESG
      READ(*,'(A20)')FN
      RETURN
      END


