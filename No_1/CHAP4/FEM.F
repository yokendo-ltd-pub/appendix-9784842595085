***************************************************************
      PROGRAM FEM
***************************************************************
      INCLUDE'FEM.H'
      CHARACTER FNM(5)*30,DUM*1
         INTEGER NZU(LA),NZL(LA),IU(LB,LA),IL(LB,LA)
         COMPLEX AU(LB,LA),AD(LA),AL(LB,LA)
         COMPLEX B(LA),X(LA)
         COMPLEX UU(LB,LA),UL(LB,LA),UD(LA)
         COMPLEX R(LA),P(LA),RU(LA),AP(LA)
         INTEGER IV(LA)

********** OPEN DATA FILES ****************
      OPEN(1,FILE='FEM.FIL',STATUS='UNKNOWN')
         READ(1,501)DUM
         READ(1,502)FNM(1)
         READ(1,501)DUM
         READ(1,502)FNM(2)
         READ(1,501)DUM
         READ(1,502)FNM(3)
         READ(1,501)DUM
         READ(1,502)FNM(4)
  501 FORMAT(A1)
  502 FORMAT(A30)
      CLOSE(1)

      WRITE(*,*)'***** EDDY CURRENT *******'
      WRITE(*,601)1,FNM(1)
      WRITE(*,*)'***** FLUX ***************'
      WRITE(*,601)2,FNM(2)
      WRITE(*,*)'***** POTENTIAL **********'
      WRITE(*,601)3,FNM(3)
      WRITE(*,*)'***** IMPEDANCE **********'
      WRITE(*,601)4,FNM(4)
      WRITE(*,*)' '
  601 FORMAT('  NO.',I1,'  DATA FILE = ',A30,)

      OPEN(10,FILE=FNM(1),STATUS='UNKNOWN')
      OPEN(11,FILE=FNM(2),STATUS='UNKNOWN')
      OPEN(12,FILE=FNM(3),STATUS='UNKNOWN')
      OPEN(13,FILE=FNM(4),STATUS='UNKNOWN')
     
      CALL INPUT (NEQ,FREQ,TURN,MODEL,RL,ZL,A0,Z0)
      CALL CHANGE(IMAX,5)
      PAI=3.141592654
      OMG=2.*PAI*FREQ
      NB=LB 
      NEQ=NPOIN*4
      EPS=1.E-3
      DO 10 KLM=1,NELEM
   10    CALL MAKEMAT (IMAX,KLM,NB,ZMU(KLM),SIGMA(KLM))
      DO 20 I=1,NCL
   20    CALL  MAKEVEC (IMAX,NNCL(I))
      CALL CONVEQ (NEQ,NB,NZU,IU,AU,AD,B)
      CALL ICDATA (NEQ,NB,NZU,IU,AU,AD) 
      CALL BOUND (NEQ,NC,B,AD,NZU,MODEL)
      CALL ICCG (NEQ,NB,NZU,NZL,IU,IL,AU,AL,AD,B,EPS,X,
     *                 UU,UL,UD,R,P,RU,AP,IV,LA,LB)
      CALL EDDY (X,0,0,0)
      CALL FLUX (X)
      CALL IMPEDANCE (1,X,A0,Z0,
     *                RL,ZL,MODEL,IMAX,NEQ,TURN)
      CALL  OUTPUT (NEQ,X)
      CALL  OUTPUTI (1)

      CLOSE (10)
      CLOSE (11)
      CLOSE (12)
      CLOSE (13)

      STOP
      END

***************************************************************
      SUBROUTINE INPUT (NEQ,FREQ,TURN,MODEL,RL,ZL,A0,Z0)
***************************************************************
      INCLUDE'FEM.H'
*****  OUTPUT TO FEA *****
      OPEN (60,FILE='MESH.DAT',STATUS='UNKNOWN')
      OPEN (65,FILE='DOMAIN.DAT',STATUS='UNKNOWN')

         READ(60,710)NPOIN,NELEM
         CALL ERR(N1,NELEM,' ERR NELEM')
         CALL ERR(N2,NPOIN,' ERR NPOIN') 
         DO 110 I=1,NPOIN
         READ(60,720) XYZ(I,1),XYZ(I,2),XYZ(I,3)
  110    CONTINUE
         DO 120 I=1,NELEM
            READ(60,710) (NOD(I,J),J=1,8)
  120    CONTINUE
      CLOSE (60,STATUS='KEEP')

  710 FORMAT (10I5)
  720 FORMAT (6E10.3)
      
      READ(65,610) NCL,NMTL,MTLP,NDEF,NBD,NEQ,MODEL
      CALL ERR(N3,NCL,' ERR NNCL')
      CALL ERR(N3,NDEF,' ERR NNDEF ')
      CALL ERR(N4,NMLT,' ERR NNMTL ')
      CALL ERR(N5,MTLP,' ERR NMTLP ')
      CALL ERR(N6,NBD,' ERR NNBD ')
      
      READ(65,600) FREQ,CUR1(1),TURN,RL,ZL,A0,Z0
****** COIL ******* 
      DO 10 I=1,NCL,12
         READ(65,620) (NNCL(I+J),J=0,11)
   10 CONTINUE
      DO 15 I=1,NCL,2
         READ(65,600) (CUR(I+J,1),CUR(I+J,2),CUR(I+J,3),J=0,1)
   15 CONTINUE
****** METAL ******
      DO 20 I=1,NMTL,12
         READ(65,620) (NNMTL(I+J),J=0,11)
   20 CONTINUE      
      DO 21 I=1,MTLP,12
         READ(65,620) (NMTLP(I+J),J=0,11)
   21 CONTINUE
      DO 22 I=1,NELEM,6
         READ(65,600) (ZMU(I+J),J=0,5)
   22 CONTINUE
      DO 23 I=1,NELEM,6
         READ(65,600) (SIGMA(I+J),J=0,5)
   23 CONTINUE
****** DEFECT ******    
      DO 30 I=1,NDEF,12
         READ(65,620) (NNDEF(I+J),J=0,11)
   30 CONTINUE
****** BOUNDARY *******
      DO 40 I=1,NBD,12
         READ(65,620) (NNBD(I+J),J=0,11)
   40 CONTINUE 
      DO 45 I=1,NBD,6
         READ(65,600) (VNBD(I+J),J=0, 5)
   45 CONTINUE
      
      CLOSE (65,STATUS='KEEP')

  600 FORMAT(6E12.5)
  610 FORMAT(10I7)
  620 FORMAT(12I5)

      RETURN
      END

**************************************************************
      SUBROUTINE ERR (N,M,MESSG)
**************************************************************
      CHARACTER*20 MESSG

      IF(M.GT.N) THEN
         WRITE(*,900)MESSG,M
         STOP
      END IF
  900 FORMAT(A20,I5)

      RETURN
      END
      
**************************************************************
      SUBROUTINE  MAKEMAT (IMAX,KLM,NB,ZMIU,SIGM)
**************************************************************
      INCLUDE'FEM.H'

      DO 200 I=1,8
  200    MND(I)=NOD(KLM,I)

      DO 100 I=1,8
      DO 100 J=1,8
         N11(I,J)=0. 
         N22(I,J)=0.
         N3A(I,J)=0.
         N3B(I,J)=0.
         N3C(I,J)=0.
         N1D(I,J)=0.
  100 CONTINUE
         
      DO 10 L=1, IMAX
      DO 10 M=1, IMAX
      DO 10 N=1, IMAX
        CALL SHAPEFUN(GAUSS(L),GAUSS(M),GAUSS(N))
        CALL CJACOB(KLM, CJACOBD)
        G=GCOEF(L)*GCOEF(M)*GCOEF(N)*CJACOBD
        DO 20 I=1, 8
          DO 20 K=1, 8
            DO 30 J=1, 3
              N11(I,K)=N11(I,K)+DERIVC(J,K)*DERIVC(J,I)*G
  30        CONTINUE
              N22(I,K)=N22(I,K)+SHAPE(I)*SHAPE(K)*G
  20    CONTINUE 
 
        DO 40 I=1, 8
        DO 40 K=1, 8
           N3A(I,K) =N3A(I,K) +SHAPE(I)*DERIVC(1,K)*G
           N3B(I,K) =N3B(I,K) +SHAPE(I)*DERIVC(2,K)*G
           N3C(I,K) =N3C(I,K) +SHAPE(I)*DERIVC(3,K)*G
  40    CONTINUE
  10  CONTINUE
      DO 60 I=1,8
      DO 60 K=1, 8
            N1D(I,K)= N11(I,K)*SIGM
            N3A(I,K)= N3A(I,K)*SIGM
            N3B(I,K)= N3B(I,K)*SIGM
            N3C(I,K)= N3C(I,K)*SIGM
            N22(I,K)= N22(I,K)*SIGM
            N11(I,K)= N11(I,K) /ZMIU
  60  CONTINUE
 
      CALL  ADDMAT (SIGM,NB)
 
      RETURN
      END


**************************************************************
      SUBROUTINE  MAKEVEC (IMAX,KLM)
**************************************************************
      INCLUDE'FEM.H'

      DO 200 I=1,8
  200    MND(I)=NOD(KLM,I)
      DO 10 L=1, IMAX
      DO 10 M=1, IMAX
      DO 10 N=1, IMAX
         CALL SHAPEFUN(GAUSS(L),GAUSS(M),GAUSS(N))
         CALL CJACOB(KLM, CJACOBD)
         G=GCOEF(L)*GCOEF(M)*GCOEF(N)*CJACOBD
         DO 20 I=1,NCL
            IF(NNCL(I).EQ.KLM) THEN
               DO 30 K=1,3
               NP=NPOIN*(K-1)
               DO 30 J=1,8
                  VEC(MND(J)+NP)=VEC(MND(J)+NP)+SHAPE(J)*G*CUR(I,K)
  30           CONTINUE
            END IF
  20  CONTINUE
  10  CONTINUE

      RETURN
      END

**************************************************************
      SUBROUTINE  ADDMAT (SIGM,NB)
**************************************************************
      INCLUDE'FEM.H'
         COMPLEX S,Z,C

      Z=( 0. , 0. )
***** N11+JW*N22      
      DO 10 I=1,8
      DO 10 J=1,8
      DO 10 L=1,3
         MNDI=MND(I)+NPOIN*(L-1)
         MNDJ=MND(J)+NPOIN*(L-1)
         S=CMPLX(N11(I,J),OMG*N22(I,J))
                  IF(MNDJ.EQ.MNDI)  MTXD(MNDI)=MTXD(MNDI)+S 
         IF(S.EQ.Z) GOTO 100
           
            IF(MNDJ.GT.MNDI) THEN
            DO 1 K=1,NB
               IF(IELM(K,MNDI).EQ.MNDJ) THEN
                  MTXU(K,MNDI)=MTXU(K,MNDI)+S
                  GOTO 100
               END IF
    1       CONTINUE
            DO 2 K=1,NB
               IF(IELM(K,MNDI).EQ.0) THEN
                  IELM(K,MNDI)=MNDJ
                  MTXU(K,MNDI)=MTXU(K,MNDI)+S
                  GOTO 100
               END IF
    2       CONTINUE
         END IF
  100    CONTINUE
   10 CONTINUE

************ CONDUCTOR *************
      IF(SIGM.EQ.0.0) GOTO 80
      IF(OMG.EQ.0.0) GOTO 80
      IF(MTLP.EQ.0) GOTO 80
      NPOIN2=NPOIN*2
      NPOIN3=NPOIN*3
***** N3A
      DO 20 I=1,8
      DO 20 J=1,8
            MNDI=MND(I)
            MNDJ=MND(J)+NPOIN3
            S=N3A(I,J) 
            IF(S.EQ.Z) GOTO 120
            DO 21 K=1,NB
               IF(IELM(K,MNDI).EQ.MNDJ) THEN
                  MTXU(K,MNDI)=MTXU(K,MNDI)+S
                  GOTO 120
               END IF
   21       CONTINUE
            DO 22 K=1,NB
               IF(IELM(K,MNDI).EQ.0) THEN
                  IELM(K,MNDI)=MNDJ
                  MTXU(K,MNDI)=MTXU(K,MNDI)+S
                  GOTO 120
               END IF
   22       CONTINUE
  120     CONTINUE
   20     CONTINUE
***** N3B
      DO 30 I=1,8
      DO 30 J=1,8
            MNDI=MND(I)+NPOIN
            MNDJ=MND(J)+NPOIN3
            S=N3B(I,J) 
            IF(S.EQ.Z) GOTO 130
            DO 31 K=1,NB
               IF(IELM(K,MNDI).EQ.MNDJ) THEN
                  MTXU(K,MNDI)=MTXU(K,MNDI)+S
                  GOTO 130
               END IF
   31       CONTINUE
            DO 32 K=1,NB
               IF(IELM(K,MNDI).EQ.0) THEN
                  IELM(K,MNDI)=MNDJ
                  MTXU(K,MNDI)=MTXU(K,MNDI)+S
                  GOTO 130
               END IF
   32       CONTINUE
  130     CONTINUE
   30     CONTINUE
***** N3C
      DO 40 I=1,8
      DO 40 J=1,8
            MNDI=MND(I)+NPOIN2
            MNDJ=MND(J)+NPOIN3
            S=N3C(I,J) 
            IF(S.EQ.Z) GOTO 140
            DO 41 K=1,NB
               IF(IELM(K,MNDI).EQ.MNDJ) THEN
                  MTXU(K,MNDI)=MTXU(K,MNDI)+S
                  GOTO 140
               END IF
   41       CONTINUE
            DO 42 K=1,NB
               IF(IELM(K,MNDI).EQ.0) THEN
                  IELM(K,MNDI)=MNDJ
                  MTXU(K,MNDI)=MTXU(K,MNDI)+S
                  GOTO 140
               END IF
   42       CONTINUE
  140     CONTINUE
   40     CONTINUE
   81 CONTINUE

***** N1D
      C=1/(CMPLX(0,1)*OMG)
      DO 50 I=1,8
      DO 50 J=1,8
         MNDI=MND(I)+NPOIN3
         MNDJ=MND(J)+NPOIN3
         S=C*N1D(I,J)
         IF(S.EQ.Z) GOTO 150
         IF(MNDJ.EQ.MNDI)  MTXD(MNDI)=MTXD(MNDI)+S            
            IF(MNDJ.GT.MNDI) THEN
            DO 51 K=1,NB
               IF(IELM(K,MNDI).EQ.MNDJ) THEN
                  MTXU(K,MNDI)=MTXU(K,MNDI)+S
                  GOTO 150
               END IF
   51       CONTINUE
            DO 52 K=1,NB
               IF(IELM(K,MNDI).EQ.0) THEN
                  IELM(K,MNDI)=MNDJ
                  MTXU(K,MNDI)=MTXU(K,MNDI)+S
                  GOTO 150
               END IF
   52       CONTINUE
         END IF
  150    CONTINUE

   50 CONTINUE
   80 CONTINUE

      RETURN
      END

**************************************************************
      SUBROUTINE SHAPEFUN(XI, ETA, ZETA)
**************************************************************
      INCLUDE'FEM.H'

      REAL XI,ETA,ZETA
      REAL A, XIP, ETAP,ZETAP,XIM,ETAM,ZETAM

      A    =0.125E0
      XIP  =1.0E0+XI
      ETAP =1.0E0+ETA
      ZETAP=1.0E0+ZETA
      XIM  =1.0E0-XI
      ETAM =1.0E0-ETA
      ZETAM=1.0E0-ZETA

      SHAPE(1)= A*XIM*ETAM*ZETAP
      SHAPE(2)= A*XIP*ETAM*ZETAP
      SHAPE(3)= A*XIP*ETAP*ZETAP
      SHAPE(4)= A*XIM*ETAP*ZETAP
      SHAPE(5)= A*XIM*ETAM*ZETAM
      SHAPE(6)= A*XIP*ETAM*ZETAM
      SHAPE(7)= A*XIP*ETAP*ZETAM
      SHAPE(8)= A*XIM*ETAP*ZETAM

      DERIV(1,1)= -A*ETAM*ZETAP
      DERIV(1,2)=  A*ETAM*ZETAP
      DERIV(1,3)=  A*ETAP*ZETAP
      DERIV(1,4)= -A*ETAP*ZETAP
      DERIV(1,5)= -A*ETAM*ZETAM
      DERIV(1,6)=  A*ETAM*ZETAM
      DERIV(1,7)=  A*ETAP*ZETAM
      DERIV(1,8)= -A*ETAP*ZETAM

      DERIV(2,1)= -A*XIM*ZETAP
      DERIV(2,2)= -A*XIP*ZETAP
      DERIV(2,3)=  A*XIP*ZETAP
      DERIV(2,4)=  A*XIM*ZETAP
      DERIV(2,5)= -A*XIM*ZETAM
      DERIV(2,6)= -A*XIP*ZETAM
      DERIV(2,7)=  A*XIP*ZETAM
      DERIV(2,8)=  A*XIM*ZETAM

      DERIV(3,1)=  A*XIM*ETAM
      DERIV(3,2)=  A*XIP*ETAM
      DERIV(3,3)=  A*XIP*ETAP
      DERIV(3,4)=  A*XIM*ETAP
      DERIV(3,5)= -A*XIM*ETAM
      DERIV(3,6)= -A*XIP*ETAM
      DERIV(3,7)= -A*XIP*ETAP
      DERIV(3,8)= -A*XIM*ETAP

      RETURN
      END

**************************************************************
      SUBROUTINE CJACOB(KLM, CJACOBD)
**************************************************************
       INCLUDE'FEM.H'

      DO 10 I= 1, 3
      DO 10 J= 1, 3
         JACOB(I,J)=0.0D0
  10  CONTINUE

      DO 20 I=1, 8
         K=NOD(KLM,I)
         DO 30 J=1, 3
            JACOB(J,1)=JACOB(J,1)+DERIV(J,I)*XYZ(K,1)
            JACOB(J,2)=JACOB(J,2)+DERIV(J,I)*XYZ(K,2)
            JACOB(J,3)=JACOB(J,3)+DERIV(J,I)*XYZ(K,3)
  30     CONTINUE
  20  CONTINUE

      CJACOBD=  JACOB(1,1)*JACOB(2,2)*JACOB(3,3)
     &        +JACOB(1,2)*JACOB(2,3)*JACOB(3,1)
     &        +JACOB(2,1)*JACOB(3,2)*JACOB(1,3)
     &        -JACOB(1,3)*JACOB(2,2)*JACOB(3,1)
     &        -JACOB(1,2)*JACOB(2,1)*JACOB(3,3)
     &        -JACOB(1,1)*JACOB(2,3)*JACOB(3,2)
      CJACOBD =ABS(CJACOBD)

      JACOBI(1,1)=JACOB(2,2)*JACOB(3,3)-JACOB(2,3)*JACOB(3,2)
      JACOBI(1,2)=JACOB(1,3)*JACOB(3,2)-JACOB(1,2)*JACOB(3,3)
      JACOBI(1,3)=JACOB(1,2)*JACOB(2,3)-JACOB(1,3)*JACOB(2,2)
      JACOBI(2,1)=JACOB(2,3)*JACOB(3,1)-JACOB(2,1)*JACOB(3,3)
      JACOBI(2,2)=JACOB(1,1)*JACOB(3,3)-JACOB(1,3)*JACOB(3,1)
      JACOBI(2,3)=JACOB(1,3)*JACOB(2,1)-JACOB(1,1)*JACOB(2,3)
      JACOBI(3,1)=JACOB(2,1)*JACOB(3,2)-JACOB(2,2)*JACOB(3,1)
      JACOBI(3,2)=JACOB(1,2)*JACOB(3,1)-JACOB(1,1)*JACOB(3,2)
      JACOBI(3,3)=JACOB(1,1)*JACOB(2,2)-JACOB(1,2)*JACOB(2,1)

      IF(CJACOBD.EQ.0.) THEN
         WRITE(6,*)'ERR   ELEMENT VOLUME = 0.'
         STOP
      END IF

      DO 40 I= 1, 3
      DO 40 J= 1, 3
         JACOBI(I,J)=JACOBI(I,J)/CJACOBD
   40 CONTINUE

      DO 50 I=1, 8
      DO 50 J=1, 3
         DERIVC(J,I)= JACOBI(J,1)*DERIV(1,I)+JACOBI(J,2)*DERIV(2,I)
     &               +JACOBI(J,3)*DERIV(3,I)
   50 CONTINUE

      RETURN
      END

**************************************************************
      SUBROUTINE CHANGE(KMAX, N)
**************************************************************
      INCLUDE'FEM.H'

      KMAX=N
      IF (N.EQ.3)  GOTO 103
      IF (N.EQ.5)  GOTO 105
      IF (N.EQ.7)  GOTO 107
      IF (N.EQ.9)  GOTO 109

      WRITE(*,*) 'ERROR'
      STOP

  103 CONTINUE
      GAUSS(1)= 0.774596669241483E0
      GAUSS(2)= 0.0E0
      GAUSS(3)=-0.774596669241483E0
      GCOEF(1)= 0.555555555555556E0
      GCOEF(2)= 0.888888888888889E0
      GCOEF(3)= 0.555555555555556E0
      RETURN

  105 CONTINUE
      GAUSS(1)= 0.906179845938664E0
      GAUSS(2)= 0.538469310105683E0
      GAUSS(3)= 0.0D0
      GAUSS(4)=-0.538469310105683E0
      GAUSS(5)=-0.906179845938664E0

      GCOEF(1)= 0.236926885056189E0
      GCOEF(2)= 0.478628670499366E0
      GCOEF(3)= 0.568888888888889E0
      GCOEF(4)= 0.478628670499366E0
      GCOEF(5)= 0.236926885056189E0
      RETURN

  107 CONTINUE
      GAUSS(1)= 0.949107912342759E0
      GAUSS(2)= 0.741531185599394E0
      GAUSS(3)= 0.405845151377397E0
      GAUSS(4)= 0.0D0
      GAUSS(5)=-0.405845151377397E0
      GAUSS(6)=-0.741531185599394E0
      GAUSS(7)=-0.949107912342759E0

      GCOEF(1)= 0.129484966168870E0
      GCOEF(2)= 0.279705391489277E0
      GCOEF(3)= 0.381830050505119E0
      GCOEF(4)= 0.417959183673469E0
      GCOEF(5)= 0.381830050505119E0
      GCOEF(6)= 0.279705391489277E0
      GCOEF(7)= 0.129484966168870E0
      RETURN

  109 CONTINUE
      GAUSS(1)= 0.968160239507626E0
      GAUSS(2)= 0.836031107326636E0
      GAUSS(3)= 0.613371432700590E0
      GAUSS(4)= 0.324253423403809E0
      GAUSS(5)= 0.0D0
      GAUSS(6)=-0.324253423403809E0
      GAUSS(7)=-0.613371432700590E0
      GAUSS(8)=-0.836031107326636E0
      GAUSS(9)=-0.968160239507626E0

      GCOEF(1)= 0.081274388361574E0
      GCOEF(2)= 0.180648160694857E0
      GCOEF(3)= 0.260610696402935E0
      GCOEF(4)= 0.312347077040003E0
      GCOEF(5)= 0.330239355001260E0
      GCOEF(6)= 0.312347077040003E0
      GCOEF(7)= 0.260610696402935E0
      GCOEF(8)= 0.180648160694857E0
      GCOEF(9)= 0.081274388361574E0
      RETURN

      END

**************************************************************
      SUBROUTINE CONVEQ (NEQ,NB,NZU,IU,AU,AD,B)
**************************************************************
      INCLUDE'FEM.H'
      INTEGER IU(LB,LA)
      COMPLEX AU(LB,LA),AD(LA),B(LA)

      DO 10 I=1,NPOIN*4
         AD(I)=MTXD(I)
         B(I)=VEC(I)
         DO 10 J=1,LB
	    AU(J,I)=MTXU(J,I)
	    IU(J,I)=IELM(J,I)
   10 CONTINUE

      RETURN
      END 

**************************************************************
      SUBROUTINE ICDATA (NEQ,NB,NZU,IU,AU,AD)
**************************************************************
      INCLUDE'FEM.H'
      INTEGER NZU(LA),IU(LB,LA)
      COMPLEX AU(LB,LA),AD(LA)
      DIMENSION LL(LB),AA(LB)
      COMPLEX Z,AA

      Z=( 0. , 0. )
      DO 70 I=1,LB 
         LL(I)=0
         AA(I)=Z
   70 CONTINUE

      NN=0
      DO 10 I=1,NEQ
         MM=0
         DO 20 J=1,NB
            IF(AU(J,I).EQ.Z) IU(J,I)=0
               IF(IU(J,I).NE.0) MM=MM+1
   20    CONTINUE
         NZU(I)=MM
         IF(MM.GT.NN) NN=MM
   10 CONTINUE

      IF(NN*2.4.GT.LB) THEN
         WRITE(*,*) 'LB IN ICDATA',REAL(NN*2.4)/LB
         STOP
      END IF

      DO 30 I=1,NEQ
         DO 40 J=1,NB
         M=100000000
         N=0
            DO 50 K=1,NB
               IF(IU(K,I).GT.0.AND.IU(K,I).LE.M)  THEN
                  M=IU(K,I)
                  N=K
               END IF
   50       CONTINUE 
            IF(N.EQ.0) GOTO 45
            LL(J)=IU(N,I)
            AA(J)=AU(N,I)
               IU(N,I)=0
               AU(N,I)=Z
   40    CONTINUE
   45    CONTINUE
         DO 60 K=1,NB
            IU(K,I)=LL(K)
               LL(K)=0
            AU(K,I)=AA(K)
               AA(K)=Z
   60    CONTINUE
   30 CONTINUE
      
      RETURN
      END

**************************************************************
      SUBROUTINE BOUND (NEQ,NC,B,AD,NZU,MODEL)
**************************************************************
      INCLUDE'FEM.H'
      INTEGER NZU(LA)
      COMPLEX AD(LA),B(LA),Z

      PENA=1.0E 03
      Z=(0. , 0. )
      DO 10 I=1,NEQ
         IF(AD(I).NE.Z) GOTO 10
         AD(I)=(1.0E 1 ,0.0 )
   10 CONTINUE

      DO 20 I=1,NBD
         AD(NNBD(I))=PENA*AD(NNBD(I))
         B (NNBD(I))=AD(NNBD(I))*VNBD(I)
   20 CONTINUE
           
      RETURN
      END
      
**************************************************************
      SUBROUTINE ICCG (NEQ,NB,NZU,NZL,IU,IL,AU,AL,AD,B,EPS,X,
     *                 UU,UL,UD,R,P,RU,AP,IV,LA,LB)
**************************************************************
      INTEGER NZU(LA),NZL(LA),IU(LB,LA),IL(LB,LA)
      COMPLEX AU(LB,LA),AD(LA),AL(LB,LA)
      COMPLEX B(LA),X(LA)
      COMPLEX UU(LB,LA),UL(LB,LA),UD(LA)
      COMPLEX R(LA),P(LA),RU(LA),AP(LA)
      INTEGER IV(LA)
      
      LAB=LA*LB
      NINT=(LA*3+LAB*2)*4
      NCMLX=(LA*8+LAB*4)*8
      WRITE(*,*)'  SINGLE PRECISION',   (NINT+NCMLX)/1.E6,' MBYTES'
      WRITE(*,*)'  DOUBLE PRECISION', 2*(NINT+NCMLX)/1.E6,' MBYTES'
         
       BALL=0.
       DO 1 I=1,NEQ
         BALL=BALL+CABS(B(I))
    1 CONTINUE   
      IF(BALL.EQ.0.) THEN
         WRITE(*,*)'  B={0} IN A*X=B'
         STOP
      END IF 

      CALL CRCONV (NEQ,NB,NZU,NZL,IU,IL,AU,AL,IV,LA,LB)
      CALL ICDCMP (NEQ,NB,NZL,IL,AL,AD,UL,UD,LA,LB)
      CALL CRCONV (NEQ,NB,NZL,NZU,IL,IU,UL,UU,IV,LA,LB)
      CALL CG     (NEQ,NB,NZU,NZL,IU,IL,AU,AL,AD,B,EPS,X,
     *             UU,UL,UD,R,P,RU,AP,ITERA,LA,LB)

      WRITE(*,*) '  ICCGM EPS          ',EPS
      WRITE(*,*) '  ICCGM CONVERGENCE  ',REAL(ITERA)/NEQ*100,'%'
      
      RETURN 
      END
      
**************************************************************
      SUBROUTINE CRCONV (NEQ,NB,NZ1,NZ2,I1,I2,A1,A2,IV,LA,LB)
**************************************************************
      INTEGER NZ1(LA),NZ2(LA),I1(LB,LA),I2(LB,LA)
      COMPLEX A1(LB,LA),A2(LB,LA)
      INTEGER IV(LA)

      DO 10 I=1,NEQ
         NZ2(I)=0
         DO 10 J=1,NB
            I2(J,I)=0
            A2(J,I)=( 0.0 , 0.0 )
   10 CONTINUE
            
      DO 20 I=1,NEQ
         DO 30 NU=1,NZ1(I)
         IV(I1(NU,I))=NZ2(I1(NU,I))+1
   30 CONTINUE
         DO 40 NU=1,NZ1(I)
            NZ2(I1(NU,I))=IV(I1(NU,I))
            I2(IV(I1(NU,I)),I1(NU,I))=I
            A2(IV(I1(NU,I)),I1(NU,I))=A1(NU,I)
   40 CONTINUE
   20 CONTINUE
      NZ1MAX=0
      NZ2MAX=0
      DO 200 I=1,NEQ
         IF(NZ1(I).GT.NZ1MAX) NZ1MAX=NZ1(I)
         IF(NZ2(I).GT.NZ2MAX) NZ2MAX=NZ2(I)
         IF(NZ1(I).GT.LB.OR.NZ2(I).GT.LB) THEN
             WRITE(*,*)I,NZ1(I),NZ2(I)
             WRITE(*,*)'  ERR IN SUBROUTINE CRCONV'
            STOP
	 END IF   
  200 CONTINUE

      RETURN
      END
      
**************************************************************
      SUBROUTINE ICDCMP (NEQ,NB,NZL,IL,AL,AD,UL,UD,LA,LB)
**************************************************************
      INTEGER NZL(LA),IL(LB,LA)
      COMPLEX AD(LA),AL(LB,LA),UL(LB,LA),UD(LA)
      COMPLEX T
      
      EPS=1.E-06
      DO 10 K=1,NEQ
         DO 20 MU=1,NZL(K)
            I=IL(MU,K)
            MUI=1
            MUK=1
            UL(MU,K)=AL(MU,K)
    1       CONTINUE
            IF(IL(MUI,I).GT.IL(MUK,K)) THEN
               MUK=MUK+1
            ELSE IF(IL(MUI,I).LT.IL(MUK,K)) THEN
               MUI=MUI+1
            ELSE
               UL(MU,K)=UL(MU,K)
     *                 -UD(IL(MUK,K))*UL(MUI,I)*UL(MUK,K)
               MUI=MUI+1
               MUK=MUK+1
            END IF
            
         IF((MUI.LE.NZL(I)) .AND. (MUK.LE.NZL(K))) GOTO 1
   
   20 CONTINUE      
   
         T=AD(K)
         DO 30 MU=1,NZL(K)
            T=T-UD(IL(MU,K))*CMPLX(UL(MU,K))**2
   30    CONTINUE
         IF(CABS(T).LE.EPS) T=EPS
         
         UD(K)=1/T
     
   10 CONTINUE
   
      RETURN
      END
      
**************************************************************
      SUBROUTINE CG (NEQ,NB,NZU,NZL,IU,IL,AU,AL,AD,B,EPS,X,
     *         UU,UL,UD,R,P,RU,AP,ITERA,LA,LB)
**************************************************************
      INTEGER NZU(LA),NZL(LA),IU(LB,LA),IL(LB,LA)
      COMPLEX AU(LB,LA),AD(LA),AL(LB,LA)
      COMPLEX B(LA),X(LA)
      COMPLEX UU(LB,LA),UL(LB,LA),UD(LA)
      COMPLEX R(LA),P(LA),RU(LA),AP(LA)
      COMPLEX RUR0,RUR1,PAP,RES2,B2
      COMPLEX ALPHA,BETA,Z
      
      KEND=NEQ
      EPS2=EPS**2
      Z=( 0. , 0. )

      CALL ICSOLV(NEQ,NB,NZU,NZL,IU,IL,UU,UL,UD, B , X ,LA,LB)            
      CALL PRODMV(NEQ,NB,NZU,NZL,IU,IL,AU,AL,AD,X,AP,LA,LB)
      
      DO 510 I=1,NEQ
         R(I)=B(I)-AP(I)
  510 CONTINUE
  
      CALL ICSOLV(NEQ,NB,NZU,NZL,IU,IL,UU,UL,UD,R,RU,LA,LB)
            
      RUR0=Z
      DO 520 I=1,NEQ
         P(I)=RU(I)
         RUR0=RUR0+R(I)*CMPLX(RU(I))
  520 CONTINUE
  
      B2=Z
      DO 530 I=1,NEQ
         B2=B2+B(I)**2
  530 CONTINUE
  
      DO 10 K=1,KEND
            CALL PRODMV(NEQ,NB,NZU,NZL,IU,IL,AU,AL,AD,P,AP,LA,LB)
         PAP=Z
         DO 570 I=1,NEQ
            PAP=PAP+P(I)*CMPLX(AP(I))
  570    CONTINUE
         
         ALPHA=RUR0/PAP
         RES2=0.0
         
         DO 550 I=1,NEQ
            X(I)=X(I)+ALPHA*P(I)
            R(I)=    R(I)-ALPHA*AP(I)
            RES2=RES2+R(I)**2
  550    CONTINUE
  
         IF(CABS(RES2/B2).LE.EPS2) THEN
	 ITERA=K
	 RETURN
	 END IF

         CALL ICSOLV(NEQ,NB,NZU,NZL,IU,IL,UU,UL,UD,R,RU,LA,LB)
         
         RUR1=Z        
         DO 560 I=1,NEQ
            RUR1=RUR1+R(I)*CMPLX(RU(I))
  560    CONTINUE
  
         BETA=RUR1/RUR0
         RUR0=RUR1
         
         DO 540 I=1,NEQ
            P(I)=RU(I)+BETA*P(I)
  540    CONTINUE
  
   10 CONTINUE
      
      RETURN
      END
      
**************************************************************
      SUBROUTINE PRODMV (NEQ,NB,NZU,NZL,IU,IL,AU,AL,AD,X,AP,LA,LB)
**************************************************************
      INTEGER NZU(LA),NZL(LA),IU(LB,LA),IL(LB,LA)
      COMPLEX AU(LB,LA),AD(LA),AL(LB,LA)
      COMPLEX X(LA),AP(LA)
      
       COMPLEX APV
      
      DO 10 I=1,NEQ
         APV=0.0D0
         DO 20 MU=1,NZL(I)
            APV=APV+AL(MU,I)*CMPLX(X(IL(MU,I)))
   20    CONTINUE
         APV=APV+AD(I)*CMPLX(X(I))
         DO 30 NU=1,NZU(I)
            APV=APV+AU(NU,I)*CMPLX(X(IU(NU,I)))
   30    CONTINUE
         AP(I)=APV
   
   10 CONTINUE
   
      RETURN
      END
      
**************************************************************
      SUBROUTINE ICSOLV(NEQ,NB,NZU,NZL,IU,IL,UU,UL,UD,B,X,LA,LB)
**************************************************************
      INTEGER NZU(LA),NZL(LA),IU(LB,LA),IL(LB,LA)
      COMPLEX B(LA),X(LA)
      COMPLEX UU(LB,LA),UL(LB,LA),UD(LA)
      
       COMPLEX T
      
      DO 10 K=1,NEQ
         T=B(K)
         DO 520 MU=1,NZL(K)
            T=T-UL(MU,K)*CMPLX(X(IL(MU,K)))
  520    CONTINUE
         X(K)=UD(K)*T
   10 CONTINUE
   
      DO 30 K=NEQ,1,-1
         T=( 0.0 , 0.0 )
         DO 570 NU=1,NZU(K)
            T=T+UU(NU,K)*CMPLX(X(IU(NU,K)))
  570    CONTINUE
         X(K)=X(K)-UD(K)*T
   30 CONTINUE
   
      RETURN
      END

**************************************************************
      SUBROUTINE EDDY (X,XP,YP,ZP)
**************************************************************
      INCLUDE'FEM.H'
      COMPLEX X(N6)
      COMPLEX C,DT,XX

      IF(OMG.EQ.0.) RETURN
      DO 90 I=1, NMTL
         JX(I)=(0.0,0.0)
         JY(I)=(0.0,0.0)
         JZ(I)=(0.0,0.0)
   90 CONTINUE

      CALL SHAPEFUN(XP,YP,ZP)
      C=(0.,1.)

      NPOIN2=NPOIN*2
      NPOIN3=NPOIN*3
      DO 10 I=1,NMTL
         M=NNMTL(I)
         CALL CJACOB(M, CJACOBD)
         DT=C*OMG*SIGMA(M)
         DO 10 J=1, 8
            K=NOD(M,J)
            XX = SIGMA(M)*X(NPOIN3+K)
            JX(I)=JX(I) -SHAPE(J)*X(       K)*DT-XX*DERIVC(1,J)
            JY(I)=JY(I) -SHAPE(J)*X(NPOIN +K)*DT-XX*DERIVC(2,J)
            JZ(I)=JZ(I) -SHAPE(J)*X(NPOIN2+K)*DT-XX*DERIVC(3,J)
            XAE(I) = XAE(I) +SHAPE(J)*XYZ(K,1)
            YAE(I) = YAE(I) +SHAPE(J)*XYZ(K,2)
            ZAE(I) = ZAE(I) +SHAPE(J)*XYZ(K,3)
  10  CONTINUE
 
      RETURN
      END

**************************************************************
      SUBROUTINE FLUX (X)
**************************************************************
      INCLUDE'FEM.H'
      COMPLEX X(N6)

      DO 100 I=1, NELEM
         FX(I)=(0.0,0.0)
         FY(I)=(0.0,0.0)
         FZ(I)=(0.0,0.0)
  100 CONTINUE 

      NP2=NPOIN*2
      NP1=NPOIN*1

      CALL SHAPEFUN(XP,YP,ZP)
      DO 10 I=1,NELEM
         CALL CJACOB(I, CJACOBD)
         DO 10 J=1, 8
            K=NOD(I,J)
            FX(I)=FX(I)+DERIVC(2,J)*X(NP2+K)-DERIVC(3,J)*X(NP1+K)
            FY(I)=FY(I)+DERIVC(3,J)*X(   +K)-DERIVC(1,J)*X(NP2+K)
            FZ(I)=FZ(I)+DERIVC(1,J)*X(NP1+K)-DERIVC(2,J)*X(    K)
  10  CONTINUE

      RETURN
      END

**************************************************************
      SUBROUTINE IMPEDANCE (ICAL,V,A0,Z0,
     *                      RL,ZL,MODEL,IMAX,NEQ,TURN)
**************************************************************
      INCLUDE 'FEM.H' 

      PARAMETER ( NRI =5, NTHI = 19,  NZI = 5 )
      PARAMETER ( NFMAX=(NRI)*(NTHI)*(NZI) )
      COMPLEX   V(N6),AV(NFMAX,3)
      COMPLEX   JX0,JY0,JZ0,DT,Y1,Y2,Y3,XX
      REAL      CXYZ(NFMAX,3)

      IF(OMG.EQ.0.) RETURN
      ZMU0 = 4.0D-7*PAI
      SS = 0.25/PAI*ZMU0
      AT = ZL*RL*CUR1(1)
      
      DZ = ZL/(NZI-1)
      DR = RL/(NRI-1)
      DTH= 2.0*PAI/(NTHI-1)
      DTH= 1.0*PAI/(NTHI-1)

      ND=0
      DO 10 I=1, NZI
         ZPP= Z0 + DZ*(I-1)
         DO 10 J=1, NRI
            APP = A0 + DR*(J-1)
            DO 10 K=1, NTHI
            ND=ND+1
            DEL=DTH*(K-1)
            CXYZ(ND,1)=APP*COS(DEL)
            CXYZ(ND,2)=APP*SIN(DEL)
            CXYZ(ND,3)=ZPP         
  10  CONTINUE
      NF=NRI*NZI*NTHI

      DO 20 I=1, NF
         AV(I,1)=CMPLX(0.0D0,0.0D0)
         AV(I,2)=CMPLX(0.0D0,0.0D0)
         AV(I,3)=CMPLX(0.0D0,0.0D0)
  20  CONTINUE
 
      DO 2000 L=1, IMAX
      DO 2000 M=1, IMAX
      DO 2000 N=1, IMAX
        CALL SHAPEFUN(GAUSS(L),GAUSS(M),GAUSS(N))
        DO 1000 I=1,NMTL
          DT=CMPLX( 0.0 , OMG*SIGMA(NNMTL(I)) )
          JX0=CMPLX(0.0D0,0.0D0)
          JY0=CMPLX(0.0D0,0.0D0)
          JZ0=CMPLX(0.0D0,0.0D0)
          XA=0.0D0
          YA=0.0D0
          ZA=0.0D0
          CALL CJACOB(NNMTL(I) , CJACOBD)
          G=GCOEF(L)*GCOEF(M)*GCOEF(N)*CJACOBD*SS

          DO 30 J=1, 8
             K=NOD(NNMTL(I),J)
             XX = SIGMA(NNMTL(I))*V(NPOIN*3+K)
             Y1 =-SHAPE(J)*V(        K)*DT
             Y2 =-SHAPE(J)*V(NPOIN  +K)*DT
             Y3 =-SHAPE(J)*V(NPOIN*2+K)*DT
             JX0 = JX0 -XX*DERIVC(1,J)+Y1
             JY0 = JY0 -XX*DERIVC(2,J)+Y2
             JZ0 = JZ0 -XX*DERIVC(3,J)+Y3
             XA = XA +SHAPE(J)*XYZ(K,1)
             YA = YA +SHAPE(J)*XYZ(K,2)
             ZA = ZA +SHAPE(J)*XYZ(K,3)
 30       CONTINUE
 
          DO 40 J=1, NF
             X = (CXYZ(J,1)-XA)**2
             Y = (CXYZ(J,2)-YA)**2
             Z = (CXYZ(J,3)-ZA)**2
             R=G/SQRT(X+Y+Z)
             AV(J,1)=AV(J,1)+JX0*R
             AV(J,2)=AV(J,2)+JY0*R
             AV(J,3)=AV(J,3)+JZ0*R

             X = (CXYZ(J,1)-XA)**2
             Y = (CXYZ(J,2)+YA)**2
             Z = (CXYZ(J,3)-ZA)**2
             R=G/SQRT(X+Y+Z)
             AV(J,1)=AV(J,1)-JX0*R
             AV(J,2)=AV(J,2)+JY0*R
             AV(J,3)=AV(J,3)-JZ0*R
  40      CONTINUE

       IF(MODEL.NE.2) GOTO 200
          DO 45 J=1, NF
             X = (CXYZ(J,1)+XA)**2
             Y = (CXYZ(J,2)-YA)**2
             Z = (CXYZ(J,3)-ZA)**2
             R=G/SQRT(X+Y+Z)
             AV(J,1)=AV(J,1)+JX0*R
             AV(J,2)=AV(J,2)-JY0*R
             AV(J,3)=AV(J,3)-JZ0*R

             X = (CXYZ(J,1)+XA)**2
             Y = (CXYZ(J,2)+YA)**2
             Z = (CXYZ(J,3)-ZA)**2
             R=G/SQRT(X+Y+Z)
             AV(J,1)=AV(J,1)-JX0*R
             AV(J,2)=AV(J,2)-JY0*R
             AV(J,3)=AV(J,3)-JZ0*R
  45      CONTINUE
  200   CONTINUE

 1000   CONTINUE
 2000  CONTINUE

      ZIMP(ICAL)=CMPLX(0.0D0,0.0D0)
 
      ND=0
      DO 50 I=1, NZI
      DO 50 J=1, NRI
        A = A0 + DR*(J-1)
        DO 50 K=1, NTHI
          ND=ND+1
          DEL=DTH*(K-1)
          IF(K.EQ.NTHI) GOTO 55
          ZIMP(ICAL)=ZIMP(ICAL)
     *              +2*(-SIN(DEL)*AV(ND,1)+COS(DEL)*AV(ND,2))*DTH*A
 55   CONTINUE
 50   CONTINUE
      ZIMP(ICAL)=ZIMP(ICAL)*(-1.0)*CMPLX(0.0,OMG)*TURN**2/(AT*NRI*NZI)

      RETURN
      END

**************************************************************
      SUBROUTINE  OUTPUTI (NCAL)
**************************************************************
      INCLUDE'FEM.H'
      COMPLEX Z

         WRITE(13,*)' *****  IMPEDANCE  *****'
         WRITE(13,*)'    REAL  ,  IMAGINARY'
         DO 20 I=1,NCAL
   20       WRITE(13,*) ZIMP(I)

      IF(NCAL.GT.1) THEN
         WRITE(13,*)' *****  IMPEDANCE CHANGE  *****'
         WRITE(13,*)'    REAL  ,  IMAGINARY  ,  ABSOLUTE'
         DO 10 I=1,NCAL
         Z=ZIMP(I)-ZIMP(NCAL)
   10       WRITE(13,*) Z,CABS(Z) 
      ELSE
      END IF

      RETURN
      END

**************************************************************
      SUBROUTINE OUTPUT (NEQ,V)
**************************************************************
      INCLUDE'FEM.H'
      COMPLEX V(N6)
 
      WRITE(10,600) NMTL
      DO 10 I=1,NMTL
         WRITE(10,610)XAE(I),YAE(I),ZAE(I)
   10 CONTINUE
      DO 15 I=1,NMTL
         WRITE(10,610)JX(I),JY(I),JZ(I)
   15 CONTINUE 
 
      WRITE(11,600) NELEM
      DO 20 I=1,NELEM
         XAF=0.
         YAF=0.
         ZAF=0.
         DO 21 J=1,8
         K=NOD(I,J)
            XAF=XAF+XYZ(K,1)
            YAF=YAF+XYZ(K,2)
            ZAF=ZAF+XYZ(K,3)
   21    CONTINUE
         XAF=XAF/8
         YAF=YAF/8
         ZAF=ZAF/8
         WRITE(11,610)XAF,YAF,ZAF
   20 CONTINUE
      DO 25 I=1,NELEM
         WRITE(11,610)FX(I),FY(I),FZ(I)
   25 CONTINUE 

      WRITE(12,600) NEQ
      DO 100 I=1,NEQ,6
         WRITE(12,610) (V(I+J),J=0,5)
  100 CONTINUE
  600 FORMAT(12I5)
  610 FORMAT(6E12.5)
 
      RETURN
      END

