C    TRANSIENT EDDY CURRENT ANALYSIS ROUTINE
      DIMENSION  D(100000)
      OPEN(9,FILE='F9')
      OPEN(8,FILE='F8')
      OPEN(12,FILE='F12')
      OPEN(30,FILE='F30',FORM='UNFORMATTED')
      OPEN(31,FILE='f31',FORM='UNFORMATTED') 
C
      NDIM = 100000
C     READ *, NP, NE, NV, NMAT
      NE   =  400
      NP   =  21 * 11
      NMAT =  1
      NV   =  3      
      WRITE(9,303)  NP, NE
 303  FORMAT(2I5)
      M1   =   1  + NP
      M2   =  M1  + NP
      M3   =  M2  + NMAT
      M4   =  M3  + NMAT
      M5   =  M4  + NP
      M6   =  M5  + NP
      M7   =  M6  + NP
      M8   =  M7  + NP
      CALL  PREPAR  ( D, D(M1), D(M2), D(M3), NP, NE, NV, NMAT, NBAND )
      CALL  PREPRO  ( D, D(M1), D(M2), D(M3), D(M4), D(M5),
     *                D(M6), D(M7), D(M8), NP, NE, NV, NMAT, NBAND )
      M9   =  M8  + NBAND*NP
      CALL  CONTRL  ( D, D(M1), D(M2), D(M3), D(M4), D(M5),
     *                D(M6), D(M7), D(M8), D(M9), NP, NE,NMAT, NBAND )
      STOP
      END
      SUBROUTINE  PREPAR ( X, Y, AMU, SIG, NP, NE, NV, NMAT, NBAND )
      DIMENSION  X(NP), Y(NP), AMU(NMAT), SIG(NMAT)
      DIMENSION  KT(10)
      DATA  KT /10*1/
      DATA  NKT, IJ, IV / 1, 0, 0 /
C
C     DO 10 I=1,NP
C10      READ *, X(I), Y(I)
      K      = 0
      DO 10 I=1,21
      DO 10 J=1,11
      K      = K + 1
      X(K)   = 0.001*(I-1)
 10   Y(K)   = 0.001*(J-1)
      DO 15 I=1,NP
 15      WRITE(9,301)  X(I), Y(I)
 301  FORMAT( 2E15.8)
C      DO 20 I=1,NMAT
C 20      READ *, AMU(I), SIG(I)
      AMU(1)  = 3.771E-4
      SIG(1)  = 8.33E5
      WRITE(12,303)  NKT, IJ, IV
 303  FORMAT(3I5)
      WRITE(12,305)  AMU(1), SIG(1)
      WRITE(12,305)  1.
 305  FORMAT(2E15.8) 
      DO 25 I=1,20
      DO 25 J=1,10
      K1     = 1 + (J-1) + 11*(I-1)
      K2     = 2 + (J-1) + 11*(I-1)
      K3     =13 + (J-1) + 11*(I-1)
      K4     =12 + (J-1) + 11*(I-1) 
      KMAT   = 1
      WRITE(31)  K1, K4, K3, KMAT
 25   WRITE(31)  K1, K2, K3, KMAT
      REWIND  31
       NBAND  = 0
      DO 30 I=1,NE
C      READ *, K1, K2, K3, KMAT
      READ(31) K1, K2, K3, KMAT
      NBAND  = MAX( NBAND, ABS(K1-K2) )
      NBAND  = MAX( NBAND, ABS(K1-K3) )
      NBAND  = MAX( NBAND, ABS(K2-K3) )
      A1     = X(K2)*Y(K3) - X(K3)*Y(K2)
      A2     = X(K3)*Y(K1) - X(K1)*Y(K3)
      A3     = X(K1)*Y(K2) - X(K2)*Y(K1)
      DET    =  A1 + A2 + A3
      IF( DET .GT. 0.0 )  THEN
      WRITE(30)  K1, K2, K3, KMAT 
      WRITE(9,308)  K1, K2, K3
 308  FORMAT(3I5)
                          ELSE
      WRITE(30)  K1, K3, K2, KMAT
      WRITE(9,308)  K1, K3, K2
                          ENDIF
 30   CONTINUE
      N  = NE / 10
      DO 31 I=1,N
 31      WRITE(9,300)  KT
 300     FORMAT(10I5)
      IF( NE .GT. 10*N )  WRITE(9,300)  KT
      NBAND  = NBAND + 1                    
      RETURN
      END
      SUBROUTINE  PREPRO  ( X, Y, AMU, SIG, AK, A0, C, AJ, S,
     *                      NP, NE, NV, NMAT, NBAND          )
      DIMENSION  X(NP), Y(NP), AMU(NMAT), SIG(NMAT)
      DIMENSION  AK(NP), A0(NP), C(NP), AJ(NP), S(NBAND,NP)
      DIMENSION  KBOU(31)
      COMMON     DT, T, TMAX, GPHI
      DATA       KBOU  / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 
     *     11, 12, 23, 34, 45,
     *     56, 67, 78, 89, 100, 111, 122, 133, 144, 155, 166, 177, 
     *    188, 199, 210, 221 /
C
C *** MATRIX & VECTOR GENERATE       
      REWIND  30
      DO 10 I=1,NP
      AJ(I)  = 0.0
      C(I)   = 0.0
      DO 10 J=1,NBAND
 10      S(J,I) = 0.0
C      READ *, DT, T, TMAX, GPHI
      DT     = 0.07 / 1000
      T      = 0.7  / 1000
      T      = 100 * T  
      TMAX   = T
      GPHI   = 1.0 
      DO 20 I=1,NE
      READ(30)  K1, K2, K3, KMAT
      ANU    = 1.0 / AMU(KMAT)
      A1     = X(K2)*Y(K3) - X(K3)*Y(K2)
      A2     = X(K3)*Y(K1) - X(K1)*Y(K3)
      A3     = X(K1)*Y(K2) - X(K2)*Y(K1)
      DET    =  A1 + A2 + A3
      B1     = Y(K2) - Y(K3)
      B2     = Y(K3) - Y(K1)
      B3     = Y(K1) - Y(K2)
      C1     = X(K3) - X(K2)
      C2     = X(K1) - X(K3)
      C3     = X(K2) - X(K1)
      W      = (ANU/(2*DET))
      S(1,K1)= S(1,K1) + W*(B1**2+C1**2)
      S(1,K2)= S(1,K2) + W*(B2**2+C2**2)
      S(1,K3)= S(1,K3) + W*(B3**2+C3**2)
      J      = MIN(K1,K2)
      K      = ABS(K1-K2) + 1
      S(K,J) = S(K,J)  + W*(B1*B2+C1*C2)
      J      = MIN(K1,K3)
      K      = ABS(K1-K3) + 1
      S(K,J) = S(K,J)  + W*(B1*B3+C1*C3)
      J      = MIN(K2,K3)
      K      = ABS(K2-K3) + 1      
      S(K,J) = S(K,J)  + W*(B2*B3+C2*C3)
      C(K1)  = C(K1)   + (1.0/DT)*(DET/6.0)*SIG(KMAT)
      C(K2)  = C(K2)   + (1.0/DT)*(DET/6.0)*SIG(KMAT)
      C(K3)  = C(K3)   + (1.0/DT)*(DET/6.0)*SIG(KMAT)
      AJ(K1) = AJ(K1)  +          (DET/6.0)*SIG(KMAT)*GPHI
      AJ(K2) = AJ(K2)  +          (DET/6.0)*SIG(KMAT)*GPHI
      AJ(K3) = AJ(K3)  +          (DET/6.0)*SIG(KMAT)*GPHI
 20   CONTINUE
C *** Q-MATRIX GENERATE
      DO 30 I=1,NP
 30      S(1,I) = S(1,I) + C(I)
C *** BOUNDARY CONDITION TREATMENT
C      READ *, NBOU
       NBOU   = 31
      DO 60 I=1,NBOU
C      READ *, J, AS
      J       = KBOU(I) 
      AS      = 0.0
      DO 40 K=2,NBAND
      L      = J + K - 1
      IF( L .LE. NP )  THEN
      AJ(L)  = AJ(L) - S(K,J)*AS
      S(K,J) = 0.0
                       ENDIF
 40                    CONTINUE
      DO 50 K=2,NBAND
      L      = J - K + 1
      IF( L .GE. 1 )   THEN
      AJ(L)  = AJ(L) - S(K,L)*AS
      S(K,L) =0.0
                       ENDIF
 50                    CONTINUE
      S(1,J) = 1.0
      C(J)   = 0.0
      AJ(J)  = AS
 60   CONTINUE
C *** INITIAL CONDITION SET
      DO 70 I=1,NP
 70      A0(I) = 0.0
C *** TRIANGULAR DECOMPOSITION
      DO 90 N=1,NP
      DO 90 L=2,NBAND
      IF( S(L,N) .EQ. 0.0 )  GO TO 90
      D      = S(L,N) / S(1,N)
      I      = N + L - 1
      J      = 0 
      DO 80 K=L,NBAND
      J      = J + 1
 80   S(J,I) = S(J,I) - D*S(K,N)
      S(L,N) = D
 90   CONTINUE
      RETURN
      END
      SUBROUTINE  CONTRL  ( X, Y, AMU, SIG, AK, A0, C, AJ, Q, AIE, 
     *                      NP, NE,     NMAT, NBAND          )
      DIMENSION  X(NP), Y(NP), AMU(NMAT), SIG(NMAT)
      DIMENSION  AK(NP), A0(NP), C(NP), AJ(NP), Q(NBAND,NP)
      DIMENSION  AIE(NP)
      COMMON     DT, T, TMAX, GPHI
C

C *** WRITE PARAMETER SET
      NT     = T / DT    + 0.5
      NLOOP  = TMAX / T  + 0.5  
C *** BEGIN STEP-BY-STEP METHOD
      DO 100 M=1,NLOOP
      DO  60 J=1,NT
      DO  10 I=1,NP
 10      AK(I) = C(I)*A0(I) + AJ(I)
C *** FORWARD/BACKWARD SUBSTITUTION 
      DO 20 N=1,NP
      DO 20 L=2,NBAND
      I      = N + L - 1
 20   AK(I)  = AK(I) - Q(L,N)*AK(N)
      DO 30 N=1,NP
 30      AK(N)  = AK(N)/Q(1,N)
      DO 40 N=NP,1,-1
      DO 40 K=2,NBAND
 40   AK(N)   = AK(N) - Q(K,N)*AK(N+K-1)
      IF( J .EQ. NT )  GO TO 60
      DO 50 N=1,NP
 50      A0(N)  = AK(N)
 60   CONTINUE
C *** EDDY CURRENT OUT PUT
      DO 70 I=1,NP
c       AIE(I)   = AJ(I) - (AK(I)-A0(I))*C(I)
 70    A0(I) = AK(I)
      WRITE(8,555)   AK, AK
 555  FORMAT(5E15.8)
 100   CONTINUE
       RETURN
       END

