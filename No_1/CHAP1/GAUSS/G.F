      PARAMETER( NNN=1000 )
      IMPLICIT DOUBLE PRECISION (A-H,O-W)
      DIMENSION C(NNN+1,NNN)
*
      WRITE(*,*)
      WRITE(*,*) '### Gauss Elimination Method (fortran) ###'
      WRITE(*,*)
      WRITE(*,*) 'N? N<=',NNN
      READ(*,*) NP
*
      DO 10 I=2,NP-1
         DO 11 J=1,NP+1
            C(J,I)=1E-10
   11    CONTINUE
         C(I-1,I)=-1.
         C(I+1,I)=-1.
         C(I,I)=2.
   10 CONTINUE
      DO 20 J=1,NP
         C(J,1)=0.
   20    C(J,NP)=0.
      C(1,1)=1.
      C(NP,NP)=1.
      C(NP+1,1)=1.
      C(NP+1,NP)=FLOAT(NP)
* 
      WRITE(*,*) 'Start ( N=',NP,')  - Input any Number & Cr - '
      READ(*,*) NDUMMY
*
      CALL GAUSS(NP,C,NNN)
*
      N2=NP/2
      WRITE(*,*) N2,C(NP+1,N2)
      END

      SUBROUTINE GAUSS(NP,CC,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-W)
      DIMENSION CC(N+1,N)
      NPP=NP+1
      DO 1000 K=1,NP-1
         DO 10 I=K+1,NP
            CI=-CC(K,I)/CC(K,K)
            DO 50 J=K+1,NPP

               CC(J,I)=CC(J,I)+CC(J,K)*CI
   50       CONTINUE
   10    CONTINUE
 1000 CONTINUE
      CC(NPP,NP)=CC(NPP,NP)/CC(NP,NP)
      DO 2000 I=NP-1,1,-1
         DO 30 J=I+1,NP
            CC(NPP,I)=CC(NPP,I)-CC(NPP,J)*CC(J,I)
   30    CONTINUE
         CC(NPP,I)=CC(NPP,I)/CC(I,I)
 2000 CONTINUE
      RETURN
      END

