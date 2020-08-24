C###########################################################################
C            NON-LINEAR STATIC MAGNETIC 2D FIELD
C###########################################################################
      PROGRAM MAIN
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      PARAMETER (NP1=6000,NP2=200000,NE1=10000,NP4=20)
      COMMON/FINITE/NP,NE,ISEN
      COMMON/PARAM/IBM1,IBM3,IBM5
      COMMON/MATTER/ NMAT,NCO,NPB
      COMMON/MIMAC1/NN,NVP,MAXCOL,IP,JPN
      DIMENSION X(NP1),Y(NP1)
      DIMENSION IPIV(NP1),IRST(NP1)
      DIMENSION IREN(NP1),SDIA(NP1)
      DIMENSION ICOL(NP2),SMAT(NP2)
      DIMENSION NSP(3,NE1),DELTA(NE1)
      DIMENSION IQ(NP1),IR(NP1),IL(NP1),IU(NP1)
      DIMENSION WORK(NP1,4)
      DIMENSION POTENTIAL(NP1),DELTAPOTEN(NP1),CURENTVECTOR(NP1),
     &          CDEN(NP1)
      DIMENSION NBR(NP1)
      DIMENSION NIRON(NP4),IRON(NE1),NOO(NP4),NO(NE1)
      DIMENSION AMPRITUDE(NP4),AMEW(NP4)
      CHARACTER*33 DOCU03,DOCU07,DOCU17

C***  OPEN STATEMENT FOR UNIX >>>>>>>>>>>>>>>>>>>>>>>>>>>
      DOCU03='INPUT FEM CONECTION INF.          '
      DOCU07='INPUT FEM TRIANGLE PART DATA INF. '
      DOCU17='OUTPUT ALL AREA POTENTIAL INF.    '
      CALL ASGUNX(03,1,'in ',DOCU03,'bcd')
      CALL ASGUNX(07,1,'in ',DOCU07,'dat')
      CALL ASGUNX(17,2,'out',DOCU17,'pot')

      IBM1=NP1
      IBM3=NP2
      IBM5=NP2
      ISEN=3

C **  SUB1: FINITE ELEMENT DATA INPUT....
      CALL TRIAIN( X,Y,NSP,NP1,NE1,NP4,NIRON,IRON,NOO,NO )
      NITR=NP*3
      EITR=1.E-6

C **  BOUNDS : BOUNDRY CONDITION.... 
      CALL BOUNDS( X,Y,NSP,NBR,NP1,NE1,ITERA,AMEW,AMPRITUDE,NP4 )

c--------------------------------------------------
      DO 1000 INUM=1,ITERA
         DO 10 I=1,NP
           IPIV(I)=0
           CURENTVECTOR(I)=0.D0
   10    CONTINUE

C **  MATRIX : CONSTRUCTION OF THE GLOBAL COEFFICIENT MATRIX
      CALL MATRIX( X,Y,NSP,DELTA,IRST,IREN,ICOL,SMAT,SDIA,
     &           NBR,NIRON,IRON,NOO,NO,AMEW,IPIV,
     &           POTENTIAL,CURENTVECTOR,NP1,NP2,NP4,NE1,INUM )

       WRITE(*,'(1H ,''NO OF NONZERO ELEMENT IN CO-EFFICIENT MATRIX='',
     &                I6)' ) MAXCOL       

C **  VECTOR : CONSTRUCT THE INPUT CURRENT DENSITY VECTOR.
      CALL VECTOR( NSP,DELTA,CURENTVECTOR,CDEN,INUM,X,Y,NP1,NP4,NE1,
     &             NBR,NIRON,IRON,NOO,NO,AMPRITUDE )


      WRITE(*,*) ' *** ICCG   START *** '
C **  ICCG : 
      CALL ICCG( SDIA,SMAT,CURENTVECTOR,DELTAPOTEN,WORK,IRST,IREN,ICOL,
     &           IQ,IR,IL,IU,NITR,EITR,np,NP1,NP2)

       DO 20  I=1,NP
 	 POTENTIAL(I)=POTENTIAL(I)+DELTAPOTEN(I)
  20   CONTINUE

          DO 30 I=1,NP
            IF( DABS(POTENTIAL(I)).GT.1.E-20 .AND.
     &        DABS( DELTAPOTEN(I)/POTENTIAL(I) ).GT.4.E-2 ) GOTO 200
 30        CONTINUE
           GO TO 300

  200 WRITE(*,'(1h  ''NEWTON-ERROR='',E15.5)') DABS( DELTAPOTEN(I)/
     &                                             POTENTIAL(I))
      WRITE(*,'(1h  ''SETTEN BANGOU='',I7)') I
      WRITE(*,'(1h  ''N.R.M. NO KAISUU='',I3)') INUM

 1000 CONTINUE

C **  OUTPUT : OUTPUT POTENTIAL DATA
  300   CALL OUTPUT( POTENTIAL,NP1 )
        STOP 'NORMAL'
        END

C********************************************************************
C  OUTPUT : OUTPUT POTENTIAL
C*******************************************************************
       SUBROUTINE OUTPUT(POTENTIAL,NP1)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       COMMON/FINITE/NP,NE,ISEN
       DIMENSION POTENTIAL(NP1)
       DO 100 I=1,NP
	  WRITE(17,*) POTENTIAL(I)
 100   CONTINUE
       RETURN
       END 
C*****************************************************************
C     INPUT : FINITE ELEMENT DIVISION DATA INPUT
C*****************************************************************
      SUBROUTINE TRIAIN (X,Y,NSP,NP1,NE1,NP4,NIRON,IRON,NOO,NO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      COMMON/FINITE/NP,NE,ISEN
      COMMON/MATTER/ NMAT,NCO,NPB
      DIMENSION X(NP1),Y(NP1),NSP(3,NE1)
      DIMENSION NIRON(NP4),IRON(NE1),NOO(NP4),NO(NE1)

C ****   READING NUMBER OF POINTS AND ELEMNTS
         READ(7,*) NP,NE

         IF(NP.GT.NP1 .OR. NE.GT.NE1) THEN
            STOP 'DIMENSION OVER: CHANGE NP1 OR NE1'
         ENDIF

C ****   READING COORDINATES OF POINTS
         DO 100 I=1,NP
           READ(7,*) X(I),Y(I)
 100     CONTINUE

C ****   READING ELEMENT NUMBER OF EACH TRIANGULAR ELEMENT
         DO 200 I=1,NE
          READ(7,*) (NSP(J,I),J=1,3)
 200     CONTINUE

C ****   READING  NUMBER OF MATERIALS
         READ(7,'(I3)') NMAT
         IF(NMAT.GT.NP4) THEN
           STOP 'DIMENSION OVER: CHANGE NP4'
         ENDIF
           READ(7,'(20I4)') (NIRON(I),I=1,NMAT)

        NIRN=0
        DO 10 I=1,NMAT
          NIRN=NIRN+NIRON(I)
   10   CONTINUE
          READ(7,'(20I4)') (IRON(I),I=1,NIRN)

C ****  READING  NUMBER OF COILS      
        READ(7,'(I3)') NCO
        READ(7,'(20I4)')(NOO(I),I=1,NCO)
        NCOIL=0
        DO 20 I=1,NCO
         NCOIL=NCOIL+NOO(I)
   20   CONTINUE
        READ(7,'(20I4)')(NO(I),I=1,NCOIL)

      DO 30 I=1,NP
        X(I)=X(I)*1.D-3
        Y(I)=Y(I)*1.D-3
   30 CONTINUE
      RETURN
      END
C     *************************************
C     BOUNDS : INPUT OF BOUNDARY CONDITION
C     *************************************

      SUBROUTINE BOUNDS(X,Y,NSP,NBR,NP1,NE1,ITERA,AMEW,AMPRITUDE,NP4)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      COMMON/FINITE/NP,NE,ISEN
      COMMON/MATTER/NMAT,NCO,NPB 
      DIMENSION NBR(NP1)
      DIMENSION X(NP1),Y(NP1),NSP(3,NE1)
      DIMENSION NS(30),NL(30),IBOD(30)
      DIMENSION AMEW(NP4),AMPRITUDE(NP4)

      CHARACTER*2 TYPE(2),ATY
      DATA TYPE/'CL','SH'/
      HANTEI(XS,YS,XE,YE,XO,YO)=(YE-YS)*XO+(XS-XE)*YO-XS*YE+YS*XE
        NC=0
        NPB=0
        READ(3,'(6X,I2,7X,I5)') NLINE,ITERA
        IF(ITERA.EQ.0) ITERA=2*NP

      DO 100 I=1,NLINE
            READ(3,'(A2,2I3)') ATY,NS(I),NL(I)

            DO 101 J=1,2
              IF(ATY.EQ.TYPE(J)) GO TO 102
  101       CONTINUE
            STOP 'UNMATCH BOUNDRY CONDITION'
  102     CONTINUE

                    IBOD(I)=J
                    IF(J.EQ.1)  THEN 
                            IF(NPB.EQ.0) THEN
                                NPB=NPB+1
                                NBR(NPB)=NS(I)

                                IF(NS(I).EQ.NL(I)) GO TO 100
                                NPB=NPB+1
                                NBR(NPB)=NL(I)

                                GO TO 100
                            ENDIF

                            DO 106 J1=1,NPB
                                IF(NBR(J1).EQ.NS(I)) GO TO 107
  106                       CONTINUE
                             NPB=NPB+1
                             NBR(NPB)=NS(I)
 
  107                        CONTINUE
                            DO 108 J1=1,NPB
                               IF(NBR(J1).EQ.NL(I)) GO TO 100
  108                       CONTINUE
                             NPB=NPB+1
                             NBR(NPB)=NL(I)
                    ENDIF

  100 CONTINUE

      KKP=0
      DO 200 I=1,NLINE
              IF(IBOD(I).NE.2) THEN 
                   I1=NS(I)
                   I2=NL(I)
                    IF(I1.EQ.I2) GO TO 200
                   XS=X(I1)
                   YS=Y(I1)
                   XE=X(I2)
                   YE=Y(I2)
                    IF(XS.GT.XE) THEN
                       XL=XE-1.E-3
                       XG=XS+1.E-3
                    ELSE
                       XL=XS-1.E-3
                       XG=XE+1.E-3
                    ENDIF

                    IF(YS.GT.YE) THEN
	               YL=YE-1.E-3
	               YG=YS+1.E-3
                    ELSE
	               YL=YS-1.E-3
	               YG=YE+1.E-3
                    ENDIF
                 XY2=SQRT((XS-XE)*(XS-XE)+(YS-YE)*(YS-YE))

                 DO 300 J=1,NP
                    N1=J
                    XO=X(N1)
                    YO=Y(N1)
                     IF(XO.LT.XL.OR.XO.GT.XG) GOTO 300
                     IF(YO.LT.YL.OR.YO.GT.YG) GOTO 300
                       HANT=HANTEI(XS,YS,XE,YE,XO,YO)
                       HANT=HANT/XY2
                        IF(DABS(HANT).GT.1.E-3) GO TO 300
                         IF(IBOD(I).EQ.1) THEN  
                              DO 50 K=1,NPB
                                 IF(NBR(K).EQ.N1) GO TO 300
   50                         CONTINUE
                               NPB=NPB+1
                               NBR(NPB)=N1
                               GO TO 300
                         ENDIF
  300            CONTINUE
             ENDIF
  200  CONTINUE

           IF(NMAT.NE.0) THEN
                 READ(3,'(5X,4E15.7)') (AMEW(I),I=1,NMAT)
           ENDIF
           DO 110 I=1,NCO
              READ(3,'(11X,F10.3)') AMPRITUDE(I)
  110      CONTINUE
       RETURN
       END

C     *****************************************************************
C      MATRIX : CONSTRUCTION OF THE GLOBAL MATRIX
C     ****************************************************************
      SUBROUTINE MATRIX(X,Y,NSP,DELTA,IRST,IREN,ICOL,SMAT,SDIA,
     &                  NBR,NIRON,IRON,NOO,NO,AMEW,
     &                  IPIV,POTENTIAL,CURENTVECTOR,NP1,NP2,NP4,
     &                  NE1,INUM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      COMMON/FINITE/NP,NE,ISEN
      COMMON/MATTER/NMAT,NCO,NPB
      COMMON/MIMAC1/NN,NVP,MAXCOL,IP,JPN
      COMMON/SAI1/IPIVN0,ICOL0,MAX,IFROW,LAROW,JCOL0,IWRT
      COMMON/LIMIT/LMAT,LCOL,LERS,LLCOL

      DIMENSION NIRON(NP4),IRON(NE1),NOO(NP4),NO(NE1),AMEW(NP4)  
      DIMENSION X(NP1),Y(NP1),NSP(3,NE1),DELTA(NE1),
     &       IRST(NP1),IREN(NP1),ICOL(NP2),SMAT(NP2),SDIA(NP1),
     &       IPIV(NP1),POTENTIAL(NP1),NBR(NP1)
      DIMENSION YMAT(3,3),CURENTVECTOR(NP1)

      IF(INUM.EQ.1) THEN
           DO 10 I=1,NP
               POTENTIAL(I)=0.d0
   10      CONTINUE
      ENDIF
      LCOL=NP2
      LMAT=NP2
      LERS=NP1
      MAXCOL=0
      LLCOL=LCOL
      LAROW=NP
      IFROW=1
      ICOL0=1
      IPIVN0=1
      N=0
      NN=NP
C **  POINTER SET OF MATRIX MAKING 
      NUM=LCOL/NP
      DO 30 I=1,NP
         IRST(I)=N+1
         IREN(I)=N
         N=N+NUM
   30 CONTINUE
      IRST(NP+1)=LCOL+1

      DO 40 I=1,NP
         SDIA(I)=0.d0
   40 CONTINUE

      DO 50 I=1,NP2
         SMAT(I)=0.d0
   50 CONTINUE

      NVP=0
      IFLUX=0

      DO 1000 I=1,NE
                  JJJ=0
                  IF(NMAT.EQ.0) GO TO 500
                  NIRN2=0
                  DO 510 J=1,NMAT
                     NIRN1=NIRN2+1
                     NIRN2=NIRN2+NIRON(J)
                    JJJ=JJJ+1
                   DO 510 K=NIRN1,NIRN2
                      IF(I.EQ.IRON(K)) THEN
                        IF(AMEW(J).LT.0) THEN
                           ISIGN=-JJJ
                        ELSE
                           ISIGN=JJJ
                        ENDIF
                           GOTO 500
                      ENDIF
  510             CONTINUE
                  ISIGN=0

  500             CONTINUE

C **  SUBMAT : CONSTRUCT ELEMENT MATRICES FOR EACH TRIANGULAR ELEMENTS **
                    CALL SUBMAT(ISIGN,X,Y,NSP,I,NP1,NP4,NE1,DELTA,
     &                       POTENTIAL,INUM,CURENTVECTOR,YMAT,AMEW)

                  CALL SUBSTI (NP1,NP2,NE1,1,ISEN,IREN,IRST,
     &                 SMAT,SDIA,NSP,ICOL,IPIV,I,YMAT,NBR)

                  IF(NVP.GT.LMAT) STOP 'LMAT DIMENSION OVER'
 1000 CONTINUE

      RETURN
      END
C     ******************************************************************
C      SUBMAT: CONSTRUCTION OF ELEMENT MATRICES FOR 
C               1ST ORDER TRIANGULAR ELEMENT
C     ******************************************************************
      SUBROUTINE SUBMAT( ISIGN,X,Y,NSP,JELEM,NP1,NP4,NE1,DELTA,
     &                    POTENTIAL,INUM,CURENTVECTOR,YMAT,AMEW )
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      COMMON/FINITE/NP,NE,ISEN
      DIMENSION  YMAT(ISEN,ISEN)
      DIMENSION  X(NP1),Y(NP1),NSP(3,NE1),DELTA(NE1)
      DIMENSION  CURENTVECTOR(NP1),AMEW(NP4)
      DIMENSION  UU(3),UUU(3,3)
      DIMENSION  POTENTIAL(NP1),BE1(3,3),CE1(3,3)              
      DIMENSION  Q(3),R(3),POT(3)
      REL=7.95774D5
      IF(ISIGN .GT. 0) THEN
         REL=REL/ ( AMEW(ISIGN) )
      ENDIF

      NSP1=NSP(1,JELEM)
      NSP2=NSP(2,JELEM)
      NSP3=NSP(3,JELEM)

      Q(1)=Y(NSP2)-Y(NSP3)
      Q(2)=Y(NSP3)-Y(NSP1)
      Q(3)=Y(NSP1)-Y(NSP2)
      R(1)=X(NSP3)-X(NSP2)
      R(2)=X(NSP1)-X(NSP3)
      R(3)=X(NSP2)-X(NSP1)
      BE1(1,1)=Q(1)*Q(1)
      BE1(1,2)=Q(1)*Q(2)
      BE1(1,3)=Q(1)*Q(3)
      BE1(2,2)=Q(2)*Q(2)
      BE1(2,3)=Q(2)*Q(3)
      BE1(3,3)=Q(3)*Q(3)
      CE1(1,1)=R(1)*R(1)
      CE1(1,2)=R(1)*R(2)
      CE1(1,3)=R(1)*R(3)
      CE1(2,2)=R(2)*R(2)
      CE1(2,3)=R(2)*R(3)
      CE1(3,3)=R(3)*R(3)
      DO 10 J=2,3
	N=J-1
	DO 10 K=1,N
	  BE1(J,K)=BE1(K,J)
	  CE1(J,K)=CE1(K,J)
   10 CONTINUE

C **  DELTA(JELEM) : THE AREA OF TRIANGLE ELEMENT ( TWO DELTA-S )
      DELTA(JELEM)=X(NSP1)*Q(1)+X(NSP2)*Q(2)+X(NSP3)*Q(3)
      IF(DELTA(JELEM).LE.0.D0) THEN
            STOP 'AREA OF ELEMENT BECOMES NEGATIVE'
      ENDIF
      DLT=DELTA(JELEM)
      POT(1)=POTENTIAL(NSP1)
      POT(2)=POTENTIAL(NSP2)
      POT(3)=POTENTIAL(NSP3)

             DO 40 L=1,3
                UU(L)=0.D0
  	        DO 50 J=1,3
	          UU(L) = UU(L)+((CE1(L,J)+BE1(L,J))*0.5/DLT)*POT(J)
   50           CONTINUE
   40        CONTINUE

        DO 90 I=1,3
  	  DO 90 J=I,3
            UUU(I,J)=UU(I)*UU(J)
   90   CONTINUE

      IF( INUM.NE.1 .AND. ISIGN.LT.0) THEN
           CALL FLUX(REL,DRDB,DLT,JELEM,Q,R,POT)
      ENDIF

      DO 200 I=1,3
         DO 300  J=I,3
                   YMAT(I,J)= REL*(BE1(I,J)+CE1(I,J))*0.5/DLT  +
     &             4.*DRDB*UUU(I,J)/DLT
 300     CONTINUE
 200  CONTINUE

      DO 110 K=2,3
          J=K-1
      DO 110 N=1,J
          YMAT(K,N)=YMAT(N,K)
 110  CONTINUE

       CURENTVECTOR(NSP1)=CURENTVECTOR(NSP1) - REL*UU(1)
       CURENTVECTOR(NSP2)=CURENTVECTOR(NSP2) - REL*UU(2)
       CURENTVECTOR(NSP3)=CURENTVECTOR(NSP3) - REL*UU(3)

       RETURN
       END
c     ***********************************************************
c     FLUX : CALCULATION OF RELUCTIVITY  FOR 1ST ORDER ELEMENT
c     ***********************************************************
      SUBROUTINE FLUX(REL,DRDB,DLT,IU,Q,R,POT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE  
      DIMENSION RE(3),QE(3)
      DIMENSION Q(3),R(3)
      DIMENSION POT(3)
      RRC=0.d0
      QQC=0.d0

      DO 20  I=1,3
        RRC=RRC+R(I)*POT(I)
        QQC=QQC-Q(I)*POT(I)
  20  CONTINUE

       FDX=RRC*RRC/(DLT*DLT)
       FDY=QQC*QQC/(DLT*DLT)
       FD = FDX+FDY

       IF(FD.LE.1.) THEN
          REL=191.3033+8.69673*FD
          DRDB=8.69673
       ELSEIF(FD.LE.2.56) THEN
          REL=199. + FD**8.69673
          DRDB=8.69673*(FD**7.69673)
       ELSE
          REL = 12063.37*FD-27132.22
          DRDB= 12063.37
       ENDIF

      RETURN
      END
C**********************************************************************
C          VECTOR : CONSTRUCTION OF THE RIGHT SIDE VECTOR
C***********************************************************************
      SUBROUTINE VECTOR(NSP,DELTA,CURENTVECTOR,CDEN,
     &       INUM,X,Y,NP1,NP4,NE1,NBR,NIRON,
     &       IRON,NOO,NO,AMPRITUDE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      COMMON/FINITE/NP,NE,ISEN
      COMMON/MATTER/NMAT,NCO,NPB
      DIMENSION NBR(NP1),NIRON(NP4)
      DIMENSION IRON(NE1),NOO(NP4),NO(NE1)
      DIMENSION NSP(3,NE1),DELTA(NE1)
      DIMENSION UU(3),AMPRITUDE(NP4)
      DIMENSION CDEN(NP1),CURENTVECTOR(NP1)
      DIMENSION X(NP1),Y(NP1)
      DIMENSION BE1(3,3),CE1(3,3),Q(3),R(3)
       IF(INUM.EQ.1) THEN

              DO 100 IP=1,NP
	        CDEN(IP)=0.d0
  100        CONTINUE

                   NO2=0
                   DO 20 I=1,NCO
                            NO1=NO2+1
                            NO2=NO2+NOO(I)
                            DEL=0.D0

                        DO 30 J=NO1,NO2
                          L=NO(J)
                          DEL=DEL+DELTA(L)
   30                   CONTINUE
                        DEL=DEL*0.5D0
                        CD=AMPRITUDE(I)/DEL
                        CONTINUE
                           DO 40 J=NO1,NO2
                             L=NO(J)
                             N1=NSP(1,L)
                             N2=NSP(2,L)
                             N3=NSP(3,L)
                             DECD=DELTA(L)*CD/6.D0
                             CDEN(N1)=CDEN(N1)+DECD
                             CDEN(N2)=CDEN(N2)+DECD
                             CDEN(N3)=CDEN(N3)+DECD
   40                      CONTINUE
   20                 CONTINUE
         ENDIF

        DO 270 I=1,NP
           CURENTVECTOR(I)=CURENTVECTOR(I) + CDEN(I)
  270   CONTINUE

      IF(NPB.NE.0) THEN
         DO 280 I=1,NPB
              L=NBR(I)
              CURENTVECTOR(L)=0.d0
  280    CONTINUE
      ENDIF
      RETURN
      END

C     *******************************************************
C     SUBROUTINE SUBSTI : SETTING ZERO BOUNDARY CONDITION
C     *******************************************************
      SUBROUTINE SUBSTI (NP1,NP2,NE1,NSET1,NSET2,IREN,IRST,
     &                  SMAT,SDIA,NSP,ICOL,IPIV,I,YMAT,NBR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      COMMON/FINITE/NP,NE,ISEN
      COMMON/MATTER/NMAT,NCO,NPB
      DIMENSION YMAT(ISEN,ISEN)
      DIMENSION SDIA(NP1),SMAT(NP2),ICOL(NP2),IRST(NP1),IREN(NP1)
      DIMENSION NSP(3,NE1),IPIV(NP1)
      DIMENSION NBR(NP1)
      DO 10 J=NSET1,NSET2
             NSPJ=NSP(J,I)
             IF(NPB.NE.0) THEN
                 DO 30 K=1,NPB
                   IF(NSPJ.EQ.NBR(K)) GO TO 70
   30            CONTINUE
             ENDIF

         DO 20 K=NSET1,NSET2
               NSPK=NSP(K,I)
               AX=YMAT(J,K)
               IF(DABS(AX).LT.1.E-10) GO TO 20

               IF(NSPK.LT.NSPJ) GO TO 20

               IF(NPB.NE.0) THEN
                 DO 90 L=1,NPB
                  IF(NBR(L).EQ.NSPK) GO TO 20
   90            CONTINUE
               ENDIF

                   CALL ARANGE(NSPJ,NSPK,AX,IRST,IREN,SDIA,
     &                        IPIV,ICOL,SMAT,NP1,NP2)
   20    CONTINUE
         GO TO 10
   70    SDIA(NSPJ)=1.d0
   10 CONTINUE
      RETURN
      END

C*********************************************************************C
C     ARRANGE OF GLOBAL MATRIX OF THE SYSTEM                          C
C*********************************************************************C
      SUBROUTINE ARANGE(NSPJ,NSPK,AX,
     &                  IRST,IREN,SDIA,IPIV,ICOL,SMAT,NP1,NP2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /MIMAC1/ NN,NVP,MAXCOL,IP,JPN
      COMMON /LIMIT/  LMAT,LCOL,LERS,LLCOL
      COMMON /SAMAT/  BMAT

      DIMENSION IRST(NP1),IREN(NP1)
      DIMENSION SDIA(NP1),IPIV(NP1)
      DIMENSION ICOL(NP2),SMAT(NP2)
      SAVE

      IF(NSPJ.EQ.NSPK) THEN
        SDIA(NSPJ)=SDIA(NSPJ)+AX
      ELSE
          I1=IRST(NSPJ)
          I2=IREN(NSPJ)
          DO 500 ITN=I1,I2
            IF(ICOL(ITN).EQ.NSPK) THEN
              SMAT(ITN)=SMAT(ITN)+AX
                                                   GOTO 100
            ELSE IF(ICOL(ITN).GT.NSPK) THEN
              ILO=ITN
                                                   GOTO 100
            ENDIF
  500     CONTINUE
          ILO=I2+1
 100     CONTINUE
            I3=IRST(NSPK)
            I4=IREN(NSPK)
         DO 600 ITN=I3,I4
            IF(ICOL(ITN).EQ.NSPJ) THEN
               SMAT(ITN)=SMAT(ITN)+AX
                                                   GOTO 300
            ELSE IF(ICOL(ITN).GT.NSPJ) THEN
               ILP=ITN
                                                   GOTO 200
            ENDIF
  600     CONTINUE
           ILP=I4+1
 200      CONTINUE
           BMAT=AX

           CALL RENEW(NSPJ,NSPK,ILO,ILP,IS,IRST,IREN,ICOL,IPIV,
     &                SMAT)

           MAXCOL=MAXCOL+2
          IF(MAXCOL.GT.LCOL) STOP 'SMAT IS OVER'
      ENDIF
  300 CONTINUE
      RETURN
      END
C****************************************************************
C
C****************************************************************
      SUBROUTINE RENEW(J1,J2,ILO,ILP,IVAR,IRST,IREN,ICOL,IPIV,AMAT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      COMMON/PARAM/IBM1,IBM3,IBM5
      COMMON/MIMAC1/NN,NVP,MAXCOL,IP,JPN
      COMMON/MIMAC5/MARGIN(2,6),LMG
      COMMON /SAI1/ IPIVN0,ICOL0,MAX,IFROW,LAROW,JCOL0,IWRT
      DIMENSION IRST(IBM1),IREN(IBM1),ICOL(IBM5),AMAT(IBM3),IPIV(IBM1)
      DIMENSION NCASE(6),NKASE(2,6)
      DATA NKASE/1,2,2,4,2,5,3,4,3,5,5,6/
      LMG=0
      DO 1 I=1,6
      NCASE(I)=30000000
    1 MARGIN(2,I)=-10000000
      CALL SERCH(J1,J2-1,1,3,IRST,IREN,IPIV)
      CALL SERCH(J2-1,J1,1,4,IRST,IREN,IPIV)
      LIFT=LMG
      LB=2
      IF(LIFT.NE.0) LB=1
      LI=LMG
      IF(J1.EQ.IFROW) GO TO 2
      CALL SERCH(J1-1,IFROW,LB,2,IRST,IREN,IPIV)
    2 LJ=LMG
      LK=LB-LJ+LI
      IF(LK.LE.0) GO TO 3
      L=2
      NP=IRST(IFROW)-1
    4 IF(NP.LE.ICOL0-1) GO TO 3
      MARGIN(1,L)=IFROW-1
      MARGIN(2,L)=NP
      LMG=LMG+1
      L=L-1
      NP=NP-1
      IF(L.LE.0) GO TO 3
      IF(LK.EQ.2) GO TO 4
    3 CALL SERCH(J2,LAROW,LB,5,IRST,IREN,IPIV)
      IF(LMG.LE.1) STOP 10
      NCASE(1)=ILP-MARGIN(2,1)-2
      NCASE(2)=ILO-MARGIN(2,2)+ILP-MARGIN(2,4)-2
      M3=MARGIN(2,3)-ILO
      M5=MARGIN(2,5)-ILP
      M6=MARGIN(2,6)-ILO-1
      IF(M5.GE.0) NCASE(3)=ILO-MARGIN(2,2)-1+M5
      IF(M3.GE.0) NCASE(4)=M3+ILP-MARGIN(2,4)-1
      IF(M3.GE.0.AND.M5.GE.0) NCASE(5)=M3+M5
      IF(M6.GE.0) NCASE(6)=M6
      IF(MARGIN(2,3).EQ.MARGIN(2,4)) NCASE(4)=30000000
C     EXECUTE OF SLIDING
      MIN=30000000
      DO 310 I=1,6
      IF(MIN.LE.NCASE(I)) GO TO 310
      MIN=NCASE(I)
      MINI=I
  310 CONTINUE
  320 CONTINUE
      MGO=NKASE(1,MINI)
      MGP=NKASE(2,MINI)
      MGO1=MARGIN(1,MGO)
      MGO2=MARGIN(2,MGO)
      MGP1=MARGIN(1,MGP)
      MGP2=MARGIN(2,MGP)
      IF(MINI.EQ.6) GO TO 470
      IF(MINI.GT.1) GO TO 480
      CALL SLIDE(MGO1,MGP1,MGO2,MGP2,1,0,IVAR,IRST,IREN,ICOL,IPIV,AMAT)
      CALL SLIDE(MGP1,J1,MGP2,ILO,2,J2,IVAR,IRST,IREN,ICOL,IPIV,AMAT)
      ILO=ILO-1
      CALL SLIDE(J1,J2,ILO,ILP,1,J1,IVAR,IRST,IREN,ICOL,IPIV,AMAT)
      GO TO 490
  470 CONTINUE
      MP=MGO2+1
      CALL SLIDE(MGP1,MGO1,MGP2,MP,1,0,IVAR,IRST,IREN,ICOL,IPIV,AMAT)
      CALL SLIDE(MGO1,J2,MGO2,ILP,2,J1,IVAR,IRST,IREN,ICOL,IPIV,AMAT)
      CALL SLIDE(J2,J1,ILP,ILO,1,J2,IVAR,IRST,IREN,ICOL,IPIV,AMAT)
      GO TO 490
  480 CONTINUE
      CALL SLIDE(MGO1,J1,MGO2,ILO,1,J2,IVAR,IRST,IREN,ICOL,IPIV,AMAT)
      CALL SLIDE(MGP1,J2,MGP2,ILP,1,J1,IVAR,IRST,IREN,ICOL,IPIV,AMAT)
  490 CONTINUE
      RETURN
      END

C***********************************************************************
C
C***********************************************************************
      SUBROUTINE SERCH(I1,I2,M,MG,IRST,IREN,IPIV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      COMMON /SAI1/ IPIVN0,ICOL0,MAX,IFROW,LAROW,JCOL0,IWRT
      COMMON/MIMAC5/MARGIN(2,6),LMG
      COMMON/PARAM/IBM1,IBM3,IBM5
      DIMENSION IRST(IBM1),IREN(IBM1),IPIV(IBM1)
      NZ=1
      IF(I1.GT.I2) NZ=-1
      I=I1-NZ
      L=M
      MAG=MG
  100 CONTINUE
      J=1
      IF(NZ.LT.0) J=0
      MZ=J
      I=I+NZ
  320 IJ=I+J
      IF(IJ.LE.LAROW) GO TO 340
      JR=IBM5+1
      GO TO 330
  340 JR=IRST(IJ)
      IF(IPIV(IJ).EQ.2) GO TO 300
  330 IF(NZ.GT.0) NP=JR-1-IREN(I)
      IF(NZ.LT.0) NP=IRST(I+1)-1-IREN(IJ)
      IF(NP.LE.0) GO TO 200
      IF(NZ.LT.0) GO TO 110
      MG2=IREN(I)+1
      GO TO 120
  110 CONTINUE
      MG2=IRST(I+1)-1
  120 CONTINUE
      IF(NP.EQ.1) GO TO 130
      LMG=LMG+1
      MARGIN(1,MAG)=I
      MARGIN(2,MAG)=MG2
      L=L-1
      IF(L.LE.0) RETURN
      LMG=LMG+1
      K=MAG+NZ
      MARGIN(1,K)=I
      MARGIN(2,K)=MG2+NZ
      RETURN
  130 CONTINUE
      LMG=LMG+1
      MARGIN(1,MAG)=I
      MARGIN(2,MAG)=MG2
      L=L-1
      MAG=MAG+NZ
      IF(L.LE.0) RETURN
  200 CONTINUE
      I=I+J-MZ
  310 IF(NZ.GT.0.AND.I.LT.I2) GO TO 100
      IF(NZ.LT.0.AND.I.GT.I2) GO TO 100
      RETURN
  300 J=J+NZ
      GO TO 320
      END

C**********************************************************************
C
C**********************************************************************
      SUBROUTINE SLIDE(JR1,JR2,I1,I2,NZ,MG,IVAR,IRST,IREN,ICOL,IPIV,
     +AMAT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      COMMON/MIMAC1/NN,NVP,MAXCOL,IP,JPN
      COMMON/PARAM/IBM1,IBM3,IBM5
      COMMON/SAMAT/BMAT
      DIMENSION IRST(IBM1),IREN(IBM1),ICOL(IBM5),IPIV(IBM1),AMAT(IBM3)
      IF(JR1.GE.JR2) GO TO 300
      I11=I1+1
      I21=I2-1
      IF(I11.GT.I21) GO TO 110
      IF(IP.GT.JR1.AND.IP.LE.JR2) IVAR=IVAR-NZ
      DO 100 I=I11,I21
      K=I-NZ
      ICOL(K)=ICOL(I)
      AMAT(K)=AMAT(I)
  100 CONTINUE
  110 CONTINUE
      JR11=JR1+1
      DO 120 I=JR11,JR2
      IF(IPIV(I).EQ.2) GO TO 120
      IRST(I)=IRST(I)-NZ
      IREN(I)=IREN(I)-NZ
  120 CONTINUE
      INZ=I2-NZ
      GO TO 250
  300 CONTINUE
      I11=I1-1
      IF(I2.GT.I11) GO TO 210
      IF(IP.GT.JR2.AND.IP.LE.JR1) IVAR=IVAR+NZ
      I2I=I2+I11
      DO 200 IJ=I2,I11
      I=I2I-IJ
      KP=I+NZ
      ICOL(KP)=ICOL(I)
      AMAT(KP)=AMAT(I)
  200 CONTINUE
  210 CONTINUE
      JR21=JR2+1
      IF(JR21.GT.JR1) GO TO 230
      DO 220 J=JR21,JR1
      IF(IPIV(J).EQ.2) GO TO 220
      IRST(J)=IRST(J)+NZ
      IREN(J)=IREN(J)+NZ
  220 CONTINUE
  230 CONTINUE
      INZ=I2+NZ-1
  250 CONTINUE
      IF(MG.LE.0) GO TO 240
      ICOL(INZ)=MG
      AMAT(INZ)=BMAT
      IREN(JR2)=IREN(JR2)+1
  240 CONTINUE
      RETURN
      END
C*********************************************************************C
C     EXECUTE ASGN FOR UNIX                                           C
C        NUNIT          :NUMBER OF UNIT                    (I)        C
C        IOST           :1 - FOR READ  2 - FOR WRITE       (I)        C
C                        3 - WORK FILE                                C
C	 DIR		:DIRECTORY			   (I)	      C
C        DOCU           :DOCUMENT TO STANDARD-UNIT         (I)        C
C        KAKU           :EXTENSION OF FILENAME             (I)        C
C*********************************************************************C

      SUBROUTINE ASGUNX(NUNIT,IOST,DIR,DOCU,KAKU)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      CHARACTER DOCU*33,DIR*3,KAKU*3,FNAME*45,FNAME1*45
      DATA FNAME1/'????                                         '/

      IF(IOST.EQ.1.OR.IOST.EQ.2) THEN
        DO 10 I=45,1,-1
          IF(FNAME1(I:I).NE.' ') GOTO 100
   10   CONTINUE
  100   CONTINUE
	DO 5 ID=3,1,-1
	  IF(DIR(ID:ID).NE.' ') GOTO 6
    5   CONTINUE
    6   CONTINUE
        WRITE(*,'(A,'' <'',A,''/'',A,''.'',A,''> ? '',$)')
     &		DOCU,DIR(1:ID),FNAME1(1:I),KAKU
        READ(*,'(A)') FNAME


        IF(FNAME.EQ.' ') THEN
	  FNAME(1:ID)=DIR(1:ID)
	  FNAME(ID+1:ID+1)='/'
          FNAME(ID+2:I+ID+1)=FNAME1(1:I)
          FNAME(I+ID+2:I+ID+2)='.'
          FNAME(I+ID+3:I+ID+5)=KAKU
        ELSE
          DO 20 I=45,1,-1
            IF(FNAME(I:I).NE.' ') GOTO 200
   20     CONTINUE
          STOP ' STANDARD-INPUTS IS WRONG'
  200     CONTINUE
C*** DETECT EXTENSION >>>>>
          NEXT=1
          DO 30 J=I,1,-1
            IF(FNAME(J:J).EQ.'/') GOTO 300
            IF(FNAME(J:J).EQ.'.') GOTO 400
   30     CONTINUE
  300     NEXT=0
  400     CONTINUE
          NEXTD=0
          DO 35 J=1,I
            IF(FNAME(J:J).EQ.'/') IDD=J
   35     CONTINUE
	  IF(IDD.NE.0) NEXTD=1
          FNAME1=FNAME
          IF(NEXT.EQ.0) THEN
            FNAME(I+1:I+1)='.'
            FNAME(I+2:I+4)=KAKU
	    I=I+4
          ELSE
            FNAME1(I-3:I)='    '
          ENDIF
          IF(NEXTD.EQ.0) THEN
	    DO 37 K=I,1,-1
	      FNAME(ID+1+K:ID+1+K)=FNAME(K:K)
   37	    CONTINUE
            FNAME(1:ID)=DIR
            FNAME(ID+1:ID+1)='/'
          ELSE
            FNAME1(1:I-4-IDD)=FNAME(IDD+1:I-4)
	    DO 38 K=I-3-IDD,45
	      FNAME1(K:K)=' '
   38	    CONTINUE
          ENDIF
        ENDIF
      ELSE
        DO 50 I=1,45
          FNAME(I:I)=' '
   50   CONTINUE
        FNAME(1:9)='WORK FILE'
      ENDIF
      IF(IOST.EQ.1)THEN
        OPEN(NUNIT,FILE=FNAME,STATUS='OLD',IOSTAT=IO)
        REWIND NUNIT
      ELSEIF(IOST.EQ.2)THEN
        OPEN(NUNIT,FILE=FNAME,STATUS='UNKNOWN',IOSTAT=IO)
        REWIND NUNIT
      ELSEIF(IOST.EQ.3)THEN
        OPEN(NUNIT,STATUS='SCRATCH',ACCESS='DIRECT',RECL=256,
     &       FORM='UNFORMATTED',IOSTAT=IO)
      ENDIF
      IF(IO.NE.0) THEN
        WRITE(*,'(A,I2,2A)') 'OPEN FAILED :UNIT-NO.=',NUNIT,' ',FNAME
        WRITE(*,*) ('',I=1,10)
        STOP 
      ELSE
        WRITE(*,'(A,I2,2A)') 'OPEN SUCCESS:UNIT-NO.=',NUNIT,' ',FNAME
      ENDIF
      RETURN
      END
C*********************************************************************C
C     MATRIX-SOLVER BY ICCG-METHOD (FOR SYMMETRICAL-SPARSE-MATRIX)    C
C                           SM*X= BB                                  C
C*********************************************************************C

      SUBROUTINE ICCG( SD,SM,BB,X,WW,IRST,IREN,ICOL,
     &                IQ,IR,IL,IU,NITR,EITR,NNDM,NP1,NP2 )
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NIC1=6000,NIC2=200000)
      SAVE
      DIMENSION SD(NP1),SM(NP2),BB(NP1),X(NP1),
     &          WW(NP1,4)
      DIMENSION IRST(NP1),IREN(NP1),ICOL(NP2)
      DIMENSION IQ(NP1),IR(NP1),IL(NP1),IU(NP1)
      DIMENSION RSD(NIC1),SL(NIC2),XX(NIC1)

      NP=NNDM

      DO 10 I=1,NP
        IL(I)=0
        IU(I)=0
        XX(I)=0.D0
   10 CONTINUE
      DO 101 I=1,NP
        WW(I,1)=0.D0
        WW(I,2)=0.D0
        WW(I,3)=0.D0
        WW(I,4)=0.D0
  101 CONTINUE

      L=1
      L0=1
      DO 20 I=1,NP
        RSD(I)=SQRT(SD(I))
        SD(I) =1.D0
        I1=IRST(I)
        I2=IREN(I)
        DO 21 J=I1,I2
          ICOL(L)=ICOL(J)
          SM(L)=SM(J)
          L=L+1
   21   CONTINUE
        IRST(I)=L0
        IREN(I)=L-1
        L0=L
   20 CONTINUE

      JD=0
      DO 200 II=1,NP
        BB(II)=BB(II)/RSD(II)
        IF(IRST(II).LE.IREN(II)) THEN
          JD=JD+1
          IR(JD)=II
          IQ(II)=JD
        ELSE
          IQ(II)=0
          BB(II)=BB(II)/RSD(II)
        ENDIF
  200 CONTINUE

      NPP=JD

      DO 210 IJ=1,NPP
        II=IR(IJ)
        DO 220 IK=IRST(II),IREN(II)
          IF(ICOL(IK).GT.II) THEN
            IU(IJ)=IU(IJ)+1
          ELSE
            IL(IJ)=IL(IJ)+1
          ENDIF
  220   CONTINUE
  210 CONTINUE
      DO 121 II=1,NP
        I1=IRST(II)
        I2=IREN(II)
        DO 122 IJ=I1,I2
          J1=ICOL(IJ)
          SM(IJ)=SM(IJ)/(RSD(II)*RSD(J1))
          SL(IJ)=SM(IJ)
  122   CONTINUE
  121 CONTINUE

       IF(NPP.NE.0) THEN
           CALL ILUCG(SD,SM,SL,BB,XX,WW,IRST,IREN,ICOL,IQ,IR,IL,IU,
     &               NPP,NITR,EITR,NP,NP1,NP2,NIC1,NIC2)
       ENDIF

      DO 30 I=1,NPP
        JK=IR(I)
        X(JK)=XX(I)/RSD(JK)
   30 CONTINUE
      RETURN
      END
C*********************************************************************C
C     MAIN PROCESSING OF ICCG                                         C
C*********************************************************************C
      SUBROUTINE ILUCG(SD,SM,SL,BB,XX,WW,IRST,IREN,ICOL,IQ,IR,IL,IU,
     &                 NPP,NITR,EITR,NP,NP1,NP2,NIC1,NIC2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      DIMENSION SD(NP1),SM(NP2),BB(NP1),WW(NP1,4)
      DIMENSION IRST(NP1),IREN(NP1),ICOL(NP2),IQ(NP1)
      DIMENSION IR(NP1),IL(NP1),IU(NP1)
      DIMENSION SL(NIC2),XX(NIC1)

      CALL DECOMP(SD,SL,NPP,WW(1,1),
     &            IL,IR,IRST,IREN,ICOL,IQ,IU,NP1,NP2,NIC2)
      N=NPP
      BNORM=0.

      DO 10 I=1,N
        JK=IR(I)
        BNORM=BNORM+ABS(BB(JK))**2
        WW(I,2)=BB(JK)-WW(I,2)
   10 CONTINUE

        IF(BNORM.EQ.0.) THEN
          DO 11 I=1,N
 11        XX(I)=0.
           NITR=1
           EITR=0.
          RETURN
        ENDIF


      CALL LUSUB(WW(1,1),SL,IRST,IREN,ICOL,IR,IL,NPP,
     &           WW(1,2),WW(1,3),IQ,IU,NP1,NP2,NIC1)

      RR1=0.D0
      DO 60 I=1,N
        RR1=RR1+WW(I,3)*WW(I,2)
   60 CONTINUE

      DO 100 L=1,NITR
           CALL AXMLT(SD,SM,NPP,WW(1,3),WW(1,4),IR,IRST,
     &               IREN,ICOL,IQ,NP1,NP2)
           RR2=0.D0
           DO 120 I=1,N
              RR2=RR2+WW(I,3)*WW(I,4)
  120      CONTINUE
           ALP=RR1/RR2
           RNORM=0.
           DO 140 I=1,N
             XX(I)=XX(I)+ALP*WW(I,3)
             WW(I,2)=WW(I,2)-ALP*WW(I,4)
             RNORM=RNORM+ABS(WW(I,2))**2
  140      CONTINUE

           ERR=SQRT(RNORM/BNORM)

             WRITE(*,'('' N='',I4,''  ERROR='',E12.4)') L,ERR

             IF(ERR.LT.EITR)  GOTO 1100

              CALL LUSUB(WW(1,1),SL,IRST,IREN,ICOL,IR,IL,NPP,
     &                 WW(1,2),WW(1,4),IQ,IU,NP1,NP2,NIC1)
              RR3=0.D0
              DO 160 I=1,N
                 RR3=RR3+WW(I,4)*WW(I,2)
  160         CONTINUE
              BETA=RR3/RR1
              DO 180 I=1,N
                 WW(I,3)=WW(I,4)+BETA*WW(I,3)
  180         CONTINUE
              RR1=RR3
  100 CONTINUE
 1100 CONTINUE
      RETURN
      END
C*********************************************************************C
C     PRE-PROCESSING (MATRIX DE-COMPOSITION)                          C
C*********************************************************************C
      SUBROUTINE DECOMP(SD,SL,NPP,DD,
     &                  IL,IR,IRST,IREN,ICOL,IQ,IU,NP1,NP2,NIC2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      DIMENSION SD(NP1),DD(NP1)
      DIMENSION ICOL(NP2)
      DIMENSION IRST(NP1),IREN(NP1),IL(NP1),IR(NP1),IQ(NP1),IU(NP1)
      DIMENSION SL(NIC2)

      DO 100 L=1,NPP
        LR=IR(L)
        DD(L)=SD(LR)
        DO 110 K=IRST(LR),IRST(LR)+IL(L)-1
          KQ=IQ(ICOL(K))
          DD(L)=DD(L)-SL(K)*SL(K)
  110   CONTINUE
        IF(IU(L).EQ.0)                             GOTO 100
        DO 200 J=IRST(LR)+IL(L),IREN(LR)
          JJ=ICOL(J)
          JJQ=IQ(JJ)
          IF(IL(L).EQ.0.OR.IL(JJQ).EQ.0)           GOTO 250
          DO 210 K=IRST(LR),IRST(LR)+IL(L)-1
            KK=ICOL(K)
            KKQ=IQ(KK)
            IF(KK.GE.LR)                           GOTO 250
            DO 220 I=IRST(JJ),IRST(JJ)+IL(JJQ)-1
              II=ICOL(I)
              IF(II.GE.LR)                         GOTO 210
              IF(II.EQ.KK) THEN
                SL(J)=SL(J)-SL(K)*SL(I)
                                                   GOTO 220
              END IF
  220       CONTINUE
  210     CONTINUE
  250     SL(J)=SL(J)/DD(L)
          DO 300 I=IRST(JJ),IRST(JJ)+IL(JJQ)-1
            IF(ICOL(I).EQ.LR)                      GOTO 310
  300     CONTINUE
          STOP 'DECOMP -2'
  310     SL(I)=SL(J)*DD(L)
  200   CONTINUE
  100 CONTINUE
      DO 400 I=1,NPP
        DD(I)=1./DD(I)
  400 CONTINUE
      RETURN
      END
C*********************************************************************C
C     MLUTYPLY MATRIX & VECTOR FOR NORM-CALCULATION                   C
C*********************************************************************C
      SUBROUTINE AXMLT(SD,SM,NPP,XX,YY,IR,IRST,
     &                 IREN,ICOL,IQ,NP1,NP2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      DIMENSION SD(NP1),SM(NP2),XX(NP1),YY(NP1),IQ(NP1)
      DIMENSION IR(NP1),IRST(NP1),IREN(NP1),ICOL(NP2)
      DO 10 I=1,NPP
        JK=IR(I)

        YY(I)=SD(JK)*XX(I)
        DO 11 J=IRST(JK),IREN(JK)
          JL=IQ(ICOL(J))
          YY(I)=YY(I)+SM(J)*XX(JL)
   11   CONTINUE
   10 CONTINUE
      RETURN
      END
C*********************************************************************C
C     NANISHITERUKA SHIRANAI                                          C
C*********************************************************************C
      SUBROUTINE LUSUB(DD,SM,IRST,IREN,ICOL,IR,IL,NPP,
     &                 BB,XX,IQ,IU,NP1,NP2,NIC1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
      DIMENSION DD(NP1),SM(NP2),BB(NP1)
      DIMENSION IRST(NP1),IREN(NP1),ICOL(NP2)
      DIMENSION IR(NP1),IL(NP1),IU(NP1),IQ(NP1)
      DIMENSION XX(NIC1)
      DO 10 L=1,NPP
        JK=IR(L)
        XX(L)=BB(L)
        DO 12 K=IRST(JK),IRST(JK)+IL(L)-1
          JL=IQ(ICOL(K))
          XX(L)=XX(L)-SM(K)*XX(JL)
   12   CONTINUE
        XX(L)=DD(L)*XX(L)
   10 CONTINUE
      DO 20 L=NPP-1,1,-1
        JK=IR(L)
        IF(IU(L).EQ.0)                             GOTO 20
        DO 24 K=IRST(JK)+IL(L),IREN(JK)
          JL=IQ(ICOL(K))
          XX(L)=XX(L)-SM(K)*XX(JL)
   24   CONTINUE
   20 CONTINUE
      RETURN
      END

