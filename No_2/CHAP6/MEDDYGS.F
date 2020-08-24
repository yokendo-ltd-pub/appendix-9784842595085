      SUBROUTINE GAUS( N )
      PARAMETER (NPPN=4000,NEEN=8000,NVVN=4000,IFMX=8)
      implicit real ( A-H, O-Z )
      COMPLEX*8 AA,AMX,DMM,T,BB,ANS,TMP
      COMMON /FLD/AA(IFMX,2*NVVN+1),NUNK,NBLK,IAMR
      COMMON /RES/ANS(NVVN+NVVN),BB(IFMX,2*NVVN+1)
      ND=N
      NLP=NBLK
      IF( IAMR .NE. 0 ) NLP = NBLK + 1
      WRITE(*,'(1H ,40X,I4,''/'',I4)')N,IFMX
      DO 10 MRD=1,NLP
         WRITE(*,'(1H+,20X,A10,I4)')'BLOCK No.',MRD
         READ(35,REC=MRD)((BB(II,JJ),II=1,IFMX),JJ=1,ND+1)
         NLP2=IFMX
         IF( MRD.GT.NBLK)NLP2=IAMR
         DO 30 MRF=1,MRD
            READ(35,REC=MRF)((AA(II,JJ),II=1,IFMX),JJ=1,ND+1)
            IF( MRF .NE. MRD )THEN
               DO 50 K=1,IFMX
                  K1=IFMX*(MRF-1)+K
                  AMX=AA(K,K1)
                  IF( CABS( AMX).LT. 1.E-20) PAUSE
                  DO 50 I=1,IFMX
                     TMP=BB(I,K1)/AMX
                     DO 50 J=K1,ND+1
                        BB(I,J)=BB(I,J)-AA(K,J)*TMP
   50          CONTINUE
            ELSE
               DO 55 K=1,NLP2
                  K1=IFMX*(MRF-1)+K
                  AMX=BB(K,K1)
                  IF( CABS( AMX).LT. 1.E-20) PAUSE
                  DO 55 I=K+1,NLP2
                     TMP=BB(I,K1)/AMX
                     DO 55 J=K1,ND+1
                        BB(I,J)=BB(I,J)-BB(K,J)*TMP
   55          CONTINUE
            ENDIF
   30    CONTINUE
         WRITE(35,REC=MRD)((BB(II,JJ),II=1,IFMX),JJ=1,ND+1)
   10 CONTINUE
      WRITE(*,'(1H ,20X,A10)') 'BACKFORD'
      DO 60 MRD=NLP,1,-1
         READ(35,REC=MRD)((AA(II,JJ),II=1,IFMX),JJ=1,ND+1)
         IF( IAMR .NE. 0 .AND. MRD.EQ.NLP) THEN
            DO 70 K=IAMR,1,-1
               K1 = ( MRD -1 )*IFMX + K
               ANS( K1 ) = AA( K, ND+1)
               DO 80 I=K1+1,ND
   80             ANS( K1 ) = ANS( K1 ) - AA(K,I)*ANS(I)
               ANS( K1 ) = ANS( K1 ) / AA( K,K1)
   70       CONTINUE
         ELSE
            DO 90 K=IFMX,1,-1
               K1 = ( MRD -1 )*IFMX + K
               ANS( K1 ) = AA( K, ND+1)
               DO 95 I=K1+1,ND
   95             ANS( K1 ) = ANS( K1 ) - AA(K,I)*ANS(I)
               ANS( K1 ) = ANS( K1 ) / AA( K,K1)
   90       CONTINUE
         ENDIF
   60 CONTINUE
      RETURN
      END
C
      SUBROUTINE TSYMCOND( JJ, JORG, ASGN, BSGN )
C
      PARAMETER (NPPN=4000,NEEN=8000,NVVN=4000,IFMX=8)
C
      implicit real ( A-H, O-Z )
      COMPLEX*8 AA,BB,ANS
      COMMON /TR1/NP,NE3,ZZ(NPPN,3),NOE3(NEEN,6),ZNP(NEEN,3)
      COMMON /VCT/ZLN(NVVN,6),NOV(NVVN,2),NV
      COMMON /TR4/BN(NEEN,3),T1(NEEN,3),T2(NEEN,3),ZI(NEEN,3)
      COMMON /B00/B0(3*NEEN,3),A0(3*NEEN,3)
      COMMON /SYM/NPOLD,NE3OLD,NVOLD,IISYM(8,3),ITBL(16,3*NEEN)
      COMMON /SBE/NXYZ,IESYM(8,3),IBSYM(8,3)
      COMMON /FLD/AA(IFMX,2*NVVN+1),NUNK,NBLK,IAMR
      COMMON /RES/ANS(NVVN+NVVN),BB(IFMX,2*NVVN+1)
C
      JORG=0
      BSGN=1.
      ASGN=1.
      DO 10 I=2,NXYZ
         IF( ITBL(I,JJ) .NE. 0 )THEN
            IORG=ITBL(I,JJ)
            JORG=IORG
            JSYM=JJ
            TX=ZLN(JSYM,1)
            TY=ZLN(JSYM,2)
            TZ=ZLN(JSYM,3)
            BORGX=ZLN(IORG,1)*IBSYM(I,1)
            BORGY=ZLN(IORG,2)*IBSYM(I,2)
            BORGZ=ZLN(IORG,3)*IBSYM(I,3)
            EORGX=ZLN(IORG,1)*IESYM(I,1)
            EORGY=ZLN(IORG,2)*IESYM(I,2)
            EORGZ=ZLN(IORG,3)*IESYM(I,3)
            FLG = BORGX*TX+BORGY*TY+BORGZ*TZ
            FLGE= EORGX*TX+EORGY*TY+EORGZ*TZ
            IF( FLGE .GT. 0 )THEN
               ASGN=1.
            ELSE
               ASGN=-1.
            ENDIF
            IF( FLG .GT. 0 )THEN
               BSGN=1.
            ELSE
               BSGN=-1.
            ENDIF
         ENDIF
   10 CONTINUE
      IF( JORG .EQ. 0 ) JORG=JJ
      RETURN
      END

