************************************************************************
*                                                                      *
      SUBROUTINE INIZIO
*                                                                      *
*     Questo sotto-programma legge i dati in ingresso e prepara le     *
*     grandezze necessarie all elaborazione                            *
*                                                                      *
************************************************************************

***** PARAMETRI
*     numero massimo di particelle
      INTEGER,PARAMETER::NPMAX=200000

***** COSTANTI MATEMATICHE
*     pi greco
      REAL*8,PARAMETER::PI=3.1415926535897932D+0

***** GRANDEZZE MICROSCOPICHE
*     posizione delle particelle
      REAL*8 X(NPMAX)
*     componenti della velocita delle particelle
      REAL*8 VX(NPMAX),VY(NPMAX),VZ(NPMAX)
*     energia interna delle particelle
      REAL*8 EI(NPMAX)
*      numero totale di particelle
      INTEGER NPAR
      COMMON /MICRO/ X,VX,VY,VZ,EI,NPAR

***** STATO DELLA PARETE
*     velocita relativa
      REAL*8 UW
*     temperatura
      REAL*8 TMPW
      COMMON /PARETE/ UW,TMPW

***** DOMINIO SPAZIALE
*     estremi del dominio
      REAL*8 XMIN,XMAX
*     passo spaziale
      REAL*8 DX
*     numero totale di nodi
      INTEGER NNOD
      COMMON /SPAZIO/ XMIN,XMAX,DX,NNOD

***** DOMINIO TEMPORALE
*     istante finale
      REAL*8 TMAX
*     passo temporale
      REAL*8 DT
*     istante d inizio campionamento grandezze macroscopiche
      REAL*8 TIM
*     passo campionamento
      REAL*8 DTM
      COMMON /TEMPO/ TMAX,DT,TIM,DTM

***** CARATTERISTICHE DEL GAS
*     rapporto calori specifici
      REAL*8 GAMMA
*     numero di gradi di liberta interni
      INTEGER ZETA
      COMMON /GAS/ GAMMA,ZETA

***** PARAMETRI COLLISIONALI
*     sezione d urto
      REAL*8 XSECT
*     frazione di urti anaelastici
      REAL*8 LAMBDA
      COMMON /COLLIS/ XSECT,LAMBDA

***** PARAMETRI SIMULAZIONE
*     numero iniziale di particelle per elemento
      INTEGER NPE0
      COMMON /SIM/ NPE0

***** VARIABILI LOCALI
*     stringa per la descrizione degli input
      CHARACTER COMMEN (58)
*     numero di Knudsen
      REAL*8 KN
*     numero di nodi lungo x
      INTEGER NX
*     seme casuale
      INTEGER SEME
*     quadrato del diametro delle particelle
      REAL*8 DIAM2
*     indice generico
      INTEGER J
*     indice particelle
      INTEGER JP

***** LETTURA DATI DI INGRESSO

*     apertura file di input
      OPEN(UNIT=1,FILE='couette.inp',ACCESS='SEQUENTIAL',STATUS='OLD')

*     apertura file di output
      OPEN(UNIT=2,FILE='couette.out',ACCESS='SEQUENTIAL',
     &     STATUS='UNKNOWN')

*     lettura file input e copiatura su output
      READ (1,100)(COMMEN (J),J=1,58)
      WRITE(2,100)(COMMEN(J),J=1,58)
      READ (1,100)(COMMEN (J),J=1,58)
      WRITE(2,100)(COMMEN(J),J=1,58)

      READ (1,300)(COMMEN (J),J=1,58),KN
      WRITE(2,300)(COMMEN(J),J=1,58),KN
      READ (1,300)(COMMEN (J),J=1,58),UW
      WRITE(2,300)(COMMEN(J),J=1,58),UW
      READ (1,300)(COMMEN (J),J=1,58),TMPW
      WRITE(2,300)(COMMEN(J),J=1,58),TMPW
      READ (1,200)(COMMEN (J),J=1,58),NPE0
      WRITE(2,200)(COMMEN(J),J=1,58),NPE0
      READ (1,300)(COMMEN (J),J=1,58),LAMBDA
      WRITE(2,300)(COMMEN(J),J=1,58),LAMBDA

      READ (1,100)(COMMEN (J),J=1,58)
      WRITE(2,100)(COMMEN(J),J=1,58)
      READ (1,100)(COMMEN (J),J=1,58)
      WRITE(2,100)(COMMEN(J),J=1,58)
      READ (1,100)(COMMEN (J),J=1,58)
      WRITE(2,100)(COMMEN(J),J=1,58)

      READ (1,200)(COMMEN (J),J=1,58),NX
      WRITE(2,200)(COMMEN(J),J=1,58),NX

      READ (1,100)(COMMEN (J),J=1,58)
      WRITE(2,100)(COMMEN(J),J=1,58)
      READ (1,100)(COMMEN (J),J=1,58)
      WRITE(2,100)(COMMEN(J),J=1,58)
      READ (1,100)(COMMEN (J),J=1,58)
      WRITE(2,100)(COMMEN(J),J=1,58)

      READ (1,300)(COMMEN (J),J=1,58),DT
      WRITE(2,300)(COMMEN(J),J=1,58),DT
      READ (1,300)(COMMEN (J),J=1,58),TMAX
      WRITE(2,300)(COMMEN(J),J=1,58),TMAX
      READ (1,300)(COMMEN (J),J=1,58),TIM
      WRITE(2,300)(COMMEN(J),J=1,58),TIM
      READ (1,300)(COMMEN (J),J=1,58),DTM
      WRITE(2,300)(COMMEN(J),J=1,58),DTM
      READ (1,200)(COMMEN (J),J=1,58),SEME
      WRITE(2,200)(COMMEN(J),J=1,58),SEME

      READ (1,100)(COMMEN (J),J=1,58)
      WRITE(2,100)(COMMEN(J),J=1,58)
      READ (1,100)(COMMEN (J),J=1,58)
      WRITE(2,100)(COMMEN(J),J=1,58)
      READ (1,100)(COMMEN (J),J=1,58)
      WRITE(2,100)(COMMEN(J),J=1,58)

*     definisco gli estremi del dominio
*     hp : libero cammino medio=1
      XMIN=-0.5D+0/ KN
      XMAX=+0.5D+0/ KN

*     definisco la griglia
      NNOD=NX
      DX=(XMAX-XMIN)/ NX

*     calcolo sezione d urto
      DIAM2=1.D+0/(DSQRT (2.D+0)*PI*NPE0/DX)
      XSECT=PI*DIAM2

*     rapporto calori specifici
      IF(LAMBDA.EQ .0.D+0) THEN
         GAMMA=5.D+0/3.D+0
         ZETA=2
      ELSE
         GAMMA=7.D+0/5.D+0
         ZETA=0
      END IF

*     numero di particelle
      NPAR=NPE0*NNOD
      WRITE(2,400)
     & '      NUMERO DI PARTICELLE INIZIALE .....................',NPAR

*     controllo che il numero di particelle non superi NPMAX
      IF(NPAR.GT.NPMAX) GOTO 1000

*     genero paritcelle con distribuzione della posizione uniforme lungo
*     x e distribuzione della velocita maxwelliana
      CALL SRAND(SEME)
      DO JP=1,NPAR
         X(JP)=XMIN +(XMAX-XMIN)*RAND()
         CALL MAXWELL(TMPW,VX(JP),VY (JP),VZ(JP),EI(JP))
      END DO

      RETURN

***** FORMATI **********************************************************

100   FORMAT (58(A1))
200   FORMAT (58(A1),I14)
300   FORMAT (58(A1),E14.4)
400   FORMAT(A58,I14)

***** ERRORE NPAR > NPMAX **********************************************

1000  WRITE(2,1100)
     & '
      '
      WRITE(2,1100)
     & '***** ERROR ********************************************* '
      WRITE(2,1100)
     & '                                                          '
      WRITE(2,1200)
     & '      NUMERO DI PARTICELLE > ',NPMAX,'                    '
      WRITE(2,1100)
     & '                                                          '
      WRITE(2,1100)
     & '***** STOP ********************************************** '
      WRITE(*,1300) 'error : NPAR > ',NPMAX
      WRITE(*,*) 'Aborted'

      STOP

***** FINE ERRORE ******************************************************

1100  FORMAT(A58)
1200  FORMAT(A29,I9.9,A20)
1300  FORMAT(A14,I9.9)

      END

***** FINE SUBROUTINE INIZIO *******************************************
