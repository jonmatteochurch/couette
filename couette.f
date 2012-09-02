************************************************************************
*                                                                *
*     Programma per la simulazione diretta del flusso di calore tra    *
*     due piani paralleli in moto relativo                       *
*                                                                *
      PROGRAM COUETTE
*                                                                *
*     Ultima modifica ore 17 del 3 Settembre 2012                *
*     comprende i file:                                          *
*     - collisioni.f                                             *
*     - inizio.f                                                 *
*     - maxwell.f                                                *
*     - velpar.f                                                 *
*                                                                *
************************************************************************

***** VARIABILI ********************************************************

***** PARAMETRI
*     numero massimo di particelle
*     (deve coincidere col valore in inizio.f e collisioni.f)
      INTEGER,PARAMETER::NPMAX=200000
*     numero massimo di nodi
*     (deve coincidere col valore in collisioni.f)
      INTEGER,PARAMETER::NNMAX=1000

***** GRANDEZZE MICROSCOPICHE
*     posizione delle particelle
      REAL*8 X(NPMAX)
*     componenti della velocita delle particelle
      REAL*8 VX(NPMAX),VY(NPMAX),VZ(NPMAX)
*     energia interna delle particelle
      REAL*8 EI(NPMAX)
*     numero totale di particelle
      INTEGER NPAR
      COMMON /MICRO/ X,VX,VY,VZ,EI,NPAR

***** STATO DELLA PARETE
*     velocita relativa pareti
      REAL*8 UW
*     temperatura
      REAL*8 TMPW
      COMMON /PARETE/ UW,TMPW

***** DOMINIO SPAZIALE
*     estremi del dominio
      REAL*8 XMIN,XMAX
*     passo spaziali
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

***** PARAMETRI SIMULAZIONE
*     numero iniziale di particelle per elemento
      INTEGER NPE0
      COMMON /SIM/ NPE0

***** INDICI PARTICELLE/ELEMENTI
*     indici delle particelle ordinate per elemento di appartenenza
      INTEGER IP(NPMAX)
*     indici delle prime particelle di ciascun elemento
      INTEGER IPIN(NNMAX)
*     numero di particelle per elemento della griglia
      INTEGER NPE(NNMAX)
      COMMON /IND/ IP,IPIN,NPE

***** VARIABILI LOCALI
*     numero totale di particelle campionate per elemento
      REAL*8NPME(NNMAX)
*     componenti della velocita molecolare media e loro prodotti
      REAL*8 VXMED(NNMAX),VYMED(NNMAX),VZMED(NNMAX)
      REAL*8 VX2MED(NNMAX),VY2MED(NNMAX),VZ2MED(NNMAX)
      REAL*8 VXVYMED(NNMAX),VXVZMED(NNMAX),VYVZMED(NNMAX)
*     modulo quadro medio medio della velocita molecolare e sue
*     combinazioni
      REAL*8 V2MED(NNMAX),VXV2MED(NNMAX),VYV2MED(NNMAX),VZV2MED(NNMAX)
*     densita dei particelle
      REAL*8 N(NNMAX)
*     temperatura complessiva,traslazionale e rotazionale
      REAL*8 TMP(NNMAX),TMPTR(NNMAX),TMPROT(NNMAX)
*     energia interna media
      REAL*8 EIMED(NNMAX)
*     tensore degli sforzi e pressione scalare
      REAL*8 PXX(NNMAX),PYY(NNMAX),PZZ(NNMAX),PXY(NNMAX),PXZ(NNMAX)
      REAL*8 PYZ(NNMAX),P(NNMAX)
*     flusso di molecole
      REAL*8 FX(NNMAX),FY(NNMAX),FZ(NNMAX)
*     flusso di calore
      REAL*8 QX(NNMAX),QY(NNMAX),QZ(NNMAX)
*     indice particelle
      INTEGER JP
*     indice nodi
      INTEGER JN
*     indice elementi
      INTEGER JE
*     istante temporale corrente
      REAL*8 T
*     prossimo istante temporale di campionamento grandezze
*     macroscopiche
      REAL*8 TM
*     posizione raggiunta dalla particella dopo DT
      REAL*8 XNEW
*     tempo trascorso fra l istante precedente e l urto con una parete
      REAL*8 DTCOLL
*     numero campionamenti
      INTEGER NM
*     indice particella ordinata per elemento
      INTEGER JPNEW
*     quadrato delle componenti della velocita media di una particella
      REAL*8 VXMED2,VYMED2,VZMED2
*     modulo quadro della velocita delle particelle
      REAL*8 V2
*     modulo quadro della velocita media del gas
      REAL*8 VMED2
*     coordinata x del punto medio di un elemento
      REAL*8 XE
*     carriage return
      CHARACTER,PARAMETER::CR=ACHAR (13)
*     numero medio particelle
      INTEGER NPMED
*     progresso della simulazione in percentuale
      REAL PROGR
*     indici degli elementi della griglia a cui appartengono le
*     particelle
      INTEGER IE(NPMAX)

***** FINE VARIABILI ***************************************************

***** INIZIO inizializza
*     - le grandezze microscopiche X,VX,VY,VZ,EI,NPAR
*     - lo stato della parete UW,TMPW
*     - il dominio spaziale XMIN,XMAX,DX,NNOD
*     - il dominio temporale TMAX,DT,TIM,DTM
*     - le caratteristiche del gas GAMMA,ZETA
*     - i parametri collisionali XSECT,LAMBDA
*     - i parametri della simulazione NPE0

      WRITE(*,100,ADVANCE='no') ' Inizializzazione... '
      CALL INIZIO
      WRITE(*,*) 'ok'

***** PASSO INIZIALE

      PROGR=0.E+1
      WRITE(*,200,ADVANCE=' no') CR,PROGR

      T=0.D+0
      TM=TIM
      NM=0
      NPMED=0

      DO JE=1,NNOD
         VXMED(JE)=0.D+0
         VYMED(JE)=0.D+0
         VZMED(JE)=0.D+0
         VX2MED(JE)=0.D+0
         VY2MED(JE)=0.D+0
         VZ2MED(JE)=0.D+0
         VXVYMED(JE)=0.D+0
         VXVZMED(JE)=0.D+0
         VYVZMED(JE)=0.D+0
         V2MED(JE)=0.D+0
         VXV2MED(JE)=0.D+0
         VYV2MED(JE)=0.D+0
         VZV2MED(JE)=0.D+0
      END DO

***** CICLO TEMPORALE

      DO WHILE(T.LT.TMAX)

         T=T+DT

*****    movimento particelle

         DO JP=1,NPAR
*           aggiorno posizione
            XNEW=X(JP)+VX(JP)*DT
            IF(XNEW.LE.XMIN) THEN
*              particella urta parete inferiore
               DTCOLL=(XMIN-X(JP))/ VX(JP)
               CALL VELPAR (-UW/2.D+0,TMPW,VX(JP),VY(JP),VZ(JP),EI(JP))
               XNEW=XMIN+VX(JP)*(DT-DTCOLL)
            END IF
            IF(XNEW.GT.XMAX) THEN
*              particella urta parete superiore
               DTCOLL=(XMAX-X(JP))/ VX(JP)
               CALL VELPAR(UW/2.D+0,TMPW,VX(JP),VY(JP),VZ(JP),EI(JP))
               VX(JP)=-VX(JP)
               XNEW=XMAX+VX(JP)*(DT-DTCOLL)
            END IF
            X(JP)=XNEW
         END DO

*****    elimino le particelle finite fuori del dominio

*        ciclo all indietro sulle particelle
         JP=NPAR
         DO WHILE(JP.GE.1)
            IF(X(JP).LT.XMIN.OR.X(JP).GE.XMAX) THEN
*              particella fuori dal dominio
*              le sostituisco l ultima in posizione NPAR e aggiorno il 
*              valore di NPAR
               X(JP)=X(NPAR)
               VX(JP)=VX(NPAR)
               VY(JP)=VY(NPAR)
               VZ(JP)=VZ(NPAR)
               EI(JP)=EI(NPAR)
               NPAR=NPAR -1
            END IF
*           aggiorno indice particella
            JP=JP -1
         END DO

*****    riordino le particelle

*        azzero il numero di particelle per elemento
         DO JN=1,NNOD
            NPE(JN)=0
         END DO

*        aggiorno gli indici degli elementi della griglia a cui
*        appartengono le particelle ed il numero di particelle per
*        elemento
         DO JP=1,NPAR
            IE(JP)=1+INT ((X(JP)-XMIN)/ DX)
            NPE(IE(JP))=NPE(IE(JP))+1
         END DO

*        aggiorno il vettore degli indici relativi alla prima particella
*        di ciascun elemento e contestualmente pongo a zero il numero di
*        particelle per elemento
         IPIN (1)=0
         DO JN=2,NNOD
            IPIN(JN)=IPIN(JN -1)+NPE(JN -1)
            NPE(JN -1)=0
         END DO
         NPE(NNOD)=0
 
*        calcolo i nuovi indici di ciascuna particella in modo che siano
*        ordinati per elemento di appartenenza e aggiorno il numero di
*        particelle per elemento
         DO JP=1,NPAR
            JE=IE(JP)
            NPE(JE)=NPE(JE)+1
            JPNEW=IPIN(JE)+NPE(JE)
            IP(JPNEW)=JP
         END DO

*****    calcolo le collisioni

         DO JE=1,NNOD
            IF(NPE(JE). GT .1) THEN
               CALL COLLISIONI(JE,T)
            END IF
         END DO

*****    campiono le quantita macroscopiche

         IF(T.GE.TM) THEN
*           aggiorno numero di particelle medio
            NPMED=NPMED+NPAR
*           aggiorno numero di campionamenti
            NM=NM +1
            DO JP=1,NPAR
*
               JE=IE(JP)
              NPME(JE)=NPME(JE)+1

               VXMED(JE)=VXMED(JE)+VX (JP)
               VYMED(JE)=VYMED(JE)+VY (JP)
               VZMED(JE)=VZMED(JE)+VZ (JP)

               VX2MED(JE)=VX2MED(JE)+VX(JP)**2
               VY2MED(JE)=VY2MED(JE)+VY(JP)**2
               VZ2MED(JE)=VZ2MED(JE)+VZ(JP)**2

               VXVYMED(JE)=VXVYMED(JE)+VX(JP)*VY(JP)
               VXVZMED(JE)=VXVZMED(JE)+VX(JP)*VZ(JP)
               VYVZMED(JE)=VYVZMED(JE)+VY(JP)*VZ(JP)

               V2=VX(JP)**2+VY(JP)**2+VZ(JP)**2
               V2MED(JE)=V2MED(JE)+V2
               VXV2MED(JE)=VXV2MED(JE)+VX(JP)*V2
               VYV2MED(JE)=VYV2MED(JE)+VY(JP)*V2
               VZV2MED(JE)=VZV2MED(JE)+VZ(JP)*V2
*
               EIMED(JE)=EIMED(JE)+EI (JP)
            END DO
*           aggiorno prossimo istante per campionamento
            TM=TM+DTM
         END IF

         IF(INT(T/TMAX *1000).GT.PROGR*10) THEN
            PROGR=PROGR +1.0E-1
            WRITE(*,200,ADVANCE='no') CR,PROGR
         END IF

      END DO

      WRITE(*,200) CR,100 E +0
      WRITE(*,300,ADVANCE=' no')' Finalizzazione ... '

      WRITE(2,400)
     & '      NUMERO DI CAMPIONAMENTI ........................... ',Nm
      WRITE(2,500)
     & '      NUMERO MEDIO DI PARTICELLE ........................ ',
     & DFLOAT(NPMED)/ NM
      WRITE(2,600)
     & '                                                          ',

      DO JE=1,NNOD
         IF(NPME(JE).NE.0.D+0) THEN
*           scalo per il numero di campionamenti effettuati su ciascun
*           elemeneto i valori campionari cumulati delle combinazioni
*           delle velocita
            VXMED(JE)=VXMED(JE)/NPME(JE)
            VYMED(JE)=VYMED(JE)/NPME(JE)
            VZMED(JE)=VZMED(JE)/NPME(JE)

            VX2MED(JE)=VX2MED(JE)/NPME(JE)
            VY2MED(JE)=VY2MED(JE)/NPME(JE)
            VZ2MED(JE)=VZ2MED(JE)/NPME(JE)

            VXVYMED(JE)=VXVYMED(JE)/NPME(JE)
            VXVZMED(JE)=VXVZMED(JE)/NPME(JE)
            VYVZMED(JE)=VYVZMED(JE)/NPME(JE)

            V2MED(JE)=V2MED(JE)/NPME(JE)
            VXV2MED(JE)=VXV2MED(JE)/NPME(JE)
            VYV2MED(JE)=VYV2MED(JE)/NPME(JE)
            VZV2MED(JE)=VZV2MED(JE)/NPME(JE)

            EIMED(JE)=EIMED(JE)/NPME(JE)
         END IF
         VXMED2=VXMED(JE)**2
         VYMED2=VYMED(JE)**2
         VZMED2=VZMED(JE)**2
         VMED2=VXMED2+VYMED2+VZMED2
*        densita di molecole
         N(JE)=NPME(JE)/(NM*NPE0)
*        flusso di molecole
         FX(JE)=N(JE)*VXMED(JE)
         FY(JE)=N(JE)*VYMED(JE)
         FZ(JE)=N(JE)*VZMED(JE)
*        tensore degli sforzi
         PXX(JE)=N(JE)*(VX2MED(JE)-VXMED2)
         PYY(JE)=N(JE)*(VY2MED(JE)-VYMED2)
         PZZ(JE)=N(JE)*(VZ2MED(JE)-VZMED2)
         PXY(JE)=N(JE)*(VXVYMED(JE)-VXMED(JE)*VYMED(JE))
         PXZ(JE)=N(JE)*(VXVZMED(JE)-VXMED(JE)*VZMED(JE))
         PYZ(JE)=N(JE)*(VYVZMED(JE)-VYMED(JE)*VZMED(JE))
*        pressione scalare
         P(JE)=(PXX(JE)+PYY(JE)+PZZ(JE))/3.D+0
*        temperatura traslazionale,rotazionale e totale
         TMPTR(JE)=(V2MED(JE)-VMED2)/3.D+0
         TMPROT(JE)=2.D+0*EIMED(JE)/ZETA
         TMP(JE)=(3.D+0*TMPTR(JE)+ZETA*TMPROT(JE))/(3.D+0+ZETA)
*        flusso di calore
         QX(JE)=N(JE)*(VXV2MED(JE)/2.D+0-VXMED(JE)*VX2MED(JE)-
     &                 VYMED(JE)*VXVYMED(JE)-VZMED(JE)*VXVZMED(JE)+
     &                 VXMED(JE)*VMED2-VXMED(JE)*V2MED(JE)/2.D+0)
         QY(JE)=N(JE)*(VYV2MED(JE)/2.D+0-VXMED(JE)*VXVYMED(JE)
     &                 -VYMED(JE)*VY2MED(JE)-VZMED(JE)*VYVZMED(JE)
     &                 +VZMED(JE)*VMED2-VZMED(JE)*V2MED(JE)/2.D+0)
         QZ(JE)=N(JE)*(VZV2MED(JE)/2.D+0-VXMED(JE)*VXVZMED(JE)
     &                 -VYMED(JE)*VYVZMED(JE)-VZMED(JE)*VZ2MED(JE)
     &                 +VZMED(JE)*VMED2-VZMED(JE)*V2MED(JE)/2.D+0)
      END DO

*     stampo le variabili macroscopiche in funzione della posizione

      WRITE(2,700)
     & '***** GRANDEZZE MACROSCOPICHE *****************************',
     & '************'
      WRITE(2,700)
     & '                                                           ',
     & '            '
      WRITE(2,800) '**** X *****','**** N *****','** VXMED ***',
     &             '** VYMED ***','** VZMED ***'
      DO JN=1,NNOD
         XE=XMIN+DX*(JN -.5D+0)
         WRITE(2,900) XE,N(JN),VXMED (JN),VYMED(JN),VZMED(JN)
      END DO
      WRITE(2,700)
     & '                                                           ',
     & '            '
      WRITE(2,800) '**** X *****','** FLUX X **','** FLUX Y **',
     &             '** FLUX Z **','*** PXX ****'
      DO JN=1,NNOD
          XE=XMIN+DX*(JN +.5D+0)
          WRITE(2,900) XE,FX(JN),FY(JN),FZ(JN),PXX(JN)
      END DO
      WRITE(2,700)
     & '                                                           ',
     & '            '
      WRITE(2,800) '**** X *****','*** PYY ****','*** PZZ ****',
      &            '*** PXY ****','*** PXZ ****'
      DO JN=1,NNOD
         XE=XMIN+DX*(JN +.5D+0)
         WRITE(2,900) XE,PYY(JN),PZZ(JN),PXY(JN),PXZ(JN)
      END DO
      WRITE(2,700)
     & '                                                           ',
     & '            '
      WRITE(2,800) '**** X *****','*** PYZ ****','**** P *****',
      &            '** TMP TR **','* TMP ROT **'
      DO JN=1,NNOD
         XE=XMIN+DX*(JN +.5D+0)
         WRITE(2,900) XE,PYZ(JN),P(JN),TMPTR(JN),TMPROT(JN)
      END DO
      WRITE(2,700)
     & '                                                           ',
     & '            '
      WRITE(2,800) '**** X *****','*** TMP ****','**** QX ****',
      &            '**** QY ****','**** QZ ****'
      DO JN=1,NNOD
         XE=XMIN+DX*(JN +.5D+0)
         WRITE(2,900) XE,TMP(JN),QX(JN),QY(JN),QZ(JN)
      END DO

      WRITE(*,*) 'ok'

      STOP

100   FORMAT(A20)
200   FORMAT(A1,"Simulazione ... ",F5.1,"% ")
300   FORMAT(A18)
400   FORMAT(A58,I14)
500   FORMAT(A58,E14.4)
600   FORMAT(A58)
700   FORMAT(A60,A12)
800   FORMAT(5(A12,3 X))
900   FORMAT(5(E12.4,3 X))

      END

***** FINE PROGRAM COUETTE *********************************************
