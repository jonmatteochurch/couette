************************************************************************
*                                                                      *
      SUBROUTINE COLLISIONI(JE)
*                                                                      *
************************************************************************

***** VARIABILI ********************************************************

***** INPUT
*     indice dell elemento della griglia all interno del quale
*     simulare le collisioni
      INTEGER JE

***** PARAMETRI
*     numero massimo di particelle
*     (deve coincidere col valore in couette.f ed inizio.f)
      INTEGER,PARAMETER::NPMAX=200000
*     numero massimo di nodi
*     (deve coincidere col valore in couette.f)
      INTEGER,PARAMETER::NNMAX=1000

***** COSTANTI MATEMATICHE
*     pi greco e suoi multipli
      REAL*8,PARAMETER::PI=3.1415926535897932D+0
      REAL*8,PARAMETER::PI2=6.2831853071795864D+0

***** GRANDEZZE MICROSCOPICHE
*     posizione delle particelle
      REAL*8 X(NPMAX)
*     componenti della velocita delle particelle
      REAL*8 VX(NPMAX),VY(NPMAX),VZ(NPMAX)
*     energia interna delle particelle
      REAL*8 EI(NPMAX)
*     numero totale di particelle
      INTEGER NPAR
      COMMON/MICRO/ X,VX,VY,VZ,EI,NPAR

***** DOMINIO SPAZIALE
*     estremi del dominio
      REAL*8 XMIN,XMAX
*     passi spaziali
      REAL*8 DX
*     numero totale di nodi
      INTEGER NNOD
      COMMON/SPAZIO/ XMIN,XMAX,DX,NNOD

***** DOMINIO TEMPORALE
*     passo temporale
      REAL*8 DT
*     istante finale
      REAL*8 TMAX
*     istante d inizio campionamento grandezze macroscopiche
      REAL*8 TIM
*     passo campionamento
      REAL*8 DTM
      COMMON/TEMPO/ TMAX,DT,TIM,DTM

***** PARAMETRI COLLISIONALI
*     sezione d urto
      REAL*8 XSECT
*     frazione di urti anaelastici
      REAL*8 LAMBDA
      COMMON/COLLIS/ XSECT,LAMBDA

***** INDICI PARTICELLE/ELEMENTI
*     indici delle particelle ordinate per elemento di appartenenza
      INTEGER IP(NPMAX)
*     indici delle prime particelle di ciascun elemento
      INTEGER IPIN(NNMAX)
*     numero di particelle per elemento della griglia
      INTEGER NPE(NNMAX)
      COMMON/IND/ IP,IPIN,NPE

***** VARIABILI LOCALI
*     componenti velocita massima
      REAL*8 VXMAX,VYMAX,VZMAX
*     componenti velocita minima
      REAL*8 VXMIN,VYMIN,VZMIN
*     velocita relativa delle particelle che collidono e suo quadrato
      REAL*8 VR,VR2
*     componenti velocita relativa dopo urto
      REAL*8 VRX,VRY,VRZ
*     modulo velocita relativa massima
      REAL*8 VRMAX
*     energia totale particelle che collidono
      REAL*8 ETOT
*     energia cinetica ed interna particelle dopo l urto
      REAL*8 ECIN,EINT
*     numero atteso di collisioni
      REAL*8 AVNCOLL
*     numero di collisioni
      INTEGER NCOLL
*     contatore per le collisioni
      INTEGER JCOLL
*     tempo trascorso fra l istante precedente e l ultima collisione
      REAL*8 DTCOLL
*     indici delle particelle che collidono
      INTEGER JP1,JP2
*     indice della particella ordinata per elemento
      INTEGER JP
*     indici prima e ultima particella dell elemento
      INTEGER JPEMIN,JPEMAX
*     indice della particella nell elemento
      INTEGER JPE
*     parametri trig method
      REAL*8 B ,C,PHI,SITETA
*     numero casuale con distribuzione uniforme
      REAL*8 R

***** FINE VARIABILI ***************************************************

***** calcolo il numero atteso di collisioni (reali+fittizie)

      JPEMIN=IPIN(JE)+1
      JPEMAX=IPIN(JE)+NPE(JE)

*     calcolo componenti velocita massima e minimadelle particelle
      VXMAX=-1.E+38
      VXMIN=1.E+38
      VYMAX=-1.E+38
      VYMIN=1.E+38
      VZMAX=-1.E+38
      VZMIN=1.E+38
      DO JPE=JPEMIN,JPEMAX
         JP=IP(JPE)
         VXMIN=DMIN1(VXMIN,VX(JP))
         VXMAX=DMAX1(VXMAX,VX(JP))
         VYMIN=DMIN1(VYMIN,VY(JP))
         VYMAX=DMAX1(VYMAX,VY(JP))
         VZMIN=DMIN1(VZMIN,VZ(JP))
         VZMAX=DMAX1(VZMAX,VZ(JP))
      END DO

*     calcolo velocita relativa massima
      VRMAX=DSQRT((VXMAX-VXMIN)**2+(VYMAX-VYMIN)**2+(VZMAX-VZMIN)**2)

*     processo di poisson per determinare il numero di collisioni
      AVNCOLL=NPE(JE)*NPE(JE)*XSECT*VRMAX/(2.D+0*DX)
      NCOLL=0
      DTCOLL=0.D+0
10    R=1.D+0-RAND()
      DTCOLL=DTCOLL-DLOG(R)/ AVNCOLL
      IF(DTCOLL.LT.DT) THEN
        NCOLL=NCOLL +1
        GO TO 10
      END IF

***** calcolo le collisioni fra le particelle

      DO JCOLL=1,NCOLL
*        seleziono due particelle a caso dell elemento
         JP1=JPEMIN+INT(NPE(JE)*RAND())
         JP1=IP(JP1)
         JP2=JPEMIN+INT(NPE(JE)*RAND())
         JP2=IP(JP2)
         calcolo velocita relativa
         VR2=(VX(JP2)-VX(JP1))**2+(VY(JP2)-VY(JP1))**2+
     &       (VZ(JP2)-VZ(JP1))**2
         VR=DSQRT(VR2)

         R=RAND()
         IF(R. LT.VR/VRMAX) THEN
*           collisione reale
            R=RAND()
            IF(R.LT.LAMBDA) THEN

*****          collisione anelastiche - modello di Larsen/Borgnakke

               ETOT=VR2/4.D0+EI(JP1)+EI(JP2)
*              genero energia cinetica dopo l urto
               R=1.D0-RAND()
               ECIN=ETOT*(DCOS(DASIN(2.0D+0*R-1.0D+0)/3.0D+0-PI/2.0D+0)+
     &                    0.5D+0)
*              spartisco energia interna fra le due particelle
               EINT=ETOT-ECIN
               EI(JP1)=EINT*RAND()
               EI(JP2)=EINT-EI(JP1)
               VR=2.0D+0*DSQRT(ECIN)
               B=1.D0-2.D0*RAND()
               SITETA=DSQRT(1.D0-B*B)
               PHI=PI2*RAND()
               VRZ=VR*B
               VRX=VR*SITETA*DCOS(PHI)
               VRY=VR*SITETA*DSIN(PHI)
               C=VX(JP1)+VX(JP2)
               VX(JP1)=(C-VRX)*0.5D0
               VX(JP2)=(C+VRX)*0.5D0
               C=VY(JP1)+VY(JP2)
               VY(JP1)=(C-VRY)*0.5D0
               VY(JP2)=(C+VRY)*0.5D0
               C=VZ(JP1)+VZ(JP2)
               VZ(JP1)=(C-VRZ)*0.5D0
               VZ(JP2)=(C+VRZ)*0.5D0
 
            ELSE

*****          collisioni elastiche

               B=1.D0 -2.D0*RAND()
               SITETA=DSQRT(1.D0-B*B)
               PHI=PI2*RAND()
               VRZ=VR*B
               VRX=VR*SITETA*DCOS(PHI)
               VRY=VR*SITETA*DSIN(PHI)
               C=VX(JP1)+VX(JP2)
               VX(JP1)=(C-VRX)*0.5D0
               VX(JP2)=(C+VRX)*0.5D0
               C=VY(JP1)+VY(JP2)
               VY(JP1)=(C-VRY)*0.5D0
               VY(JP2)=(C+VRY)*0.5D0
               C=VZ(JP1)+VZ(JP2)
               VZ(JP1)=(C-VRZ)*0.5D0
               VZ(JP2)=(C+VRZ)*0.5D0

            END IF
         END IF
      END DO

      RETURN
      END

***** FINE SUBROUTINE COLLIS *******************************************
