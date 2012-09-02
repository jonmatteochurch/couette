
************************************************************************
*                                                                      *
      SUBROUTINE MAXWELL(TMP,VX,VY,VZ,EI)
*                                                                      *
*     Subroutine per la generazione di velocita ed energie secondo     *
*     una distribuzione di Maxwell tramite l argoritmo di Box-Mueller  *
*                                                                      *
************************************************************************

***** VARIABILI ********************************************************

***** INPUT
*     temperatura
      REAL*8 TMP

***** OUTPUT
*     componenti della velocita
      REAL*8 VX,VY,VZ
*     energia interna
      REAL*8 EI

***** VARIABILI LOCALI
*     due volte pi greco
      REAL*8,PARAMETER::PI2=6.2831853071795864
*     radice quadrata della temperatura
      REAL*8 SIGMA
*     numero casuale da distribuzione uniforme
      REAL*8 R
*     coordinate polari
      REAL*8 RHO,THETA

***** FINE VARIABILI ***************************************************

      SIGMA=DSQRT(TMP)

*     genero modulo velocita nel piano yz
      R=1.D0-RAND()
      RHO=DSQRT(-2.D0*DLOG(R))
*     genero coordinata angolare
      R=RAND()
      THETA=PI2*R
*     calcolo coordinate yz della velocita
      VY=RHO*DCOS(THETA)*SIGMA
      VZ=RHO*DSIN(THETA)*SIGMA

*     genero modulo per x
      R=1.D0-RAND()
      RHO=DSQRT(-2.D0*DLOG(R))
*     genero angolo per x
      R=RAND()
      THETA=PI2*R
      VX=RHO*DCOS(THETA)*SIGMA

*     genero valore energia interna
      R=1.D0-RAND()
      EI=-DLOG(R)*TMP

      RETURN

      END

***** FINE SUBROUTINEMAXWELL *******************************************
