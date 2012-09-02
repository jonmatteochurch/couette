************************************************************************
*                                                                      *
      SUBROUTINE VELPAR(UW,TMPW,VX,VY,VZ,EI)
*                                                                      *
*                                                                      *
************************************************************************

***** VARIABILI ********************************************************

***** INPUT
*     velocita parete lungo z
      REAL*8 UW
*     temperatura parete
      REAL*8 TMPW

***** OUTPUT
*     componenti della velocita dopo l urto
      REAL*8 VX,VY,VZ
*     energia interna
      REAL*8 EI

***** VARIABILI LOCALI
*     due volte pi greco
      REAL*8,PARAMETER::PI2=6.2831853071795864D+0
*     radice quadrata della temperatura
      REAL*8 SIGMA
*     numero casuale da distribuzione uniforme
      REAL*8 R
*     coordinate polari
      REAL*8 RHO,THETA

***** FINE VARIABILI ***************************************************

      SIGMA=DSQRT(TMPW)

*     genero modulo velocita nel piano yz
      R=1.D0-RAND()
      RHO=DSQRT(-2.D0*DLOG(R))
*     genero coordinata angolare
      R=RAND()
      THETA=PI2*R
*     calcolo coordinate yz della velocita
      VY=RHO*DCOS(THETA)*SIGMA
      VZ=RHO*DSIN(THETA)*SIGMA+UW

*     genero velocita lungo x
10    R=1.D0-RAND()
      IF (R.EQ .1.0 D0) GOTO 10
      VX=DSQRT(-2.D0*TMPW*DLOG(R))

*     genero valore energia interna
      R=1.D0-RAND()
      EI=-DLOG(R)*TMPW

      RETURN

      END

***** FINE SUBROUTINE VELPAR *******************************************
