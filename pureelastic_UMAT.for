      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*8 CMNAME
C
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2 PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DROT(3, 3),
     3 DFGRD0(3, 3), DFGRD1(3, 3),TIME(2)
 
      DIMENSION EELAS(6), EPLAS(6), FLOW(6)
C     User defined variable
      DIMENSION MSTRESS(6,1),VOLUME(1)
      DIMENSION OLD_DSTRAIN(6), OLD_DSTRESS(6), OLD_STRESS(NTENS)
      DIMENSION STRAN_total(NTENS)      
C
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0,
     1 ENUMAX=.4999D0, NEWTON=10, TOLER=1.0D-6)
C
C ----------------------------------------------------------------
C UMAT FOR ISOTROPIC ELASTICITY
C CANNOT BE USED FOR PLANE STRESS
C ----------------------------------------------------------------
C PROPS(1) - E
C PROPS(2) - NU
C ----------------------------------------------------------------
C ELASTIC PROPERTIES
       EMOD=PROPS(1)
       ENU=PROPS(2)
       EBULK3=EMOD/(ONE-TWO*ENU)
       EG2=EMOD/(ONE+ENU)
       EG=EG2/TWO
       EG3=THREE*EG
       ELAM=(EBULK3-EG2)/THREE
C	   
C     ELASTIC STIFFNESS
C
      DO K1=1, NDI
       DO K2=1, NDI
          DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
      END DO
      DO K1=NDI+1, NTENS
        DDSDDE(K1 ,K1)=EG
      END DO
C
      STRAN_total = STRAN+DSTRAN
      STRESS = MATMUL(DDSDDE,STRAN_total)
      STATEV(1:6) = STRAN_total
      write(*,*) 'STRAN_total'
      write(*,*) STRAN_total
      write(*,*) "        "
      RETURN
      END
