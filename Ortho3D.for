      ! Orthotropic Elastic Main File

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      real*8::eps(6)
      
C     Why is eps =STRAN + DSTRAN? why not STRAN alone?      
      eps=STRAN + DSTRAN
      
      call  orthoelas(NSTATV, NPROPS,NOEL, 
     *     TIME(1), TIME(2), DTIME, CELENT, 
     *     PROPS, eps, STRESS, STATEV, DDSDDE, 1)


      RETURN
      END
      
      include 'ortho_elastic.for'
      include 'matrixInverse.for'
      include 'myStandardSupport.for'