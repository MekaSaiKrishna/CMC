! DIANYUN ZHANG
! DIANYUN.ZHANG@UCONN.EDU
! LAST DEBUDDED: 07/17/2018

!
! USER MATERIAL SUBROUTINE TO MODEL THE FIBER TOW RESPONSE 
! THE FIBER TOW IS MODELED BASED ON COMPOSITE CYLINDER MODEL AND GENERALIZED 
! SELF-CONSISTENT MDTHOD (EXTENDED TO THE NONLINEAR REGIME BY SENCANT METHOD)
! FAILURE CRITERIA FOR ORTHOTROPIC COMPOSITE ARE INCLUDED
! 11 -- FIBER BREAKAGE
! 22/33 -- MATRIX CRACKING
! 12/13/23 -- MATRIX SHEAR CRACKING



! FIBER:  TRANSVERSELY ISOTROPIC (E1f,E2f,NU12f,G12f,G23f)
! MATRIX: ISOTROPIC DAMADING SOLID USING DEFORMATION THEORY (Em, NUm, YIELD, K1, K2)


! ORDERING SCHEME [EPS11 EPS22 EPS33 GAM12 GAM13 GAM23]
!                 [SIG11 SIG22 SIG33 SIG12 SIG13 SIG23]

! REMAINING ISSUES: 1) PRECISION (SINGLE/DOUBLE)
!                   2) IMPROVE THE EFFICIENCY      

      SUBROUTINE CCM3D_ANISO(NSTATV, NPROPS, NOEL, stepTIME, totalTIME, 
     1 DT,L,PROPS, EPS, SIG, STATEV,DDSDDE, isIMPLICIT,DEPS)

C USER INPUT PROPERTIES
C PROPS(1)  = E1f
C PROPS(2)  = E2f
C PROPS(3)  = NU12f
C PROPS(4)  = G12f
C PROPS(5)  = G23f
C PROPS(6)  = Em
C PROPS(7)  = NUm
C PROPS(8)  = YIELD (YIELD STRESS OF THE MATRIX)
C PROPS(9)  = K1 (HARDENING PARAMETER FOR MATRIX)
C PROPS(10) = K2 (HARDENING PARAMETER FOR MATRIX)
C PROPS(11) = MATRIX CRITICAL STRENGTH
C PROPS(12) = ADDITIONAL HARDENING
C PROPS(13) = FIBER VOLUME FRACTION
C PROPS(14) = FIBER TENSILE FAILURE STRESS
C PROPS(15) = FIBER COMPRESSIVE FAILURE STRESS
C PROPS(16) = 90 DEGREEE TENSILE FAILURE STRESS
C PROPS(17) = 90 DEGREEE COMPRESSIVE FAILURE STRESS
C PROPS(18) = IN-PLANE SHEAR FAILURE STRESS
C PROPS(19) = FIBER BREAKAGE FRACTURE TOUGHNESS
C PROPS(20) = FRACTURE TOUGHNESS FOR FIBER COMPRESSIVE FAILURE
C PROPS(21) = MATRIX NORMAL CRACKING FRACTURE TOUGHNESS
C PROPS(22) = MATRIX SHEAR CRACKING FRACTURE TOUGHNESS
C PROPS(23) = DAMPING
C PROPS(24) = E1f COMPRESSION



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      IMPLICIT NONE
      REAL*8::SIG(6),EPS(6),DEPS(6)
      REAL*8::DDSDDE(6,6)
      REAL*8::STATEV
      REAL*8::PROPS
      REAL*8::StepTime,TotalTime,DT,L
      INTEGER::NSTATV,NPROPS,NOEL,isIMPLICIT
    
      
      REAL*8::E1f,E2f,NU12f,G12f,G23f
      REAL*8::Em,NUm,YIELD,K1,K2,Vf
      REAL*8::isYIELD,isCRACKED,isKINK
      REAL*8::EPSEQ_M,SIGEQ_M
      REAL*8::EmS_OLD,NUmS_OLD,EPSEQ_M_OLD
      REAL*8::EmS,NUmS
      REAL*8::E1c,E2c,NU12c,G12c,G23c
      REAL*8::Q(6,6)
      REAL*8::G1C_F,G2C_F,G1C_M,G2C_M
      REAL*8::DAMP
      REAL*8::SIGcr11T,SIGcr11C,SIGcr22T,SIGcr22C
      REAL*8::SIGcr33T,SIGcr33C
      REAL*8::TAUcr12,TAUcr13,TAUcr23
      REAL*8::EPScr11T,EPScr11C,EPScr22T,EPScr22C
      REAL*8::EPScr33T,EPScr33C
      REAL*8::GAMMAcr12,GAMMAcr13,GAMMAcr23      
      REAL*8::HARD,EPSCR_M
      REAL*8::E1c_LIM, E2c_LIM, NU12c_LIM, G12c_LIM, G23c_LIM
      REAL*8::Q_LIM(6,6)
      
      INTEGER::MODE      
      DIMENSION STATEV(NSTATV),PROPS(NPROPS)
      
      
C READ IN PROPERTIES
      E1f = PROPS(1)
      IF (EPS(1) .LT. 0.0D0) E1f = PROPS(24)
      E2f = PROPS(2)
      NU12f = PROPS(3)
      G12f = PROPS(4)
      G23f = PROPS(5)
      Em = PROPS(6)
      NUm = PROPS(7)
      YIELD = PROPS(8)
      K1 = PROPS(9)
      K2 = PROPS(10)
      EPSCR_M = PROPS(11)
      HARD = PROPS(12)
      Vf = PROPS(13)
      SIGcr11T = PROPS(14)
      SIGcr11C = PROPS(15)
      SIGcr22T = PROPS(16)
      SIGcr22C = PROPS(17)
      TAUcr12 = PROPS(18)
      G1C_F = PROPS(19)
      G2C_F = PROPS(20)
      G1C_M = PROPS(21)
      G2C_M = PROPS(22)
      DAMP = PROPS(23)*L
      

      isYIELD = STATEV(1)
      isCRACKED = STATEV(10)
   
      IF ((isYIELD .LT. 0.5) .AND. (isCRACKED .LT. 0.5)) THEN
          CALL MODULI_CALC_GSCM(E1f,E2f,NU12f,G12f,G23f,Em,NUm,Vf,
     1     E1c,E2c,NU12c,G12c,G23c)  
             
C CHECK MATRIX DAMAGE MODE      
          CALL MATRIX_STRAIN(E1f,E2f,NU12f,G12f,G23f,Em,NUm,Vf,
     1       EPS,1,EPSEQ_M)  ! FIND MATRIX STRAINS
               
          IF (EPSEQ_M .GE. (YIELD/Em)) THEN
              isYIELD = 1.0D0
              STATEV(1) = 1.0D0
              STATEV(2)  = EPSEQ_M
              STATEV(3) = Em
          ELSE    
              CALL STIFFNESS_MATRIX(E1c,E2c,NU12c,G12c,G23c,Q)
              SIG = MATMUL(Q,EPS)
              STATEV(3) = Em
              IF(isIMPLICIT .EQ.1) DDSDDE = Q
          END IF
          
          STATEV(4) = E1c
          STATEV(5) = E2c
          STATEV(6) = NU12c
          STATEV(7) = G12c
          STATEV(8) = G23c
          
         IF (SIG(1) .GT. SIGcr11T) THEN
             MODE = 1
             isCRACKED = 1.0D0
         ELSE IF (SIG(1) .LE. SIGcr11C) THEN
             MODE = 2
             isCRACKED = 2.0D0
         ELSE IF (SIG(2) .GE. SIGcr22T) THEN
             MODE = 3
             isCRACKED = 3.0D0
         ELSE IF (SIG(2) .LE. SIGcr22C) THEN
             MODE = 4
             isCRACKED = 4.0D0
         ELSE IF (SIG(3) .GE. SIGcr22T) THEN
             MODE = 5
             isCRACKED = 5.0D0
         ELSE IF (SIG(3) .LE. SIGcr22C) THEN
             MODE = 6
             isCRACKED = 6.0D0
         ELSE IF (DABS(SIG(4)) .GE. TAUcr12) THEN
             MODE = 7
             isCRACKED = 7.0D0
         ELSE IF (DABS(SIG(5)) .GE. TAUcr12) THEN
             MODE = 8
             isCRACKED = 8.0D0
         END IF
         
         IF (isCRACKED .GT. 0.5) THEN
            STATEV(10) = isCRACKED
            STATEV(11) = DABS(SIG(1))
            STATEV(12) = DABS(SIG(2))
            STATEV(13) = DABS(SIG(3))
            STATEV(14) = DABS(SIG(4))
            STATEV(15) = DABS(SIG(5))
            STATEV(16) = DABS(SIG(6))
            STATEV(18:20)=1.0D-8  
         END IF        
      END IF
    
            
      IF ((isYIELD .GE. 0.5) .AND. (isCRACKED .LT. 0.5)) THEN
          
          EPSEQ_M_OLD = STATEV(2)
          EmS_OLD = STATEV(3)
          NUmS_OLD = dmin1(0.5+EmS_OLD/Em*(NUm-0.5),4.900d-1)
          CALL MATRIX_STRAIN(E1f,E2f,NU12f,G12f,G23f,EmS_OLD,
     1       NUmS_OLD,Vf,EPS,1,EPSEQ_M)
          !statev(9)= EPSEQ_M-EPSEQ_M_OLD
          
          
          IF (EPSEQ_M-EPSEQ_M_OLD.lt.-1.0d-3) THEN
              CALL MODULI_CALC_GSCM(E1f,E2f,NU12f,G12f,G23f,Em,NUm,
     1          Vf,E1c,E2c,NU12c,G12c,G23c)
          
              CALL STIFFNESS_MATRIX(E1c,E2c,NU12c,G12c,G23c,Q)
              SIG = SIG+MATMUL(Q,DEPS)            
              EmS = EmS_OLD
              NUmS = NUmS_OLD

             ! CALL MODULI_CALC_GSCM(E1f,E2f,NU12f,G12f,G23f,EmS_OLD,     
     1       !     NUmS_OLD,Vf,E1c,E2c,NU12c,G12c,G23c)
          
              !CALL STIFFNESS_MATRIX(E1c,E2c,NU12c,G12c,G23c,Q)              
          ELSE
            CALL SIGEQ_CALC(YIELD,K1,K2,Em,EPSCR_M,HARD,EPSEQ_M,SIGEQ_M)
            EmS = DMIN1(SIGEQ_M/EPSEQ_M,EmS_OLD) ! THE MODULOUS KEEP DECREASING
            NUmS =dmin1( 0.5+EmS/Em*(NUm-0.5) , 4.900d-1)
          
            CALL MODULI_CALC_GSCM(E1f,E2f,NU12f,G12f,G23f,EmS,NUmS,
     1       Vf,E1c,E2c,NU12c,G12c,G23c)
          
            CALL STIFFNESS_MATRIX(E1c,E2c,NU12c,G12c,G23c,Q)
            SIG = MATMUL(Q,EPS)
          END IF
          
          STATEV(2) = DMAX1(EPSEQ_M,EPSEQ_M_OLD)
          STATEV(3) = EmS
          STATEV(4) = E1c
          STATEV(5) = E2c
          STATEV(6) = NU12c
          STATEV(7) = G12c
          STATEV(8) = G23c
          statev(9) = NUmS

          IF(isIMPLICIT .EQ.1) DDSDDE = Q
 
C    CHECK IF THERE EXISTS MATRIX CRACKING IN THE FIBER TOW             
C          IF ( (EPSEQ_M .GT. 0.05) .OR. 
C     1         (SIG(1) .GT. SIGcr11T(1)) .OR.
C     2         (SIG(1) .LT. SIGcr11C(1)) ) THEN  

         IF (SIG(1) .GT. SIGcr11T) THEN
             MODE = 1
             isCRACKED = 1.0D0
         ELSE IF (SIG(1) .LE. SIGcr11C) THEN
             MODE = 2
             isCRACKED = 2.0D0
         ELSE IF (SIG(2) .GE. SIGcr22T) THEN
             MODE = 3
             isCRACKED = 3.0D0
         ELSE IF (SIG(2) .LE. SIGcr22C) THEN
             MODE = 4
             isCRACKED = 4.0D0
         ELSE IF (SIG(3) .GE. SIGcr22T) THEN
             MODE = 5
             isCRACKED = 5.0D0
         ELSE IF (SIG(3) .LE. SIGcr22C) THEN
             MODE = 6
             isCRACKED = 6.0D0
         ELSE IF (DABS(SIG(4)) .GE. TAUcr12) THEN
             MODE = 7
             isCRACKED = 7.0D0
         ELSE IF (DABS(SIG(5)) .GE. TAUcr12) THEN
             MODE = 8
             isCRACKED = 8.0D0
         END IF
         
         IF (isCRACKED .GT. 0.5) THEN
            STATEV(10) = isCRACKED
            STATEV(11) = DABS(SIG(1))
            STATEV(12) = DABS(SIG(2))
            STATEV(13) = DABS(SIG(3))
            STATEV(14) = DABS(SIG(4))
            STATEV(15) = DABS(SIG(5))
            STATEV(16) = DABS(SIG(6))
            STATEV(18:20)=1.0D-8  
         END IF           
      END IF   
          
      IF (isCRACKED .GT. 0.5) THEN
         CALL SCA3Daniso(NSTATV, NPROPS, NOEL, stepTIME, totalTIME, 
     1 DT, L, PROPS, EPS, SIG, STATEV, DDSDDE, isIMPLICIT)   
         
      END IF  
      
      RETURN 
      END
      

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC          

     


