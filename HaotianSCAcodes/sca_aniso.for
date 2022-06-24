! THIS IS THE NEW ORTHOTROPIC MODEL!!!!!

!DIANYUN@UMICH.EDU

! EXTENDED SMEARED CRACK APPROACH

! ordering scheme [eps11 eps22 eps33 gam12 gam13 gam23]
!                 [sig11 sig22 sig33 sig12 sig13 sig23]



! remaining issues:
! - include characteristic length dependent on angle, not on Abaqus
!*********************************************************************!
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

! STATE VARIABLE LIST
C STATEV(11:16) = MAXIMUM CRACK STRAIN 11 22 33 12 13 23 (ENGINEERING)
C STATEV(17:22) = OLD CRACK STRAIN 11 22 33 12 13 23 (ENGINEERING)  

!*********************************************************************!
      SUBROUTINE SCA3Daniso(NSTATV, NPROPS, NOEL, stepTIME, totalTIME, 
     1 DT, L, PROPS, EPS, SIG, STATEV, DDSDDE, isIMPLICIT)
          
      IMPLICIT NONE
      
      REAL*8::SIG(6),EPS(6)
      REAL*8::DDSDDE(6,6)
      REAL*8::Dco(6,6),Dcr(3,3),Dcocr(6,6)
      REAL*8::STATEV,PROPS
      REAL*8::stepTIME,totalTIME,DT,L,DAMP
      REAL*8::Dda(3,3)
      INTEGER::NSTATV,NPROPS,NOEL,isIMPLICIT
     
      REAL*8::SIGcr11T,SIGcr11C,SIGcr22T,SIGcr22C
      REAL*8::SIGcr11, SIGcr22, SIGcr33, TAUcr12, TAUcr13, TAUcr23
c      REAL*8::EPScr11T,EPScr11C,EPScr22T,EPScr22C
c      REAL*8::GAMMAcr12,GAMMAcr13,GAMMAcr23 
      REAL*8::G1C_F,G2C_F,G1C_M,G2C_M
      REAL*8::isCRACKED
      INTEGER::MODE
      REAL*8::A(3,3)   !MATERIAL ORIENTATION / CRACK ORIENTATION
      REAL*8::N1(6,3)
      REAL*8::EPSCR(3),EPSCR_OLD(3),EPSCR_MAX(3)
      REAL*8::TERM63(6,3),TERM33(3,3),invTERM33(3,3)
      REAL*8::TERM66(6,6),invTERM66(6,6)
      REAL*8::E1c,E2c,NU12c,G12c,G23c
      REAL*8::EPSEQ_M
      
      
      DIMENSION STATEV(NSTATV),PROPS(NPROPS)
      

      MODE = INT(STATEV(10))
           
      G1C_F = PROPS(19)
      G2C_F = PROPS(20)
      G1C_M = PROPS(21)
      G2C_M = PROPS(22)
      DAMP = PROPS(23)*L
      
      E1c = STATEV(4)
      E2c = STATEV(5)
      NU12c = STATEV(6)
      G12c = STATEV(7)
      G23c = STATEV(8)
      
      CALL STIFFNESS_MATRIX(E1c,E2c,NU12c,G12c,G23c,Dco)
     
 
      
      EPSCR_MAX = STATEV(18:20)    
      EPSCR_OLD = STATEV(21:23)
      
         IF( (MODE .EQ. 1) .OR. (MODE .EQ. 2) )  THEN
            SIGcr11 = STATEV(11) 
            TAUcr12 = STATEV(14)
            TAUcr13 = STATEV(15)
         END IF

         IF( (MODE .EQ. 3) .OR. (MODE .EQ. 4) .OR. (MODE .EQ. 7)) THEN
            SIGcr22 = STATEV(12)
            TAUcr12 = STATEV(14)
            TAUcr23 = STATEV(16)
         END IF
       
         IF(( MODE .EQ. 5) .OR. (MODE .EQ. 6) .OR. (MODE .EQ. 8)) THEN
            SIGcr22 = STATEV(13)
            TAUcr12 = STATEV(15)
            TAUcr23 = STATEV(16)
         END IF
      
         Dda = 0.0D0
         IF (MODE .EQ. 1) THEN 
             DAMP = DAMP*400.D0
         ELSE IF(MODE .EQ. 2) THEN
             DAMP = DAMP*4.0D0
         END IF
      
      ! TOW CRACK ORIENTATION IS ALIGNED WITH THE MATERIAL COORDINATE

      A = 0.0D0
      A(1,1) = 1.0D0
      A(2,2) = 1.0D0
      A(3,3) = 1.0D0

      CALL getN_UD(A,N1,MODE)
      
      CALL QcalcEPScr_lamina(EPSCR,EPSCR_MAX,EPSCR_OLD,EPS,Dco,
     1   Dda,Dcr,SIGcr11,SIGcr22,SIGcr33,TAUcr12,TAUcr13,
     2   TAUcr23,G1C_F,G2C_F,G1C_M,G2C_M,L,DT,MODE,N1)
     
      
      TERM33 = Dcr+MATMUL(TRANSPOSE(N1),MATMUL(Dco,N1))+Dda/DT
        CALL matrixInverse(TERM33,invTERM33,3,3)
        TERM63 = MATMUL(Dco,MATMUL(N1,invTERM33))
        TERM66 = MATMUL(TERM63,MATMUL(TRANSPOSE(N1),Dco))        
        Dcocr = Dco-TERM66
        SIG = MATMUL(Dco,EPS-MATMUL(N1,EPSCR))
c        IF ( (SIG(1) .GT. 0.0) .AND. (SIG(1) .LT.500.)) THEN
c            STATEV(43) = 0.0    ! ELEMENT DELETE 
c        ENDIF
        
      STATEV(18) = DMAX1(STATEV(18),DABS(EPSCR(1)))
      STATEV(19) = DMAX1(STATEV(19),DABS(EPSCR(2)))
      STATEV(20) = DMAX1(STATEV(20),DABS(EPSCR(3)))
                    
      STATEV(21:23) = EPSCR      
      IF (isIMPLICIT .EQ. 1) DDSDDE = Dcocr
        

      
      RETURN
      END
     
     

