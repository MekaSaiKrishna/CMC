C Dianyun Zhang     
C dianyun.zhang@uconn.edu
C Last Debugged: 07/17/2018



C EXTENDED SMEARED CRACK APPROACH FOR LAMINAR FAILURE

C LAMINA IS LINEAR ELASTIC UP TO FAILURE
C SEPARATE FAILURE MODE

C ordering scheme [eps11 eps22 eps33 gam12 gam13 gam23]
C                 [sig11 sig22 sig33 sig12 sig13 sig23]



C remaining issues:
C     - include characteristic length dependent on angle, not on Abaqus
C     - use deformation gradient for large deformation
C     - mix-mode failure
      
C*********************************************************************!
C USER INPUT PROPERTIES
C PROPS(1)  = E1T (LAMINAR 0 DEGREE TENSION)
C PROPS(2)  = E1C (LMAINAR 0 DEGREE COMPRESSION)
C PROPS(3) = E2  (LAMINAR 90 DEGREE MODULUS)      
C PROPS(4)  = NU12 (LAMINAR IN-PLANE POISSON'S RATIO)
C PROPS(5)  = G12  (LAMINAR IN-PLANE SHEAR)
C PROPS(6)  = G23  (LAMINAR TRANSVERSE SHEAR)
C PROPS(7) = FIBER TENSILE FAILURE STRESS
C PROPS(8) = FIBER COMPRESSIVE FAILURE STRESS
C PROPS(9) = 90 DEGREEE TENSILE FAILURE STRESS
C PROPS(10) = 90 DEGREEE COMPRESSIVE FAILURE STRESS
C PROPS(11) = IN-PLANE SHEAR FAILURE STRESS (1-2)
C PROPS(12) = TRANSVERSE SHEAR FAILURE STRESS (2-3)
C PROPS(13) = FIBER BREAKAGE FRACTURE TOUGHNESS
C PROPS(14) = FRACTURE TOUGHNESS FOR FIBER COMPRESSIVE FAILURE
C PROPS(15) = MATRIX NORMAL CRACKING FRACTURE TOUGHNESS
C PROPS(16) = MATRIX SHEAR CRACKING FRACTURE TOUGHNESS
C PROPS(17) = DAMPING

      
C STATE VARIABLE LIST
C STATEV(31:33) = Thermal Strain
C STATEV(1) = FAILURE INDEX (1-9)

   
C STATEV(11:16) = MAXIMUM CRACK STRAIN 11 22 33 12 13 23 (ENGINEERING)
C STATEV(17:22) = OLD CRACK STRAIN 11 22 33 12 13 23 (ENGINEERING)  

!*********************************************************************!
      SUBROUTINE SCA_LAMINA(NSTATV, NPROPS, NOEL, stepTIME, totalTIME, 
     1 DT, L, PROPS, EPS, SIG, STATEV, DDSDDE, DTEMP, isIMPLICIT)
          
      IMPLICIT NONE
      
      REAL*8::SIG(6),EPS(6)
      REAL*8::DDSDDE(6,6)
      REAL*8::Dco(6,6),Dcr(3,3),Dcocr(6,6)
      REAL*8::STATEV,PROPS
      REAL*8::stepTIME,totalTIME,DT,L,DAMP,DTEMP
      REAL*8::Dda(3,3)
      INTEGER::NSTATV,NPROPS,NOEL,isIMPLICIT
     
      REAL*8::SIGcr11T,SIGcr11C,SIGcr22T,SIGcr22C,SIGcr33T
      REAL*8::SIGcr11, SIGcr22, SIGcr33, TAUcr12, TAUcr13, TAUcr23
c      REAL*8::EPScr11T,EPScr11C,EPScr22T,EPScr22C
c      REAL*8::GAMMAcr12,GAMMAcr13,GAMMAcr23 
      REAL*8::G1C_F,G2C_F,G1C_M,G2C_M
      REAL*8::isCRACKED
      INTEGER::MODE
      REAL*8::A(3,3)   !MATERIAL ORIENTATION / CRACK ORIENTATION
      REAL*8::B(3,3)
      REAL*8::N1(6,3)
      REAL*8::EPSCR(3),EPSCR_OLD(3),EPSCR_MAX(3)
      REAL*8::TERM63(6,3),TERM33(3,3),invTERM33(3,3)
      REAL*8::TERM66(6,6),invTERM66(6,6)
      REAL*8::E1,E2,NU12,G12,G23
      REAL*8::FI_FT,FI_FC,FII_12,FII_13
      REAL*8::temp_SIG(6)
      REAL*8::prc_temp_SIG(3),prc_temp_angle(3,3)
      REAL*8::matrixs_shear_ROT(3,3),matrixs_crack_ori(3,3)
      REAL*8::brocken
      REAL*8::alpha1,alpha2,alpha3
      INTEGER::max_index,min_index

      real*8::sigCr(3)


      DIMENSION STATEV(NSTATV),PROPS(NPROPS)

C READ IN PROPERTIES
      E1 = PROPS(1)
      IF (EPS(1) .LT. 0.0D0) E1 = PROPS(2)
      E2 = PROPS(3)
      NU12 = PROPS(4)
      G12 = PROPS(5)
      G23 = PROPS(6)
      SIGcr11T = PROPS(7)
      SIGcr11C = PROPS(8)
      SIGcr22T = PROPS(9)
      SIGcr33T = PROPS(10)
      TAUcr12 = PROPS(11)
      TAUcr23 = PROPS(12)
      G1C_F = PROPS(13)
      G2C_F = PROPS(14)
      G1C_M = PROPS(15)
      G2C_M = PROPS(16)
      DAMP = PROPS(17)*L
      alpha1 = PROPS(18)
      alpha2 = PROPS(19)
      alpha3 = PROPS(20)
      
      brocken = 0.0D0
      isCRACKED = STATEV(10)

      STATEV(31) = STATEV(31)+alpha1*DTEMP
      STATEV(32) = STATEV(32)+alpha2*DTEMP
      STATEV(33) = STATEV(33)+alpha3*DTEMP
      STATEV(34) = STATEV(32)

      EPS(1) = EPS(1)-STATEV(31)
      EPS(2) = EPS(2)-STATEV(32)
      EPS(3) = EPS(3)-STATEV(33)

C     BEFORE FAILURE -- LINEAR ELASTIC RESPONSE
      IF (isCRACKED .LT. 0.5) THEN

         CALL STIFFNESS_MATRIX(E1,E2,NU12,G12,G23,Dco)

         SIG = MATMUL(Dco,EPS)
         IF(isIMPLICIT .EQ.1) DDSDDE = Dco

         FI_FT = (SIG(1)/SIGcr11T)**2
         FI_FC = (SIG(1)/SIGcr11C)**2
         FII_12 = SIG(4)**2/TAUcr12**2
         FII_13 = SIG(5)**2/TAUcr12**2


         IF ((SIG(1) .GT. 0.0D0) .AND. (FI_FT .GE. 1.0D0)) THEN 
            MODE = 1
            isCRACKED = 1.0D0
         END IF
 
         IF ((SIG(1) .LT. 0.0D0) .AND. (FI_FC .GE. 1.0D0)) THEN 
            MODE = 2
            isCRACKED = 2.0D0
         END IF

         IF (FII_12 .GE. 1.0D0) THEN
             MODE = 7
             isCRACKED = 7.0D0
         END IF

         IF (FII_13 .GE. 1.0D0) THEN            
             MODE = 8
             isCRACKED = 8.0D0
         END IF

         IF (((MODE .NE. 1) .or. (MODE .NE. 2) .or. (MODE .NE. 7) 
     1     .or. (MODE .NE. 8)) .and. (SIG(2)+SIG(3) .GT. 0.0D0))THEN 
            temp_SIG = SIG
            temp_SIG(1) = 0.0D0
            temp_SIG(4) = 0.0D0
            temp_SIG(5) = 0.0D0

            call sprind(temp_SIG,prc_temp_SIG,prc_temp_angle,1,3,3)
            STATEV(13) = DMAX1(prc_temp_SIG(1),
     1         prc_temp_SIG(2),prc_temp_SIG(3))
            STATEV(14) = DMIN1(prc_temp_SIG(1),
     1         prc_temp_SIG(2),prc_temp_SIG(3))
            STATEV(15) = sqrt((temp_SIG(2)-temp_SIG(3))**2.0D0+
     1         4.0D0*temp_SIG(6)**2.0D0)/2

            IF (SIG(2) .GT. SIG(3)) THEN
              IF (STATEV(13) .GT. SIGcr22T) THEN
                MODE = 3
                isCRACKED = 3.0D0
              ELSE IF (STATEV(15) .GT. TAUcr23) THEN
                MODE = 9
                isCRACKED = 9.0D0
              END IF
            ELSE
              IF (STATEV(13) .GT. SIGcr33T) THEN
                MODE = 4
                isCRACKED = 4.0D0
              ELSE IF (STATEV(15) .GT. TAUcr23) THEN
                MODE = 9
                isCRACKED = 9.0D0
              END IF
            END IF
         END IF
          
         IF (((MODE .NE. 1) .or. (MODE .NE. 2) .or. (MODE .NE. 7) 
     1     .or. (MODE .NE. 8)) .and. (SIG(2)+SIG(3) .LT. 0.0D0))THEN  
            temp_SIG = SIG
            temp_SIG(1) = 0.0D0
            temp_SIG(4) = 0.0D0
            temp_SIG(5) = 0.0D0

            call sprind(temp_SIG,prc_temp_SIG,prc_temp_angle,1,3,3)
            STATEV(13) = DMAX1(prc_temp_SIG(1),
     1         prc_temp_SIG(2),prc_temp_SIG(3))
            STATEV(14) = DMIN1(prc_temp_SIG(1),
     1         prc_temp_SIG(2),prc_temp_SIG(3))
            STATEV(15) = sqrt((temp_SIG(2)-temp_SIG(3))**2.0D0+
     1         4.0D0*temp_SIG(6)**2.0D0)/2

C             IF (STATEV(14) .LT. -1.0D0*SIGcr22C) THEN
C               MODE = 4
C               isCRACKED = 4.0D0
C             ELSE IF (STATEV(15) .GT. TAUcr23) THEN
C               MODE = 9
C               isCRACKED = 9.0D0
C             END IF
         END IF

         IF (isCRACKED .GT. 0.5) THEN
            IF ((MODE .EQ. 1) .OR. (MODE .EQ. 2)) THEN

              STATEV(11) = DABS(SIG(1))
              STATEV(12) = DABS(SIG(1))

C               write(*,*) "DABS(SIG(1))"
C               write(*,*) DABS(SIG(1))
              
              STATEV(1:3) = (/1,0,0/)
              STATEV(4:6) = (/0,1,0/)
              STATEV(7:9) = (/0,0,1/)
            END IF

            IF (MODE .EQ. 7) THEN
              STATEV(11) = DABS(SIG(1))
              STATEV(12) = DABS(SIG(4))
              
              STATEV(1:3) = (/0,1,0/)
              STATEV(4:6) = (/1,0,0/)
              STATEV(7:9) = (/0,0,1/)
            END IF

            IF (MODE .EQ. 8) THEN
              STATEV(11) = DABS(SIG(1))
              STATEV(12) = DABS(SIG(5))
              
              STATEV(1:3) = (/0,0,1/)
              STATEV(4:6) = (/1,0,0/)
              STATEV(7:9) = (/0,1,0/)
            END IF

            IF (MODE .EQ. 3) THEN
              temp_SIG = SIG
              temp_SIG(1) = 0.0D0
              temp_SIG(4) = 0.0D0
              temp_SIG(5) = 0.0D0

              call sprind(temp_SIG,prc_temp_SIG,prc_temp_angle,1,3,3)
              STATEV(13) = DMAX1(prc_temp_SIG(1),
     1           prc_temp_SIG(2),prc_temp_SIG(3))
              STATEV(14) = DMIN1(prc_temp_SIG(1),
     1           prc_temp_SIG(2),prc_temp_SIG(3))
              STATEV(15) = sqrt((temp_SIG(2)-temp_SIG(3))**2.0D0+
     1           4.0D0*temp_SIG(6)**2.0D0)/2

              IF (STATEV(13).EQ.prc_temp_SIG(1)) THEN
                max_index = 1
              ELSE IF (STATEV(13).EQ.prc_temp_SIG(2))THEN
                max_index = 2
              ELSE
                max_index = 3
              END IF

              matrixs_crack_ori(:,1) = prc_temp_angle(max_index,:)
              matrixs_crack_ori(:,2) = (/1,0,0/)

              matrixs_crack_ori(1,3) = 
     1        matrixs_crack_ori(2,1)*matrixs_crack_ori(3,2)-
     2        matrixs_crack_ori(3,1)*matrixs_crack_ori(2,2)

              matrixs_crack_ori(2,3) = 
     1        matrixs_crack_ori(1,2)*matrixs_crack_ori(3,1)-
     2        matrixs_crack_ori(1,1)*matrixs_crack_ori(3,2)

              matrixs_crack_ori(3,3) = 
     1        matrixs_crack_ori(1,1)*matrixs_crack_ori(2,2)-
     2        matrixs_crack_ori(2,1)*matrixs_crack_ori(1,2)

              STATEV(1:3) = matrixs_crack_ori(:,1)
              STATEV(4:6) = matrixs_crack_ori(:,2)
              STATEV(7:9) = matrixs_crack_ori(:,3)
            END IF

            IF (MODE .EQ. 4) THEN
              temp_SIG = SIG
              temp_SIG(1) = 0.0D0
              temp_SIG(4) = 0.0D0
              temp_SIG(5) = 0.0D0

              call sprind(temp_SIG,prc_temp_SIG,prc_temp_angle,1,3,3)
              STATEV(13) = DMAX1(prc_temp_SIG(1),
     1           prc_temp_SIG(2),prc_temp_SIG(3))
              STATEV(14) = DMIN1(prc_temp_SIG(1),
     1           prc_temp_SIG(2),prc_temp_SIG(3))
              STATEV(15) = sqrt((temp_SIG(2)-temp_SIG(3))**2.0D0+
     1           4.0D0*temp_SIG(6)**2.0D0)/2

              IF (STATEV(14).EQ.prc_temp_SIG(1)) THEN
                min_index = 1
              ELSE IF (STATEV(13).EQ.prc_temp_SIG(2))THEN
                min_index = 2
              ELSE
                min_index = 3
              END IF

              matrixs_crack_ori(:,1) = prc_temp_angle(min_index,:)
              matrixs_crack_ori(:,2) = (/1,0,0/)

              matrixs_crack_ori(1,3) = 
     1        matrixs_crack_ori(2,1)*matrixs_crack_ori(3,2)-
     2        matrixs_crack_ori(3,1)*matrixs_crack_ori(2,2)

              matrixs_crack_ori(2,3) = 
     1        matrixs_crack_ori(1,2)*matrixs_crack_ori(3,1)-
     2        matrixs_crack_ori(1,1)*matrixs_crack_ori(3,2)

              matrixs_crack_ori(3,3) = 
     1        matrixs_crack_ori(1,1)*matrixs_crack_ori(2,2)-
     2        matrixs_crack_ori(2,1)*matrixs_crack_ori(1,2)

              STATEV(1:3) = matrixs_crack_ori(:,1)
              STATEV(4:6) = matrixs_crack_ori(:,2)
              STATEV(7:9) = matrixs_crack_ori(:,3)
            END IF

            IF (MODE .EQ. 9) THEN

              temp_SIG = SIG
              temp_SIG(1) = 0.0D0
              temp_SIG(4) = 0.0D0
              temp_SIG(5) = 0.0D0

              call sprind(temp_SIG,prc_temp_SIG,prc_temp_angle,1,3,3)
              STATEV(13) = DMAX1(DABS(prc_temp_SIG(1)),
     1           DABS(prc_temp_SIG(2)),DABS(prc_temp_SIG(3)))
              STATEV(14) = DMIN1(DABS(prc_temp_SIG(1)),
     1           DABS(prc_temp_SIG(2)),DABS(prc_temp_SIG(3)))
              IF (STATEV(13).EQ.DABS(prc_temp_SIG(1))) THEN
                max_index = 1
              ELSE IF (STATEV(13).EQ.DABS(prc_temp_SIG(2)))THEN
                max_index = 2
              ELSE
                max_index = 3
              END IF

              IF (STATEV(14).EQ.DABS(prc_temp_SIG(1))) THEN
                min_index = 1
              ELSE IF (STATEV(14).EQ.DABS(prc_temp_SIG(2)))THEN
                min_index = 2
              ELSE
                min_index = 3
              END IF

              matrixs_crack_ori(:,1) = prc_temp_angle(max_index,:)
              matrixs_crack_ori(:,2) = (/1,0,0/)

              matrixs_crack_ori(1,3) = 
     1        matrixs_crack_ori(2,1)*matrixs_crack_ori(3,2)-
     2        matrixs_crack_ori(3,1)*matrixs_crack_ori(2,2)

              matrixs_crack_ori(2,3) = 
     1        matrixs_crack_ori(1,2)*matrixs_crack_ori(3,1)-
     2        matrixs_crack_ori(1,1)*matrixs_crack_ori(3,2)

              matrixs_crack_ori(3,3) = 
     1        matrixs_crack_ori(1,1)*matrixs_crack_ori(2,2)-
     2        matrixs_crack_ori(2,1)*matrixs_crack_ori(1,2)

              matrixs_shear_ROT = 0.0D0

              matrixs_shear_ROT(1,1) =1.0D0
              matrixs_shear_ROT(2,2) = sqrt(2.0D0)/2.0D0
              matrixs_shear_ROT(3,3) = sqrt(2.0D0)/2.0D0
              matrixs_shear_ROT(2,3) = -sqrt(2.0D0)/2.0D0
              matrixs_shear_ROT(3,2) = -matrixs_shear_ROT(2,3)

              matrixs_crack_ori(:,1) = matmul(matrixs_shear_ROT,
     1         matrixs_crack_ori(:,1))  
              matrixs_crack_ori(:,2) = matmul(matrixs_shear_ROT,
     1         matrixs_crack_ori(:,2)) 
              matrixs_crack_ori(:,3) = matmul(matrixs_shear_ROT,
     1         matrixs_crack_ori(:,3)) 

              STATEV(1:3) = matrixs_crack_ori(:,1)
              STATEV(4:6) = matrixs_crack_ori(:,2)
              STATEV(7:9) = matrixs_crack_ori(:,3)
            END IF

            STATEV(10) = isCRACKED
            STATEV(18:20)=1.0D-8  
         END IF          

C     DEFINE THE CRITICAL STRESS
      ELSE IF (isCRACKED .GE. 0.5) THEN

         MODE = INT(STATEV(10))
         EPSCR_MAX = STATEV(18:20)    
         EPSCR_OLD = STATEV(21:23)
         
         CALL STIFFNESS_MATRIX(E1,E2,NU12,G12,G23,Dco)

         IF (MODE .EQ. 1) THEN
            SIGcr11 = STATEV(11) 
         END IF

         IF (MODE .EQ. 2) THEN
            SIGcr11 = STATEV(12)
C             write(*,*) "SIGcr11_top"
C             write(*,*) SIGcr11 
         END IF


         IF((MODE .EQ. 3) .OR. (MODE .EQ. 4)) THEN
            SIGcr22 = STATEV(13)
         END IF
      
         IF((MODE .EQ. 7) .OR. (MODE .EQ. 8))THEN
            TAUcr12 = STATEV(12)
            TAUcr13 = STATEV(12)
         END IF

         IF (MODE .EQ. 9) THEN
            TAUcr23 = STATEV(15) 
         END IF
      
         Dda = 0.0D0
         IF (MODE .EQ. 1) THEN 
             DAMP = DAMP*400.D0
         ELSE IF(MODE .EQ. 2) THEN
             DAMP = DAMP*4.0D0
         END IF
         
         Dda(1,1) = DAMP
         Dda(2,2) = DAMP
         Dda(3,3) = DAMP

      
      ! TOW CRACK ORIENTATION IS ALIGNED WITH THE MATERIAL COORDINATE

         A = 0.0D0
         A(1,:) = STATEV(1:3)
         A(2,:) = STATEV(4:6)
         A(3,:) = STATEV(7:9)

         CALL getN_UD(A,N1,MODE)
      
         CALL QcalcEPScr_lamina(EPSCR,EPSCR_MAX,EPSCR_OLD,EPS,Dco,
     1   Dda,Dcr,SIGcr11,SIGcr22,SIGcr33,TAUcr12,TAUcr13,
     2   TAUcr23,G1C_F,G2C_F,G1C_M,G2C_M,L,DT,MODE,N1,brocken)
      
         TERM33 = Dcr+MATMUL(TRANSPOSE(N1),MATMUL(Dco,N1))+Dda/DT
         CALL matrixInverse(TERM33,invTERM33,3,3)
         TERM63 = MATMUL(Dco,MATMUL(N1,invTERM33))
         TERM66 = MATMUL(TERM63,MATMUL(TRANSPOSE(N1),Dco))        
         Dcocr = Dco-TERM66
         SIG = MATMUL(Dco,EPS-MATMUL(N1,EPSCR))
c        IF ( (SIG(1) .GT. 0.0) .AND. (SIG(1) .LT.500.)) THEN
c            STATEV(43) = 0.0    ! ELEMENT DELETE 
c        ENDIF

C          sigCr = matmul(Dcr,EPSCR)
C          STATEV(31:33)=sigCr
        
         STATEV(18) = DMAX1(STATEV(18),DABS(EPSCR(1)))
         STATEV(19) = DMAX1(STATEV(19),DABS(EPSCR(2)))
         STATEV(20) = DMAX1(STATEV(20),DABS(EPSCR(3)))
                    
         STATEV(21:23) = EPSCR
         STATEV(24) = STATEV(22)
         STATEV(25) = brocken

         STATEV(26) = Dcr(1,1)
         STATEV(27) = Dcr(2,2)
         STATEV(28) = Dcr(3,3)

         STATEV(29) = L
         IF (isIMPLICIT .EQ. 1) DDSDDE = Dcocr
        
      END IF 
     
      RETURN
      END
     
     
!------------------------------------------------------------------ 
      subroutine getN_UD(o,N,MODE)
      IMPLICIT NONE
      
      real*8::o(3,3)
      real*8::N(6,3)
      INTEGER::MODE
      
      N(1,1) = o(1,1) ** 2
      N(1,2) = o(1,1) * o(2,1)
      N(1,3) = o(3,1) * o(1,1)
      N(2,1) = o(1,2) ** 2
      N(2,2) = o(1,2) * o(2,2)
      N(2,3) = o(3,2) * o(1,2)
      N(3,1) = o(1,3) ** 2
      N(3,2) = o(1,3) * o(2,3)
      N(3,3) = o(1,3) * o(3,3)
      N(4,1) = 2.0D0 * o(1,1) * o(1,2)
      N(4,2) = o(1,1) * o(2,2) + o(1,2) * o(2,1)
      N(4,3) = o(3,1) * o(1,2) + o(3,2) * o(1,1)
      N(5,1) = 2.0D0 * o(1,3) * o(1,1)
      N(5,2) = o(1,3) * o(2,1) + o(1,1) * o(2,3)
      N(5,3) = o(3,3) * o(1,1) + o(3,1) * o(1,3)
      N(6,1) = 2.0D0 * o(1,2) * o(1,3)
      N(6,2) = o(1,2) * o(2,3) + o(1,3) * o(2,2)
      N(6,3) = o(3,2) * o(1,3) + o(3,3) * o(1,2)

C       write(*,*) "o:  ", o
C       write(*,*) "N:  ", N
      RETURN
      END      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C SUBROUTINE TO FORM THE STIFFNESS MATRIX FOR TRANSVERSELY ISOTROPIC MATERIAL
!      SUBROUTINE STIFFNESS_MATRIX(E1,E2,NU12,G12,G23,Q)
!      
!      IMPLICIT NONE
!      REAL*8::E1,E2,NU12,G12,G23,NU23
!      REAL*8::DELTA 
!      REAL*8::Q(6,6)
!      
!      NU23 = E2/(2.0*G23)-1.0
!      DELTA = 2.0*E2*NU12**2+E1*(NU23-1.0)
!      Q = 0.0
!      Q(1,1) = E1**2*(NU23-1.0)/DELTA
!      Q(1,2) = -E1*E2*NU12/DELTA
!      Q(1,3) = Q(1,2)
!      Q(2,1) = Q(1,2)
!      Q(2,2) = E2*(E2*NU12**2-E1)/((1.0+NU23)*DELTA)
!      Q(2,3) = -(E2*(E2*NU12**2+E1*NU23))/((1.0+NU23)*DELTA)
!      Q(3,1) = Q(1,3)
!      Q(3,2) = Q(2,3)
!      Q(3,3) = Q(2,2)
!      Q(4,4) = G12
!      Q(5,5) = G12
!      Q(6,6) = G23 
!      
!      RETURN
!      END 
      

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC          
      SUBROUTINE  QcalcEPScr_lamina(EPSCR,EPSCR_MAX,EPSCR_OLD,EPS,Dco,
     1 Dda,Dcr,SIGcr11,SIGcr22,SIGcr33,TAUcr12,TAUcr13,
     2 TAUcr23,G1C_F,G2C_F,G1C_M,G2C_M,L,DT,MODE,N,brocken)    
   
      IMPLICIT NONE
      REAL*8::EPSCR(3),EPSCR_MAX(3),EPSCR_OLD(3),EPS(6)
      REAL*8::Dco(6,6),Dda(3,3),Dcr(3,3)
      REAL*8::SIGcr11,SIGcr22,SIGcr33
      REAL*8::TAUcr12,TAUcr13,TAUcr23
      REAL*8::G1C_F,G2C_F,G1C_M,G2C_M
      REAL*8::L,DT
      REAL*8::N(6,3)
      
      REAL*16::QEPSCR_MAX(3),QEPSCR_OLD(3),QEPS(6)
      REAL*16::QDco(6,6),QDda(3,3),QDcr(3,3)
      REAL*16::QSIGcr11,QSIGcr22,QSIGcr33
      REAL*16::QTAUcr12,QTAUcr13,QTAUcr23
      REAL*16::QG1C_F,QG2C_F,QG1C_M,QG2C_M
      REAL*16::QL,QDT
      
      REAL*16::QN(6,3),QNT(3,6),QNT_DCO_N(3,3)
      
      REAL*16::F_TOL,DELTA
      REAL*16::X(3),DX(3),FO(3),FN(3),XplusDX(3),X0(3)
      REAL*16::J(3,3),invJ(3,3)
      REAL*8::brocken

      INTEGER::IMAX,I,K,flag
      
      INTEGER::MODE
      
      flag = 1
      QEPSCR_MAX = EPSCR_MAX
      QEPSCR_OLD = EPSCR_OLD
      QEPS = EPS
      QDco = Dco
      QDda = Dda
      QSIGcr11 = SIGcr11
      QSIGcr22 = SIGcr22
      QSIGcr33 = SIGcr33
      QTAUcr12 = TAUcr12
      QTAUcr13 = TAUcr13
      QTAUcr23 = TAUcr23
      QG1C_F = G1C_F
      QG2C_F = G2C_F
      QG1C_M = G1C_M
      QG2C_M = G2C_M
      QL = L
      QDT = DT
      QNT = TRANSPOSE(N)
      QNT_DCO_N = matmul(transpose(N),matmul(Dco,N))
     
      IMAX = 20
      F_TOL = 1.0D-4
      X = 0.0D0
      DO I = 1, IMAX
        CALL QcalcDcr_lamina(X,QEPSCR_MAX,QSIGcr11,QSIGcr22,QSIGcr33,
     1     QTAUcr12,QTAUcr13,QTAUcr23, QG1C_F, QG2C_F, QG1C_M,QG2C_M,
     2     QL,QDcr,MODE,brocken)

        FO=matmul(QDcr+QNT_Dco_N+QDda/QDT,X)
     1  -matmul(QNT,matmul(QDco,QEPS))-matmul(QDda/QDT,QEPSCR_OLD)

        IF( QABS(FO(1))+QABS(FO(2))+QABS(FO(3)) .LT. F_TOL )EXIT 
     
C         write(*,*) ' '
C         write(*,*) '*****************'
C         write(*,*) I
        DO K=1,3
          DELTA = QSIGN(1.0Q-16,X(K))
          XplusDX = X
          XplusDX(K)=X(K)+DELTA
          CALL QcalcDcr_lamina(XplusDX,QEPSCR_MAX,QSIGcr11,
     1     QSIGcr22, QSIGcr33,QTAUcr12,QTAUcr13,QTAUcr23,
     2     QG1C_F, QG2C_F, QG1C_M,QG2C_M,QL,QDcr,MODE,brocken)
          
          FN=MATMUL(QDcr+QNT_Dco_N+QDda/QDT,XplusDX)
     1       -MATMUL(QNT,matmul(QDco,QEPS))-MATMUL(QDda/QDT,QEPSCR_OLD)
          
          J(:,K)=(FN-FO)/DELTA    

C           write(*,*) '*****************'
C           write(*,*) 'K:  ',K

C C           write(*,*) '***********'
C C           write(*,*) 'FO:'
C C           write(*,*) FO
  
C C           write(*,*) '***********'
C C           write(*,*) 'FN:'
C C           write(*,*) FN

C           write(*,*) '***********'
C           write(*,*) 'FN-FO:'
C           write(*,*) FN-FO

C           write(*,*) '***********'
C           write(*,*) 'J:'
C           write(*,*) J

        END DO
        DX = FO
        CALL qSolveLinSysLU(J,DX,3,flag)
        IF (flag. EQ. 0) EXIT       
        X = X-DX
      END DO
           
      IF ((I.GE.IMAX).OR.(flag .EQ. 0)) THEN
      !A little bit of diagnostics
c$$$      write(*,*)
c$$$      write(*,*) 'check tow'
C           write(*,*) 'Newtons`s method failed to converge in max number'
          write(*,*) 'I'
          write(*,*) I
          write(*,*) 'flag'
          write(*,*) flag
c$$$      !write(*,*) 'of increments'
c$$$      write(*,*) 'Dco='
c$$$      write(*,*) Dco
c$$$      write(*,*) 'Dcr='
c$$$      write(*,*) DBLE(QDcr(1,1)), DBLE(QDcr(2,2)), DBLE(QDcr(3,3))
c$$$      write(*,*) 'MODE',MODE
c$$$      write(*,*) DBLE(QTAUcr12),dble(QSIGcr11T),dble(QSIGcr11C)
c$$$      !write(*,*) 'eps='
c$$$      !write(*,*) EPS
c$$$      !write(*,*) 'epscrmax='
c$$$      !write(*,*) EPSCR_MAX
c$$$      !write(*,*) 'epscrold'
c$$$      !write(*,*) EPSCR_OLD
c$$$      !write(*,*) 'GICm,GICf,GIICm,sigcr11,sigcr22,taucr12,taucr23'
c$$$      !write(*,*) G1C_F,G1C_M,G2C_M,SIGcr11,SIGcr22,TAUcr12,TAUcr23
c$$$      !write(*,*) 'The final values of the iteration'
c$$$      write(*,*) 'epscr='
c$$$      write(*,*) DBLE(X)
c$$$      write(*,*) 'F'
c$$$      write(*,*) DBLE(FO)

      !call myExit()
      
      endif 
      
      EPSCR = DBLE(X)
      Dcr = DBLE(QDcr)
    
      RETURN 
      END
      

      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      SUBROUTINE QcalcDcr_lamina(EPSCR,EPSCR_MAX,SIGcr11,
     1 SIGcr22,SIGcr33,TAUcr12,TAUcr13,TAUcr23,
     2 G1C_F,G2C_F,G1C_M,G2C_M,L,Dcr,MODE,brocken)
      
      IMPLICIT NONE
      REAL*16::EPSCR(3),EPSCR_MAX(3),Dcr(3,3)
      REAL*16::Dcr11_T,Dcr11_C
      REAL*16::SIGcr11,SIGcr22,SIGcr33
      REAL*16::TAUcr12,TAUcr13,TAUcr23
      REAL*16::G1C_F,G2C_F,G1C_M,G2C_M
      REAL*16::L
      REAL*16::SIGcr11_R,SIGcr22_R
      REAL*16::TAUcr12_R,TAUcr13_R,TAUcr23_R
      REAL*16::DNORM, DSHEAR1,DSHEAR2
      REAL*8::brocken
      INTEGER::MODE
      
      
      DNORM = QMAX1(QABS(EPSCR(1)),QABS(EPSCR_MAX(1)))
      DSHEAR1 = QMAX1(QABS(EPSCR(2)),QABS(EPSCR_MAX(2)))
      DSHEAR2 = QMAX1(QABS(EPSCR(3)),QABS(EPSCR_MAX(3)))
      
      Dcr = 0.0Q0
      
      IF (MODE .EQ. 1) SIGcr11_R = SIGcr11*1.0Q-5
      IF (MODE .EQ. 2) THEN
          SIGcr11_R = SIGcr11*0.5Q0
          G1C_F = G2C_F
      END IF 
      
      IF (MODE .EQ. 3) SIGcr22_R = SIGcr22*1.0Q-8
      IF (MODE .EQ. 4) SIGcr22_R = SIGcr22*1.0Q-8

      TAUcr13_R = TAUcr13*1.0Q-5
      TAUcr12_R = TAUcr12*1.0Q-5
      TAUcr23_R = TAUcr23*1.0Q-5


      IF ((MODE .EQ. 1) .OR. (MODE .EQ. 2)) THEN
          IF (DNORM .LT. (2.0Q0*G1C_F/L/(SIGcr11-SIGcr11_R)) ) THEN
             Dcr(1,1) = (-(SIGcr11_R-SIGcr11)**2/(2.0Q0*G1C_F/L)*
     1       DNORM+SIGcr11)/DNORM
          ELSE
            Dcr(1,1) = SIGcr11_R/DNORM
            brocken = 1.1D-1
          END IF            

C           IF (DNORM .LT. (2.0Q0*G1C_F/L/(SIGcr11-SIGcr11_R)) ) THEN
C              Dcr(2,2) = 1.0Q16
C              Dcr(3,3) = 1.0Q16
C           ELSE
C              Dcr(2,2) = 1.0Q-5
C              Dcr(3,3) = 1.0Q-5
C           END IF   
      END IF


      IF ((MODE .EQ. 3) .OR. (MODE .EQ. 4)) THEN
          IF (DNORM .LT. (2.0Q0*G1C_M/L/(SIGcr22-SIGcr22_R)) ) THEN
             Dcr(1,1) = (-(SIGcr22_R-SIGcr22)**2/(2.0Q0*G1C_M/L)*
     1       DNORM+SIGcr22)/DNORM
          ELSE
             Dcr(1,1) = SIGcr22_R/DNORM
             brocken = 1.1D-1
          END IF
          Dcr(2,2) = Dcr(1,1)
          Dcr(3,3) = Dcr(1,1)
      END IF

      IF ((MODE .EQ. 7) .OR. (MODE .EQ. 8)) THEN
C           WRITE(*,*) "G2C_M:"
C           write(*,*) G2C_M
C           WRITE(*,*) "TAUcr12:"
C           write(*,*) TAUcr12
C           WRITE(*,*) "TAUcr12_R:"
C           write(*,*) TAUcr12_R
C           WRITE(*,*) "DSHEAR1:"
C           write(*,*) DSHEAR1
          IF (DSHEAR1 .LT. (2.0Q0*G2C_M/L/(TAUcr12-TAUcr12_R)) ) THEN
              Dcr(2,2) = (-(TAUcr12_R-TAUcr12)**2/(2.0Q0*G2C_M/L)*
     1        DSHEAR1+TAUcr12)/DSHEAR1
          ELSE
              Dcr(2,2) = TAUcr12_R/DSHEAR1
              brocken = 1.1D-1
          END IF
          Dcr(1,1) = Dcr(2,2)
          Dcr(3,3) = Dcr(2,2)
      END IF

      IF (MODE .EQ. 9) THEN
          IF (DSHEAR2 .LT. (2.0Q0*G2C_M/L/(TAUcr23-TAUcr23_R)) ) THEN
             Dcr(3,3) = (-(TAUcr23_R-TAUcr23)**2/(2.0Q0*G2C_M/L)*
     1       DSHEAR2+TAUcr23)/DSHEAR2
          ELSE
            Dcr(3,3) = TAUcr23_R/DSHEAR2
            brocken = 1.1D-1
          END IF
          Dcr(1,1) = Dcr(3,3)
          Dcr(2,2) = Dcr(3,3)
      END IF

      write(*,*) "SIGcr22"
      write(*,*) SIGcr22  


      RETURN 
      END  
      
      


