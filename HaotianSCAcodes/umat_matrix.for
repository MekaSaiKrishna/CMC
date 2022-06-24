! DIANYUN ZHANG, UNIVERSITY OF MICHIGAN
! dianyun@umich.edu
!
! INTERFACE TO USE J2 DEFORMATION THEORY AND SMEARED CRACK BAND APPROACH WITHIN ABAQUS/STANDARD
! PLASTICITY: ISOTROPIC NONLINEAR HARDENING USING DEFORMATION THEORY
! SMEARED CRACK BAND: EXPONENTIAL SOFTERNING


! ORDERING SCHEME [EPS11 EPS22 EPS33 GAM12 GAM13 GAM23]
!                 [SIG11 SIG22 SIG33 SIG12 SIG13 SIG23]
!                 [PE11 PE22 PE33 PE12 PE13 PE 23]

      SUBROUTINE UMAT_MATRIX(NSTATV, NPROPS, NOEL, StepTime, TotalTime,
     1 DT, L, PROPS, EPS, SIG, STATEV, DDSDDE, isIMPLICIT)
      



C USER INPUT PROPERTIES
C PROPS(1) = ELASTIC MODULUS
C PROPS(2) = POISSON'S RATIO
C PROPS(3) = YIELD (YIELD STRESS)
C PROPS(4) = K1 (HARDENING PARAMETER)
C PROPS(5) = K2 (HARDENING PARAMETER)
C PROPS(6) = MODE I CRITICAL STRESS
C PROPS(7) = MODE I FRACTURE TOUGHNESS
C PROPS(8) = DAMPING


C UPDATED STATE VARIABLES
C STATEV(1) = YIELD FLAG
C STATEV(2) = EQUIVALENT STRAIN
C STATEV(3) = EQUIVALENT STRESS
C STATEV(4) = SECANT MODULUS
C STATEV(5) = POISSON'S RATIO
C STATEV(6) = J2 MAX IN THE PREVIOUS STEPS
C STATEV(9) = PLASTIC ENERGY DISSIPATION
C STATEV(10) = CRACK FLAG
C STATEV(11) = MAX. CRACK STRAIN 1
C STATEV(12) = MAX. CRACK STRAIN 2
C STATEV(13) = MAX. CRACK STRAIN 3
C STATEV(14:16) = NORMAL TO CRACK STRAIN 1
C STATEV(17:19) = NORMAL TO CRACK STRAIN 2
C STATEV(20:22) = NORMAL TO CRACK STRAIN 3
C STATEV(23:31) = OLD CRACK STRAIN


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
      
      REAL*8::SIG(6),EPS(6)
      REAL*8::DDSDDE(6,6)
      REAL*8::STATEV
      REAL*8::PROPS
      REAL*8::StepTime,TotalTime,DT,L
      INTEGER::NSTATV,NPROPS,NOEL,isIMPLICIT
      
      REAL*8:: E,NU,YIELD,K1,K2,SIGCR0,GIC,DAMP
      REAL*8:: MU,LAMBDA,BULK
      REAL*8:: ES,NUS,NUS_NEW,MUS,LAMBDAS
      REAL*8:: isYIELD,isCRACKED,TOL
      REAL*8:: J2,J2P,J2P_MAX
      REAL*8:: EPSEQ,SIGEQ
   
      REAL*8:: DE(6,6),DS(6,6)
      REAL*8:: S(6)
      
      real*8::princsig(3)
      real*8::T(3,3)
      
      DIMENSION STATEV(NSTATV),PROPS(NPROPS)
      
      
      E = PROPS(1)
      NU = PROPS(2)
      YIELD = PROPS(3)
      K1 = PROPS(4)
      K2 = PROPS(5)
      SIGCR0 = PROPS(6)
      GIC = PROPS(7)
      DAMP = PROPS(8)
  
C CHECK THE ELEMENT SIZE FOR CRACK BAND ELEMENTS
!      IF (CELENT .GT. 2.0*GIC*E/SIGCR0**2) THEN
!        WRITE(*,*) '  '
!        WRITE(*,*) 'ELEMENT SIZE IS TOO LARGE'
!        WRITE(*,*) 'CELENT - (2.0*GIC*E/SIGCRO**2)'
!        WRITE(*,*) CELENT, 2.0*GIC*E/SIGCRO**2
!        CALL XIT
!      END IF

C COMPUTE ELASTIC CONSTANS AND STIFFNESS MATRIX      
      MU = 0.5*E/(1.0+NU)
      LAMBDA = E*NU/((1.0+NU)*(1.0-2.0*NU))
      BULK = E/(3.0*(1.0-2.0*NU))
      
      DE = 0.0
      DE(1,1) = LAMBDA+2.0*MU
      DE(1,2) = LAMBDA
      DE(1,3) = LAMBDA
      DE(2,1) = LAMBDA
      DE(2,2) = LAMBDA+2.0*MU
      DE(2,3) = LAMBDA
      DE(3,1) = LAMBDA
      DE(3,2) = LAMBDA
      DE(3,3) = LAMBDA+2.0*MU
      DE(4,4) = MU
      DE(5,5) = MU
      DE(6,6) = MU
      
      isYIELD = STATEV(1)
      isCRACKED = STATEV(10)      
      
      IF ((isYIELD .LT. 0.5) .AND. (isCRACKED .LT. 0.5)) THEN
C LINEAR ELASTIC UNTIL YIELDING
        DDSDDE = DE
        SIG = MATMUL(DE,EPS)     
      
C CHECK WHETHER YIELDING OCCURS (J2 THEORY)
        
        CALL StressDeviator(SIG,S)
        
        CALL vonMises(S,J2)
        
        J2P = 1.0/6.0*((EPS(1)-EPS(2))**2+
     1  (EPS(2)-EPS(3))**2+(EPS(3)-EPS(1))**2)+
     2   1.0/4.0*(EPS(4)**2+EPS(5)**2+EPS(6)**2)
      
        IF (SQRT(3.0*J2) .GE. YIELD) THEN
          isYIELD = 1.0D0
          STATEV(1) = 1.0D0
          STATEV(4) = E
          STATEV(5) = NU
          STATEV(6) = J2P
        END IF        

        
C     CHECK IF THE CRACK EXISTS
        CALL mysprind(SIG,princsig,T,1,3,3)
      
        IF (dmax1(princsig(1),princsig(2),princsig(3)).GT.SIGCR0) THEN
        isCRACKED     =1.0D0
        STATEV(10)    =1.0D0
        STATEV(11)    =1.0D-8
        STATEV(12)    =1.0D-8
        STATEV(13)    =1.0D-8
        STATEV(14:16)  =T(1,:)
        STATEV(17:19) =T(2,:)
        STATEV(20:22)=T(3,:)
        END IF


      END IF

      IF ((isYIELD .GT. 0.5) .AND.(isCRACKED .LT. 0.5)) THEN
              
        ES = STATEV(4)
        NUS = STATEV(5)
        J2P_MAX = STATEV(6)
        

        J2P = 1.0/6.0*((EPS(1)-EPS(2))**2+
     1  (EPS(2)-EPS(3))**2+(EPS(3)-EPS(1))**2)+
     2   1.0/4.0*(EPS(4)**2+EPS(5)**2+EPS(6)**2)
     
C CONDITION FOR LOADING     
        IF (J2P .GT. J2P_MAX) THEN  
            
            J2P_MAX = J2P
            
C COMPUTE THE EQUIVALENT STRAIN
            TOL = 1.0D-4
9           CONTINUE
            EPSEQ = 1.0/(1.0+NUS)*SQRT(0.5*((EPS(1)-EPS(2))**2+
     1       (EPS(2)-EPS(3))**2+(EPS(3)-EPS(1))**2)+
     2        3.0/4.0*(EPS(4)**2+EPS(5)**2+EPS(6)**2))        
            
C COMPUTE THE EQUIVALENT STRESS
            SIGEQ = YIELD - K1/K2*(EXP(-K2*EPSEQ)-EXP(-K2*(YIELD/E)))
C COMPUTE THE NEW SECANT MODULUS        
            ES = SIGEQ/EPSEQ
C UPDATE THE NEW SECANT POISSON'S RATIO
            NUS_NEW = 0.5+ES/E*(NU-0.5)       
            
            IF(ABS(NUS_NEW-NUS) .GT. TOL) THEN
                NUS = NUS_NEW
                GOTO 9
            END IF
            
            NUS = NUS_NEW
            STATEV(2) = EPSEQ
            STATEV(3) = SIGEQ
            
            
C CALCULATE UPDATED SECANT STIFFNESS        
            MUS = 0.5*ES/(1.0+NUS)
            LAMBDAS = ES*NUS/((1.0+NUS)*(1.0-2.0*NUS))
            DS = 0.0
            DS(1,1) = LAMBDAS+2.0*MUS
            DS(1,2) = LAMBDAS
            DS(1,3) = LAMBDAS
            DS(2,1) = LAMBDAS
            DS(2,2) = LAMBDAS+2.0*MUS
            DS(2,3) = LAMBDAS
            DS(3,1) = LAMBDAS
            DS(3,2) = LAMBDAS
            DS(3,3) = LAMBDAS+2.0*MUS
            DS(4,4) = MUS
            DS(5,5) = MUS
            DS(6,6) = MUS
            
            DDSDDE = DS
            SIG = MATMUL(DS,EPS)
            
            ELSE
C UNDER THE YIELDING SURFACE, NO UPDATE ON THE SECANT STIFFNESS
                MUS = 0.5*ES/(1.0+NUS)
                LAMBDAS = ES*NUS/((1.0+NUS)*(1.0-2.0*NUS))
                DS = 0.0
                DS(1,1) = LAMBDAS+2.0*MUS
                DS(1,2) = LAMBDAS
                DS(1,3) = LAMBDAS
                DS(2,1) = LAMBDAS
                DS(2,2) = LAMBDAS+2.0*MUS
                DS(2,3) = LAMBDAS
                DS(3,1) = LAMBDAS
                DS(3,2) = LAMBDAS
                DS(3,3) = LAMBDAS+2.0*MUS
                DS(4,4) = MUS
                DS(5,5) = MUS
                DS(6,6) = MUS
                
                DDSDDE = DS
                SIG = MATMUL(DS,EPS)
            END IF
        
            
            STATEV(4) = ES
            STATEV(5) = NUS
            STATEV(6) = J2P_MAX

        
C CHECK IF THE CRACK EXISTS
        CALL mysprind(SIG,princsig,T,1,3,3)
      
        IF (dmax1(princsig(1),princsig(2),princsig(3)).GT.SIGCR0) THEN
        isCRACKED     =1.0D0
        STATEV(10)    =1.0D0
        STATEV(11)    =1.0D-8
        STATEV(12)    =1.0D-8
        STATEV(13)    =1.0D-8
        STATEV(14:16)  =T(1,:)
        STATEV(17:19) =T(2,:)
        STATEV(20:22)=T(3,:)
        END IF
        
      END IF  
      
      IF (isCRACKED .GT. 0.5) THEN     
          CALL sca3diso(NSTATV, NPROPS,NOEL,StepTime,TotalTime,
     1         DT, L,PROPS, EPS, SIG, STATEV, DDSDDE,isIMPLICIT)
      END IF
       
      RETURN
      END
      
      
    

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE StressDeviator(STRESS,S)
      IMPLICIT NONE
      REAL*8:: PRESS
      REAL*8:: STRESS(6), S(6)
      
      PRESS = (STRESS(1)+STRESS(2)+STRESS(3))/3.0
      S(1)=STRESS(1)-PRESS
      S(2)=STRESS(2)-PRESS
      S(3)=STRESS(3)-PRESS
      S(4)=STRESS(4)
      S(5)=STRESS(5)
      S(6)=STRESS(6)
      
      RETURN
      END

     
      SUBROUTINE StrainDeviator(STRAN,EPS_MEAN,EPSD)
      IMPLICIT NONE
      REAL*8:: EPS_MEAN
      REAL*8:: STRAN(6),EPSD(6)
      
      EPS_MEAN = (STRAN(1)+STRAN(2)+STRAN(3))/3.0
      EPSD(1) = STRAN(1)-EPS_MEAN
      EPSD(2) = STRAN(2)-EPS_MEAN
      EPSD(3) = STRAN(3)-EPS_MEAN
      EPSD(4) = STRAN(4)/2.0
      EPSD(5) = STRAN(5)/2.0
      EPSD(6) = STRAN(6)/2.0   
      
      RETURN      
      END
      
      
      SUBROUTINE vonMises(S,J2) 
      IMPLICIT NONE
      REAL*8:: J2
      REAL*8:: S(6)     
      
      J2 = 0.5*(S(1)**2+S(2)**2+S(3)**2+2.0*S(4)**2+2.0*S(5)**2
     1          +2.0*S(6)**2)
      
      RETURN
      END
      
      
      
