! DIANYUN ZHANG, UNIVERSITY OF MICHIGAN
! dianyun@umich.edu
! LAST DEBUDDED: 12/03/2012
!
! SUBROUTINES FOR CCM3D
! REMAINING TASK:
!     1) BETTER ALGORITHM TO COMPUTE MAXIMUM AND SORTING
!     2) MORE EFFICIENT WAY TO SOLVE THE LINEAR EQUATIONS

! BETA MODEL -- APPROXIMATION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                            C 
C                                 SUBROUTINES                                C
C                                                                            C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C SUBROUTINE TO CALCUALTE THE COMPOSITE CONSTANTS BASED ON COMPOSITE CYLINDER
C MODEL AND SELF-CONSISTENT METHOD
C E11,NU12,G12 -> CCM
C G23, E22 -> GSCM

      SUBROUTINE MODULI_CALC_GSCM(E1f,E2f,NU12f,G12f,G23f,Em,NUm,Vf,
     1   E1c,E2c,NU12c,G12c,G23c)
     
      IMPLICIT NONE
      REAL*8::Em,NUm,Vf
      REAL*8::E1f,E2f,NU12f,G12f,G23f
      REAL*8::E1c,E2c,NU12c,G12c,G23c,K23c
      REAL*8::NU21f,NU23f,K23f,Gm,K23m
      REAL*8::TERM,GAMMA,DELTA
      REAL*8::TERM1,TERM2,PSI      
      
      
      NU21f = E2f*NU12f/E1f
      NU23f = E2f/(2.0*G23f)-1.0
      K23f = E2f/(2.0*(1.0-2.0*NU12f*NU21f-NU23f))
      
      Gm = Em/(2.0*(1.0+NUm))
      K23m = Em/(2.0*(1.0-NUm-2.0*NUm**2))

C LONGITUDINAL YOUNG'S MODULUS E1C   
      TERM = E2f*(1.0+NUm)*(1.0+Vf*(1.0-2.0*NUm))+Em*
     1   (1.0-NU23f-2.0*NU12f*NU21f)*(1.0-Vf)
      GAMMA = (2.0*Em*NU21f*(1.0-Vf)*(NU12f-NUm))/TERM
      DELTA = (2.0*E2f*NUm*Vf*(NUm-NU12f))/TERM
      E1c = E1f*(1.0+GAMMA)*Vf+Em*(1.0+DELTA)*(1.0-Vf)
      
C MAJOR POISSON'S RATIO NU12C
      TERM1 = (1.0-Vf)*(1.0-NU23f-2.0*NU12f*NU21f)*NUm*Em+(NUm+Vf*
     1  (2.0*NU12f-NUm)+(NUm**2*(1.0-2.0*Vf*NU12f-Vf)))*E2f
      TERM2 = (1.0-Vf)*(1.0-NU23f-2.0*NU12f*NU21f)*Em+(1.0+Vf+(1.0-Vf)*
     1   NUm-2.0*Vf*NUm**2)*E2f
      NU12c = TERM1/TERM2
      
C IN-PLANE SHEAR MODULUS G12C
      TERM1 = G12f*(1.0+Vf)+Gm*(1.0-Vf);
      TERM2 = G12f*(1.0-Vf)+Gm*(1.0+Vf);
      G12c = Gm*(TERM1/TERM2);

C CONSTANTS G23 AND K23(PLANE-STRAIN BULK MODULUS)
      CALL CALC23_GSCM(K23f,G23f,K23m,Gm,Vf,K23c,G23c)

C TRANSVERSE YOUNG'S MODULUS E2c
      PSI = 1.0+(4.0*K23c*NU12c**2)/E1c
      E2c = (4.0*G23c*K23c)/(K23c+PSI*G23c);       
      
      
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
C SUBROUTINE TO COMPUTE 2-3 PLANE PROPERTIES: G23 & K23
      SUBROUTINE CALC23_GSCM(K23f,G23f,K23m,Gm,Vf,K23c,G23c)
      IMPLICIT NONE 
      REAL*8::K23f,G23f,K23m,Gm,Vf,K23c,G23c
      REAL*8::A,B,C      
      
C OUT-OF-PLANE SHEAR MODULUS G23
C SOLVE THE QUADRATIC EQUATION: A*G23^2+B*G23+C=0
      A = 2.0*(-(2.0*G23f*Gm +(G23f+Gm)*K23f)*(2.0*Gm+K23m)*(2.0*G23f*Gm
     1 +(G23f+Gm)*K23m) + 4.0*(G23f-Gm)*(2.0*G23f*Gm +(G23f+Gm)*K23f)*
     2 (Gm**2 + Gm*K23m + K23m**2)*Vf - 6.0*(G23f - Gm)*(2.0*G23f*Gm + 
     3 (G23f+Gm)*K23f)*K23m**2*Vf**2+4.0*(G23f**2*Gm**2*K23f+G23f**2*Gm*
     4 (-Gm + K23f)*K23m +(G23f*(G23f-2.0*Gm)*Gm +(G23f-Gm)*(G23f+Gm)*
     5  K23f)*K23m**2)*Vf**3 +(G23f-Gm)*(2.0*Gm + K23m)*
     6 (Gm*K23f*K23m - G23f*(2.0*Gm*(K23f - K23m) + K23f*K23m))*Vf**4)
      
      B = 4.0*Gm*(Gm*(2.0*G23f*Gm + (G23f + Gm)*K23f)*(2.0*G23f*Gm + 
     1 (G23f + Gm)*K23m)+2.0*(G23f-Gm)*(2.0*G23f*Gm +(G23f + Gm)*K23f)*
     2 (Gm - K23m)*K23m*Vf + 6.0*(G23f - Gm)*(2.0*G23f*Gm + (G23f + Gm)*
     3  K23f)*K23m**2*Vf**2 - 4.0*(G23f**2*Gm**2*K23f + G23f**2*Gm*
     4(-Gm + K23f)*K23m +(G23f*(G23f-2.0*Gm)*Gm +(G23f -Gm)*(G23f+Gm)*
     5 K23f)*K23m**2)*Vf**3 + Gm*(-G23f + Gm)*(Gm*K23f*K23m -
     6 G23f*(2.0*Gm*(K23f - K23m)+ K23f*K23m))*Vf**4)
      
      C = 2.0*Gm**2*(Gm**2*K23f*K23m**2*(-1.0+Vf)**4+2.0*G23f*Gm*K23m*
     1 (-K23f*K23m*(-1.0+Vf**4)+Gm*(K23f + K23m*(-1.0 + Vf)**4 - 
     2 K23f*Vf**4))+G23f**2*(4.0*Gm**2*(K23m +K23f*Vf**3-K23m*Vf**3)
     3 + 2.0*Gm*K23m*(K23f + K23m + 4.0*K23m*Vf - 6.0*K23m*Vf**2 + 
     4 2.0*(K23f + K23m)*Vf**3 +(K23f-K23m)*Vf**4) + K23f*K23m**2*
     5 (1.0 + Vf*(4.0 + Vf*(-6.0 + Vf*(4.0 + Vf))))))

      G23c = (-B-DSQRT(B**2-4.0*A*C))/(2.0*A);      
      
C PLANE-STRAIN BULK MODULUS K23
      K23c = K23m+Vf/(1.0d0/(K23f-K23m)+(1.0d0-Vf)/(K23m+Gm))
     
      RETURN
      END
     
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C SUBROUTINE TO FORM THE STIFFNESS MATRIX FOR TRANSVERSELY ISOTROPIC MATERIAL
      SUBROUTINE STIFFNESS_MATRIX(E1,E2,NU12,G12,G23,Q)
      
      IMPLICIT NONE
      REAL*8::E1,E2,NU12,G12,G23,NU23
      REAL*8::DELTA 
      REAL*8::Q(6,6)
      
      NU23 = E2/(2.0*G23)-1.0
      DELTA = 2.0*E2*NU12**2+E1*(NU23-1.0)
      Q = 0.0
      Q(1,1) = E1**2*(NU23-1.0)/DELTA
      Q(1,2) = -E1*E2*NU12/DELTA
      Q(1,3) = Q(1,2)
      Q(2,1) = Q(1,2)
      Q(2,2) = E2*(E2*NU12**2-E1)/((1.0+NU23)*DELTA)
      Q(2,3) = -(E2*(E2*NU12**2+E1*NU23))/((1.0+NU23)*DELTA)
      Q(3,1) = Q(1,3)
      Q(3,2) = Q(2,3)
      Q(3,3) = Q(2,2)
      Q(4,4) = G12
      Q(5,5) = G12
      Q(6,6) = G23 
      
      RETURN
      END 
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C SUBROUTINES TO COMPUTE EQUIVALENT MATRIX STRAIN
      SUBROUTINE MATRIX_STRAIN(E1f,E2f,NU12f,G12f,G23f,Em,NUm,Vf,
     1     EPS,ISV,EPSEQ)   
C ISV:1,3 (MODE 2 CALCULATED INDIVIDUALLY)      
      IMPLICIT NONE
      INTEGER::ISV
      REAL*8::E1f,E2f,NU12f,G12f,G23f
      REAL*8::Em,NUm,Vf
      REAL*8::NU21f,NU23f,K23f,K23m,Gm
      REAL*8::EPS(6)
      REAL*8::G23c,K23c
      REAL*8::M2,N2,A2,B2,C2,D2
      REAL*8::P11,P41,P42,P21,P22,P23,P24,P25
      REAL*8::G11,G22,G33,G66,G12,G13,G16,G23
      REAL*8::G26,G36
      REAL*8::AVG45,VNORMAL,VSHEAR,VMAX
      REAL*8::ANGLE,ANGLE3,V3,RATIO
      REAL*8::EPSEQ
      REAL*8::PI
      PARAMETER(PI=3.1415927D0)  
      
        
      NU21f = E2f*NU12f/E1f
      NU23f = E2f/(2.0*G23f)-1.0
      K23f = E2f/(2.0*(1.0-2.0*NU12f*NU21f-NU23f))
      
      Gm = Em/(2.0*(1.0+NUm))
      K23m = Em/(2.0*(1.0-NUm-2.0*NUm**2))

      
C COMPUTE F22, F32, F23, F33, F62, F63, F26, F36, F66 BASED UPON 2-3 PLANE STRAIN PROBLEM
      CALL CALC23_GSCM(K23f,G23f,K23m,Gm,Vf,K23c,G23c) 
            
      CALL COEFF_23(K23f,G23f,K23m,Gm,Vf,M2,N2,A2,B2,C2,D2)
C DEFINE SOME CONSTATNS
      P11 = (Vf*(-NUm/K23f+NU12f/K23m))/
     1   (Vf/K23f+(1.0-2.0*NUm)/K23f+(1.0-Vf)/K23m)
      
      P41 = Vf*(G12f-Gm)/(G12f+Gm-Vf*(G12f-Gm))
      P42 = (G12f+Gm)/(G12f+Gm-Vf*(G12f-Gm))
      
      P21 = 0.5*K23c/K23m*N2
      P22 = G23c/(2.0*Gm)*(A2+3.0*B2*Vf)
      P23 = K23c/(4.0*Gm)*M2/Vf
      P24 = G23c/(2.0*K23m)*(3.0*B2*Vf-D2/Vf)
      P25 = G23c/(4.0*Gm)*(3.0*C2/(Vf**2)+2.0*D2/Vf)
      
C COMPUTE GIJ
      G11 = (P11-1.0)**2
      G22 = P21**2+3.0*(P22**2+P23**2-2.0*P22*P25+P25**2)
      G33 = G22
      G66 = 3.0*(P22-P25)**2
      G12 = 2.0*(P11-1.0)*P21
      G13 = G12
      G23 = 2.0*(P21**2-3.0*(P22**2-P23**2-2.0*P22*P25+P25**2))
      AVG45 = 0.75*((P41/Vf)**2+P42**2)*(EPS(4)**2+EPS(5)**2)
      VNORMAL = G11*EPS(1)**2+G22*(EPS(2)**2+EPS(3)**2)+
     1  G12*EPS(1)*(EPS(2)+EPS(3))+G23*EPS(2)*EPS(3)
      VSHEAR = G66*EPS(6)**2
      VMAX = VNORMAL+VSHEAR+AVG45
      !write(*,*)'test',vnormal,P11,NUm,K23m,Em
      

      
C COMPUTE EQUIVALENT MATRIX STRAIN FOR 3 CASES
      SELECT CASE(ISV)
      CASE(1)
        EPSEQ = 1.0/(1.0+NUm)*DSQRT(VMAX)
        
      CASE(3)
        RATIO = VSHEAR/(G11*EPS(1)**2+G22*(EPS(2)**2+EPS(3)**2))
        ANGLE = DATAN(RATIO)
        ANGLE3 = DMIN1(ABS(ANGLE),ABS(ANGLE-0.25*PI),ABS(ANGLE-0.5*PI),
     1      ABS(ANGLE+0.25*PI),ABS(ANGLE+0.5*PI))
    
        V3 = VMAX*DCOS(2.0*ANGLE3)
        EPSEQ = 1.0/(1.0+NUm)*DSQRT(V3)  ! NEED MODIFIED LATER ON
      
      END SELECT
      
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SUBROUTINES TO COMPUTE F-MATRIX THAT RELATES MATRIX STRAINS AND COMPOSITE STRAINS

!      SUBROUTINE CALC_V1(Ef,NUf,Em,NUm,Vf,TH,EPS,V1)
!     
!      IMPLICIT NONE
!      REAL*8::Ef,NUf,Em,NUm,Vf,TH
!      REAL*8::EPS(6)
!      REAL*8::F11,F12,F22,F32,F62,F13,F23,F33,F63
!      REAL*8::F44,F54,F45,F55,F26,F36,F66
!      REAL*8::Gf,Kf,K23f,Gm,Km,K23m
!      REAL*8::K23c,G23c
!      REAL*8::M2,N2,A2,B2,C2,D2
!      REAL*8::P11,P41,P42,P21,P22,P23,P24,P25
!      REAL*8::TERM1,TERM2,TERM3,TERM4,TERM5,TERM6
!      REAL*8::EPS_M_11,EPS_M_22,EPS_M_33
!      REAL*8::GAMMA_M_12,GAMMA_M_13,GAMMA_M_23
!      REAL*8::AVG45,V1
!           
!      Gf = Ef/(2.0*(1.0+NUf))
!      Kf = Ef/(3.0*(1.0-2.0*NUf))   ! BULK MODULUS
!      K23f = Kf + Gf/3.0;           ! PLANE-STRAIN BULK MODULUS
!      
!      Gm = Em/(2.0*(1.0+NUm))
!      Km = Em/(3.0*(1.0-2.0*NUm))  ! BULK MODULUS
!      K23m = Km + Gm/3.0           ! PLANE-STRAIN BULK MODULUS
!
!      
!C COMPUTE F22, F32, F23, F33, F62, F63, F26, F36, F66 BASED UPON 2-3 PLANE STRAIN PROBLEM
!      CALL CALC23_GSCM(Ef,NUf,Em,NUm,Vf,K23c,G23c) 
!            
!      CALL COEFF_23(Ef,NUf,Em,NUm,Vf,M2,N2,A2,B2,C2,D2)
!     
!C DEFINE SOME CONSTATNS
!      P11 = (Vf*(-NUm/K23f+NUf/K23m))/
!     1   (Vf/K23f+(1.0-2.0*NUm)/K23f+(1.0-Vf)/K23m)
!      
!      P41 = Vf*(Gf-Gm)/(Gf+Gm-Vf*(Gf-Gm))
!      P42 = (Gf+Gm)/(Gf+Gm-Vf*(Gf-Gm))
!      
!      P21 = 0.5*K23c/K23m*N2
!      P22 = G23c/(2.0*Gm)*(A2+3.0*B2*Vf)
!      P23 = K23c/(4.0*Gm)*M2/Vf
!      P24 = G23c/(2.0*K23m)*(3.0*B2*Vf-D2/Vf)
!      P25 = G23c/(4.0*Gm)*(3.0*C2/(Vf**2)+2.0*D2/Vf)
!      
!      F11 = 1.0
!      F12 = P11*(1.0-1.0/Vf*DCOS(2.0*TH))
!      F13 = P11*(1.0+1.0/Vf*DCOS(2.0*TH))
!      F44 = P42-P41/Vf*DCOS(2.0*TH)
!      F55 = P42+P41/Vf*DCOS(2.0*TH)
!      F45 = P41/Vf*DSIN(2.0*TH)
!      F54 = F45
!      F22 = (P21-P22)-(P23+P24)*DCOS(2.0*TH)-P25*DCOS(4.0*TH)
!      F23 = (P21+P22)+(P23-P24)*DCOS(2.0*TH)+P25*DCOS(4.0*TH)
!      F26 = 2.0*P23*DSIN(2.0*TH)+2.0*P25*DSIN(4.0*TH)
!      F33 = (P21-P22)+(P23+P24)*DCOS(2.0*TH)-P25*DCOS(4.0*TH)
!      F32 = (P21+P22)-(P23-P24)*DCOS(2.0*TH)+P25*DCOS(4.0*TH)
!      F36 = 2.0*P23*DSIN(2.0*TH)-2.0*P25*DSIN(4.0*TH)
!      F62 = P25*DSIN(4.0*TH)+P24*DSIN(2.0*TH)
!      F63 = -P25*DSIN(4.0*TH)+P24*DSIN(2.0*TH)
!      F66 = -2.0*P22+2.0*P25*DCOS(4.0*TH)
!      
!    
!C COMPUTE MATRIX STRAINS
!      EPS_M_11 = F11*EPS(1)
!      EPS_M_22 = F12*EPS(1)+F22*EPS(2)+F32*EPS(3)+F62*EPS(6)
!      EPS_M_33 = F13*EPS(1)+F23*EPS(2)+F33*EPS(3)+F63*EPS(6)
!C      GAMMA_M_12 = F44*EPS(4)+F54*EPS(5)
!C      GAMMA_M_13 = F45*EPS(4)+F55*EPS(5)
!      GAMMA_M_23 = F26*EPS(2)+F36*EPS(3)+F66*EPS(6)
!      
!      TERM1 = (EPS_M_11 - EPS_M_22)**2
!      TERM2 = (EPS_M_22 - EPS_M_33)**2
!      TERM3 = (EPS_M_33 - EPS_M_11)**2
!C      TERM4 = GAMMA_M_12**2
!C      TERM5 = GAMMA_M_13**2
!      TERM6 = GAMMA_M_23**2     
!      
!      AVG45 = 0.75*((P41/Vf)**2+P42**2)*(EPS(4)**2+EPS(5)**2)
!      V1 = 0.5*(TERM1+TERM2+TERM3)+0.75*TERM6 +AVG45  
!      
!      if(TH .eq. 0.0d0) then
!        !write(*,*) P11,P21,P22,P23,P24,P25
!        
!      
!      end if      
!      
!           
!      RETURN
!      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SUBROUTINES TO CALCULATE V2
      SUBROUTINE MATRIX_STRAIN2(E1f,E2f,NU12f,G12f,G23f,Em,NUm,Vf,
     1     EPS,EPSEQ)
      IMPLICIT NONE
      REAL*8::E1f,E2f,NU12f,G12f,G23f
      REAL*8::Em,NUm,Vf
      REAL*8::EPS(6)
      REAL*8::NU21f,NU23f,K23f,K23m,Gm
      REAL*8::G23c,K23c
      REAL*8::M2,N2,A2,B2,C2,D2
      REAL*8::P11,P41,P42,P21,P22,P23,P24,P25
      REAL*8::G11AVG,G22AVG,G66AVG,G12AVG,G23AVG
      REAL*8::AVG45,V2,EPSEQ
           
      NU21f = E2f*NU12f/E1f
      NU23f = E2f/(2.0*G23f)-1.0
      K23f = E2f/(2.0*(1.0-2.0*NU12f*NU21f-NU23f))
      
      Gm = Em/(2.0*(1.0+NUm))
      K23m = Em/(2.0*(1.0-NUm-2.0*NUm**2))

      
C COMPUTE F22, F32, F23, F33, F62, F63, F26, F36, F66 BASED UPON 2-3 PLANE STRAIN PROBLEM
      CALL CALC23_GSCM(K23f,G23f,K23m,Gm,Vf,K23c,G23c) 
            
      CALL COEFF_23(K23f,G23f,K23m,Gm,Vf,M2,N2,A2,B2,C2,D2)
      
C DEFINE SOME CONSTATNS
      P11 = (Vf*(-NUm/K23f+NU12f/K23m))/
     1   (Vf/K23f+(1.0-2.0*NUm)/K23f+(1.0-Vf)/K23m)
      
      P41 = Vf*(G12f-Gm)/(G12f+Gm-Vf*(G12f-Gm))
      P42 = (G12f+Gm)/(G12f+Gm-Vf*(G12f-Gm))
      
      P21 = 0.5*K23c/K23m*N2
      P22 = G23c/(2.0*Gm)*(A2+3.0*B2*Vf)
      P23 = K23c/(4.0*Gm)*M2/Vf
      P24 = G23c/(2.0*K23m)*(3.0*B2*Vf-D2/Vf)
      P25 = G23c/(4.0*Gm)*(3.0*C2/(Vf**2)+2.0*D2/Vf)
       
      G11AVG = 1.0-2.0*P11+P11**2*(1.0+1.5/Vf**2)
      G22AVG = P21**2+3.0*P22**2+3.0*P23**2+0.5*P24**2+3.0*P25**2
!      G33AVG = G22AVG
      G66AVG = 0.5*(6.0*P22**2+P24**2+6.0*P25**2)
      G12AVG = 2.0*(-1.0+P11)*P21+3.0*P11*P23/Vf
!      G13AVG = G12AVG
      G23AVG = 2.0*P21**2-6.0*P22**2+6.0*P23**2-P24**2-6.0*P25**2

      AVG45 = 0.75*((P41/Vf)**2+P42**2)*(EPS(4)**2+EPS(5)**2) 
      
      V2 = G11AVG*EPS(1)**2+G22AVG*(EPS(2)**2+ EPS(3)**2)+
     1  G12AVG*EPS(1)*(EPS(2)+ EPS(3))+G23AVG*EPS(2)*EPS(3)+
     2  AVG45 + G66AVG*EPS(6)**2
     
      EPSEQ = 1.0/(1.0+NUm)*DSQRT(V2)
     
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SUBTROUTINE TO COMPUTE COEFFICIENTS IN 2-3 PLANE STRAIN PROBLEM
      SUBROUTINE COEFF_23(K23f,G23f,K23m,Gm,Vf,M2,N2,A2,B2,C2,D2)
      
      IMPLICIT NONE
      
      REAL*8::K23f,G23f,K23m,Gm,Vf
      REAL*8::K23c,G23c
      REAL*8::N1,M2,N2,M3
      REAL*8::A1,B1,A2,B2,C2,D2,C3,D3
      REAL*8::ABCD(8,8),ABCD_INV(8,8)
      REAL*8::MN(4,4),MN_INV(4,4)
      REAL*8::F1(8,1),F2(4,1)
      REAL*8::SOL1(8,1),SOL2(4,1)

      CALL CALC23_GSCM(K23f,G23f,K23m,Gm,Vf,K23c,G23c)
      
      ABCD = RESHAPE(
     1 (/-1.0D0, 0.0D0, 1.0D0, 0.0D0, 1.5/(Vf**2), 2.0D0/Vf, 0.0D0, 
     1   0.0D0,
     2   1.0D0, 3.0D0*Vf, -1.0D0, -3.0D0*Vf, 1.5/(Vf**2), 1.0D0/Vf, 
     2   0.0D0, 0.0D0,
     3   -2.0D0, 2.0D0*Vf*(G23f/K23f-1.0D0), 2.0D0*G23f/Gm, 
     3   -2.0*Vf*(Gm/K23m-1.0D0)*G23f/Gm, -G23f/Gm*(1.0D0/Vf)**2,
     3   -2.0D0*(Gm/K23m+1.0D0)*G23f/Gm/Vf, 0.0D0, 0.0D0,
     4   2.0D0, 2.0D0*Vf*(G23f/K23f+2.0D0), -2.0D0*G23f/Gm, 
     4   -2.0*Vf*(Gm/K23m+2.0D0)*G23f/Gm, -G23f/Gm*(1.0D0/Vf)**2,
     4   2.0D0*G23f/K23m/Vf, 0.0D0, 0.0D0,
     5   0.0D0, 0.0D0, -1.0D0, 0.0D0, -1.5D0, -2.0D0, 1.5D0, 2.0D0,
     6   0.0D0, 0.0D0, 1.0D0, 3.0D0, -1.5D0, -1.0D0, 1.5D0, 1.0D0,
     7   0.0D0, 0.0D0, -2.0D0*G23c/Gm, 2.0*(Gm/K23m-1.0D0)*G23c/Gm,
     7   G23c/Gm, 2.0D0*(Gm/K23m+1.0D0)*G23c/Gm, -1.0D0, 
     7   -2.0D0*(1.0D0+G23c/K23c),
     8   0.0D0, 0.0D0, 2.0D0*G23c/Gm, 2.0D0*(Gm/K23m+2.0D0)*G23c/Gm,
     8   G23c/Gm, -2.0D0*G23c/K23m, -1.0D0, 2.0D0*G23c/K23c   
     5 /),(/8,8/))
      
      F1 = 0.0D0
      F1(5,1) = 1.0D0
      F1(6,1) = -1.0D0
      F1(7,1) = 2.0D0
      F1(8,1) = -2.0D0
      CALL matrixInverse(TRANSPOSE(ABCD),ABCD_INV,8,8)
      SOL1 = MATMUL(ABCD_INV,F1)
C      A1 = SOL1(1,1)
C      B1 = SOL1(2,1)
      A2 = SOL1(3,1)
      B2 = SOL1(4,1)
      C2 = SOL1(5,1)
      D2 = SOL1(6,1)
C      C3 = SOL1(7,1)
C      D3 = SOL1(8,1)

      MN = RESHAPE(
     1 (/1.0D0, -0.5/(Vf), -1.0D0, 0.0D0,
     2   2.0*G23f/K23f, G23f/(Gm*Vf), -2.0*G23f/K23m, 0.0D0,
     3   0.0D0, 0.5D0, 1.0D0, -0.5D0,
     4   0.0D0, -G23c/Gm, 2.0*G23c/K23m, 1.0D0
     5 /),(/4,4/))
      
      F2 = 0.0
      F2(3,1) = 1.0
      F2(4,1) = 2.0*G23c/K23c
      CALL matrixInverse(TRANSPOSE(MN),MN_INV,4,4)
      SOL2 = MATMUL(MN_INV,F2)
C      N1 = SOL2(1,1)
      M2 = SOL2(2,1)
      N2 = SOL2(3,1)
C      M3 = SOL2(4,1)
      
      RETURN
      END
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C SUBROUTINES TO TRANSFER THE STRAINS IN POLAR TO CARTISIAN COORDIANTE
      SUBROUTINE Pol2CartStrain(EPSr,EPSt,GAMMArt,THETA,
     1 EPSx,EPSy,GAMMAxy)

      IMPLICIT NONE
      REAL*8::EPSr,EPSt,GAMMArt,THETA
      REAL*8::EPSx,EPSy,GAMMAxy
      
      EPSx = EPSr*(DCOS(THETA))**2 + EPSt*(DSIN(THETA))**2 - 
     1  GAMMArt*DSIN(THETA)*DCOS(THETA)
     
      EPSy = EPSr*(DSIN(THETA))**2 + EPSt**(DCOS(THETA))**2 + 
     1  GAMMArt*DSIN(THETA)*DCOS(THETA)
     
      GAMMAxy = 2.0*((EPSr-EPSt)*DSIN(THETA)*DCOS(THETA) + 
     1  0.5*GAMMArt*(DCOS(THETA)**2-DSIN(THETA)**2))
     
      RETURN
      END

     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C SUBROUTINE TO COMPUTE THE AVERAGE OF AN ARRAY
      SUBROUTINE MEAN(A,A_AVG)
      
      IMPLICIT NONE 
      INTEGER::NUM,II
C      REAL*8,ALLOCATABLE::A(:)
      REAL*8::A(1000)
      REAL*8::SUM
      REAL*8::A_AVG
C      ALLOCATE(A(NUM))
      
      NUM = 1000
      SUM = 0.0
      DO II = 1,NUM
          SUM = SUM + A(II)
      END DO
      
      A_AVG = SUM/NUM
        
      RETURN
      END
      
      

      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C SUBROUTINES TO CHECK THE YIELDING CONDITION          
      
      SUBROUTINE YIELD_CHECK(E,NU,EPS,YIELD,isYIELD) 
      
      IMPLICIT NONE 
      REAL*8::E,NU,MU,YIELD, LAMBDA
      REAL*8::D(6,6)
      REAL*8::EPS(6),SIG(6)
      REAL*8::J2
      REAL*8::isYIELD
      
      
      isYIELD = 0.0 
      MU = 0.5*E/(1.0+NU)
      LAMBDA = E*NU/((1.0+NU)*(1.0-2.0*NU))
      D = 0.0
      D(1,1) = LAMBDA+2.0*MU
      D(1,2) = LAMBDA
      D(1,3) = LAMBDA
      D(2,1) = LAMBDA
      D(2,2) = LAMBDA+2.0*MU
      D(2,3) = LAMBDA
      D(3,1) = LAMBDA
      D(3,2) = LAMBDA
      D(3,3) = LAMBDA+2.0*MU
      D(4,4) = MU
      D(5,5) = MU
      D(6,6) = MU
      SIG = MATMUL(D,EPS)
      CALL J2_CALC(SIG,J2)
        IF (DSQRT(3.0*J2) .GE.YIELD) isYIELD = 1.0
      
      RETURN
      END  
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C SUBROUTINES TO CALCULATE J2 (STRESS)
      SUBROUTINE J2_CALC(SIG,J2)
      
      IMPLICIT NONE
      REAL*8::SIG(6)
      REAL*8::J2
      
      J2 = 1.0/6.0*((SIG(1)-SIG(2))**2+(SIG(2)-SIG(3))**2+
     1 (SIG(3)-SIG(1))**2)+SIG(4)**2+SIG(5)**2+SIG(6)**2
      
      RETURN
      END
           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C SUBROUTINES TO CALCULATE J2_PRIME (STRAIN)
      SUBROUTINE J2P_CALC(EPS,J2P)
      
      IMPLICIT NONE
      REAL*8::EPS(6)
      REAL*8::J2P
      
      J2P = 1.0/6.0*((EPS(1)-EPS(2))**2+(EPS(2)-EPS(3))**2+
     1 (EPS(3)-EPS(1))**2)+1.0/4.0*(EPS(4)**2+EPS(5)**2+EPS(6)**2)
           
      RETURN
      END
    
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C SUBROUTINES TO CALCULATE THE EQUIVALENT STRAIN      
      
      SUBROUTINE EPSEQ_CALC(NU,EPS,EPSEQ)
      
      IMPLICIT NONE 
      REAL*8::NU,EPSEQ
      REAL*8::EPS(6)
      
      EPSEQ = 1.0/(1.0+NU)*DSQRT(0.5*((EPS(1)-EPS(2))**2+
     1 (EPS(2)-EPS(3))**2+(EPS(3)-EPS(1))**2)+
     2  3.0/4.0*(EPS(4)**2+EPS(5)**2+EPS(6)**2))
     
      RETURN
      END
         
         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C SUBROUTINES TO CALCULATE THE EQUIVALENT STRESS
         
      SUBROUTINE SIGEQ_CALC(YIELD,K1,K2,E,EPSCR,HARD,EPSEQ,SIGEQ) 
      
      IMPLICIT NONE
      REAL*8::YIELD,K1,K2,E,EPSEQ,SIGEQ
      REAL*8::EPSCR,SIGCR,HARD

      !HARD = 400 
      !EPSCR = 0.04
      !SIGCR = YIELD - K1/K2*(DEXP(-K2*EPSCR)-DEXP(-K2*(YIELD/E)))
      ! 
      !IF (EPSEQ .LT. EPSCR) THEN
      !  SIGEQ = YIELD - K1/K2*(DEXP(-K2*EPSEQ)-DEXP(-K2*(YIELD/E)))
      !ELSE
      !  SIGEQ = SIGCR + HARD*(EPSEQ-EPSCR) 
      !END IF
      
C      EPSCR = -(DLOG(DEXP(-K2*YIELD/E)+K2/K1*(YIELD-SIGCR)))/K2
      SIGEQ = YIELD - K1/K2*(DEXP(-K2*EPSEQ)-DEXP(-K2*(YIELD/E)))
c      IF (EPSEQ .GT. EPSCR) THEN
c          SIGCR = YIELD - K1/K2*(DEXP(-K2*EPSCR)-DEXP(-K2*(YIELD/E)))
c          SIGEQ = SIGCR + HARD*(EPSEQ-EPSCR)
c      END IF    

      RETURN
      END  

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

