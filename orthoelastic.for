 !     ORTHOTROPIC ElASTIC
!     Analytic solution used, so no time dependence


      !! Note: Double Exlcamation == Notes to myself

! ordering scheme [eps11 eps22 eps33 gam12 gam13 gam23]
!                 [sig11 sig22 sig33 sig12 sig13 sig23]

! remaining issues:
! - include characteristic length dependent on angle, not on Abaqus
! - compression failure
! - crack closing in compression     
      
      subroutine orthoelas(nstatv, nprops, noel, stepTime, totalTime,
     1 dt, props, eps, sig, statev, ddsdde)

      implicit none
      
      real*8::sig(6)
      real*8::eps(6)
      real*8::ddsdde(6,6)
      real*8::statev
      real*8::props
      real*8::stepTime,totalTime,dt
      integer::nstatv,nprops,noel

      real*8::E1, E2, E3, nu12, nu13, nu23, nu21, nu31, nu32, G12, G13, G23
      real*8::Sco(6,6)       !Compliance Matrix of Continuum
      real*8::Dco(6,6)       !Stiffness Matrix of Continuum     


      dimension statev(nstatv),props(nprops)
      
      !Property list
      !(1)  - Elastic modulus (E1)
      !(2)  - Elastic modulus (E2)
      !(3)  - Elastic modulus (E3)
      !(4)  - Shear modulus   (G12)
      !(5)  - Shear modulus   (G13)
      !(6)  - Shear modulus   (G23)
      !(7)  - Poisson's ratio (nu12)
      !(8)  - Poisson's ratio (nu13)
      !(9)  - Poisson's ratio (nu23)
      
      !State variable list
      !(1:3)  - stress
      !(4:6)  - strain
 
      E1      = props(1)
      E2      = props(2)
      E3      = props(3)
      G12     = props(4)
      G13     = props(5)
      G23     = props(6)
      nu12    = props(7)
      nu13    = props(8)
      nu23    = props(9)

      nu21 = (nu12/E1)*E2
      nu32 = (nu23/E2)*E3
      nu31 = (nu13/E1)*E3
      
     ! Compliance Matrix of the Continuum
      Sco = reshape(
     1 (/1.0D0/E1, -nu21/E2, -nu31/E3, 0.0D0, 0.0D0, 0.0D0;
     2   -nu12/E1, 1.0D0/E2, -nu32/E3, 0.0D0, 0.0D0, 0.0D0;
     3   -nu13/E1, -nu23/E2, 1.0D0/E3, 0.0D0, 0.0D0, 0.0D0;
     4    0.0D0,    0.0D0,    0.0D0, 0.5D0/G12, 0.0D0, 0.0D0;
     5    0.0D0,    0.0D0,    0.0D0,  0.0D0, 0.5D0/G13,0.0D0;
     6    0.0D0,    0.0D0,    0.0D0,  0.0D0, 0.0D0, 0.5D0/G23
     7 /),(/6,6/))
      
     ! Stiffness Matrix of the Continuum
      call matrixInverse(Sco,Dco,6,6)

     ! Finding Global Stress  
      sig = matmul(Dco,eps)

      statev(1:3) = sig
      statev(4:6) = eps
      
      return
      end
