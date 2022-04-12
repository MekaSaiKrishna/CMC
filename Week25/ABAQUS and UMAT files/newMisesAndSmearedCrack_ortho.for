C---------------------------------------------------------
C Problem: Orthotropic Material + Orthotropic SCA
C---------------------------------------------------------
C Crack Planes: Normal along = [0,  1,    0];
C---------------------------------------------------------
C UNIAXIAL TENSILE LOADING CONDITIONS ONLY CONSIDERED
C
C Formulation assumes GIIc=GIIIc and max stress criteria as the failure
C criteria

!     ordering scheme [eps11 eps22 eps33 gam12 gam13 gam23]
!                     [sig11 sig22 sig33 sig12 sig13 sig23]

!     remaining issues:
!    - include characteristic length dependent on angle, not on Abaqus
!    - compression failure
!    - crack closing in compression     
      
      subroutine sca3dortho(nstatv, nprops, noel, stepTime, totalTime, 
     1 dt, L, props, eps, sig, statev, ddsdde, isImplicit)
          
      implicit none
      
      real*8::sig(6)
      real*8::eps(6)
      real*8::ddsdde(6,6)
      real*8::statev
      real*8::props
      real*8::stepTime,totalTime,dt,L
      integer::nstatv,nprops,noel,isImplicit
      
      real*8::damp
      real*8::E1,E2,G12,G23,nu12,nu23
      real*8::sigcr0,tau1cr0,tau2cr0
      real*8::GIc,GIIc,GIIIc  
      real*8::XT,XC,YT,YC,ZT,ZC,S12,S13,S23    
      real*8::peeq
      real*8::isCracked
      real*8::twomu
      real*8::Dco(6,6)
      real*8::Sco(6,6)
C     real*8::Dcr(9,9)
      real*8::Dcr(3,3)
C     real*8::Dda(9,9)
      real*8::Dda(3,3)
      real*8::Dcocr(6,6)
C     real*8::princsig(3)
C     real*8::T(3,3)
      real*8::T(1,3)
C     real*8::epscrmax(3)
      real*8::epscrmax
C     real*8::epscr(9)
      real*8::epscr(3)
C     real*8::epscrold(9)
      real*8::epscrold(3)
C     real*8::N(6,9)
      real*8::N(6,3)
C     real*8::NT(9,6) 
      real*8::NT(3,6)
C     real*8::term99(9,9)
      real*8::term33(3,3)
C     real*8::invTerm99(9,9)
      real*8::invTerm33(3,3)
      real*8::term66(6,6)
C     real*8::term69(6,9)
      real*8::term63(6,3)
      real*8::term
C     real*8::sigCr(9)
      real*8::sigCr(3)

      real*8::PI
      real*8::minExp2
      real*8::minExp3
      
      dimension statev(nstatv),props(nprops)
      
C     !Property list
C     !(1) - Elastic modulus             (E1)
C     !(2) - Elastic modulus             (E2)
C     !(3) - Shear modulus               (G12)
C     !(4) - Poisson's ratio             (nu12)
C     !(5) - Poisson's ratio             (nu23)
C     !(6) - Mode   I critical stress-1  (sigcr0)
C     !(7) - Mode  II critical stress-2  (tau1cr0)
C     !(8) - Mode III critical stress-3  (tau2cr0)
C     !(9) - Mode   I fracture toughness (GIc)
C     !(10)- Mode  II fracture toughness (GIIc)
C     !(11)- Mode III fracture toughness (GIIc)
C     !(12)- Long. Strength-Tensile      (XT)
C     !(13)- Long. Strength-Compressive  (XC)
C     !(14)- Trans.Strength-Tensile      (YT)
C     !(15)- Trans.Strength-Compressive  (YC)
C     !(16)- Trans.Strength-Tensile      (ZT)
C     !(17)- Trans.Strength-Compressive  (ZC)
C     !(18)- Shear-Strength 1-2          (S12)
C     !(19)- Shear-Strength 1-3          (S13)
C     !(20)- Shear-Strength 2-3          (S23)
C     !(21)- damping 
      
C     !State variable list
C     !(1)     - equivalent plastic strain
C     !(2)     - crack flag
C     !(3)     - max. crack strain 1
C     !(4:6)   - normal to crack 1
C     !(7:9)   - old crack strain
C     !(10:15) - total strain (this is necessary for explicit)
C     !(16)    - characteristic length (element length)
C     !(17)    - energy dissipated so far [Crack-1]
C     !(18:20) - Dcr diagonal components  [Crack-1]
      
C     statev(40) = L
      statev(16) = L
      E1       = props(1)  !Elastic modulus
      E2       = props(2)  !Elastic modulus
      G12      = props(3)  !Shear   modulus
      nu12     = props(4)  !Poisson's ratio
      nu23     = props(5)  !Poisson's ratio
      sigcr0   = props(6)  !Mode I   critical stress
      tau1cr0  = props(7)  !Mode II  critical stress
      tau2cr0  = props(8)  !Mode III critical stress
      GIc      = props(9)  !Mode I  fracture toughness
      GIIc     = props(10) !Mode II fracture toughness
      GIIIc    = props(11) !Mode II fracture toughness
      XT       = props(12) !Long. Strength-Tensile
      XC       = props(13) !Long. Strength-Compressive      
      YT       = props(14) !Trans.Strength-Tensile
      YC       = props(15) !Trans.Strength-Compressive
      ZT       = props(16) !Trans.Strength-Tensile
      ZC       = props(17) !Trans.Strength-Compressive
      S12      = props(18) !Shear-Strength 1-2
      S13      = props(19) !Shear-Strength 1-3
      S23      = props(20) !Shear-Strength 2-3
      damp     = props(21) !damping 
      
C    ! Compliance Matrix of the Continuum
      
      G23 = E2/(2.0D0*(1.0D0+nu23))

      Sco = reshape(
     1 (/1.0D0/E1, -nu12/E1, -nu12/E1, 0.0D0, 0.0D0, 0.0D0,
     2   -nu12/E1, 1.0D0/E2, -nu23/E2, 0.0D0, 0.0D0, 0.0D0,
     3   -nu12/E1, -nu23/E2, 1.0D0/E2, 0.0D0, 0.0D0, 0.0D0,
     4    0.0D0,    0.0D0,    0.0D0, 1.0D0/G12, 0.0D0, 0.0D0,
     5    0.0D0,    0.0D0,    0.0D0,  0.0D0, 1.0D0/G12,0.0D0,
     6    0.0D0,    0.0D0,    0.0D0,  0.0D0, 0.0D0, 2.0D0*(1.0D0+nu23)/E2
     7 /),(/6,6/))
      
C    ! Stiffness Matrix of the Continuum
      call matrixInverse(Sco,Dco,6,6)

      peeq      =statev(1)
      isCracked =statev(2)
      
C#######################################################################      
C     Defining PI 
      PI=4.D0*DATAN(1.D0)
C#######################################################################

C       write(*,*) ''
C       write(*,*) 'SQRT(1/2)'
C       write(*,*) DSQRT(1.D0/2.D0)
C       write(*,*) ''
C       write(*,*) 'COS(PI/4)'
C       write(*,*) DCOS(PI/4.D0)

C####################################################################### 
C Defining the T matrix (normals of cracks stacked vertically)
C We have 1 crack: normal +90 degrees to '1'      

      T(1,:) = (/DCOS(PI/2.D0), DSIN(PI/2.D0), 0.0d0/)
C     T(2,:) = (/DCOS(PI/4.D0),-DSIN(PI/4.D0), 0.0d0/)
C     T(3,:) = (/0.0d0, 0.0d0, 1.0d0/)

C####################################################################### 
C Check if element size is small enough
C Here instead of l1,l2,l3 we are looking at L alone i.e. cubic element

      !check if element size is small enough (fibre)
C     if ((2.0d0*GIc*E/sigcr0**2).lt.L) then

      if ((2.0d0*GIc*E1/XT**2).lt.L) then  
      write(*,*) ' '
      write(*,*) 'Element Size too large - 1'
      write(*,*) '(2.0d0*GIc*E1/XT**2) - L'
C     write(*,*) 2.0d0*GIc*E/XT**2,L
      write(*,*) 2.0d0*GIc*E1/XT**2,L
      call myExit()      
      end if 

C       !check if element size is small enough (matrix)
C       minExp2 = dmin1(2.0d0*GIc*E2/YT**2, 2.0d0*GIIc*G12/S12**2, 2.0d0*GIIIc*G23/S23**2)

C       if ((minExp2).lt.L) then  
C       write(*,*) ' '
C       write(*,*) 'Element Size too large - 2'
C C     write(*,*) 'term-1'
C C     write(*,*) 2.0d0*GIc*E2/YT**2
C C     write(*,*) 'term-2'
C C     write(*,*) 2.0d0*GIIc*G12/S12**2
C C     write(*,*) 'term-3'
C C     write(*,*) 2.0d0*GIIIc*G23/S23**2
C       write(*,*) '(minExp2) - L'
C       write(*,*) minExp2,L
C       call myExit()      
C       end if      
      
C       minExp3 = dmin1(2.0d0*GIc*E2/ZT**2, 2.0d0*GIIc*G12/S13**2, 2.0d0*GIIIc*G23/S23**2)

C       if ((minExp3).lt.L) then  
C       write(*,*) ' '
C       write(*,*) 'Element Size too large - 3'
C C     write(*,*) 'term-1'
C C     write(*,*) 2.0d0*GIc*E2/ZT**2
C C     write(*,*) 'term-2'
C C     write(*,*) 2.0d0*GIIc*G12/S12**2
C C     write(*,*) 'term-3'
C C     write(*,*) 2.0d0*GIIIc*G23/S23**2
C       write(*,*) '(minExp3) - L'
C       write(*,*) minExp3,L
C       call myExit()      
C       end if        

      if (isCracked.lt.0.5) then
        ! linear elasticity until crack      
        sig=matmul(Dco,eps)
        ddsdde=Dco
      
        ! check if crack exists

C###########################################################################
C     Initiation of failure i.e. Transition criteria as Max Stress Criteria

C     We know sig(2)=sigma0; 
C     and failure criteria: (sigma0/2).gt.dmin1(XT,YT,S12) 

C     if (dmax1(princsig(1),princsig(2),princsig(3)).gt.sigcr0) then

        if (sig(2).gt.dmin1(2.0d0*XT,2.0d0*YT,2.0d0*S12)) then
        isCracked    =1.0d0
        statev(2)    =1.0d0
        statev(3)    =1.0d-8
C       statev(4)    =1.0d-8
C       statev(5)    =1.0d-8
        statev(4:6)  =T(1,:)
C       statev(9:11) =T(2,:)
C       statev(12:14)=T(3,:)
        endif
C#######################################################################
      
      endif

      if (isCracked.gt.0.5) then
C     epscrmax = statev( 3: 5)
      epscrmax = statev( 3)
C     T(1,;)   = statev( 6: 8)
      T(1,:)   = statev( 4: 6)
C     T(2,:)   = statev( 9:11)
C     T(3,:)   = statev(12:14)
C     epscrold = statev(15:23)
      epscrold = statev( 9:11)
      
      !form the damping matrix
C     term=damp/(2.0d0*(1.0d0+nu))
      term=damp/(2.0d0*(1.0d0+nu12))
      Dda=0.0d0
      Dda(1,1)=damp
      Dda(2,2)=term
      Dda(3,3)=term
C     Dda(4,4)=damp
C     Dda(5,5)=term
C     Dda(6,6)=term
C     Dda(7,7)=damp
C     Dda(8,8)=term
C     Dda(9,9)=term
      
      !find the transformation matrix
      call getN(T,N)

     
      !Determine crack strain consistently
      call qCalcEpscr(epscr,epscrmax,epscrold,eps,Dco,Dda,Dcr,
     # GIC,GIIC,sigcr0,tau1cr0,tau2cr0,L,dt,N,nu12)
C     statev(31) = Dcr(1,1)
C     statev(32) = Dcr(2,2)
C     statev(33) = Dcr(3,3)

      statev(18) = Dcr(1,1)
      statev(19) = Dcr(2,2)
      statev(20) = Dcr(3,3)
      
C     statev(34) = Dcr(4,4)
C     statev(35) = Dcr(5,5)
C     statev(36) = Dcr(6,6)
C     statev(37) = Dcr(7,7)
C     statev(38) = Dcr(8,8)
C     statev(39) = Dcr(9,9)

      !calculate the stress
      !if no Jacobian needed, this is faster: stress=matmul(Dco,eps-matmul(N,epscr))
      !build the jacobian

      !call calcDcr(epscr,epscrmax,GIC,sigcr0,L,Dcr,nu)
C     term99=Dcr+matmul(transpose(N),matmul(Dco,N))+1.0d0/dt*Dda
      term33=Dcr+matmul(transpose(N),matmul(Dco,N))+1.0d0/dt*Dda
C     call matrixInverse(term99,invTerm99,9,9)
      call matrixInverse(term33,invTerm33,3,3)
C     term69=matmul(Dco,matmul(N,invTerm99))
      term63=matmul(Dco,matmul(N,invTerm33))
C     term66=matmul(term69,matmul(transpose(N),Dco))
      term66=matmul(term63,matmul(transpose(N),Dco))
      
      Dcocr=Dco-term66
      
C      sig=matmul(Dcocr,eps)
C     #     -1.0d0/dt*matmul(term69,matmul(Dda,epscrold))
      
      sig=matmul(Dcocr,eps)
     #     -1.0d0/dt*matmul(term63,matmul(Dda,epscrold))

      statev(3) = dmax1(statev(3) , epscr(1))
C     statev(4) = dmax1(statev(4) , epscr(4))
C     statev(5) = dmax1(statev(5) , epscr(7))
      
C     OLD CRACK STRAIN     
      statev(7:9)=epscr 
    
      if (isImplicit.eq.1) ddsdde=Dcocr
      
      sigCr=matmul(Dcr,epscr)

C#######################################################################
C     Energy Calculation: 
C#######################################################################

      !calculate the energy dissipated during fracture
      !this is not correct for two reasons:
      !-the stress should include damping effects
      !-this is a Rieman-Sum, not a trapezoidal rule


C     WILL THIS STILL BE VALID? 
C       statev(24)=statev(24)+dmax1(0.0d0,sigCr(1)*(epscr(1)-epscrold(1)))
C      # *L/GIc
C       statev(25)=statev(25)+dmax1(0.0d0,sigCr(4)*(epscr(4)-epscrold(4)))
C      # *L/GIc
C       statev(26)=statev(26)+dmax1(0.0d0,sigCr(7)*(epscr(7)-epscrold(7)))
C      # *L/GIc

      endif
      
      return
      end

C#######################################################################
C     Defining other SUBROUTINES:
C#######################################################################
!----------------------------------------------------------------------   

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     Subroutine:: getN(o,N):= Finding 'N' matrix 
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Turning a coordinate transformation based on a tensor
! into a coordinate transformation based on a matrix (Voigt Notation)
! Was: sig_ijkl=o_iq*o_jr*o_ks*o_lt*sig_qrst
! Is : sig_i=N_ij*sig_j
! See: TING's book on Anisotropy
! Ordering of the crack strains:
! [epscr11,gamcr12,gamcr12,epscr22,...]


      subroutine getN(o,N)
      implicit none
      
      real*8::o(3,3)
C     real*8::N(6,9)
      real*8::N(6,3)
      
          
      N(1,1) = o(1,1) ** 2
      N(1,2) = 0.1D1 * o(1,1) * o(2,1)
      N(1,3) = 0.1D1 * o(3,1) * o(1,1)
C     N(1,4) = o(2,1) ** 2
C     N(1,5) = 0.1D1 * o(1,1) * o(2,1)
C     N(1,6) = 0.1D1 * o(2,1) * o(3,1)
C     N(1,7) = o(3,1) ** 2
C     N(1,8) = 0.1D1 * o(3,1) * o(1,1)
C     N(1,9) = 0.1D1 * o(2,1) * o(3,1)
      N(2,1) = o(1,2) ** 2
      N(2,2) = 0.1D1 * o(1,2) * o(2,2)
      N(2,3) = 0.1D1 * o(3,2) * o(1,2)
C     N(2,4) = o(2,2) ** 2
C     N(2,5) = 0.1D1 * o(1,2) * o(2,2)
C     N(2,6) = 0.1D1 * o(2,2) * o(3,2)
C     N(2,7) = o(3,2) ** 2
C     N(2,8) = 0.1D1 * o(3,2) * o(1,2)
C     N(2,9) = 0.1D1 * o(2,2) * o(3,2)
      N(3,1) = o(1,3) ** 2
      N(3,2) = 0.1D1 * o(1,3) * o(2,3)
      N(3,3) = 0.1D1 * o(3,3) * o(1,3)
C     N(3,4) = o(2,3) ** 2
C     N(3,5) = 0.1D1 * o(1,3) * o(2,3)
C     N(3,6) = 0.1D1 * o(2,3) * o(3,3)
C     N(3,7) = o(3,3) ** 2
C     N(3,8) = 0.1D1 * o(3,3) * o(1,3)
C     N(3,9) = 0.1D1 * o(2,3) * o(3,3)
      N(4,1) = 0.2D1 * o(1,1) * o(1,2)
      N(4,2) = o(1,1) * o(2,2) + o(1,2) * o(2,1)
      N(4,3) = o(3,1) * o(1,2) + o(3,2) * o(1,1)
C     N(4,4) = 0.2D1 * o(2,1) * o(2,2)
C     N(4,5) = o(1,1) * o(2,2) + o(1,2) * o(2,1)
C     N(4,6) = o(2,1) * o(3,2) + o(2,2) * o(3,1)
C     N(4,7) = 0.2D1 * o(3,1) * o(3,2)
C     N(4,8) = o(3,1) * o(1,2) + o(3,2) * o(1,1)
C     N(4,9) = o(2,1) * o(3,2) + o(2,2) * o(3,1)
      N(5,1) = 0.2D1 * o(1,3) * o(1,1)
      N(5,2) = o(1,3) * o(2,1) + o(1,1) * o(2,3)
      N(5,3) = o(3,3) * o(1,1) + o(3,1) * o(1,3)
C     N(5,4) = 0.2D1 * o(2,3) * o(2,1)
C     N(5,5) = o(1,3) * o(2,1) + o(1,1) * o(2,3)
C     N(5,6) = o(2,3) * o(3,1) + o(2,1) * o(3,3)
C     N(5,7) = 0.2D1 * o(3,3) * o(3,1)
C     N(5,8) = o(3,3) * o(1,1) + o(3,1) * o(1,3)
C     N(5,9) = o(2,3) * o(3,1) + o(2,1) * o(3,3)
      N(6,1) = 0.2D1 * o(1,2) * o(1,3)
      N(6,2) = o(1,2) * o(2,3) + o(1,3) * o(2,2)
      N(6,3) = o(3,2) * o(1,3) + o(3,3) * o(1,2)
C     N(6,4) = 0.2D1 * o(2,2) * o(2,3)
C     N(6,5) = o(1,2) * o(2,3) + o(1,3) * o(2,2)
C     N(6,6) = o(2,2) * o(3,3) + o(2,3) * o(3,2)
C     N(6,7) = 0.2D1 * o(3,2) * o(3,3)
C     N(6,8) = o(3,2) * o(1,3) + o(3,3) * o(1,2)
C     N(6,9) = o(2,2) * o(3,3) + o(2,3) * o(3,2)

      
      return
      end


C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     Subroutine:: qCalcEpscr:= Finding 'epscr' with quad precision
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      
C     ! Same as above, just with quad precision      
C     subroutine qCalcEpscr(epscr,epscrmax,epscrold,eps,Dco,Dda,Dcr,
C    # GIC,sigcr0,L,dt,N,nu)

      subroutine qCalcEpscr(epscr,epscrmax,epscrold,eps,Dco,Dda,Dcr,
     # GIC,GIIC,sigcr0,tau1cr0,tau2cr0,L,dt,N,nu12)
      
      implicit none

C     WHAT IS EPSCRMAX? WHAT IS THE DIMENSION?      
C     real*8::epscr(9),epscrmax(3),eps(6),epscrold(9)
      real*8::epscr(3),epscrmax,eps(6),epscrold(3)
C     real*8::Dco(6,6),Dda(9,9),Dcr(9,9)
      real*8::Dco(6,6),Dda(3,3),Dcr(3,3)
C     real*8::GIC,sigcr0,L,dt,N(6,9),nu12
      real*8::GIC,sigcr0,L,dt,N(6,3),nu12
      real*8::GIIC,tau1cr0,tau2cr0  
      integer::imax,i,k,flag


      
      ! everything internal is handled as quad precision due to 
      ! very small initial cracks
C     real*16::J(9,9),invJ(9,9),Ftol,delta
      real*16::J(3,3),invJ(3,3),Ftol,delta
C     real*16::x(9),dx(9),Fo(9),Fn(9),xPlusDx(9),x0(9) 
      real*16::x(3),dx(3),Fo(3),Fn(3),xPlusDx(3),x0(3) 
C     real*16::qepscrmax(3),qDco(6,6),qDcr(9,9), qeps(6),qepscrold(9)
      real*16::qepscrmax,qDco(6,6),qDcr(3,3), qeps(6),qepscrold(3)
C     real*16::qDda(9,9),qN(6,9),qNT(9,6)
      real*16::qDda(3,3),qN(6,3),qNT(3,6)
      real*16::qGIC,qsigcr0,qL
      real*16::qGIIC,qtau1cr0,qtau2cr0
C     real*16::qdti,qnu
      real*16::qdti,qnu12
C     real*16::qNT_Dco_N(9,9)
      real*16::qNT_Dco_N(3,3)

      flag       = 1            
      qepscrmax  = epscrmax
      qepscrold  = epscrold
      qeps       = eps
      qDco       = Dco
      qDda       = Dda
      qGIC       = GIC
      qGIIC      = qGIIC
      qsigcr0    = sigcr0
      qtau1cr0   = tau1cr0
      qtau2cr0   = tau2cr0
      qL         = L
      qdti       = 1.0d0/dt
      qN         = N
      qNT        = transpose(N)
C     qnu        = nu
      qnu12      = nu12
           
      imax=20  !maximum number of iterations
      Ftol=1.0d-4 !tolerance on the force
      
      !precalculate some stuff
      !N^T*D_co*N
      qNT_Dco_N=matmul(transpose(N),matmul(Dco,N))
      
      !initial guess of crackstrain
      x=0.0d0

      !iterate:
      do i=1,imax
      !calc current force
       call qcalcDcr(x,qepscrmax,qGIC,qGIIC,qsigcr0,qtau1cr0,qtau2cr0,qL,qDcr,qnu12)
      Fo=matmul(qDcr+qNT_Dco_N+qdti*qDda,x)
     #  -matmul(qNT,matmul(qDco,qeps))-matmul(qdti*qDda,qepscrold)  
      !check for convergence 
C      if (qabs(Fo(1))+qabs(Fo(2))+qabs(Fo(3))+qabs(Fo(4))+qabs(Fo(5))
C     #   +qabs(Fo(6))+qabs(Fo(7))+qabs(Fo(8))+qabs(Fo(9)).lt.Ftol) exit 
      if (qabs(Fo(1))+qabs(Fo(2))+qabs(Fo(3)).lt.Ftol) exit 

      !calc jacobian
C         do k=1,9
          do k=1,3
          delta=qsign(1.0q-16,x(k))
          xPlusDx=x
          xPlusDx(k)=x(k)+delta
C      call qcalcDcr(xPlusDx,qepscrmax,qGIC,qsigcr0,qL,qDcr,qnu)
       call qcalcDcr(xPlusDx,qepscrmax,qGIC,qGIIC,qsigcr0,qtau1cr0,qtau2cr0,qL,qDcr,qnu12)
          Fn=matmul(qDcr+qNT_Dco_N+qdti*qDda,xPlusDx)
     #      -matmul(qNt,matmul(qDco,qeps))-matmul(qdti*qDda,qepscrold)
          J(:,k)=(Fn-Fo)/delta
          end do
      !dx=-matmul(J^-1,F)
      dx=Fo !Fo gets replaced by dx in the process
C     call qSolveLinSysLU(J,dx,9,flag)
      call qSolveLinSysLU(J,dx,3,flag)
      if (flag.eq.0) exit
      !x=x+dx <-the minus sign has not been enforced before
      x=x-dx
      end do

      if ((i.ge.imax).or.(flag.eq.0)) then
      !A little bit of diagnostics
C       write(*,*)
      write(*,*) 'Newtons`s method failed to converge in max number'
C       write(*,*) 'of increments'
C       write(*,*) 'Dco='
C       write(*,*) Dco
C       write(*,*) 'Dda='
C       write(*,*) Dda
C       write(*,*) 'eps='
C       write(*,*) eps
C       write(*,*) 'epscrmax='
C       write(*,*) epscrmax
C       write(*,*) 'epscrold'
C       write(*,*) epscrold
C       write(*,*) 'GIC,sigcr0,L,dt'
C       write(*,*) GIC,sigcr0,L,dt
C       write(*,*) 'Transformation matrix'
C       write(*,*) N
C       write(*,*) 'The final values of the iteration'
C       write(*,*) 'epscr='
C       write(*,*) dble(x)
C       write(*,*) 'F'
C       write(*,*) dble(Fo)
 
      
C       call myExit()
C       try to comment the call exit() function. 
      endif 
      
      epscr=dble(x)
      Dcr=dble(qDcr)

      
      return
      end
      

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     Subroutine:: qcalcDcr:= Finding 'Dcr' with quad precision 
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
      
C     ! Same subroutine using quad precision (just change real*8 --> real*16,
C     ! and intrinsic functions if available)
C     subroutine qcalcDcr(epscr,epscrmax,GIC,sigcr0,L,Dcr,nu)
      subroutine qcalcDcr(epscr,epscrmax,GIC,GIIC,sigcr0,tau1cr0,tau2cr0,L,Dcr,nu12)
      implicit none
      
C     DIMENSION OF EPSCRMAX? WHAT IS IT? 
C     real*16::GIC,sigcr0,L,epscr(9),epscrmax(3),Dcr(9,9),nu12,term
      real*16::GIC,sigcr0,L,epscr(3),epscrmax,Dcr(3,3),nu12,term
      real*16::GIIC,tau1cr0,tau2cr0
C     real*16::epscr_temp(9),sigcr0_r
      real*16::epscr_temp(3),sigcr0_r
      real*16::Dcr_min1,Dcr_min2,Dcr_min3
      real*16::tau1cr0_r,tau2cr0_r

      Dcr=0.0Q0
      
C     term=1.0Q0/(2.0Q0*(1.0Q0+nu))
C     term=1.0Q0/(2.0Q0*(1.0Q0+nu12))

C     WHAT DOES EPSCRMAX MEAN ?   
C-----------------------------------------------------------------------
C Is this correct? can this be modified?
      epscr_temp(1) = max(qabs(epscr(1)),qabs(epscrmax))
      epscr_temp(2) = max(qabs(epscr(2)),qabs(epscrmax))
      epscr_temp(3) = max(qabs(epscr(3)),qabs(epscrmax))

C     epscr_temp(4) = max(qabs(epscr(4)),qabs(epscrmax(4)))
C     epscr_temp(5) = max(qabs(epscr(5)),qabs(epscrmax(5)))
C     epscr_temp(6) = max(qabs(epscr(6)),qabs(epscrmax(6)))

C     epscr_temp(7) = max(qabs(epscr(7)),qabs(epscrmax(7)))
C     epscr_temp(8) = max(qabs(epscr(8)),qabs(epscrmax(8)))
C     epscr_temp(9) = max(qabs(epscr(9)),qabs(epscrmax(9)))

      sigcr0_r  = sigcr0*5Q-10
      tau1cr0_r = tau1cr0*5Q-10
      tau2cr0_r = tau2cr0*5Q-10

      !1st direct crack strain
C       if (epscr(1).lt.epscrmax(1)) then
C         Dcr(1,1)=sigcr0*qexp(-sigcr0*L/GIC*epscrmax(1))/epscrmax(1)
C       else
C         Dcr(1,1)=sigcr0*qexp(-sigcr0*L/GIC*epscr(1))/epscr(1)
C       end if
      IF (epscr_temp(1).LT.(2.0Q0*GIC/L/(sigcr0-sigcr0_r))) THEN
       Dcr(1,1) = (-(sigcr0_r-sigcr0)**2/(2.0Q0*GIC/L)*
     1 epscr_temp(1)+sigcr0)/epscr_temp(1)
      ELSE
       Dcr(1,1) = sigcr0_r/epscr_temp(1)
      END IF 

      IF (epscr_temp(2).LT.(2.0Q0*GIIC/L/(tau1cr0-tau1cr0_r))) THEN
       Dcr(2,2) = (-(tau1cr0_r-tau1cr0)**2/(2.0Q0*GIIC/L)*
     1 epscr_temp(2)+tau1cr0)/epscr_temp(2)
      ELSE
       Dcr(2,2) = tau1cr0_r/epscr_temp(2)
      END IF 

      IF (epscr_temp(3).LT.(2.0Q0*GIIC/L/(tau2cr0-tau2cr0_r))) THEN
       Dcr(3,3) = (-(tau2cr0_r-tau2cr0)**2/(2.0Q0*GIIC/L)*
     1 epscr_temp(3)+tau2cr0)/epscr_temp(3)
      ELSE
       Dcr(3,3) = tau2cr0_r/epscr_temp(3)
      END IF

C     !2nd direct crack strain
C C   ###########################
C C       if (epscr(4).lt.epscrmax(2)) then
C C         Dcr(4,4)=sigcr0*qexp(-sigcr0*L/GIC*epscrmax(2))/epscrmax(2)
C C       else
C C         Dcr(4,4)=sigcr0*qexp(-sigcr0*L/GIC*epscr(4))/epscr(4)
C C       end if

C C       if (epscr(4).lt.epscrmax(2)) then
C C         Dcr(4,4)=(-(-sigcr0)**2/(2.0d0*GIC/L)*epscrmax(2)+
C C      1     sigcr0)/epscrmax(2)
C C       else
C C         Dcr(4,4)=(-(-sigcr0)**2/(2.0d0*GIC/L)*epscr(4)+
C C      1     sigcr0)/epscr(4)
C C       end if  

C       IF (epscr_temp(4).LT.(2.0Q0*GIC/L/(sigcr0-sigcr0_r))) THEN
C        Dcr(4,4) = (-(sigcr0_r-sigcr0)**2/(2.0Q0*GIC/L)*
C      1 epscr_temp(4)+sigcr0)/epscr_temp(4)
C       ELSE
C        Dcr(4,4) = sigcr0_r/epscr_temp(4)
C       END IF 

C       IF (epscr_temp(5).LT.(2.0Q0*GIIC/L/(tau1cr0-tau1cr0_r))) THEN
C        Dcr(5,5) = (-(tau1cr0_r-tau1cr0)**2/(2.0Q0*GIIC/L)*
C      1 epscr_temp(5)+tau1cr0)/epscr_temp(5)
C       ELSE
C        Dcr(5,5) = tau1cr0_r/epscr_temp(5)
C       END IF 

C       IF (epscr_temp(6).LT.(2.0Q0*GIIC/L/(tau2cr0-tau2cr0_r))) THEN
C        Dcr(6,6) = (-(tau2cr0_r-tau2cr0)**2/(2.0Q0*GIIC/L)*
C      1 epscr_temp(6)+tau2cr0)/epscr_temp(6)
C       ELSE
C        Dcr(6,6) = tau2cr0_r/epscr_temp(6)
C       END IF 

C       !3rd direct crack strain
C C       if (epscr(7).lt.epscrmax(3)) then
C C         Dcr(7,7)=sigcr0*qexp(-sigcr0*L/GIC*epscrmax(3))/epscrmax(3)
C C       else
C C         Dcr(7,7)=sigcr0*qexp(-sigcr0*L/GIC*epscr(7))/epscr(7)
C C       end if
C       IF (epscr_temp(7).LT.(2.0Q0*GIC/L/(sigcr0-sigcr0_r))) THEN
C        Dcr(7,7) = (-(sigcr0_r-sigcr0)**2/(2.0Q0*GIC/L)*
C      1 epscr_temp(7)+sigcr0)/epscr_temp(7)
C       ELSE
C        Dcr(7,7) = sigcr0_r/epscr_temp(7)
C       END IF

C       IF (epscr_temp(8).LT.(2.0Q0*GIIC/L/(tau1cr0-tau1cr0_r))) THEN
C        Dcr(8,8) = (-(tau1cr0_r-tau1cr0)**2/(2.0Q0*GIIC/L)*
C      1 epscr_temp(8)+tau1cr0)/epscr_temp(8)
C       ELSE
C        Dcr(8,8) = tau1cr0_r/epscr_temp(8)
C       END IF 

C       IF (epscr_temp(9).LT.(2.0Q0*GIIC/L/(tau2cr0-tau2cr0_r))) THEN
C        Dcr(9,9) = (-(tau2cr0_r-tau2cr0)**2/(2.0Q0*GIIC/L)*
C      1 epscr_temp(9)+tau2cr0)/epscr_temp(9)
C       ELSE
C        Dcr(9,9) = tau2cr0_r/epscr_temp(9)
C       END IF  
           
      !construct the rest of the matrix
C       Dcr_min1 = qmin1(Dcr(1,1), Dcr(4,4), Dcr(7,7))
C       Dcr(1,1) = Dcr_min1
C       Dcr(4,4) = Dcr_min1
C       Dcr(7,7) = Dcr_min1

C       Dcr_min2 = qmin1(Dcr(2,2), Dcr(5,5), Dcr(8,8))
C       Dcr(2,2) = Dcr_min2
C       Dcr(5,5) = Dcr_min2
C       Dcr(8,8) = Dcr_min2

C       Dcr_min3 = qmin1(Dcr(3,3), Dcr(6,6), Dcr(9,9))
C       Dcr(3,3) = Dcr_min3
C       Dcr(6,6) = Dcr_min3
C       Dcr(9,9) = Dcr_min3


C      Dcr(2,2)=Dcr(1,1)*term*1.0D-5
C      Dcr(3,3)=Dcr(1,1)*term*1.0D-5
C      Dcr(5,5)=Dcr(4,4)*term*1.0D-5
C      Dcr(6,6)=Dcr(4,4)*term*1.0D-5
C      Dcr(8,8)=Dcr(7,7)*term*1.0D-5
C      Dcr(9,9)=Dcr(7,7)*term*1.0D-5

      
      return
      end
