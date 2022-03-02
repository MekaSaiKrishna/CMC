! ordering scheme [eps11 eps22 eps33 gam12 gam13 gam23]
!                 [sig11 sig22 sig33 sig12 sig13 sig23]

! remaining issues:
! - include characteristic length dependent on angle, not on Abaqus
! - compression failure
! - crack closing in compression     
      
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
      
C     real*8::E,nu,lambda,mu,damp
      real*8::damp
      real*8::E1,E2,G12,nu12,nu23
      real*8::sigcr0,GIc
C     GIIC, tau1cr0, and tau1cr0
      real*8::tau1cr0,tau2cr0,GIIC      
      real*8::peeq
      real*8::isCracked
      real*8::twomu
      real*8::Dco(6,6)
      real*8::Sco(6,6)
      real*8::Dcr(9,9)
      real*8::Dda(9,9)
      real*8::Dcocr(6,6)
      real*8::princsig(3)
      real*8::T(3,3)
      real*8::epscrmax(3)
      real*8::epscr(9)
      real*8::epscrold(9)
      real*8::N(6,9)
      real*8::NT(9,6)
      real*8::term99(9,9)
      real*8::invTerm99(9,9)
      real*8::term66(6,6)
      real*8::term69(6,9)
      real*8::term
      real*8::sigCr(9)
      
      
      dimension statev(nstatv),props(nprops)
      
C     !Property list
C     !(1) - Elastic modulus            (E1)
C     !(2) - Elastic modulus            (E2)
C     !(3) - Shear modulus              (G12)
C     !(4) - Poisson's ratio            (nu12)
C     !(5) - Poisson's ratio            (nu23)
C     !(6) - Mode I critical stress     (sigcr0)
C     !(7) - Mode I fracture toughness  (GIC)
C     !(8) - Mode II fracture toughness (GIIC)
C     !(9) - Mode II critical stress-1  (tau1cr0)
C     !(10)- Mode II critical stress-2  (tau2cr0)
C     !(11)- damping 
      
C     !State variable list
C     !(1)    - equivalent plastic strain
C     !(2)    - crack flag
C     !(3)    - max. crack strain 1
C     !(4)    - max. crack strain 2
C     !(5)    - max. crack strain 3
C     !(6:8)  - normal to crack strain 1
C     !(9:11) - normal to crack strain 2
C     !(12:14)- normal to crack strain 2
C     !(15:23)- old crack strain
C     !(24:26)- energy dissipated so far
C     !(27:32)- total strain (this is necessary for explicit)
C     !(33)   - status for failure !!!
      
      statev(40) = L
      E1      = props(1)  !Elastic modulus
      E2      = props(2)  !Elastic modulus
      G12     = props(3)  !Shear   modulus
      nu12    = props(4)  !Poisson's ratio
      nu23    = props(5)  !Poisson's ratio
      sigcr0  = props(6)  !Mode I critical stress
      GIc     = props(7)  !Mode I fracture toughness
      GIIc    = props(8)  !Mode II fracture toughness
      tau1cr0 = props(9)  !Mode II critical stress-1
      tau2cr0 = props(10) !Mode II critical stress-2
      damp    = props(11) !damping 
      

C      mu= 0.5*E/( 1.0d0 + nu )
C      twomu=2.0d0*mu
C      lambda=E*nu/((1.0d0+nu)*(1.0d0-2.0d0*nu))
      
C     Dco=reshape(
C    1 (/lambda+ twomu,  lambda, lambda, 0.0d0, 0.0d0, 0.0d0, 
C    2   lambda, lambda+twomu,  lambda, 0.0d0, 0.0d0, 0.0d0, 
C    3   lambda, lambda, lambda+twomu,  0.0d0, 0.0d0, 0.0d0, 
C    4   0.0d0,  0.0d0,  0.0d0,  mu,    0.0d0, 0.0d0,
C    5   0.0d0,  0.0d0,  0.0d0,  0.0d0, mu,    0.0d0,
C    6   0.0d0,  0.0d0,  0.0d0,  0.0d0, 0.0d0, mu 
C    7 /),(/6,6/)) 

C     ! Compliance Matrix of the Continuum
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
      
      !check if element size is small enough
C     if ((2.0d0*GIc*E/sigcr0**2).lt.L) then
      if ((2.0d0*GIc*E1/sigcr0**2).lt.L) then  
      write(*,*) ' '
      write(*,*) 'Element Size too large'
      write(*,*) '(2.0d0*GIc*E1/sigcr0**2) - L'
C     write(*,*) 2.0d0*GIc*E/sigcr0**2,L
      write(*,*) 2.0d0*GIc*E1/sigcr0**2,L
      call myExit()      
      end if      
      
      if (isCracked.lt.0.5) then
        ! linear elasticity until crack      
        sig=matmul(Dco,eps)
        ddsdde=Dco
      
        ! check if crack exists
        call mysprind(sig,princsig,T,1,3,3)
C ##################################################      
C         T(1,:) = (/1.0d0, 0.0d0, 0.0d0/)
C         T(2,:) = (/0.0d0, 1.0d0, 0.0d0/)
C         T(3,:) = (/1.0d0, 0.0d0, 1.0d0/)

        if (dmax1(princsig(1),princsig(2),princsig(3)).gt.sigcr0) then
        isCracked    =1.0d0
        statev(2)    =1.0d0
C       Why are we setting these(3,4,5) values as max strain?   
        statev(3)    =1.0d-8
        statev(4)    =1.0d-8
        statev(5)    =1.0d-8
        statev(6:8)  =T(1,:)
        statev(9:11) =T(2,:)
        statev(12:14)=T(3,:)
        endif
        
      !ordering of principal direction: T(index of eval,1:3)
      
      endif

      if (isCracked.gt.0.5) then
      epscrmax = statev( 3: 5)
      T(1,:)   = statev( 6: 8)
      T(2,:)   = statev( 9:11)
      T(3,:)   = statev(12:14)
      epscrold = statev(15:23)
      
      !form the damping matrix
C     term=damp/(2.0d0*(1.0d0+nu))
      term=damp/(2.0d0*(1.0d0+nu12))
      Dda=0.0d0
      Dda(1,1)=damp
      Dda(2,2)=term
      Dda(3,3)=term
      Dda(4,4)=damp
      Dda(5,5)=term
      Dda(6,6)=term
      Dda(7,7)=damp
      Dda(8,8)=term
      Dda(9,9)=term
      
      !find the transformation matrix
      call getN(T,N)

     
      !Determine crack strain consistently
C     call qCalcEpscr(epscr,epscrmax,epscrold,eps,Dco,Dda,Dcr,
C    # GIC,sigcr0,L,dt,N,nu)
      call qCalcEpscr(epscr,epscrmax,epscrold,eps,Dco,Dda,Dcr,
     # GIC,GIIC,sigcr0,tau1cr0,tau2cr0,L,dt,N,nu12)
      statev(31) = Dcr(1,1)
      statev(32) = Dcr(2,2)
      statev(33) = Dcr(3,3)
      statev(34) = Dcr(4,4)
      statev(35) = Dcr(5,5)
      statev(36) = Dcr(6,6)
      statev(37) = Dcr(7,7)
      statev(38) = Dcr(8,8)
      statev(39) = Dcr(9,9)

      !calculate the stress
      !if no Jacobian needed, this is faster: stress=matmul(Dco,eps-matmul(N,epscr))
      !build the jacobian
      !call calcDcr(epscr,epscrmax,GIC,sigcr0,L,Dcr,nu)
      term99=Dcr+matmul(transpose(N),matmul(Dco,N))+1.0d0/dt*Dda
      call matrixInverse(term99,invTerm99,9,9)
      term69=matmul(Dco,matmul(N,invTerm99))
      term66=matmul(term69,matmul(transpose(N),Dco))
      
      Dcocr=Dco-term66
      
      sig=matmul(Dcocr,eps)
     #     -1.0d0/dt*matmul(term69,matmul(Dda,epscrold))
      

      statev(3) = dmax1(statev(3) , epscr(1))
      statev(4) = dmax1(statev(4) , epscr(4))
      statev(5) = dmax1(statev(5) , epscr(7))
      
      statev(15:23)=epscr
    
      if (isImplicit.eq.1) ddsdde=Dcocr
      
      !calculate the energy dissipated during fracture
      !this is not correct for two reasons:
      !-the stress should include damping effects
      !-this is a Rieman-Sum, not a trapezoidal rule
      sigCr=matmul(Dcr,epscr)

C     WILL THIS STILL BE VALID? 
      statev(24)=statev(24)+dmax1(0.0d0,sigCr(1)*(epscr(1)-epscrold(1)))
     # *L/GIc
      statev(25)=statev(25)+dmax1(0.0d0,sigCr(4)*(epscr(4)-epscrold(4)))
     # *L/GIc
      statev(26)=statev(26)+dmax1(0.0d0,sigCr(7)*(epscr(7)-epscrold(7)))
     # *L/GIc
      endif
      
      



      return
      end


!----------------------------------------------------------------------      
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
      real*8::N(6,9)
      
          
      N(1,1) = o(1,1) ** 2
      N(1,2) = 0.1D1 * o(1,1) * o(2,1)
      N(1,3) = 0.1D1 * o(3,1) * o(1,1)
      N(1,4) = o(2,1) ** 2
      N(1,5) = 0.1D1 * o(1,1) * o(2,1)
      N(1,6) = 0.1D1 * o(2,1) * o(3,1)
      N(1,7) = o(3,1) ** 2
      N(1,8) = 0.1D1 * o(3,1) * o(1,1)
      N(1,9) = 0.1D1 * o(2,1) * o(3,1)
      N(2,1) = o(1,2) ** 2
      N(2,2) = 0.1D1 * o(1,2) * o(2,2)
      N(2,3) = 0.1D1 * o(3,2) * o(1,2)
      N(2,4) = o(2,2) ** 2
      N(2,5) = 0.1D1 * o(1,2) * o(2,2)
      N(2,6) = 0.1D1 * o(2,2) * o(3,2)
      N(2,7) = o(3,2) ** 2
      N(2,8) = 0.1D1 * o(3,2) * o(1,2)
      N(2,9) = 0.1D1 * o(2,2) * o(3,2)
      N(3,1) = o(1,3) ** 2
      N(3,2) = 0.1D1 * o(1,3) * o(2,3)
      N(3,3) = 0.1D1 * o(3,3) * o(1,3)
      N(3,4) = o(2,3) ** 2
      N(3,5) = 0.1D1 * o(1,3) * o(2,3)
      N(3,6) = 0.1D1 * o(2,3) * o(3,3)
      N(3,7) = o(3,3) ** 2
      N(3,8) = 0.1D1 * o(3,3) * o(1,3)
      N(3,9) = 0.1D1 * o(2,3) * o(3,3)
      N(4,1) = 0.2D1 * o(1,1) * o(1,2)
      N(4,2) = o(1,1) * o(2,2) + o(1,2) * o(2,1)
      N(4,3) = o(3,1) * o(1,2) + o(3,2) * o(1,1)
      N(4,4) = 0.2D1 * o(2,1) * o(2,2)
      N(4,5) = o(1,1) * o(2,2) + o(1,2) * o(2,1)
      N(4,6) = o(2,1) * o(3,2) + o(2,2) * o(3,1)
      N(4,7) = 0.2D1 * o(3,1) * o(3,2)
      N(4,8) = o(3,1) * o(1,2) + o(3,2) * o(1,1)
      N(4,9) = o(2,1) * o(3,2) + o(2,2) * o(3,1)
      N(5,1) = 0.2D1 * o(1,3) * o(1,1)
      N(5,2) = o(1,3) * o(2,1) + o(1,1) * o(2,3)
      N(5,3) = o(3,3) * o(1,1) + o(3,1) * o(1,3)
      N(5,4) = 0.2D1 * o(2,3) * o(2,1)
      N(5,5) = o(1,3) * o(2,1) + o(1,1) * o(2,3)
      N(5,6) = o(2,3) * o(3,1) + o(2,1) * o(3,3)
      N(5,7) = 0.2D1 * o(3,3) * o(3,1)
      N(5,8) = o(3,3) * o(1,1) + o(3,1) * o(1,3)
      N(5,9) = o(2,3) * o(3,1) + o(2,1) * o(3,3)
      N(6,1) = 0.2D1 * o(1,2) * o(1,3)
      N(6,2) = o(1,2) * o(2,3) + o(1,3) * o(2,2)
      N(6,3) = o(3,2) * o(1,3) + o(3,3) * o(1,2)
      N(6,4) = 0.2D1 * o(2,2) * o(2,3)
      N(6,5) = o(1,2) * o(2,3) + o(1,3) * o(2,2)
      N(6,6) = o(2,2) * o(3,3) + o(2,3) * o(3,2)
      N(6,7) = 0.2D1 * o(3,2) * o(3,3)
      N(6,8) = o(3,2) * o(1,3) + o(3,3) * o(1,2)
      N(6,9) = o(2,2) * o(3,3) + o(2,3) * o(3,2)

      
      return
      end

      
!-----------------------------------------------------------------------
! subroutine to find the crack strain at the end of the increment
! using Newton's method (Broyden's method does not work due to too large
! variation in the Jacobian)

C      subroutine calcEpscr(epscr,epscrmax,epscrold,eps,Dco,Dda,Dcr,
C    # GIC,sigcr0,L,dt,N,nu)

       subroutine calcEpscr(epscr,epscrmax,epscrold,eps,Dco,Dda,Dcr,
     # GIC,GIIC,sigcr0,tau1cr0,tau2cr0,L,dt,N,nu12)
      
      implicit none
      real*8::epscr(9),epscrmax(3),eps(6),epscrold(9)
      real*8::Dco(6,6),Dda(9,9),Dcr(9,9)
C     real*8::GIC,sigcr0,L,dt,N(6,9),nu
      real*8::GIC,sigcr0,L,dt,N(6,9),nu12 
      real*8::GIIC,tau1cr0,tau2cr0
      integer::imax,i,k,flag


      
      ! everything internal is handled as quad precision due to 
      ! very small initial cracks
      real*8::J(9,9),invJ(9,9),Ftol,delta
      real*8::x(9),dx(9),Fo(9),Fn(9),xPlusDx(9),x0(9) 
      real*8::qepscrmax(3),qDco(6,6),qDcr(9,9), qeps(6),qepscrold(9)
      real*8::qDda(9,9),qN(6,9),qNT(9,6)
      real*8::qGIC,qsigcr0,qL
      real*8::qGIIC,qtau1cr0,qtau2cr0
C     real*8::qdti,qnu
      real*8::qdti,qnu12
      real*8::qNT_Dco_N(9,9)
      
      flag       = 1      
      qepscrmax  = epscrmax
      qepscrold  = epscrold
      qeps       = eps
      qDco       = Dco
      qDda       = Dda
      qGIC       = GIC
      qGIIC      = GIIC
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
C      call calcDcr(x,qepscrmax,qGIC,qsigcr0,qL,qDcr,qnu)
       call calcDcr(x,qepscrmax,qGIC,qGIIC,qsigcr0,qtau1cr0,qtau2cr0,qL,qDcr,qnu12)
C     !check for convergence
      Fo=matmul(qDcr+qNT_Dco_N+qdti*qDda,x)
     #  -matmul(qNT,matmul(qDco,qeps))-matmul(qdti*qDda,qepscrold)   
      if (dabs(Fo(1))+dabs(Fo(2))+dabs(Fo(3))+dabs(Fo(4))+dabs(Fo(5))
     #   +dabs(Fo(6))+dabs(Fo(7))+dabs(Fo(8))+dabs(Fo(9)).lt.Ftol) exit 
      !calc jacobian
          do k=1,9
          delta=dsign(1.0d-12,x(k))
          xPlusDx=x
          xPlusDx(k)=x(k)+delta
C      call calcDcr(xPlusDx,qepscrmax,qGIC,qsigcr0,qL,qDcr,qnu)
       call calcDcr(xPlusDx,qepscrmax,qGIC,qGIIC,qsigcr0,qtau1cr0,qtau2cr0,qL,qDcr,qnu12)
          Fn=matmul(qDcr+qNT_Dco_N+qdti*qDda,xPlusDx)
     #      -matmul(qNt,matmul(qDco,qeps))-matmul(qdti*qDda,qepscrold)
          J(:,k)=(Fn-Fo)/delta
          end do

C     Why are we not finding 'dx' using dx=-matmul(J^-1,F)          
      !dx=-matmul(J^-1,F)
      dx=Fo !Fo gets replaced by dx in the process

C     Why is it not solveLinSysLU(J,-F,n,flag) 
      call solveLinSysLU(J,dx,9,flag) 
      if (flag.eq.0) exit
      !x=x+dx <-the minus sign has not been enforced before      
      x=x-dx
      end do

 
      if ((i.ge.imax).or.(flag.eq.0)) then
!      !A little bit of diagnostics
!      write(*,*)
!      write(*,*) 'Newtons`s method failed to converge in max number'
!      write(*,*) 'of increments'
!      write(*,*) 'Dco='
!      write(*,*) Dco
!      write(*,*) 'Dda='
!      write(*,*) Dda
!      write(*,*) 'eps='
!      write(*,*) eps
!      write(*,*) 'epscrmax='
!      write(*,*) epscrmax
!      write(*,*) 'epscrold'
!      write(*,*) epscrold
!      write(*,*) 'GIC,sigcr0,L,dt'
!      write(*,*) GIC,sigcr0,L,dt
!      write(*,*) 'Transformation matrix'
!      write(*,*) N
!      write(*,*) 'The final values of the iteration'
!      write(*,*) 'epscr='
!      write(*,*) dble(x)
!      write(*,*) 'F'
!      write(*,*) dble(Fo)
 
      
!      call myExit()
!      write(*,*) 'double precision did not work'
      !write(*,*) 'Entering '
C      call qCalcEpscr(epscr,epscrmax,epscrold,eps,Dco,Dda,Dcr,C
C     # GIC,sigcr0,L,dt,N,nu)
      call qCalcEpscr(epscr,epscrmax,epscrold,eps,Dco,Dda,Dcr,
     # GIC,sigcr0,L,dt,N,nu12)
      !write(*,*) 'Exitting'
      endif 
      
      epscr=dble(x)
      Dcr=dble(qDcr)

      
      return
      end
      
      
C     ! Same as above, just with quad precision      
C     subroutine qCalcEpscr(epscr,epscrmax,epscrold,eps,Dco,Dda,Dcr,
C    # GIC,sigcr0,L,dt,N,nu)

      subroutine qCalcEpscr(epscr,epscrmax,epscrold,eps,Dco,Dda,Dcr,
     # GIC,GIIC,sigcr0,tau1cr0,tau2cr0,L,dt,N,nu12)
      
      implicit none

C     WHAT IS EPSCRMAX? WHAT IS THE DIMENSION?      
      real*8::epscr(9),epscrmax(3),eps(6),epscrold(9)
      real*8::Dco(6,6),Dda(9,9),Dcr(9,9)
C     real*8::GIC,sigcr0,L,dt,N(6,9),nu
      real*8::GIC,sigcr0,L,dt,N(6,9),nu12
      real*8::GIIC,tau1cr0,tau2cr0  
      integer::imax,i,k,flag


      
      ! everything internal is handled as quad precision due to 
      ! very small initial cracks
      real*16::J(9,9),invJ(9,9),Ftol,delta
      real*16::x(9),dx(9),Fo(9),Fn(9),xPlusDx(9),x0(9) 
      real*16::qepscrmax(3),qDco(6,6),qDcr(9,9), qeps(6),qepscrold(9)
      real*16::qDda(9,9),qN(6,9),qNT(9,6)
      real*16::qGIC,qsigcr0,qL
      real*16::qGIIC,qtau1cr0,qtau2cr0
C     real*16::qdti,qnu
      real*16::qdti,qnu12
      real*16::qNT_Dco_N(9,9)

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
C      call qcalcDcr(x,qepscrmax,qGIC,qsigcr0,qL,qDcr,qnu)
       call qcalcDcr(x,qepscrmax,qGIC,qGIIC,qsigcr0,qtau1cr0,qtau2cr0,qL,qDcr,qnu12)
      Fo=matmul(qDcr+qNT_Dco_N+qdti*qDda,x)
     #  -matmul(qNT,matmul(qDco,qeps))-matmul(qdti*qDda,qepscrold)  
      !check for convergence 
      if (qabs(Fo(1))+qabs(Fo(2))+qabs(Fo(3))+qabs(Fo(4))+qabs(Fo(5))
     #   +qabs(Fo(6))+qabs(Fo(7))+qabs(Fo(8))+qabs(Fo(9)).lt.Ftol) exit 
      !calc jacobian
          do k=1,9
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
      call qSolveLinSysLU(J,dx,9,flag)
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
      
      
      
C     ! ----------------------------------------------------------------      
C     ! function to find the crack slope
C     ! using exponentials 

C     subroutine calcDcr(epscr,epscrmax,GIC,sigcr0,L,Dcr,nu12)
      subroutine calcDcr(epscr,epscrmax,GIC,GIIC,sigcr0,tau1cr0,tau2cr0,L,Dcr,nu12)
      implicit none
      
C     real*8::GIC,sigcr0,L,epscr(9),epscrmax(9),Dcr(9,9),nu,term
      real*8::GIC,sigcr0,L,epscr(9),epscrmax(9),Dcr(9,9),nu12,term
      real*8::GIIC,tau1cr0,tau2cr0
      real*8::epscr_temp(9),sigcr0_r
      real*8::tau1cr0_r,tau2cr0_r
      
      Dcr=0.0d0
      
C     term=1.0d0/(2.0d0*(1.0d0+nu))
C     term=1.0d0/(2.0d0*(1.0d0+nu12))

C     WHAT DOES EPSCRMAX MEAN ?     
      epscr_temp(4) = max(dabs(epscr(4)),dabs(epscrmax(4)))
      epscr_temp(5) = max(dabs(epscr(5)),dabs(epscrmax(5)))
      epscr_temp(6) = max(dabs(epscr(6)),dabs(epscrmax(6)))

      sigcr0_r  = sigcr0*5d-3
      tau1cr0_r = tau1cr0*5d-3
      tau2cr0_r = tau2cr0*5d-3

      !1st direct crack strain
C       if (epscr(1).lt.epscrmax(1)) then
C         Dcr(1,1)=sigcr0*dexp(-sigcr0*L/GIC*epscrmax(1))/epscrmax(1)
C       else
C         Dcr(1,1)=sigcr0*dexp(-sigcr0*L/GIC*epscr(1))/epscr(1)
C       end if
      !2nd direct crack strain
C ########################################
C       if (epscr(4).lt.epscrmax(2)) then
C         Dcr(4,4)=sigcr0*dexp(-sigcr0*L/GIC*epscrmax(2))/epscrmax(2)
C       else
C         Dcr(4,4)=sigcr0*dexp(-sigcr0*L/GIC*epscr(4))/epscr(4)
C       end if      

C       if (epscr(4).lt.epscrmax(2)) then
C         Dcr(4,4)=(-(-sigcr0)**2/(2.0d0*GIC/L)*epscrmax(2)+
C      1     sigcr0)/epscrmax(2)
C       else
C         Dcr(4,4)=(-(-sigcr0)**2/(2.0d0*GIC/L)*epscr(4)+
C      1     sigcr0)/epscr(4)
C       end if  

      IF (epscr_temp(4).LT.(2.0d0*GIC/L/(sigcr0-sigcr0_r))) THEN
       Dcr(4,4) = (-(sigcr0_r-sigcr0)**2/(2.0d0*GIC/L)*
     1 epscr_temp(4)+sigcr0)/epscr_temp(4)
      ELSE
       Dcr(4,4) = sigcr0_r/epscr_temp(4)
      END IF 

      IF (epscr_temp(5).LT.(2.0d0*GIIC/L/(tau1cr0-tau1cr0_r))) THEN
       Dcr(5,5) = (-(tau1cr0_r-tau1cr0)**2/(2.0d0*GIIC/L)*
     1 epscr_temp(5)+tau1cr0)/epscr_temp(5)
      ELSE
       Dcr(5,5) = tau1cr0_r/epscr_temp(5)
      END IF

      IF (epscr_temp(6).LT.(2.0d0*GIIC/L/(tau2cr0-tau2cr0_r))) THEN
       Dcr(6,6) = (-(tau2cr0_r-tau2cr0)**2/(2.0d0*GIIC/L)*
     1 epscr_temp(6)+tau2cr0)/epscr_temp(6)
      ELSE
       Dcr(6,6) = tau2cr0_r/epscr_temp(6)
      END IF

      !3rd direct crack strain
C       if (epscr(7).lt.epscrmax(3)) then
C         Dcr(7,7)=sigcr0*dexp(-sigcr0*L/GIC*epscrmax(3))/epscrmax(3)
C       else
C         Dcr(7,7)=sigcr0*dexp(-sigcr0*L/GIC*epscr(7))/epscr(7)
C       end if
      
      !construct the rest of the matrix
      Dcr(1,1) = Dcr(4,4)
      Dcr(7,7) = Dcr(4,4)

      Dcr(2,2) = Dcr(5,5)
      Dcr(3,3) = Dcr(6,6)
      Dcr(8,8) = Dcr(5,5)
      Dcr(9,9) = Dcr(6,6)

C      Dcr(2,2)=Dcr(1,1)*term
C      Dcr(3,3)=Dcr(1,1)*term
C      Dcr(5,5)=Dcr(2,2)*term
C      Dcr(6,6)=Dcr(2,2)*term
C      Dcr(8,8)=Dcr(3,3)*term
C      Dcr(9,9)=Dcr(3,3)*term

      
      return
      end
      
C     ! Same subroutine using quad precision (just change real*8 --> real*16,
C     ! and intrinsic functions if available)
C     subroutine qcalcDcr(epscr,epscrmax,GIC,sigcr0,L,Dcr,nu)
      subroutine qcalcDcr(epscr,epscrmax,GIC,GIIC,sigcr0,tau1cr0,tau2cr0,L,Dcr,nu12)
      implicit none
      
C     DIMENSION OF EPSCRMAX? WHAT IS IT? 
C     real*16::GIC,sigcr0,L,epscr(9),epscrmax(3),Dcr(9,9),nu12,term
      real*16::GIC,sigcr0,L,epscr(9),epscrmax(9),Dcr(9,9),nu12,term
      real*16::GIIC,tau1cr0,tau2cr0
      real*16::epscr_temp(9),sigcr0_r
      real*16::Dcr_min1,Dcr_min2,Dcr_min3
      real*16::tau1cr0_r,tau2cr0_r

      Dcr=0.0Q0
      
C     term=1.0Q0/(2.0Q0*(1.0Q0+nu))
C     term=1.0Q0/(2.0Q0*(1.0Q0+nu12))

C     WHAT DOES EPSCRMAX MEAN ?   
      epscr_temp(1) = max(qabs(epscr(1)),qabs(epscrmax(1)))
      epscr_temp(2) = max(qabs(epscr(2)),qabs(epscrmax(2)))
      epscr_temp(3) = max(qabs(epscr(3)),qabs(epscrmax(3)))

      epscr_temp(4) = max(qabs(epscr(4)),qabs(epscrmax(4)))
      epscr_temp(5) = max(qabs(epscr(5)),qabs(epscrmax(5)))
      epscr_temp(6) = max(qabs(epscr(6)),qabs(epscrmax(6)))

      epscr_temp(7) = max(qabs(epscr(7)),qabs(epscrmax(7)))
      epscr_temp(8) = max(qabs(epscr(8)),qabs(epscrmax(8)))
      epscr_temp(9) = max(qabs(epscr(9)),qabs(epscrmax(9)))

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
      !2nd direct crack strain
C###########################
C       if (epscr(4).lt.epscrmax(2)) then
C         Dcr(4,4)=sigcr0*qexp(-sigcr0*L/GIC*epscrmax(2))/epscrmax(2)
C       else
C         Dcr(4,4)=sigcr0*qexp(-sigcr0*L/GIC*epscr(4))/epscr(4)
C       end if

C       if (epscr(4).lt.epscrmax(2)) then
C         Dcr(4,4)=(-(-sigcr0)**2/(2.0d0*GIC/L)*epscrmax(2)+
C      1     sigcr0)/epscrmax(2)
C       else
C         Dcr(4,4)=(-(-sigcr0)**2/(2.0d0*GIC/L)*epscr(4)+
C      1     sigcr0)/epscr(4)
C       end if  

      IF (epscr_temp(4).LT.(2.0Q0*GIC/L/(sigcr0-sigcr0_r))) THEN
       Dcr(4,4) = (-(sigcr0_r-sigcr0)**2/(2.0Q0*GIC/L)*
     1 epscr_temp(4)+sigcr0)/epscr_temp(4)
      ELSE
       Dcr(4,4) = sigcr0_r/epscr_temp(4)
      END IF 

      IF (epscr_temp(5).LT.(2.0Q0*GIIC/L/(tau1cr0-tau1cr0_r))) THEN
       Dcr(5,5) = (-(tau1cr0_r-tau1cr0)**2/(2.0Q0*GIIC/L)*
     1 epscr_temp(5)+tau1cr0)/epscr_temp(5)
      ELSE
       Dcr(5,5) = tau1cr0_r/epscr_temp(5)
      END IF 

      IF (epscr_temp(6).LT.(2.0Q0*GIIC/L/(tau2cr0-tau2cr0_r))) THEN
       Dcr(6,6) = (-(tau2cr0_r-tau2cr0)**2/(2.0Q0*GIIC/L)*
     1 epscr_temp(6)+tau2cr0)/epscr_temp(6)
      ELSE
       Dcr(6,6) = tau2cr0_r/epscr_temp(6)
      END IF 

      !3rd direct crack strain
C       if (epscr(7).lt.epscrmax(3)) then
C         Dcr(7,7)=sigcr0*qexp(-sigcr0*L/GIC*epscrmax(3))/epscrmax(3)
C       else
C         Dcr(7,7)=sigcr0*qexp(-sigcr0*L/GIC*epscr(7))/epscr(7)
C       end if
      IF (epscr_temp(7).LT.(2.0Q0*GIC/L/(sigcr0-sigcr0_r))) THEN
       Dcr(7,7) = (-(sigcr0_r-sigcr0)**2/(2.0Q0*GIC/L)*
     1 epscr_temp(7)+sigcr0)/epscr_temp(7)
      ELSE
       Dcr(7,7) = sigcr0_r/epscr_temp(7)
      END IF

      IF (epscr_temp(8).LT.(2.0Q0*GIIC/L/(tau1cr0-tau1cr0_r))) THEN
       Dcr(8,8) = (-(tau1cr0_r-tau1cr0)**2/(2.0Q0*GIIC/L)*
     1 epscr_temp(8)+tau1cr0)/epscr_temp(8)
      ELSE
       Dcr(8,8) = tau1cr0_r/epscr_temp(8)
      END IF 

      IF (epscr_temp(9).LT.(2.0Q0*GIIC/L/(tau2cr0-tau2cr0_r))) THEN
       Dcr(9,9) = (-(tau2cr0_r-tau2cr0)**2/(2.0Q0*GIIC/L)*
     1 epscr_temp(9)+tau2cr0)/epscr_temp(9)
      ELSE
       Dcr(9,9) = tau2cr0_r/epscr_temp(9)
      END IF  
           
      !construct the rest of the matrix
      Dcr_min1 = qmin1(Dcr(1,1), Dcr(4,4), Dcr(7,7))
      Dcr(1,1) = Dcr_min1
      Dcr(4,4) = Dcr_min1
      Dcr(7,7) = Dcr_min1

      Dcr_min2 = qmin1(Dcr(2,2), Dcr(5,5), Dcr(8,8))
      Dcr(2,2) = Dcr_min2
      Dcr(5,5) = Dcr_min2
      Dcr(8,8) = Dcr_min2

      Dcr_min3 = qmin1(Dcr(3,3), Dcr(6,6), Dcr(9,9))
      Dcr(3,3) = Dcr_min3
      Dcr(6,6) = Dcr_min3
      Dcr(9,9) = Dcr_min3


C      Dcr(2,w2)=Dcr(1,1)*term*1.0D-5
C      Dcr(3,3)=Dcr(1,1)*term*1.0D-5
C      Dcr(5,5)=Dcr(4,4)*term*1.0D-5
C      Dcr(6,6)=Dcr(4,4)*term*1.0D-5
C      Dcr(8,8)=Dcr(7,7)*term*1.0D-5
C      Dcr(9,9)=Dcr(7,7)*term*1.0D-5

      
      return
      end
