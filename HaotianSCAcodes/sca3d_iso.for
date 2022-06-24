! ordering scheme [eps11 eps22 eps33 gam12 gam13 gam23]
!                 [sig11 sig22 sig33 sig12 sig13 sig23]

! remaining issues:
! - include characteristic length dependent on angle, not on Abaqus
! - compression failure
! - crack closing in compression     
      
      subroutine sca3diso(nstatv, nprops, noel, stepTime, totalTime, 
     1 dt, L, props, eps, sig, statev, ddsdde, isImplicit)
          
      implicit none
      
      real*8::sig(6)
      real*8::eps(6)
      real*8::ddsdde(6,6)
      real*8::statev
      real*8::props
      real*8::stepTime,totalTime,dt,L
      integer::nstatv,nprops,noel,isImplicit
      
      real*8::E,nu,lambda,mu,damp
      real*8::sigcr0,GIc
      real*8::peeq
      real*8::isCracked, failure
      real*8::twomu
      real*8::Dco(6,6)
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
      
      !Property list
      !(1) - Elastic modulus
      !(2) - Poisson's ratio
      !(3) - Mode I critical streess
      !(4) - Mode I fracture toughness
      !(5) - damping 
      
      !State variable list
      !(1)    - equivalent plastic strain
      !(2)    - crack flag
      !(3)    - max. crack strain 1
      !(4)    - max. crack strain 2
      !(5)    - max. crack strain 3
      !(6:8)  - normal to crack strain 1
      !(9:11) - normal to crack strain 2
      !(12:14)- normal to crack strain 2
      !(15:23)- old crack strain
      !(24:26)- energy dissipated so far
      !(27:32)- total strain (this is necessary for explicit)
      !(33)   - status for failure !!!
      
      isCracked =statev(10)
      
      E      = statev(4)
      nu     = statev(5)
      sigcr0 = props(6)
      GIc    = props(7)
      damp   = props(8)*L
      

      mu= 0.5*E/( 1.0 + nu )
      twomu=2.0d0*mu
      lambda=E*nu/((1.0+nu)*(1.0-2.0*nu))
      
      Dco=reshape(
     1 (/lambda+ twomu,  lambda, lambda, 0.0d0, 0.0d0, 0.0d0, 
     2   lambda, lambda+twomu,  lambda, 0.0d0, 0.0d0, 0.0d0, 
     3   lambda, lambda, lambda+twomu,  0.0d0, 0.0d0, 0.0d0, 
     4   0.0d0,  0.0d0,  0.0d0,  mu,    0.0d0, 0.0d0,
     5   0.0d0,  0.0d0,  0.0d0,  0.0d0, mu,    0.0d0,
     6   0.0d0,  0.0d0,  0.0d0,  0.0d0, 0.0d0, mu 
     7 /),(/6,6/)) 
      
      if (isCracked.lt.0.5) then
        ! check if crack exists
        call mysprind(sig,princsig,T,1,3,3)
      
        if (dmax1(princsig(1),princsig(2),princsig(3)).gt.sigcr0) then
           isCracked    =1.0d0
           statev(10)    =1.0d0
           statev(11)    =1.0d-8
           statev(12)    =1.0d-8
           statev(13)    =1.0d-8
           statev(14:16)  =T(1,:)
           statev(17:19) =T(2,:)
           statev(20:22)=T(3,:)    
        endif
        !ordering of principal direction: T(index of eval,1:3)
      endif

      if (isCracked.gt.0.5) then
         epscrmax = statev( 11: 13)
         T(1,:)   = statev( 14: 16)
         T(2,:)   = statev( 17:19)
         T(3,:)   = statev(20:22)
         epscrold = statev(23:31)
      
         !form the damping matrix
         term=damp/(2.0d0*(1.0d0+nu))
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
         call qcalcEpscr(epscr,epscrmax,epscrold,eps,Dco,Dda,Dcr,
     1 GIC,sigcr0,L,dt,N,nu,failure)      
      
      !calculate the stress
      !if no Jacobian needed, this is faster: stress=matmul(Dco,eps-matmul(N,epscr))
      !build the jacobian
      !call calcDcr(epscr,epscrmax,GIC,sigcr0,L,Dcr,nu)

         term99=Dcr+matmul(transpose(N),matmul(Dco,N))+1.0d0/dt*Dda
         call matrixInverse(term99,invTerm99,9,9)
         term69=matmul(Dco,matmul(N,invTerm99))
         term66=matmul(term69,matmul(transpose(N),Dco))
      
         Dcocr=Dco-term66
      
         sig=matmul(Dco,eps-matmul(N,epscr))

         statev(11) = dmax1(statev(11) , epscr(1))
         statev(12) = dmax1(statev(12) , epscr(4))
         statev(13) = dmax1(statev(13) , epscr(7))
      
         statev(23:31)=epscr
    
         if (isImplicit.eq.1) ddsdde=Dcocr
      
      !calculate the energy dissipated during fracture
      !this is not correct for two reasons:
      !-the stress should include damping effects
      !-this is a Rieman-Sum, not a trapezoidal rule
      ! sigCr=matmul(Dcr,epscr)
      ! statev(32)=statev(32)+dmax1(0.0d0,sigCr(1)*(epscr(1)-epscrold(1)))
      !# *L/GIc
      ! statev(33)=statev(33)+dmax1(0.0d0,sigCr(4)*(epscr(4)-epscrold(4)))
      !# *L/GIc
      ! statev(34)=statev(34)+dmax1(0.0d0,sigCr(7)*(epscr(7)-epscrold(7)))
      !# *L/GIc
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

      
      
! Same as above, just with quad precision      
      subroutine qCalcEpscr(epscr,epscrmax,epscrold,eps,Dco,Dda,Dcr,
     1 GIC,sigcr0,L,dt,N,nu,failure)
      
      implicit none
      real*8::epscr(9),epscrmax(3),eps(6),epscrold(9)
      real*8::Dco(6,6),Dda(9,9),Dcr(9,9)
      real*8::GIC,sigcr0,L,dt,N(6,9),nu
      real*8::failure
      integer::imax,i,k,flag


      
      ! everything internal is handled as quad precision due to 
      ! very small initial cracks
      real*16::J(9,9),invJ(9,9),Ftol,delta
      real*16::x(9),dx(9),Fo(9),Fn(9),xPlusDx(9),x0(9) 
      real*16::qepscrmax(3),qDco(6,6),qDcr(9,9), qeps(6),qepscrold(9)
      real*16::qDda(9,9),qN(6,9),qNT(9,6)
      real*16::qGIC,qsigcr0,qL
      real*16::qdti,qnu
      real*16::qNT_Dco_N(9,9)

      flag       = 1            
      qepscrmax  = epscrmax
      qepscrold  = epscrold
      qeps       = eps
      qDco       = Dco
      qDda       = Dda
      qGIC       = GIC
      qsigcr0    = sigcr0
      qL         = L
      qdti       = 1.0d0/dt
      qN         = N
      qNT        = transpose(N)
      qnu        = nu
           
      imax=20  !maximum number of iterations
      Ftol=1.0q-4 !tolerance on the force
      
      !precalculate some stuff
      !N^T*D_co*N
      qNT_Dco_N=matmul(transpose(N),matmul(Dco,N))
      
      !initial guess of crackstrain
      x=0.0d0

      !iterate:
      do i=1,imax
      !calc current force
       call qcalcDcr(x,qepscrmax,qGIC,qsigcr0,qL,qDcr,qnu,failure)
      Fo=matmul(qDcr+qNT_Dco_N+qdti*qDda,x)
     1  -matmul(qNT,matmul(qDco,qeps))-matmul(qdti*qDda,qepscrold)  
      !check for convergence 
      if (qabs(Fo(1))+qabs(Fo(2))+qabs(Fo(3))+qabs(Fo(4))+qabs(Fo(5))
     1   +qabs(Fo(6))+qabs(Fo(7))+qabs(Fo(8))+qabs(Fo(9)).lt.Ftol) exit 
      !calc jacobian
          do k=1,9
          delta=qsign(1.0q-16,x(k))
          xPlusDx=x
          xPlusDx(k)=x(k)+delta
       call qcalcDcr(xPlusDx,qepscrmax,qGIC,qsigcr0,qL,qDcr,qnu,failure)
          Fn=matmul(qDcr+qNT_Dco_N+qdti*qDda,xPlusDx)
     1      -matmul(qNt,matmul(qDco,qeps))-matmul(qdti*qDda,qepscrold)
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
      !write(*,*)
      !write(*,*) 'Newtons`s method failed to converge in max number'
      !write(*,*) 'of increments'
      !write(*,*) 'Dco='
      !write(*,*) Dco
      !write(*,*) 'Dda='
      !write(*,*) Dda
      !write(*,*) 'eps='
      !write(*,*) eps
      !write(*,*) 'epscrmax='
      !write(*,*) epscrmax
      !write(*,*) 'epscrold'
      !write(*,*) epscrold
      !write(*,*) 'GIC,sigcr0,L,dt'
      !write(*,*) GIC,sigcr0,L,dt
      !write(*,*) 'Transformation matrix'
      !write(*,*) N
      !write(*,*) 'The final values of the iteration'
      !write(*,*) 'epscr='
      !write(*,*) dble(x)
      !write(*,*) 'F'
      !write(*,*) dble(Fo)
      !
      !
      !call myExit()
        write(*,*) 'Check matrix failure!'  
      
      endif 
      
      epscr=dble(x)
      Dcr=dble(qDcr)

      
      return
      end
      
      
      
! ----------------------------------------------------------------------------      
! function to find the crack slope
! using exponentials 

      
! Same subroutine using quad precision (just change real*8 --> real*16,
! and intrinsic functions if available)
      subroutine qcalcDcr(epscr,epscrmax,GIC,sigcr0,L,Dcr,nu,failure)
      implicit none
      
      real*16::GIC,sigcr0,L,epscr(9),epscrmax(3),Dcr(9,9),nu,term
      real*8::failure
      
      Dcr=0.0q0
      
      term=1.0q0/(2.0q0*(1.0q0+nu))
      

      !1st direct crack strain
      if (epscr(1).lt.epscrmax(1)) then
        Dcr(1,1)=sigcr0*qexp(-sigcr0*L/GIC*epscrmax(1))/epscrmax(1)
      else
        Dcr(1,1)=sigcr0*qexp(-sigcr0*L/GIC*epscr(1))/epscr(1)
      end if
      !2nd direct crack strain
      if (epscr(4).lt.epscrmax(2)) then
        Dcr(4,4)=sigcr0*qexp(-sigcr0*L/GIC*epscrmax(2))/epscrmax(2)
      else
        Dcr(4,4)=sigcr0*qexp(-sigcr0*L/GIC*epscr(4))/epscr(4)
      end if
      !3rd direct crack strain
      if (epscr(7).lt.epscrmax(3)) then
        Dcr(7,7)=sigcr0*qexp(-sigcr0*L/GIC*epscrmax(3))/epscrmax(3)
      else
        Dcr(7,7)=sigcr0*qexp(-sigcr0*L/GIC*epscr(7))/epscr(7)
      end if

c$$$        if (qmax1(qabs(epscr(1)),qabs(epscrmax(1))) .gt. 
c$$$     1         0.01/L ) then
c$$$          Dcr(1,1)=sigcr0*QEXP(-sigcr0/GIC*0.01)
c$$$     1     /qmax1(qabs(epscr(1)),qabs(epscrmax(1)))
c$$$           failure = 1.0
c$$$        end if         
c$$$      
c$$$         if (qmax1(qabs(epscr(4)),qabs(epscrmax(4))) .gt. 
c$$$     1         0.01/L ) then
c$$$          Dcr(4,4)=sigcr0*QEXP(-sigcr0/GIC*0.01)
c$$$     1     /qmax1(qabs(epscr(4)),qabs(epscrmax(4)))
c$$$         failure = 1.0
c$$$        end if 
c$$$     
c$$$        if (qmax1(qabs(epscr(7)),qabs(epscrmax(7))) .gt. 
c$$$     1         0.01/L ) then
c$$$          Dcr(7,7)=sigcr0*QEXP(-sigcr0/GIC*0.01)
c$$$     1     /qmax1(qabs(epscr(7)),qabs(epscrmax(7)))
c$$$            failure = 1.0
c$$$        end if       
      !construct the rest of the matrix
      Dcr(2,2)=Dcr(1,1)*term
      Dcr(3,3)=Dcr(1,1)*term
      Dcr(5,5)=Dcr(4,4)*term
      Dcr(6,6)=Dcr(4,4)*term
      Dcr(8,8)=Dcr(7,7)*term
      Dcr(9,9)=Dcr(7,7)*term
      
      
      return
      end
