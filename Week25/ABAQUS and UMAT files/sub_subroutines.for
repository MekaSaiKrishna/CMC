C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     Subroutine:: calcEpscr:= Finding 'epscr'
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

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
      Fo=matmul(qDcr+qNT_Dco_N+qdti*qDda,x)
      !check for convergence
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
      !dx=-matmul(J^-1,F)
      dx=Fo !Fo gets replaced by dx in the process
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

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     Subroutine:: calcDcr:= Finding 'Dcr' 
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
      
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
