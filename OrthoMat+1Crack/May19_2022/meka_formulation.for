C***********************************************************************
C     ! Orthotropic SCA 3D 
C     ! ----------------------------------------------------------------
C     ! Contains the code concerned with the formulation for SCA in 
C     ! orthotropic composite material, this is used in the 
C     ! 'meka_orthoSCA3D.for' file
C     ! ----------------------------------------------------------------
C     ! Last Debugged: 05/19/2022 
C     ! meka1@purdue.edu
C***********************************************************************
C ! Description: There is only one crack whose crack plane normal is 
C !              perpendicular to the fiber direction

C ! ordering scheme [eps11 eps22 eps33 gam12 gam13 gam23]
C !                 [sig11 sig22 sig33 sig12 sig13 sig23]
C !
C ! Formulation assumes GIIc=GIIIc and max stress criteria as the failure
C ! criteria
C !
C ! remaining issues:
C ! - include characteristic length dependent on angle, not on Abaqus
C ! - compression failure
C ! - crack closing in compression
C***********************************************************************
C ! Comments:
C
C ! (1) For Transverse Tension(Mode-3) and Compression(Mode-4) we assume 
C       the crack plane like we did for Longitudinal case
C
C ! (2) Shear23 case is not present in this code
C***********************************************************************
      SUBROUTINE sca3dortho(nstatv, nprops, noel, stepTime, totalTime,
     1 dt, L, props, eps, sig, statev, ddsdde, isImplicit)

      implicit none 

      real*8::sig(6),eps(6)
C       real*8::stress(6)
      real*8::Dco(6,6),Sco(6,6),Dcr(3,3),Dcocr(6,6),ddsdde(6,6),Dda(3,3)
      real*8::statev,props 
      real*8::stepTime,totalTime,dt,L
      integer::nstatv,nprops,noel,isImplicit 

      real*8::E1,E2,G12,G23,nu12,nu23
      real*8::sigcr0,TAUcr12,TAUcr13,TAUcr23
      real*8::XT,XC,YT,YC,ZT
      real*8::GIC_F,GIIC_F,GIC_M,GIIC_M 
      real*8::isCracked

      INTEGER::MODE
      real*8::T(3,3),B(3,3)
      real*8::N(6,3),NT(3,6),N1(6,3) 
      real*8::epscr(3),epscr_old(3),epscr_max(3)
      real*8::term63(6,3),term33(3,3),invTerm33(3,3),term66(6,6),invTerm66(6,6)
      real*8::sigCr(3)

      real*8::temp_sig(6),prc_temp_sig(3),prc_temp_angle(3,3)
      real*8::matrix_crack_ori(3,3)
      integer::max_index,min_index

      real*8::damp,term  

      real*8::FI_FT,FI_FC,FII_12,FII_13

      dimension statev(nstatv),props(nprops)

C***********************************************************************
C       !State Variables List
C-----------------------------------------------------------------------
C       (1)     crack flag
C       (2:4)   old crack strain [local]
C       (5:10)  total strain
C       (11:13) Dcr components
C       (14)    MODE
C       (15)    characteristic length
C       (16)    sigcr0
C       (17:19) 1st row of 'a' matrix
C       (20:22) 2nd row of 'a' matrix
C       (23:25) 3rd row of 'a' matrix
C       (26:28) epscr_max
C       (29:31) sigCr [local crack stress]

C       !Property Variables List
C-----------------------------------------------------------------------
C       (1) Elastic Modulus-1
C       (2) Elastic Modulus-2
C       (3) Shear Modulus-12
C       (4) Poisson's ratio-12
C       (5) Shear Modulus-23
C       (6) Fracture Toughness-Fibre  (Tension)
C       (7) Fracture Toughness-Fibre  (Compression)
C       (8) Fracture Toughness-Matrix (Tension)
C       (9) Fracture Toughness-Matrix (Compression)
C       (10) Tensile     Strength     (Longitudinal)
C       (11) Compressive Strength     (Longitudinal)
C       (12) Tensile     Strength     (Transverse)
C       (13) Compressive Strength     (Transverse)
C       (14) TAUcr12
C       (15) TAUcr23
C       (16) damping
C***********************************************************************
C     MODE = 1: Longitudinal Tension
C     MODE = 2: Longitudinal Compression
C     MODE = 3: Transverse Tension
C     MODE = 4: Transverse Compression
C     MODE = 7: Shear12
C     MODE = 8: Shear13
C-----------------------------------------------------------------------
C     !Property List
      E1      = props(1)
      E2      = props(2)
      G12     = props(3)
      nu12    = props(4)
      G23     = props(5)
      GIC_F   = props(6)
      GIIC_F  = props(7)
      GIC_M   = props(8)
      GIIC_M  = props(9)
      XT      = props(10)
      XC      = props(11)
      YT      = props(12)
      YC      = props(13)
      TAUcr12 = props(14)
      TAUcr23 = props(15)
      damp    = props(16)

C-----------------------------------------------------------------------
      isCracked = statev(1) !Flag, tells us when the transition from linear
C regime to non-linear regime occurs i.e. from damage phase to failure phase
C-----------------------------------------------------------------------
C     !Compliance Matrix of Continuum (Sco)

C       G23=E2/(2.0D0*(1.0D0+nu23))
      
      nu23 = ((E2/G23)/2.0d0)-(1.0d0)

      Sco = reshape(
     1 (/1.0d0/E1, -nu12/E1, -nu12/E1, 0.0d0, 0.0d0, 0.0d0,
     2   -nu12/E1, 1.0d0/E2, -nu23/E2, 0.0d0, 0.0d0, 0.0d0,
     3   -nu12/E1, -nu23/E2, 1.0d0/E2, 0.0d0, 0.0d0, 0.0d0,
     4      0.0d0,  0.0d0,    0.0d0, 1.0d0/G12,0.0d0,0.0d0,
     5      0.0d0,  0.0d0,    0.0d0, 0.0d0,1.0d0/G12,0.0d0,
     6      0.0d0,  0.0d0,    0.0d0 ,0.0d0,0.0d0,1.0d0/G23
     7/),(/6,6/))

C       write(*,*) 'Sco'
C       write(*,*) Sco

C     !Stiffness Matric of the Continuum (Dco)
      call matrixInverse(Sco,Dco,6,6)
C-----------------------------------------------------------------------
C     Bazant Limit (element size is small enough?):
C
C     Here instead of L1,L2,L3 we are looking at L alone i.e. cubic element
C
C     Check if the element size is small enough
C
C       if ((2.0d0*GIC_F*E1/XT**2).lt.L) then
C             write(*,*) ' '
C             write(*,*) 'Element Size TOO large - 1'
C             write(*,*) '(2.0d0*GIC_F*E1/XT**2)-L'
C             write(*,*) (2.0d0*GIC_F*E1/(XT**2)),L
C             call myExit()
C       end if
C
C ! NOTE: In future we will have to add extra bazant limit conditions when
C      !  element isn't cubic, refer eq(22-24) in 
C      !  http://dx.doi.org/10.1016/j.jcomc.2020.100073
C-----------------------------------------------------------------------
C     ! Failure Criteria: MAX STRESS
C
C     (1) Longitudinal Tensile Failure:       sigma11>XT
C     (2) Longitudinal Compressive Failure: |sigma11|>XC
C     (3) Transverse Tensile Failure:         sigma22>YT
C     (4) Transverse Compressive Failure:   |sigma22|>YC
C     (5) Shear Failure (12 and 13):        (sigma12)**2>TAUcr12**2 (or) 
C                                           (sigma13)**2>TAUcr13**2
C-----------------------------------------------------------------------
C     Linear Elastic Regime:(Damage Phase i.e. Pre-Peak)

      if (isCracked .lt. 0.5) then
            sig = matmul(Dco,eps)
            if (isImplicit .eq. 1) ddsdde=Dco 

            FI_FT  = (sig(1)/XT)**2
            FI_FC  = (sig(1)/XC)**2
            FII_12 = (sig(4)/TAUcr12)**2
            FII_13 = (sig(5)/TAUcr12)**2

            !Check if Crack Exists
            !Longitudinal Tensile Failure:
            if ((sig(1) .ge. 0.0d0) .and. (sig(1) .ge. XT)) then
                  isCracked = 1.0d0
                  MODE = 1.0d0 
            end if
            !Longitudinal Compressive Failure:
            if ((sig(1) .le. 0.0d0) .and. (dabs(sig(1)) .ge. XC)) then
                  isCracked = 2.0d0
                  MODE = 2.0d0
            end if

            !Transverse Tensile Failure:
            if ((sig(2) .ge. 0.0d0) .and. (sig(2) .ge. YT)) then
                  isCracked = 3.0d0
                  MODE = 3.0d0 
            end if
            !Transverse Compressive Failure:
            if ((sig(2) .le. 0.0d0) .and. (dabs(sig(2)) .ge. YC)) then
                  isCracked = 4.0d0
                  MODE = 4.0d0

            end if
            ! Shear-12 Failure:
            if (FII_12 .GE. 1.0D0) then
                  isCracked = 7.0d0
                  MODE = 7.0d0
            end if 
            ! Shear-13 Failure:
            if (FII_13 .GE. 1.0D0) then
                  isCracked = 8.0d0
                  MODE = 8.0d0
            end if
            
            statev(14)=MODE

            if (isCracked .gt. 0.5) then
                  if ((MODE .eq. 1) .OR. (MODE .eq. 2)) then
                        statev(16)=dabs(sig(1))

                        statev(17:19)=(/ 1.0d0, 0.0d0, 0.0d0/)
                        statev(20:22)=(/ 0.0d0, 1.0d0, 0.0d0/)
                        statev(23:25)=(/ 0.0d0, 0.0d0, 1.0d0/)
                  end if
                  if ((MODE .eq. 3) .OR. (MODE .eq. 4)) then
                        statev(16)=dabs(sig(2))

                        statev(17:19)=(/ 0.0d0, 1.0d0, 0.0d0/)
                        statev(20:22)=(/-1.0d0, 0.0d0, 0.0d0/)
                        statev(23:25)=(/ 0.0d0, 0.0d0, 1.0d0/)
                  end if
                  if (MODE .eq. 7) then
                        statev(32)=dabs(sig(1))
                        statev(33)=dabs(sig(4))

                        statev(17:19)=(/ 0.0d0, 1.0d0, 0.0d0/)
                        statev(20:22)=(/ 1.0d0, 0.0d0, 0.0d0/)
                        statev(23:25)=(/ 0.0d0, 0.0d0,-1.0d0/)
                  end if 
                  if (MODE .eq. 8) then
                        statev(32)=dabs(sig(1))
                        statev(33)=dabs(sig(5))
                        statev(17:19)=(/ 0.0d0, 0.0d0, 1.0d0/)
                        statev(20:22)=(/ 1.0d0, 0.0d0, 0.0d0/)
                        statev(23:25)=(/ 0.0d0, 1.0d0, 0.0d0/)
                  end if
                  statev(1)=isCracked
                  statev(26:28)=1.0d-8
            end if

C     SPRIND subroutine: https://searchcode.com/codesearch/view/86783693/            
      end if

C     Post-peak Regime (Failure Phase)
      if (isCRACKED .gt. 0.5) then 
            MODE      = INT(STATEV(14))
            epscr_max = statev(26:28)
            epscr_old = statev(2:4)

            if ((MODE .eq. 1) .or. (MODE .eq. 2)) then
                  sigcr0 = statev(16) 
            end if
            if ((MODE .eq. 3) .or. (MODE .eq. 4)) then
                  sigcr0 = statev(16) 
            end if
            if ((MODE .eq. 7) .or. (MODE .eq. 8)) then
                  TAUcr12 = statev(33)
                  TAUcr13 = statev(33)
            end if   
C     !Damping Matrix
            damp = 0.0d0
         
            Dda(1,1) = damp
            Dda(2,2) = damp
            Dda(3,3) = damp
C-----------------------------------------------------------------------
C     !Transformation Matrix (N)
            B      = 0.0d0
            B(1,:) = statev(17:19)
            B(2,:) = statev(20:22)
            B(3,:) = statev(23:25)

            call getN(B,N)
C-----------------------------------------------------------------------
C     !Determine Crack Strain 
            call qcalcEpscr(epscr,epscr_max,epscr_old,eps,Dco,Dda,Dcr,
     1                      GIC_F,GIIC_F,GIC_M,GIIC_M,sigcr0,
     2                      TAUcr12,TAUcr13,TAUcr23,L,dt,N,MODE)

            term33 = Dcr + matmul(transpose(N),matmul(Dco,N)) + Dda/dt 
            call matrixInverse(term33,invTerm33,3,3)
            term63 = matmul(Dco,matmul(N,invTerm33))
            term66 = matmul(term63,matmul(transpose(N),Dco))
            Dcocr = Dco-term66

C           sig = matmul(Dcocr,eps)-matmul(term63,matmul(Dda,epscr_old))*(1.0d0/dt)
            sig = matmul(Dco,eps-matmul(N,epscr)) 

            sigCr = matmul(Dcr,epscr)

            statev(26) = DMAX1(statev(26),DABS(epscr(1)))
            statev(27) = DMAX1(statev(27),DABS(epscr(2)))
            statev(28) = DMAX1(statev(28),DABS(epscr(3)))

            statev(2:4)=epscr 

            statev(11)=Dcr(1,1)
            statev(12)=Dcr(2,2)
            statev(13)=Dcr(3,3)

            statev(29:31)=sigCr

            statev(15)=L
            if (isImplicit .eq. 1) ddsdde = Dcocr

      end if

      return
      end

C !----------------------------------------------------------------------- 
C !     Defining other SUBROUTINES:
C !-----------------------------------------------------------------------   
C !     Subroutine:: getN(o,N):= Finding 'N' matrix 
C !----------------------------------------------------------------------- 

C ! Turning a coordinate transformation based on a tensor
C ! into a coordinate transformation based on a matrix (Voigt Notation)
C ! Was: sig_ijkl=o_iq*o_jr*o_ks*o_lt*sig_qrst
C ! Is : sig_i=N_ij*sig_j
C ! See: TING's book on Anisotropy
C ! Ordering of the crack strains:
C ! [epscr11,gamcr12,gamcr12,epscr22,...]


      subroutine getN(o,N)
      implicit none
      
      real*8::o(3,3)
      real*8::N(6,3)
      
      N(1,1) = o(1,1) ** 2
      N(1,2) = o(1,1) * o(2,1)
      N(1,3) = o(3,1) * o(1,1)
      N(2,1) = o(1,2) ** 2
      N(2,2) = o(1,2) * o(2,2)
      N(2,3) = o(3,2) * o(1,2)
      N(3,1) = o(1,3) ** 2
      N(3,2) = o(1,3) * o(2,3)
      N(3,3) = o(3,3) * o(1,3)
      N(4,1) = 2.0d0 * o(1,1) * o(1,2)
      N(4,2) = o(1,1) * o(2,2) + o(1,2) * o(2,1)
      N(4,3) = o(3,1) * o(1,2) + o(3,2) * o(1,1)
      N(5,1) = 2.0d0 * o(1,3) * o(1,1)
      N(5,2) = o(1,3) * o(2,1) + o(1,1) * o(2,3)
      N(5,3) = o(3,3) * o(1,1) + o(3,1) * o(1,3)
      N(6,1) = 2.0d0 * o(1,2) * o(1,3)
      N(6,2) = o(1,2) * o(2,3) + o(1,3) * o(2,2)
      N(6,3) = o(3,2) * o(1,3) + o(3,3) * o(1,2)

      return
      end

C !----------------------------------------------------------------------- 
C !     Subroutine:: qCalcEpscr:= Finding 'epscr' with quad precision
C !----------------------------------------------------------------------- 
      subroutine qcalcEpscr(epscr,epscr_max,epscr_old,eps,Dco,Dda,Dcr,
     1                      GIC_F,GIIC_F,GIC_M,GIIC_M,sigcr0,
     2                      TAUcr12,TAUcr13,TAUcr23,L,dt,N,MODE)

      implicit none
      real*8::epscr(3),epscr_old(3),eps(6),epscr_max(3)
      real*8::Dco(6,6),Dda(3,3),Dcr(3,3)
      real*8::sigcr0,TAUcr12,TAUcr13,TAUcr23
      real*8::GIC_F,GIIC_F,GIC_M,GIIC_M
      real*8::L,dt
      real*8::N(6,3) 

      real*16::qepscr_old(3),qeps(6),qepscr_max(3)
      real*16::qDco(6,6),qDda(3,3),qDcr(3,3)
      real*16::qsigcr0,qTAUcr12,qTAUcr13,qTAUcr23
      real*16::qGIC_F,qGIIC_F,qGIC_M,qGIIC_M
      real*16::qL,qdt

      real*16::qNT_Dco_N(3,3),qNT(3,6),qN(6,3)

      real*16::Ftol,delta
      real*16::J(3,3),invJ(3,3)
      real*16::x(3),dx(3),Fo(3),Fn(3),xPlusDx(3),x0(3)

      integer::imax,i,k,flag

      integer::MODE 

      flag       = 1
      qepscr_old = epscr_old
      qeps       = eps 
      qepscr_max = epscr_max
      qDco       = Dco 
      qDda       = Dda 
      qGIC_F     = GIC_F
      qGIIC_F    = GIIC_F
      qGIC_M     = GIC_M
      qGIIC_M    = GIIC_M 
      qsigcr0    = sigcr0
      qTAUcr12   = TAUcr12
      qTAUcr13   = TAUcr13
      qTAUcr23   = TAUcr23 
      qL         = L 
      qdt        = dt 
      qN         = N 
      qNT        = transpose(N)

      qNT_Dco_N = matmul(transpose(N),matmul(Dco,N))

      imax = 20        !max number of iterations
      Ftol = 1.0d-4    !tolerance on f(x)

C     !initial guess of crack strain
      x = 0.0d0

C     !iterate:
      do i=1, imax
            call qcalcDcr(x,qepscr_max,qGIC_F,qGIIC_F,qGIC_M,qGIIC_M,
     1                   qsigcr0,qTAUcr12,qTAUcr13,qTAUcr23,qL,qDcr,MODE)

            Fo=matmul(qDcr + qNT_Dco_N + qDda/qdt,x) - matmul(qNT,matmul(qDco,qeps))
     1        -matmul(qDda/qdt,qepscr_old)

            if (qabs(Fo(1))+qabs(Fo(2))+qabs(Fo(3)) .lt. Ftol) exit 

            !Calculate Jacobian
            do k=1,3
                  delta   = qsign(1.0q-16,x(k))
                  xPlusDx = x
                  xPlusDx(k) = x(k)+delta
                  call qcalcDcr(xPlusDx,qepscr_max,qGIC_F,qGIIC_F,qGIC_M,qGIIC_M,
     1                          qsigcr0,qTAUcr12,qTAUcr13,qTAUcr23,qL,qDcr,MODE)

                  Fn=matmul(qDcr+qNT_Dco_N+(qDda/qdt),xPlusDx)
     1              -matmul(qNT,matmul(qDco,qeps))-matmul(qDda/qdt,qepscr_old)

                  J(:,k)=(Fn-Fo)/delta
            end do 
            dx=Fo 
            call qSolveLinSysLU(J,dx,3,flag)
            if (flag. eq. 0) exit
            x = x - dx 
      end do 
      
      if ((i .ge. imax) .or. (flag .eq. 0)) then
        write(*,*) 'Newton`s method fails to converge in max number'
        write(*,*) 'i:'
        write(*,*) i
        write(*,*) 'flag:'
        write(*,*) flag
C       write(*,*) 'of increments'
C       write(*,*) i
C       write(*,*) 'Dco='
C       write(*,*) Dco
C       write(*,*) 'Dda='
C       write(*,*) Dda
C       write(*,*) 'eps='
C       write(*,*) eps
C       write(*,*) 'epscr_old'
C       write(*,*) epscr_old
C       write(*,*) 'GIC_F,sigcr0,L,dt'
C       write(*,*) GIC_F,sigcr0,L,dt
C       write(*,*) 'Transformation matrix'
C       write(*,*) N
C       write(*,*) 'The final values of the iteration'
C       write(*,*) 'epscr='
C       write(*,*) dble(x)
C       write(*,*) 'F'
C       write(*,*) dble(Fo)
      end if 

      epscr = dble(x)
      Dcr   = dble(qDcr)

      return 
      end

C !----------------------------------------------------------------------- 
C !     Subroutine:: qcalcDcr:= Finding 'Dcr' with quad precision 
C !-----------------------------------------------------------------------       
      
C     ! Same subroutine using quad precision (change real*8 --> real*16,
C     ! and intrinsic functions if available)

      subroutine qcalcDcr(epscr,epscr_max,GIC_F,GIIC_F,GIC_M,GIIC_M,
     1                    sigcr0,TAUcr12,TAUcr13,TAUcr23,L,Dcr,MODE)

      implicit none 

      real*16::epscr(3),epscr_max(3),Dcr(3,3)
      real*16::GIC_F,GIIC_F,GIC_M,GIIC_M 
      real*16::sigcr0,sigcr0_r
      real*16::TAUcr12,TAUcr13,TAUcr23
      real*16::TAUcr12_R,TAUcr13_R,TAUcr23_R
      real*16::L
      real*16::DNORM,DSHEAR1,DSHEAR2

      integer::MODE 

      DNORM   = qmax1(qabs(epscr(1)),qabs(epscr_max(1)))
      DSHEAR1 = qmax1(qabs(epscr(2)),qabs(epscr_max(2)))
      DSHEAR2 = qmax1(qabs(epscr(3)),qabs(epscr_max(3)))

      Dcr=0.0Q0

      TAUcr12_R = TAUcr12*1.0Q-5
      TAUcr13_R = TAUcr13*1.0Q-5
      TAUcr23_R = TAUcr23*1.0Q-5

      if (MODE .eq. 1) sigcr0_r = sigcr0*1.0Q-5
      if (MODE .eq. 2) then 
            sigcr0_r = sigcr0*0.5Q0
            GIC_F=GIIC_F
      end if

      if (MODE .eq. 3) sigcr0_r = sigcr0*1.0Q-8
C       if (MODE .eq. 4) sigcr0_r = sigcr0*1.0Q-8
      if (MODE .eq. 4) then 
            sigcr0_r = sigcr0*1.0Q-5
            GIC_M=GIIC_M
      end if

      if ((MODE .eq. 1) .OR. (MODE .eq. 2)) then
            if (DNORM .LT. (2.0q0*GIC_F/L/(sigcr0-sigcr0_r))) then
                  Dcr(1,1) = (-(sigcr0_r-sigcr0)**2/(2.0q0*GIC_F/L)*
     1            DNORM+sigcr0)/DNORM
            else 
                  Dcr(1,1)=sigcr0_r/DNORM
            end if 
      end if

      if ((MODE .eq. 3) .OR. (MODE .eq. 4)) then
            if (DNORM .LT. (2.0q0*GIC_M/L/(sigcr0-sigcr0_r))) then
                  Dcr(1,1) = (-(sigcr0_r-sigcr0)**2/(2.0q0*GIC_M/L)*
     1            DNORM+sigcr0)/DNORM
            else 
                  Dcr(1,1)=sigcr0_r/DNORM
            end if 
      end if

      if ((MODE .eq. 7) .or. (MODE .eq. 8)) then
            if (DSHEAR1 .LT. (2.0q0*GIIC_M/L/(TAUcr12-TAUcr12_R))) then
                  Dcr(2,2) = (-(TAUcr12_R-TAUcr12)**2/(2.0q0*GIIC_M/L)*
     1            DSHEAR1+TAUcr12)/DSHEAR1
            else
                  Dcr(2,2) = TAUcr12_R/DSHEAR1
            end if 
            Dcr(1,1) = Dcr(2,2)
            Dcr(3,3) = Dcr(2,2)
      end if

      return
      end
