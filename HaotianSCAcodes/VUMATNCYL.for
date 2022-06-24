! Last Debugged: 07/18/2013
! VUMAT
! dianyun@umich.edu


      
c User subroutine VUMAT
      subroutine vumat (
c Read only -
     1     jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
c Write only -
     7     stressNew, stateNew, enerInternNew, enerInelasNew )
c
      include 'vaba_param.inc'
c
      dimension jblock(*), props(nprops),density(*), coordMp(*),
     1     charLength(*), strainInc(*),
     2     relSpinInc(*), tempOld(*),
     3     stretchOld(*),
     4     defgradOld(*),
     5     fieldOld(*), stressOld(*),
     6     stateOld(*), enerInternOld(*),
     7     enerInelasOld(*), tempNew(*),
     8     stretchNew(*),
     9     defgradNew(*),
     1     fieldNew(*),
     2     stressNew(*), stateNew(*),
     3     enerInternNew(*), enerInelasNew(*)
c
      character*80 cmname

      parameter (
     1     i_umt_nblock = 1,
     2     i_umt_npt    = 2,
     3     i_umt_layer  = 3,
     4     i_umt_kspt   = 4,
     5     i_umt_noel   = 5 )

      call  vumatXtrArg ( jblock(i_umt_nblock),
     1     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
     7     stressNew, stateNew, enerInternNew, enerInelasNew,
     8     jblock(i_umt_noel), jblock(i_umt_npt),
     9     jblock(i_umt_layer), jblock(i_umt_kspt))

      return
      end

! The VUMAT that has information about element number, integration point number ...
c
      subroutine vumatXtrArg (
c Read only -
     1     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
c Write only -
     7     stressNew, stateNew, enerInternNew, enerInelasNew,
c Read only extra arguments -
     8     nElement, nMatPoint, nLayer, nSecPoint )
C
      include 'vaba_param.inc'
C
C begin: no vaba_param.inc
!      implicit none
!
!      real(4)::stepTime, totalTime, dt, coordMp, charLength,
!     1  props, density, strainInc, relSpinInc,
!     2  tempOld, stretchOld, defgradOld, fieldOld,
!     3  stressOld, stateOld, enerInternOld, enerInelasOld,
!     4  tempNew, stretchNew, defgradNew, fieldNew,
!C Write only (modifiable) variables -
!     5  stressNew, stateNew, enerInternNew, enerInelasNew
!
!      INTEGER:: nblock, ndir, nshr, nstatev, nfieldv, nprops,
!     1  lanneal,nprecd, nElement, nMatPoint, nLayer, nSecPoint
!
!      parameter (j_sys_Dimension = 2)
!      parameter( n_vec_Length = 544  )
!      parameter( maxblk = n_vec_Length  )
!      parameter(i_ipm_sta = -6)
!      character*5 j_ipm_Error
!      parameter(j_ipm_Error = "Error")
!      parameter(j_ipm_Aborted = 20)
!      parameter(r_MaxVal = 1.d+30)
!
C end: no vaba_param.inc
C
c
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1     charLength(nblock), strainInc(nblock,ndir+nshr),
     2     relSpinInc(nblock,nshr), tempOld(nblock),
     3     stretchOld(nblock,ndir+nshr),
     4     defgradOld(nblock,ndir+nshr+nshr),
     5     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6     stateOld(nblock,nstatev), enerInternOld(nblock),
     7     enerInelasOld(nblock), tempNew(nblock),
     8     stretchNew(nblock,ndir+nshr),
     9     defgradNew(nblock,ndir+nshr+nshr),
     1     fieldNew(nblock,nfieldv),
     2     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3     enerInternNew(nblock), enerInelasNew(nblock)
c
c Documentation of extra arguments:
c  nElement: Array of internal element numbers
      dimension nElement(nblock)
c  nMatPoint: Integration point number
c  nLayer   : Layer number for composite shells and layered solids
c  nSecPoint: Section point number within the current layer
c
      character*80 cmname

      real*8::eps(6)
      real*8::sig(6)
      real*8::statev(nstatev)
      real*8::ddsdde(6,6)

      do i=1,nblock

        sig(1)=stressOld(i,1)
        sig(2)=stressOld(i,2)
        sig(3)=stressOld(i,3)
        sig(4)=stressOld(i,4)
        sig(5)=stressOld(i,6)
        sig(6)=stressOld(i,5)
        
        
        statev=stateOld(i,:)
        
      
c$$$        if (cmname(1:2) .eq. 'MA') then
c$$$          
c$$$          statev(32:34)=statev(32:34)+strainInc(i,1:3)
c$$$          statev(35)=statev(35)+strainInc(i,4)*2.0d0
c$$$          statev(36)=statev(36)+strainInc(i,6)*2.0d0
c$$$          statev(37)=statev(37)+strainInc(i,5)*2.0d0
c$$$          eps=statev(32:37)
c$$$          
c$$$          CALL VUMAT_MATRIX(nstatev, nprops, nElement(i), 
c$$$     1      stepTime, totalTime, dt, charLength(i),
c$$$     2      props, eps, sig, statev, ddsdde, 0)

        if (cmname(1:2) .eq. 'LA') then
          statev(24:26)=statev(24:26)+strainInc(i,1:3)
          statev(27)=statev(27)+strainInc(i,4)*2.0d0
          statev(28)=statev(28)+strainInc(i,6)*2.0d0
          statev(29)=statev(29)+strainInc(i,5)*2.0d0
          eps=statev(24:29)
          
          call SCA_LAMINA(nstatev, nprops, nElement(i), 
     1      stepTime, totalTime, dt, charLength(i),
     2      props, eps, sig, statev, ddsdde, 0)
          
        
        else if (cmname(1:2) .eq. 'CA') then
          statev(24:26)=statev(24:26)+strainInc(i,1:3)
          statev(27)=statev(27)+strainInc(i,4)*2.0d0
          statev(28)=statev(28)+strainInc(i,6)*2.0d0
          statev(29)=statev(29)+strainInc(i,5)*2.0d0
          eps=statev(24:29)       
          
          call CCM3D_ANISO(nstatev, nprops, nElement(i), 
     1      stepTime, totalTime, dt, charLength(i),
     2         props, eps, sig, statev, ddsdde, 0)

          
        end if
      !convert to engineering strain
      !using small strain assumption
      !eps approx U-I
!      eps(1)=stretchNew(i,1)-1.0d0
!      eps(2)=stretchNew(i,2)-1.0d0
!      eps(3)=stretchNew(i,3)-1.0d0
!      eps(4)=stretchNew(i,4)*2.0d0
!      eps(5)=stretchNew(i,5)*2.0d0
!      eps(6)=stretchNew(i,6)*2.0d0

           

      !characteristic length is currently taken as provided by ABAQUS
      !it would require to store on direction and all coordinates of
      !each element
      !L=charLength(i)

         

!            !Abaqus remembers the values from the packaging step
!      if (totalTime.lt.1.0d-16) then
!      stateNew(:,27:32)=0.0d0
!      endif

        stressNew(i,1)=sig(1)
        stressNew(i,2)=sig(2)
        stressNew(i,3)=sig(3)
        stressNew(i,4)=sig(4)
        stressNew(i,5)=sig(6)
        stressNew(i,6)=sig(5)
        
c        statev(39)= charLength(i)
c        statev(43)=relSpinInc(i,1)/dt
c        statev(44)=relSpinInc(i,2)/dt
c        statev(45)=relSpinInc(i,3)/dt
c		if (dabs(statev(43)).gt.1.0D3) statev(46)=0.0d0
c		if (dabs(statev(44)).gt.1.0D3) statev(46)=0.0d0
c		if (dabs(statev(45)).gt.1.0D3) statev(46)=0.0d0
        
        stateNew(i,:)=statev

      enddo

      return
      end      
      
      

      INCLUDE 'matrix_inverse.for'
      INCLUDE 'myaux.for'            
      INCLUDE 'ccm3d_aniso.for'
      INCLUDE 'subroutines_for_ccm3d_aniso.for'
      INCLUDE 'my_explicit_support.for'
      INCLUDE 'sca_aniso.for'
      INCLUDE 'sca_lamina.for'
