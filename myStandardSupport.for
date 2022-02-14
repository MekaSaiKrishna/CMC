!other subroutines needed in Abaqus/Standard when running SCA
!
!These are utility routines that are differently named in Standard
!and Explicit
!

! Finding principle directions 
      subroutine mysprind(S,PS,AN,LSTR,NDI,NSHR)
      implicit none
      integer::NDI,NSHR,LSTR
      real*8::S(NDI+NSHR)
      real*8::PS(3)
      real*8::AN(3,3)
      call SPRIND(S,PS,AN,LSTR,NDI,NSHR)
      return 
      end
      
      
      
      
! Exitting
      subroutine myExit()
      call XIT
      return
      end