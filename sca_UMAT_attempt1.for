C-----------------------------------------------------------------------
C
C     Smeared Crack Approach -  Attempt1 (Jan31,2022)
C     
C-----------------------------------------------------------------------
      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)
C     
      INCLUDE 'ABA_PARAM.INC'
C     
      CHARACTER*8 CMNAME
C     
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2 PREDEF(1), DPRED(2), PROPS(NPROPS), COORDS(3), DROT(3, 3),
     3 DFGRD0(3, 3), DFGRD1(3, 3), TIME(2)
C
      DIMENSION EELAS(6), EPLAS(6), FLOW(6)
C     USER DEFINED VARIABLE
      DIMENSION MSTRESS(6, 1), VOLUME(1)
      DIMENSION OLD_DSTRIAN(6), OLD_DSTRESS(6), OLD_STRES(NTENS)
      DIMENSION STRAN_total(NTENS)
C     
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0,
     1 ENUMAX=.4999D0, NEWTON=10, TOLER=1.0D-6)
C
C-----------------------------------------------------------------------
C     Matrix Inverse Subroutine "matrixInverse()"
C-----------------------------------------------------------------------
      SUBROUTINE matrixInverse(a,y,n,np)
      ! a is the input matrix
      ! y is the output i.e. inverse
      
      implicit none

      INTEGER np, indx(np)
      REAL a(np,np), y(np,np)

      INTEGER i,j

      ! Set up Identity Matrix as 'y' to calculate 'a' inverse using
      ! "lubksb" subroutine

      do i=1,n
            do j=1,n
                  y(i,j)=0.D0 
            end do 
            y(i,i)=1.D0
      end do 

      ! Perform LU decomposition once
      call ludcmp(a,n,np,indx,d)

      ! Find inverse by columns
      do j=1,n 
            call lubksb(a,n,np,indx,y(1,j))
            !FORTRAN stores 2D matrices by column, so y(1,j) is the
            !address of the j'th column of y
      end do
C     
C     lubksb(a,n,np,indx,b): Solves the set of n linear equations A·X=B. 
C     Here 'a' is input, not as the matrix A but rather as its LU decom-
C     -position,determined by the routine ludcmp. indx is input as the 
C     permutation vector returned by ludcmp. b(1:n) is input as the  
C     right-hand side vector B, and returns with the solution vector X
C
      return
      end

C-----------------------------------------------------------------------      
C     System of Linear Equations Solving Subroutine "solveLinSysLU()"
C-----------------------------------------------------------------------
      SUBROUTINE solveLinSysLU(A,b,n,flag)
      ! solves Ax=b via LU decomposition
      ! double precisiom - Why is it needed and how is it ensured?
      ! INPUT:: A,b
      ! OUTPUT:: x in place of b
      ! flag=1 if successful, flag=0 if not

      implicit none 
      !Feature of Fortran that by default treats all variables that start with 
      !the letters i, j, k, l, m and n as integers and all other variables
      !as real arguments. Implicit None should always be used. It prevents
      !potential confusion in variable types, and makes detection of typographic
      !errors easier.

      INTEGER n, indx(n), flag
      REAL*8 A(n,n), b(n), d 

      call ludcmp(A,n,n,indx,d,flag)

      ! A(n,n):: input,real values, matrix to be decomposed
      !
      !     n :: input,integer,dimension of the square matrix
      !
      !indx(n):: output,integer values, vector that records for the row
      !          permutation effected by the partial pivoting
      !
      !     d :: output, integer, output as 1 or -1 depending on whether
      !          the number of row interchanges was even or odd
      !
      !  flag :: output, integer, error message, 1=success, 0=singular matrix

      if (flag.eq.1) then
            call lubksb(A,n,n,indx,b)
      end if

      return
      end

C-----------------------------------------------------------------------      
C     LUDCMP SUBROUTINE
C-----------------------------------------------------------------------
      SUBROUTINE ludcmp(a,n,np,indx,d,flag)

      implicit none

      INTEGER n,np,indx(n),NMAX,flag
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=10, TINY=1.0E-20) !Largest expected 'n' and a small
      ! number. Given a matirx 'a(1:n,1:n)' with physical dimension np by np,
      ! this routine replaces it by the LU decomposition of a rowwise permutation
      ! of itself. a and n are input. a is output, arranged as in equation
      ! (2.3.14); indx(1:n) is an output vector that records the row permutation
      ! effected by the partial pivoting; 'd' is output as +1 or -1 
      ! depending on whether the number of row interchanges was even or odd
      ! respectively. This routine is used in combination with "lubksb"
      ! to solve linear equations or invert a matrix

      INTEGER i,imax,j,k 
      REAL*8 aamax,dum,sum,vv(NMAX)   !vv stores the implicit scaling of each row
      flag=1                          !So far successful...
      d=1.                            !No row interchanges yet.
      do i=1,n 
      do j=1,n 
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
      end do 
      if (aamax.eq.0) then
            !'Singular Matrix in LUDCMP'
            flag=0
            exit
      end if

      vv(i)=1./aamax                  !Save the scaling 
      end do 

      do j=1,n                        !loop over columns of Crout's Method
            do i=1,j-1                !Eq:2.3.12 except for i=j
                  sum=a(i,j)          
                  do k=1,i-1
                        sum=sum-a(i,k)*a(k,j)
                  end do 
                  a(i,j)=sum
            end do 
            aamax=0.                        !Initialize for the search of the largest pivot element.
            do i=j,n                        !This is i=j of equation(2.3.12)and i=j+1,,,N of equation (2.3.13)
                  sum=a(i,j) 
                  do k=1,j-1 
                        sum=sum-a(i,k)*a(k,j)
                  end do 
                  a(i,j)=sum 
                  dum=vv(i)*abs(sum)        !Figure of merit for the pivot
                  if (dum.ge.aamax) then    !is it better than the best so far
                        imax=i 
                        aamax=dum 
                  end if 
            end do 
            if (j.ne.imax) then             !Do we need to interchange rows?
                  do k=1,n                  !Yes, do so...
                        dum=a(imax,k)
                        a(imax,k)=a(j,k)
                        a(j,k)=dum
                  end do 
                  d=-d                      !... and change the parsity of d
                  vv(imax)=vv(j)            !Also interchange the scale factor
            end if 
            indx(j)=imax
            if (a(j,j).eq.0)a(j,j)=TINY
            if(j.ne.n)then                  !Now finally divide the pivot element
                  dum=1./a(j,j)
                  do i=j+1,n 
                        a(i,j)=a(i,j)*dum
                  end do 
            end if 
      end do                                !Go back for the next column in reduction
      return 
      END

C-----------------------------------------------------------------------      
C     LUBKSB SUBROUTINE: for fwd substitution and backsubstitution
C-----------------------------------------------------------------------
      SUBROUTINE lubksb(a,n,np,indx,b)

      implicit none

      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      ! Solves the set of linear equations A*X=B. Here 'a' is input not
      ! as the matrix A, but rather as its LU decomp., determined by the
      ! routine ludcmp. 'indx' is input as the permutation vector returned
      ! by ludcmp. 'b(1:n)' is input as the right-hand side vector B, and
      ! returns with the solution vector X. a,n,np, and indx are not 
      ! modified by this routine and can be left in place for succesive 
      ! calls with different right hand sides b. This routine takes into
      ! account the possibility that b will begin with many zero elements
      ! so it is efficient for use in matrix inversion

      INTEGER i, ii, j, ll 
      REAL sum 
      ii=0                    !When ii is set to a positive value, it will become the
      do i=1,n                !index of the first non-vanishing element of b. We now do
            ll=indx(i)        !the forward substitution, equation(2.3.6). The only new
            sum=b(ll)         !wrinkle is to unscramble the permutation as we go
            b(ll)=b(i) 
            if(ii.ne.0)then 
                  do j=ii,i-1
                        sum=sum-a(i,j)*b(j)
                  end do
            else if (sum.ne.0)then
                  ii=i        !A non-zero element was encountered, so from now on we will
            end if            !have to do the sums in the loop above
            b(i)=sum
      end do
      do 1=n,1,-1             !Now we do back substitution, equation (2.3.7)
            sum=b(i)
            do j=i+1,n
                  sum=sum-a(i,j)*b(j)
            end do
            b(i)=sum/a(i,i)   !Store a component of the solution vector X
      end do
      return
      end
C-----------------------------------------------------------------------      
C     References:
C     ----------
C     Fortran Numerical Recipes: https://websites.pmc.ucsc.edu/~fnimmo/eart290c_17/NumericalRecipesinF77.pdf
C     SCA Tutorial BOX Folder: https://app.box.com/s/dq3vg9t7gix3nh0v9z39nh3eorfel734 
C-----------------------------------------------------------------------

