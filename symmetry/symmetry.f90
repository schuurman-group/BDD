!#######################################################################
! symmetry: routines for the determination of whether a given coupling
!           coefficient is zero by symmetry
!#######################################################################

module symmetry

  use constants
  
  implicit none

  !----------------------------------------------------------------------
  ! Symmetry information
  !----------------------------------------------------------------------
  integer, allocatable           :: kappa_mask(:,:)
  integer, allocatable           :: lambda_mask(:,:,:)
  integer, allocatable           :: gamma_mask(:,:,:)
  integer, allocatable           :: mu_mask(:,:,:,:)
  integer, allocatable           :: iota_mask(:,:)
  integer, allocatable           :: tau_mask(:,:,:)
  integer, allocatable           :: epsilon_mask(:,:)
  integer, allocatable           :: xi_mask(:,:,:)
  integer, allocatable           :: dip1_mask(:,:,:,:)
  integer, allocatable           :: dip2_mask(:,:,:,:,:)
  integer, allocatable           :: dip3_mask(:,:,:,:)
  integer, allocatable           :: dip4_mask(:,:,:,:)
  integer, dimension(2)          :: nkappa
  integer, dimension(2)          :: nlambda
  integer, dimension(2)          :: ngamma
  integer, dimension(2)          :: nmu
  integer, dimension(2)          :: niota
  integer, dimension(2)          :: ntau
  integer, dimension(2)          :: nepsilon
  integer, dimension(2)          :: nxi
  integer, dimension(2)          :: ntot
  integer, dimension(2)          :: ndip1
  integer, dimension(2)          :: ndip2
  integer, dimension(2)          :: ndip3
  integer, dimension(2)          :: ndip4
  integer, dimension(2)          :: ndiptot
  integer, allocatable           :: cut_mask(:,:)
  character(len=3), allocatable  :: stalab(:)
  character(len=3), dimension(3) :: diplab
  
!----------------------------------------------------------------------
! Point group information
!----------------------------------------------------------------------
! pntgrp: name of the point group
! irreplab: irrep labels
! chdim: order of the point group
! chel: elements of the character table
!----------------------------------------------------------------------
  character(len=16)                     :: pntgrp
  character(len=16), dimension(8), save :: irreplab
  integer, save                         :: chdim
  real(dp), dimension(8,8),save         :: chel
  
contains

!######################################################################

  subroutine create_mask(dipolein)

    use constants
    use sysinfo
    
    implicit none

    integer                    :: s1,s2,m1,m2,c
    integer, dimension(3)      :: dipchk
    integer, dimension(nmodes) :: nmchk
    integer, dimension(nsta)   :: stachk
    logical, optional          :: dipolein
    logical                    :: dipole

!----------------------------------------------------------------------
! Optional fitting of diabatic dipole matrix elements
!----------------------------------------------------------------------
    dipole=.false.
    if (present(dipolein)) then
       dipole=dipolein
    endif
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! First-order, intrastate
    allocate(kappa_mask(nmodes,nsta))
    kappa_mask=0

    ! First-order, interstate
    allocate(lambda_mask(nmodes,nsta,nsta))
    lambda_mask=0

    ! Second-order, intrastate
    allocate(gamma_mask(nmodes,nmodes,nsta))
    gamma_mask=0

    ! Second-order, interstate
    allocate(mu_mask(nmodes,nmodes,nsta,nsta))
    mu_mask=0

    ! Third-order, cubic-only, intrastate
    allocate(iota_mask(nmodes,nsta))
    iota_mask=0

    ! Third-order, cubic-only, interstate
    allocate(tau_mask(nmodes,nsta,nsta))
    tau_mask=0

    ! Fourth-order, quartic-only, intrastate
    allocate(epsilon_mask(nmodes,nsta))
    epsilon_mask=0

    ! Fourth-order, quartic-only, interstate
    allocate(xi_mask(nmodes,nsta,nsta))
    xi_mask=0

    ! Diabatic dipole matrix expansion coefficients
    if (dipole) then
       allocate(dip1_mask(nmodes,nsta,nsta,3))
       dip1_mask=0
       allocate(dip2_mask(nmodes,nmodes,nsta,nsta,3))
       dip2_mask=0
       allocate(dip3_mask(nmodes,nsta,nsta,3))
       dip3_mask=0
       allocate(dip4_mask(nmodes,nsta,nsta,3))
       dip4_mask=0
    endif
    
!----------------------------------------------------------------------
! First, convert the user specified symmetry labels to the format used
! in the code, then get the characters.
!----------------------------------------------------------------------
    ! Clean up the normal mode symmetry labels
    do m1=1,nmodes
       call cleanirreps(nmlab(m1),len(nmlab))
    enddo

    ! Clean up the state symmetry labels
    do s1=1,nsta
       call cleanirreps(stalab(s1),len(stalab))
    enddo

    ! Clean up the dipole symmetry labels
    if (dipole) then
       do c=1,3
          call cleanirreps(diplab(c),len(diplab(c)))
       enddo
    endif

    ! Get the characters
    call getcharacters

!----------------------------------------------------------------------
! Determine the mask values from the calculation of the number of
! times that the totally symmetric irreducible representation enters
! into the decomposition of the direct product representations
! corresponding to the integrands of the parameters of the diabatic
! potential, i.e., phi_a* d^n/dq1...dqn phi_b
!----------------------------------------------------------------------
    ! Set dipchk to zero until we reach the dipole expansion terms
    dipchk=0
    
    ! kappa
    do s1=1,nsta
       do m1=1,nmodes
          nmchk=0
          stachk=0
          nmchk(m1)=1
          stachk(s1)=2
          kappa_mask(m1,s1)=integralsym(nmchk,stachk,dipchk)                
       enddo
    enddo

    ! lambda
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m1=1,nmodes
             nmchk=0
             stachk=0
             nmchk(m1)=1
             stachk(s1)=1
             stachk(s2)=1
             lambda_mask(m1,s1,s2)=integralsym(nmchk,stachk,dipchk)
             lambda_mask(m1,s2,s1)=lambda_mask(m1,s1,s2)
          enddo
       enddo
    enddo

    ! gamma
    do s1=1,nsta
       do m1=1,nmodes
          do m2=m1,nmodes
             nmchk=0
             stachk=0
             nmchk(m1)=nmchk(m1)+1
             nmchk(m2)=nmchk(m2)+1
             stachk(s1)=2
             gamma_mask(m1,m2,s1)=integralsym(nmchk,stachk,dipchk)
             gamma_mask(m2,m1,s1)=gamma_mask(m1,m2,s1)
          enddo
       enddo
    enddo

    ! mu
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m1=1,nmodes
             do m2=m1,nmodes
                nmchk=0
                stachk=0
                nmchk(m1)=nmchk(m1)+1
                nmchk(m2)=nmchk(m2)+1
                stachk(s1)=stachk(s1)+1
                stachk(s2)=stachk(s2)+1
                mu_mask(m1,m2,s1,s2)=integralsym(nmchk,stachk,dipchk)
                mu_mask(m2,m1,s1,s2)=mu_mask(m1,m2,s1,s2)
                mu_mask(m1,m2,s1,s2)=mu_mask(m1,m2,s1,s2)
                mu_mask(m2,m1,s2,s1)=mu_mask(m1,m2,s1,s2)
             enddo
          enddo
       enddo
    enddo

    ! iota
    do s1=1,nsta
       do m1=1,nmodes
          nmchk=0
          stachk=0
          nmchk(m1)=3
          stachk(s1)=2
          iota_mask(m1,s1)=integralsym(nmchk,stachk,dipchk)
       enddo
    enddo
    
    ! tau
    do s1=1,nsta
       do s2=s1+1,nsta
          do m1=1,nmodes
             nmchk=0
             stachk=0
             nmchk(m1)=3
             stachk(s1)=stachk(s1)+1
             stachk(s2)=stachk(s2)+1
             tau_mask(m1,s1,s2)=integralsym(nmchk,stachk,dipchk)
             tau_mask(m1,s2,s1)=tau_mask(m1,s1,s2)
          enddo
       enddo
    enddo

    ! epsilon
    do s1=1,nsta
       do m1=1,nmodes
          nmchk=0
          stachk=0
          nmchk(m1)=4
          stachk(s1)=2
          epsilon_mask(m1,s1)=integralsym(nmchk,stachk,dipchk)
       enddo
    enddo

    ! xi
    do s1=1,nsta
       do s2=s1+1,nsta
          do m1=1,nmodes
             nmchk=0
             stachk=0
             nmchk(m1)=4
             stachk(s1)=stachk(s1)+1
             stachk(s2)=stachk(s2)+1
             xi_mask(m1,s1,s2)=integralsym(nmchk,stachk,dipchk)
             xi_mask(m1,s2,s1)=xi_mask(m1,s1,s2)
          enddo
       enddo
    enddo

    ! Diabatic dipole matrix expansion coefficients
    if (dipole) then

       ! 1st-order terms
       do c=1,3
          do s1=1,nsta
             do s2=s1,nsta
                do m1=1,nmodes
                   nmchk=0
                   stachk=0
                   dipchk=0
                   nmchk(m1)=1
                   stachk(s1)=stachk(s1)+1
                   stachk(s2)=stachk(s2)+1
                   dipchk(c)=1
                   dip1_mask(m1,s1,s2,c)=&
                        integralsym(nmchk,stachk,dipchk)
                   dip1_mask(m1,s2,s1,c)=dip1_mask(m1,s1,s2,c)
                enddo
             enddo
          enddo
       enddo

       ! 2nd-order terms
       do c=1,3
          do s1=1,nsta
             do s2=s1,nsta
                do m1=1,nmodes
                   do m2=m1,nmodes
                      nmchk=0
                      stachk=0
                      dipchk=0
                      nmchk(m1)=nmchk(m1)+1
                      nmchk(m2)=nmchk(m2)+1
                      stachk(s1)=stachk(s1)+1
                      stachk(s2)=stachk(s2)+1
                      dipchk(c)=1
                      dip2_mask(m1,m2,s1,s2,c)=&
                           integralsym(nmchk,stachk,dipchk)
                      dip2_mask(m1,m2,s2,s1,c)=dip2_mask(m1,m2,s1,s2,c)
                      dip2_mask(m2,m1,s2,s1,c)=dip2_mask(m1,m2,s1,s2,c)
                   enddo
               enddo
             enddo
          enddo
       enddo

       ! 3rd-order terms
       do c=1,3
          do s1=1,nsta
             do s2=s1,nsta
                do m1=1,nmodes
                   nmchk=0
                   stachk=0
                   dipchk=0
                   nmchk(m1)=3
                   stachk(s1)=stachk(s1)+1
                   stachk(s2)=stachk(s2)+1
                   dipchk(c)=1
                   dip3_mask(m1,s1,s2,c)=&
                        integralsym(nmchk,stachk,dipchk)
                   dip3_mask(m1,s2,s1,c)=dip3_mask(m1,s1,s2,c)
                enddo
             enddo
          enddo
       enddo

       ! 4th-order terms
       do c=1,3
          do s1=1,nsta
             do s2=s1,nsta
                do m1=1,nmodes
                   nmchk=0
                   stachk=0
                   dipchk=0
                   nmchk(m1)=4
                   stachk(s1)=stachk(s1)+1
                   stachk(s2)=stachk(s2)+1
                   dipchk(c)=1
                   dip4_mask(m1,s1,s2,c)=&
                        integralsym(nmchk,stachk,dipchk)
                   dip4_mask(m1,s2,s1,c)=dip3_mask(m1,s1,s2,c)
                enddo
             enddo
          enddo
       enddo
       
    endif
    
    return
    
  end subroutine create_mask
  
!######################################################################
  
  function integralsym (nmchk,stachk,dipchk) result(check)

    use constants
    use sysinfo
    
    implicit none

    integer                                :: i,j,k,l,r,numsym,prodnum
    integer, dimension(nmodes), intent(in) :: nmchk
    integer, dimension(nsta), intent(in)   :: stachk
    integer, dimension(3), intent(in)      :: dipchk
    integer, dimension(10)                 :: prodchar
    integer                                :: check
    real(dp)                               :: junk
    character(len=16)                      :: dummy

    dummy=''
    check=0
    prodnum=0

!----------------------------------------------------------------------
! N.B. In the matching of the mode/state/axis symmetry labels to those
! of the irreps of the point group, all letters are converted to 
! uppercase, and all instances of ` are converted to ' in order to
! allow for a free(ish) input format
!
! In the following,
!
! chdim: order of the point group
!
! numsym: no. of times that the totally symmetric irrep enters into
!         the decomposition of the direct product representation
!
! chel(i,j): character of the ith irrep of the jth symmetry operation
!----------------------------------------------------------------------    
    ! Contributions from differentiation wrt the normal modes
    do i=1,nmodes
       if (nmchk(i).gt.0) then            
          do j=1,nmchk(i)
             prodnum=prodnum+1
             do k=1,chdim
                if (nmlab(i).eq.irreplab(k)) prodchar(prodnum)=k
             enddo
          enddo
       endif
    enddo

    ! Contributions from the electronic wavefunctions
    do i=1,nsta
       if (stachk(i).gt.0) then            
          do j=1,stachk(i)
             prodnum=prodnum+1
             do k=1,chdim
                if (stalab(i).eq.irreplab(k)) prodchar(prodnum)=k
             enddo
          enddo
       endif
    enddo

    ! Contributions from the dipole operator
    do i=1,3
       if (dipchk(i).gt.0) then            
          do j=1,dipchk(i)
             prodnum=prodnum+1
             do k=1,chdim
                if (diplab(i).eq.irreplab(k)) prodchar(prodnum)=k
             enddo
          enddo
       endif
    enddo

    numsym=0
    ! Loop over operations
    do r=1,chdim
       ! Loop over irreps intering into direct product
       junk=1
       do i=1,prodnum
          junk=junk*chel(prodchar(i),r)
       enddo
       numsym=numsym+junk
    enddo
    
    numsym=numsym/chdim

    if (numsym.ne.0) check=1
    
    return

  end function integralsym

!######################################################################
!
! getcharacters determines the symmetry labels and characters of the 
! irreps of the point group of interest
!
!######################################################################

  subroutine getcharacters

    use iomod
    use sysinfo
    
    implicit none

    integer           :: i,j
    character(len=80) :: dummy

!----------------------------------------------------------------------
! BEFORE ADDING CHARACTERS FOR A NEW POINT GROUP, READ AND OBEY THE 
! FOLLOWING (unless you can think of a better way of doing it) 
!----------------------------------------------------------------------
!
! 1.) All pointgroup labels added here must not contain any lowercase 
!     letters.
!
! 2.) All irrep symmetry labels must not contain any lowercase letters,
!     and also must not contain any ` characters (only ' characters are
!     allowed).
!
! 3.) No " characters are allowed, only ''.
!
! 4.) All numbers must appear before instances of ' and ''
!
! chdim:     order of the point group
! chel(i,j): character of the ith irrep of the jth symmetry operation
!----------------------------------------------------------------------
    irreplab=''

    dummy=pntgrp 
    call uppercase(dummy)
    
    if (dummy.eq.'C1') then
       
       chdim=1
       
       irreplab(1)='A'
       
       chel(1,1)=1.0d0
       
    else if (dummy.eq.'C2V') then
       
       chdim=4
       
       irreplab(1)='A1'
       irreplab(2)='A2'
       irreplab(3)='B1'
       irreplab(4)='B2'
           
       chel(1,1)=1.0d0
       chel(1,2)=1.0d0
       chel(1,3)=1.0d0
       chel(1,4)=1.0d0
           
       chel(2,1)=1.0d0
       chel(2,2)=1.0d0
       chel(2,3)=-1.0d0
       chel(2,4)=-1.0d0
           
       chel(3,1)=1.0d0
       chel(3,2)=-1.0d0
       chel(3,3)=1.0d0
       chel(3,4)=-1.0d0

       chel(4,1)=1.0d0
       chel(4,2)=-1.0d0
       chel(4,3)=-1.0d0
       chel(4,4)=1.0d0
       
    else if (dummy.eq.'CS') then
       
       chdim=2
           
       irreplab(1)='A'''
       irreplab(2)='A'''''

       chel(1,1)=1.0d0
       chel(1,2)=1.0d0

       chel(2,1)=1.0d0
       chel(2,2)=-1.0d0

    else if (dummy.eq.'D2H') then

       chdim=8
       
       irreplab(1)='AG'
       irreplab(2)='B1G'
       irreplab(3)='B2G'
       irreplab(4)='B3G'
       irreplab(5)='AU'
       irreplab(6)='B1U'
       irreplab(7)='B2U'
       irreplab(8)='B3U'

       chel(1,1)=1.0d0
       chel(1,2)=1.0d0
       chel(1,3)=1.0d0
       chel(1,4)=1.0d0
       chel(1,5)=1.0d0
       chel(1,6)=1.0d0
       chel(1,7)=1.0d0
       chel(1,8)=1.0d0

       chel(2,1)=1.0d0
       chel(2,2)=1.0d0
       chel(2,3)=-1.0d0
       chel(2,4)=-1.0d0
       chel(2,5)=1.0d0
       chel(2,6)=1.0d0
       chel(2,7)=-1.0d0
       chel(2,8)=-1.0d0
       
       chel(3,1)=1.0d0
       chel(3,2)=-1.0d0
       chel(3,3)=1.0d0
       chel(3,4)=-1.0d0
       chel(3,5)=1.0d0
       chel(3,6)=-1.0d0
       chel(3,7)=1.0d0
       chel(3,8)=-1.0d0

       chel(4,1)=1.0d0
       chel(4,2)=-1.0d0
       chel(4,3)=-1.0d0
       chel(4,4)=1.0d0
       chel(4,5)=1.0d0
       chel(4,6)=-1.0d0
       chel(4,7)=-1.0d0
       chel(4,8)=1.0d0

       chel(5,1)=1.0d0
       chel(5,2)=1.0d0
       chel(5,3)=1.0d0
       chel(5,4)=1.0d0
       chel(5,5)=-1.0d0
       chel(5,6)=-1.0d0
       chel(5,7)=-1.0d0
       chel(5,8)=-1.0d0

       chel(6,1)=1.0d0
       chel(6,2)=1.0d0
       chel(6,3)=-1.0d0
       chel(6,4)=-1.0d0
       chel(6,5)=-1.0d0
       chel(6,6)=-1.0d0
       chel(6,7)=1.0d0
       chel(6,8)=1.0d0

       chel(7,1)=1.0d0
       chel(7,2)=-1.0d0
       chel(7,3)=1.0d0
       chel(7,4)=-1.0d0
       chel(7,5)=-1.0d0
       chel(7,6)=1.0d0
       chel(7,7)=-1.0d0
       chel(7,8)=1.0d0
       
       chel(8,1)=1.0d0
       chel(8,2)=-1.0d0
       chel(8,3)=-1.0d0
       chel(8,4)=1.0d0
       chel(8,5)=-1.0d0
       chel(8,6)=1.0d0
       chel(8,7)=1.0d0
       chel(8,8)=-1.0d0
           
    else if (dummy.eq.'CI') then
           
       chdim=2

       irreplab(1)='AG'
       irreplab(2)='AU'
       
       chel(1,1)=1.0d0
       chel(1,2)=1.0d0
       
       chel(2,1)=1.0d0
       chel(2,2)=-1.0d0

    else if (dummy.eq.'C2') then
       
       chdim=2

       irreplab(1)='A'
       irreplab(2)='B'
       
       chel(1,1)=1.0d0
       chel(1,2)=1.0d0
       
       chel(2,1)=1.0d0
       chel(2,2)=-1.0d0

    else if (dummy.eq.'D2') then
       
       chdim=4

       irreplab(1)='A'
       irreplab(2)='B1'
       irreplab(3)='B2'
       irreplab(4)='B3'
           
       chel(1,1)=1.0d0
       chel(1,2)=1.0d0
       chel(1,3)=1.0d0
       chel(1,4)=1.0d0
       
       chel(2,1)=1.0d0
       chel(2,2)=1.0d0
       chel(2,3)=-1.0d0
       chel(2,4)=-1.0d0
           
       chel(3,1)=1.0d0
       chel(3,2)=-1.0d0
       chel(3,3)=1.0d0
       chel(3,4)=-1.0d0

       chel(4,1)=1.0d0
       chel(4,2)=-1.0d0
       chel(4,3)=-1.0d0
       chel(4,4)=1.0d0

    else if (dummy.eq.'C2H') then

       chdim=4

       irreplab(1)='AG'
       irreplab(2)='BG'
       irreplab(3)='AU'
       irreplab(4)='BU'
           
       chel(1,1)=1.0d0
       chel(1,2)=1.0d0
       chel(1,3)=1.0d0
       chel(1,4)=1.0d0

       chel(2,1)=1.0d0
       chel(2,2)=-1.0d0
       chel(2,3)=1.0d0
       chel(2,4)=-1.0d0
           
       chel(3,1)=1.0d0
       chel(3,2)=1.0d0
       chel(3,3)=-1.0d0
       chel(3,4)=-1.0d0

       chel(4,1)=1.0d0
       chel(4,2)=-1.0d0
       chel(4,3)=-1.0d0
       chel(4,4)=1.0d0
       
    else

       errmsg='The point group'//trim(pntgrp)//' is not supported.'
       call error_control

    endif
        
    return

  end subroutine getcharacters
      
!######################################################################

  subroutine cleanirreps(string,dim)

    implicit none

    integer            :: dim,j,i
    character(len=dim) :: string,tmp
    logical(kind=4)    :: check

    tmp=''
    check=.false.

    ! Convert to upper case letters 
    call uppercase(string)

    ! Convert all instances of ` to '
    do j=1,len_trim(string)
       if (string(j:j).eq.'`') string(j:j)=''''
    enddo

    ! Convert all instances of " to ''
    do j=1,len_trim(string)
       if (string(j:j).eq.'"') then
          check=.true.
          string(j:j)=''''
          do i=1,j
             tmp(i:i)=string(i:i)
          enddo
          do i=j+1,len_trim(string)+1
             tmp(i:i)=string(i-1:i-1)
          enddo
       endif
    enddo
    if (check) then 
       string=tmp
       tmp=''
    endif
    
    ! Place all numbers before primes
    do j=1,len_trim(string)
       if (string(j:j).ge.'1'.and.string(j:j).le.'5') then
          if (string(j-1:j-1).eq.'''') then
             ! Double prime
             if (string(j-2:j-2).eq.'''') then
                tmp(1:1)=string(j:j)
                tmp(2:2)=string(j-1:j-1)
                tmp(3:3)=string(j-2:j-2)
                string(j-2:j-2)=tmp(1:1)
                string(j-1:j-1)=tmp(2:2)
                string(j:j)=tmp(3:3)
                ! Single prime
             else
                tmp(1:1)=string(j:j)
                tmp(2:2)=string(j-1:j-1)
                string(j-1:j-1)=tmp(1:1)
                string(j:j)=tmp(2:2)              
             endif
          endif
       endif
    enddo

    return
    
  end subroutine cleanirreps

!######################################################################
! 
!   uppercase: replaces all lowercase letters in a given character
!              string with the corresponding uppercase letters
!
!######################################################################

  subroutine uppercase(string)
      
    implicit none
      
    integer*8    ::  i,j,length
    character(*) :: string
      
    length=len(string)

    do j = 1,length
       if(string(j:j).ge."a".and.string(j:j).le."z")&
            string(j:j)=achar(iachar(string(j:j))-32)
    enddo
    
    return

  end subroutine uppercase

!######################################################################

  subroutine getnpar(dipolein)

    use constants
    use sysinfo
    
    implicit none

    integer           :: m,m1,m2,s,s1,s2,c
    logical, optional :: dipolein
    logical           :: dipole
    
!----------------------------------------------------------------------
! Optional fitting of diabatic dipole matrix elements
!----------------------------------------------------------------------
    dipole=.false.
    if (present(dipolein)) then
       dipole=dipolein
    endif

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    nkappa=0
    nlambda=0
    ngamma=0
    nmu=0
    niota=0
    ntau=0
    nepsilon=0
    nxi=0
    ndip1=0
    ndip2=0
    ndip3=0
    ndip4=0
    
!----------------------------------------------------------------------
! Coefficients of the vibronic coupling Hamiltonian
!----------------------------------------------------------------------
    ! kappa
    do s=1,nsta
       do m=1,nmodes
          ! Total number
          nkappa(1)=nkappa(1)+1
          ! Symmetry allowed
          if (kappa_mask(m,s).eq.1) nkappa(2)=nkappa(2)+1
       enddo
    enddo

    ! lambda
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m=1,nmodes
             ! Total number
             nlambda(1)=nlambda(1)+1
             ! Symmetry allowed
             if (lambda_mask(m,s1,s2).eq.1) nlambda(2)=nlambda(2)+1
          enddo
       enddo
    enddo

    ! gamma
    do s=1,nsta
       do m1=1,nmodes
          do m2=m1,nmodes
             ! Total number
             ngamma(1)=ngamma(1)+1
             ! Symmetry allowed
             if (gamma_mask(m1,m2,s).eq.1) ngamma(2)=ngamma(2)+1
          enddo
       enddo
    enddo

    ! mu
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m1=1,nmodes
             do m2=m1,nmodes
                ! Total number
                nmu(1)=nmu(1)+1
                ! Symmetry allowed
                if (mu_mask(m1,m2,s1,s2).eq.1) nmu(2)=nmu(2)+1
             enddo
          enddo
       enddo
    enddo

    ! iota
    do s1=1,nsta-1
       do m1=1,nmodes
          ! Total number
          niota(1)=niota(1)+1
          ! Symmetry allowed
          if (iota_mask(m1,s1).eq.1) niota(2)=niota(2)+1
       enddo
    enddo

    ! tau
    do s1=1,nsta
       do s2=s1+1,nsta
          do m1=1,nmodes
             ! Total number
             ntau(1)=ntau(1)+1
             ! Symmetry allowed
             if (tau_mask(m1,s1,s2).eq.1) ntau(2)=ntau(2)+1
          enddo
       enddo
    enddo

    ! epsilon
    do s1=1,nsta-1
       do m1=1,nmodes
          ! Total number
          nepsilon(1)=nepsilon(1)+1
          ! Symmetry allowed
          if (epsilon_mask(m1,s1).eq.1) nepsilon(2)=nepsilon(2)+1
       enddo
    enddo

    ! xi
    do s1=1,nsta
       do s2=s1+1,nsta
          do m1=1,nmodes
             ! Total number
             nxi(1)=nxi(1)+1
             ! Symmetry allowed
             if (xi_mask(m1,s1,s2).eq.1) nxi(2)=nxi(2)+1
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! Coefficients of the dipole matrix expansion
!----------------------------------------------------------------------
    if (dipole) then

       ! 1st-order terms
       do c=1,3
          do s1=1,nsta
             do s2=s1,nsta
                do m1=1,nmodes
                   ! Total number
                   ndip1(1)=ndip1(1)+1
                   ! Symmetry allowed
                   if (dip1_mask(m1,s1,s2,c).eq.1) ndip1(2)=ndip1(2)+1
                enddo
             enddo
          enddo
       enddo

       ! 2nd-order terms
       do c=1,3
          do s1=1,nsta
             do s2=s1,nsta
                do m1=1,nmodes
                   do m2=m1,nmodes
                      ! Total number
                      ndip2(1)=ndip2(1)+1
                      ! Symmetry allowed
                      if (dip2_mask(m1,m2,s1,s2,c).eq.1) &
                           ndip2(2)=ndip2(2)+1
                   enddo
                enddo
             enddo
          enddo
       enddo

       ! 3rd-order terms
       do c=1,3
          do s1=1,nsta
             do s2=s1,nsta
                do m1=1,nmodes
                   ! Total number
                   ndip3(1)=ndip3(1)+1
                   ! Symmetry allowed
                   if (dip3_mask(m1,s1,s2,c).eq.1) ndip3(2)=ndip3(2)+1
                enddo
             enddo
          enddo
       enddo
       
       ! 4th-order terms
       do c=1,3
          do s1=1,nsta
             do s2=s1,nsta
                do m1=1,nmodes
                   ! Total number
                   ndip4(1)=ndip4(1)+1
                   ! Symmetry allowed
                   if (dip4_mask(m1,s1,s2,c).eq.1) ndip4(2)=ndip4(2)+1
                enddo
             enddo
          enddo
       enddo
       
    endif

!----------------------------------------------------------------------
! Total
!----------------------------------------------------------------------
    ntot=nkappa+nlambda+ngamma+nmu+niota+ntau+nepsilon+nxi &
         +ndip1+ndip2+ndip3+ndip4
    
    return
    
  end subroutine getnpar

!######################################################################

  subroutine get_cutmask

    use constants
    use sysinfo
    
    implicit none

    integer :: m1,m2,s,s1,s2
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(cut_mask(nmodes,nmodes))
    cut_mask=0

!----------------------------------------------------------------------
! Determine which 2D cuts to make based on the symmetry-allowed
! coupling coefficients
!----------------------------------------------------------------------
    ! gamma
    do s=1,nsta
       do m1=1,nmodes-1
          do m2=m1+1,nmodes
             if (gamma_mask(m1,m2,s).eq.1) then
                cut_mask(m1,m2)=1
                cut_mask(m2,m1)=1
             endif
          enddo
       enddo
    enddo
    
    ! mu
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m1=1,nmodes-1
             do m2=m1+1,nmodes
                if (mu_mask(m1,m2,s1,s2).eq.1) then
                   cut_mask(m1,m2)=1
                   cut_mask(m2,m1)=1
                endif
             enddo
          enddo
       enddo
    enddo

    return
    
  end subroutine get_cutmask
  
!######################################################################
  
end module symmetry
