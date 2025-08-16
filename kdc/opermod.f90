!######################################################################
! opermod: routines for the writing of MCTDH operator files
!######################################################################

module opermod

  implicit none

contains

!######################################################################

  subroutine wroper_mctdh

    use constants
    use channels
    use iomod
    use sysinfo
    use parameters
    use symmetry
    use kdcglobal
    
    implicit none

    integer                        :: unit,n,m,m1,m2,s,s1,s2,i,j,k,c,&
                                      nl,ncurr,ndof,fel
    integer                        :: nzeta
    integer                        :: nzdip0
    integer                        :: nzdip1
    integer                        :: nzdip2
    integer                        :: nzdip3
    integer                        :: nzdip4
    real(dp), parameter            :: thrsh=5e-4_dp
    integer                        :: fac
    character(len=3)               :: an,am,am1,am2,as,as1,as2,afel,aft
    character(len=5)               :: aunit
    character(len=90)              :: string
    character(len=3)               :: amode
    character(len=8)               :: atmp
    character(len=1), dimension(3) :: acomp

!----------------------------------------------------------------------
! Determine the no. non-zero coupling coefficients in each class
!----------------------------------------------------------------------
    call get_nzpar(nzeta,nzdip0,nzdip1,nzdip2,nzdip3,nzdip4,thrsh)
    
!----------------------------------------------------------------------
! Unit label
!----------------------------------------------------------------------
    aunit=' , ev'

!----------------------------------------------------------------------
! Index of the electronic DOF
!----------------------------------------------------------------------
    fel=nmodes+1
    write(afel,'(i3)') fel

!----------------------------------------------------------------------
! Index of the Time DOF
!----------------------------------------------------------------------
    write(aft,'(i3)') fel+1
    
!----------------------------------------------------------------------
! Cartesian axis labels
!----------------------------------------------------------------------
    acomp=(/ 'x' , 'y' , 'z' /)
    
!----------------------------------------------------------------------
! Write the op_define section
!----------------------------------------------------------------------
    write(iop,'(a)') 'op_define-section'
    write(iop,'(a)') 'title'
    write(iop,'(a)') 'MCTDH operator file created by KDC'
    write(iop,'(a)') 'end-title'
    write(iop,'(a)') 'end-op_define-section'

!----------------------------------------------------------------------
! Write the parameter section
!----------------------------------------------------------------------
    ! Starting line
    write(iop,'(/,a)') 'parameter-section'

    ! Frequencies
    write(iop,'(/,a)') '# Frequencies'
    do m=1,nmodes
       write(am,'(i3)') m
       write(iop,'(a,F9.6,a)') &
            'omega_'//adjustl(am)//' =',freq(m),aunit
    enddo

    ! Energies
    write(iop,'(/,a)') '# Energies'
    do s=1,nsta
       write(as,'(i3)') s
       write(iop,'(a,F9.6,a)') &
            'E'//adjustl(as)//' = ',e0(s),aunit
    enddo

    ! On-diagonal one-mode coupling coefficients
    do n=1,order1
       write(iop,'(/,a,i0,x,a)') '# Order-',n,&
            'on-diagonal one-mode coupling coefficients'
       do s=1,nsta
          do m=1,nmodes
             if (coeff1_mask(m,s,s,n) == 0) cycle
             if (abs(coeff1(m,s,s,n)) < thrsh) cycle
             write(iop,'(a,4(i0,a),F9.6,a)') &
                  'tau',n,'_',m,'_',s,'_',s,' = ',&
                  coeff1(m,s,s,n),aunit
          enddo
       enddo
    enddo

    ! Off-diagonal one-mode coupling coefficients
    do n=1,order1
       write(iop,'(/,a,i0,x,a)') '# Order-',n,&
            'off-diagonal one-mode coupling coefficients'
       do s2=1,nsta-1
          do s1=s2+1,nsta
             do m=1,nmodes
                if (coeff1_mask(m,s1,s2,n) == 0) cycle
                if (abs(coeff1(m,s1,s2,n)) < thrsh) cycle
                write(iop,'(a,4(i0,a),F9.6,a)') &
                     'tau',n,'_',m,'_',s2,'_',s1,' = ',&
                     coeff1(m,s1,s2,n),aunit
             enddo
          enddo
       enddo
    enddo

    ! Two-mode coupling coefficients
    if (nzeta > 0) then
       
       ! On-diagonal 
       write(iop,'(/,a)') &
            '# Order-2 on-diagonal two-mode coupling coefficients'
       do s=1,nsta
          do m2=1,nmodes-1
             do m1=m2+1,nmodes
                if (coeff2_mask(m1,m2,s,s) == 0) cycle
                if (abs(coeff2(m1,m2,s,s)) < thrsh) cycle
                write(iop,'(a,4(i0,a),F9.6,a)') &
                     'eta_',m2,'_',m1,'_',s,'_',s,' = ',&
                     coeff2(m1,m2,s,s),aunit
             enddo
          enddo
       enddo

       ! Off-diagonal 
       write(iop,'(/,a)') &
            '# Order-2 off-diagonal two-mode coupling coefficients'
       do s2=1,nsta-1
          do s1=s2+1,nsta
             do m2=1,nmodes-1
                do m1=m2+1,nmodes
                   if (coeff2_mask(m1,m2,s1,s2) == 0) cycle
                   if (abs(coeff2(m1,m2,s1,s2)) < thrsh) cycle
                   write(iop,'(a,4(i0,a),F9.6,a)') &
                        'eta_',m2,'_',m1,'_',s2,'_',s1,' = ',&
                        coeff2(m1,m2,s1,s2),aunit
                enddo
             enddo
          enddo
       enddo
       
    endif
       
    ! Dipole matrix elements
    if (ldipfit) then

       ! Zeroth-order terms
       if (nzdip0.gt.0) then
          write(iop,'(/,a)') '# Dipole matrix elements: 0th-order &
               terms'
          do c=1,3
             do s1=1,nsta
                write(as1,'(i3)') s1
                do s2=s1,nsta
                   write(as2,'(i3)') s2
                   if (abs(dip0(s1,s2,c)).lt.thrsh) cycle
                   write(iop,'(a,F10.7)') 'dip0'//acomp(c)&
                        //trim(adjustl(as1))//'_'//trim(adjustl(as2))&
                        //' = ',dip0(s1,s2,c)
                enddo
             enddo
          enddo
       endif
          
       ! First-order terms
       if (nzdip1.gt.0) then
          write(iop,'(/,a)') '# Dipole matrix elements: 1st-order &
               terms'
          do c=1,3
             do s1=1,nsta
                write(as1,'(i3)') s1
                do s2=s1,nsta
                   write(as2,'(i3)') s2
                   do m=1,nmodes
                      if (abs(dip1(m,s1,s2,c)).lt.thrsh) cycle
                      if (dip1_mask(m,s1,s2,c).eq.0) cycle
                      write(am,'(i3)') m
                      write(iop,'(a,F10.7)') 'dip1'//acomp(c)&
                           //trim(adjustl(as1))//'_'&
                           //trim(adjustl(as2))//'_'&
                           //trim(adjustl(am))//' = ',dip1(m,s1,s2,c)
                   enddo
                enddo
             enddo
          enddo
       endif
          
       ! Second-order terms
       if (nzdip2.gt.0) then
          write(iop,'(/,a)') '# Dipole matrix elements: 2nd-order &
               terms'
          do c=1,3
             do s1=1,nsta
                write(as1,'(i3)') s1
                do s2=s1,nsta
                   write(as2,'(i3)') s2
                   do m1=1,nmodes
                      do m2=m1,nmodes
                         if (abs(dip2(m1,m2,s1,s2,c)).lt.thrsh) cycle
                         if (dip2_mask(m1,m2,s1,s2,c).eq.0) cycle
                         write(am1,'(i3)') m1
                         write(am2,'(i3)') m2
                         write(iop,'(a,F10.7)') 'dip2'//acomp(c)&
                              //trim(adjustl(as1))//'_'&
                              //trim(adjustl(as2))//'_'&
                              //trim(adjustl(am1))//'_'&
                              //trim(adjustl(am2))&
                              //' = ',dip2(m1,m2,s1,s2,c)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       endif

       ! Third-order terms
       if (nzdip3.gt.0) then
          write(iop,'(/,a)') '# Dipole matrix elements: 3rd-order &
               terms'
          do c=1,3
             do s1=1,nsta
                write(as1,'(i3)') s1
                do s2=s1,nsta
                   write(as2,'(i3)') s2
                   do m=1,nmodes
                      if (abs(dip3(m,s1,s2,c)).lt.thrsh) cycle
                      if (dip3_mask(m,s1,s2,c).eq.0) cycle
                      write(am,'(i3)') m
                      write(iop,'(a,F10.7)') 'dip3'//acomp(c)&
                           //trim(adjustl(as1))//'_'&
                           //trim(adjustl(as2))//'_'&
                           //trim(adjustl(am))//' = ',dip3(m,s1,s2,c)
                   enddo
                enddo
             enddo
          enddo
       endif

       ! Fourth-order terms
       if (nzdip4.gt.0) then
          write(iop,'(/,a)') '# Dipole matrix elements: 4th-order &
               terms'
          do c=1,3
             do s1=1,nsta
                write(as1,'(i3)') s1
                do s2=s1,nsta
                   write(as2,'(i3)') s2
                   do m=1,nmodes
                      if (abs(dip4(m,s1,s2,c)).lt.thrsh) cycle
                      if (dip4_mask(m,s1,s2,c).eq.0) cycle
                      write(am,'(i3)') m
                      write(iop,'(a,F10.7)') 'dip4'//acomp(c)&
                           //trim(adjustl(as1))//'_'&
                           //trim(adjustl(as2))//'_'&
                           //trim(adjustl(am))//' = ',dip4(m,s1,s2,c)
                   enddo
                enddo
             enddo
          enddo
       endif

       ! Pulse parameters
       write(iop,'(/,a)') '# Pulse parameters'
       write(iop,'(a)') 'A = 2.7726'
       write(iop,'(a)') 'B = A/PI'
       write(iop,'(a)') 'C = B^0.5'
       write(iop,'(a)') 'width = 30 , fs'
       write(iop,'(a)') 'freq = 5.0 , ev'
       write(iop,'(a)') 't0 = 0.0 , fs'
       write(iop,'(a,F20.3,a)') 'I0 = 0.00014286'
       write(iop,'(a)') 's = I0*width/C'
       
    endif
    
    ! Finishing line
    write(iop,'(/,a)') 'end-parameter-section'

!----------------------------------------------------------------------
! Write the labels sections
!----------------------------------------------------------------------
    if (ldipfit) then
       
       ! Starting line
       write(iop,'(/,a)') 'LABELS-SECTION'

       ! Pulse functions
       write(iop,'(/,a)') 'pulse = gauss[A/width^2,t0]'
       write(iop,'(a)') 'cosom = cos[freq,t0]'
       write(iop,'(a)') 'stepf = step[t0-1.25*width]'
       write(iop,'(a)') 'stepr = rstep[t0+1.25*width]'

       ! Finishing line
       write(iop,'(/,a)') 'end-labels-section'

    endif
    
!----------------------------------------------------------------------
! Write the Hamiltonian section
!----------------------------------------------------------------------
    ! Starting line
    write(iop,'(/,a)') 'hamiltonian-section'

    ! Modes section
    write(iop,'(/,38a)') ('-',i=1,38)
    ndof=nmodes+1
    if (ldipfit) ndof=ndof+1
    m=0
    nl=ceiling((real(ndof))/10.0d0)
    do i=1,nl
       ncurr=min(10,ndof-10*(i-1))
       string='modes|'
       do k=1,ncurr
          m=m+1
          if (m.lt.nmodes+1) then
             write(amode,'(i3)') m
             write(atmp,'(a)') ' Q'//adjustl(amode)//' '
             if (m.lt.10) then
                string=trim(string)//trim(atmp)//'  |'
             else
                string=trim(string)//trim(atmp)//' |'
             endif
          else
             if (m.eq.nmodes+1) then
                string=trim(string)//' el'
             else
                string=trim(string)//'  | Time'
             endif
          endif
       enddo
       write(iop,'(a)') trim(string)
    enddo
    write(iop,'(38a)') ('-',i=1,38)

    ! Kinetic energy operator
    write(iop,'(/,a)') '# Kinetic energy'
    do m=1,nmodes
       write(am,'(i3)') m
       write(iop,'(a)') &
            'omega_'//adjustl(am)//'  |'//adjustl(am)//'  KE'
    enddo

    ! Zeroth-order potential: VEEs
    write(iop,'(/,a)') '# Zeroth-order potential: VEEs'
    do s=1,nsta
       write(as,'(i3)') s
       write(iop,'(a)') &
            'E'//adjustl(as)//'  |'//adjustl(afel)&
            //'  S'//trim(adjustl(as))//'&'//trim(adjustl(as))
    enddo

    ! Zeroth-order potential: Harmonic oscillators
    write(iop,'(/,a)') '# Zeroth-order potential: &
         Harmonic oscillators'
    do m=1,nmodes
       write(am,'(i3)') m
       write(iop,'(a)') &
            '0.5*omega_'//adjustl(am)//'  |'//adjustl(am)//'  q^2'
    enddo

    ! On-diagonal one-mode coupling coefficients
    fac=1
    do n=1,order1
       fac=n*fac
       write(iop,'(/,a,i0,x,a)') '# Order-',n,&
            'on-diagonal one-mode coupling coefficients'
       do s=1,nsta
          do m=1,nmodes
             if (coeff1_mask(m,s,s,n) == 0) cycle
             if (abs(coeff1(m,s,s,n)) < thrsh) cycle
             write(iop,'(a,i0,a,8(a,i0))') '1.0/',fac,'.0',&
                  '*tau',n,'_',m,'_',s,'_',s,&
                  '  |',m,'  q^',n,&
                  '  |'//adjustl(afel)//'  S',s,&
                  '&',s
          enddo
       enddo
    enddo

    ! Off-diagonal one-mode coupling coefficients
    fac=1
    do n=1,order1
       fac=n*fac
       write(iop,'(/,a,i0,x,a)') '# Order-',n,&
            'off-diagonal one-mode coupling coefficients'
       do s2=1,nsta-1
          do s1=s2+1,nsta
             do m=1,nmodes
                if (coeff1_mask(m,s1,s2,n) == 0) cycle
                if (abs(coeff1(m,s1,s2,n)) < thrsh) cycle
                write(iop,'(a,i0,a,8(a,i0))') '1.0/',fac,'.0',&
                     '*tau',n,'_',m,'_',s2,'_',s1,&
                     '  |',m,'  q^',n,&
                     '  |'//adjustl(afel)//'  S',s2,&
                     '&',s1
             enddo
          enddo
       enddo
    enddo

    ! Two-mode coupling coefficients
    if (nzeta > 0) then
    
       ! On-diagonal
       write(iop,'(/,a)') &
            '# Order-2 on-diagonal two-mode coupling coefficients'
       do s=1,nsta
          do m2=1,nmodes-1
             do m1=m2+1,nmodes
                if (coeff2_mask(m1,m2,s,s) == 0) cycle
                if (abs(coeff2(m1,m2,s,s)) < thrsh) cycle
                write(iop,'(8(a,i0))') &
                     'eta_',m2,&
                     '_',m1,&
                     '_',s,&
                     '_',s,&
                     '  |',m2,&
                     '  q  |',m1,&
                     '  q  |'//adjustl(afel)//'  S',s,&
                     '&',s
             enddo
          enddo
       enddo

       ! Off-diagonal
       write(iop,'(/,a)') &
            '# Order-2 off-diagonal two-mode coupling coefficients'
       do s2=1,nsta-1
          do s1=s2+1,nsta
             do m2=1,nmodes-1
                do m1=m2+1,nmodes
                   if (coeff2_mask(m1,m2,s1,s2) == 0) cycle
                   if (abs(coeff2(m1,m2,s1,s2)) < thrsh) cycle
                   write(iop,'(8(a,i0))') &
                        'eta_',m2,&
                        '_',m1,&
                        '_',s2,&
                        '_',s1,&
                        '  |',m2,&
                        '  q  |',m1,&
                        '  q  |'//adjustl(afel)//'  S',s2,&
                        '&',s1
                enddo
             enddo
          enddo
       enddo
          
    endif
          
    ! Light-molecule interaction terms
    if (ldipfit) then
       
       ! Zeroth-order terms
       if (nzdip0.gt.0) then
          write(iop,'(/,a)') '# Light-molecule interaction: &
               0th-order terms'
          do c=1,3
             do s1=1,nsta
                write(as1,'(i3)') s1
                do s2=s1,nsta
                   write(as2,'(i3)') s2
                   if (abs(dip0(s1,s2,c)).lt.thrsh) cycle
                   write(iop,'(a)') '-dip0'//acomp(c)&
                        //trim(adjustl(as1))//'_'//trim(adjustl(as2))&
                        //'*s*C/width'//'  |'//adjustl(afel)&
                        //'  S'//trim(adjustl(as1))//'&'&
                        //trim(adjustl(as2))//'  |'//adjustl(aft)&
                        //'  cosom*pulse*stepf*stepr'
                enddo
             enddo
          enddo
       endif

       ! First-order terms
       if (nzdip1.gt.0) then
          write(iop,'(/,a)') '# Light-molecule interaction: &
               1st-order terms'
          do c=1,3
             do m=1,nmodes
                write(am,'(i3)') m
                do s1=1,nsta
                   write(as1,'(i3)') s1
                   do s2=s1,nsta
                      write(as2,'(i3)') s2
                      if (abs(dip1(m,s1,s2,c)).lt.thrsh) cycle
                      if (dip1_mask(m,s1,s2,c).eq.0) cycle
                      write(iop,'(a,F10.7)') '-dip1'//acomp(c)&
                           //trim(adjustl(as1))//'_'&
                           //trim(adjustl(as2))//'_'&
                           //trim(adjustl(am))//'*s*C/width'&
                           //'  |'//adjustl(am)//'  q'&
                           //'  |'//adjustl(afel)&
                           //'  S'//trim(adjustl(as1))//'&'&
                           //trim(adjustl(as2))//'  |'//adjustl(aft)&
                           //'  cosom*pulse*stepf*stepr'
                   enddo
                enddo
             enddo
          enddo
       endif
          
       ! Second-order terms
       if (nzdip2.gt.0) then
          write(iop,'(/,a)') '# Light-molecule interaction: &
               2nd-order terms'
          ! Quadratic terms
          do c=1,3
             do m=1,nmodes
                write(am,'(i3)') m
                do s1=1,nsta
                   write(as1,'(i3)') s1
                   do s2=s1,nsta
                      write(as2,'(i3)') s2
                      if (abs(dip2(m,m,s1,s2,c)).lt.thrsh) cycle
                      if (dip2_mask(m,m,s1,s2,c).eq.0) cycle
                      write(iop,'(a,F10.7)') '-0.5*dip2'//acomp(c)&
                           //trim(adjustl(as1))//'_'&
                           //trim(adjustl(as2))//'_'&
                           //trim(adjustl(am))//'_'&
                           //trim(adjustl(am))//'*s*C/width'&
                           //'  |'//adjustl(am)//'  q^2'&
                           //'  |'//adjustl(afel)&
                           //'  S'//trim(adjustl(as1))//'&'&
                           //trim(adjustl(as2))//'  |'//adjustl(aft)&
                           //'  cosom*pulse*stepf*stepr'
                   enddo
                enddo
             enddo
          enddo
          ! Bi-linear terms
          do c=1,3
             do m1=1,nmodes
                write(am1,'(i3)') m1
                do m2=m1+1,nmodes
                   write(am1,'(i3)') m1
                   do s1=1,nsta
                      write(as1,'(i3)') s1
                      do s2=s1,nsta
                         write(as2,'(i3)') s2
                         if (abs(dip2(m1,m2,s1,s2,c)).lt.thrsh) cycle
                         if (dip2_mask(m1,m2,s1,s2,c).eq.0) cycle
                         write(iop,'(a,F10.7)') '-dip2'//acomp(c)&
                              //trim(adjustl(as1))//'_'&
                              //trim(adjustl(as2))//'_'&
                              //trim(adjustl(am1))//'_'&
                              //trim(adjustl(am2))//'*s*C/width'&
                              //'  |'//adjustl(am1)//'  q'&
                              //'  |'//adjustl(am2)//'  q'&
                              //'  |'//adjustl(afel)&
                              //'  S'//trim(adjustl(as1))//'&'&
                              //trim(adjustl(as2))//'  |'//adjustl(aft)&
                              //'  cosom*pulse*stepf*stepr'
                      enddo
                   enddo
                enddo
             enddo
          enddo
       endif
       
       ! Third-order terms
       if (nzdip3.gt.0) then
          write(iop,'(/,a)') '# Light-molecule interaction: &
               3rd-order terms'
          do c=1,3
             do m=1,nmodes
                write(am,'(i3)') m
                do s1=1,nsta
                   write(as1,'(i3)') s1
                   do s2=s1,nsta
                      write(as2,'(i3)') s2
                      if (abs(dip3(m,s1,s2,c)).lt.thrsh) cycle
                      if (dip3_mask(m,s1,s2,c).eq.0) cycle
                      write(iop,'(a,F10.7)') '-0.166667*dip3'//acomp(c)&
                           //trim(adjustl(as1))//'_'&
                           //trim(adjustl(as2))//'_'&
                           //trim(adjustl(am))//'*s*C/width'&
                           //'  |'//adjustl(am)//'  q^3'&
                           //'  |'//adjustl(afel)&
                           //'  S'//trim(adjustl(as1))//'&'&
                           //trim(adjustl(as2))//'  |'//adjustl(aft)&
                           //'  cosom*pulse*stepf*stepr'
                   enddo
                enddo
             enddo
          enddo
       endif

       ! Fourth-order terms
       if (nzdip4.gt.0) then
          write(iop,'(/,a)') '# Light-molecule interaction: &
               4th-order terms'
          do c=1,3
             do m=1,nmodes
                write(am,'(i3)') m
                do s1=1,nsta
                   write(as1,'(i3)') s1
                   do s2=s1,nsta
                      write(as2,'(i3)') s2
                      if (abs(dip4(m,s1,s2,c)).lt.thrsh) cycle
                      if (dip4_mask(m,s1,s2,c).eq.0) cycle
                      write(iop,'(a,F10.7)') '-0.041666*dip4'//acomp(c)&
                           //trim(adjustl(as1))//'_'&
                           //trim(adjustl(as2))//'_'&
                           //trim(adjustl(am))//'*s*C/width'&
                           //'  |'//adjustl(am)//'  q^4'&
                           //'  |'//adjustl(afel)&
                           //'  S'//trim(adjustl(as1))//'&'&
                           //trim(adjustl(as2))//'  |'//adjustl(aft)&
                           //'  cosom*pulse*stepf*stepr'
                   enddo
                enddo
             enddo
          enddo
       endif
       
    endif
    
    ! Finishing line
    write(iop,'(/,a)') 'end-hamiltonian-section'

!----------------------------------------------------------------------
! Vertical excitation operators
!----------------------------------------------------------------------
    if (nsta.gt.1) then
       do s=2,nsta

          write(as,'(i3)') s

          ! Starting line
          write(iop,'(/,a)') 'hamiltonian-section_excite'&
               //trim(adjustl(as))

          ! Modes
          write(iop,'(38a)') ('-',i=1,38)
          m=0    
          nl=ceiling((real(nmodes+1))/10.0d0)
          do i=1,nl
             ncurr=min(10,nmodes+1-10*(i-1))
             string='modes|'
             do k=1,ncurr
                m=m+1
                if (m.lt.nmodes+1) then
                   write(amode,'(i3)') m
                   write(atmp,'(a)') ' Q'//adjustl(amode)//' '
                   if (m.lt.10) then
                      string=trim(string)//trim(atmp)//'  |'
                   else
                      string=trim(string)//trim(atmp)//' |'
                   endif
                else
                   if (m.eq.nmodes+1) then
                      string=trim(string)//' el'
                   else
                      string=trim(string)//'  | Time'
                   endif
                endif
             enddo
             write(iop,'(a)') trim(string)
          enddo
          write(iop,'(38a)') ('-',i=1,38)

          ! Vertical excitation operator
          write(iop,'(a)') '1.0 |'//adjustl(afel)//'  S1'&
               //'&'//trim(adjustl(as))

          ! Finishing line
          write(iop,'(a)') 'end-hamiltonian-section'
          
       enddo
    endif
    
!----------------------------------------------------------------------
! End line
!----------------------------------------------------------------------
    write(iop,'(/,a)') 'end-operator'
    
!----------------------------------------------------------------------
! Close the MCTDH operator file
!----------------------------------------------------------------------
    close(iop)
    
    return
    
  end subroutine wroper_mctdh

!######################################################################

  subroutine wroper_multiqd

    use constants
    use channels
    use iomod
    use sysinfo
    use parameters
    use symmetry
    use kdcglobal

    integer                        :: unit,n,m,m1,m2,s,s1,s2,i,j,k,c,&
                                      nl,ncurr,ndof,fel
    integer                        :: nzeta
    integer                        :: nzdip0
    integer                        :: nzdip1
    integer                        :: nzdip2
    integer                        :: nzdip3
    integer                        :: nzdip4
    real(dp), parameter            :: thrsh=5e-4_dp
    integer                        :: fac
    character(len=3)               :: an,am,am1,am2,as,as1,as2
    character(len=5)               :: aunit
    character(len=9)               :: apre
    
!----------------------------------------------------------------------
! Determine the no. non-zero coupling coefficients in each class
!----------------------------------------------------------------------
    call get_nzpar(nzeta,nzdip0,nzdip1,nzdip2,nzdip3,nzdip4,thrsh)
    
!----------------------------------------------------------------------
! Unit label
!----------------------------------------------------------------------
    aunit=' ev'

!----------------------------------------------------------------------
! Write the parameters values
!----------------------------------------------------------------------
    ! Frequencies
    write(iop,'(a)') '# Frequencies'
    do m=1,nmodes
       write(am,'(i3)') m
       write(iop,'(a,F9.6,a)') &
            'omega_'//adjustl(am)//' =',freq(m),aunit
    enddo

    ! Energies
    write(iop,'(/,a)') '# Energies'
    do s=1,nsta
       write(as,'(i3)') s
       write(iop,'(a,F9.6,a)') &
            'E'//adjustl(as)//' = ',e0(s),aunit
    enddo

    ! On-diagonal one-mode coupling coefficients
    do n=1,order1
       write(iop,'(/,a,i0,x,a)') '# Order-',n,&
            'on-diagonal one-mode coupling coefficients'
       do s=1,nsta
          do m=1,nmodes
             if (coeff1_mask(m,s,s,n) == 0) cycle
             if (abs(coeff1(m,s,s,n)) < thrsh) cycle
             write(iop,'(a,4(i0,a),F9.6,a)') &
                  'tau',n,'_',m,'_',s,'_',s,' = ',&
                  coeff1(m,s,s,n),aunit
          enddo
       enddo
    enddo

    ! Off-diagonal one-mode coupling coefficients
    do n=1,order1
       write(iop,'(/,a,i0,x,a)') '# Order-',n,&
            'off-diagonal one-mode coupling coefficients'
       do s2=1,nsta-1
          do s1=s2+1,nsta
             do m=1,nmodes
                if (coeff1_mask(m,s1,s2,n) == 0) cycle
                if (abs(coeff1(m,s1,s2,n)) < thrsh) cycle
                write(iop,'(a,4(i0,a),F9.6,a)') &
                     'tau',n,'_',m,'_',s2,'_',s1,' = ',&
                     coeff1(m,s1,s2,n),aunit
             enddo
          enddo
       enddo
    enddo

    ! Two-mode coupling coefficients
    if (nzeta > 0) then
       
       ! On-diagonal 
       write(iop,'(/,a)') &
            '# Order-2 on-diagonal two-mode coupling coefficients'
       do s=1,nsta
          do m2=1,nmodes-1
             do m1=m2+1,nmodes
                if (coeff2_mask(m1,m2,s,s) == 0) cycle
                if (abs(coeff2(m1,m2,s,s)) < thrsh) cycle
                write(iop,'(a,4(i0,a),F9.6,a)') &
                     'eta_',m2,'_',m1,'_',s,'_',s,' = ',&
                     coeff2(m1,m2,s,s),aunit
             enddo
          enddo
       enddo

       ! Off-diagonal 
       write(iop,'(/,a)') &
            '# Order-2 off-diagonal two-mode coupling coefficients'
       do s2=1,nsta-1
          do s1=s2+1,nsta
             do m2=1,nmodes-1
                do m1=m2+1,nmodes
                   if (coeff2_mask(m1,m2,s1,s2) == 0) cycle
                   if (abs(coeff2(m1,m2,s1,s2)) < thrsh) cycle
                   write(iop,'(a,4(i0,a),F9.6,a)') &
                        'eta_',m2,'_',m1,'_',s2,'_',s1,' = ',&
                        coeff2(m1,m2,s1,s2),aunit
                enddo
             enddo
          enddo
       enddo
       
    endif

!----------------------------------------------------------------------
! Write the operator terms
!----------------------------------------------------------------------
    ! Kinetic energy operator
    write(iop,'(/,a)') '# Kinetic energy'
    do m=1,nmodes
       write(am,'(i3)') m
       write(iop,'(a)') &
            '-0.5*omega_'//adjustl(am)//' * d2q_'//adjustl(am) &
            //' @ sum |i><i|'
    enddo

    ! Zeroth-order potential: VEEs
    write(iop,'(/,a)') '# Zeroth-order potential: VEEs'
    do s=1,nsta
       write(as,'(i3)') s
       write(iop,'(a)') &
            'E'//adjustl(as)//' *'&
            //' |'//trim(adjustl(as))//'><'//trim(adjustl(as))//'|'
    enddo

    ! Zeroth-order potential: Harmonic oscillators
    write(iop,'(/,a)') '# Zeroth-order potential: &
         Harmonic oscillators'
    do m=1,nmodes
       write(am,'(i3)') m
       write(iop,'(a)') &
            '0.5*omega_'//adjustl(am)//' * q_'//&
            trim(adjustl(am))//'^2'//' @ sum |i><i|'
    enddo

    ! On-diagonal one-mode coupling coefficients
    fac=1
    do n=1,order1
       fac=n*fac
       write(iop,'(/,a,i0,x,a)') '# Order-',n,&
            'on-diagonal one-mode coupling coefficients'
       do s=1,nsta
          do m=1,nmodes
             if (coeff1_mask(m,s,s,n) == 0) cycle
             if (abs(coeff1(m,s,s,n)) < thrsh) cycle
             write(apre,'(F9.6)') 1.0/fac
             write(iop,'(a,8(a,i0),a)') apre,'*tau',&
                  n,'_',m,'_',s,'_',s,&
                  ' * q_',m,'^',n,&
                  ' @ |',s,'><',s,'|'
          enddo
       enddo
    enddo

    ! Off-diagonal one-mode coupling coefficients
    fac=1
    do n=1,order1
       fac=n*fac
       write(iop,'(/,a,i0,x,a)') '# Order-',n,&
            'off-diagonal one-mode coupling coefficients'
       do s2=1,nsta-1
          do s1=s2+1,nsta
             do m=1,nmodes
                if (coeff1_mask(m,s1,s2,n) == 0) cycle
                if (abs(coeff1(m,s1,s2,n)) < thrsh) cycle

                write(apre,'(F9.6)') 1.0/fac
                write(iop,'(a,8(a,i0),a)') apre,'*tau',&
                     n,'_',m,'_',s,'_',s,&
                     ' * q_',m,'^',n,&
                     ' @ |',s2,'><',s1,'| + hc'
                
             enddo
          enddo
       enddo
    enddo

    ! Two-mode coupling coefficients
    if (nzeta > 0) then
    
       ! On-diagonal
       write(iop,'(/,a)') &
            '# Order-2 on-diagonal two-mode coupling coefficients'
       do s=1,nsta
          do m2=1,nmodes-1
             do m1=m2+1,nmodes
                if (coeff2_mask(m1,m2,s,s) == 0) cycle
                if (abs(coeff2(m1,m2,s,s)) < thrsh) cycle
                write(iop,'(8(a,i0),a)') &
                     'eta_',m2,&
                     '_',m1,&
                     '_',s,&
                     '_',s,&
                     '  * q_',m2,&
                     ' @ q_',m1,&
                     ' @ |',s,&
                     '><',s,&
                     '|'
             enddo
          enddo
       enddo

       write(iop,'(/,a)') &
            '# Order-2 off-diagonal two-mode coupling coefficients'
       do s2=1,nsta-1
          do s1=s2+1,nsta
             do m2=1,nmodes-1
                do m1=m2+1,nmodes
                   if (coeff2_mask(m1,m2,s1,s2) == 0) cycle
                   if (abs(coeff2(m1,m2,s1,s2)) < thrsh) cycle
                   write(iop,'(8(a,i0),a)') &
                     'eta_',m2,&
                     '_',m1,&
                     '_',s,&
                     '_',s,&
                     '  * q_',m2,&
                     ' @ q_',m1,&
                     ' @ |',s2,&
                     '><',s1,&
                     '| + hc'
                enddo
             enddo
          enddo
       enddo
                   
    endif
       
    return
    
  end subroutine wroper_multiqd
  
!######################################################################

  subroutine get_nzpar(nzeta,nzdip0,nzdip1,nzdip2,nzdip3,nzdip4,thrsh)

    use constants
    use sysinfo
    use symmetry
    use parameters
    use kdcglobal
    
    implicit none

    integer              :: nzeta
    integer              :: nzdip0
    integer              :: nzdip1
    integer              :: nzdip2
    integer              :: nzdip3
    integer              :: nzdip4
    integer              :: m,s,m1,m2,s1,s2,c
    real(dp), intent(in) :: thrsh

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    nzeta=0
    nzdip0=0
    nzdip1=0
    nzdip2=0
    nzdip3=0
    nzdip4=0
    
!----------------------------------------------------------------------
! Coupling coefficients of the vibronic coupling Hamiltonian
!----------------------------------------------------------------------
    ! Two-mode terms
    do s2=1,nsta
       do s1=s2,nsta
          do m1=1,nmodes-1
             do m2=m1+1,nmodes
                if (coeff2_mask(m1,m2,s1,s2) == 0) cycle
                if (abs(coeff2(m1,m2,s1,s2)) < thrsh) cycle
                nzeta=nzeta+1
             enddo
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! Coefficients of the expansion of the diabatic dipole matrix
!----------------------------------------------------------------------
    if (ldipfit) then

       ! 0th-order terms
       do c=1,3
          do s1=1,nsta
             do s2=s1,nsta
                if (abs(dip0(s1,s2,c)).lt.thrsh) cycle
                nzdip0=nzdip0+1
             enddo
          enddo
       enddo
                
       ! 1st-order terms
       do c=1,3
          do s1=1,nsta
             do s2=s1,nsta
                do m=1,nmodes
                   if (dip1_mask(m,s1,s2,c).eq.0) cycle
                   if (abs(dip1(m,s1,s2,c)).lt.thrsh) cycle
                   nzdip1=nzdip1+1
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
                      if (dip2_mask(m1,m2,s1,s2,c).eq.0) cycle
                      if (abs(dip2(m1,m2,s1,s2,c)).lt.thrsh) cycle
                      nzdip2=nzdip2+1
                   enddo
                enddo
             enddo
          enddo
       enddo

       ! 3rd-order terms
       do c=1,3
          do s1=1,nsta
             do s2=s1,nsta
                do m=1,nmodes
                   if (dip3_mask(m,s1,s2,c).eq.0) cycle
                   if (abs(dip3(m,s1,s2,c)).lt.thrsh) cycle
                   nzdip3=nzdip3+1
                enddo
             enddo
          enddo
       enddo

       ! 4th-order terms
       do c=1,3
          do s1=1,nsta
             do s2=s1,nsta
                do m=1,nmodes
                   if (dip4_mask(m,s1,s2,c).eq.0) cycle
                   if (abs(dip4(m,s1,s2,c)).lt.thrsh) cycle
                   nzdip4=nzdip4+1
                enddo
             enddo
          enddo
       enddo
       
    endif
    
    return
    
  end subroutine get_nzpar
    
!######################################################################

end module opermod
