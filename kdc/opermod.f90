!######################################################################
! opermod: routines for the writing of MCTDH operator files
!######################################################################

module opermod

  implicit none

contains

!######################################################################

  subroutine wroper

    use constants
    use channels
    use iomod
    use sysinfo
    use parameters
    use symmetry
    use kdcglobal
    
    implicit none

    integer                        :: unit,m,m1,m2,s,s1,s2,i,j,k,c,nl,&
                                      ncurr,ndof,fel
    integer                        :: nzkappa
    integer                        :: nzlambda
    integer                        :: nzgamma1,nzgamma2
    integer                        :: nzmu1,nzmu2
    integer                        :: nziota
    integer                        :: nztau
    integer                        :: nzepsilon
    integer                        :: nzxi
    integer                        :: nzdip0
    integer                        :: nzdip1
    integer                        :: nzdip2
    integer                        :: nzdip3
    integer                        :: nzdip4
    real(dp), parameter            :: thrsh=5e-4_dp
    character(len=3)               :: am,am1,am2,as,as1,as2,afel,aft
    character(len=5)               :: aunit
    character(len=90)              :: string
    character(len=3)               :: amode
    character(len=8)               :: atmp
    character(len=1), dimension(3) :: acomp

    
    !! TEST
    !real(dp) :: Fij(ncoo)
    !real(dp) :: DeltaE
    !
    !! (1) Transform to the Cartesian coordinate basis
    !Fij = matmul(transpose(coonm), lambda(:,2,6))
    !
    !! (2) Convert to a.u.
    !Fij = Fij / eh2ev / ang2bohr
    !
    !! (3) Divide through by the energy difference
    !DeltaE = (e0(6)-e0(2)) / eh2ev
    !Fij = Fij / DeltaE
    !
    !do i=1,natm
    !   write(6,'(i0,3(2x,F12.7))') i,(Fij(j),j=i*3-2,i*3)
    !enddo
    !
    !stop
    !! TEST
    
    
!----------------------------------------------------------------------
! Determine the no. non-zero coupling coefficients in each class
!----------------------------------------------------------------------
    call get_nzpar(nzkappa,nzlambda,nzgamma1,nzgamma2,nzmu1,nzmu2,&
         nziota,nztau,nzepsilon,nzxi,nzdip0,nzdip1,nzdip2,nzdip3,&
         nzdip4,thrsh)
    
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

    ! 1st-order intrastate coupling constants (kappa)
    if (nzkappa.gt.0) then
       write(iop,'(/,a)') '# 1st-order intrastate coupling &
            constants (kappa)'
       do s=1,nsta
          write(as,'(i3)') s
          do m=1,nmodes
             if (kappa_mask(m,s).eq.0) cycle
             if (abs(kappa(m,s)).lt.thrsh) cycle
             write(am,'(i3)') m
             write(iop,'(a,F9.6,a)') 'kappa'//trim(adjustl(as))&
                  //'_'//adjustl(am)//' = ',kappa(m,s),aunit
          enddo
       enddo
    endif
       
    ! 1st-order intrastate coupling constants (lambda)
    if (nzlambda.gt.0) then
       write(iop,'(/,a)') '# 1st-order interstate coupling &
            constants (lambda)'
       do s1=1,nsta-1
          write(as1,'(i3)') s1
          do s2=s1+1,nsta
             write(as2,'(i3)') s2
             do m=1,nmodes
                if (lambda_mask(m,s1,s2).eq.0) cycle
                if (abs(lambda(m,s1,s2)).lt.thrsh) cycle
                write(am,'(i3)') m
                write(iop,'(a,F9.6,a)') 'lambda'//trim(adjustl(as1))&
                     //'_'//trim(adjustl(as2))//'_'//adjustl(am)//&
                     ' = ',lambda(m,s1,s2),aunit
             enddo
          enddo
       enddo
    endif

    ! Quadratic (aa) 2nd-order intrastate coupling constants (gamma)
    if (nzgamma1.gt.0) then
       write(iop,'(/,a)') '# Quadratic 2nd-order intrastate &
            coupling constants (gamma)'
       do s=1,nsta
          write(as,'(i3)') s
          do m=1,nmodes
             if (gamma_mask(m,m,s).eq.0) cycle
             if (abs(gamma(m,m,s)).lt.thrsh) cycle
             write(am,'(i3)') m
             write(iop,'(a,F9.6,a)') 'gamma'//trim(adjustl(as))&
                  //'_'//trim(adjustl(am))//'_'//adjustl(am)&
                  //' = ',gamma(m,m,s),aunit
          enddo
       enddo
    endif
       
    ! Bi-linear (ab) 2nd-order intrastate coupling constants (gamma)
    if (nzgamma2.gt.0) then
       write(iop,'(/,a)') '# Bi-linear 2nd-order intrastate &
            coupling constants (gamma)'
       do s=1,nsta
          write(as,'(i3)') s
          do m1=1,nmodes-1
             do m2=m1+1,nmodes
                if (gamma_mask(m1,m2,s).eq.0) cycle
                if (abs(gamma(m1,m2,s)).lt.thrsh) cycle
                write(am1,'(i3)') m1
                write(am2,'(i3)') m2
                write(iop,'(a,F9.6,a)') 'gamma'//trim(adjustl(as))&
                     //'_'//trim(adjustl(am1))//'_'//adjustl(am2)&
                     //' = ',gamma(m1,m2,s),aunit
             enddo
          enddo
       enddo
    endif
       
    ! Quadratic (aa) 2nd-order interstate coupling constants (mu)
    if (nzmu1.gt.0) then
       write(iop,'(/,a)') '# Quadratic 2nd-order interstate &
            coupling constants (mu)'
       do s1=1,nsta-1
          write(as1,'(i3)') s1
          do s2=s1+1,nsta
             write(as2,'(i3)') s2          
             do m=1,nmodes
                if (mu_mask(m,m,s1,s2).eq.0) cycle
                if (abs(mu(m,m,s1,s2)).lt.thrsh) cycle
                write(am,'(i3)') m
                write(iop,'(a,F9.6,a)') 'mu'//trim(adjustl(as1))&
                     //'_'//trim(adjustl(as2))//'_'//&
                     trim(adjustl(am))//'_'//adjustl(am)//' = ',&
                     mu(m,m,s1,s2),aunit
             enddo
          enddo
       enddo
    endif
       
    ! Bi-linear (ab) 2nd-order interstate coupling constants (mu)
    if (nzmu2.gt.0) then
       write(iop,'(/,a)') '# Bi-linear 2nd-order interstate &
            coupling constants (mu)'
       do s1=1,nsta-1
          write(as1,'(i3)') s1
          do s2=s1+1,nsta
             write(as2,'(i3)') s2          
             do m1=1,nmodes-1
                do m2=m1+1,nmodes
                   if (mu_mask(m1,m2,s1,s2).eq.0) cycle
                   if (abs(mu(m1,m2,s1,s2)).lt.thrsh) cycle
                   write(am1,'(i3)') m1
                   write(am2,'(i3)') m2
                   write(iop,'(a,F9.6,a)') 'mu'//trim(adjustl(as1))&
                        //'_'//trim(adjustl(as2))//'_'//&
                        trim(adjustl(am1))//'_'//adjustl(am2)//&
                        ' = ',mu(m1,m2,s1,s2),aunit
                enddo
             enddo
          enddo
       enddo
    endif
       
    ! Cubic (aaa) 3rd-order intrastate coupling constants (iota)
    if (nziota.gt.0) then
       write(iop,'(/,a)') '# Cubic 3rd-order intrastate coupling &
            constants (iota)'
       do s=1,nsta
          write(as,'(i3)') s
          do m=1,nmodes
             if (iota_mask(m,s).eq.0) cycle
             if (abs(iota(m,s)).lt.thrsh) cycle
             write(am,'(i3)') m
             write(iop,'(a,F9.6,a)') 'iota'//trim(adjustl(as))&
                  //'_'//trim(adjustl(am))//' = ',iota(m,s),aunit
          enddo
       enddo
    endif
       
    ! Cubic (aaa) 3rd-order interstate coupling constants (tau)
    if (nztau.gt.0) then
       write(iop,'(/,a)') '# Cubic 3rd-order interstate coupling &
            constants (tau)'
       do s1=1,nsta-1
          write(as1,'(i3)') s1
          do s2=s1+1,nsta
             write(as2,'(i3)') s2          
             do m=1,nmodes
                if (tau_mask(m,s1,s2).eq.0) cycle
                if (abs(tau(m,s1,s2)).lt.thrsh) cycle
                write(am,'(i3)') m
                write(iop,'(a,F9.6,a)') 'tau'//trim(adjustl(as1))&
                     //'_'//trim(adjustl(as2))//'_'//&
                     trim(adjustl(am))//' = ',&
                     tau(m,s1,s2),aunit
             enddo
          enddo
       enddo
    endif
       
    ! Quartic (aaaa) 4rd-order intrastate coupling constants (epsilon)
    if (nzepsilon.gt.0) then
       write(iop,'(/,a)') '# Cubic 4rd-order intrastate coupling &
            constants (epsilon)'
       do s=1,nsta
          write(as,'(i3)') s
          do m=1,nmodes
             if (epsilon_mask(m,s).eq.0) cycle
             if (abs(epsilon(m,s)).lt.thrsh) cycle
             write(am,'(i3)') m
             write(iop,'(a,F9.6,a)') 'epsilon'//trim(adjustl(as))&
                  //'_'//trim(adjustl(am))//' = ',epsilon(m,s),aunit
          enddo
       enddo
    endif
       
    ! Quartic (aaaa) 4rd-order interstate coupling constants (xi)
    if (nzxi.gt.0) then
       write(iop,'(/,a)') '# Cubic 3rd-order interstate coupling &
            constants (xi)'
       do s1=1,nsta-1
          write(as1,'(i3)') s1
          do s2=s1+1,nsta
             write(as2,'(i3)') s2          
             do m=1,nmodes
                if (xi_mask(m,s1,s2).eq.0) cycle
                if (abs(xi(m,s1,s2)).lt.thrsh) cycle
                write(am,'(i3)') m
                write(iop,'(a,F9.6,a)') 'xi'//trim(adjustl(as1))&
                     //'_'//trim(adjustl(as2))//'_'//&
                     trim(adjustl(am))//' = ',&
                     xi(m,s1,s2),aunit
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

    ! 1st-order intrastate coupling terms (kappa)
    if (nzkappa.gt.0) then
       write(iop,'(/,a)') '# 1st-order intrastate coupling &
            constants (kappa)'
       do m=1,nmodes
          write(am,'(i3)') m
          do s=1,nsta
             if (kappa_mask(m,s).eq.0) cycle
             if (abs(kappa(m,s)).lt.thrsh) cycle
             write(as,'(i3)') s
             write(iop,'(a)') 'kappa'//trim(adjustl(as))&
                  //'_'//adjustl(am)&
                  //'  |'//adjustl(am)//'  q'//'  |'//adjustl(afel)&
                  //'  S'//trim(adjustl(as))//'&'//trim(adjustl(as))
          enddo
       enddo
    endif
       
    ! 1st-order interstate coupling terms (lambda)
    if (nzlambda.gt.0) then
       write(iop,'(/,a)') '# 1st-order interstate coupling &
            constants (lambda)'
       do m=1,nmodes
          write(am,'(i3)') m
          do s1=1,nsta-1
             write(as1,'(i3)') s1
             do s2=s1+1,nsta
                if (lambda_mask(m,s1,s2).eq.0) cycle
                if (abs(lambda(m,s1,s2)).lt.thrsh) cycle
                write(as2,'(i3)') s2
                write(iop,'(a)') 'lambda'//trim(adjustl(as1)) &
                     //'_'//trim(adjustl(as2))//'_'//adjustl(am) &
                     //'  |'//adjustl(am)//'  q'//'  |'//adjustl(afel)&
                     //'  S'//trim(adjustl(as1))//'&'&
                     //trim(adjustl(as2))
             enddo
          enddo
       enddo
    endif
       
    ! Quadratic 2nd-order intrastate coupling constants (gamma)
    if (nzgamma1.gt.0) then
       write(iop,'(/,a)') '# Quadratic 2nd-order intrastate &
            coupling constants (gamma)'
       do m=1,nmodes
          write(am,'(i3)') m
          do s=1,nsta
             if (gamma_mask(m,m,s).eq.0) cycle
             if (abs(gamma(m,m,s)).lt.thrsh) cycle
             write(as,'(i3)') s
             write(iop,'(a)') '0.5*gamma'//trim(adjustl(as))&
                  //'_'//trim(adjustl(am))//'_'//adjustl(am)&
                  //'  |'//adjustl(am)//'  q^2'//'  |'//adjustl(afel)&
                  //'  S'//trim(adjustl(as))//'&'//trim(adjustl(as))
          enddo
       enddo
    endif
       
    ! Bi-linear 2nd-order intrastate coupling constants (gamma)
    if (nzgamma2.gt.0) then
       write(iop,'(/,a)') '# Bi-linear 2nd-order intrastate &
            coupling constants (gamma)'
       do m1=1,nmodes-1
          do m2=m1+1,nmodes
             do s=1,nsta
                if (gamma_mask(m1,m2,s).eq.0) cycle
                if (abs(gamma(m1,m2,s)).lt.thrsh) cycle
                write(am1,'(i3)') m1
                write(am2,'(i3)') m2
                write(as,'(i3)') s
                write(iop,'(a)') 'gamma'//trim(adjustl(as))&
                     //'_'//trim(adjustl(am1))//'_'//adjustl(am2)&
                     //'  |'//adjustl(am1)//'  q'&
                     //'  |'//adjustl(am2)//'  q'&
                     //'  |'//adjustl(afel)&
                     //'  S'//trim(adjustl(as))//'&'&
                     //trim(adjustl(as))
             enddo
          enddo
       enddo
    endif
       
    ! Quadratic 2nd-order interstate coupling constants (mu)
    if (nzmu1.gt.0) then
       write(iop,'(/,a)') '# Quadratic 2nd-order interstate &
            coupling constants (mu)'
       do m=1,nmodes
          write(am,'(i3)') m
          do s1=1,nsta-1
             do s2=s1+1,nsta
                if (mu_mask(m,m,s1,s2).eq.0) cycle
                if (abs(mu(m,m,s1,s2)).lt.thrsh) cycle
                write(as1,'(i3)') s1
                write(as2,'(i3)') s2
                write(iop,'(a)') '0.5*mu'//trim(adjustl(as1))&
                     //'_'//trim(adjustl(as2))&
                     //'_'//trim(adjustl(am))//'_'//adjustl(am)&
                     //'  |'//adjustl(am)//'  q^2'//'  |'&
                     //adjustl(afel)//'  S'//trim(adjustl(as1))&
                     //'&'//trim(adjustl(as2))
             enddo
          enddo
       enddo
    endif
       
    ! Bi-linear 2nd-order interstate coupling constants (mu)
    if (nzmu2.gt.0) then
       write(iop,'(/,a)') '# Bi-linear 2nd-order interstate &
            coupling constants (mu)'
       do m1=1,nmodes-1
          write(am1,'(i3)') m1
          do m2=m1+1,nmodes
             write(am2,'(i3)') m2
             do s1=1,nsta-1
                write(as1,'(i3)') s1
                do s2=s1+1,nsta
                   write(as2,'(i3)') s2
                   if (mu_mask(m1,m2,s1,s2).eq.0) cycle
                   if (abs(mu(m1,m2,s1,s2)).lt.thrsh) cycle
                   write(iop,'(a)') 'mu'//trim(adjustl(as1))&
                        //'_'//trim(adjustl(as2))&
                        //'_'//trim(adjustl(am1))//'_'//adjustl(am2)&
                        //'  |'//adjustl(am1)//'  q'&
                        //'  |'//adjustl(am2)//'  q'&
                        //'  |'//adjustl(afel)&
                        //'  S'//trim(adjustl(as1))//'&'&
                        //trim(adjustl(as2))
                enddo
             enddo
          enddo
       enddo
    endif

    ! Cubic (aaa) 3rd-order intrastate coupling constants (iota)
    if (nziota.gt.0) then
       write(iop,'(/,a)') '# Cubic 3rd-order intrastate coupling &
            constants (iota)'
       do m=1,nmodes
          write(am,'(i3)') m
          do s=1,nsta
             if (iota_mask(m,s).eq.0) cycle
             if (abs(iota(m,s)).lt.thrsh) cycle
             write(as,'(i3)') s
             write(iop,'(a)') '0.166667*iota'//trim(adjustl(as))&
                  //'_'//trim(adjustl(am))&
                  //'  |'//adjustl(am)//'  q^3'//'  |'//adjustl(afel)&
                  //'  S'//trim(adjustl(as))//'&'//trim(adjustl(as))
          enddo
       enddo
    endif

    ! Cubic (aaa) 3rd-order interstate coupling constants (tau)
    if (nztau.gt.0) then
       write(iop,'(/,a)') '# Cubic 3rd-order interstate coupling &
            constants (tau)'
       do m=1,nmodes
          write(am,'(i3)') m
          do s1=1,nsta-1
             do s2=s1+1,nsta
                if (tau_mask(m,s1,s2).eq.0) cycle
                if (abs(tau(m,s1,s2)).lt.thrsh) cycle
                write(as1,'(i3)') s1
                write(as2,'(i3)') s2
                write(iop,'(a)') '0.166667*tau'//trim(adjustl(as1))&
                     //'_'//trim(adjustl(as2))&
                     //'_'//trim(adjustl(am))&
                     //'  |'//adjustl(am)//'  q^3'//'  |'&
                     //adjustl(afel)//'  S'//trim(adjustl(as1))&
                     //'&'//trim(adjustl(as2))
             enddo
          enddo
       enddo
    endif

     ! Quartic (aaaa) 4rd-order intrastate coupling constants (epsilon)
    if (nzepsilon.gt.0) then
       write(iop,'(/,a)') '# Cubic 4rd-order intrastate coupling &
            constants (epsilon)'
       do m=1,nmodes
          write(am,'(i3)') m
          do s=1,nsta
             if (epsilon_mask(m,s).eq.0) cycle
             if (abs(epsilon(m,s)).lt.thrsh) cycle
             write(as,'(i3)') s
             write(iop,'(a)') '0.041666*epsilon'//trim(adjustl(as))&
                  //'_'//trim(adjustl(am))&
                  //'  |'//adjustl(am)//'  q^4'//'  |'//adjustl(afel)&
                  //'  S'//trim(adjustl(as))//'&'//trim(adjustl(as))
          enddo
       enddo
    endif

    ! Quartic (aaaa) 4rd-order interstate coupling constants (xi)
    if (nzxi.gt.0) then
       write(iop,'(/,a)') '# Cubic 3rd-order interstate coupling &
            constants (xi)'
       do m=1,nmodes
          write(am,'(i3)') m
          do s1=1,nsta-1
             do s2=s1+1,nsta
                if (xi_mask(m,s1,s2).eq.0) cycle
                if (abs(xi(m,s1,s2)).lt.thrsh) cycle
                write(as1,'(i3)') s1
                write(as2,'(i3)') s2
                write(iop,'(a)') '0.041666*xi'//trim(adjustl(as1))&
                     //'_'//trim(adjustl(as2))&
                     //'_'//trim(adjustl(am))&
                     //'  |'//adjustl(am)//'  q^4'//'  |'&
                     //adjustl(afel)//'  S'//trim(adjustl(as1))&
                     //'&'//trim(adjustl(as2))
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
    
  end subroutine wroper

!######################################################################

  subroutine get_nzpar(nzkappa,nzlambda,nzgamma1,nzgamma2,nzmu1,&
       nzmu2,nziota,nztau,nzepsilon,nzxi,nzdip0,nzdip1,nzdip2,nzdip3,&
       nzdip4,thrsh)

    use constants
    use sysinfo
    use symmetry
    use parameters
    use kdcglobal
    
    implicit none

    integer              :: nzkappa
    integer              :: nzlambda
    integer              :: nzgamma1,nzgamma2
    integer              :: nzmu1,nzmu2
    integer              :: nziota
    integer              :: nztau
    integer              :: nzepsilon
    integer              :: nzxi
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
    nzkappa=0
    nzlambda=0
    nzgamma1=0
    nzgamma2=0
    nzmu1=0
    nzmu2=0
    nziota=0
    nztau=0
    nzepsilon=0
    nzxi=0
    nzdip0=0
    nzdip1=0
    nzdip2=0
    nzdip3=0
    nzdip4=0
    
!----------------------------------------------------------------------
! Coupling coefficients of the vibronic coupling Hamiltonian
!----------------------------------------------------------------------
    ! 1st-order intrastate coupling constants (kappa)
    do s=1,nsta
       do m=1,nmodes
          if (kappa_mask(m,s).eq.0) cycle
          if (abs(kappa(m,s)).lt.thrsh) cycle
          nzkappa=nzkappa+1
       enddo
    enddo

    ! 1st-order intrastate coupling constants (lambda)
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m=1,nmodes
             if (lambda_mask(m,s1,s2).eq.0) cycle
             if (abs(lambda(m,s1,s2)).lt.thrsh) cycle
             nzlambda=nzlambda+1
          enddo
       enddo
    enddo

    ! Quadratic (aa) 2nd-order intrastate coupling constants (gamma)
    do s=1,nsta
       do m=1,nmodes
          if (gamma_mask(m,m,s).eq.0) cycle
          if (abs(gamma(m,m,s)).lt.thrsh) cycle
          nzgamma1=nzgamma1+1
       enddo
    enddo

    ! Bi-linear (ab) 2nd-order intrastate coupling constants (gamma)
    do s=1,nsta
       do m1=1,nmodes-1
          do m2=m1+1,nmodes
             if (gamma_mask(m1,m2,s).eq.0) cycle
             if (abs(gamma(m1,m2,s)).lt.thrsh) cycle
             nzgamma2=nzgamma2+1
          enddo
       enddo
    enddo

    ! Quadratic (aa) 2nd-order interstate coupling constants (mu)
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m=1,nmodes
             if (mu_mask(m,m,s1,s2).eq.0) cycle
             if (abs(mu(m,m,s1,s2)).lt.thrsh) cycle
             nzmu1=nzmu1+1
          enddo
       enddo
    enddo

    ! Bi-linear (ab) 2nd-order interstate coupling constants (mu)
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m1=1,nmodes-1
             do m2=m1+1,nmodes
                if (mu_mask(m1,m2,s1,s2).eq.0) cycle
                if (abs(mu(m1,m2,s1,s2)).lt.thrsh) cycle
                nzmu2=nzmu2+1
             enddo
          enddo
       enddo
    enddo

    ! Cubic (aaa) 3rd-order intrastate coupling constants (iota)
    do s=1,nsta
       do m=1,nmodes
          if (iota_mask(m,s).eq.0) cycle
          if (abs(iota(m,s)).lt.thrsh) cycle
          nziota=nziota+1
       enddo
    enddo
    
    ! Cubic (aaa) 3rd-order interstate coupling constants (tau)
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m=1,nmodes
             if (tau_mask(m,s1,s2).eq.0) cycle
             if (abs(tau(m,s1,s2)).lt.thrsh) cycle
             nztau=nztau+1
          enddo
       enddo
    enddo

    ! Quartic (aaaa) 4rd-order intrastate coupling constants (epsilon)
    do s=1,nsta
       do m=1,nmodes
          if (epsilon_mask(m,s).eq.0) cycle
          if (abs(epsilon(m,s)).lt.thrsh) cycle
          nzepsilon=nzepsilon+1
       enddo
    enddo

    ! Quartic (aaaa) 4rd-order interstate coupling constants (xi)
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m=1,nmodes
             if (xi_mask(m,s1,s2).eq.0) cycle
             if (abs(xi(m,s1,s2)).lt.thrsh) cycle
             nzxi=nzxi+1
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
