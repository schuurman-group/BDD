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
    use symmetry
    use kdcglobal
    
    implicit none

    integer             :: unit,m,m1,m2,s,s1,s2,i,j,k,c,nl,ncurr,fel
    real(dp), parameter :: thrsh=5e-4
    character(len=2)    :: am,am1,am2,as,as1,as2,afel
    character(len=5)    :: aunit
    character(len=90)   :: string
    character(len=3)    :: amode
    character(len=8)    :: atmp

!----------------------------------------------------------------------
! nit label
!----------------------------------------------------------------------
    aunit=' , ev'

!----------------------------------------------------------------------
! Index of the electronic DOF
!----------------------------------------------------------------------
    fel=nmodes+1
    write(afel,'(i2)') fel

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
       write(am,'(i2)') m
       write(iop,'(a,F9.6,a)') &
            'omega_'//adjustl(am)//' =',freq(m),aunit
    enddo

    ! Energies
    write(iop,'(/,a)') '# Energies'
    do s=1,nsta
       write(as,'(i2)') s
       write(iop,'(a,F9.6,a)') &
            'E'//adjustl(as)//' = ',(q0pot(s)-q0pot(1))*eh2ev,aunit
    enddo

    ! 1st-order intrastate coupling constants (kappa)
    write(iop,'(/,a)') '# 1st-order intrastate coupling &
         constants (kappa)'
    do s=1,nsta
       write(as,'(i2)') s
       do m=1,nmodes
          if (kappa_mask(m,s).eq.0) cycle
          if (abs(kappa(m,s)).lt.thrsh) cycle
          write(am,'(i2)') m
          write(iop,'(a,F9.6,a)') 'kappa'//trim(adjustl(as))&
               //'_'//adjustl(am)//' = ',kappa(m,s),aunit
       enddo
    enddo

    ! 1st-order intrastate coupling constants (lambda)
    write(iop,'(/,a)') '# 1st-order interstate coupling &
         constants (lambda)'
    do s1=1,nsta-1
       write(as1,'(i2)') s1
       do s2=s1+1,nsta
          write(as2,'(i2)') s2
          do m=1,nmodes
             if (lambda_mask(m,s1,s2).eq.0) cycle
             if (abs(lambda(m,s1,s2)).lt.thrsh) cycle
             write(am,'(i2)') m
             write(iop,'(a,F9.6,a)') 'lambda'//trim(adjustl(as1))&
               //'_'//trim(adjustl(as2))//'_'//adjustl(am)//&
               ' = ',lambda(m,s1,s2),aunit
          enddo
       enddo
    enddo

    ! Quadratic 2nd-order intrastate coupling constants (gamma)
    write(iop,'(/,a)') '# Quadratic 2nd-order intrastate &
         coupling constants (gamma)'
    do s=1,nsta
       write(as,'(i2)') s
       do m=1,nmodes
          if (gamma_mask(m,m,s).eq.0) cycle
          if (abs(gamma(m,m,s)).lt.thrsh) cycle
          write(am,'(i2)') m
          write(iop,'(a,F9.6,a)') 'gamma'//trim(adjustl(as))&
               //'_'//trim(adjustl(am))//'_'//adjustl(am)&
               //' = ',gamma(m,m,s),aunit
       enddo
    enddo

    ! Bi-linear 2nd-order intrastate coupling constants (gamma)
    write(iop,'(/,a)') '# Bi-linear 2nd-order intrastate &
         coupling constants (gamma)'
    do s=1,nsta
       write(as,'(i2)') s
       do m1=1,nmodes-1
          do m2=m1+1,nmodes
             if (gamma_mask(m1,m2,s).eq.0) cycle
             if (abs(gamma(m1,m2,s)).lt.thrsh) cycle
             write(am1,'(i2)') m1
             write(am2,'(i2)') m2
             write(iop,'(a,F9.6,a)') 'gamma'//trim(adjustl(as))&
                  //'_'//trim(adjustl(am1))//'_'//adjustl(am2)&
                  //' = ',gamma(m1,m2,s),aunit
          enddo
       enddo
    enddo

    ! Quadratic 2nd-order interstate coupling constants (mu)
    write(iop,'(/,a)') '# Quadratic 2nd-order interstate &
         coupling constants (mu)'
    do s1=1,nsta-1
       write(as1,'(i2)') s1
       do s2=s1+1,nsta
          write(as2,'(i2)') s2          
          do m=1,nmodes
             if (mu_mask(m,m,s1,s2).eq.0) cycle
             if (abs(mu(m,m,s1,s2)).lt.thrsh) cycle
             write(am,'(i2)') m
             write(iop,'(a,F9.6,a)') 'mu'//trim(adjustl(as1))&
                  //'_'//trim(adjustl(as2))//'_'//&
                  trim(adjustl(am))//'_'//adjustl(am)//' = ',&
                  mu(m,m,s1,s2),aunit
          enddo
       enddo
    enddo

    ! Bi-linear 2nd-order interstate coupling constants (mu)
    write(iop,'(/,a)') '# Bi-linear 2nd-order interstate &
         coupling constants (mu)'
    do s1=1,nsta-1
       write(as1,'(i2)') s1
       do s2=s1+1,nsta
          write(as2,'(i2)') s2          
          do m1=1,nmodes-1
             do m2=m1+1,nmodes
                if (mu_mask(m1,m2,s1,s2).eq.0) cycle
                if (abs(mu(m1,m2,s1,s2)).lt.thrsh) cycle
                write(am1,'(i2)') m1
                write(am2,'(i2)') m2
                write(iop,'(a,F9.6,a)') 'mu'//trim(adjustl(as1))&
                     //'_'//trim(adjustl(as2))//'_'//&
                     trim(adjustl(am1))//'_'//adjustl(am2)//&
                     ' = ',mu(m1,m2,s1,s2),aunit
             enddo
          enddo
       enddo
    enddo

    ! Finishing line
    write(iop,'(/,a)') 'end-parameter-section'

!----------------------------------------------------------------------
! Write the Hamiltonian section
!----------------------------------------------------------------------
    ! Starting line
    write(iop,'(/,a)') 'hamiltonian-section'

    ! Modes section
    write(iop,'(/,38a)') ('-',i=1,38)
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

    ! Kinetic energy operator
    write(iop,'(/,a)') '# Kinetic energy'
    do m=1,nmodes
       write(am,'(i2)') m
       write(iop,'(a)') &
            'omega_'//adjustl(am)//'  |'//adjustl(am)//'  KE'
    enddo

    ! Zeroth-order potential: VEEs
    write(iop,'(/,a)') '# Zeroth-order potential: VEEs'
    do s=1,nsta
       write(as,'(i2)') s
       write(iop,'(a)') &
            'E'//adjustl(as)//'  |'//adjustl(afel)&
            //'  S'//trim(adjustl(as))//'&'//trim(adjustl(as))
    enddo

    ! Zeroth-order potential: Harmonic oscillators
    write(iop,'(/,a)') '# Zeroth-order potential: &
         Harmonic oscillators'
    do m=1,nmodes
       write(am,'(i2)') m
       write(iop,'(a)') &
            '0.5*omega_'//adjustl(am)//'  |'//adjustl(am)//'  q^2'
    enddo

    ! 1st-order intrastate coupling terms (kappa)
    write(iop,'(/,a)') '# 1st-order intrastate coupling &
         constants (kappa)'
    do m=1,nmodes
       write(am,'(i2)') m
       do s=1,nsta
          if (kappa_mask(m,s).eq.0) cycle
          if (abs(kappa(m,s)).lt.thrsh) cycle
          write(as,'(i2)') s
          write(iop,'(a)') 'kappa'//trim(adjustl(as))&
               //'_'//adjustl(am)&
               //'  |'//adjustl(am)//'  q'//'  |'//adjustl(afel)&
               //'  S'//trim(adjustl(as))//'&'//trim(adjustl(as))
       enddo
    enddo

    ! 1st-order interstate coupling terms (lambda)
    write(iop,'(/,a)') '# 1st-order interstate coupling &
         constants (lambda)'
    do m=1,nmodes
       write(am,'(i2)') m
       do s1=1,nsta-1
          write(as1,'(i2)') s1
          do s2=s1+1,nsta
             if (lambda_mask(m,s1,s2).eq.0) cycle
             if (abs(lambda(m,s1,s2)).lt.thrsh) cycle
             write(as2,'(i2)') s2
             write(iop,'(a)') 'lambda'//trim(adjustl(as1)) &
                  //'_'//trim(adjustl(as2))//'_'//adjustl(am) &
                  //'  |'//adjustl(am)//'  q'//'  |'//adjustl(afel)&
                  //'  S'//trim(adjustl(as1))//'&'//trim(adjustl(as2))
          enddo
       enddo
    enddo

    ! Quadratic 2nd-order intrastate coupling constants (gamma)
    write(iop,'(/,a)') '# Quadratic 2nd-order intrastate &
         coupling constants (gamma)'
    do m=1,nmodes
       write(am,'(i2)') m
       do s=1,nsta
          if (gamma_mask(m,m,s).eq.0) cycle
          if (abs(gamma(m,m,s)).lt.thrsh) cycle
          write(as,'(i2)') s
          write(iop,'(a)') '0.5*gamma'//trim(adjustl(as))&
               //'_'//trim(adjustl(am))//'_'//adjustl(am)&
               //'  |'//adjustl(am)//'  q^2'//'  |'//adjustl(afel)&
               //'  S'//trim(adjustl(as))//'&'//trim(adjustl(as))
       enddo
    enddo

    ! Bi-linear 2nd-order intrastate coupling constants (gamma)
    write(iop,'(/,a)') '# Bi-linear 2nd-order intrastate &
         coupling constants (gamma)'
    do m1=1,nmodes-1
       do m2=m1+1,nmodes
          do s=1,nsta
             if (gamma_mask(m1,m2,s).eq.0) cycle
             if (abs(gamma(m1,m2,s)).lt.thrsh) cycle
             write(am1,'(i2)') m1
             write(am2,'(i2)') m2
             write(as,'(i2)') s
             write(iop,'(a)') 'gamma'//trim(adjustl(as))&
                  //'_'//trim(adjustl(am1))//'_'//adjustl(am2)&
                  //'  |'//adjustl(am1)//'  q'&
                  //'  |'//adjustl(am2)//'  q'&
                  //'  |'//adjustl(afel)&
                  //'  S'//trim(adjustl(as))//'&'//trim(adjustl(as))
          enddo
       enddo
    enddo

    ! Quadratic 2nd-order interstate coupling constants (mu)
    write(iop,'(/,a)') '# Bi-linear 2nd-order intrastate &
         coupling constants (gamma)'
    do m=1,nmodes
       write(am,'(i2)') m
       do s1=1,nsta-1
          do s2=s1+1,nsta
             if (mu_mask(m,m,s1,s2).eq.0) cycle
             if (abs(mu(m,m,s1,s2)).lt.thrsh) cycle
             write(as1,'(i2)') s1
             write(as2,'(i2)') s2
             write(iop,'(a)') '0.5*mu'//trim(adjustl(as1))&
                  //'_'//trim(adjustl(as2))&
                  //'_'//trim(adjustl(am))//'_'//adjustl(am)&
                  //'  |'//adjustl(am)//'  q^2'//'  |'//adjustl(afel)&
                  //'  S'//trim(adjustl(as1))//'&'//trim(adjustl(as2))
          enddo
       enddo
    enddo

    ! Bi-linear 2nd-order interstate coupling constants (mu)
    write(iop,'(/,a)') '# Bi-linear 2nd-order intrastate &
         coupling constants (gamma)'
    do m1=1,nmodes-1
       write(am1,'(i2)') m1
       do m2=m1+1,nmodes
          write(am2,'(i2)') m2
          do s1=1,nsta-1
             write(as1,'(i2)') s1
             do s2=s1+1,nsta
                write(as1,'(i2)') s1
                if (mu_mask(m1,m2,s1,s2).eq.0) cycle
                if (abs(mu(m1,m2,s1,s2)).lt.thrsh) cycle
                write(iop,'(a)') 'mu'//trim(adjustl(as1))&
                     //'_'//trim(adjustl(as2))&
                     //'_'//trim(adjustl(am1))//'_'//adjustl(am2)&
                     //'  |'//adjustl(am1)//'  q'&
                     //'  |'//adjustl(am2)//'  q'&
                     //'  |'//adjustl(afel)&
                     //'  S'//trim(adjustl(as1))//'&'//trim(adjustl(as2))
             enddo
          enddo
       enddo
    enddo
       
    ! Finishing line
    write(iop,'(/,a)') 'end-hamiltonian-section'

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

end module opermod
