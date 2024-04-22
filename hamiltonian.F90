!------------------------------------------------------------------------------!
!                      POTENTIAL ENERGY CALCULATION                            !
!------------------------------------------------------------------------------!

double precision function potential(n_el, n_nu, r_el, r_nu, Z)

  implicit none

  double precision, intent(in)  :: r_el(n_el,3), r_nu(n_nu,3), Z(n_nu)
  integer,          intent(in)  :: n_el, n_nu

  double precision              :: dist, dr(3)
  integer                       :: i, j

  potential = 0.d0

  ! Calculate electron-nucleus interaction

  do i = 1, n_nu
    do j = 1, n_el

      dr(:) = r_el(j,:) - r_nu(i,:)
      dist = dsqrt(dr(1) * dr(1) + dr(2) * dr(2) + dr(3) * dr(3))

      if (dist > 1.d-7) then
        
        potential = potential - Z(i) / dist

      else
        
        write(*,*) "ERROR potential electron-nucleus diverges"
        potential = potential - huge(1.d0)

      end if

    end do
  end do

  ! Calculate electron-electron interaction

  do i = 1, n_el - 1
    do j = i + 1, n_el

      dr(:) = r_el(j,:) - r_el(i,:)
      dist = dsqrt(dr(1) * dr(1) + dr(2) * dr(2) + dr(3) * dr(3))

      if (dist > 1.d-7) then
        
        potential = potential + 1.d0 / dist

      else
        
        write(*,*) "ERROR potential electron-electron diverges"
        potential = potential + huge(1.d0)

      end if

    end do
  end do

  ! Calculate nucleus-nucleus interaction 

  do i = 1, n_nu - 1
    do j = i + 1, n_nu

      dr(:) = r_nu(j,:) - r_nu(i,:)
      dist = dsqrt(dr(1) * dr(1) + dr(2) * dr(2) + dr(3) * dr(3))

      if (dist > 1.d-7) then
        
        potential = potential + Z(i) * Z(j) / dist

      else
        
        write(*,*) "ERROR potential nucleus-nucleus diverges"
        potential = potential + huge(1.d0)

      end if

    end do
  end do

end function potential



!------------------------------------------------------------------------------!
!                        KINETIC ENERGY CALCULATION                            !
!------------------------------------------------------------------------------!

double precision function kinetic(a, n_el, n_nu, r_el, r_nu)

  implicit none

  double precision, intent(in)  :: a, r_el(n_el,3), r_nu(n_nu,3)
  integer,          intent(in)  :: n_el, n_nu

  double precision              :: dist, dr(3), k, local_k, total_psi
  integer                       :: i, j
  double precision, external    :: psi_nucleus, psi

  k = 0.d0

  ! Calculate the wavefunction dependent kinetic term

  total_psi = psi(a, n_el, n_nu, r_el, r_nu)

  do i = 1, n_nu
    
    local_k = psi_nucleus(a, n_el, r_el, r_nu(i,:)) / total_psi

    do j = 1, n_el
      
      dr(:) = r_el(j,:) - r_nu(i,:)
      dist = dsqrt(dr(1) * dr(1) + dr(2) * dr(2) + dr(3) * dr(3))

      if (dist > 1.d-7) then
        
        k = k + a / dist * local_k

      else
        
        write(*,*) "ERROR kinetic energy diverges"
        k = k + huge(1.d0)

      end if

    end do
  end do

  ! sum the dependent term with the constant term to fin the kinetic energy

  kinetic = - (a * a * 0.5d0 * dble(n_el)) + k

end function kinetic



!------------------------------------------------------------------------------!
!                          LOCAL ENERGY CALCULATION                            !
!------------------------------------------------------------------------------!

double precision function e_loc(a, n_el, n_nu, r_el, r_nu, Z)
  
  implicit none

  double precision, intent(in)  :: a, r_el(n_el,3), r_nu(n_nu,3), Z(n_nu)
  integer,          intent(in)  :: n_el, n_nu

  double precision, external    :: kinetic
  double precision, external    :: potential

  ! Sum potential and kinetic contribution to get the local energy

  e_loc = kinetic(a, n_el, n_nu, r_el, r_nu) + potential(n_el, n_nu, r_el, r_nu, Z) 

end function e_loc


!------------------------------------------------------------------------------!
!                                 WAVEFUNCTION                                 !
!------------------------------------------------------------------------------!

! Calculate the Hartree-product for nuclei

double precision function psi_nucleus(a, n_el, r_el, r_nu)

  implicit none

  double precision, intent(in)  :: a, r_el(n_el,3), r_nu(3)
  integer,          intent(in)  :: n_el

  double precision              :: dist, dr(3)
  integer                       :: i

  dist = 0.d0

  do i = 1, n_el
      
    dr(:) = r_el(i,:) - r_nu(:)
    dist = dist + dsqrt(dr(1) * dr(1) + dr(2) * dr(2) + dr(3) * dr(3))

  end do

  psi_nucleus = dexp(-a * dist)

end function psi_nucleus



! Calculate the wavefunction

double precision function psi(a, n_el, n_nu, r_el, r_nu)

  implicit none

  double precision, intent(in)  :: a, r_el(n_el,3), r_nu(n_nu,3)
  integer,          intent(in)  :: n_el, n_nu

  integer                       :: i
  double precision, external    :: psi_nucleus
  
  psi = 0.d0

  do i = 1, n_nu

    psi = psi + psi_nucleus(a, n_el, r_el, r_nu(i,:))

  end do

end function psi




