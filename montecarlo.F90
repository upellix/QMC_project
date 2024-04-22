!------------------------------------------------------------------------------!
!                         DRIFT VECTOR CALCULATION                             !
!------------------------------------------------------------------------------!

subroutine drift(a, n_nu, r_el, r_nu, d)

  implicit none

  double precision, intent(in)   :: a, r_el(3), r_nu(n_nu,3)
  integer,          intent(in)   :: n_nu
  double precision, intent(out)  :: d(3)

  double precision               :: ar_inv, dr(3)
  integer                        :: i

  d(:) = 0.d0

  do i = 1, n_nu

    dr(:) = r_el(:) - r_nu(i,:)
    ar_inv = -a / dsqrt(dr(1) * dr(1) + dr(2) * dr(2) + dr(3) * dr(3))

    d(:) = d(:) - dr(:) * ar_inv

  end do

end subroutine drift



!------------------------------------------------------------------------------!
!                   VARIATIONAL MONTE CARLO CALCULATION                        !
!------------------------------------------------------------------------------!

subroutine variational_montecarlo(a, nmax, dt, energy ,accep, r_nu, n_el, n_nu, Z)

  implicit none

  double precision, intent(in)   :: a, dt, r_nu(n_nu,3), Z(n_nu)
  integer,          intent(in)   :: n_nu, n_el
  integer*8,        intent(in)   :: nmax
  double precision, intent(out)  :: energy, accep

  integer                        :: i
  integer*8                      :: istep, n_accep
  double precision               :: sq_dt, chi(n_el,3), d2_old(n_el), d2_new(n_el) 
  double precision               :: prod(n_el), u, psi_old, psi_new, argexpo, q
  double precision               :: r_old(n_el,3), r_new(n_el,3), d_old(n_el,3), d_new(n_el,3)
  double precision, external     :: e_loc, psi

  sq_dt = dsqrt(dt)

  energy  = 0.d0
  n_accep = 0_8

  ! Initialization for each electron of the system

  do i = 1, n_el

    call random_gauss(r_old(i,:), 3)

    call drift(a, n_nu,  r_old(i,:), r_nu, d_old(i,:))
 
    ! Calculate the square of drift vector

    d2_old(i) = sum( d_old(i,:) * d_old(i,:) )

  end do

  psi_old = psi(a, n_el, n_nu, r_old, r_nu)

  ! Propagation

  do istep = 1, nmax
  
    energy = energy + e_loc(a, n_el, n_nu, r_old, r_nu, Z)

    ! Calculate the new move of every electron

    do i = 1, n_el

      call random_gauss(chi(i,:), 3)
      r_new(i,:) = r_old(i,:) + dt * d_old(i,:) + chi(i,:) * sq_dt

      call drift(a, n_nu, r_new(i,:), r_nu, d_new(i,:))

      ! Recalculate the square of drift vector with the new electrons position
      
      d2_new(i) = sum( d_new(i,:) * d_new(i,:) )

    end do

    ! Calculate the new wavefunction

    psi_new = psi(a, n_el, n_nu, r_new, r_nu)

    ! Perform the metropolis calculation

    prod(:) = (d_new(:,1) + d_old(:,1))*(r_new(:,1) - r_old(:,1)) + &
              (d_new(:,2) + d_old(:,2))*(r_new(:,2) - r_old(:,2)) + &
              (d_new(:,3) + d_old(:,3))*(r_new(:,3) - r_old(:,3))

    argexpo = 0.5d0 * sum( (d2_new(:) - d2_old(:)) * dt + prod(:) )

    q = psi_new / psi_old
    q = dexp(-argexpo) * q * q

    ! Accept or reject the new move

    call random_number(u)

    if (u <= q) then

      n_accep = n_accep + 1_8

      r_old(:,:) = r_new(:,:)
      d_old(:,:) = d_new(:,:)
      d2_old(:)  = d2_new(:)
      psi_old    = psi_new

    end if

  end do

  ! Calculate the energy and the acceptance

  energy = energy / dble(nmax)
  accep = dble(n_accep) / dble(nmax)

end subroutine variational_montecarlo



!------------------------------------------------------------------------------!
!                   PURE DIFFUSION MONTE CARLO CALCULATION                     !
!------------------------------------------------------------------------------!

subroutine pd_montecarlo(a, nmax, dt, energy ,accep, tau ,E_ref, r_nu, n_el, n_nu, Z)

  implicit none

  double precision, intent(in)   :: a, dt, tau, E_ref, r_nu(n_nu,3), Z(n_nu)
  integer,          intent(in)   :: n_el, n_nu
  integer*8,        intent(in)   :: nmax
  double precision, intent(out)  :: energy, accep

  integer                        :: i
  integer*8                      :: istep, n_accep
  double precision               :: sq_dt, chi(n_el,3), d2_old(n_el), d2_new(n_el)
  double precision               :: prod(n_el), u, psi_old, psi_new, argexpo, q
  double precision               :: r_old(n_el,3), r_new(n_el,3), d_old(n_el,3), d_new(n_el,3)
  double precision               :: e, w, norm, tau_current
  double precision, external     :: e_loc, psi

  sq_dt = dsqrt(dt)

  ! Inizialization

  energy      = 0.d0
  n_accep     = 0_8
  norm        = 0.d0

  w           = 1.d0
  tau_current = 0.d0

  ! Inizialization for each electron in the system

  do i = 1, n_el

    call random_gauss(r_old(i,:), 3)

    call drift(a, n_nu, r_old(i,:), r_nu, d_old(i,:))

    ! Calculate the square of drift vector

    d2_old(i) = sum( d_old(i, :) * d_old(i,:) )
  
  end do

  psi_old = psi(a, n_el, n_nu, r_old, r_nu)

  ! Propagation

  do istep = 1, nmax
  
    e = e_loc(a, n_el, n_nu, r_old, r_nu, Z)

    ! Calculate the weight for pd monte carlo

    w = w * dexp(-dt * (e - E_ref))

    ! Normalize the energy and update its value

    norm = norm + w
    energy = energy + w * e

    ! Update the projection time

    tau_current = tau_current + dt

    ! Reset when tau is reached

    if (tau_current > tau) then
      
      w           = 1.d0
      tau_current = 0.d0
  
    end if

    ! Calculate the new move of every electron

    do i = 1, n_el

      call random_gauss(chi(i,:), 3)
      r_new(i,:) = r_old(i,:) + dt * d_old(i,:) + chi(i,:) * sq_dt

      call drift(a, n_nu, r_new(i,:), r_nu, d_new(i,:))

      ! Recalculate the square of drift vector with the new electron position
      
      d2_new(i) = sum( d_new(i, :) * d_new(i,:) )
    
    end do

    ! Calculate the new wavefunction

    psi_new = psi(a, n_el, n_nu, r_new, r_nu)

    ! Perform the metropolis calculation

    prod(:) = (d_new(:,1) + d_old(:,1))*(r_new(:,1) - r_old(:,1)) + &
              (d_new(:,2) + d_old(:,2))*(r_new(:,2) - r_old(:,2)) + &
              (d_new(:,3) + d_old(:,3))*(r_new(:,3) - r_old(:,3))

    argexpo = 0.5d0 * sum( (d2_new(:) - d2_old(:)) * dt + prod(:) )

    q = psi_new / psi_old
    q = dexp(-argexpo) * q * q

    ! Accept or reject the new move

    call random_number(u)

    if (u <= q) then

      n_accep = n_accep + 1_8

      r_old(:,:) = r_new(:,:)
      d_old(:,:) = d_new(:,:)
      d2_old(:)  = d2_new(:)
      psi_old    = psi_new

    end if

  end do

  ! Calculate the energy and the acceptance

  energy = energy / norm
  accep = dble(n_accep) / dble(nmax)

end subroutine pd_montecarlo
