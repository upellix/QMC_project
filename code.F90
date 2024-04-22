!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                              !
!              PROGRAM TO PERFORM QUANTUM MONTE CARLO SIMULATION               !
!                          Autor : Pelliccia Umberto                           !
!                                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program monte_carlo

  implicit none

  ! Define variables and parameters

  double precision, parameter    :: A_to_Bohr = 1.8897259886d0

  double precision               :: ave, err, a, E_ref, dt, tau
  double precision, allocatable  :: energy(:), accep(:), r_nu(:,:), Z(:)
  integer*8                      :: nmax
  integer                        :: n_el, n_nu, nruns, irun, i
  character*6                    :: geom_file
  character*4                    :: method
  character                      :: atom

!------------------------------------------------------------------------------!
!                             INIZIALIZATION                                   !
!------------------------------------------------------------------------------!
  
  ! Read information from the input file

  open(10, file='mc_input.inp', status='old')

  read(10,*) 
  read(10,*)
  read(10,*) geom_file           ! Geometry file
  read(10,*) n_el                ! Number of electrons
  read(10,*) a                   ! Wavefunction parameter
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*) method              ! Method for monte-carlo calculation
  read(10,*) nmax                ! Number of steps
  read(10,*) nruns               ! Number of walkers
  read(10,*) dt                  ! Time steps
  
  if (method == 'dif') then

    read(10,*)
    read(10,*) E_ref             ! Reference energy 
    read(10,*) tau               ! Projection time

  end if

  close(10)

  ! Read information from geometry file

  open(20, file=geom_file, status='old')

  read(20,*) n_nu                ! Number of atoms
  read(20,*)

  allocate(r_nu(n_nu,3), Z(n_nu), energy(nruns), accep(nruns))

  do i = 1, n_nu
  
    read(20,*) atom, r_nu(i,:)   ! Atom type and position
    
    if (atom == 'H') then 

      Z(i) = 1.d0

    else if (atom == 'He') then

      Z(i) = 2.d0

    else

      stop 'ERROR : Invalid atom, please use only "H" and "He" '

    end if

  end do

  close (20)

  ! Conversion of nuclei's position from Angstrom to a.u.

  r_nu(:,:) = r_nu(:,:) * A_to_Bohr

!------------------------------------------------------------------------------!
!                 MAIN LOOP FOR MONTE CARLO CALCULATION                        !
!------------------------------------------------------------------------------!

  ! Run variational monte carlo calculation

  if (method == 'vmc') then

    do irun = 1, nruns

      call variational_montecarlo(a, nmax, dt, energy(irun), accep(irun), r_nu, n_el, n_nu, Z)
  
    end do

    write(*,*) 'Results for Variational Monte Carlo method:'

  ! Run pure diffusion monte carlo simulation

  else if (method == 'pdmc') then

    do irun = 1, nruns

      call pd_montecarlo(a, nmax, dt, energy(irun), accep(irun), tau, E_ref, r_nu, n_el, n_nu, Z)

    end do

    write(*,*) 'Results for Pure Diffusion Monte Carlo method:'
  
  else
    
    stop 'ERROR : Invalid method, please chose between "vmc" and "pdmc" '
 
  end if

!------------------------------------------------------------------------------!
!                                 OUTPUT                                       !
!------------------------------------------------------------------------------!

  ! Write the results for energy and acceptance

  call ave_error(energy, nruns ,ave, err)
  write(*,80) 'ENERGY = ', ave, '+/-', err
    
  call ave_error(accep, nruns, ave, err)
  write(*,80) 'ACCEPTANCE = ', ave, '+/-', err

  80 format(a13, f20.12, a7, f20.6)

  ! Generate the output file for the plot of energy

  call plot_file(nruns, energy(:))

end program monte_carlo
