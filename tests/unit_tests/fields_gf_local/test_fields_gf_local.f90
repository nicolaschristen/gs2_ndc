
!> A program that runs unit tests on the fields_gf_local module.
!! Currently success is defined as giving the same results as
!! the old implicit module...could be improved by an absolute test
!!
!! This is free software released under the MIT license
!!   Written by: Adrian Jackson (adrianj@epcc.ed.ac.uk)
!!   Based on the fields_local test
program test_fields_gf_local
  use unit_tests
  use init_g, only: ginit
  use fields_gf_local, only: fields_gf_local_unit_test_init_fields_matrixlocal, fields_gf_local_functional
  use fields_implicit, only: fields_implicit_unit_test_init_fields_implicit
  use fields, only: fieldopt_switch,fieldopt_implicit,fieldopt_gf_local
  use fields, only: fields_pre_init
  use egrid
  use mp, only: init_mp, finish_mp, broadcast, mp_comm
  use file_utils, only: init_file_utils
  use species, only: init_species, nspec
  use constants, only: pi
  use dist_fn, only: init_dist_fn,  gf_lo_integrate
  use dist_fn_arrays, only: g,gnew
  use kt_grids, only: naky, ntheta0
  use theta_grid, only: ntgrid
  use gs2_layouts, only: g_lo
  use gs2_main, only: initialize_gs2, gs2_program_state_type, initialize_equations, initialize_diagnostics
  use gs2_main, only: finalize_gs2, finalize_equations, finalize_diagnostics
  use gs2_init, only: init, init_level_list
  implicit none
  real :: eps
  type(gs2_program_state_type) :: state
  complex, dimension (:,:,:), allocatable :: gbak
  complex, dimension (:,:,:), allocatable :: phi_imp, apar_imp, bpar_imp
  complex, dimension (:,:,:), allocatable :: phi_loc, apar_loc, bpar_loc
  logical :: restarted

  ! General config
  eps = 1.0e-10

  if (precision(eps).lt. 11) eps = eps * 1000.0

  ! Set up depenencies
  call init_mp

  state%init%opt_ov%override_field_option = .true.
  state%init%opt_ov%field_option = 'gf_local'
  state%init%opt_ov%override_gf_lo_integrate = .true.
  state%init%opt_ov%gf_lo_integrate = .true.
  state%init%opt_ov%override_gf_local_fields = .true.
  state%init%opt_ov%gf_local_fields = .true.
  state%init%opt_ov%override_simple_gf_decomposition = .true.
  state%init%opt_ov%simple_gf_decomposition = .true.
  state%mp_comm_external = .true.
  state%mp_comm = mp_comm

  call initialize_gs2(state)
  call initialize_equations(state)
  call initialize_diagnostics(state)

  state%print_times = .false.

  call announce_module_test('fields_gf_local')

  allocate(gbak(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
  allocate(phi_imp(-ntgrid:ntgrid,ntheta0,naky))
  allocate(apar_imp(-ntgrid:ntgrid,ntheta0,naky))
  allocate(bpar_imp(-ntgrid:ntgrid,ntheta0,naky))
  allocate(phi_loc(-ntgrid:ntgrid,ntheta0,naky))
  allocate(apar_loc(-ntgrid:ntgrid,ntheta0,naky))
  allocate(bpar_loc(-ntgrid:ntgrid,ntheta0,naky))

  phi_imp = 0; apar_imp = 0; bpar_imp =0;  
  phi_loc = 0; apar_loc = 0; bpar_loc =0;


  call ginit(restarted)
  gbak = g

  if (fields_gf_local_functional()) then

    call announce_test('init_fields_matrixlocal')
    call process_test(fields_gf_local_unit_test_init_fields_matrixlocal(), 'init_fields_matrixlocal')
    call test_local_advance(phi_loc,apar_loc,bpar_loc)

  else 

    write (*,*) "WARNING: fields_gf_local is non-functional in your build. &
      & Skipping the fields_gf_local unit test. &
      & If you are using the PGI compilers this is to be expected. "
  end if 
  
  
  call finalize_diagnostics(state)
  call finalize_equations(state)
  call finalize_gs2(state)

  state%init%opt_ov%override_field_option = .true.
  state%init%opt_ov%field_option = 'local'
  state%init%opt_ov%override_gf_lo_integrate = .true.
  state%init%opt_ov%gf_lo_integrate = .false.
  state%init%opt_ov%override_gf_local_fields = .true.
  state%init%opt_ov%gf_local_fields = .false.
  state%mp_comm_external = .true.
  state%mp_comm = mp_comm
  call initialize_gs2(state)
  call initialize_equations(state)
  call initialize_diagnostics(state)

  g = gbak
  gnew = g
  
  call announce_test('init_fields_implicit')
  call process_test(fields_implicit_unit_test_init_fields_implicit(), 'init_fields_implicit')

  call test_implicit_advance(phi_imp,apar_imp,bpar_imp)

  call announce_test('advance')
  call process_test(check_results(eps,phi_imp,phi_loc,apar_imp,apar_loc,bpar_imp,bpar_loc),'advance')

  call finalize_diagnostics(state)
  call finalize_equations(state)
  call finalize_gs2(state)

  deallocate(gbak,phi_imp,apar_imp,bpar_imp,phi_loc,apar_loc,bpar_loc)

  call close_module_test('fields_gf_local')
  call finish_mp
contains

  subroutine test_implicit_advance(phi_i, apar_i, bpar_i)
    use fields_implicit, only: init_allfields_implicit, advance_implicit
    use fields_arrays
    use run_parameters, only: nstep
    use fields, only: remove_zonal_flows_switch
    implicit none
    complex, dimension (:,:,:) :: phi_i, apar_i, bpar_i
    integer :: istep

    !Now we setup the initial fields to be consistent with g
    call init_allfields_implicit()

    !Now we can do a timestep (or lots)
    do istep=1,nstep
        call advance_implicit(istep, remove_zonal_flows_switch)
    enddo

    !!Now we store the results
    phi_i=phinew
    apar_i=aparnew
    bpar_i=bparnew

  end subroutine test_implicit_advance

  subroutine test_local_advance(phi_l, apar_l, bpar_l)
    use fields_arrays
    use run_parameters, only: nstep
    use fields, only: remove_zonal_flows_switch
    use fields_gf_local, only:  advance_gf_local, init_allfields_gf_local

    implicit none
    complex, dimension (:,:,:) :: phi_l, apar_l, bpar_l
    integer :: istep

    !Now we setup the initial fields to be consistent with g
    call init_allfields_gf_local()
    !Now we can do a timestep (or lots)
    do istep=1,nstep
        call advance_gf_local(istep, remove_zonal_flows_switch)
    enddo

    !Now we store the results
    phi_l=phinew
    apar_l=aparnew
    bpar_l=bparnew


  end subroutine test_local_advance


   function check_results(eps,phi_imp,phi_loc,apar_imp,apar_loc,bpar_imp,bpar_loc)

    use kt_grids, only: kwork_filter
    use mp, only: iproc,proc0
    use gs2_layouts, only: g_lo
    implicit none

    real, intent(in) :: eps
    complex, dimension (:,:,:) :: phi_imp, apar_imp, bpar_imp
    complex, dimension (:,:,:) :: phi_loc, apar_loc, bpar_loc
    integer :: ik,it,ikit
    character(len=29) :: message
    logical :: check_results
    
    check_results = .true.

    if(proc0) then
       do ik = 1,g_lo%naky
          do it = 1,g_lo%ntheta0
             write(message, fmt="(A19, I2, A6, I2)") 'talue of phi,  it =', it, ' ik = ', ik
             call announce_check(message)
             call process_check(check_results, agrees_with(phi_imp(:,it,ik), phi_loc(:,it,ik), eps), message)
             write(message, fmt="(A19, I2, A6, I2)") 'talue of apar, it =', it, ' ik = ', ik
             call announce_check(message)
             call process_check(check_results, agrees_with(apar_imp(:,it,ik), apar_loc(:,it,ik), eps), message)
             write(message, fmt="(A19, I2, A6, I2)") 'talue of bpar, it =', it, ' ik = ', ik
             call announce_check(message)
             call process_check(check_results, agrees_with(bpar_imp(:,it,ik), bpar_loc(:,it,ik), eps), message)       
          end do
       end do
    end if

    do ikit=1,g_lo%ikitrange
       ik = g_lo%local_ikit_points(ikit)%ik 
       it = g_lo%local_ikit_points(ikit)%it
       if(kwork_filter(it,ik)) cycle
       write(message, fmt="(A19, I2, A6, I2)") 'value of phi,  it =', it, ' ik = ', ik
       call announce_check(message)
       call process_check(check_results, agrees_with(phi_imp(:, it, ik), phi_loc(:, it, ik), eps), message)
       write(message, fmt="(A19, I2, A6, I2)") 'value of apar, it =', it, ' ik = ', ik
       call announce_check(message)
       call process_check(check_results, agrees_with(apar_imp(:, it, ik), apar_loc(:, it, ik), eps), message)
       write(message, fmt="(A19, I2, A6, I2)") 'value of bpar, it =', it, ' ik = ', ik
       call announce_check(message)
       call process_check(check_results, agrees_with(bpar_imp(:, it, ik), bpar_loc(:, it, ik), eps), message)       
    end do

  end function check_results


end program test_fields_gf_local
