module fields
  use text_options, only: text_option

  implicit none

  private 

  public :: init_fields, finish_fields
  public :: read_parameters, wnml_fields, check_fields
  public :: advance, force_maxwell_reinit
  public :: reset_init, set_init_fields
  public :: fields_init_response, set_dump_and_read_response
  public :: dump_response_to_file
  public :: init_fields_parameters
  public :: finish_fields_parameters
  public :: set_overrides
  public :: init_fields_level_1, init_fields_level_2
  public :: finish_fields_level_1, finish_fields_level_2
  public :: fieldopt_switch, fieldopt_implicit, fieldopt_test, fieldopt_local, fieldopt_gf_local

  !> Made public for unit tests
  public :: fields_pre_init
  public :: remove_zonal_flows_switch
  !> Made public for replay
  public :: allocate_arrays

  interface fieldlineavgphi
     module procedure fieldlineavgphi_loc
     module procedure fieldlineavgphi_tot
  end interface

  ! knobs
  integer :: fieldopt_switch
  logical :: remove_zonal_flows_switch
  logical :: force_maxwell_reinit
  integer, parameter :: fieldopt_implicit = 1, fieldopt_test = 2, fieldopt_local = 3, fieldopt_gf_local = 4
  logical :: dump_response, read_response
  logical :: initialized = .false.
  logical :: parameters_read = .false.
  logical :: exist

  !EGH made fieldopts module level for overrides
  type (text_option), dimension (6), parameter :: fieldopts = &
       (/ text_option('default', fieldopt_implicit), &
          text_option('implicit', fieldopt_implicit), &
          text_option('test', fieldopt_test),&
          text_option('local', fieldopt_local),&
          text_option('gf_local', fieldopt_gf_local),&
          text_option('implicit_local', fieldopt_local)/)
contains

  subroutine check_fields(report_unit)
    use fields_local, only: do_smart_update, minNRow
    implicit none
    integer, intent(in) :: report_unit
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       write (report_unit, fmt="('The field equations will be advanced in time implicitly.')")
       if(dump_response) write (report_unit, fmt="('The response matrix will be dumped to file.')")
       if(read_response) write (report_unit, fmt="('The response matrix will be read from file.')")
    case (fieldopt_test)
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The field equations will only be tested.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    case (fieldopt_local)
       write (report_unit, fmt="('The field equations will be advanced in time implicitly with decomposition respecting g_lo layout.')")
       if(dump_response) write (report_unit, fmt="('The response matrix will be dumped to file.')")
       if(read_response) write (report_unit, fmt="('The response matrix will be read from file.')")
       write(report_unit, fmt="('Using a min block size of ',I0)") minNrow
       if(do_smart_update) write(report_unit, fmt="('Using optimised field update.')")
    case (fieldopt_gf_local)
       write (report_unit, fmt="('The field equations will be advanced in time implicitly with decomposition respecting gf_lo layout.')")
    end select
  end subroutine check_fields

  subroutine wnml_fields(unit)
    use fields_local, only: do_smart_update, minNRow
    implicit none
    integer, intent(in) :: unit
    if (.not. exist) return 
    write (unit, *)
    write (unit, fmt="(' &',a)") "fields_knobs"
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       write (unit, fmt="(' field_option = ',a)") '"implicit"'
    case (fieldopt_test)
       write (unit, fmt="(' field_option = ',a)") '"test"'
    case (fieldopt_local)
       write (unit, fmt="(' field_option = ',a)") '"local"'
       write (unit, fmt="(' minNrow = ',I0)") minNrow
       write (unit, fmt="(' do_smart_update = ',L1)") do_smart_update
    case (fieldopt_gf_local)
       write (unit, fmt="(' field_option = ',a)") '"local gf"'
    end select
    if(dump_response) write (unit, fmt="(' dump_response = ',L1)") dump_response
    if(read_response) write (unit, fmt="(' read_response = ',L1)") read_response
    write (unit, fmt="(' /')")
  end subroutine wnml_fields

  subroutine set_overrides(opt_ov)
    use overrides, only: optimisations_overrides_type
    use file_utils, only: error_unit
    use text_options, only: get_option_value
    use fields_local, only: minnrow
    use fields_implicit, only: field_subgath
    use fields_local, only: do_smart_update
    use fields_local, only: field_local_allreduce
    use fields_local, only: field_local_allreduce_sub
    type(optimisations_overrides_type), intent(in) :: opt_ov
    integer :: ierr
    if (opt_ov%override_field_option) then
       ierr = error_unit()
       call get_option_value &
            (opt_ov%field_option, fieldopts, fieldopt_switch, &
            ierr, "field_option in set_overrides",.true.)
    end if
    if (opt_ov%override_minnrow) minnrow = opt_ov%minnrow
    if (opt_ov%override_field_subgath) field_subgath = opt_ov%field_subgath
    if (opt_ov%override_do_smart_update) do_smart_update = & 
      opt_ov%do_smart_update
    if (opt_ov%override_field_local_allreduce) field_local_allreduce = & 
      opt_ov%field_local_allreduce
    if (opt_ov%override_field_local_allreduce_sub) field_local_allreduce_sub =& 
      opt_ov%field_local_allreduce_sub
  end subroutine set_overrides

  !> Calls all initialisations required for init_fields_implicit/local, 
  !! reads parameters and allocates field arrays
  subroutine fields_pre_init
    use theta_grid, only: init_theta_grid
    use run_parameters, only: init_run_parameters
    use dist_fn, only: init_dist_fn, gf_lo_integrate
    use init_g, only: ginit, init_init_g
    use antenna, only: init_antenna
    use unit_tests, only: debug_message
    use kt_grids, only: naky, ntheta0
    use mp, only: nproc, proc0
    implicit none
    integer, parameter :: verb=3
    
    call init_fields_parameters

    call debug_message(verb, "init_fields: init_theta_grid")
    call init_theta_grid
    
!CMR,30/3/2009:
! call init_init_g before init_run_parameters to read delt from restart file

    call debug_message(verb, "init_fields: init_init_g")
    call init_init_g
    call debug_message(verb, "init_fields: init_run_parameters")
    call init_run_parameters
    call debug_message(verb, "init_fields: init_dist_fn")
    call init_dist_fn

    if(nproc .lt. ntheta0*naky .and. fieldopt_switch .eq. fieldopt_gf_local) then
       fieldopt_switch = fieldopt_local
       gf_lo_integrate = .false.
       if(proc0) then
          write(*,*) 'gf local fields cannot be used as you are running less MPI processes than there are'
          write(*,*) 'naky*ntheta0 points.  Defaulting to local fields.  You need to use at least',naky*ntheta0
          write(*,*) 'MPI processes for this simulation'
       end if
    end if

    !call debug_message(verb, "init_fields: init_parameter_scan")
    !call init_parameter_scan
    call debug_message(verb, "init_fields: init_antenna")
    call init_antenna !Must come before allocate_arrays so we know if we need apar_ext
    !call debug_message(verb, "init_fields: read_parameters")
    !call read_parameters
    call debug_message(verb, "init_fields: allocate_arrays")
    call allocate_arrays
  end subroutine fields_pre_init
  
  subroutine init_fields_parameters
    use unit_tests, only: debug_message
    implicit none
    integer, parameter :: verb=3
    if (parameters_read) return
    call debug_message(verb, "init_fields: read_parameters")
    call read_parameters
    parameters_read = .true.
  end subroutine init_fields_parameters

  subroutine finish_fields_parameters
    parameters_read = .false.
  end subroutine finish_fields_parameters

  subroutine init_fields_level_1
    use unit_tests, only: debug_message
    implicit none
    integer, parameter :: verb=3
    call debug_message(verb, "init_fields: allocate_arrays")
    call allocate_arrays
  end subroutine init_fields_level_1

  subroutine finish_fields_level_1
    call finish_fields
  end subroutine finish_fields_level_1

  subroutine init_fields_level_2
    call init_fields
  end subroutine init_fields_level_2

  subroutine finish_fields_level_2
    call reset_init
  end subroutine finish_fields_level_2

  subroutine fields_init_response
    use fields_implicit, only: init_fields_implicit
    use fields_test, only: init_fields_test
    use fields_local, only: init_fields_local
    use fields_gf_local, only: init_fields_gf_local
    use unit_tests, only: debug_message
    implicit none
    integer, parameter :: verb=3
    logical, parameter :: debug = .false.
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call debug_message(verb, &
         "fields::fields_init_response init_fields_implicit")
       call init_fields_implicit
    case (fieldopt_test)
       call debug_message(verb, "fields::fields_init_response init_fields_test")
       call init_fields_test
    case (fieldopt_local)
       call debug_message(verb, &
         "fields::fields_init_response init_fields_local")
       call init_fields_local
    case (fieldopt_gf_local)
       call debug_message(verb, &
         "fields::fields_init_response init_fields_gf_local")
       call init_fields_gf_local
    case default
       !Silently ignore unsupported field options
    end select
  end subroutine fields_init_response

  subroutine init_fields
!CMR,18/2/2011:
! add optional logical arg "noalloc" to avoid array allocations in ingen
    use theta_grid, only: init_theta_grid
    use run_parameters, only: init_run_parameters
    use dist_fn, only: init_dist_fn
    use init_g, only: ginit, init_init_g
    use nonlinear_terms, only: nl_finish_init => finish_init
    use antenna, only: init_antenna
    use kt_grids, only: gridopt_switch, gridopt_box, kwork_filter
    implicit none
    logical :: restarted
    logical, parameter :: debug=.false.

    if (initialized) return
    initialized = .true.
    
    call fields_pre_init

    call fields_init_response

! Turn on nonlinear terms.
    if (debug) write(6,*) "init_fields: nl_finish_init"
    call nl_finish_init

    ! EGH Commented out the following lines as they are now
    ! handled by gs2_init

    !if (debug) write(6,*) "init_fields: ginit"
    !call ginit (restarted)
    !if (restarted .and. .not. force_maxwell_reinit) return
    !if (debug) write(6,*) "init_fields: init_antenna"
    !call init_antenna

    !!Set the initial fields
    !call set_init_fields

    !If running in flux tube disable evolution of ky=kx=0 mode
    if(gridopt_switch.eq.gridopt_box) kwork_filter(1,1)=.true.
  end subroutine init_fields

  !>Force the current 
  subroutine dump_response_to_file(suffix)
    use fields_implicit, only: dump_response_to_file_imp
    use fields_local, only: dump_response_to_file_local
    use fields_gf_local, only: dump_response_to_file_gf_local
    implicit none
    character(len=*), intent(in), optional :: suffix 
    !Note can pass optional straight through as long as also optional
    !in called routine (and not different routines combined in interface)
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call dump_response_to_file_imp(suffix)
    case (fieldopt_local)
       call dump_response_to_file_local(suffix)
    case (fieldopt_gf_local)
       call dump_response_to_file_gf_local(suffix)
    case default
       !Silently ignore unsupported field options
    end select
  end subroutine dump_response_to_file

  subroutine set_init_fields
    use fields_implicit, only: init_allfields_implicit
    use fields_test, only: init_phi_test
    use mp, only: proc0
    use fields_local, only: init_allfields_local
    use fields_gf_local, only: init_allfields_gf_local
    use dist_fn, only: gf_lo_integrate
    use gs2_layouts, only: gf_local_fields
    use kt_grids, only: naky, ntheta0
    implicit none
    logical, parameter :: debug=.false.
    ! New flow shear algorithm only supports implicit fieldopt_implicit. -- NDC 08/18
    if(proc0.and.debug) write(6,*) "Syncing fields with g."
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       if (debug) write(6,*) "init_fields: init_allfields_implicit"
       call init_allfields_implicit
    case (fieldopt_test)
       if (debug) write(6,*) "init_fields: init_phi_test"
       call init_phi_test
    case (fieldopt_local)
       if (debug) write(6,*) "init_fields: init_allfields_local"
       call init_allfields_local
    case (fieldopt_gf_local)
       if (debug) write(6,*) "init_fields: init_allfields_gf_local"
       if(.not. gf_lo_integrate .or. .not. gf_local_fields) then
          if(proc0) then
             write(*,*) 'gf_lo_integrate',gf_lo_integrate,'gf_local_fields',gf_local_fields
             write(*,*) 'gf local fields cannot be used by gf_lo_integrate'
             write(*,*) 'defaulting to local fields'
             write(*,*) 'if you want to use gf local fields then set gf_lo_integrate to true in dist_fn_knobs'
             write(*,*) 'and gf_local_fields to true in layout_knobs and make sure you are running on more MPI'
             write(*,*) 'mpi processes than naky*ntheta0.  For this simulation that is: ',naky*ntheta0
          end if
          fieldopt_switch = fieldopt_local
          call init_allfields_local
       else
          call init_allfields_gf_local
       end if
    end select
  end subroutine set_init_fields

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast, mp_abort
    use fields_implicit, only: field_subgath
    use fields_local, only: minNRow
    use fields_local, only: do_smart_update, field_local_allreduce, field_local_allreduce_sub
    use fields_local, only: field_local_tuneminnrow,  field_local_nonblocking_collectives
    use fields_arrays, only: response_file
    use file_utils, only: run_name
    implicit none
    character(20) :: field_option
    character(len=256) :: response_dir
    namelist /fields_knobs/ field_option, remove_zonal_flows_switch, field_subgath, force_maxwell_reinit,&
         dump_response, read_response, minNRow, do_smart_update, field_local_allreduce, field_local_allreduce_sub,&
         response_dir, field_local_tuneminnrow,  field_local_nonblocking_collectives
    integer :: ierr, in_file

    if (proc0) then
       field_option = 'default'
       remove_zonal_flows_switch = .false.
       field_subgath=.false.
       force_maxwell_reinit=.true. ! NDCQUEST: why do we not rely on saved phi ?
       dump_response=.false.
       read_response=.false.
       do_smart_update=.false.
       field_local_allreduce=.false.
       field_local_allreduce_sub=.false.
       field_local_tuneminnrow=.false.
       field_local_nonblocking_collectives=.false.

       response_dir=''
       in_file = input_unit_exist ("fields_knobs", exist)
       if (exist) read (unit=in_file, nml=fields_knobs)

       ierr = error_unit()
       call get_option_value &
            (field_option, fieldopts, fieldopt_switch, &
            ierr, "field_option in fields_knobs",.true.)

       if(trim(response_dir).eq.'')then
          write(response_file,'(A)') trim(run_name)
       else
          write(response_file,'(A,"/",A)') trim(response_dir),trim(run_name)
       endif

    end if



    call broadcast (fieldopt_switch)
    call broadcast (remove_zonal_flows_switch)
    call broadcast (field_subgath)
    call broadcast (force_maxwell_reinit)
    call broadcast (dump_response)
    call broadcast (read_response)

    !Setup response file location
    call broadcast(response_dir)
    call broadcast(response_file)

    !Set the solve type specific flags
    call set_dump_and_read_response(dump_response, read_response)
    !select case (fieldopt_switch)
    !case (fieldopt_implicit)
    !case (fieldopt_test)
    !case (fieldopt_local)
    ! EGH commented it out as we want to be able to 
    ! change field option dynamically and so we want
    ! these parameters to be broadcast in all circumstances.
       call broadcast (minNrow)
       call broadcast (do_smart_update)
       call broadcast (field_local_allreduce)
       call broadcast (field_local_tuneminnrow)
       call broadcast (field_local_nonblocking_collectives)
       call broadcast (field_local_allreduce_sub)
    !end select

  end subroutine read_parameters

  subroutine set_dump_and_read_response(dump_flag, read_flag)
    use fields_implicit, only: dump_response_imp => dump_response, read_response_imp=>read_response
    use fields_local, only: dump_response_loc => dump_response, read_response_loc=>read_response
    use fields_gf_local, only: dump_response_gf => dump_response, read_response_gf=>read_response
    implicit none
    logical, intent(in) :: dump_flag, read_flag
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       dump_response_imp=dump_flag
       read_response_imp=read_flag
    case (fieldopt_local)
       dump_response_loc=dump_flag
       read_response_loc=read_flag
    case (fieldopt_gf_local)
       dump_response_gf=dump_flag
       read_response_gf=read_flag
    case default
       !Silently ignore unsupported field types
    end select
  end subroutine set_dump_and_read_response

  subroutine allocate_arrays
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, &
        mixed_flowshear, &
        explicit_flowshear ! NDCTESTmichaelnew
    use antenna, only: no_driver
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew, apar_ext, &
        aparold, &
        phistar_old, phistar_new ! NDCTESTmichaelnew
    use fields_arrays, only: gf_phi, gf_apar, gf_bpar, gf_phinew, gf_aparnew, gf_bparnew
    use unit_tests, only: debug_message
    implicit none
    integer, parameter :: verb=3
    logical :: michael_exp = .true. ! NDCTESTswitchexp

    if (.not. allocated(phi)) then
       call debug_message(verb, 'fields::allocate_arrays allocating')
       allocate (     phi (-ntgrid:ntgrid,ntheta0,naky))
       allocate (    apar (-ntgrid:ntgrid,ntheta0,naky))
       allocate (   bpar (-ntgrid:ntgrid,ntheta0,naky))
       allocate (  phinew (-ntgrid:ntgrid,ntheta0,naky))
       allocate ( aparnew (-ntgrid:ntgrid,ntheta0,naky))
       allocate (bparnew (-ntgrid:ntgrid,ntheta0,naky))
       if(mixed_flowshear .or. explicit_flowshear) then
           allocate(aparold(-ntgrid:ntgrid,ntheta0,naky))
           if(explicit_flowshear .and. michael_exp) then
               allocate ( phistar_old (-ntgrid:ntgrid,ntheta0,naky)) ! NDCTESTmichaelnew
               allocate ( phistar_new (-ntgrid:ntgrid,ntheta0,naky)) ! NDCTESTmichaelnew
           end if
       end if
       if(fieldopt_switch .eq. fieldopt_gf_local) then
          !AJ It should be possible to reduce the size of these by only allocating them
          !AJ the extend of it and ik in gf_lo.  However, it would need to be done carefully 
          !AJ to ensure the proc0 has the full space if we are still reducing data to proc0 
          !AJ for diagnostics.
          !AJ With some careful thought it should also be possible to remove this all together and 
          !AJ simply use the arrays above, but that would need checking
          allocate (    gf_phi (-ntgrid:ntgrid,ntheta0,naky))
          allocate (   gf_apar (-ntgrid:ntgrid,ntheta0,naky))
          allocate (   gf_bpar (-ntgrid:ntgrid,ntheta0,naky))          
          allocate ( gf_phinew (-ntgrid:ntgrid,ntheta0,naky))
          allocate (gf_aparnew (-ntgrid:ntgrid,ntheta0,naky))
          allocate (gf_bparnew (-ntgrid:ntgrid,ntheta0,naky))          
       end if
    endif
    !AJ This shouldn't be necessary as it is done in the fields code....
    phi = 0.; phinew = 0.
    apar = 0.; aparnew = 0.
    bpar = 0.; bparnew = 0.
    if(explicit_flowshear .and. michael_exp) then
        phistar_old = 0. ! NDCTESTmichaelnew
        phistar_new = 0. ! NDCTESTmichaelnew
    end if
    if(fieldopt_switch .eq. fieldopt_gf_local) then
       gf_phi = 0.; gf_phinew = 0.
       gf_apar = 0.; gf_aparnew = 0.
       gf_bpar = 0.; gf_bparnew = 0.
    end if
    if(.not.allocated(apar_ext).and.(.not.no_driver))then
       allocate (apar_ext (-ntgrid:ntgrid,ntheta0,naky))
       apar_ext = 0.
    endif
  end subroutine allocate_arrays

  subroutine advance (istep)
    use fields_implicit, only: advance_implicit
    use fields_test, only: advance_test
    use fields_local, only: advance_local
    use fields_gf_local, only: advance_gf_local

    implicit none
    integer, intent (in) :: istep

    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call advance_implicit (istep, remove_zonal_flows_switch)
    case (fieldopt_test)
       call advance_test (istep)
    case (fieldopt_local)
       call advance_local (istep, remove_zonal_flows_switch)
    case (fieldopt_gf_local)
       call advance_gf_local (istep, remove_zonal_flows_switch)
    end select
  end subroutine advance

  !This routine has a potentially misleading name and isn't used anywhere
  subroutine kperp (ntgrid_output, akperp)
    use theta_grid, only: delthet
    use kt_grids, only: naky, aky, ntheta0, kperp2
    use run_parameters, only: fphi, fapar, fbpar
    use fields_arrays, only: phinew, aparnew, bparnew
    implicit none
    integer, intent (in) :: ntgrid_output
    real, dimension (:,:), intent (out) :: akperp
    real :: anorm
    integer :: ik, it

    do ik = 1, naky
       do it = 1, ntheta0
          anorm = sum(abs(phinew(-ntgrid_output:ntgrid_output,it,ik)*fphi &
                         + aparnew(-ntgrid_output:ntgrid_output,it,ik)*fapar &
                         + bparnew(-ntgrid_output:ntgrid_output,it,ik)*fbpar)**2 &
                      *delthet(-ntgrid_output:ntgrid_output))
          if (anorm < 2.0*epsilon(0.0) .or. aky(ik) == 0.0) then
             akperp(it,ik) = 0.0
          else
             akperp(it,ik) &
                  = sqrt(sum(kperp2(-ntgrid_output:ntgrid_output,it,ik) &
                     *abs(phinew(-ntgrid_output:ntgrid_output,it,ik)*fphi &
                          + aparnew(-ntgrid_output:ntgrid_output,it,ik)*fapar &
                          + bparnew(-ntgrid_output:ntgrid_output,it,ik)*fbpar)**2 &
                     *delthet(-ntgrid_output:ntgrid_output))/anorm)
          end if
       end do
    end do
  end subroutine kperp

  !This routine isn't used anywhere
  subroutine fieldlineavgphi_loc (ntgrid_output, it, ik, phiavg)
    use theta_grid, only: ntgrid, drhodpsi, gradpar, bmag, delthet
    use fields_arrays, only: phi
    implicit none
    integer, intent (in) :: ntgrid_output, ik, it
    complex, intent (out) :: phiavg
    real, dimension (-ntgrid:ntgrid) :: jac

    jac = 1.0/abs(drhodpsi*gradpar*bmag)
    phiavg = sum(phi(-ntgrid_output:ntgrid_output,it,ik) &
                 *jac(-ntgrid_output:ntgrid_output) &
                 *delthet(-ntgrid_output:ntgrid_output)) &
            /sum(delthet(-ntgrid_output:ntgrid_output) &
                 *jac(-ntgrid_output:ntgrid_output))
  end subroutine fieldlineavgphi_loc

  !This doesn't look like a useful routine and isn't used anywhere
  subroutine fieldlineavgphi_tot (phiavg)
    use theta_grid, only: ntgrid, drhodpsi, gradpar, bmag
!    use theta_grid, only: delthet
    implicit none
    complex, dimension (:,:), intent (out) :: phiavg
    real, dimension (-ntgrid:ntgrid) :: jac
!    integer :: ntg

!    ntg = ntgrid
    jac = 1.0/abs(drhodpsi*gradpar*bmag)
!    phiavg = sum(phi(-ntg:ntg,:,:)*jac(-ntg:ntg)*delthet(-ntg:ntg)) &
!         /sum(delthet(-ntg:ntg)*jac(-ntg:ntg))
    phiavg = 0.
    write(*,*) 'error in fields'
  end subroutine fieldlineavgphi_tot


  !!> This generates a flux surface average of phi. 

  !subroutine flux_surface_average_phi (phi_in, phi_average)
    !use theta_grid, only: ntgrid, drhodpsi, gradpar, bmag, delthet
    !use kt_grids, only: ntheta0, naky

    !implicit none
    !complex, intent (in) :: phi_in
    !complex, intent (out) :: phi_average
    !complex, dimension(-ntgrid:ntgrid,1:ntheta0,1:naky) :: phi_fieldline_avg
    !integer it, ig

    !call fieldline_average_phi(phi_in, phi_fieldline_avg)
    !do it = 1,ntheta0
      !do ig = -ntgrid,ntgrid
        !phi_average(ig, it, :) = sum(phi_fieldline_avg(ig, it, :))/real(naky)
      !end do
    !end do

  !end subroutine fieldline_average_phi

  subroutine reset_init
    use fields_implicit, only: fi_reset => reset_init
    use fields_test, only: ft_reset => reset_init
    use fields_local, only: fl_reset => reset_fields_local
    use fields_gf_local, only: flgf_reset => reset_fields_gf_local
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use fields_arrays, only: gf_phi, gf_apar, gf_bpar, gf_phinew, gf_aparnew, gf_bparnew
    implicit none
    initialized  = .false.
    phi = 0.
    phinew = 0.
    apar = 0.
    aparnew = 0.
    bpar = .0
    bparnew = 0.
    if(fieldopt_switch .eq. fieldopt_gf_local) then
       gf_phi = 0.
       gf_apar = 0.
       gf_bpar = 0.
       gf_phinew = 0.
       gf_aparnew = 0.
       gf_bparnew = 0.
    end if
    !What about apar_ext?
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       call fi_reset
    case (fieldopt_test)
       call ft_reset
    case (fieldopt_local)
       call fl_reset
    case (fieldopt_gf_local)
       call flgf_reset
    end select
  end subroutine reset_init

  subroutine finish_fields

    use fields_implicit, only: implicit_reset => reset_init
    use fields_test, only: test_reset => reset_init
    use fields_local, only: finish_fields_local
    use fields_gf_local, only: finish_fields_gf_local
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew, &
        aparold, &
        phistar_old, phistar_new ! NDCTESTmichaelnew
    use fields_arrays, only: apar_ext, gf_phi, gf_apar, gf_bpar
    use fields_arrays, only: apar_ext, gf_phinew, gf_aparnew, gf_bparnew
    use unit_tests, only: debug_message

    implicit none
    integer, parameter :: verbosity = 3


    initialized  = .false.
!AJ Does these need zero'd if they are to be deallocated below?
    phi = 0.
    phinew = 0.
    apar = 0.
    aparnew = 0.
    bpar = .0
    bparnew = 0.
    call debug_message(verbosity, &
      'fields::finish_fields zeroed fields')
    
    select case (fieldopt_switch)
    case (fieldopt_implicit)
       ! This line is no longer necessary
       ! because fields_implicit::reset_init is
       ! called by finish_fields_level_2
       !call implicit_reset
    case (fieldopt_test)
       ! This line is no longer necessary
       ! because fields_test::reset_init is
       ! called by finish_fields_level_2
       !call test_reset
    case (fieldopt_local)
       call finish_fields_local
    case (fieldopt_gf_local)
       call finish_fields_gf_local
    end select

    call debug_message(verbosity, &
      'fields::finish_fields called subroutines')

    if (allocated(phi)) deallocate (phi, apar, bpar, phinew, aparnew, bparnew)
    if (allocated(aparold)) deallocate (aparold)
    if (allocated(phistar_old)) deallocate (phistar_old, phistar_new) ! NDCTESTmichaelnew
    if (allocated(gf_phi)) deallocate(gf_phi, gf_apar, gf_bpar, gf_phinew, gf_aparnew, gf_bparnew)
    if (allocated(apar_ext)) deallocate (apar_ext)
    call debug_message(verbosity, &
      'fields::finish_fields deallocated fields')

  end subroutine finish_fields

end module fields
