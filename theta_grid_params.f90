module theta_grid_params
  implicit none

  private

  public :: init_theta_grid_params
  public :: finish_theta_grid_params
  public :: wnml_theta_grid_params, write_trinity_parameters
  public :: rhoc, rmaj, r_geo, eps, epsl, qinp, shat, alpmhd
  public :: pk, geoType, aSurf, shift, shiftVert, mMode, nMode, deltam, deltan
  public :: deltampri, deltanpri, thetam, thetan
  public :: btor_slab, betaprim, ntheta, nperiod
  public :: set_overrides

  real :: rhoc, rmaj, r_geo, eps, epsl
  real :: qinp, shat, alpmhd, pk
  integer :: geoType, mMode, nMode
  real :: aSurf, shift, shiftVert, raxis, zaxis, deltam, deltan
  real :: akappa, tri, delta2, delta3, deltampri, deltanpri
  real :: akappri, tripri, thetam, thetan, thetak, thetad, theta2, theta3
  real :: btor_slab, betaprim
  real :: asym, asympri ! JRB - these two geometry parameters are NOT functional

  integer :: ntheta, nperiod

  logical :: initialized = .false.
  real :: kp = -1.
  logical :: exist

contains
  subroutine init_theta_grid_params
    use unit_tests, only: debug_message
    implicit none
    integer, parameter :: verb=3
!    logical, save :: initialized = .false.

    
    call debug_message(verb, "theta_grid_params::init_theta_grid_params start")
    if (initialized) return
    initialized = .true.

    call debug_message(verb, &
      "theta_grid_params::init_theta_grid_params call read_parameters")
    call read_parameters

    call debug_message(verb, "theta_grid_params::init_theta_grid_params end")
  end subroutine init_theta_grid_params

  subroutine finish_theta_grid_params
    implicit none
    initialized = .false.
  end subroutine finish_theta_grid_params

  subroutine read_parameters
    use file_utils, only: input_unit, input_unit_exist
    use unit_tests, only: debug_message
    use constants, only: pi=>pi
    implicit none
    integer, parameter :: verb=3
    integer :: in_file
    character(4) :: ntheta_char

!CMR,2/2/2011: add btor_slab
! btor_slab = btor/bpol defines direction of a flow relative to B in slab 
! geometry, where flow is by definition in the toroidal direction.

    namelist /theta_grid_parameters/ rhoc, rmaj, r_geo, eps, epsl, &
         qinp, shat, alpmhd, pk, geoType, aSurf, shift, shiftVert, raxis, zaxis, &
         mMode, nMode, deltam, deltan, akappa, tri, delta2, delta3, deltampri, deltanpri, &
         akappri, tripri, thetam, thetan, thetak, thetad, theta2, theta3, &
         ntheta, nperiod, kp, btor_slab, asym, asympri
    
       call debug_message(verb, "theta_grid_params::read_parameters start")


    rhoc = 0.5
    rmaj = 3.0
    r_geo = 3.0
    eps = 0.3
    epsl = 0.3
    qinp = 1.5
    shat = 0.75
    pk = 0.3
    geoType = 0 ! NB The default for this MUST be zero otherwise the Trinity interface will break
    aSurf = 1.0
    shift = 0.0
    shiftVert = 0.0
    raxis = 3.0
    zaxis = 0.0
    mMode = 2
    nMode = 3
    deltam = 1.0
    deltan = 1.0
    akappa = 1.0
    tri = 0.0
    delta2 = 1.0
    delta3 = 1.0
    deltampri = 0.0
    deltanpri = 0.0
    akappri = 0.0
    tripri = 0.0
    thetam = 0.0
    thetan = 0.0
    thetak = 0.0
    thetad = 0.0
    theta2 = 0.0
    theta3 = 0.0
    asym = 0.0
    asympri = 0.0
    btor_slab = 0.0
    ntheta = 24
    nperiod = 2
    in_file = input_unit_exist("theta_grid_parameters", exist)
    if (exist) read (unit=in_file, nml=theta_grid_parameters)

    if (kp > 0.) pk = 2.*kp

    ! Print warning if the user has specified non-default values for the
    ! geometrical input parameters of other geometry types ! JRB
    select case (geoType)
       case (0)
          if (mMode/=2) write (*,*) "WARNING: ignoring value of mMode, value not needed"
          if (nMode/=3) write (*,*) "WARNING: ignoring value of nMode, value not needed"
          if (deltam/=1.0) write (*,*) "WARNING: ignoring value of deltam, did you mean akappa?"
          if (deltan/=1.0) write (*,*) "WARNING: ignoring value of deltan, did you mean tri?"
          if (delta2/=1.0) write (*,*) "WARNING: ignoring value of delta2, did you mean akappa?"
          if (delta3/=1.0) write (*,*) "WARNING: ignoring value of delta3, did you mean tri?"
          if (deltampri/=0.0) write (*,*) "WARNING: ignoring value of deltampri, did you mean akappri?"
          if (deltanpri/=0.0) write (*,*) "WARNING: ignoring value of deltanpri, did you mean tripri?"
          if (thetam/=0.0) write (*,*) "WARNING: ignoring value of thetam, did you mean thetak?"
          if (thetan/=0.0) write (*,*) "WARNING: ignoring value of thetan, did you mean thetad?"
          if (theta2/=0.0) write (*,*) "WARNING: ignoring value of theta2, did you mean thetak?"
          if (theta3/=0.0) write (*,*) "WARNING: ignoring value of theta3, did you mean thetad?"
          if (aSurf/=1.0) write (*,*) "WARNING: ignoring value of aSurf, value not needed"
          if (raxis/=3.0) write (*,*) "WARNING: ignoring value of raxis, did you mean shift?"
          if (zaxis/=0.0) write (*,*) "WARNING: ignoring value of zaxis, did you mean shiftVert?"
          deltam = akappa
          deltan = tri
          deltampri = akappri
          deltanpri = tripri
          thetam = thetak
          thetan = thetad
       case (1)
          if (mMode/=2) write (*,*) "WARNING: ignoring value of mMode, value not needed"
          if (nMode/=3) write (*,*) "WARNING: ignoring value of nMode, value not needed"
          if (deltam/=1.0) write (*,*) "WARNING: ignoring value of deltam, did you mean delta2?"
          if (deltan/=1.0) write (*,*) "WARNING: ignoring value of deltan, did you mean delta3?"
          if (akappa/=1.0) write (*,*) "WARNING: ignoring value of akappa, did you mean delta2?"
          if (tri/=0.0) write (*,*) "WARNING: ignoring value of tri, did you mean delta3?"
          if (deltampri/=0.0) write (*,*) "WARNING: ignoring value of deltampri, value not needed"
          if (deltanpri/=0.0) write (*,*) "WARNING: ignoring value of deltanpri, value not needed"
          if (akappri/=0.0) write (*,*) "WARNING: ignoring value of akappri, value not needed"
          if (tripri/=0.0) write (*,*) "WARNING: ignoring value of tripri, value not needed"
          if (thetam/=0.0) write (*,*) "WARNING: ignoring value of thetam, did you mean theta2?"
          if (thetan/=0.0) write (*,*) "WARNING: ignoring value of thetan, did you mean theta3?"
          if (thetak/=0.0) write (*,*) "WARNING: ignoring value of thetak, did you mean theta2?"
          if (thetad/=0.0) write (*,*) "WARNING: ignoring value of thetad, did you mean theta3?"
          if (shift/=0.0) write (*,*) "WARNING: ignoring value of shift, did you mean Raxis?"
          if (shiftVert/=0.0) write (*,*) "WARNING: ignoring value of shiftVert, did you mean Zaxis?"

          if (raxis>=Rmaj+aSurf) write (*,*) "WARNING: Raxis>=Rmag+aSurf, may specify non-nested flux surfaces"
          if (raxis<=Rmaj-aSurf) write (*,*) "WARNING: Raxis<=Rmag-aSurf, may specify non-nested flux surfaces"
          if (zaxis>=aSurf) write (*,*) "WARNING: Zaxis>=aSurf, may specify non-nested flux surfaces"
          if (zaxis<=-aSurf) write (*,*) "WARNING: Zaxis<=-aSurf, may specify non-nested flux surfaces"

          shift = raxis
          shiftVert = zaxis
          deltam = delta2
          deltan = delta3
          thetam = theta2
          thetan = theta3
       case (2)
          if (delta2/=1.0) write (*,*) "WARNING: ignoring value of delta2, did you mean deltam?"
          if (delta3/=1.0) write (*,*) "WARNING: ignoring value of delta3, did you mean deltan?"
          if (akappa/=1.0) write (*,*) "WARNING: ignoring value of akappa, did you mean deltam?"
          if (tri/=0.0) write (*,*) "WARNING: ignoring value of tri, did you mean deltan?"
          if (akappri/=0.0) write (*,*) "WARNING: ignoring value of akappri, did you mean deltampri?"
          if (tripri/=0.0) write (*,*) "WARNING: ignoring value of tripri, did you mean deltanpri?"
          if (theta2/=0.0) write (*,*) "WARNING: ignoring value of theta2, did you mean thetam?"
          if (theta3/=0.0) write (*,*) "WARNING: ignoring value of theta3, did you mean thetan?"
          if (thetak/=0.0) write (*,*) "WARNING: ignoring value of thetak, did you mean thetam?"
          if (thetad/=0.0) write (*,*) "WARNING: ignoring value of thetad, did you mean thetan?"
          if (aSurf/=1.0) write (*,*) "WARNING: ignoring value of aSurf, value not needed"
          if (raxis/=3.0) write (*,*) "WARNING: ignoring value of raxis, did you mean shift?"
          if (zaxis/=0.0) write (*,*) "WARNING: ignoring value of zaxis, did you mean shiftVert?"
          
       case (3)
          if (delta2/=1.0) write (*,*) "WARNING: ignoring value of delta2, did you mean deltam?"
          if (delta3/=1.0) write (*,*) "WARNING: ignoring value of delta3, did you mean deltan?"
          if (akappa/=1.0) write (*,*) "WARNING: ignoring value of akappa, did you mean deltam?"
          if (tri/=0.0) write (*,*) "WARNING: ignoring value of tri, did you mean deltan?"
          if (akappri/=0.0) write (*,*) "WARNING: ignoring value of akappri, did you mean deltampri?"
          if (tripri/=0.0) write (*,*) "WARNING: ignoring value of tripri, did you mean deltanpri?"
          if (theta2/=0.0) write (*,*) "WARNING: ignoring value of theta2, did you mean thetam?"
          if (theta3/=0.0) write (*,*) "WARNING: ignoring value of theta3, did you mean thetan?"
          if (thetak/=0.0) write (*,*) "WARNING: ignoring value of thetak, did you mean thetam?"
          if (thetad/=0.0) write (*,*) "WARNING: ignoring value of thetad, did you mean thetan?"
          if (aSurf/=1.0) write (*,*) "WARNING: ignoring value of aSurf, value not needed"
          if (raxis/=3.0) write (*,*) "WARNING: ignoring value of raxis, did you mean shift?"
          if (zaxis/=0.0) write (*,*) "WARNING: ignoring value of zaxis, did you mean shiftVert?"
          
       case default
          write (*,*) "ERROR: invalid analytic geometry specification"
    end select

    if (asym/=0.0) write (*,*) "WARNING: ignoring value of asym (not functional)"
    if (asympri/=0.0) write (*,*) "WARNING: ignoring value of asympri (not functional)"

  end subroutine read_parameters

  subroutine wnml_theta_grid_params(unit)
    implicit none
    integer, intent(in) :: unit
    if (.not.exist) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "theta_grid_parameters"
    write (unit, fmt="(' ntheta =   ',i4)") ntheta
    write (unit, fmt="(' nperiod =  ',i4)") nperiod
    write (unit, fmt="(' rhoc =     ',f7.4)") rhoc
    write (unit, fmt="(' Rmaj =     ',f7.4)") rmaj
    write (unit, fmt="(' R_geo =    ',f7.4)") r_geo
    write (unit, fmt="(' eps =      ',f7.4)") eps
    write (unit, fmt="(' epsl =     ',f7.4)") epsl
    write (unit, fmt="(' qinp =     ',f7.4)") qinp
    write (unit, fmt="(' shat =     ',f7.4)") shat
    write (unit, fmt="(' alpmhd =   ',f7.4)") alpmhd
    write (unit, fmt="(' pk =       ',f7.4)") pk
    write (unit, fmt="(' kp =       ',f7.4)") kp
    write (unit, fmt="(' geoType =  ',i4)") geoType
    write (unit, fmt="(' aSurf =    ',f7.4)") aSurf
    write (unit, fmt="(' shift =    ',f7.4)") shift
    write (unit, fmt="(' shiftVert =',f7.4)") shiftVert
    write (unit, fmt="(' mMode =    ',i4)") mMode
    write (unit, fmt="(' nMode =    ',i4)") nMode
    write (unit, fmt="(' deltam =   ',f7.4)") deltam
    write (unit, fmt="(' deltan =   ',f7.4)") deltan
    write (unit, fmt="(' deltampri =',f7.4)") deltampri
    write (unit, fmt="(' deltanpri =',f7.4)") deltanpri
    write (unit, fmt="(' thetam =   ',f7.4)") thetam
    write (unit, fmt="(' thetan =   ',f7.4)") thetan
    write (unit, fmt="(' btor_slab =',f7.4)") btor_slab
    write (unit, fmt="(' /')")
  end subroutine wnml_theta_grid_params

  subroutine write_trinity_parameters(trinpars_unit)
    integer, intent(in) :: trinpars_unit
    write (trinpars_unit, "(A22)") '&theta_grid_parameters'
    write (trinpars_unit, *) ' rhoc = ', rhoc
    write (trinpars_unit, *) ' qinp = ', qinp
    write (trinpars_unit, *) ' shat = ', shat
    write (trinpars_unit, *) ' rmaj = ', rmaj
    write (trinpars_unit, *) ' r_geo = ', r_geo
    write (trinpars_unit, *) ' geoType = ', geoType
    write (trinpars_unit, *) ' aSurf = ', aSurf
    write (trinpars_unit, *) ' shift = ', shift
    write (trinpars_unit, *) ' shiftVert = ', shiftVert
    write (trinpars_unit, *) ' mMode = ', mMode
    write (trinpars_unit, *) ' nMode = ', nMode
    write (trinpars_unit, *) ' deltam = ', deltam
    write (trinpars_unit, *) ' deltan = ', deltan
    write (trinpars_unit, *) ' deltampri = ', deltampri
    write (trinpars_unit, *) ' deltanpri = ', deltanpri
    write (trinpars_unit, *) ' thetam = ', thetam
    write (trinpars_unit, *) ' thetan = ', thetan
    write (trinpars_unit, *) ' betaprim = ', betaprim
    write (trinpars_unit, "(A1)") '/'
  end subroutine write_trinity_parameters


  subroutine set_overrides(mgeo_ov)
    use overrides, only: miller_geometry_overrides_type
    type(miller_geometry_overrides_type), intent(in) :: mgeo_ov
          !write (*,*) 'Calling tgpso'
    if (mgeo_ov%override_rhoc) rhoc = mgeo_ov%rhoc
    if (mgeo_ov%override_qinp) qinp = mgeo_ov%qinp
    if (mgeo_ov%override_shat) shat = mgeo_ov%shat
    if (mgeo_ov%override_rgeo_lcfs) r_geo = mgeo_ov%rgeo_lcfs
    if (mgeo_ov%override_rgeo_local) rmaj = mgeo_ov%rgeo_local
    if (mgeo_ov%override_geoType) geoType = mgeo_ov%geoType
    if (mgeo_ov%override_aSurf) aSurf = mgeo_ov%aSurf
    if (mgeo_ov%override_shift) shift = mgeo_ov%shift
    if (mgeo_ov%override_shiftVert) shiftVert = mgeo_ov%shiftVert
    if (mgeo_ov%override_mMode) mMode = mgeo_ov%mMode
    if (mgeo_ov%override_nMode) nMode = mgeo_ov%nMode
    if (mgeo_ov%override_deltam) then
      deltam = mgeo_ov%deltam
    else if (mgeo_ov%override_akappa) then
      deltam = mgeo_ov%akappa
    end if
    if (mgeo_ov%override_deltan) then
      deltan = mgeo_ov%deltan
    else if (mgeo_ov%override_tri) then
      deltan = mgeo_ov%tri
    end if
    if (mgeo_ov%override_deltampri) then 
      deltampri = mgeo_ov%deltampri
    else if (mgeo_ov%override_akappri) then 
      deltampri = mgeo_ov%akappri
    end if
    if (mgeo_ov%override_deltanpri) then 
      deltanpri = mgeo_ov%deltanpri
    else if (mgeo_ov%override_tripri) then 
      deltanpri = mgeo_ov%tripri
    end if
    if (mgeo_ov%override_thetam) thetam = mgeo_ov%thetam
    if (mgeo_ov%override_thetan) thetan = mgeo_ov%thetan
    if (mgeo_ov%override_betaprim) betaprim = mgeo_ov%betaprim
  end subroutine set_overrides


end module theta_grid_params
