!> A container for the arrays that are used to store the distribution function among other things.
!!  These need to be accessible at a lower dependency level than the dist_fn module itself.
!! These arrays are allocated in the function dist_fn::allocate_arrays. 

module dist_fn_arrays

  implicit none

  private

  public :: g, gnew, g_restart_tmp, kx_shift, theta0_shift, vpa, vpac
  public :: vperp2, vpar, ittp, aj0, aj1
!AJ
  public :: vpa_gf, vperp2_gf, aj0_gf, aj1_gf
  public :: c_rate
  public :: g_adjust, check_g_bouncepoints
#ifdef LOWFLOW
  public :: hneoc, vparterm, wdfac, wstarfac, wdttpfac
#endif

  public :: kx_shift_old

  public :: bessel_tdep_ptr_type
  public :: aj0_tdep, aj1_tdep
  public :: aj0_ptr, aj1_ptr
  public :: aj0_ptr_src, aj1_ptr_src
  public :: aj0_left, aj1_left
  public :: aj0_right, aj1_right

  public :: gamtots_tdep_ptr_type
  public :: gamtot_tdep, gamtot1_tdep, gamtot2_tdep, gamtot3_tdep
  public :: gamtot_ptr, gamtot1_ptr, gamtot2_ptr, gamtot3_ptr
  public :: gamtot_left, gamtot1_left, gamtot2_left, gamtot3_left
  public :: gamtot_right, gamtot1_right, gamtot2_right, gamtot3_right

  public :: a, b, r, ainv
  public :: r_ptr, ainv_ptr
  public :: r_left, ainv_left
  public :: r_right, ainv_right

  public :: jump
  
  ! dist fn
  complex, dimension (:,:,:), allocatable :: g, gnew, g_restart_tmp
  ! (-ntgrid:ntgrid,2, -g-layout-)

  real, dimension(:), allocatable :: kx_shift, theta0_shift
  ! (naky)

  real, dimension (:,:,:), allocatable :: vpa, vpac, vpar
!AJ
  real, dimension (:,:,:,:), allocatable :: vpa_gf
  real, dimension (:,:,:), allocatable :: vperp2_gf
  ! (-ntgrid:ntgrid,2, -g-layout-)

  real, dimension (:,:), allocatable :: vperp2
  real, dimension (:,:), allocatable, target :: aj0, aj1
  ! (-ntgrid:ntgrid, -g-layout-)

!AJ
  real, dimension (:,:,:,:,:), allocatable :: aj0_gf, aj1_gf

  integer, dimension (:), allocatable :: ittp
  ! (-ntgrid:ntgrid)

  ! collisional diagnostic of heating rate
  complex, dimension (:,:,:,:,:), allocatable :: c_rate
  ! (-ntgrid:ntgrid,ntheta0,naky,nspecies,2) replicated

#ifdef LOWFLOW
  ! v-dependent factors in low-flow terms
  real, dimension (:,:,:), allocatable :: hneoc, vparterm, wdfac, wstarfac
  ! (-ntgrid:ntgrid,2, -g-layout-)
  real, dimension (:,:,:,:,:,:), allocatable :: wdttpfac
  ! (-ntgrid:ntgrid,ntheta0,naky,negrid,nspecies,2)
#endif

  real, dimension(:), allocatable :: kx_shift_old

  integer, dimension(:), allocatable :: jump

  complex, dimension (:,:,:), allocatable :: a, b
  complex, dimension(:,:,:), allocatable, target :: r
  complex, dimension(:,:,:), allocatable, target :: r_left, r_right
  complex, dimension(:,:,:), pointer :: r_ptr
  complex, dimension(:,:,:), allocatable, target :: ainv
  complex, dimension(:,:,:), allocatable, target :: ainv_left, ainv_right
  complex, dimension(:,:,:), pointer :: ainv_ptr
  
  ! Time-dependent versions of bessel functions and gamtots, containing values for present and next
  ! time-steps. Used for flow-shear cases. --NDC 11/2018
  type :: bessel_tdep_type
      real, dimension(:,:), allocatable :: old
      real, dimension(:,:), allocatable :: new
  end type bessel_tdep_type

  type :: bessel_tdep_ptr_type
      real, dimension(:,:), pointer :: old
      real, dimension(:,:), pointer :: new
  end type bessel_tdep_ptr_type

  type(bessel_tdep_type), target :: aj0_tdep, aj1_tdep
  type(bessel_tdep_ptr_type) :: aj0_ptr, aj1_ptr
  type(bessel_tdep_ptr_type) :: aj0_ptr_src, aj1_ptr_src
  real, dimension(:,:), allocatable, target :: aj0_left, aj0_right
  real, dimension(:,:), allocatable, target :: aj1_left, aj1_right
  
  type :: gamtots_tdep_type
      real, dimension(:,:,:), allocatable :: old
      real, dimension(:,:,:), allocatable :: new
  end type gamtots_tdep_type
  
  type :: gamtots_tdep_ptr_type
      real, dimension(:,:,:), pointer :: old
      real, dimension(:,:,:), pointer :: new
  end type gamtots_tdep_ptr_type

  type(gamtots_tdep_type), target :: gamtot_tdep, gamtot1_tdep, gamtot2_tdep, gamtot3_tdep
  type(gamtots_tdep_ptr_type) :: gamtot_ptr, gamtot1_ptr, gamtot2_ptr, gamtot3_ptr
  real, dimension(:,:,:), allocatable, target :: gamtot_left, gamtot_right
  real, dimension(:,:,:), allocatable, target :: gamtot1_left, gamtot1_right
  real, dimension(:,:,:), allocatable, target :: gamtot2_left, gamtot2_right
  real, dimension(:,:,:), allocatable, target :: gamtot3_left, gamtot3_right
  
contains

  subroutine g_adjust (g, phi, bpar, facphi, facbpar, adjust_old)
    !CMR, 17/4/2012: 
    ! g_adjust transforms between representations of perturbed dist'n func'n.
    !    <delta_f> = g_wesson J0(Z) - q phi/T F_m  where <> = gyroaverage
    !        g_gs2 = g_wesson - q phi/T J0(Z) F_m - m v_||^2/T B_||/B J1(Z)/Z F_m
    ! For numerical convenience the GS2 variable g uses the form g_gs2.
    ! g_wesson (see Wesson's book, Tokamaks) is often a more convenient form:
    !     e.g. for handling collisions, calculating v-space moments in real space.
    !
    ! To transform gnew from g_gs2 to g_wesson form:
    !    call g_adjust(gnew,phinew,bparnew,fphi,fbpar)
    ! or transform from gnew from g_wesson to g_gs2 form:
    !    call g_adjust(gnew,phinew,bparnew,-fphi,-fbpar)
    !
    use species, only: spec
    use theta_grid, only: ntgrid
    use le_grids, only: anon, jend, ng2
    use gs2_layouts, only: g_lo, ik_idx, it_idx, ie_idx, is_idx, il_idx
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    real, intent (in) :: facphi, facbpar
    logical, intent(in), optional :: adjust_old

    integer :: iglo, ig, ik, it, ie, is, il
    complex :: adj
    logical :: trapped=.false.
    real, dimension(:,:), pointer :: my_aj0, my_aj1

    if(present(adjust_old)) then
        my_aj0 => aj0_ptr%old
        my_aj1 => aj1_ptr%old
    else
        my_aj0 => aj0_ptr%new
        my_aj1 => aj1_ptr%new
    end if

    if (minval(jend) .gt. ng2) trapped=.true.

    ! Flow shear cases: we can pass time dependent bessel funcs as arguments -- NDC 11/18
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)

       ! BD:  bpar == delta B_parallel / B_0(theta) so no extra factor of 
       ! 1/bmag is needed here.
       do ig = -ntgrid, ntgrid
          !<DD>Don't adjust in the forbidden region as we don't set g/h here
!CMR, 13/10/2014: attempting minor reduction in memory bandwidth by replacing 
!         if( forbid(ig,il) ) cycle with 
!         if ( trapped .and. il > jend(ig))
          if ( trapped .and. il > jend(ig)) cycle 

          adj = anon(ie)*2.0*vperp2(ig,iglo)*my_aj1(ig,iglo) &
               *bpar(ig,it,ik)*facbpar &
               + spec(is)%z*anon(ie)*phi(ig,it,ik)*my_aj0(ig,iglo) &
               /spec(is)%temp*facphi
          g(ig,1,iglo) = g(ig,1,iglo) + adj
          g(ig,2,iglo) = g(ig,2,iglo) + adj
       end do
    end do

    nullify(my_aj0, my_aj1)

  end subroutine g_adjust

  subroutine check_g_bouncepoints(g, ik,it,il,ie,is,err,tol)
    ! CMR, 3/10/2013: 
    !  This routine checks trapped particle bounce conditions: 
    !     g(thetab,1:ik,it,il,ie,is)=g(thetab,2:ik,it,il,ie,is)
    !  and flags fractional errors exceeding a threshold tolerance, tol. 
    !
    use theta_grid, only: ntgrid, bmag
    use gs2_layouts, only: g_lo, idx
    use le_grids, only: ng2, jend, al
    use mp, only: mp_abort
    
    implicit none
    integer, intent(in) :: ik, it, il, ie, is
    complex, dimension(-ntgrid:,:,g_lo%llim_proc:), intent(in) :: g
    real, intent(out):: err
    real, optional, intent(in):: tol
    real :: tolerance, dg
    integer :: iglo, ig
    logical :: started

    started=.false.
    if (present(tol)) then
       tolerance=tol 
    else 
       tolerance=1.0e-6
    endif
    iglo=idx(g_lo,ik,it,il,ie,is)
    if (iglo.lt.g_lo%llim_proc .or. iglo.gt.g_lo%ulim_proc .or. il.le.ng2) then
       return
    endif
    err=0.0
    do ig=-ntgrid,ntgrid
       ! if at a bounce point, check error on g
       if (il.eq.jend(ig) .and.  al(il)*bmag(ig).gt.0.999999) then
          dg=abs(g(ig,1,iglo)-g(ig,2,iglo))/max(abs(g(ig,1,iglo)),abs(g(ig,2,iglo)))
          if ( dg .gt. tolerance) then
             if (.not. started) then
                write(6,fmt='(T7,"ig",T17,"g(ig,1,iglo)",T43,"g(ig,2,iglo)",T63,"FracBP Error" )')
                started=.true.
             endif
             write(6,fmt='(i8, "  (",1pe11.4,",",e11.4,"),(",e11.4,",",e11.4,"), ", e11.4)') ig, g(ig,1,iglo),g(ig,2,iglo), dg
             err=max(err,dg)
          endif
       endif
    enddo
    write(6,fmt='(t5,"ik",t11,"it",t17,"il",t23,"ie",t29,"is",t33,"MaxFracBP Error")')
    write(6,fmt='(5i6,1pe12.4)')ik,it,il,ie,is,err
    write(6,*) "-----"
  end subroutine check_g_bouncepoints
end module dist_fn_arrays
