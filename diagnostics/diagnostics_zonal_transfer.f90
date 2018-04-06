module diagnostics_zonal_transfer

!> calculate tau(kx,ky), the transfer of free energy by the non-linearity in fourier space
  implicit none

  private

!> allocate free energy transfer arrays
  public :: init_diagnostics_transfer
  public :: calculate_zonal_transfer, write_zonal_transfer

contains

  subroutine init_diagnostics_transfer(gnostics)

    use kt_grids, only: naky, ntheta0, nx, ny, akx, aky
    use theta_grid, only: ntgrid
    use species, only: nspec
    use le_grids, only: init_le_grids, nlambda, negrid
    use gs2_transforms, only: init_transforms
    use gs2_layouts, only: init_dist_fn_layouts
    use diagnostics_config, only: diagnostics_type

    implicit none

    type(diagnostics_type), intent(inout) :: gnostics

    logical :: accel
    
    accel=.false.

    !> initialize zonal_transfer to zero
    gnostics%current_results%zonal_transfer = 0.0
    
    !> call every initiation routine we need:
    !! le_grids for integration, init_gs2 transforms to initiate fourier transforms routines, etc.
    
    !> calls init_theta_grid, init_kt_grids, init_gs2_layouts, and init_species for later use
    call init_le_grids(accel,accel)
    !> will need to use fft routines
    call init_transforms(ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accel)
    call init_dist_fn_layouts(naky, ntheta0, nlambda, negrid, nspec)
!!!! accel_x =false : useless
    
  end subroutine init_diagnostics_transfer
  
  !> calculate nonlinear term appearing in GKE
  !! for now assume electrostatic
  subroutine non_linear_term (g1)

    use theta_grid, only: ntgrid, kxfac
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use gs2_layouts, only: yxf_lo
    use dist_fn_arrays, only: g, g_adjust
    use species, only: spec
    use gs2_transforms, only: transform2, inverse2
    use run_parameters, only: fphi, fapar, fbpar
    use kt_grids, only: aky, akx
    use fields_arrays, only: phi, bpar
    use constants, only: zi
    use species, only: nspec, spec
    
    implicit none
    
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent(out) :: g1
    integer :: i, j, k, isgn
!    complex, dimension(-ntgrid:ntgrid,nspec,negrid,nlambda,2) :: test_sum
    integer :: iglo, ik, it, ig, is
    
    real, dimension(:,:), allocatable :: banew, gbnew, bracketnew

    allocate (banew(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc)) ; banew = 0.
    allocate (gbnew(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc)) ; gbnew = 0.
    allocate (bracketnew(yxf_lo%ny,yxf_lo%llim_proc:yxf_lo%ulim_alloc)) ; bracketnew = 0.

    !> Form g1=i*kx*phi
    call load_kx_phi
    
    !> Transform to real space
    !! 2D fft of g1 is returned as banew
    call transform2 (g1, banew)
    
!    g1=g   !!!! g is an originally empty array from dist_fn_array, assigned with g somewhere
!    call g_adjust(g1,phi,bpar,fphi,fbpar)   !!! transforms g into h (minus boltzmann response), anyway fbpar=0
    
    !> Form g1=i*ky*g
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
!       g1(:,:,iglo)=g1(:,:,iglo)*zi*aky(ik) !!!!for all indices in velocity space, multiply by i*ky 
       g1(:,:,iglo)=g(:,:,iglo)*zi*aky(ik)
    end do
    
    !> Transform to real space    
    call transform2 (g1, gbnew) 
    
    !> Calculate (d phi /dx).(dg/dy)
    do j = yxf_lo%llim_proc, yxf_lo%ulim_proc
       do i = 1, yxf_lo%ny
          bracketnew(i,j) = banew(i,j)*gbnew(i,j)*kxfac
       end do
    end do
    
    !> Form g1=i*ky*chi
    call load_ky_phi
    
    !>Transform to real space    
    call transform2 (g1, banew)

    !Form g1=i*kx*h
!     g1=g

!     call g_adjust(g1,phi,bpar,fphi,fbpar)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
!       g1(:,:,iglo)=g1(:,:,iglo)*zi*akx(it)
       g1(:,:,iglo)=g(:,:,iglo)*zi*akx(it)
    enddo
       
    !> Transform to real space
    call transform2 (g1, gbnew)

    !> Calculate (d phi /dy).(dg/dx) and subtract from (d phi /dx).(dg/dy)
    do j = yxf_lo%llim_proc, yxf_lo%ulim_proc
       do i = 1, yxf_lo%ny
          bracketnew(i,j) = bracketnew(i,j) - banew(i,j)*gbnew(i,j)*kxfac
       end do
    end do
    
    !> Transform nonlinearity back to spectral space
    !! g1 contains Fourier coefficients associated with nonlinearity
    call inverse2 (bracketnew, g1)

!    g2=g
!    call g_adjust(g2,phi,bpar,fphi,fbpar)
    !> change variables from g to h
    call g_adjust(g,phi,bpar,fphi,fbpar)
    
    !> we have to multiply nonlinearity by hk
!    do iglo = g_lo%llim_proc, g_lo%ulim_proc
!       do isgn=1,2
!          do ig=-ntgrid,ntgrid
!             g1(ig,isgn,iglo)=g1(ig,isgn,iglo)*conjg(g2(ig,isgn,iglo))
!          end do
!       end do
!    end do

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       g1(:,:,iglo)=g1(:,:,iglo)*conjg(g(:,:,iglo))*spec(is)%temp
    end do
    
    call g_adjust(g,phi,bpar,-fphi,-fbpar)

!!!!! the sum over kx,ky of g1 should be 0! we calculate it
  
  !write(*,*) g3(10,1,g_lo%ulim_proc-10)
  
!   test_sum=0.
!   do iglo = g_lo%llim_proc, g_lo%ulim_proc
!      it = it_idx(g_lo,iglo)
!      ik = ik_idx(g_lo,iglo)
!      is = is_idx(g_lo,iglo)
!      ie = ie_idx(g_lo,iglo)
!      il = il_idx(g_lo,iglo)
!      do ig = -ntgrid,ntgrid
!         do isgn=1,2
!            test_sum(ig,is,ie,il,isgn) = test_sum(ig,is,ie,il,isgn) + g1(ig,isgn,iglo)
!         end do
!      end do
!   end do
!		write(*,*) abs(test_sum(1,1,2,20,1))/maxval(abs(g3))   !!!! should be zero

!     do iglo = g_lo%llim_proc, g_lo%ulim_proc
!        is = is_idx(g_lo,iglo)
!        do isgn=1,2
!           do ig=-ntgrid,ntgrid
!              g1(ig,isgn,iglo)=spec(is)%temp*g3(ig,isgn,iglo)
!           end do
!        end do
!     end do

    deallocate (banew, gbnew, bracketnew)
    
  contains
    
    subroutine load_kx_phi
      
      use dist_fn_arrays, only: aj0  !!!bessel function from the fourier coeff of average
      complex :: fac

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         it = it_idx(g_lo,iglo)
         ik = ik_idx(g_lo,iglo)
         do ig = -ntgrid, ntgrid
            fac = zi*akx(it)*aj0(ig,iglo)*phi(ig,it,ik)
            g1(ig,1,iglo) = fac
            g1(ig,2,iglo) = fac
         end do
      end do

    end subroutine load_kx_phi

    subroutine load_ky_phi

      use dist_fn_arrays, only: aj0
      complex :: fac

      do iglo = g_lo%llim_proc, g_lo%ulim_proc
         it = it_idx(g_lo,iglo)
         ik = ik_idx(g_lo,iglo)
         do ig = -ntgrid, ntgrid
            fac = zi*aky(ik)*aj0(ig,iglo)*phi(ig,it,ik)
            g1(ig,1,iglo) = fac
            g1(ig,2,iglo) = fac
         end do
      end do

    end subroutine load_ky_phi

  end subroutine non_linear_term

!!!!!!!!!!!!!! integration over velocity space
  
!   subroutine integrate (g, total)  !!!!! integrates g with maxwellian weights, output in total array
    
!     use le_grids, only: integrate_moment  !!!! Not public, until told so in le_grids!
!     use theta_grid, only: ntgrid
!     use gs2_layouts, only: g_lo
!     use species, only:nspec
    
!     implicit none
    
!     complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g
!     complex, dimension (-ntgrid:,:,:,:), intent (out) :: total
    
!     call integrate_moment(g, total, all=.true.)  !!!! can only be called once w and wl are assigned, after init_le_grids
    
!   end subroutine integrate

  !> average over theta
  subroutine average_theta(a, axb)
    
    use theta_grid, only: ntgrid, delthet, jacob
    
    implicit none
    
    complex, dimension (-ntgrid:), intent (in) :: a
    complex, intent (out) :: axb
    
    real, dimension (-ntgrid:ntgrid) :: wgt
    real :: anorm
    
    wgt = delthet(-ntgrid:ntgrid)*jacob(-ntgrid:ntgrid)
    anorm = sum(wgt)
    
    axb = sum(a(-ntgrid:ntgrid)*wgt)/anorm
    
  end subroutine average_theta

  !> returns tau, the transfer of free energy as a function of (kx,ky)
  subroutine total_term(tau)
    
    use mp, only: proc0, broadcast
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, it_idx
    use species, only: nspec, spec
    use kt_grids, only: ntheta0, naky
    use le_grids, only: integrate_moment

    implicit none
    
    complex, dimension(:,:), intent(out) :: tau

    complex, dimension (:,:,:), allocatable :: g0
    complex, dimension (:,:,:), allocatable :: inter
    complex, dimension (:,:,:,:), allocatable :: inter_s

    integer :: ik, it, is, iglo
    
    allocate (g0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (inter(-ntgrid:ntgrid,ntheta0,naky)) ; inter = 0.
    allocate (inter_s(-ntgrid:ntgrid,ntheta0,naky,nspec)) ; inter_s = 0.

    ! obtain the nonlinear term
    call non_linear_term(g0)
    
    !> integrate over velocities
    !! can only be called after w and wl are assigned in init_le_grids
!    call integrate(g0,inter_s)
    call integrate_moment(g0, inter_s)

    if (proc0) then
       !> sum over species
       do is=1,nspec
          inter = inter + inter_s(:,:,:,is)
       end do
    
       !> integrate over theta, must be done after velocity space integration
       do ik = 1, naky
          do it = 1, ntheta0
             call average_theta(inter(:,it,ik),tau(it,ik))
          end do
       end do
    end if
    
    do ik = 1, naky
       call broadcast (tau(:,ik))
    end do

    deallocate (g0, inter, inter_s)
    
  end subroutine total_term
  
  subroutine write_zonal_transfer(gnostics)
    
!    use diagnostics_create_and_write, only: create_and_write_distributed_fieldlike_variable
    use diagnostics_create_and_write, only: create_and_write_variable
    use diagnostics_dimensions, only: dim_string
    use fields_parallelization, only: field_k_local
    use kt_grids, only: ntheta0, naky
    use diagnostics_config, only: diagnostics_type
    
    implicit none
    
    type(diagnostics_type), intent(in) :: gnostics
    
!    call create_and_write_distributed_fieldlike_variable(gnostics, gnostics%rtype, "zonal_transfer", &
!         dim_string([gnostics%dims%ri,gnostics%dims%kx,gnostics%dims%ky,gnostics%dims%time]), &
!         "Time rate of change of free energy due to nonlinear transfer, as a function of kx, ky and time", &
!         "Tr2*rhor*c*ns*vth6 / Ba2*e*a3", gnostics%current_results%zonal_transfer)
! check gnostics%rtype
    call create_and_write_variable(gnostics, gnostics%rtype, "zonal_transfer", &
         dim_string([gnostics%dims%ri,gnostics%dims%kx,gnostics%dims%ky,gnostics%dims%time]), &
         "Time rate of change of free energy due to nonlinear transfer, as a function of kx, ky and time", &
         "Tr2*rhor*c*ns*vth6 / Ba2*e*a3", gnostics%current_results%zonal_transfer)
    
  end subroutine write_zonal_transfer
  
  subroutine calculate_zonal_transfer(gnostics)

    use diagnostics_config, only: diagnostics_type
    use kt_grids, only: ntheta0, naky

    implicit none
    
    type(diagnostics_type), intent(inout) :: gnostics
    
    complex, dimension (:,:), allocatable ::  zonal_transfer

    allocate (zonal_transfer(ntheta0,naky)) ; zonal_transfer = 0.

    call total_term(zonal_transfer)
    gnostics%current_results%zonal_transfer = zonal_transfer

    deallocate (zonal_transfer)
    
  end subroutine calculate_zonal_transfer
  
end module diagnostics_zonal_transfer













