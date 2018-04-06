module collisional_heating

  implicit none
  private
  
  public :: init_collisional, finish_collisional
  public :: write_collisional, calculate_collisional
  
  complex, dimension (:,:,:), allocatable :: g0   !!!! g0 will be used in calculate_collisional,
  complex, dimension (:,:,:,:), allocatable :: tot              !!!! they must be dummy arguments...
  
  complex, dimension(:,:), allocatable :: coll_heating, coll_heating_2 !!! we will use the hk variable in diagnostics_heating
  
contains
  
  subroutine init_collisional(gnostics) 
    use kt_grids, only: naky, ntheta0, nx, ny, akx, aky
    use theta_grid, only: ntgrid
    use species, only: nspec
    use diagnostics_config, only: diagnostics_type
    use le_grids, only: init_le_grids, nlambda, negrid
    use gs2_transforms, only: init_transforms
    use gs2_layouts, only: g_lo, init_dist_fn_layouts
    use gs2_layouts, only: yxf_lo
    use hyper, only: init_hyper
    
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    logical :: accel
    
    accel=.false.
    
!!!calls every initiation routine we need: le_grids for integration, init_gs2 transforms to initiate fourier transforms routines...
    
    call init_le_grids(accel,accel)  !!!! benefit: calls init_theta_grid, init_kt_grids, init_gs2_layouts and init_species for later use
    call init_transforms(ntgrid, naky, ntheta0, nlambda, negrid, nspec, nx, ny, accel)  !!! in order to use fft transform routines
    call init_dist_fn_layouts(naky, ntheta0, nlambda, negrid, nspec)
!!!! accel_x =false : useless
    
!!!!! now need to initiate and assign hypervisc_filter
    call init_hyper
    
!!!!! allocate and initialize arrays
    allocate(tot(-ntgrid:ntgrid,ntheta0,naky,nspec))
    allocate(g0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    
    allocate (coll_heating(ntheta0,naky))
    coll_heating=0.0
    
    allocate (coll_heating_2(ntheta0,naky))
    coll_heating_2=0.0
    
  end subroutine init_collisional
  
  subroutine finish_collisional
    deallocate(coll_heating)
    deallocate(coll_heating_2)
    deallocate(tot)
    deallocate(g0)
  end subroutine finish_collisional
  
  subroutine calculate_collisional(gnostics)
    
    use fields_arrays, only: phinew
    use dist_fn_arrays, only: g, gnew
    use diagnostics_config, only: diagnostics_type
    use dist_fn_arrays, only: c_rate, g_adjust
    use theta_grid, only: ntgrid, delthet, jacob
    use kt_grids, only: naky, ntheta0, akx, aky
    use species, only: spec, nspec
    use le_grids, only: integrate_moment
    use gs2_layouts, only: g_lo, it_idx, ik_idx, is_idx
    use gs2_time, only: code_dt
    use hyper, only: hyper_diff, hypervisc_filter
    use le_derivatives, only: vspace_derivatives
    use fields_arrays, only: phi, bpar, phinew, bparnew
    use run_parameters, only: fphi, fbpar

    implicit none

    type(diagnostics_type), intent(inout) :: gnostics
    real, dimension(-ntgrid:ntgrid) :: wgt !! weights for theta integration
    complex, dimension(ntheta0,naky) :: hyper_part, hyper_part_2, colls_part   !!! two parts of c_rate array: collisions and hypervisocity
    !complex, dimension(-ntgrid:ntgrid,ntheta0,naky) :: tot
    integer :: is, ik, it, ig, isgn, iglo
    complex :: dgdt_hypervisc, havg
    real :: fac, fac2
    real, dimension(1:nspec) :: weights 

    hyper_part = 0.
    hyper_part_2 = 0.
    colls_part = 0.

    wgt = delthet*jacob
    wgt = wgt/sum(wgt)
    g0=g
    weights=1.
  
    call g_adjust(g0,phi,bpar,fphi,fbpar)   !!! transforms g into h (minus boltzmann response), anyway fbpar=0

    !call vspace_derivatives(gnew,g,g0,phi, bpar, phinew, bparnew)  !!!!! assigns values to c_rate
    call hyper_diff(g0,phi) !!!! assigns values to hypervisc_filter
       
    do is = 1, nspec
       do ik = 1, naky
          fac = 0.5
          if (aky(ik) < epsilon(0.)) fac = 1.0
          do it = 1, ntheta0
             do ig= -ntgrid,ntgrid
                colls_part(it,ik) = colls_part(it,ik) + real(c_rate(ig,it,ik,is,3))*fac*wgt(ig)*spec(is)%temp*spec(is)%dens
                hyper_part_2(it,ik) = hyper_part_2(it,ik) + real(c_rate(ig,it,ik,is,2))*fac*wgt(ig)*spec(is)%temp*spec(is)%dens
             end do
          end do
       end do
    end do

    do iglo=g_lo%llim_proc, g_lo%ulim_proc
      is = is_idx(g_lo, iglo)
      it = it_idx(g_lo, iglo)
      ik = ik_idx(g_lo, iglo)
      do isgn=1,2
        do ig=-ntgrid, ntgrid-1
                
          havg = 0.25*(g(ig,isgn,iglo)+ g(ig+1,isgn,iglo)+ gnew(ig,isgn,iglo)+ gnew(ig+1,isgn,iglo)) 

          dgdt_hypervisc = 0.5*((1.0-1./hypervisc_filter(ig,it,ik))*gnew(ig,isgn,iglo) &
                     + (1.0-1./hypervisc_filter(ig+1,it,ik))*gnew(ig+1,isgn,iglo))/code_dt

          g0(ig,isgn,iglo) = spec(is)%dens*spec(is)%temp*conjg(havg)*dgdt_hypervisc
					
        end do
      end do
    end do

    !call integrate_species_original(g0, weights, tot)    !!!!! integrate species_original is dubious...
    call integrate_moment(g0, tot)
    !write(*,*) tot
    do ik = 1, naky
       fac2 = 0.5
       if (aky(ik) < epsilon(0.0)) fac2 = 1.0
       do it = 1, ntheta0
          do is = 1, nspec
             do ig = -ntgrid, ntgrid-1
                hyper_part(it,ik) = hyper_part(it,ik) + real(tot(ig,it,ik,is))*wgt(ig)*fac2
             end do
          end do
       end do
    end do
    
    coll_heating= colls_part + hyper_part
    coll_heating_2 = colls_part + hyper_part_2
    
  end subroutine calculate_collisional
  
  subroutine write_collisional(gnostics)
    
    use diagnostics_create_and_write, only: create_and_write_distributed_fieldlike_variable
    use diagnostics_dimensions, only: dim_string
    use fields_parallelization, only: field_k_local
    use kt_grids, only: ntheta0, naky
    use diagnostics_config, only: diagnostics_type
    
    implicit none
    
    type(diagnostics_type), intent(in) :: gnostics
    
    call create_and_write_distributed_fieldlike_variable(gnostics, gnostics%rtype, "collisional_heating", &
         dim_string([gnostics%dims%ri,gnostics%dims%kx,gnostics%dims%ky,gnostics%dims%time]), &
         "Collisional and hyper viscous rate of loss of free energy for each mode and each t", &
         "unit", coll_heating)
    
    call create_and_write_distributed_fieldlike_variable(gnostics, gnostics%rtype, "collisional_heating_2", &
         dim_string([gnostics%dims%ri,gnostics%dims%kx,gnostics%dims%ky,gnostics%dims%time]), &
         "Collisional and hyper viscous rate of loss of free energy for each mode and each t, other method", &
         "unit", coll_heating_2)
    
  end subroutine write_collisional

end module collisional_heating
