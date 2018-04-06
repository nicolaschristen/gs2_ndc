module flowshear_continuous
    implicit none

    private

    public :: update_kperp2_tdep
    public :: update_aj0_tdep
    public :: update_gamtot_tdep

    contains

        ! Updating time-dependent kperp2 for cases with flow-shear --NDC 11/2017
        subroutine update_kperp2_tdep
            use kt_grids, only: aky, naky, theta0, ntheta0, kperp2, kperp2_tdep
            use theta_grid, only: ntgrid, gds21, gds22, shat
            use dist_fn_arrays, only: kx_shift
            
            implicit none

            integer :: ik, it

            do ik = 1, naky
                ! kperp is constant if ky==0
                if (aky(ik) /= 0.0) then
                    kperp2_tdep%old(:,:,ik) = kperp2_tdep%new(:,:,ik)
                    do it = 1, ntheta0
                        kperp2_tdep%new(:,it,ik) = kperp2(:,it,ik) + &
                            kx_shift(ik)*kx_shift(ik)*gds22/(shat*shat) + &
                            2.0*kx_shift(ik)*aky(ik)*gds21/shat + &
                            2.0*kx_shift(ik)*theta0(it,ik)*aky(ik)*gds22/shat
                    end do
                end if
            end do
        end subroutine update_kperp2_tdep

        ! Updating time-dependent aj0 for cases with flow-shear --NDC 11/2017
        subroutine update_aj0_tdep
            use kt_grids, only: kperp2_tdep, aky
            use species, only: spec
            use theta_grid, only: ntgrid, bmag
            use le_grids, only: energy, al
            use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
            use spfunc, only: j0
            use dist_fn_arrays, only: aj0_tdep
            
            implicit none

            integer :: iglo, ik, it, il, ie, is, ig
            real :: arg

            do iglo = g_lo%llim_proc, g_lo%ulim_proc
                ik = ik_idx(g_lo,iglo)
                ! aj0 is constant if ky==0
                if (aky(ik) /= 0.0) then
                    it = it_idx(g_lo,iglo)
                    il = il_idx(g_lo,iglo)
                    ie = ie_idx(g_lo,iglo)
                    is = is_idx(g_lo,iglo)
                    aj0_tdep%old(:,iglo) = aj0_tdep%new(:,iglo)
                    do ig = -ntgrid, ntgrid
                        arg = spec(is)%bess_fac*spec(is)%smz*sqrt(energy(ie)*al(il)/bmag(ig)*kperp2_tdep%new(ig,it,ik))
                        aj0_tdep%new(ig,iglo) = j0(arg)
                    end do
                end if
            end do
        end subroutine update_aj0_tdep

        ! Updating time-dependent gamtot for cases with flow-shear --NDC 11/2017
        subroutine update_gamtot_tdep
            use dist_fn_arrays, only: aj0_tdep
            use species, only: nspec, spec
            use theta_grid, only: ntgrid
            use kt_grids, only: ntheta0, naky, kperp2_tdep
            use le_grids, only: anon, integrate_species
            use gs2_layouts, only: g_lo, ie_idx, is_idx
            use dist_fn, only: poisfac, g0
            use dist_fn_arrays, only: gamtot_tdep
            
            implicit none
            
            complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: tot
            real, dimension (nspec) :: wgt
            integer :: iglo, ie, is, isgn
            
            gamtot_tdep%old = gamtot_tdep%new

            do iglo = g_lo%llim_proc, g_lo%ulim_proc ! why is not ulim_alloc used for g0 in dist_fn ? --NDC 11/2017
                ie = ie_idx(g_lo,iglo)
                is = is_idx(g_lo,iglo)
                do isgn = 1,2
                    g0(:,isgn,iglo) = (1.0 - aj0_tdep%new(:,iglo)**2)*anon(ie)
                end do
            end do
            wgt = spec%z*spec%z*spec%dens/spec%temp
            call integrate_species (g0, wgt, tot)
      
            gamtot_tdep%new = real(tot) + kperp2_tdep%new*poisfac
        end subroutine update_gamtot_tdep

end module flowshear_continuous
