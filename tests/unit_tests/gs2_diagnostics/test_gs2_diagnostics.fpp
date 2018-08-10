#define CONCAT //

!> A regression test for the old diagnostics module. This runs 
!! a linear cyclone test case and then checks the value of diffusivity
!! against the value produced by svn revision 4242.
!!
!! This is free software released under the MIT license
!!   Written by: Edmund Highcock (edmundhighcock@users.sourceforge.net)

!module checks_mod
!  use unit_tests
!  use functional_tests
!  public checks
!  contains
!    function checks()
!      logical :: checks
!      checks = .true.
!    end function checks
!end module checks_mod
!
program test_gs2_diagnostics
  !use functional_tests
  !use checks_mod
  !call test_gs2('Linear CBC (unit test) to test new diagnostics', checks)
    use gs2_main, only: run_gs2, finish_gs2
    use unit_tests
    use file_utils, only: run_name, error_unit
    use mp, only: init_mp, mp_comm, proc0, test_driver_flag, finish_mp, mp_abort
    use gs2_diagnostics
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: finish_gs2_diagnostics_new
#endif 
    implicit none
    real :: eps

    eps = 1.0e-7
    if (precision(eps).lt. 11) eps = eps * 1000.0

    call init_mp

    test_driver_flag = .true.
    functional_test_flag = .true.

    call announce_module_test("gs2_diagnostics")

    call run_gs2(mp_comm)

    call announce_test('diffusivity')

    if( trim(run_name) .eq. 'test_gs2_diagnostics' ) then
      call process_test(diagnostics_unit_test_diffusivity(19.144575471638724, eps), 'diffusivity')
    else if( trim(run_name) .eq. 'test_gs2_diagnostics_nspec_2' ) then
      call process_test(diagnostics_unit_test_diffusivity(6.5580345754760385, eps), 'diffusivity')
    else
      call mp_abort('test_gs2_diagnostics calling invalid test name')
    end if

    call finish_gs2_diagnostics(ilast_step)
#ifdef NEW_DIAG
    call finish_gs2_diagnostics_new
#endif 

    call finish_gs2



   call close_module_test("gs2_diagnostics")

    call finish_mp

contains
  


end program test_gs2_diagnostics
