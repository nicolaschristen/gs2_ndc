
!> A program that runs unit tests on the generalized elongation geometry specification.
!!
!! Current tests:
!!   1. Calculates several geometric coefficients from the generalized
!!      elongation flux surface specification to verify that they have not changed.
!!
!!
!! This is free software released under the MIT license
!!   Written by: Edmund Highcock (edmundhighcock@users.sourceforge.net)
!!   Adapted by: Justin Ball (Justin.Ball@epfl.ch)
program test_asym_geo_genElong
  use unit_tests
  use theta_grid, only: init_theta_grid, finish_theta_grid
  use theta_grid, only: initialized, theta
  use theta_grid, only: bmag, grho, gbdrift, cvdrift, gds2, gds21
  use egrid
  use mp, only: init_mp, finish_mp, proc0
  use file_utils, only: init_file_utils, finish_file_utils
  use species, only: init_species, nspec, spec
  use constants, only: pi
  use kt_grids, only: naky, ntheta0, init_kt_grids
  use gs2_layouts, only: init_gs2_layouts, g_lo, ie_idx
  implicit none
  real :: eps
  logical :: dummy
  real, dimension(33, 7) :: new_rslt
  real, dimension(33, 7) :: prev_rslt


  eps = 1.0e-7

  if (precision(eps) .lt. 11) eps = eps * 1000.0

  call init_mp
  call announce_module_test('asym_geo_genElong')

  new_rslt = 0.

  prev_rslt = prev_result()
  if (proc0) call init_file_utils(dummy, .true., .true., .false., 'local')
  call init_theta_grid
  call get_new_results()
  if (proc0) call finish_file_utils
  call announce_test(' that generalized elongation has not changed')
  call process_test( &
    test_results(new_rslt(:,:), prev_rslt, 2e-2, 3.0e-1), &
    ' that generalized elongation has not changed')



  call close_module_test('asym_geo_genElong')
  call finish_mp

contains
  function test_results(rslt1, rslt2, epsin, bmagerror)
    real, dimension(:,:) :: rslt1, rslt2
    logical :: test_results
    real :: epsin, bmagerror
    test_results = .true.
    call announce_check('bmag')
    call process_check(test_results, &
      agrees_with(rslt1(:,1), rslt2(:,1), bmagerror), 'bmag')
    call announce_check('grho')
    call process_check(test_results, &
      agrees_with(rslt1(:,2), rslt2(:,2), bmagerror), 'grho')
    call announce_check('gbdrift')
    call process_check(test_results, &
      agrees_with(rslt1(:,3), rslt2(:,3), epsin), 'gbdrift')
    call announce_check('cvdrift')
    call process_check(test_results, &
      agrees_with(rslt1(:,4), rslt2(:,4), epsin), 'cvdrift')

  end function test_results
  subroutine get_new_results()
    call init_theta_grid
    new_rslt(:,1) = spl_fit(bmag)
    new_rslt(:,2) = spl_fit(grho)
    new_rslt(:,3) = spl_fit(gbdrift)
    new_rslt(:,4) = spl_fit(cvdrift)
    new_rslt(:,5) = spl_fit(gds2)
    new_rslt(:,6) = spl_fit(gds21)
    call finish_theta_grid
  end subroutine get_new_results

  function spl_fit(var)
    use splines
    type(spline) :: spl
    real, dimension(:), intent(in) :: var
    real, dimension(33) :: spl_fit
    integer :: i
    call new_spline(size(var), theta, var, spl)
    do i = 1,33
      spl_fit(i) = splint(prev_rslt(i,7), spl)
    end do
    call delete_spline(spl)
  end function spl_fit

  function prev_result()
    real, dimension (33,7) :: prev_result

prev_result = &
reshape([&
! Bmag
1.32368869672564 , &
1.31798827414809 , &
1.30825120718884 , &
1.27532532879558 , &
1.22777557223888 , &
1.17390715913496 , &
1.11981110822868 , &
1.09280576438369 , &
1.04432442416983 , &
1.01773590055107 , &
0.976508960236828 , &
0.947706555899248 , &
0.924530413332328 , &
0.902362529026807 , &
0.880390117388593 , &
0.862352694057515 , &
0.856413201749074 , &
0.862352694057515 , &
0.880390117388593 , &
0.902362529026807 , &
0.924530413332328 , &
0.947706555899248 , &
0.976508960236828 , &
1.01773590055107 , &
1.04432442416983 , &
1.09280576438369 , &
1.11981110822868 , &
1.17390715913496 , &
1.22777557223888 , &
1.27532532879558 , &
1.30825120718884 , &
1.31798827414809 , &
1.32368869672564 , &
! grho
0.695752932359141 , &
0.679743546942874 , &
0.665290849224824 , &
0.638215466787107 , &
0.625951141739045 , &
0.6317372217417 , &
0.654923343907442 , &
0.674840731589945 , &
0.734063150575245 , &
0.786935864471537 , &
0.913854803982456 , &
1.028697853857 , &
1.09090158671961 , &
1.09466759496339 , &
1.04086707553688 , &
0.932965441588362 , &
0.834045047750794 , &
0.773667526343015 , &
0.758286965001089 , &
0.776515307346315 , &
0.805646664935207 , &
0.838906615136337 , &
0.878176362317213 , &
0.922426047027385 , &
0.939915373712938 , &
0.94466211576692 , &
0.931564049901916 , &
0.876952988329893 , &
0.80271164404367 , &
0.743827886276401 , &
0.717346444734948 , &
0.711408170694626 , &
0.695752932359142 , &
! gbdrift
-1.04371159230805 , &
-0.855378824568088 , &
-0.7058894873864 , &
-0.437880913229861 , &
-0.232122890214455 , &
-0.0733761775546126 , &
0.0579139987910111 , &
0.117921122608413 , &
0.212863975448319 , &
0.253466496387732 , &
0.29318729862726 , &
0.300683488002944 , &
0.316776871421963 , &
0.382863671050543 , &
0.486200181764614 , &
0.592104488965358 , &
0.672274461993619 , &
0.732602412541132 , &
0.775234874469865 , &
0.799287983132661 , &
0.813349163942332 , &
0.818689957822904 , &
0.808431190230252 , &
0.755345610674448 , &
0.697334406301491 , &
0.559021597330803 , &
0.473867585390355 , &
0.291759282835162 , &
0.0531671515411274 , &
-0.278086984866032 , &
-0.659858742748837 , &
-0.85074933001604 , &
-1.2161678888773 , &
! cvdrift
-1.01012083058685 , &
-0.821176230730599 , &
-0.671215534880397 , &
-0.401775641427978 , &
-0.193234298762184 , &
-0.0307457878137321 , &
0.104898775479441 , &
0.167337281948511 , &
0.267253000045124 , &
0.311042491745759 , &
0.356621190985714 , &
0.368895273778157 , &
0.38869404584708 , &
0.457885690159346 , &
0.564264450723993 , &
0.672366451431 , &
0.752389937201757 , &
0.81080552513501 , &
0.850129620530861 , &
0.8708495980906 , &
0.881871560134288 , &
0.884282189093507 , &
0.8706741186046 , &
0.813285959563862 , &
0.752744298298723 , &
0.610090469368697 , &
0.522556174429508 , &
0.335714589219273 , &
0.0927792891032188 , &
-0.241638119310802 , &
-0.62503518629388 , &
-0.816262830188617 , &
-1.1825771271561 , &
! gds2
11.61053442938 , &
10.1463705432301 , &
8.90774729066408 , &
6.5750397256193 , &
4.8714722841454 , &
3.7286790176349 , &
2.91505642150912 , &
2.565011573902 , &
1.98172353570059 , &
1.68251163265815 , &
1.25493595646217 , &
0.981380858458517 , &
0.795871931816881 , &
0.708437083918482 , &
0.728243869463019 , &
0.845888159143225 , &
1.01752601778955 , &
1.22521306242834 , &
1.47733376987732 , &
1.71990881594994 , &
1.95648174099949 , &
2.20579867659025 , &
2.50482172274685 , &
2.87070481848318 , &
3.04160385480529 , &
3.2174729183743 , &
3.28586609101937 , &
3.5993185428902 , &
4.48807927912119 , &
6.12784985487837 , &
8.21970646300642 , &
9.2524611193917 , &
11.0381340245024 , &
! gds21
2.54071339052167 , &
2.21330915023552 , &
1.92873791239689 , &
1.35525802987599 , &
0.877422149335528 , &
0.523426058332652 , &
0.259911908210243 , &
0.145146044000706 , &
-0.0559278863699604 , &
-0.173125417532216 , &
-0.35972568870552 , &
-0.409350340944689 , &
-0.285775764717725 , &
-0.089377530711544 , &
0.0377511402705118 , &
0.0245513184831693 , &
-0.0831956490056685 , &
-0.238103195986338 , &
-0.450831487526542 , &
-0.664023349653258 , &
-0.864393577096791 , &
-1.06121843559464 , &
-1.28009096651728 , &
-1.51904982661843 , &
-1.61049487463941 , &
-1.6390438059111 , &
-1.59628479000698 , &
-1.49572904423835 , &
-1.51462490810205 , &
-1.71936658689947 , &
-2.04179859339055 , &
-2.20032197185845 , &
-2.4496130366576 , &
! theta
-3.14159265358979 , &
-2.99845364795748 , &
-2.88549614323487 , &
-2.66870648128798 , &
-2.46106474089249 , &
-2.26364597731672 , &
-2.07026763099018 , &
-1.96905719664962 , &
-1.76714586764426 , &
-1.63624617374468 , &
-1.37444678594553 , &
-1.11264739814639 , &
-0.850848010347236 , &
-0.589048622548086 , &
-0.327249234748937 , &
-0.0654498469497874 , &
0.149793533673363 , &
0.340678915178904 , &
0.528442809259481 , &
0.674334391009427 , &
0.793864599963748 , &
0.905618023099863 , &
1.03655091951204 , &
1.2205957491219 , &
1.34167062621101 , &
1.5707963267949 , &
1.70169602069447 , &
1.96349540849362 , &
2.22529479629277 , &
2.48709418409192 , &
2.74889357189107 , &
2.87979326579064 , &
3.14159265358979 &
], (/33,7/))
end function prev_result

end program test_asym_geo_genElong
