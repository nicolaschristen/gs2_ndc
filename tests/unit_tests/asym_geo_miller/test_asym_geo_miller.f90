
!> A program that runs unit tests on the generalized Miller geometry specification.
!!
!! Current tests:
!!   1. Calculates several geometric coefficients from the generalized Miller
!!      flux surface specification to verify that they have not changed.
!!
!!
!! This is free software released under the MIT license
!!   Written by: Edmund Highcock (edmundhighcock@users.sourceforge.net)
!!   Adapted by: Justin Ball (Justin.Ball@epfl.ch)
program test_asym_geo_miller
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
  call announce_module_test('asym_geo_miller')

  new_rslt = 0.

  prev_rslt = prev_result()
  if (proc0) call init_file_utils(dummy, .true., .true., .false., 'local')
  call init_theta_grid
  call get_new_results()
  if (proc0) call finish_file_utils
  call announce_test(' that generalized Miller has not changed')
  call process_test( &
    test_results(new_rslt(:,:), prev_rslt, 2e-2, 3.0e-1), &
    ' that generalized Miller has not changed')



  call close_module_test('asym_geo_miller')
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
1.33635293795397 , &
1.33099298605246 , &
1.32083826048736 , &
1.30611183977309 , &
1.28747670894275 , &
1.24134318228424 , &
1.18813663296922 , &
1.13242253596807 , &
1.07791588926715 , &
1.02730267223936 , &
0.993144652098344 , &
0.962912689286433 , &
0.929288388168257 , &
0.903719992907863 , &
0.886650771022432 , &
0.879908821295329 , &
0.878602935681074 , &
0.879908821295329 , &
0.886650771022432 , &
0.903719992907863 , &
0.929288388168257 , &
0.962912689286433 , &
0.993144652098344 , &
1.02730267223936 , &
1.07791588926715 , &
1.13242253596807 , &
1.18813663296922 , &
1.24134318228424 , &
1.28747670894275 , &
1.30611183977309 , &
1.32083826048736 , &
1.33099298605246 , &
1.33635293795397 , &
! grho
0.667160314205449 , &
0.622453410657498 , &
0.602899024411932 , &
0.599142857131351 , &
0.611574280834602 , &
0.672518190265384 , &
0.758016862439057 , &
0.846773772079592 , &
0.927113987252141 , &
0.993589597835576 , &
1.03276991427162 , &
1.06254403390402 , &
1.08915559841928 , &
1.10451464927511 , &
1.11239887596248 , &
1.11250687629734 , &
1.10731572604846 , &
1.09620267112108 , &
1.06430465152906 , &
0.995225952153163 , &
0.912842971369641 , &
0.846435645844334 , &
0.824093762206988 , &
0.828268665219625 , &
0.859902010868401 , &
0.893632206020437 , &
0.907480262096228 , &
0.892513933879149 , &
0.847424295358233 , &
0.814398370817557 , &
0.775431910758813 , &
0.732696031393377 , &
0.667160314205449 , &
! gbdrift
-1.22026576022393 , &
-0.903575231182115 , &
-0.715986054037962 , &
-0.567381993978907 , &
-0.439741221931316 , &
-0.255744871661724 , &
-0.124708414671471 , &
-0.000582491179644575 , &
0.123608710922522 , &
0.240518645458633 , &
0.318060592907738 , &
0.383609438669541 , &
0.449547484612142 , &
0.490415189625785 , &
0.509141238388929 , &
0.512174058841872 , &
0.510695941276762 , &
0.507563768070829 , &
0.501666383275845 , &
0.496378611835425 , &
0.496654206809191 , &
0.494154068044408 , &
0.475761744403759 , &
0.432180428301176 , &
0.32640137154253 , &
0.167631308915994 , &
-0.0293297254386237 , &
-0.250499996166318 , &
-0.49412729371684 , &
-0.628683358940139 , &
-0.768678318970532 , &
-0.932256763916757 , &
-1.27379990614093 , &
! cvdrift
-1.18536683564224 , &
-0.869923288583155 , &
-0.682084387502849 , &
-0.532544851605555 , &
-0.403847604859617 , &
-0.217162192271672 , &
-0.0826642473409987 , &
0.0456852732459802 , &
0.174737581937518 , &
0.296934311836471 , &
0.378545526109623 , &
0.448092082838868 , &
0.518955527557727 , &
0.563890422179597 , &
0.585414837300507 , &
0.589520590774497 , &
0.588208335794617 , &
0.584808440501307 , &
0.577675433597027 , &
0.5694188058172 , &
0.565739479682056 , &
0.558670802517906 , &
0.536537762078167 , &
0.488989798180806 , &
0.377880215653248 , &
0.214167933201222 , &
0.0128878382119454 , &
-0.211852683090303 , &
-0.458210460088527 , &
-0.593747492967421 , &
-0.734620521378639 , &
-0.898636216162353 , &
-1.23890098155924 , &
! gds2
14.6622031752583 , &
10.8127164989928 , &
8.53346330670264 , &
6.79130175992707 , &
5.40912406372689 , &
3.51112002926702 , &
2.41545202039911 , &
1.74924318137651 , &
1.32640260823814 , &
1.06051005261671 , &
0.930431128403171 , &
0.838444251569858 , &
0.750231510712446 , &
0.684927220648682 , &
0.63999108702712 , &
0.624560076602478 , &
0.626459273255338 , &
0.640175550331027 , &
0.688370233068644 , &
0.815677312595317 , &
1.03775018073447 , &
1.38814490543868 , &
1.74322890827708 , &
2.1561337867837 , &
2.70954090244443 , &
3.12969183361711 , &
3.35954589888537 , &
3.49820732874776 , &
3.74474098629792 , &
4.00201179157527 , &
4.35535153947453 , &
4.96360102477642 , &
6.76079200396937 , &
! gds21
2.44236017138935 , &
1.77974713495525 , &
1.34904900141465 , &
0.995579802561283 , &
0.71327093531472 , &
0.311117444895333 , &
0.0920218878086749 , &
0.0171273876293668 , &
0.0370579437321012 , &
0.100623772863893 , &
0.149674688180077 , &
0.183450087633031 , &
0.191119417450097 , &
0.155549817028141 , &
0.0937499641579587 , &
0.0474497484950457 , &
0.0237179999240448 , &
0.00858178266280109 , &
-0.00304635627853939 , &
-0.0285080244269559 , &
-0.11894187448759 , &
-0.306855334356314 , &
-0.508135577805063 , &
-0.742719995445783 , &
-1.04925035380732 , &
-1.26337909556127 , &
-1.33491363123388 , &
-1.28254183201944 , &
-1.17094143759241 , &
-1.12059091765804 , &
-1.08876836578566 , &
-1.09175040324748 , &
-1.18481152607719 , &
! theta
-3.14159265358979 , &
-2.97594459719933 , &
-2.85073726519902 , &
-2.72853682453103 , &
-2.60980530472059 , &
-2.38010497482605 , &
-2.15984494934298 , &
-1.94567889369262 , &
-1.73538518182181 , &
-1.52684198548675 , &
-1.37095504545554 , &
-1.2149120143528 , &
-1.00530512671691 , &
-0.790644314389682 , &
-0.565215038408503 , &
-0.384918896303878 , &
-0.257895879552839 , &
-0.130899693899575 , &
0.0654498469497874 , &
0.327249234748937 , &
0.589048622548086 , &
0.850848010347236 , &
1.0471975511966 , &
1.24354709204596 , &
1.50534647984511 , &
1.76714586764426 , &
2.0291349586737 , &
2.29074464324256 , &
2.55254403104171 , &
2.68344372494128 , &
2.81434341884086 , &
2.94524311274043 , &
3.14159265358979 &
], (/33,7/))
end function prev_result

end program test_asym_geo_miller
