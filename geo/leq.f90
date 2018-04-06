module leq
  implicit none

  private
  
  public :: leq_init, leqin, gradient, eqitem, bgradient, leq_finish
  public :: invR, Rpos, Zpos, diameter, btori, dbtori, qfun, pfun
  public :: dpfun, betafun, psi, rcenter, dpdrhofun

  integer :: nr, nt

  real, allocatable, dimension (:)     :: eqpsi, fp, beta, pressure
  real, allocatable, dimension (:,:)   :: R_psi, Z_psi
  real, allocatable, dimension (:,:,:) :: drm, dzm, dbtm, dpm, dtm
  real, allocatable, dimension (:,:,:) :: dpcart, dtcart, dbtcart
  real, allocatable, dimension (:,:,:) :: dpbish, dtbish, dbtbish

  real :: beta_0
  real :: thetaShift
 
  type :: flux_surface
     real :: Rmaj, Rgeo, r, dr, aSurf, sHorz, sVert, delm, deln, delmp, delnp, thm, thn, q, shat
     integer :: nt, geoType, mMode, nMode
  end type flux_surface

  type (flux_surface) :: surf

contains
  subroutine leqin(geoType, Rmaj, Rgeo, r, dr, aSurf, sHorz, sVert, mMode, nMode, delm, deln, delmp, delnp, thm, thn, q, shat, nt_used, thShift)
    implicit none
    real, intent(in) :: Rmaj, Rgeo, r, dr, aSurf, sHorz, sVert, delm, deln, delmp, delnp, thm, thn, q, shat
    real, intent(in) :: thShift
    integer, intent(in) :: nt_used, geoType, mMode, nMode

!R0   3.1601153788426086        3.0818707511680938        1.3431053933541011
!      0.10381317967192064        3.3317843820181277E-002   8.0267210436786265E-002
!  0.37777777777777777        1.0000000000000000E-003  -4.6851042901036040E-002
!   1.0642542534450539       0.32875578050587551        0.0000000000000000        0.0000000000000000                4
!f R0   3.1601153788426086        3.0818707511680938        1.3431053933541011
!        2.5126356633383937E-316   3.3317843820181277E-002   4.9406564584124654E-324
!   2.5126356633383937E-316   1.0000000000000000E-003   9.1922790858218116E-317
!   1.0642542534450539        2.3564531823208357E-310   0.0000000000000000        0.0000000000000000                4
!f R0   3.1601153788426086        3.0818707511680938        1.3431053933541011
!        2.2724233178454796E-316   3.3317843820181277E-002   4.9406564584124654E-324
!   2.2724233178454796E-316   1.0000000000000000E-003   9.1922790858218116E-317   
!   1.0642542534450539        2.3203788153401778E-310   0.0000000000000000        0.0000000000000000                4

! kp, dp, r, s, qs
    !write (*,*) 'leq inputs', R0,Ra, k, kp, d, dp, r, dr, s, qq, qs, a, ap, nt_used
    surf%geoType = geoType
    surf%Rmaj = Rmaj
    surf%Rgeo = Rgeo
    surf%r = r
    surf%dr = dr 
    surf%aSurf = aSurf
    surf%sHorz = sHorz
    surf%sVert = sVert
    surf%mMode = mMode
    surf%nMode = nMode
    surf%delm = delm
    surf%deln = deln
    surf%delmp = delmp
    surf%delnp = delnp
    surf%thm = thm
    surf%thn = thn
    surf%q = q
    surf%shat = shat
 
    beta_0 = 1.

    thetaShift = thShift;

    nr = 3
    nt = nt_used
    if(.not.allocated(beta)) call alloc_arrays(3, nt)
    surf%nt = nt
    call leq_init
  end subroutine leqin

  subroutine alloc_arrays(nr, nt)
    integer, intent(in) :: nr, nt

    allocate(eqpsi(nr), fp(nr), beta(nr), pressure(nr))
    allocate(R_psi(nr, -nt:nt), Z_psi(nr, -nt:nt))
    allocate(drm(nr, -nt:nt, 2), dzm(nr, -nt:nt, 2), dbtm(nr, -nt:nt, 2), &
         dpm(nr, -nt:nt, 2), dtm(nr, -nt:nt, 2))
    allocate(dpcart(nr, -nt:nt, 2), dtcart(nr, -nt:nt, 2), dbtcart(nr, -nt:nt,2))
    allocate(dpbish(nr, -nt:nt, 2), dtbish(nr, -nt:nt, 2), dbtbish(nr, -nt:nt,2))

  end subroutine alloc_arrays

  subroutine dealloc_arrays
    implicit none
    if(allocated(eqpsi)) deallocate(eqpsi, fp, beta, pressure)
    if(allocated(R_psi)) deallocate(R_psi, Z_psi)
    if(allocated(drm)) deallocate(drm,dzm,dbtm,dpm,dtm)
    if(allocated(dpcart)) deallocate(dpcart,dtcart,dbtcart)
    if(allocated(dpbish)) deallocate(dpbish,dtbish,dbtbish)
  end subroutine dealloc_arrays

  subroutine leq_finish
    implicit none
    call dealloc_arrays
  end subroutine leq_finish

  subroutine leq_init
    use constants, only: pi=>pi
    implicit none
    real, dimension(nr, -nt:nt) :: eqpsi1, eqth, eqbtor
    real :: dr(3)
    real :: t, r
    integer :: i, j
   
    dr(1) = -surf%dr
    dr(2) = 0.
    dr(3) = surf%dr
   
    ! All calculation must be extended over the full [-pi,pi] range in order to
    ! support up-down asymmetric geometries (see section 3.2.2 of Ball MIT
    ! Masters thesis or section 2.2 of "GS2 analytic geometry specification") 
    do j=-nt,nt
       do i=1,nr
          r = surf%r + dr(i)
          t = (j)*pi/real(nt-1)
          R_psi(i,j) = Rpos(r, t) 
          Z_psi(i,j) = Zpos(r, t)
          eqth(i,j) = t
          eqpsi1(i,j) = 1 + dr(i)
          eqbtor(i,j) = surf%Rgeo/R_psi(i,j)
       enddo
    enddo

    do i=1,nr
       pressure(i) = -dr(i)
    enddo

    eqpsi(:) = eqpsi1(:,1)

    call derm(eqth,   dtm,  'T')
    call derm(R_psi,  drm,  'E')
    call derm(Z_psi,  dzm,  'O')
    call derm(eqbtor, dbtm, 'E')
    call derm(eqpsi1, dpm,  'E')
    
! diagnostics
!      do j=1,nt
!         do i=1,nr
!            write(*,*) i,j
!            write(*,100) drm(i,j,1),drm(i,j,2),R_psi(i,j)
!            write(*,101) dzm(i,j,1),dzm(i,j,2),Z_psi(i,j)
!            write(*,102) dtm(i,j,1),dtm(i,j,2),eqth(i,j)
!         enddo
!      enddo
! 100  format('(gr R)1 ',g10.4,' (gr R)2 ',g10.4,' R ',g10.4)
! 101  format('(gr Z)1 ',g10.4,' (gr Z)2 ',g10.4,' Z ',g10.4)
! 102  format('(gr t)1 ',g10.4,' (gr t)2 ',g10.4,' t ',g10.4)
!      write(*,*) nr, nt
!      stop

! below is actually grad(rho) instead of grad(psi),
! and 'cartesian' refers to (R,Z) coordinates -- MAB
! grad(psi) in cartesian form 
    call eqdcart(dpm, dpcart)
! grad(psi) in Bishop form 
    call eqdbish(dpcart, dpbish)

! grad(BT) in cartesian form
    call eqdcart(dbtm, dbtcart)
! grad(BT) in Bishop form
    call eqdbish(dbtcart, dbtbish)

! grad(theta) in cartesian form
    call eqdcart(dtm, dtcart)
! grad(theta) in Bishop form
    call eqdbish(dtcart, dtbish)

  end subroutine leq_init

  subroutine derm(f, dfm, char)
    use constants, only: pi=>pi
    implicit none
    real, dimension(:,:), intent(in) :: f(:,-nt:)
    real, dimension(:,:,:), intent(out) :: dfm(:,-nt:,:)
    character(1), intent(in) :: char
    integer :: i, j
    
    i=1
    dfm(i,:,1) = -0.5*(3*f(i,:)-4*f(i+1,:)+f(i+2,:))         
    
    i=nr
    dfm(i,:,1) = 0.5*(3*f(i,:)-4*f(i-1,:)+f(i-2,:))

    select case (char)
    case ('E')
       j=-nt
       dfm(:,j,2) = 0.5*(f(:,j+1)-f(:,j+1))

       j=nt
       dfm(:,j,2) = -0.5*(f(:,j-1)-f(:,j-1))
    case ('O')
       j=-nt
       dfm(:,j,2) = 0.5*(f(:,j+1)+f(:,j+1))

       j=nt
       dfm(:,j,2) = -0.5*(f(:,j-1)+f(:,j-1))
    case ('T')
       j=-nt
       dfm(:,j,2) = f(:,j+1) + pi

       j=nt
       dfm(:,j,2) = pi - f(:,j-1)
    end select

    do i=2,nr-1
       dfm(i,:,1)=0.5*(f(i+1,:)-f(i-1,:))
    enddo

    do j=-nt+1,nt-1
       dfm(:,j,2)=0.5*(f(:,j+1)-f(:,j-1))
    enddo

  end subroutine derm

  subroutine gradient(theta, grad, char, rp, nth_used, ntm)
    use splines, only: inter_d_cspl
    implicit none
    integer, intent(in) :: nth_used, ntm
    character(1), intent(in) :: char
    real, dimension(-ntm:), intent(in) :: theta
    real, dimension(-ntm:,:), intent(out) :: grad
    real, intent(in) :: rp
    real :: tmp(2), aa(1), daa(1), rpt(1)
    real, dimension(nr,-nt:nt,2) :: dcart
    integer :: i
    
    select case(char)
    case('D')  ! diagnostic 
       dcart = dbtcart
    case('P') 
       dcart = dpcart
    case('R') 
       dcart = dpcart  ! dpcart is correct for 'R'
    case('T')
       dcart = dtcart
    end select
    
    do i=-nth_used,-1
       call eqitem(theta(i), dcart(:,:,1), tmp(1), 'R')
       call eqitem(theta(i), dcart(:,:,2), tmp(2), 'Z')
       if(char == 'T') then
          grad(i,1)=-tmp(1)
          grad(i,2)=-tmp(2)
       else
          grad(i,1)=tmp(1)
          grad(i,2)=tmp(2)
       endif
    enddo

    do i=0,nth_used
       call eqitem(theta(i), dcart(:,:,1), tmp(1), 'R')
       call eqitem(theta(i), dcart(:,:,2), tmp(2), 'Z')
       grad(i,1)=tmp(1)
       grad(i,2)=tmp(2)
    enddo

!     to get grad(pressure), multiply grad(psi) by dpressure/dpsi

    if(char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(nr, eqpsi, pressure, 1, rpt, aa, daa)
       do i=-nth_used, nth_used
          grad(i,1)=grad(i,1)*daa(1)*0.5*beta_0
          grad(i,2)=grad(i,2)*daa(1)*0.5*beta_0
       enddo
    endif
  end subroutine gradient

  subroutine bgradient(theta, grad, char, rp, nth_used, ntm)
    use splines, only: inter_d_cspl
    implicit none
    integer, intent(in) :: nth_used, ntm
    character(1), intent(in) :: char
    real, dimension(-ntm:), intent(in) :: theta
    real, dimension(-ntm:,:), intent(out) :: grad
    real, intent(in) :: rp
    real :: aa(1), daa(1), rpt(1)
    real, dimension(nr,-nt:nt,2) ::  dbish
    integer :: i

    select case(char)
    case('D')  ! diagnostic
       dbish = dbtbish
    case('P') 
       dbish = dpbish
    case('R') 
       dbish = dpbish  ! dpcart is correct for 'R'
    case('T')
       dbish = dtbish
    end select

    do i=-nth_used,nth_used
       call eqitem(theta(i), dbish(:,:,1), grad(i,1), 'R')
       call eqitem(theta(i), dbish(:,:,2), grad(i,2), 'Z')
    enddo

!     to get grad(pressure), multiply grad(psi) by dpressure/dpsi

    if(char == 'R') then
       rpt(1) = rp
       call inter_d_cspl(nr, eqpsi, pressure, 1, rpt, aa, daa)
       do i=-nth_used, nth_used
          grad(i,1)=grad(i,1)*daa(1) * 0.5*beta_0
          grad(i,2)=grad(i,2)*daa(1) * 0.5*beta_0
       enddo
    endif
  end subroutine bgradient

  subroutine eqitem(theta_in, f, fstar, char)
    use constants, only: pi=>pi
    implicit none
    real, intent(in) :: theta_in
    real, dimension(:,-nt:), intent(in) :: f
    real, intent(out) :: fstar
    character(1), intent(in) :: char
    integer :: j, istar, jstar
    real :: tp
    real :: st, dt
    real, dimension(-nt:nt) :: mtheta
! find r on psi mesh
    
    istar = 2

! Now do theta direction

    tp = mod2pi(theta_in)

! get thet on theta mesh

    mtheta = (/ ( real(j)*pi/real(nt-1), j=-nt,nt) /)
  
! note that theta(1)=0 for local_eq theta 

    jstar=-nt-1
    do j=-nt,nt
       if(tp < mtheta(j)) then
          dt = tp - mtheta(j-1)
          st = mtheta(j) - tp
          jstar=j-1
          exit
       endif
       if(jstar /= -nt-1) write(*,*) 'exit error j'
    enddo
      
! treat theta = pi separately
  
    if(jstar == -nt-1) then
       jstar=nt-1
       dt=mtheta(jstar+1)-mtheta(jstar)
       st=0.
    endif

! use opposite area stencil to interpolate

    fstar=f(istar    , jstar    )  * st &
         +f(istar    , jstar + 1)  * dt
    fstar=fstar &
         /(mtheta(jstar+1)-mtheta(jstar))
!     write(*,*) i, dr, dt, sr, st
!     write(*,*) f(istar,jstar+1),f(istar+1,jstar+1)
!     write(*,*) f(istar,jstar),f(istar+1,jstar)
!     write(*,*) eqpsi(istar),eqpsi(istar+1)
!     write(*,*) mtheta(jstar),mtheta(jstar+1)
  end subroutine eqitem

  subroutine eqdcart(dfm, dfcart)
    implicit none
    real, dimension (:,-nt:,:), intent(in)  :: dfm
    real, dimension (:,-nt:,:), intent(out) :: dfcart
    real, dimension (nr,-nt:nt) :: denom

    integer :: i, j
      
    denom(:,:) = drm(:,:,1)*dzm(:,:,2) - drm(:,:,2)*dzm(:,:,1)

    dfcart = 0.
    
    dfcart(:,:,1) =   dfm(:,:,1)*dzm(:,:,2) - dzm(:,:,1)*dfm(:,:,2)
    dfcart(:,:,2) = - dfm(:,:,1)*drm(:,:,2) + drm(:,:,1)*dfm(:,:,2)
    
    do j=-nt,nt
       do i=2,nr
          dfcart(i,j,:)=dfcart(i,j,:)/denom(i,j)
       enddo
    enddo
  end subroutine eqdcart

  subroutine eqdbish(dcart, dbish)
    implicit none
    real, dimension(:,-nt:,:), intent (in) :: dcart
    real, dimension(:,-nt:,:), intent(out) :: dbish
    real, dimension(nr,-nt:nt) :: denom

    integer :: i, j

    denom(:,:) = sqrt(dpcart(:,:,1)**2 + dpcart(:,:,2)**2)

    dbish(:,:,1) = dcart(:,:,1)*dpcart(:,:,1) + dcart(:,:,2)*dpcart(:,:,2)
    dbish(:,:,2) =-dcart(:,:,1)*dpcart(:,:,2) + dcart(:,:,2)*dpcart(:,:,1)
    
    do j=-nt,nt
       do i=2,nr
          dbish(i,j,:) = dbish(i,j,:)/denom(i,j)
       enddo
    enddo
  end subroutine eqdbish

  function invR (r, theta)
    real, intent (in) :: r, theta
    real :: invR
    
    invR=1./Rpos(r, theta)
  end function invR

  function rcenter()
    real :: rcenter
    
    ! If using the global geometry, Rmaj is the center of the shaped flux surface
    ! (i.e. the surface at rhoc=aSurf), so here we calculate the major radial
    ! location of the center of the flux surface of interest ! JRB
    if (surf%geoType==1) then
       rcenter = surf%Rmaj+(1-(surf%r/surf%aSurf)**2)*surf%sHorz
    else
       rcenter = surf%Rmaj
    end if
  end function rcenter

  function Rpos (r, theta)
    implicit none
    real, intent (in) :: r, theta
    real :: Rpos, Zpos

    ! Choose one of four different analytic geometry specifications
    ! (see section 2.1 of "GS2 analytic geometry specification") ! JRB
    select case (surf%geoType)
       case (0)
          call geo_Miller(r, theta, Rpos, Zpos)
       case (1)
          call geo_global(r, theta, Rpos, Zpos)
       case (2)
          call geo_generalizedEllipticity(r, theta, Rpos, Zpos)
       case (3)
          call geo_FourierSeries(r, theta, Rpos, Zpos)
       case default
          write (*,*) "ERROR: invalid analytic geometry specification"
    end select
  end function Rpos

  function Zpos (r, theta)
    implicit none
    real, intent (in) :: r, theta
    real :: Rpos, Zpos

    ! Choose one of four different analytic geometry specifications
    ! (see section 2.1 of "GS2 analytic geometry specification") ! JRB
    select case (surf%geoType)
       case (0)
          call geo_Miller(r, theta, Rpos, Zpos)
       case (1)
          call geo_global(r, theta, Rpos, Zpos)
       case (2)
          call geo_generalizedEllipticity(r, theta, Rpos, Zpos)
       case (3)
          call geo_FourierSeries(r, theta, Rpos, Zpos)
       case default
          write (*,*) "ERROR: invalid analytic geometry specification"
    end select
  end function Zpos

  subroutine geo_Miller(r, theta, Rpos, Zpos)
    ! see section 3.2.1 of Ball MIT Masters thesis or section 2.1.1 of
    ! "GS2 analytic geometry specification"
    real, intent (in) :: r, theta
    real, intent (out) :: Rpos, Zpos
    real :: dr, thAdj 
    real :: Rcirc, Rcircp, Relong, Relongp, RelongTilt, RelongTiltp, Rtri, Rtrip, RtriTilt, RtriTiltp, Rfinal, Rfinalp
    real :: Zcirc, Zcircp, Zelong, Zelongp, ZelongTilt, ZelongTiltp, Ztri, Ztrip, ZtriTilt, ZtriTiltp, Zfinal, Zfinalp

    dr = r - surf%r

    thAdj = theta + thetaShift

    Rcirc=surf%r*cos(surf%thn-surf%thm-thAdj)
    Zcirc=-surf%r*sin(surf%thn-surf%thm-thAdj)
    Rcircp=cos(surf%thn-surf%thm-thAdj)
    Zcircp=-sin(surf%thn-surf%thm-thAdj)

    Relong=Rcirc
    Zelong=-surf%r*(surf%delm-1)*sin(surf%thn-surf%thm-thAdj)+Zcirc
    Relongp=Rcircp
    Zelongp=-((surf%delm-1)+surf%r*surf%delmp)*sin(surf%thn-surf%thm-thAdj)+Zcircp

    RelongTilt=cos(surf%thn-surf%thm)*Relong-sin(surf%thn-surf%thm)*Zelong
    ZelongTilt=sin(surf%thn-surf%thm)*Relong+cos(surf%thn-surf%thm)*Zelong
    RelongTiltp=cos(surf%thn-surf%thm)*Relongp-sin(surf%thn-surf%thm)*Zelongp
    ZelongTiltp=sin(surf%thn-surf%thm)*Relongp+cos(surf%thn-surf%thm)*Zelongp

    Rtri=-surf%r*cos(thAdj)+surf%r*cos(sin(thAdj)*surf%deln+thAdj)+RelongTilt
    Ztri=ZelongTilt
    Rtrip=-cos(thAdj)+cos(sin(thAdj)*surf%deln+thAdj)-surf%r*sin(sin(thAdj)*surf%deln+thAdj)*sin(thAdj)*surf%delnp+RelongTiltp
    Ztrip=ZelongTiltp

    RtriTilt=cos(surf%thn)*Rtri+sin(surf%thn)*Ztri
    ZtriTilt=-sin(surf%thn)*Rtri+cos(surf%thn)*Ztri
    RtriTiltp=cos(surf%thn)*Rtrip+sin(surf%thn)*Ztrip
    ZtriTiltp=-sin(surf%thn)*Rtrip+cos(surf%thn)*Ztrip

    Rfinal=RtriTilt+surf%Rmaj
    Zfinal=ZtriTilt
    Rfinalp=RtriTiltp+surf%sHorz
    Zfinalp=ZtriTiltp+surf%sVert

    Rpos=Rfinal+dr*Rfinalp
    Zpos=Zfinal+dr*Zfinalp

  end subroutine geo_Miller

  subroutine geo_global(r, theta, Rpos, Zpos)
    ! see section 2.1.2 of
    ! "GS2 analytic geometry specification"
    real, intent (in) :: r, theta
    real, intent (out) :: Rpos, Zpos
    real :: rho0, thAdj, rCyl, Rc, Zc
    real :: Cm, Cn, psiN, xAng, yAng, rAng   

    rho0 = r/surf%aSurf

    thAdj = theta + thetaShift

    Cm=(surf%delm**2-1)/(surf%delm**2+1)
    Cn=(surf%deln**2-1)/(surf%deln**3+1)

    psiN=rho0**2+Cm*rho0**2+Cn*rho0**3

    if (abs(Cn*cos(3*(thAdj+surf%thn)))<10e-6) then
      rCyl=surf%aSurf*sqrt(psiN/(1+Cm*cos(2*(thAdj+surf%thm))))
    else
      xAng=2*(1+Cm*cos(2*(thAdj+surf%thm)))**3-27*Cn**2*psiN*cos(3*(thAdj+surf%thn))**2
      yAng=3*Cn*sqrt(3*psiN)*cos(3*(thAdj+surf%thn))*sqrt(xAng+2*(1+Cm*cos(2*(thAdj+surf%thm)))**3)

      rAng=(1/3.0)*atan2(yAng,xAng)

      rCyl=surf%aSurf*(1+Cm*cos(2*(thAdj+surf%thm)))/(3*Cn*cos(3*(thAdj+surf%thn))) &
           *(cos(rAng)+sqrt(3.0)*sin(rAng)-1)
    end if

    Rc=surf%sHorz-(surf%sHorz-surf%Rmaj)*rho0**2
    Zc=surf%sVert-(surf%sVert)*rho0**2

    Rpos=Rc+rCyl*cos(thAdj)
    Zpos=Zc+rCyl*sin(thAdj)

  end subroutine geo_global

  subroutine geo_generalizedEllipticity(r, theta, Rpos, Zpos)
    ! see section 5.1.2 of Ball Oxford PhD thesis or section 2.1.3 of
    ! "GS2 analytic geometry specification"
    real, intent (in) :: r, theta
    real, intent (out) :: Rpos, Zpos
    real :: dr, thAdj, rCyl, rCylp, Rfinal, Zfinal, Rfinalp, Zfinalp

    dr = r - surf%r

    thAdj = theta + thetaShift

    rCyl=surf%r * (1 &
        + (-1+surf%delm/sqrt(1+(surf%delm**2-1)*cos(surf%mMode*(thAdj+surf%thm)/2.0)**2)) &
        + (-1+surf%deln/sqrt(1+(surf%deln**2-1)*cos(surf%nMode*(thAdj+surf%thn)/2.0)**2)))
    rCylp=1 &
     +(-1+(1+(surf%delm**2-1)*cos(surf%mMode*(thAdj+surf%thm)/2.0)**2)**(-0.5)*(surf%delm+surf%r*surf%delmp &
      *(1-(surf%delm*cos(surf%mMode*(thAdj+surf%thm)/2.0))**2/(1+(surf%delm**2-1)*cos(surf%mMode*(thAdj+surf%thm)/2.0)**2)))) &
     +(-1+(1+(surf%deln**2-1)*cos(surf%nMode*(thAdj+surf%thn)/2)**2)**(-0.5)*(surf%deln+surf%r*surf%delnp &
      *(1-(surf%deln*cos(surf%nMode*(thAdj+surf%thn)/2.0))**2/(1+(surf%deln**2-1)*cos(surf%nMode*(thAdj+surf%thn)/2.0)**2))))

    Rfinal=rCyl*cos(thAdj)+surf%Rmaj
    Zfinal=rCyl*sin(thAdj)
    Rfinalp=rCylp*cos(thAdj)+surf%sHorz
    Zfinalp=rCylp*sin(thAdj)+surf%sVert

    Rpos=Rfinal+dr*Rfinalp
    Zpos=Zfinal+dr*Zfinalp

  end subroutine geo_generalizedEllipticity

  subroutine geo_FourierSeries(r, theta, Rpos, Zpos)
    ! see section 5.1.1 of Ball Oxford PhD thesis or section 2.1.4 of
    ! "GS2 analytic geometry specification"
    real, intent (in) :: r, theta
    real, intent (out) :: Rpos, Zpos
    real :: dr, thAdj, rCyl, rCylp, Rfinal, Zfinal, Rfinalp, Zfinalp

    dr = r - surf%r

    thAdj = theta + thetaShift

    rCyl=surf%r * (1+((surf%delm-1)/2) * (1-cos(surf%mMode*(thAdj+surf%thm))) + ((surf%deln-1)/2) * (1-cos(surf%nMode*(thAdj+surf%thn))))
    rCylp=1 + (((surf%delm-1)/2) + (surf%r*surf%delmp/2)) * (1-cos(surf%mMode*(thAdj+surf%thm))) + (((surf%deln-1)/2) + (surf%r*surf%delnp/2)) * (1-cos(surf%nMode*(thAdj+surf%thn)))

    Rfinal=rCyl*cos(thAdj)+surf%Rmaj
    Zfinal=rCyl*sin(thAdj)
    Rfinalp=rCylp*cos(thAdj)+surf%sHorz
    Zfinalp=rCylp*sin(thAdj)+surf%sVert

    Rpos=Rfinal+dr*Rfinalp
    Zpos=Zfinal+dr*Zfinalp

  end subroutine geo_FourierSeries

  function psi (r)
    real, intent (in) :: r
    real :: psi

    psi = r - surf%r
  end function psi

  function mod2pi (theta)
    use constants, only: pi=>pi
    real, intent(in) :: theta
    real :: th, mod2pi
    real, parameter :: theta_tol = 1.e-6
    logical :: out
    if(theta <= pi .and. theta >= -pi) then
       mod2pi = theta
       return
    endif
    
    if(theta - theta_tol <= pi .and. theta >= -pi) then
       mod2pi = pi
       return
    endif

    if(theta <= pi .and. theta + theta_tol >= -pi) then
       mod2pi = -pi
       return
    endif

    th=theta
    out=.true.
    do while(out)
       if(th > pi) th = th - 2.*pi
       if(th <-pi) th = th + 2.*pi
       if(th <= pi .and. th >= -pi) out=.false.
    enddo
    mod2pi=th
  end function mod2pi
   
  function diameter (rp)
    real, intent(in) :: rp
    real :: diameter

    diameter = 2.*rp
  end function diameter

  function dbtori ()
    real :: dbtori
    dbtori = 1.
  end function dbtori

  function btori ()
    real :: btori
    btori = surf%Rgeo
  end function btori

  function qfun ()
    real :: qfun
    qfun = surf%q
  end function qfun

  function pfun ()
    real :: pfun
    pfun = 0.5*beta_0
  end function pfun
  
  function dpfun ()  
    real :: dpfun    
    dpfun = -1.
  end function dpfun

  function dpdrhofun()
    real :: dpdrhofun
    dpdrhofun = 0 ! surf%pp
  end function dpdrhofun
  
  function betafun ()  
    real :: betafun
    betafun = beta_0
  end function betafun
end module leq
