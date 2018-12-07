





module command_line
  ! (c) Copyright 1991 to 1998 by Michael A. Beer, William D. Dorland, 
  ! P. B. Snyder, Q. P. Liu, and Gregory W. Hammett. ALL RIGHTS RESERVED.
  !
  ! <doc>
  !  A wrapper module for handling command line arguments.
  !  This module provides subroutine cl_getarg and integer function cl_iargc.
  !  Most of the compilers have getarg and iargc as their extensions.
  !  If not, one can use POSIX pxfgetarg and ipxfargc.
  !
  !  Note that Fortran 2003 includes get_command_argument and 
  !  command_argument_count which will replace getarg and iargc.
  !
  !  RN> 20180726
  !  Added preprocessor macro F200X_INTRINSICS to invoke Fortran 2003/2008 intrinsics
  ! </doc>

  implicit none

  private

  public :: cl_getarg, cl_iargc
  public :: cl_getenv

contains

  function cl_iargc()

    ! <doc>
    !  returns the number of arguments
    !  using intrinsic iargc or POSIX ipxfargc
    ! </doc>


    implicit none
    integer :: cl_iargc
    integer :: iargc

    cl_iargc = iargc()
  end function cl_iargc

  subroutine cl_getarg (k, arg, len, ierr)
    
    ! <doc>
    !  gets k-th argument string and its length
    !  using intrinsic getarg or POSIX pxfgetarg
    ! </doc>


    implicit none
    integer,           intent (in)  :: k
    character (len=*), intent (out) :: arg
    integer,           intent (out) :: len
    integer,           intent (out) :: ierr

    call getarg (k, arg)
    len = len_trim(arg)
    ierr = 0

  end subroutine cl_getarg

  subroutine cl_getenv(name,val,len,ierr)
    implicit none
    character (len=*), intent(in) :: name
    character (len=*), intent(out) :: val
    integer, intent(out), optional :: len
    integer, intent(out), optional :: ierr
    
    call getenv (name,val)
    if (present(len)) len = len_trim(val)
    if (present(ierr)) ierr = 0
    
  end subroutine cl_getenv
  
end module command_line
