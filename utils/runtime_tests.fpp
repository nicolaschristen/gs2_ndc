#include "define.inc"
!>  This module is intended to be used for runtime tests
!!  which interrogate what is functional/what compile time
!!  options were enabled/disabled, and also any enviroment variables
!!  which affect what happens at runtime.
!!
module runtime_tests

  implicit none

  private

  public :: compiler_pgi, verbosity, build_identifier
  public :: get_svn_rev, get_compiler_name
  public :: is_release, release
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Tests for compilers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function compiler_pgi()
    logical :: compiler_pgi
    compiler_pgi = .false.
#if FCOMPILER == _PGI_
    compiler_pgi = .true.
#endif
  end function compiler_pgi

  function get_compiler_name()
    character(len=9) :: get_compiler_name
    get_compiler_name='unknown'
#if FCOMPILER == _PGI_
    get_compiler_name='pgi'
#elif FCOMPILER == _INTEL_
    get_compiler_name='intel'
#elif FCOMPILER == _GFORTRAN_
    get_compiler_name='gfortran'
#elif FCOMPILER == _XL_
    get_compiler_name='xl'
#elif FCOMPILER == _NAG_
    get_compiler_name='nag'
#elif FCOMPILER == _CRAY_
    get_compiler_name='cray'
#elif FCOMPILER == _G95_
    get_compiler_name='g95'
#elif FCOMPILER == _PATHSCALE_
    get_compiler_name='pathscale'
#elif FCOMPILER == _LAHEY_
    get_compiler_name='lahey'
#elif FCOMPILER == _ABSOFT_
    get_compiler_name='absoft'
#elif FCOMPILER == _ALPHA_
    get_compiler_name='alpha'
#elif FCOMPILER == _SUN_
    get_compiler_name='sun'
#endif
  end function get_compiler_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Tests for svn info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>This function returns the output of svnversion 
  !It would be nice if we could strip of trailing text
  !so that we're just left with the integer revision.
  function get_svn_rev()
    character(len=10) :: get_svn_rev
#ifndef SVN_REV
#define SVN_REV "unknown"
#endif
    get_svn_rev=SVN_REV
  end function get_svn_rev

  !>This function returns true if the source code has
  !been modified relative to repo
  function get_svn_modified()
    logical :: get_svn_modified
    integer :: indx
#ifndef SVN_REV
#define SVN_REV "unknown"
#endif
    indx=index(SVN_REV,"M")
    if(indx.eq.0)then
       get_svn_modified=.false.
    else
       get_svn_modified=.true.
    endif
  end function get_svn_modified
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! System info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function get_gk_system()
    character(len=20) :: get_gk_system
#ifndef GK_SYSTEM
#define GK_SYSTEM "unknown"
#endif
    get_gk_system=GK_SYSTEM
  end function get_gk_system

  !> This function returns an identifier of the system and build:
  !! "system.compiler.svnrevision". 
  function build_identifier()
    character(len=50) :: build_identifier
    character(len=10) :: svnrev
    !Do we want a call to get_svn_modified here so that we can highlight
    !if this revision is modified or not? 
    svnrev = get_svn_rev()
    build_identifier = trim(get_gk_system())//"."//trim(get_compiler_name())//'.'//svnrev
  end function build_identifier

  function is_release()
    logical :: is_release
#ifdef RELEASE
    is_release = .true.
#else
    is_release = .false.
#endif
  end function is_release

  function release()
    character(len=10) :: release
#ifdef RELEASE
    release = RELEASE
#else
    release = ''
#endif
  end function release

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Testing the runtime environment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !> This function interrogates the environment variable
  !! GK_VERBOSITY and returns its integer value. This is used
  !! to control the level of debug output (not diagnostic/physics output).
  !! Normal levels range from 0 to 5, with output getting 
  !! heavier as the value increases. Values higher than 5 can be used for 
  !! specialised/very heavy output.
  function verbosity()
    integer :: verbosity
    character(len=10) :: verbosity_char
    verbosity_char = ''
    !call getenv("VERBOSITY", verbosity_char)
    !For fortran 2003 standard should replace above with
    call get_environment_variable("GK_VERBOSITY", verbosity_char)
    read (verbosity_char,'(I10)') verbosity
    !write (*,*) 'verbosity is ', verbosity
  end function verbosity
end module runtime_tests
