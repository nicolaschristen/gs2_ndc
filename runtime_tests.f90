





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
  end function compiler_pgi

  function get_compiler_name()
    character(len=9) :: get_compiler_name
    get_compiler_name='unknown'
    get_compiler_name='gfortran'
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
    get_svn_rev="Unversioned directory"
  end function get_svn_rev

  !>This function returns true if the source code has
  !been modified relative to repo
  function get_svn_modified()
    logical :: get_svn_modified
    integer :: indx
    indx=index("Unversioned directory","M")
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
    get_gk_system="gnu_ubuntu"
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
    is_release = .false.
  end function is_release

  function release()
    character(len=10) :: release
    release = ''
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
