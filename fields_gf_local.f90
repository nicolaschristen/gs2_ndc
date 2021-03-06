!AJ This module is based on the fields_local module (which is in term a re-implementation of the 
!AJ fields implicit functionality).  fields_gf_local implements the implicit field calculation 
!AJ using a gf_lo based data decomposition to improve performance on large process counts.  It
!AJ utilises a smaller set of processes to reduce communication costs (at the expense of significant 
!AJ load imbalance in the fields calculations).  It assumes gf_lo velocity space integration and 
!AJ stores field data in gf_lo format (i.e. ik,it are decomposed over processes, everthing else is 
!AJ local).
!AJ Implemented by Adrian Jackson, EPCC, The University of Edinburgh; 2014-2016.  Funded by an 
!AJ EPSRC ARCHER eCSE project.
module fields_gf_local
  use mp, only: comm_type
  implicit none
  private

  !> Module level routines
  public :: init_fields_gf_local, init_allfields_gf_local, finish_fields_gf_local
  public :: advance_gf_local, reset_fields_gf_local, dump_response_to_file_gf_local
  public :: fields_gf_local_functional, read_response_from_file_gf_local
  public :: field_local_allreduce_sub

  !> Unit tests
  public :: fields_gf_local_unit_test_init_fields_matrixlocal

  !> User set parameters
  public :: dump_response, read_response

  !> Made public for testing
  public :: fieldmat

  !//////////////////////////////
  !// CUSTOM TYPES
  !//////////////////////////////

  !>This is the lowest level of data object and will hold the actual data
  !involved with the field matrix
  type, private :: rowblock_type
     complex, dimension(:,:), allocatable :: data !This is the field data
     integer :: row_llim, row_ulim !The row limits this block corresponds to
     integer :: col_llim, col_ulim !The col limits this block corresponds to
     integer :: nrow, ncol !The size of the row block
     complex, dimension(:),allocatable :: tmp_sum !Work space for field update
     !/Row and parent properties. Mostly for debug printing
     integer :: ir_ind !Property of rowblock
     integer :: it_ind, ik_ind !Details of parents
     integer :: is_ind, ic_ind !Details of parents so cell can find itself
   contains
     private
     procedure :: deallocate => rb_deallocate
     procedure :: allocate => rb_allocate
     procedure :: debug_print => rb_debug_print
     procedure :: mv_mult_fun => rb_mv_mult_fun
     procedure :: mv_mult => rb_mv_mult
     procedure :: irex_to_ir => rb_irex_to_irloc
     procedure :: reset => rb_reset
     procedure :: set_nrow => rb_set_nrow
  end type rowblock_type

  !>This is the next level up of data and represents the cell type.
  !A cell has a unique it,ik pair and represents all the fields/field eq
  !on the extended domain for the given it/ik
  type, private :: cell_type
     type(rowblock_type), dimension(:), allocatable :: rb !These are the row blocks, currently one for each field equation
     integer :: nrb !How many row blocks, equal to nfield
     !//NOTE: probably don't want these as this data should be stored by the row blocks only
     integer :: row_llim, row_ulim !The row limits, always 1 to nfield*nextend
     integer, dimension(:), allocatable :: col_llim, col_ulim !These are the column limits for each fdq
     !//////
     integer :: nrow, ncol !The number of rows and columns for each field equation.
     integer :: ncol_tot !The total number of columns (nfdq*ncol)
     logical :: is_local !Does this cell have any data on this proc?
     logical :: is_empty !Have we got any data for this cell on this proc?
     logical :: is_all_local !Is all of this cells data on this proc?
     logical :: ignore_boundary !Do we ignore the boundary point (column) in this cell?
     complex, dimension(:),allocatable :: tmp_sum !Work space for field update
     type(comm_type) :: parent_sub !Sub communicator involving all processors in parent
     !Cell and parent properties. Mostly for debug printing.
     integer :: ic_ind, it_ind !Cell properties
     integer :: ik_ind, is_ind !Parent properties
   contains
     private
     procedure :: deallocate => c_deallocate
     procedure :: allocate => c_allocate
     procedure :: debug_print => c_debug_print
     procedure :: mv_mult_rb => c_mv_mult_rb
     procedure :: get_field_update => c_get_field_update
     procedure :: reset => c_reset
     procedure :: set_locality => c_set_locality
     procedure :: has_row => c_has_row
  end type cell_type

  !>This is the next level up of data and represents the supercell type.
  !A supercell represents a collection of connected cells
  type, private :: supercell_type
     type(cell_type), dimension(:), allocatable :: cells !These are the cells
     integer :: ncell
     integer :: nextend !Length of the extended domain
     integer :: nrow, ncol !The number of rows and columns. Equal to nextend*nfield
     integer :: it_left !It index of leftmost cell
     integer :: head_iproc_world !The proc id (in comm world) of the head proc
     integer :: head_iproc !The proc id (in sc_sub_all) of the head proc
     integer :: head_iproc_pd !The proc id (in sc_sub_pd) of the head proc
     logical :: is_local !Does this supercell have any data on this proc?
     logical :: is_empty !Have we got any data for this supercell on this proc?
     logical :: is_all_local !Is all of this supercells data on this proc?
     complex, dimension(:),allocatable :: tmp_sum !Work space for field update
     type(comm_type) :: sc_sub_all !Sub communicator involving all processors with this supercell
     type(comm_type) :: sc_sub_pd !Sub communicator for all procs with some data but not all of it
     type(comm_type) :: sc_sub_not_full !Sub communicator involving all processors with this supercell
     type(comm_type) :: parent_sub !Sub communicator involving all processors in parent
     integer, dimension(:),allocatable :: nb_req_hand !For non-blocking broadcast request handle storage
     logical :: initdone !Have we finished initialising this block?
     logical, dimension(:), allocatable :: initialised !Have we initialised each point?
     logical :: is_head=.false. !Are we the head of this supercell?
     !Cell and parent properties. Mostly for debug printing.
     integer :: is_ind
     integer :: ik_ind !Parent properties
   contains
     private
     procedure :: deallocate => sc_deallocate
     procedure :: allocate => sc_allocate
     procedure :: debug_print => sc_debug_print
     procedure :: get_field_update => sc_get_field_update
     procedure :: reduce_tmpsum => sc_reduce_tmpsum
     procedure :: iex_to_dims => sc_iex_to_dims
     procedure :: iex_to_ifl => sc_iex_to_ifl
     procedure :: iex_to_ic => sc_iex_to_ic
     procedure :: iex_to_ig => sc_iex_to_ig
     procedure :: iex_to_it => sc_iex_to_it
     procedure :: has_it => sc_has_it
     procedure :: reset => sc_reset
     procedure :: set_locality => sc_set_locality
     procedure :: get_left_it => sc_get_left_it
     procedure :: store_fq => sc_store_fq
     procedure :: pull_rows_to_arr => sc_pull_rows_to_arr
     procedure :: push_arr_to_rows => sc_push_arr_to_rows
     procedure :: prepare => sc_prepare
     procedure :: invert => sc_invert
     procedure :: invert_local => sc_invert_local
     procedure :: invert_mpi => sc_invert_mpi
     procedure :: dump => sc_dump
     procedure :: make_subcom_1 => sc_make_subcom_1
     procedure :: make_subcom_2 => sc_make_subcom_2
  end type supercell_type

  !>This is the next level up of data and represents the ky type.
  !A ky block consists of a number of supercells
  type, private :: ky_type
     type(supercell_type), dimension(:), allocatable :: supercells !These are the supercells
     integer :: nsupercell
     logical :: is_local !Does this supercell have any data on this proc?
     logical :: is_empty !Have we got any data for this supercell on this proc?
     logical :: is_all_local !Is all of this supercells data on this proc?
     logical :: initdone !Have we finished initialising this block?
     type(comm_type) :: ky_sub_all !Sub communicator involving all processors with this ky
     type(comm_type) :: parent_sub !Sub communicator involving all processors in parent
     !Cell and parent properties. Mostly for debug printing.
     integer :: ik_ind
   contains
     private
     procedure :: deallocate => ky_deallocate
     procedure :: allocate => ky_allocate
     procedure :: debug_print => ky_debug_print
     procedure :: get_field_update => ky_get_field_update
     procedure :: reset => ky_reset
     procedure :: set_locality => ky_set_locality
     procedure :: is_from_it => ky_is_from_it
     procedure :: store_fq => ky_store_fq
     procedure :: prepare => ky_prepare
     procedure :: make_subcom_1 => ky_make_subcom_1
     procedure :: make_subcom_2 => ky_make_subcom_2
  end type ky_type

  !>This is the top level object, consisting of a collection of
  !ky blocks.
  type, private :: fieldmat_type
     type(ky_type), dimension(:), allocatable :: kyb !The ky blocks
     integer :: naky !Number of ky blocks
     integer :: ntheta0 
     integer :: npts !Total number of theta grid points on all extended domains
     integer :: nbound !Number of ignored boundary points
     logical :: is_local !Does this supercell have any data on this proc?
     logical :: is_empty !Have we got any data for the fieldmat on this proc?
     logical :: is_all_local !Is all of this supercells data on this proc?
     type(comm_type) :: fm_sub_all !Sub communicator involving all processors with fieldmat
     type(comm_type) :: fm_sub_headsc_p0 !Sub communicator involving the supercell heads and proc0
     integer :: prepare_type=0 !What sort of preparation do we do to the field matrix
     integer, dimension(:,:), allocatable :: heads ! An array that holds the iproc (from comm world) of the supercell head for a given ik,it point
     !Currently only support:
     !0   : Invert the field matrix
     !In the future may wish to do things like LU decomposition etc.
     !Could also do things like eigenvalue analysis etc.
   contains
     private
     procedure :: deallocate => fm_deallocate
     procedure :: allocate => fm_allocate
     procedure :: debug_print => fm_debug_print
     procedure :: get_field_update => fm_get_field_update
     procedure :: init => fm_init
     procedure :: reset => fm_reset
     procedure :: populate => fm_populate
     procedure :: init_next_field_points => fm_init_next_field_points
     procedure :: prepare => fm_prepare
     procedure :: make_subcom_1 => fm_make_subcom_1
     procedure :: make_subcom_2 => fm_make_subcom_2
     procedure :: calc_sc_heads => fm_calc_sc_heads
     procedure :: gather_fields => fm_gather_fields
     procedure :: gather_init_fields => fm_gather_init_fields
     procedure :: scatter_fields => fm_scatter_fields
     procedure :: unpack_to_field => fm_unpack_to_field
     procedure :: write_debug_data => fm_write_debug_data
     procedure :: set_is_local => fm_set_is_local
     procedure :: count_subcom => fm_count_subcom
     procedure :: getfieldeq_nogath => fm_getfieldeq_nogath
     procedure :: getfieldeq1_nogath => fm_getfieldeq1_nogath
     procedure :: update_fields => fm_update_fields
     procedure :: update_fields_newstep => fm_update_fields_newstep
  end type fieldmat_type

  !>This is the type which controls/organises the parallel
  !data decomposition etc.
  type, private :: pc_type
     integer, dimension(:,:), allocatable :: is_local !Is the given it,ik on this proc?
     !NOTE: We have is_local as integer rather than logical so we can reduce across procs
     !to count how many procs have a given cell
     type(comm_type), dimension(:,:), allocatable :: itik_subcom !Cell level sub-communicators
     integer, dimension(:,:), allocatable :: nproc_per_cell !How many procs have this cell
     integer, dimension(:,:), allocatable :: nresp_per_cell !How many rows are we personally responsible for in each cell
     integer, dimension(:,:), allocatable :: navail_per_cell !How many rows could we be responsible for in each cell
     integer :: nresp_tot, navail_tot !Total number of rows available/potentially available
     logical :: force_local_invert=.false. !If true then we force local inversion
     !the local fields. This may impact diagnostics/parallel writing as it means only the processors gf_lo local section
     !of the fields are evolved on this processor.
     !NOTE: If we're doing force local invert really when we populate the supercell we should store all of the 
     !data locally so we don't need to fetch it later.
   contains
     private
     procedure :: deallocate => pc_deallocate
     procedure :: allocate => pc_allocate
     procedure :: debug_print => pc_debug_print
     procedure :: current_nresp => pc_current_nresp
     procedure :: find_locality => pc_find_locality
     procedure :: count_avail => pc_count_avail
     procedure :: has_ik => pc_has_ik
     procedure :: has_it => pc_has_it
     procedure :: make_decomp => pc_make_decomp
     procedure :: init => pc_init
     procedure :: reset => pc_reset
     !//The different decomposition routines
     procedure :: decomp => pc_decomp
  end type pc_type

  !//////////////////////////////
  !// MODULE VARIABLES
  !//////////////////////////////
  integer :: lun=6 !File unit for printing
  integer :: dlun=6 !File unit for debug printing
  !Tuning this can improve advance/init time (usually one at the expense of other)
  logical,parameter :: debug=.false. !Do we do debug stuff?
  logical :: initialised=.false. !Have we initialised yet?
  logical :: reinit=.false. !Are we reinitialising?
  logical :: dump_response=.false. !Do we dump the response matrix?
  logical :: read_response=.false. !Do we read the response matrix from dump?
  logical :: field_local_allreduce_sub = .false. !If true and field_local_allreduce true then do two sub comm all reduces rather than 1
  integer :: nfield !How many fields
  integer :: nfq !How many field equations (always equal to nfield)
  type(pc_type),save :: pc !This is the parallel control object
  type(fieldmat_type),save :: fieldmat !This is the top level field matrix object

contains

!//////////////////////////////
!// CUSTOM TYPES
!//////////////////////////////
!
!================
!ROWBLOCK
!================
  !>Allocate storage space
  subroutine rb_allocate(self)
    implicit none
    class(rowblock_type), intent(inout) :: self
    allocate(self%data(self%ncol,self%nrow))
    allocate(self%tmp_sum(self%nrow))
  end subroutine rb_allocate

  !>Deallocate storage space
  subroutine rb_deallocate(self)
    implicit none
    class(rowblock_type), intent(inout) :: self
    if(allocated(self%data)) deallocate(self%data)
    if(allocated(self%tmp_sum)) deallocate(self%tmp_sum)
  end subroutine rb_deallocate

  !>Debug printing
  subroutine rb_debug_print(self)
    implicit none
    class(rowblock_type), intent(in) :: self
    write(dlun,'("Rowblock debug print. Index ir=",I0," nrow=",I0," ncol=",I0," rl,ru=",I0," ",I0, " and cl,cu=",I0," ",I0)') &
         self%ir_ind, self%nrow, self%ncol, &
         self%row_llim, self%row_ulim, &
         self%col_llim, self%col_ulim
  end subroutine rb_debug_print

  !>Matrix vector multiplication
  function rb_mv_mult_fun(self,vect)
    implicit none
    class(rowblock_type), intent(in) :: self
    complex,dimension(self%ncol), intent(in) :: vect
    complex,dimension(self%nrow) :: rb_mv_mult_fun
!    integer :: ir
    !As a 1d vector in fortran is a "row vector" (as the
    !first index represents columns) matmul is a bit annoying
    !when trying to do a matrix-vector product
    rb_mv_mult_fun=-matmul(transpose(self%data),vect)
  end function rb_mv_mult_fun

  !>Matrix vector multiplication stored in tmp_sum
  subroutine rb_mv_mult(self,vect)
    implicit none
    class(rowblock_type), intent(inout) :: self
    complex,dimension(self%ncol), intent(in) :: vect
    self%tmp_sum=self%mv_mult_fun(vect)
  end subroutine rb_mv_mult
  
  !>Convert extended row index to local row index
  function rb_irex_to_irloc(self,irex)
    implicit none
    class(rowblock_type), intent(in) :: self
    integer, intent(in) :: irex
    integer :: rb_irex_to_irloc
    rb_irex_to_irloc=1+(irex-self%row_llim)
  end function rb_irex_to_irloc

  !>Set the nrow variables
  subroutine rb_set_nrow(self)
    implicit none
    class(rowblock_type), intent(inout) :: self
    if(self%row_ulim.le.0)then
       self%nrow=0
    else
       self%nrow=1+self%row_ulim-self%row_llim
    endif
  end subroutine rb_set_nrow

  !>A routine to reset the object
  subroutine rb_reset(self)
    implicit none
    class(rowblock_type), intent(inout) :: self

    !deallocate
    call self%deallocate

    !Could zero out variables but no real need
    self%row_llim=0
    self%row_ulim=0
  end subroutine rb_reset

!------------------------------------------------------------------------

!================
!CELL
!================
  !>Allocate storage space
  subroutine c_allocate(self)
    implicit none
    class(cell_type), intent(inout) :: self
    integer :: ir
    allocate(self%tmp_sum(self%nrow))
    do ir=1,self%nrb
       call self%rb(ir)%allocate
    enddo
  end subroutine c_allocate

  !>Deallocate storage space
  subroutine c_deallocate(self)
    implicit none
    class(cell_type), intent(inout) :: self
    if(allocated(self%tmp_sum)) deallocate(self%tmp_sum)
    if(allocated(self%rb)) deallocate(self%rb)
  end subroutine c_deallocate

  !>Debug printing
  subroutine c_debug_print(self)
    implicit none
    class(cell_type), intent(in) :: self
    write(dlun,'("Cell debug print. Index ic=",I0," nrow=",I0," ncol=",I0)') &
         self%ic_ind, self%nrow, self%ncol
  end subroutine c_debug_print

  !>Do matrix vector multiplication at rowblock level
  subroutine c_mv_mult_rb(self,vect,ifq)
    implicit none
    class(cell_type), intent(inout) :: self
    complex, dimension(self%ncol), intent(in) :: vect
    integer, intent(in) :: ifq
    call self%rb(ifq)%mv_mult(vect)
  end subroutine c_mv_mult_rb

  !>Get the field update for this cells data
  !Note still need to reduce across other cells in this
  !supercell.
  subroutine c_get_field_update(self,fq,fqa,fqp)
    use run_parameters, only: fphi, fapar, fbpar
    implicit none
    class(cell_type), intent(inout) :: self
    complex, dimension(self%ncol), intent(in) :: fq,fqa,fqp
    integer :: ifq

    !If we don't have any data then exit
    if(self%is_empty) return

    !First do the multiplication at rowblock level
    ifq=0
    if(fphi.gt.epsilon(0.0)) then
       ifq=ifq+1
       call self%mv_mult_rb(fq,ifq)
    endif
    if(fapar.gt.epsilon(0.0)) then
       ifq=ifq+1
       call self%mv_mult_rb(fqa,ifq)
    endif
    if(fbpar.gt.epsilon(0.0)) then
       ifq=ifq+1
       call self%mv_mult_rb(fqp,ifq)
    endif

    !Now store at cell level
    self%tmp_sum=0
    do ifq=1,nfq
       self%tmp_sum(self%rb(ifq)%row_llim:self%rb(ifq)%row_ulim)=&
            self%tmp_sum(self%rb(ifq)%row_llim:self%rb(ifq)%row_ulim)&
            +self%rb(ifq)%tmp_sum
    enddo

  end subroutine c_get_field_update

  !>A routine to reset the object
  subroutine c_reset(self)
    implicit none
    class(cell_type), intent(inout) :: self
    integer :: ir

    !Call reset on all children
    do ir=1,self%nrb
       call self%rb(ir)%reset
    enddo

    !deallocate
    call self%deallocate

    !Could zero out variables but no real need
  end subroutine c_reset

  !>Set the locality of each object
  subroutine c_set_locality(self)
    use kt_grids, only: aky, akx
    implicit none
    class(cell_type), intent(inout) :: self
    
    !Set our locality. Note here we assume all rb have same size
    if(size(self%rb).gt.0)then
       self%is_empty=self%rb(1)%nrow.eq.0
       self%is_all_local=self%rb(1)%nrow.eq.self%nrow
    else
       self%is_empty=.true.
       self%is_all_local=.false.
    endif

    !/Also ignore the ky=kx=0 mode
    self%is_empty=(self%is_empty.or.(abs(aky(self%ik_ind)).lt.epsilon(0.0).and.abs(akx(self%it_ind)).lt.epsilon(0.0)))
    self%is_local=pc%is_local(self%it_ind,self%ik_ind).eq.1

  end subroutine c_set_locality

  !>Test if a given row belongs to the current cell
  function c_has_row(self,irow)
    implicit none
    class(cell_type), intent(in) :: self
    integer, intent(in) :: irow
    logical :: c_has_row
    !NOTE: Here we assume all row blocks in the cell have the same
    !row limits
    if(size(self%rb).gt.0)then
       c_has_row=(irow.ge.self%rb(1)%row_llim.and.irow.le.self%rb(1)%row_ulim)
    else
       c_has_row=.false.
    endif
  end function c_has_row

!------------------------------------------------------------------------

!================
!SUPERCELL
!================
  !>Allocate storage space
  subroutine sc_allocate(self)
    implicit none
    class(supercell_type), intent(inout) :: self
    integer :: ic
    allocate(self%tmp_sum(self%nrow))
    do ic=1,self%ncell
       call self%cells(ic)%allocate
    enddo
  end subroutine sc_allocate

  !>Deallocate storage space
  subroutine sc_deallocate(self)
    implicit none
    class(supercell_type), intent(inout) :: self
    if(allocated(self%tmp_sum)) deallocate(self%tmp_sum)
    if(allocated(self%cells)) deallocate(self%cells)
  end subroutine sc_deallocate

  !>Debug printing
  subroutine sc_debug_print(self)
    implicit none
    class(supercell_type), intent(in) :: self
    write(dlun,'("Supercell debug print. Index is=",I0," nrow=",I0," ncol=",I0)') &
         self%is_ind, self%nrow, self%ncol
  end subroutine sc_debug_print

  !>Get the field update
  subroutine sc_get_field_update(self,fq,fqa,fqp)
    implicit none
    class(supercell_type), intent(inout) :: self
    complex, dimension(:,:), intent(in) :: fq,fqa,fqp
    integer :: ic, it_ind, ulim

    !Exit if we don't have any of the data locally
    if(.not.self%is_local) return

    if(.not.self%is_empty)then
       !Now loop over cells and trigger field update
       do ic=1,self%ncell
          !Get cell properties
          it_ind=self%cells(ic)%it_ind
          ulim=self%cells(ic)%ncol
          !Do update
          call self%cells(ic)%get_field_update(fq(:ulim,it_ind),fqa(:ulim,it_ind),fqp(:ulim,it_ind))
       enddo
    endif

    !Now we need to reduce across cells
    call self%reduce_tmpsum

  end subroutine sc_get_field_update

  !>Reduce the field update across cells to give the final answer
  subroutine sc_reduce_tmpsum(self)
    use mp, only: sum_allreduce_sub, sum_reduce_sub
    implicit none
    class(supercell_type), intent(inout) :: self
    integer :: ic

    !NOTE: Here we rely on the cells having not reduced/
    !gathered their data yet so there are no repeated rows.
    !First do local sums
    do ic=1,self%ncell
       if(self%cells(ic)%is_empty) cycle
       self%tmp_sum=self%tmp_sum+self%cells(ic)%tmp_sum
    enddo

    if(.not.(self%is_empty.or.self%is_all_local)) then 
       call sum_reduce_sub(self%tmp_sum,0,self%sc_sub_pd)
    end if

  end subroutine sc_reduce_tmpsum
    
!<DD>The iex_to routines all need testing!
  !>Convert the extended domain index to ig, it and ifl
  subroutine sc_iex_to_dims(self,iex,ig,ic,it,ifl)
    use theta_grid, only: ntgrid
    implicit none
    class(supercell_type), intent(in) :: self
    integer, intent(in) :: iex
    integer, intent(out) :: ig,ic,it,ifl

    !Get field index
    ifl=1+int((iex-1)/self%nextend)
    
    !Get cell and it index
    ic=min(self%ncell,1+int((iex-(ifl-1)*self%nextend-1)/(2*ntgrid)))
    it=self%cells(ic)%it_ind

    !Get ig index
    ig=-ntgrid+iex-(ifl-1)*self%nextend-1-2*ntgrid*(ic-1)
  end subroutine sc_iex_to_dims
    
  !>Convert the extended domain index to ifl
  subroutine sc_iex_to_ifl(self,iex,ifl)
    implicit none
    class(supercell_type), intent(in) :: self
    integer, intent(in) :: iex
    integer, intent(out) :: ifl
    ifl=1+(iex-1)/self%nextend
  end subroutine sc_iex_to_ifl

  !>Convert the extended domain index to ic
  subroutine sc_iex_to_ic(self,iex,ic)
    use theta_grid, only: ntgrid
    implicit none
    class(supercell_type), intent(in) :: self
    integer, intent(in) :: iex
    integer, intent(out) :: ic
    integer :: ifl, iex_tmp
    !First work out how many fields there are
    call self%iex_to_ifl(iex,ifl)
    !Ensure we're looking between 1 and self%nextend
    iex_tmp=iex-(ifl-1)*self%nextend
    ic=MIN(int((iex_tmp-1)/(2*ntgrid))+1,self%ncell)
  end subroutine sc_iex_to_ic

  !>Convert the extended domain index to it
  subroutine sc_iex_to_it(self,iex,it)
    use theta_grid, only: ntgrid
    implicit none
    class(supercell_type), intent(in) :: self
    integer, intent(in) :: iex
    integer, intent(out) :: it
    integer :: ifl, iex_tmp, ic
    !First work out how many fields there are
    call self%iex_to_ifl(iex,ifl)
    !Ensure we're looking between 1 and self%nextend
    iex_tmp=iex-(ifl-1)*self%nextend
    ic=MIN(int((iex_tmp-1)/(2*ntgrid))+1,self%ncell)
    it=self%cells(ic)%it_ind
  end subroutine sc_iex_to_it

  !>Convert the extended domain index to ig
  subroutine sc_iex_to_ig(self,iex,ig)
    use theta_grid, only: ntgrid
    implicit none
    class(supercell_type), intent(in) :: self
    integer, intent(in) :: iex
    integer, intent(out) :: ig
    integer :: iex_tmp, ifl
    !First get ifl
    call self%iex_to_ifl(iex,ifl)
    !Limit to 1 and self%nextend
    iex_tmp=iex-(ifl-1)*self%nextend
    
    if(iex_tmp.eq.self%nextend) then
       ig=ntgrid
    else
       ig=-ntgrid+mod(iex_tmp-1,2*ntgrid)
    endif
  end subroutine sc_iex_to_ig

  !>Is the passed it a member of this supercell
  function sc_has_it(self,it)
    implicit none
    class(supercell_type), intent(in) :: self
    integer, intent(in) :: it
    logical :: sc_has_it
    sc_has_it=any(self%cells%it_ind.eq.it)
  end function sc_has_it

  !>A routine to reset the object
  subroutine sc_reset(self)
    use mp, only: free_comm
    implicit none
    class(supercell_type), intent(inout) :: self
    integer :: ic

    !Call reset on all children
    do ic=1,self%ncell
       call self%cells(ic)%reset
    enddo

    !AJ Free up communicators created for this supercell
    if(self%sc_sub_not_full%nproc.gt.0) then
       call free_comm(self%sc_sub_not_full)
    end if
    if(self%sc_sub_all%nproc.gt.0 .and. self%sc_sub_all%id.ne.self%parent_sub%id) then
       call free_comm(self%sc_sub_all)
    end if
    if(self%sc_sub_pd%nproc.gt.0) then
       call free_comm(self%sc_sub_pd)
    end if

    !deallocate
    call self%deallocate

    !Could zero out variables but no real need
  end subroutine sc_reset

  !>Set the locality of each object
  subroutine sc_set_locality(self)
    implicit none
    class(supercell_type), intent(inout) :: self
    integer :: ic

    !Set the locality of children
    do ic=1,self%ncell
       call self%cells(ic)%set_locality
    enddo
    
    !Set our locality
    self%is_empty=all(self%cells%is_empty)
    self%is_local=any(self%cells%is_local)
    self%is_all_local=all(self%cells%is_all_local)

  end subroutine sc_set_locality

  !>Given an it value get the it of the left connected cell
  function sc_get_left_it(self,it)
    implicit none
    class(supercell_type),intent(in) :: self
    integer,intent(in)::it
    integer :: sc_get_left_it, ic, tmp
    !Default to no connection
    tmp=-1

    !Look for early exit
    if((.not.self%has_it(it)).or.(self%ncell.eq.1)) then
       tmp=-1
    else
       !Now see which cell has this it
       do ic=1,self%ncell
          if(self%cells(ic)%it_ind.eq.it) then
             tmp=ic-1
             exit
          endif
       enddo

       !Now set appropriate value
       if(tmp.le.0) then
          tmp=-1
       else
          tmp=self%cells(tmp)%it_ind
       endif

    endif

    !Update result
    sc_get_left_it=tmp
  end function sc_get_left_it

  !>Debug routine to dump the current supercell
  subroutine sc_dump(self,prefix)
    use mp, only: proc0
    implicit none
    class(supercell_type), intent(inout) :: self
    character(len=*), optional, intent(in) :: prefix
    character(len=80) :: fname
    integer :: lun=24
    complex, dimension(:,:), allocatable :: tmp

    !Make filename
    write(fname,'("sc_is_",I0,"_ik_",I0,".dat")') self%is_ind, self%ik_ind
    if(present(prefix)) fname=trim(prefix)//trim(fname)
    
    !Allocate array and fetch data
    allocate(tmp(self%ncol,self%nrow))

    !Fetch data
    call self%pull_rows_to_arr(tmp)

    !Now we've fetched data we can write
    if(proc0) then
       !Write
       open(unit=lun,file=fname,status='replace',form="unformatted")
       write(lun) tmp
       close(lun)
    endif

    !Free memory
    deallocate(tmp)
  end subroutine sc_dump

  !>Store the field equations at row level
  subroutine sc_store_fq(self,fq,fqa,fqp,ifl_in,it_in,ig_in)
    use theta_grid, only: ntgrid
    use run_parameters, only: fphi, fapar, fbpar
    implicit none
    class(supercell_type), intent(inout) :: self
    complex, dimension(:,:,:), intent(in) :: fq, fqa, fqp
    integer, intent(in) :: ifl_in, it_in, ig_in
    integer :: ic, ic_in, irow, ifq, ulim, ir

    !If we don't have anything in this supercell then exit
    if(self%is_empty) return

    !Find out which cell has our it
    ic_in=0
    do ic=1,self%ncell
       if(self%cells(ic)%it_ind.eq.it_in) then
          ic_in=ic
          exit
       endif
    enddo

    !Work out the row we want to put data in
    irow=(ig_in+ntgrid+1)+(ic_in-1)*(2*ntgrid)+(ifl_in-1)*self%nextend

    !Loop over cells
    do ic=1,self%ncell
       !If we don't have this cell cycle
       if(self%cells(ic)%is_empty) then
          cycle
       endif

       !Check if we have this row, if not cycle
       if(.not.self%cells(ic)%has_row(irow)) then
          cycle
       endif

       !Find upper limit of column
       ulim=self%cells(ic)%ncol
       
       ifq=0
       if(fphi.gt.epsilon(0.0))then
          !Increment counter
          ifq=ifq+1

          !Convert extended irow to local
          ir=self%cells(ic)%rb(ifq)%irex_to_ir(irow)

          !Store data
          self%cells(ic)%rb(ifq)%data(:,ir)=fq(:ulim,self%cells(ic)%it_ind,self%cells(ic)%ik_ind)
       endif

       if(fapar.gt.epsilon(0.0))then
          !Increment counter
          ifq=ifq+1

          !Convert extended irow to local
          ir=self%cells(ic)%rb(ifq)%irex_to_ir(irow)

          !Store data
          self%cells(ic)%rb(ifq)%data(:,ir)=fqa(:ulim,self%cells(ic)%it_ind,self%cells(ic)%ik_ind)
       endif

       if(fbpar.gt.epsilon(0.0))then
          !Increment counter
          ifq=ifq+1

          !Convert extended irow to local
          ir=self%cells(ic)%rb(ifq)%irex_to_ir(irow)

          !Store data
          self%cells(ic)%rb(ifq)%data(:,ir)=fqp(:ulim,self%cells(ic)%it_ind,self%cells(ic)%ik_ind)
       endif
    enddo

  end subroutine sc_store_fq

  !>Prepare the field matrix for calculating field updates
  subroutine sc_prepare(self,prepare_type)
    use mp, only: mp_abort
    implicit none
    class(supercell_type), intent(inout) :: self
    integer, intent(in) :: prepare_type
    character (len=40) :: errmesg

    !Exit early if we're empty
    if(self%is_empty) return
    if(.not.self%is_local) return

    !Need to decide what method we're using to prepare the
    !field matrix
    select case(prepare_type)
    case(0) !Invert
       call self%invert
    case default
       write(errmesg,'("ERROR: Invalid prepare_type : ",I0)') prepare_type
       call mp_abort(trim(errmesg))
    end select
  end subroutine sc_prepare

  !>A routine to invert the field matrix
  subroutine sc_invert(self)
    implicit none
    class(supercell_type), intent(inout) :: self

    if(self%is_empty) return
    if(.not.self%is_local) return
    
    !If we have all of the supercell local then use an explictly local
    !matrix inversion
    !Note could have a variable which forces local inversion
    if(self%is_all_local.or.pc%force_local_invert) then
       call self%invert_local
    else
       call self%invert_mpi
    endif

  end subroutine sc_invert

  !>A routine to invert the field matrix locally
  subroutine sc_invert_local(self)
    use mat_inv, only: inverse_gj
    implicit none
    class(supercell_type), intent(inout) :: self
    complex, dimension(:,:), allocatable :: tmp_arr

    !Exit early if possible
    if(self%is_empty) return
    if(.not.self%is_local) return

    !Allocate temporary array
    allocate(tmp_arr(self%ncol,self%nrow))

    !First have to pull all row level data
    call self%pull_rows_to_arr(tmp_arr)

    !Now invert in place
    call inverse_gj(tmp_arr,self%nrow)
    tmp_arr=transpose(tmp_arr)

    !Now push back to row level
    call self%push_arr_to_rows(tmp_arr)

    !Free memory
    deallocate(tmp_arr)

  end subroutine sc_invert_local

  !>A routine to invert the field matrix using mpi
  subroutine sc_invert_mpi(self)
    use mp, only: mp_abort
    implicit none
    class(supercell_type), intent(inout) :: self
    self%is_empty=.true. !Until implemented make everything be treated as empty.
    call mp_abort("ERROR: Invert with mpi not yet implemented.")
  end subroutine sc_invert_mpi

  !>A routine to collect all the row level data and store in passed array

  !>Gather the row blocks up for this cell to fill an array
  subroutine sc_pull_rows_to_arr(self,arr)
    !AJ Nonblocking?
    use mp, only: sum_allreduce_sub
    implicit none
    class(supercell_type), intent(in) :: self
    complex, dimension(:,:), intent(out) :: arr !Note: Dimensions should be self%ncol, self%nrow
    integer :: ic, rl,ru, cl,cu, irb

    !Now see if we want to exit
    if(self%is_empty) return
    if(.not.self%is_local) return

    !If we're all local then just do things simply as we don't need comms
    if(self%is_all_local) then
       do ic=1,self%ncell
          do irb=1,self%cells(ic)%nrb
             rl=self%cells(ic)%rb(irb)%row_llim
             ru=self%cells(ic)%rb(irb)%row_ulim
             cl=self%cells(ic)%rb(irb)%col_llim
             cu=self%cells(ic)%rb(irb)%col_ulim
             arr(cl:cu,rl:ru)=self%cells(ic)%rb(irb)%data
          enddo
       enddo
    else
       !Zero out arr |--> Not needed if gathering
       arr=0

       !Place local piece in arr
       do ic=1,self%ncell
          if(self%cells(ic)%is_empty) cycle
          do irb=1,self%cells(ic)%nrb
             rl=self%cells(ic)%rb(irb)%row_llim
             ru=self%cells(ic)%rb(irb)%row_ulim
             cl=self%cells(ic)%rb(irb)%col_llim
             cu=self%cells(ic)%rb(irb)%col_ulim
             arr(cl:cu,rl:ru)=self%cells(ic)%rb(irb)%data
          enddo
       enddo

       !Now use MPI to get the rest of the data from other procs
       !<DD>FOR NOW USE ALL_REDUCE AS EASIER, BUT SHOULD BE ABLE
       !TO DO THIS WITH ALL_GATHERV WHICH SHOULD BE FASTER
       if(self%sc_sub_pd%nproc.gt.0) call sum_allreduce_sub(arr,self%sc_sub_pd%id)
    endif

  end subroutine sc_pull_rows_to_arr

  !>A routine to distribute an array to appropriate row blocks
  subroutine sc_push_arr_to_rows(self,arr)
    implicit none
    class(supercell_type), intent(inout) :: self
    complex, dimension(self%nrow,self%ncol), intent(in) :: arr
    integer :: ic, ir
    integer :: rl,ru,cl,cu
    !If empty don't have anywhere to store data so exit
    if(self%is_empty) return

    !Store all local sections
    do ic=1,self%ncell
       if(self%cells(ic)%is_empty) cycle

       do ir=1,self%cells(ic)%nrb
          !Get row/col limits
          rl=self%cells(ic)%rb(ir)%row_llim
          ru=self%cells(ic)%rb(ir)%row_ulim
          cl=self%cells(ic)%rb(ir)%col_llim
          cu=self%cells(ic)%rb(ir)%col_ulim
          
          !Store data
          self%cells(ic)%rb(ir)%data=arr(cl:cu,rl:ru)
       enddo
    enddo
  end subroutine sc_push_arr_to_rows

  !>Create primary (top level) sub communicators
  subroutine sc_make_subcom_1(self)
    use mp, only: comm_type, split, free_comm
    implicit none
    class(supercell_type), intent(inout) :: self
    integer :: ic, colour
    type(comm_type) :: tmp

    !/First make a subcommunicator involving any procs
    !which can see this supercell
    colour=0
    if(self%is_local)colour=1
    call split(colour,tmp,self%parent_sub%id)
    !Note we only store the subcom for those procs which will use it
    !this ensures an error will occur if we try to use the subcom in
    !the wrong place
    if(self%is_local)then
       !Here we check if the new subcom has same number of members as
       !parent, if so we can carry on using the parent instead
       if(tmp%nproc.eq.self%parent_sub%nproc) then
          self%sc_sub_all=self%parent_sub
          call free_comm(tmp)
       else
          self%sc_sub_all=tmp
       endif
    else
       !If we can't see any of this object then 
       !Destroy communicator
       call free_comm(tmp)
       !and return
       return
    endif

    !Now make the supercell sub communicators
    do ic=1,self%ncell
       !Set parent subcom
       self%cells(ic)%parent_sub=self%sc_sub_all
    enddo

  end subroutine sc_make_subcom_1

  !>Create the secondary (intraobject) subcommunicators
  subroutine sc_make_subcom_2(self)
    use mp, only: comm_type, split, free_comm, sum_allreduce_sub, proc0, iproc
    implicit none
    class(supercell_type), intent(inout) :: self
    integer :: colour, key, cc
    type(comm_type) :: tmp
    logical :: some_empty, some_pd, some_all, head_notfull
    !Exit if we don't have this supercell
    if(.not.self%is_local)return

    !/Work out if we have any empty procs
    some_empty=.false.
    colour=0
    if(self%is_empty)colour=1
    call sum_allreduce_sub(colour,self%sc_sub_all%id)
    if(colour.ne.0) some_empty=.true.

    !/Work out if we have any partial data procs
    some_pd=.false.
    colour=0
    if(.not.(self%is_empty.or.self%is_all_local))colour=1
    call sum_allreduce_sub(colour,self%sc_sub_all%id)
    if(colour.ne.0) some_pd=.true.

    !/Work out if we have any all local procs
    some_all=.false.
    colour=0
    if(self%is_all_local)colour=1
    call sum_allreduce_sub(colour,self%sc_sub_all%id)
    if(colour.ne.0) some_all=.true.

    !/Now make a subcommunicator involving all procs
    !which have part of the supercell data range (but not all of it)
    if(some_pd)then
       colour=0
       if((.not.(self%is_empty.or.self%is_all_local))) colour=1
       call split(colour,tmp,self%sc_sub_all%id)
       if(colour.eq.0)then
          !/Destroy communicator
          call free_comm(tmp)
       else
          !/Store communicator
          self%sc_sub_pd=tmp
       endif
    endif

    !Decide if we're the head proc
    !Note: We also decide which proc is in charge of filling the empty
    !procs with data. By allowing this to be different from the head proc
    !we have a better chance of being able to overlap comms and computations
    !Note: We try to make proc0 the head proc if it is a member of the supercell
    !which has data as it means we don't need to send this supercell to anyone.
    if(some_all)then
       !Here we make a temporary communicator for the all_local procs
       cc=0
       if(self%is_all_local) cc=1
       !/Note we don't destroy this until later
       call split(cc,tmp,self%sc_sub_all%id)

       !Now find out if we have proc0 in our new comm
       if(self%is_all_local)then
          cc=0
          if(proc0) cc=1
          call sum_allreduce_sub(cc,tmp%id)

          !Decide if we're the head, prefer proc0 to be head
          !if it's a member of our group
          if(cc.ne.0)then
             self%is_head=proc0
          else
             self%is_head=tmp%proc0
          endif

          !Here we decide which proc should be in charge of
          !informing the empty procs of the result
          head_notfull=(min(1,tmp%nproc-1).eq.tmp%iproc)
       else
          head_notfull=.false.
       endif

       !Free the communicator
       call free_comm(tmp)
    else
       if(self%sc_sub_pd%nproc.gt.0)then
          cc=0
          if(proc0) cc=1
          call sum_allreduce_sub(cc,self%sc_sub_pd%id)

          !Decide if we're the head, prefer proc0 to be head
          !if it's a member of our group
          if(cc.ne.0)then
             self%is_head=proc0
          else
             self%is_head=self%sc_sub_pd%proc0
          endif
       endif

       !Here we decide which proc should be in charge of
       !informing the empty procs of the result
       head_notfull=(min(1,self%sc_sub_pd%nproc-1).eq.self%sc_sub_pd%iproc)
    endif

    !Store the head_iproc
    self%head_iproc_world=0
    if(self%is_head) self%head_iproc_world=iproc
    call sum_allreduce_sub(self%head_iproc_world,self%sc_sub_all%id)

    self%head_iproc=0
    if(self%is_head) self%head_iproc=self%sc_sub_all%iproc
    call sum_allreduce_sub(self%head_iproc,self%sc_sub_all%id)

    self%head_iproc_pd=0
    if(self%is_head) self%head_iproc_pd=self%sc_sub_pd%iproc
    call sum_allreduce_sub(self%head_iproc_pd,self%sc_sub_all%id)

    !/Now make a subcommunicator involving all empty procs and
    !head pd proc
    if(some_empty)then
       !/Default values
       colour=0
       key=2

       !/Join if we're empty and NOT proc0
       !This is done as we handle the filling of proc0
       !separately and by actively excluding it we aim to
       !prevent stalls waiting for proc0 to catch up
       if(self%is_empty.and.(.not.proc0)) colour=1

       if(head_notfull) then
          colour=1
          key=0
       endif

       !Here we use key to make sure the proc with data is proc0
       call split(colour,key,tmp,self%sc_sub_all%id)
       if(colour.eq.0)then
          !/Free communicator
          call free_comm(tmp)
       else
          if(tmp%nproc.gt.1)then
             !/Store communicator
             self%sc_sub_not_full=tmp
             
             !/Allocate array used to store non blocking
             !comm handles
             if(tmp%proc0) then
                allocate(self%nb_req_hand(tmp%nproc-1))
             else
                allocate(self%nb_req_hand(1))
             endif
          else
             call free_comm(tmp)
          endif
       endif

    endif

  end subroutine sc_make_subcom_2

!------------------------------------------------------------------------

!================
!KYBLOCK
!================
  !>Allocate storage space
  subroutine ky_allocate(self)
    implicit none
    class(ky_type), intent(inout) :: self
    integer :: is
    do is=1,self%nsupercell
       call self%supercells(is)%allocate
    enddo
  end subroutine ky_allocate

  !>Deallocate storage space
  subroutine ky_deallocate(self)
    implicit none
    class(ky_type), intent(inout) :: self
    if(allocated(self%supercells)) deallocate(self%supercells)
  end subroutine ky_deallocate

  !>Debug printing
  subroutine ky_debug_print(self)
    implicit none
    class(ky_type), intent(in) :: self
    write(dlun,'("Ky block debug print. Index ik=",I0," nsupercell=",I0)') &
         self%ik_ind, self%nsupercell
  end subroutine ky_debug_print

  !>Get the field update for this ik
  subroutine ky_get_field_update(self,fq,fqa,fqp)
    implicit none
    class(ky_type), intent(inout) :: self
    complex, dimension(:,:), intent(in) :: fq,fqa,fqp
    integer :: is

    !First trigger calculation of field update.
    !Note: Empty supercells don't do anything
    do is=1,self%nsupercell
       !/Calculate the field update
       call self%supercells(is)%get_field_update(fq,fqa,fqp)
    enddo
  end subroutine ky_get_field_update

  !>Store the field equation at row level
  subroutine ky_store_fq(self,fq,fqa,fqp,ifl_in,it_in,ig_in)
    implicit none
    class(ky_type), intent(inout) :: self
    complex, dimension(:,:,:), intent(in) :: fq, fqa, fqp
    integer, intent(in) :: ifl_in, it_in, ig_in
    integer :: is

    !If we don't have anything in this supercell then exit
    if(self%is_empty) return

    !Delegate work to supercell
    do is=1,self%nsupercell
       if(.not.self%supercells(is)%has_it(it_in)) cycle
       call self%supercells(is)%store_fq(fq,fqa,fqp,ifl_in,it_in,ig_in)
    enddo
  end subroutine ky_store_fq

  !>Given it say what the supercell id is
  function ky_is_from_it(self,it)
    implicit none
    class(ky_type), intent(inout) :: self
    integer, intent(in) :: it
    integer :: ky_is_from_it
    integer :: is
    ky_is_from_it=-1
    do is=1,self%nsupercell
       if(self%supercells(is)%has_it(it))then
          ky_is_from_it=is
          exit
       endif
    enddo
  end function ky_is_from_it

  !>A routine to reset the object
  subroutine ky_reset(self)
    use mp, only: free_comm, mp_comm_null
    implicit none
    class(ky_type), intent(inout) :: self
    integer :: is

    !Call reset on all children
    do is=1,self%nsupercell
       call self%supercells(is)%reset
    enddo


    !AJ Free communicators associated with this ky block
    if(self%ky_sub_all%id .ne. self%parent_sub%id .and. self%ky_sub_all%nproc .gt. 0 .and. self%ky_sub_all%id .ne. mp_comm_null) then
       call free_comm(self%ky_sub_all)
    end if

    !deallocate
    call self%deallocate

    !Could zero out variables but no real need
  end subroutine ky_reset

  !>Set the locality of each object
  subroutine ky_set_locality(self)
    implicit none
    class(ky_type), intent(inout) :: self
    integer :: is

    !Set the locality of children
    do is=1,self%nsupercell
       call self%supercells(is)%set_locality
    enddo
    
    !Set our locality
    self%is_empty=all(self%supercells%is_empty)
    self%is_local=any(self%supercells%is_local)
    self%is_all_local=all(self%supercells%is_all_local)

  end subroutine ky_set_locality

  !>Prepare the field matrix for calculating field updates
  subroutine ky_prepare(self,prepare_type)
    implicit none
    class(ky_type), intent(inout) :: self
    integer, intent(in) :: prepare_type
    integer :: is

    !Exit early if we're empty
    if(self%is_empty) return
    if(.not.self%is_local) return

    !Tell each supercell to prepare
    do is=1,self%nsupercell
       call self%supercells(is)%prepare(prepare_type)
    enddo
  end subroutine ky_prepare

  !>Create the primary subcommunicators
  subroutine ky_make_subcom_1(self)
    use mp, only: comm_type, split, free_comm
    implicit none
    class(ky_type), intent(inout) :: self
    integer :: is, colour
    type(comm_type) :: tmp

    !/First make a subcommunicator involving any procs
    !which can see this ky block
    colour=0
    if(self%is_local)colour=1
    call split(colour,tmp,self%parent_sub%id)
    !Note we only store the subcom for those procs which will use it
    !this ensures an error will occur if we try to use the subcom in
    !the wrong place
    if(self%is_local)then
       if(tmp%nproc.eq.self%parent_sub%nproc) then
          self%ky_sub_all=self%parent_sub
          call free_comm(tmp)
       else
          self%ky_sub_all=tmp
       endif
    else
       !If we can't see any of this object then
       !/Destroy communicator
       call free_comm(tmp)
       !and return
       return
    endif

    !Now make the supercell sub communicators
    do is=1,self%nsupercell
       !Set parent subcom
       self%supercells(is)%parent_sub=self%ky_sub_all

       !Now make child subcom
       call self%supercells(is)%make_subcom_1
    enddo
  end subroutine ky_make_subcom_1

  !>Create the secondary subcommunicators
  subroutine ky_make_subcom_2(self)
    implicit none
    class(ky_type), intent(inout) :: self
    integer :: is

    !Exit if we don't have this kyb
    if(.not.self%is_local)return

    !Now make the supercell sub communicators
    do is=1,self%nsupercell
       !Now make child subcom
       call self%supercells(is)%make_subcom_2
    enddo
  end subroutine ky_make_subcom_2

!------------------------------------------------------------------------

!================
!FIELDMAT
!================
  !>Allocate storage space
  subroutine fm_allocate(self)
    implicit none
    class(fieldmat_type), intent(inout) :: self
    integer :: ik
    allocate(self%heads(self%ntheta0, self%naky))
    if(.not.self%is_local) return
    do ik=1,self%naky
       call self%kyb(ik)%allocate
    enddo
  end subroutine fm_allocate

  !>Deallocate storage space
  subroutine fm_deallocate(self)
    implicit none
    class(fieldmat_type), intent(inout) :: self
    if(allocated(self%kyb)) deallocate(self%kyb)
    if(allocated(self%heads)) deallocate(self%heads)
  end subroutine fm_deallocate

  !>Debug printing
  subroutine fm_debug_print(self)
    implicit none
    class(fieldmat_type), intent(in) :: self
    write(dlun,'("Field debug print. naky=",I0)') &
         self%naky
  end subroutine fm_debug_print

  !>Just work out the locality (also sets is_empty etc. but is not intended for this)
  subroutine fm_set_is_local(self)
    implicit none
    class(fieldmat_type), intent(inout) :: self
    integer :: ik

    !Set the locality of children
    do ik=1,self%naky
       call self%kyb(ik)%set_locality
    enddo
    
    !Set our locality
    self%is_empty=all(self%kyb%is_empty)
    self%is_local=any(self%kyb%is_local)
    self%is_all_local=all(self%kyb%is_local)
  end subroutine fm_set_is_local

  !>Initialise the field objects
  subroutine fm_init(self)
    use kt_grids, only: naky, ntheta0
    use dist_fn, only: itright, get_leftmost_it
    use run_parameters, only: fphi, fapar, fbpar
    use theta_grid, only: ntgrid
    implicit none
    class(fieldmat_type), intent(inout) :: self
    integer :: it, cur_id, ifq
    integer :: it_start, is_count
    integer :: ik, is, ic
    integer, dimension(:), allocatable :: nsuper
    integer, dimension(:), allocatable :: tmp_supercell_ids
    integer, dimension(:,:), allocatable :: supercell_ids
    integer, dimension(:,:), allocatable :: ncells
    integer, dimension(:,:), allocatable :: is_to_itmin

    !Count the fields
    nfield=0
    if(fphi.gt.epsilon(0.0)) nfield=nfield+1
    if(fapar.gt.epsilon(0.0)) nfield=nfield+1
    if(fbpar.gt.epsilon(0.0)) nfield=nfield+1
    nfq=nfield

    !Allocate temp arrays
    allocate(nsuper(naky))
    allocate(tmp_supercell_ids(ntheta0))
    allocate(supercell_ids(ntheta0,naky))
    allocate(ncells(ntheta0,naky))
    allocate(is_to_itmin(ntheta0,naky))

    !Initialise arrays
    nsuper=0
    ncells=0
    is_to_itmin=0
    supercell_ids=-1

    !Find the leftmost it for each supercell
    do ik=1,naky
       do it=1,ntheta0
          supercell_ids(it,ik)=get_leftmost_it(it,ik)
       enddo
    enddo

    !Now work out how many supercells there are for each
    !ik and what their properties are.
    do ik=1,naky
       !Store the leftmost it values for current ik
       tmp_supercell_ids=supercell_ids(:,ik)
       
       !Now loop over all it values and count how many
       !have which it as the leftmost one.
       it_start=0
       is_count=0
       do while(sum(tmp_supercell_ids).ne.-1*ntheta0)
          !Increase the it value
          it_start=it_start+1

          !Get the leftmost it for supercell with it_start in it
          cur_id=tmp_supercell_ids(it_start)

          !If we've seen this supercell then cycle
          if(cur_id.eq.-1) cycle

          !Increase supercell counter
          nsuper(ik)=nsuper(ik)+1
          is_count=is_count+1

          !Store the leftmost it value
          is_to_itmin(is_count,ik)=cur_id

          !Set all matching entries to -1 so we know to skip
          !this supercell again
          do it=it_start,ntheta0
             if(tmp_supercell_ids(it).eq.cur_id) then
                tmp_supercell_ids(it)=-1

                !Increase cell counter
                ncells(cur_id,ik)=ncells(cur_id,ik)+1
             endif
          enddo
       enddo
    enddo
    
    !Free memory
    deallocate(tmp_supercell_ids)

    !Now we want to allocate the basic field structure
    !Store the basic properties
    self%naky=naky
    self%ntheta0=ntheta0
    self%npts=naky*ntheta0*(2*ntgrid+1)
    self%nbound=0 !Calculated whilst we create objects
    allocate(self%kyb(naky))
    !/Make the kyb objects
    do ik=1,naky
       !Store the basic properties
       self%kyb(ik)%ik_ind=ik
       self%kyb(ik)%nsupercell=nsuper(ik)

       !Allocate the supercell list
       allocate(self%kyb(ik)%supercells(nsuper(ik)))

       !/Now make the supercell objects
       do is=1,nsuper(ik)

          !Store the basic properties
          self%kyb(ik)%supercells(is)%is_ind=is
          self%kyb(ik)%supercells(is)%ik_ind=ik
          self%kyb(ik)%supercells(is)%it_left=is_to_itmin(is,ik)
          self%kyb(ik)%supercells(is)%ncell=ncells(is_to_itmin(is,ik),ik)
          self%kyb(ik)%supercells(is)%nextend=1+2*ntgrid*self%kyb(ik)%supercells(is)%ncell
          self%kyb(ik)%supercells(is)%nrow=self%kyb(ik)%supercells(is)%nextend*nfield
          self%kyb(ik)%supercells(is)%ncol=self%kyb(ik)%supercells(is)%nrow

          !Now allocate the cell list
          allocate(self%kyb(ik)%supercells(is)%cells(self%kyb(ik)%supercells(is)%ncell))
          allocate(self%kyb(ik)%supercells(is)%initialised(self%kyb(ik)%supercells(is)%nextend))
          self%kyb(ik)%supercells(is)%initialised=.false.

          !Get the first it
          it=supercell_ids(is_to_itmin(is,ik),ik)

          !/Now make the cell objects
          do ic=1,ncells(is_to_itmin(is,ik),ik)
             
             !Store the basic properties
             self%kyb(ik)%supercells(is)%cells(ic)%nrb=nfq
             self%kyb(ik)%supercells(is)%cells(ic)%ic_ind=ic
             self%kyb(ik)%supercells(is)%cells(ic)%it_ind=it
             self%kyb(ik)%supercells(is)%cells(ic)%ik_ind=ik
             self%kyb(ik)%supercells(is)%cells(ic)%is_ind=is
             self%kyb(ik)%supercells(is)%cells(ic)%ignore_boundary=(ic.lt.ncells(is_to_itmin(is,ik),ik)) !We ignore the boundary for all cells except the rightmost
             self%kyb(ik)%supercells(is)%cells(ic)%row_llim=1
             self%kyb(ik)%supercells(is)%cells(ic)%row_ulim=self%kyb(ik)%supercells(is)%nrow
             self%kyb(ik)%supercells(is)%cells(ic)%nrow=self%kyb(ik)%supercells(is)%nrow
             allocate(self%kyb(ik)%supercells(is)%cells(ic)%col_llim(nfq))
             allocate(self%kyb(ik)%supercells(is)%cells(ic)%col_ulim(nfq))
             if(.not.self%kyb(ik)%supercells(is)%cells(ic)%ignore_boundary)then
                self%kyb(ik)%supercells(is)%cells(ic)%ncol=2*ntgrid+1
             else
                self%kyb(ik)%supercells(is)%cells(ic)%ncol=2*ntgrid
                self%nbound=self%nbound+1
             endif
             do ifq=1,nfq
                self%kyb(ik)%supercells(is)%cells(ic)%col_llim(ifq)=1+(ic-1)*(2*ntgrid)+(ifq-1)*self%kyb(ik)%supercells(is)%nrow
                self%kyb(ik)%supercells(is)%cells(ic)%col_ulim(ifq)=&
                     self%kyb(ik)%supercells(is)%cells(ic)%col_llim(ifq)+&
                     self%kyb(ik)%supercells(is)%cells(ic)%ncol-1
             enddo
             self%kyb(ik)%supercells(is)%cells(ic)%ncol_tot=self%kyb(ik)%supercells(is)%cells(ic)%ncol*nfq

             !Now allocate row blocks
             allocate(self%kyb(ik)%supercells(is)%cells(ic)%rb(self%kyb(ik)%supercells(is)%cells(ic)%nrb))

             !/Now make row block objects
             do ifq=1,nfq
                !Store basic properties
                self%kyb(ik)%supercells(is)%cells(ic)%rb(ifq)%ir_ind=ifq
                self%kyb(ik)%supercells(is)%cells(ic)%rb(ifq)%is_ind=is
                self%kyb(ik)%supercells(is)%cells(ic)%rb(ifq)%ik_ind=ik
                self%kyb(ik)%supercells(is)%cells(ic)%rb(ifq)%ic_ind=ic
                self%kyb(ik)%supercells(is)%cells(ic)%rb(ifq)%it_ind=it
                self%kyb(ik)%supercells(is)%cells(ic)%rb(ifq)%ncol=self%kyb(ik)%supercells(is)%cells(ic)%ncol
                self%kyb(ik)%supercells(is)%cells(ic)%rb(ifq)%col_llim=1+&
                     (ic-1)*(2*ntgrid)+&
                     (ifq-1)*self%kyb(ik)%supercells(is)%nextend
                self%kyb(ik)%supercells(is)%cells(ic)%rb(ifq)%col_ulim=&
                     self%kyb(ik)%supercells(is)%cells(ic)%rb(ifq)%col_llim+&
                     self%kyb(ik)%supercells(is)%cells(ic)%rb(ifq)%ncol-1
             enddo

             !Update it
             if(ic.lt.ncells(is_to_itmin(is,ik),ik)) it=itright(ik,it)

          enddo
       enddo
    enddo

    !Free memory
    deallocate(nsuper,supercell_ids,ncells,is_to_itmin)
  end subroutine fm_init

  !>Find the response of g to delta-fn field perturbations 
  !and store at the row level
  subroutine fm_populate(self)
    use run_parameters, only: fphi, fapar, fbpar
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use dist_fn_arrays, only: g
    use kt_grids, only: kwork_filter
    implicit none
    class(fieldmat_type), intent(inout) :: self
    integer :: ifl, pts_remain, ik, is

    !First initialise everything to 0
    !AJ This is also done in fields.f90 when arrays are created so we should be 
    !AJ able to remove from here.
    g=0
    phi=0.0
    phinew=0.0
    apar=0.0
    aparnew=0.0
    bpar=0.0
    bparnew=0.0

    !Initialise field counter
    ifl=0

    !Do phi
    if(fphi.gt.epsilon(0.0)) then
       !Set the number of points left
       pts_remain=self%npts-self%nbound

       !Increment the field counter
       ifl=ifl+1

       !Reset the filter array
       kwork_filter=.false.

       !Reset the init state arrays
       self%kyb%initdone=.false.
       do ik=1,self%naky
          self%kyb(ik)%supercells%initdone=.false.
          do is=1,self%kyb(ik)%nsupercell
             self%kyb(ik)%supercells(is)%initialised=.false.
          enddo
       enddo

       !Now loop over all points and initialise
       do while(pts_remain.gt.0)
          call self%init_next_field_points(phinew,pts_remain,kwork_filter,ifl)
       enddo
    endif

    !Do apar
    if(fapar.gt.epsilon(0.0)) then
       !Set the number of points left
       pts_remain=self%npts-self%nbound

       !Increment the field counter
       ifl=ifl+1

       !Reset the filter array
       kwork_filter=.false.

       !Reset the init state arrays
       self%kyb%initdone=.false.
       do ik=1,self%naky
          self%kyb(ik)%supercells%initdone=.false.
          do is=1,self%kyb(ik)%nsupercell
             self%kyb(ik)%supercells(is)%initialised=.false.
          enddo
       enddo

       !Now loop over all points and initialise
       do while(pts_remain.gt.0)
          call self%init_next_field_points(aparnew,pts_remain,kwork_filter,ifl)
       enddo
    endif

    !Do bpar
    if(fbpar.gt.epsilon(0.0)) then
       !Set the number of points left
       pts_remain=self%npts-self%nbound

       !Increment the field counter
       ifl=ifl+1

       !Reset the filter array
       kwork_filter=.false.

       !Reset the init state arrays
       self%kyb%initdone=.false.
       do ik=1,self%naky
          self%kyb(ik)%supercells%initdone=.false.
          do is=1,self%kyb(ik)%nsupercell
             self%kyb(ik)%supercells(is)%initialised=.false.
          enddo
       enddo

       !Now loop over all points and initialise
       do while(pts_remain.gt.0)
          call self%init_next_field_points(bparnew,pts_remain,kwork_filter,ifl)
       enddo
    endif
    
    !Reset the filter array
    kwork_filter=.false.
  end subroutine fm_populate

  !>Initialise the next set of delta functions, find the response
  !and store in the appropriate rows.
  subroutine fm_init_next_field_points(self,field,pts_remain,kwork_filter,ifl)
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use dist_fn, only: getfieldeq, getfieldeq_nogath, timeadv
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    use mp, only: barrier
    implicit none
    class(fieldmat_type), intent(inout) :: self
    complex, dimension(-ntgrid:ntgrid,ntheta0,naky), intent(inout) :: field
    integer, intent(inout) :: pts_remain
    integer, intent(in) :: ifl
    logical, dimension(ntheta0, naky), intent(inout) :: kwork_filter
    complex, dimension(:,:,:), allocatable :: fq, fqa, fqp
    integer :: ik, is, iex, ic, iex_init, ig, it
    integer :: it_cur, ig_cur, found_it, found_ig
    !Allocate the fieldeq arrays, should look at just allocating
    !these things once (maybe store at fieldmat level)
    allocate(fq(-ntgrid:ntgrid, ntheta0, naky))
    allocate(fqa(-ntgrid:ntgrid, ntheta0, naky))
    allocate(fqp(-ntgrid:ntgrid, ntheta0, naky))

    !AJ This is a nasty barrier, but it is in here as on ARCHER we have a problem
    !AJ that we get deadlock on large core counts without it.  This has been investigated 
    !AJ and looks like a problem with the Cray MPI library, but it's hard to get to the bottom of.
    !AJ The problem does not occur on small core counts (we see it between 2000 and 3000 cores)
    call barrier()
    
    !Loop over all ik
    do ik=1,self%naky
       !Technically should be able to skip ik values which are not local (though could be empty)
       !/would have to remember to decrease pts_remain

       !If we've finished this ky set the 
       !filter and move on to next ik
       if(self%kyb(ik)%initdone)then
          kwork_filter(:,ik)=.true.
          cycle
       endif

       !Now look at all supercells with this ik
       do is=1,self%kyb(ik)%nsupercell
          !If we've finished with this supercell
          !then set the filter and move on to next
          if(self%kyb(ik)%supercells(is)%initdone)then
             !Set filter for all it on supercell
             do ic=1,self%kyb(ik)%supercells(is)%ncell
                kwork_filter(self%kyb(ik)%supercells(is)%cells(ic)%it_ind,ik)=.true.
             enddo
             cycle
          endif

          !//If we get to this point then we know that we have at least
          !one grid point on supercell for which we've not made a delta-fn

          !Now look for the first point that we've not initialised
          iex_init=0
          do iex=1,self%kyb(ik)%supercells(is)%nextend
             !If we've already done this point then cycle
             if(self%kyb(ik)%supercells(is)%initialised(iex)) cycle
             !Else this is the first point so store it and exit
             iex_init=iex
             exit
          enddo

          !If this is the last point then mark supercell as done
          !but don't set filter yet, this will be automatically
          !applied on next loop
          if(iex_init.eq.self%kyb(ik)%supercells(is)%nextend) then
             self%kyb(ik)%supercells(is)%initdone=.true.
          endif

          !Mark this point as initialised
          self%kyb(ik)%supercells(is)%initialised(iex_init)=.true.

          !Decrease the points remaining count
          pts_remain=pts_remain-1

          !Get the ig and it indices corresponding to this point
          call self%kyb(ik)%supercells(is)%iex_to_ig(iex_init,ig)
          call self%kyb(ik)%supercells(is)%iex_to_it(iex_init,it)

          !Initialise field
          field(ig,it,ik)=1.0

          !Set the duplicate boundary point
          if((ig.eq.-ntgrid).and.(it.ne.self%kyb(ik)%supercells(is)%it_left)) then
             !Get the left connection
             it=self%kyb(ik)%supercells(is)%get_left_it(it)
             field(ntgrid,it,ik)=1.0
          endif
       enddo !Ends loop over supercells

       !Decide if we've done this kyb
       self%kyb(ik)%initdone=all(self%kyb(ik)%supercells%initdone)

    enddo !Ends loop over ik

    !Now do the time advance
    call timeadv(phi,apar,bpar,phinew,aparnew,bparnew,0)
    !Now get the field equations
    call self%getfieldeq_nogath(phinew,aparnew,bparnew,fq,fqa,fqp)

    !Now we need to store the field equations in the appropriate places
    do ik=1,self%naky
       !If no data then skip ik
       if(self%kyb(ik)%is_empty) cycle

       !Loop until we've done all supercells
       do while(sum(field(:,:,ik)).ne.0)
          !Reset variables
          ig_cur=-ntgrid-1
          it_cur=0
          found_ig=0
          found_it=0
          
          !Look for first entry which was set to a delta-fn
          do it=1,ntheta0
             do ig=-ntgrid,ntgrid
                !If this point wasn't touched ignore
                if(field(ig,it,ik).eq.0) cycle
                ig_cur=ig
                found_ig=1
                exit
             enddo
             if(found_ig.ne.0) then
                it_cur=it
                found_it=1
                exit
             endif
          enddo

          !Now check we've found a point, else cycle
          !NOTE: If this triggers we're probably in trouble
          !as nothing will be different on the next iteration
          !and hence we will be stuck in an infinite loop.
          !As such this should really trigger an abort.
          if(found_it.eq.0) cycle

          !Zero out the field at this point so we ignore
          !it on future iterations
          field(ig_cur,it_cur,ik)=0.0

          !Find the appropriate supercell
          is=self%kyb(ik)%is_from_it(it_cur)

          !Also zero out the boundary point if needed
          if((ig_cur.eq.-ntgrid).and.(it_cur.ne.self%kyb(ik)%supercells(is)%it_left))then
             it=self%kyb(ik)%supercells(is)%get_left_it(it_cur)
             field(ntgrid,it,ik)=0.0
          endif

          !Now store data
          call self%kyb(ik)%store_fq(fq,fqa,fqp,ifl,it_cur,ig_cur)
       enddo !Ends loop over points
    enddo !Ends loop over ky

    !Make sure field has been reset
    field=0

    !Deallocate field arrays
    deallocate(fq,fqa,fqp)
  end subroutine fm_init_next_field_points

  !>Prepare the field matrix for calculating field updates
  subroutine fm_prepare(self)
    implicit none
    class(fieldmat_type), intent(inout) :: self
    integer :: ik

    !Exit early if we're empty
    if(self%is_empty) return
    if(.not.self%is_local) return

    !Tell each ky to prepare
    do ik=1,self%naky
       call self%kyb(ik)%prepare(self%prepare_type)
    enddo
  end subroutine fm_prepare

  !>A routine to calculate the update to the fields
  subroutine fm_get_field_update(self, phi, apar, bpar)
    use kt_grids, only: ntheta0, naky
    use theta_grid, only: ntgrid
    use dist_fn, only: getfieldeq, getfieldeq_nogath
    implicit none
    class(fieldmat_type), intent(inout) :: self
    complex, dimension(:,:,:), intent(inout) :: phi,apar,bpar
    complex, dimension(:,:,:), allocatable :: fq, fqa, fqp
    integer :: ik, is

    do ik=1,self%naky
       if(.not.self%kyb(ik)%is_local) cycle
       do is=1,self%kyb(ik)%nsupercell
          if(.not.self%kyb(ik)%supercells(is)%is_local) cycle
          self%kyb(ik)%supercells(is)%tmp_sum=0.0
       enddo
    enddo

    !Allocate storage, may want to put these inside fieldmat
    !object to avoid continual allocation/deallocation
    allocate(fq(-ntgrid:ntgrid,ntheta0,naky))
    allocate(fqa(-ntgrid:ntgrid,ntheta0,naky))
    allocate(fqp(-ntgrid:ntgrid,ntheta0,naky))

    !First get the field eq
    call self%getfieldeq_nogath(phi,apar,bpar,fq,fqa,fqp)
    !Now get the field update for each ik
    do ik=1,self%naky
       !Skip non-local ky
       if(.not.self%kyb(ik)%is_local) cycle
       !Trigger the supercell field updates
       !AJ Non-blocking can be done here
       call self%kyb(ik)%get_field_update(fq(:,:,ik),fqa(:,:,ik),fqp(:,:,ik))
    enddo

    !Free memory
    deallocate(fq,fqa,fqp)

  end subroutine fm_get_field_update

  !>A routine to unpack the supercell tmp_sum vectors to
  !full field arrays
  subroutine fm_unpack_to_field(self,ph,ap,bp)
    use theta_grid, only: ntgrid
    use run_parameters, only: fphi,fapar,fbpar
    use mp, only: proc0
    implicit none
    class(fieldmat_type), intent(inout) :: self
    complex,dimension(-ntgrid:,:,:), intent(inout) :: ph,ap,bp
    integer :: ik,is,ic,iex,ig,it,bnd,ifl

    !We don't need to unpack if we're proc0 as we've done this already
    if(.not.proc0)then
       do ik=1,self%naky
          if(.not.self%kyb(ik)%is_local) cycle
          do is=1,self%kyb(ik)%nsupercell
             if(.not.self%kyb(ik)%supercells(is)%is_local) cycle
             
             !Finalise the receives
!             call self%kyb(ik)%supercells(is)%gfu_finish_recv
             
             iex=0
             ifl=0
             if(fphi.gt.0)then
                ifl=ifl+1
                do ic=1,self%kyb(ik)%supercells(is)%ncell
                   it=self%kyb(ik)%supercells(is)%cells(ic)%it_ind
                   bnd=0
                   if(ic.ne.self%kyb(ik)%supercells(is)%ncell) bnd=1
                   do ig=-ntgrid,ntgrid-bnd
                      iex=iex+1
                      ph(ig,it,ik)=self%kyb(ik)%supercells(is)%tmp_sum(iex)
                   enddo
                enddo
                
                !/Fix boundary points
                if(self%kyb(ik)%supercells(is)%ncell.gt.1) then
                   do ic=1,self%kyb(ik)%supercells(is)%ncell-1
                      ph(ntgrid,self%kyb(ik)%supercells(is)%cells(ic)%it_ind,ik)=&
                           ph(-ntgrid,self%kyb(ik)%supercells(is)%cells(ic+1)%it_ind,ik)
                   enddo
                endif
                
             endif
             
             if(fapar.gt.0)then
                ifl=ifl+1
                do ic=1,self%kyb(ik)%supercells(is)%ncell
                   it=self%kyb(ik)%supercells(is)%cells(ic)%it_ind
                   bnd=0
                   if(ic.ne.self%kyb(ik)%supercells(is)%ncell) bnd=1
                   do ig=-ntgrid,ntgrid-bnd
                      iex=iex+1
                      ap(ig,it,ik)=self%kyb(ik)%supercells(is)%tmp_sum(iex)
                   enddo
                enddo
                
                !/Fix boundary points
                if(self%kyb(ik)%supercells(is)%ncell.gt.1) then
                   do ic=1,self%kyb(ik)%supercells(is)%ncell-1
                      ap(ntgrid,self%kyb(ik)%supercells(is)%cells(ic)%it_ind,ik)=&
                           ap(-ntgrid,self%kyb(ik)%supercells(is)%cells(ic+1)%it_ind,ik)
                   enddo
                endif
                
             endif
             
             if(fbpar.gt.0)then
                ifl=ifl+1
                do ic=1,self%kyb(ik)%supercells(is)%ncell
                   it=self%kyb(ik)%supercells(is)%cells(ic)%it_ind
                   bnd=0
                   if(ic.ne.self%kyb(ik)%supercells(is)%ncell) bnd=1
                   do ig=-ntgrid,ntgrid-bnd
                      iex=iex+1
                      bp(ig,it,ik)=self%kyb(ik)%supercells(is)%tmp_sum(iex)
                   enddo
                enddo
                
                !/Fix boundary points
                if(self%kyb(ik)%supercells(is)%ncell.gt.1) then
                   do ic=1,self%kyb(ik)%supercells(is)%ncell-1
                      bp(ntgrid,self%kyb(ik)%supercells(is)%cells(ic)%it_ind,ik)=&
                           bp(-ntgrid,self%kyb(ik)%supercells(is)%cells(ic+1)%it_ind,ik)
                   enddo
                endif
                
             endif
          enddo
       enddo
    endif

  end subroutine fm_unpack_to_field

  !>A routine to reset the object
  subroutine fm_reset(self)
    use mp, only: free_comm
    implicit none
    class(fieldmat_type), intent(inout) :: self
    integer :: ik

    !Call reset on all children
    do ik=1,self%naky
       call self%kyb(ik)%reset
    enddo

    !deallocate
    call self%deallocate
    !AJ Free the communicators associated with this field matrix object.
    if(self%fm_sub_headsc_p0%nproc.gt.0) then
       call free_comm(self%fm_sub_headsc_p0)
    end if
    !Could zero out variables but no real need
  end subroutine fm_reset

  !>Count how many subcommunicators will be created
  subroutine fm_count_subcom(self)
    use mp, only: iproc, barrier, nproc
    implicit none
    class(fieldmat_type), intent(inout) :: self
    integer :: ik, is, ip
    integer :: nsub_tot, nsub_proc

    !Initialise
    nsub_tot=2 !Subcom belonging to fieldmat object
    nsub_proc=nsub_tot

    !Now get children
    do ik=1,self%naky
       nsub_tot=nsub_tot+2
       if(self%kyb(ik)%is_local) nsub_proc=nsub_proc+2

       do is=1,self%kyb(ik)%nsupercell
          nsub_tot=nsub_tot+2
          if(self%kyb(ik)%supercells(is)%is_local) nsub_proc=nsub_proc+2
       enddo
    enddo

    !Now report
    if(iproc.eq.0) write(dlun,'("Total number of subcommunicators is ",I0)') nsub_tot
    do ip=0,nproc-1
       if(ip.eq.iproc) write(dlun,'("Number of subcommunicators on ",I0," is ",I0)') iproc, nsub_proc
       call barrier
    enddo
    if(iproc.eq.0)write(dlun,'()')
  end subroutine fm_count_subcom

  !>Create all the necessary subcommunicators
  subroutine fm_make_subcom_1(self)
    use mp, only: iproc, nproc, mp_comm, comm_type, split, free_comm, proc0
    implicit none
    class(fieldmat_type), intent(inout) :: self
    integer :: ik

    !Set the all comm to be mp_comm
    self%fm_sub_all%id=mp_comm
    self%fm_sub_all%iproc=iproc
    self%fm_sub_all%nproc=nproc
    self%fm_sub_all%proc0=proc0

    !Set the parents and make their subcommunicators
    do ik=1,self%naky
       !Set parent subcom
       self%kyb(ik)%parent_sub=self%fm_sub_all

       !Make individual subcom
       call self%kyb(ik)%make_subcom_1
    enddo

  end subroutine fm_make_subcom_1

  !>Create the secondary subcommunicators
  subroutine fm_make_subcom_2(self)
    use mp, only: split, free_comm
    implicit none
    class(fieldmat_type), intent(inout) :: self
    integer :: ik, colour
    type(comm_type) :: tmp
    logical :: has_head_sc

    !Make child subcommunicators
    do ik=1,self%naky
       !Make individual subcom
       call self%kyb(ik)%make_subcom_2
    enddo

    !Now make the supercell head communicator used for field gathers
    !/First work out if this proc is the head of any supercells
    has_head_sc=.false.
    do ik=1,self%naky
       has_head_sc=has_head_sc.or.any(self%kyb(ik)%supercells%is_head)
    enddo
    !/Now split to make a subcom
    colour=0
    if(has_head_sc.or.self%fm_sub_all%proc0) colour=1
    call split(colour,tmp,self%fm_sub_all%id)
    if(has_head_sc.or.self%fm_sub_all%proc0) then
       !/Store subcom
       self%fm_sub_headsc_p0=tmp
    else
       self%fm_sub_headsc_p0%nproc=0
       !/Destroy subcom
       call free_comm(tmp)
    endif
       
  end subroutine fm_make_subcom_2

  !AJ Record the process that is the supercell head associated with a given ik,it point
  !AJ Supercells have a single ik but may have many it points so a given supercell head 
  !AJ may be responsible for a range of points
  subroutine fm_calc_sc_heads(self)
    use mp, only: iproc, sum_allreduce
    implicit none

    class(fieldmat_type), intent(inout) :: self
    integer :: ik, is, ic, it
    !AJ There may be a problem here is any supercell has no assigned head.
    !AJ This is the case (I think) for 1,1; but it shouldn't actually cause 
    !AJ a problem as 1,1 shouldn't be communicated anyway as only one process
    !AJ is involved in it.
    self%heads = 0
    do ik=1,self%naky
       do is=1,self%kyb(ik)%nsupercell
          if(self%kyb(ik)%supercells(is)%is_head) then
             do ic=1,self%kyb(ik)%supercells(is)%ncell               
                it = self%kyb(ik)%supercells(is)%cells(ic)%it_ind
                self%heads(it,ik) = iproc
             end do
          end if
       end do
    end do

    call sum_allreduce(self%heads)

  end subroutine fm_calc_sc_heads

  !AJ Redistribute the fields from g_lo to gf_lo
  !AJ This is used to take the initial fields that are calculated 
  !AJ in initialisation using the fields_implicit functionality, 
  !AJ and therefore in g_lo, and distribute them to gf_lo layout as 
  !AJ well so that we can start the advance_gf_local routine with the 
  !AJ correct field data in both formats.
  subroutine fm_scatter_fields(self,gf_ph,gf_ap,gf_bp,ph,ap,bp)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, kwork_filter
    use mp, only: nbsend, nbrecv, waitall, broadcast_sub
    use mp, only: rank_comm, comm_type, proc0, iproc, initialise_requests
    use run_parameters, only: fphi, fapar, fbpar
    use gs2_layouts, only: gf_lo, g_lo, proc_id, idx, ik_idx, it_idx
    implicit none
    class(fieldmat_type), intent(inout) :: self
    complex, dimension(-ntgrid:,:,:), intent(inout) :: ph,ap,bp,gf_ph,gf_ap,gf_bp
    integer :: igf, ik, it, sendreqnum, recvreqnum, sender, receiver, tag
    integer, dimension (:), allocatable :: recv_handles, send_handles
    integer :: tempfield
    complex, dimension(:,:,:,:), allocatable :: tempdata

    allocate(tempdata(-ntgrid:ntgrid,nfield,gf_lo%ntheta0,gf_lo%naky))

    !AJ This section is for sending data from one of the g_lo owner to the 
    !AJ process that owns this it,ik point in gf_lo.


    recvreqnum = 1

    allocate(recv_handles(naky*ntheta0))
    allocate(send_handles(naky*ntheta0))

    call initialise_requests(send_handles)
    call initialise_requests(recv_handles)

    if(proc0) then
       
       do ik = 1,self%naky
          do it = 1,self%ntheta0
             
             if(kwork_filter(it,ik)) cycle
             !AJ Calculate which of the processes that own this it,ik point in g_lo space will send it to proc0
             !AJ This can be optimised by ensuring we choose proc0 to send this if proc0 is one of the owners of this
             !AJ in g_lo (rather than just choosing the last process as we do here and elsewhere for load balance).
             sender = g_lo%ikit_procs_assignment(it,ik)%proc_list(g_lo%ikit_procs_assignment(it,ik)%num_procs)       
             
             if(sender .ne. iproc) then
                
                tag = ik*self%ntheta0 + it
                
                call nbrecv(tempdata(:,:,it,ik),sender,tag,recv_handles(recvreqnum))
                recvreqnum = recvreqnum + 1
                
             end if
          end do
       end do

    end if
    
    sendreqnum = 1
    
    !AJ Iterate through all the points I own in gf_lo and send that part of the field 
    !AJ arrays to proc0
    receiver = 0

    do ik = 1,g_lo%naky
       do it = 1,g_lo%ntheta0
             
          if(kwork_filter(it,ik)) cycle          
          
          !AJ This can be optimised by ensuring we choose proc0 to send this if proc0 is one of the owners of this
          !AJ in g_lo (rather than just choosing the last process as we do here and elsewhere for load balance).
          if(g_lo%ikit_procs_assignment(it,ik)%proc_list(g_lo%ikit_procs_assignment(it,ik)%num_procs) .eq. iproc) then
             
             if(iproc .ne. receiver) then
                
                tag = ik*self%ntheta0 + it
                
                tempfield = 1
                if(fphi>epsilon(0.0)) then
                   tempdata(:,tempfield,it,ik) = ph(:,it,ik)
                   tempfield = tempfield + 1
                end if
                
                if(fapar>epsilon(0.0)) then
                   tempdata(:,tempfield,it,ik) = ap(:,it,ik)
                   tempfield = tempfield + 1
                end if
                   
                if(fbpar>epsilon(0.0)) then
                   tempdata(:,tempfield,it,ik) = bp(:,it,ik)
                end if
                
                call nbsend(tempdata(:,:,it,ik),receiver,tag,send_handles(sendreqnum))
                sendreqnum = sendreqnum + 1
                
             end if
             
          end if
          
       end do
    end do
    
    
    call waitall(recvreqnum-1,recv_handles)
    
    !AJ Now copy the received data out of the temporary array to the data arrays on proc0
    if(proc0) then
       
       do ik = 1,self%naky
          do it = 1,self%ntheta0
             
             if(kwork_filter(it,ik)) cycle
             !AJ Calculate which of the processes that own this it,ik point in g_lo space will send it to proc0
             !AJ This can be optimised by ensuring we choose proc0 to send this if proc0 is one of the owners of this
             !AJ in g_lo (rather than just choosing the last process as we do here and elsewhere as a very 
             !AJ basic attempt at load balancing).
             sender = g_lo%ikit_procs_assignment(it,ik)%proc_list(g_lo%ikit_procs_assignment(it,ik)%num_procs)       
                
             if(sender .ne. iproc) then
                
                tempfield = 1
                if(fphi>epsilon(0.0)) then
                   gf_ph(:,it,ik) = tempdata(:,tempfield,it,ik)
                   tempfield = tempfield + 1
                end if
                
                if(fapar>epsilon(0.0)) then
                   gf_ap(:,it,ik) = tempdata(:,tempfield,it,ik)
                   tempfield = tempfield + 1
                end if
                
                if(fbpar>epsilon(0.0)) then
                   gf_bp(:,it,ik) = tempdata(:,tempfield,it,ik)
                end if
                
             else
                
                if(fphi>epsilon(0.0)) then
                   gf_ph(:,it,ik) = ph(:,it,ik)
                end if
                
                if(fapar>epsilon(0.0)) then
                   gf_ap(:,it,ik) = ap(:,it,ik)
                end if
                
                if(fbpar>epsilon(0.0)) then
                   gf_bp(:,it,ik) = bp(:,it,ik)
                end if
                
             end if
          end do
       end do
    end if      

    call waitall(sendreqnum-1,send_handles)

    !AJ Send the data from one of the owners in g_lo to the owner in gf_lo. 
    recvreqnum = 1

    do igf = gf_lo%llim_proc, gf_lo%ulim_proc
       ik = ik_idx(gf_lo,igf)
       it = it_idx(gf_lo,igf)
       
       if(kwork_filter(it,ik)) cycle
       
       tag = ik*self%ntheta0 + it
       
       !AJ Calculate which process of the processes that owns this point in g_lo space will send this data to 
       !AJ the owner in gf_lo
       sender = g_lo%ikit_procs_assignment(it,ik)%proc_list(g_lo%ikit_procs_assignment(it,ik)%num_procs)       
       
       if(sender .ne. iproc) then

          call nbrecv(tempdata(:,:,it,ik),sender,tag,recv_handles(recvreqnum))
          recvreqnum = recvreqnum + 1
          
       end if
    end do

    sendreqnum = 1
    
    do ik = 1,g_lo%naky
       do it = 1,g_lo%ntheta0

          !AJ Choose the sender as the LAST process listed in the process list for a given ik,it point in g_lo
          !AJ We chose the last process here to try and ensure that proc0 isn't overloaded.  Any of the processes 
          !AJ in the proc_list could have been chosen
          if(g_lo%ikit_procs_assignment(it,ik)%proc_list(g_lo%ikit_procs_assignment(it,ik)%num_procs) .eq. iproc) then
             
             if(kwork_filter(it,ik)) cycle

             tag = ik*self%ntheta0 + it

             receiver = proc_id(gf_lo,idx(gf_lo,ik,it))

             if(iproc .ne. receiver) then

                tempfield = 1

                if(fphi>epsilon(0.0)) then
                   tempdata(:,tempfield,it,ik) = ph(:,it,ik)
                   tempfield = tempfield + 1
                end if
          
                if(fapar>epsilon(0.0)) then
                   tempdata(:,tempfield,it,ik) = ap(:,it,ik)
                   tempfield = tempfield + 1
                end if
          
                if(fbpar>epsilon(0.0)) then
                   tempdata(:,tempfield,it,ik) = bp(:,it,ik)
                end if

                call nbsend(tempdata(:,:,it,ik),receiver,tag,send_handles(sendreqnum))          
                sendreqnum = sendreqnum + 1

             end if

          end if

       end do
    end do


    call waitall(recvreqnum-1,recv_handles)

    do igf = gf_lo%llim_proc, gf_lo%ulim_proc
       ik = ik_idx(gf_lo,igf)
       it = it_idx(gf_lo,igf)

       if(kwork_filter(it,ik)) cycle

       sender = g_lo%ikit_procs_assignment(it,ik)%proc_list(g_lo%ikit_procs_assignment(it,ik)%num_procs)       

       if(sender .ne. iproc) then

          tempfield = 1
          if(fphi>epsilon(0.0)) then
             gf_ph(:,it,ik) = tempdata(:,tempfield,it,ik)
             tempfield = tempfield + 1
          end if
          
          if(fapar>epsilon(0.0)) then
             gf_ap(:,it,ik) = tempdata(:,tempfield,it,ik)
             tempfield = tempfield + 1
          end if
          
          
          if(fbpar>epsilon(0.0)) then
             gf_bp(:,it,ik) = tempdata(:,tempfield,it,ik)
          end if

       else

          if(fphi>epsilon(0.0)) then
             gf_ph(:,it,ik) = ph(:,it,ik)
          end if
          
          if(fapar>epsilon(0.0)) then
             gf_ap(:,it,ik) = ap(:,it,ik)
          end if
          
          
          if(fbpar>epsilon(0.0)) then
             gf_bp(:,it,ik) = bp(:,it,ik)
          end if
          
          
       end if

    end do

    call waitall(sendreqnum-1,send_handles)

    deallocate(tempdata)
    deallocate(send_handles, recv_handles)

  end subroutine fm_scatter_fields

  !AJ Redistribute initial fields from gf_lo to g_lo.  These are fields that have not been reduced to the 
  !AJ supercell heads, just calculated on each process that owns gf_lo points
  subroutine fm_gather_init_fields(self,gf_phnew,gf_apnew,gf_bpnew,gf_ph,gf_ap,gf_bp,phnew,apnew,bpnew,ph,ap,bp)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, kwork_filter
    use mp, only: nbsend, nbrecv, waitall, broadcast_sub, sum_reduce_sub, proc0
    use mp, only: rank_comm, comm_type, initialise_requests, iproc
    use run_parameters, only: fphi, fapar, fbpar
    use gs2_layouts, only: gf_lo, g_lo, proc_id, idx, ik_idx, it_idx
    implicit none
    class(fieldmat_type), intent(inout) :: self
    !AJ We have two sets of field arrays; ph,ap,bp store the data we hold in g_lo (i.e. 
    !AJ the field points we have in the g_lo decomposition), and gf_ph,gf_ap,gf_bp hold 
    !AJ the field points we have in the gf_lo decomposition.  This was done for safety (to 
    !AJ ensure that updating fields in gf_lo or g_lo did not interfere with each other), however
    !AJ I THINK it should be possible to do this in a single array.
    complex, dimension(-ntgrid:,:,:), intent(inout) :: ph,ap,bp,gf_ph,gf_ap,gf_bp
    complex, dimension(-ntgrid:,:,:), intent(inout) :: phnew,apnew,bpnew,gf_phnew,gf_apnew,gf_bpnew
    integer :: is, ic, iex, ig, bnd
    integer :: ik, it, ikit, sendreqnum, recvreqnum, sender, receiver, tag
    integer, dimension (:), allocatable :: recv_handles, send_handles
    integer :: tempfield, head, sub, local_iproc
    complex, dimension(:,:,:,:), allocatable :: tempdata

    allocate(tempdata(-ntgrid:ntgrid,nfield,ntheta0,naky))
    !AJ At this point each process has the correct initialised fields in gf_lo (i.e. the points it owns in gf_lo).
    !AJ what we want to do in this routine is distribute them so that all those processes who are part of the supercell
    !AJ in gf_lo have the data, and all the processes that own the points in g_lo also have the data.

    allocate(recv_handles(naky*ntheta0))
    allocate(send_handles(naky*ntheta0))

    call initialise_requests(send_handles)
    call initialise_requests(recv_handles)

    !AJ The first thing we do is send the data from the owner in gf_lo to the supercell head in gf_lo.  This 
    !AJ will enable us to distribute it to the whole supercell in gf_lo, and then after that to the processes 
    !AJ that need it in g_lo

    tempdata = 0.0

    recvreqnum = 1
    sendreqnum = 1
    do ik=1,self%naky
       if(self%kyb(ik)%is_empty) cycle
       do is=1,self%kyb(ik)%nsupercell
          if(self%kyb(ik)%supercells(is)%is_empty) cycle
          receiver = self%kyb(ik)%supercells(is)%head_iproc_world
          do ic=1,self%kyb(ik)%supercells(is)%ncell                
             it=self%kyb(ik)%supercells(is)%cells(ic)%it_ind
             sender = proc_id(gf_lo,idx(gf_lo,ik,it))
             tag = ik*self%ntheta0 + it
             if(sender .eq. iproc .and. receiver .ne. iproc) then
                tempfield = 1
                if(fphi>epsilon(0.0)) then
                   tempdata(:,tempfield,it,ik) = phnew(:,it,ik)
                   !AJ This copy ensures that the field data is also up to date in the gf_lo data structures 
                   !AJ (the implicit field initialise calculates it in the g_lo data structures.  This could be
                   !AJ fixed in fields_implicit init_allfields_implicit making this copy redundant)
                   gf_ph(:,it,ik) = phnew(:,it,ik)
                   gf_phnew(:,it,ik) = phnew(:,it,ik)
                   tempfield = tempfield + 1
                end if
                if(fapar>epsilon(0.0)) then
                   tempdata(:,tempfield,it,ik) = apnew(:,it,ik)
                   !AJ This copy ensures that the field data is also up to date in the gf_lo data structures 
                   !AJ (the implicit field initialise calculates it in the g_lo data structures.  This could be
                   !AJ fixed in fields_implicit init_allfields_implicit making this copy redundant)
                   gf_ap(:,it,ik) = apnew(:,it,ik)
                   gf_apnew(:,it,ik) = apnew(:,it,ik)
                   tempfield = tempfield + 1
                end if
                if(fbpar>epsilon(0.0)) then
                   tempdata(:,tempfield,it,ik) = bpnew(:,it,ik)
                   !AJ This copy ensures that the field data is also up to date in the gf_lo data structures 
                   !AJ (the implicit field initialise calculates it in the g_lo data structures.  This could be
                   !AJ fixed in fields_implicit init_allfields_implicit making this copy redundant)
                   gf_bp(:,it,ik) = bpnew(:,it,ik)
                   gf_bpnew(:,it,ik) = bpnew(:,it,ik)
                end if
                call nbsend(tempdata(:,:,it,ik),receiver,tag,send_handles(sendreqnum))
                sendreqnum = sendreqnum + 1
             else if(receiver .eq. iproc .and. sender .ne. iproc) then
                call nbrecv(tempdata(:,:,it,ik),sender,tag,recv_handles(recvreqnum))                
                recvreqnum = recvreqnum + 1
                !AJ Both the sender and receiver are the same process, so this process 
                !AJ already has the correct gf_lo field data, but it is currently stored 
                !AJ in the g_lo arrays (phi,apar,bpar) instead of the gf_lo arrays 
                !AJ (gf_phi, gf_apar, gf_bpar).  Copy it across so it is in the correct array.
             else if(receiver .eq. iproc .and. sender .eq. iproc) then
                tempfield = 1
                if(fphi>epsilon(0.0)) then
                   gf_phnew(:,it,ik) = phnew(:,it,ik)
                   gf_ph(:,it,ik) = phnew(:,it,ik)
                   !AJ This copy is necessary so the data is in the correct place for the next 
                   !AJ communication (the sending from supercell head to the first proc in the 
                   !AJ g_lo list).
                   tempdata(:,tempfield,it,ik) = phnew(:,it,ik)
                   tempfield = tempfield + 1
                end if
                if(fapar>epsilon(0.0)) then
                   gf_apnew(:,it,ik) = apnew(:,it,ik)
                   gf_ap(:,it,ik) = apnew(:,it,ik)
                   !AJ This copy is necessary so the data is in the correct place for the next 
                   !AJ communication (the sending from supercell head to the first proc in the 
                   !AJ g_lo list).
                   tempdata(:,tempfield,it,ik) = apnew(:,it,ik)
                   tempfield = tempfield + 1
                end if
                if(fbpar>epsilon(0.0)) then
                   gf_bpnew(:,it,ik) = bpnew(:,it,ik)
                   gf_bp(:,it,ik) = bpnew(:,it,ik)
                   !AJ This copy is necessary so the data is in the correct place for the next 
                   !AJ communication (the sending from supercell head to the first proc in the 
                   !AJ g_lo list).
                   tempdata(:,tempfield,it,ik) = bpnew(:,it,ik)
                end if
             end if
          end do
       end do
    end do


    call waitall(recvreqnum-1,recv_handles)
    
    do ik=1,self%naky
       if(self%kyb(ik)%is_empty) cycle
       do is=1,self%kyb(ik)%nsupercell
          if(self%kyb(ik)%supercells(is)%is_empty) cycle
          if(.not.self%kyb(ik)%supercells(is)%is_head) cycle
          do ic=1,self%kyb(ik)%supercells(is)%ncell                
             it=self%kyb(ik)%supercells(is)%cells(ic)%it_ind
             sender = proc_id(gf_lo,idx(gf_lo,ik,it))
             if(sender .ne. iproc) then
                tempfield = 1
                if(fphi>epsilon(0.0)) then
                   gf_phnew(:,it,ik) = tempdata(:,tempfield,it,ik)
                   gf_ph(:,it,ik) = gf_phnew(:,it,ik)
                   tempfield = tempfield + 1
                end if
                if(fapar>epsilon(0.0)) then
                   gf_apnew(:,it,ik) = tempdata(:,tempfield,it,ik)
                   gf_ap(:,it,ik) = gf_apnew(:,it,ik)
                   tempfield = tempfield + 1
                end if
                if(fbpar>epsilon(0.0)) then
                   gf_bpnew(:,it,ik) = tempdata(:,tempfield,it,ik)
                   gf_bp(:,it,ik) = gf_bpnew(:,it,ik)
                end if
             end if
          end do
       end do
    end do
    
    call waitall(sendreqnum-1,send_handles)


    !AJ Now we collect all the data from the supercell heads to proc0 for diagnostics.  If
    !AJ we implemented distributed diagnostics then this would not be necessary.    
    if(self%fm_sub_headsc_p0%nproc.gt.0)then
       call sum_reduce_sub(tempdata,0,self%fm_sub_headsc_p0)
       
       !AJ proc0 unpacks the received field data and updates it's fields from it.
       !AJ Because proc0 now has all the field data it can update both the gf_lo
       !AJ fields and the g_lo fields.
       if(proc0) then
          
          tempfield = 1
          if(fphi>epsilon(0.0)) then
             ph = tempdata(:,tempfield,:,:)
             phnew = ph
             gf_ph = ph
             gf_phnew = ph
             tempfield = tempfield + 1
          end if
             
          if(fapar>epsilon(0.0)) then
             ap = tempdata(:,tempfield,:,:)
             apnew = ap 
             gf_ap = ap
             gf_apnew = ap
             tempfield = tempfield + 1
          end if
          
          if(fbpar>epsilon(0.0)) then
             bp = tempdata(:,tempfield,:,:)
             bpnew = bp
             gf_bp = bp
             gf_bpnew = bp
          end if
          
       end if
       
    end if

       
    call initialise_requests(send_handles)
    call initialise_requests(recv_handles)
    
    !AJ Now send the data from the supercell head in gf_lo to the first process in the list of owners for
    !AJ the same point in g_lo
    recvreqnum = 1
    do ikit=1, g_lo%ikitrange
       ik=g_lo%local_ikit_points(ikit)%ik
       it=g_lo%local_ikit_points(ikit)%it
       if(kwork_filter(it,ik)) cycle
       sender = self%heads(it,ik)
       if(sender .ne. iproc .and. .not. proc0) then
          if(g_lo%ikit_procs_assignment(it,ik)%proc_list(1) .eq. iproc) then
             tag = ik*self%ntheta0 + it          
             call nbrecv(tempdata(:,:,it,ik),sender,tag,recv_handles(recvreqnum))             
             recvreqnum = recvreqnum + 1
          end if
       end if
    end do
    
    sendreqnum = 1
    
    do ik=1,self%naky
       if(self%kyb(ik)%is_empty) cycle
       do is=1,self%kyb(ik)%nsupercell
          if(self%kyb(ik)%supercells(is)%is_empty) cycle
          if(self%kyb(ik)%supercells(is)%is_head) then
             do ic=1,self%kyb(ik)%supercells(is)%ncell                
                it=self%kyb(ik)%supercells(is)%cells(ic)%it_ind
                if(kwork_filter(it,ik)) cycle
                receiver = g_lo%ikit_procs_assignment(it,ik)%proc_list(1) 
                if(receiver .ne. iproc .and. receiver .ne. 0) then
                   tag = ik*self%ntheta0 + it    
                   call nbsend(tempdata(:,:,it,ik),receiver,tag,send_handles(sendreqnum))
                   sendreqnum = sendreqnum + 1
                !AJ If I'm to be sending to myself then just copy the data from gf_lo arrays (where 
                !AJ is was put in the communication to get it to supercell heads) to the correct g_lo arrays.
                else if(receiver .eq. iproc .and. .not. proc0) then
                   if(fphi>epsilon(0.0)) then
                      phnew(:,it,ik) = gf_ph(:,it,ik)
                      ph(:,it,ik) = phnew(:,it,ik)
                   end if
                   if(fapar>epsilon(0.0)) then
                      apnew(:,it,ik) = gf_ap(:,it,ik)
                      ap(:,it,ik) = apnew(:,it,ik)
                   end if
                   if(fbpar>epsilon(0.0)) then
                      bpnew(:,it,ik) = gf_bp(:,it,ik)
                      bp(:,it,ik) = bpnew(:,it,ik)
                   end if
                end if
             end do
          end if
       end do
    end do
    
    call waitall(recvreqnum-1,recv_handles)
    
    do ikit=1, g_lo%ikitrange
       ik=g_lo%local_ikit_points(ikit)%ik
       it=g_lo%local_ikit_points(ikit)%it
       if(kwork_filter(it,ik)) cycle
       sender = self%heads(it,ik)
       receiver = g_lo%ikit_procs_assignment(it,ik)%proc_list(1) 
       if(sender .ne. iproc .and. receiver .eq. iproc .and. .not. proc0) then
          tempfield = 1
          if(fphi>epsilon(0.0)) then
             phnew(:,it,ik) = tempdata(:,tempfield,it,ik)
             ph(:,it,ik) = phnew(:,it,ik)
             tempfield = tempfield + 1
          end if
          if(fapar>epsilon(0.0)) then
             apnew(:,it,ik) = tempdata(:,tempfield,it,ik)
             ap(:,it,ik) = apnew(:,it,ik)
             tempfield = tempfield + 1
          end if
          if(fbpar>epsilon(0.0)) then
             bpnew(:,it,ik) = tempdata(:,tempfield,it,ik)
             bp(:,it,ik) = bpnew(:,it,ik)
          end if
       end if
    end do
    
    call waitall(sendreqnum-1,send_handles)

    deallocate(send_handles, recv_handles)
    
    !AJ Now get the data from the process in g_lo that has received a given ik,it point to the rest of the 
    !AJ processes that own it.
    !AJ This could be optimised if the xyblock_comm is available, but as a first start the code below works
    !AJ and doesn't involve global collectives.
    do ik = 1,self%naky
       do it = 1,self%ntheta0       
          
          if(kwork_filter(it,ik)) cycle

          if(g_lo%ikit_procs_assignment(it,ik)%mine) then

             if(g_lo%ikit_procs_assignment(it,ik)%num_procs .gt. 1) then
                sender = g_lo%ikit_procs_assignment(it,ik)%sub_proc_list(1)
                local_iproc = g_lo%ikit_procs_assignment(it,ik)%comm%iproc
                
                !AJ This could be avoided on some proceses as tempdata will already have the correct data 
                !AJ (providing sender .ne. receiver in the above calls)
                if(sender .eq. local_iproc .and. .not. proc0) then
                   tempfield = 1
                   if(fphi>epsilon(0.0)) then
                      tempdata(:,tempfield,it,ik) = phnew(:,it,ik)
                      tempfield = tempfield + 1
                   end if
                   
                   if(fapar>epsilon(0.0)) then
                      tempdata(:,tempfield,it,ik) = apnew(:,it,ik)
                      tempfield = tempfield + 1
                   end if
                   
                   if(fbpar>epsilon(0.0)) then
                      tempdata(:,tempfield,it,ik) = bpnew(:,it,ik)
                   end if
                   
                end if
                call broadcast_sub(tempdata(:,:,it,ik),sender,g_lo%ikit_procs_assignment(it,ik)%comm%id)
                
                if(sender .ne. local_iproc .and. .not. proc0) then
                   tempfield = 1
                   if(fphi>epsilon(0.0)) then
                      phnew(:,it,ik) = tempdata(:,tempfield,it,ik)
                      ph(:,it,ik) = phnew(:,it,ik)
                      tempfield = tempfield + 1
                   end if
                   
                   if(fapar>epsilon(0.0)) then
                      apnew(:,it,ik) = tempdata(:,tempfield,it,ik)
                      ap(:,it,ik) = apnew(:,it,ik)
                      tempfield = tempfield + 1
                   end if
                   
                   if(fbpar>epsilon(0.0)) then
                      bpnew(:,it,ik) = tempdata(:,tempfield,it,ik)
                      bp(:,it,ik) = bpnew(:,it,ik)
                   end if
                end if
             end if
          end if
       end do
    end do
    
    deallocate(tempdata)

  end subroutine fm_gather_init_fields
  

  !AJ Redistribute fields from gf_lo to g_lo.  This first involves unpacking the calculated 
  !AJ fields from tmpsum on the supercell heads to the field arrays and then sending them back 
  !AJ to the processes that calculated them in gf_lo and the processes that need them in g_lo.
  subroutine fm_gather_fields(self,gf_ph,gf_ap,gf_bp,ph,ap,bp)
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, kwork_filter
    use mp, only: nbsend, nbrecv, waitall, broadcast_sub, sum_reduce_sub
    use mp, only: rank_comm, comm_type, proc0, iproc, initialise_requests
    use run_parameters, only: fphi, fapar, fbpar
    use gs2_layouts, only: gf_lo, g_lo, proc_id, idx, ik_idx, it_idx
    implicit none
    class(fieldmat_type), intent(inout) :: self
    !AJ We have two sets of field arrays; ph,ap,bp store the data we hold in g_lo (i.e. 
    !AJ the field points we have in the g_lo decomposition), and gf_ph,gf_ap,gf_bp hold 
    !AJ the field points we have in the gf_lo decomposition.  This was done for safety (to 
    !AJ ensure that updating fields in gf_lo or g_lo did not interfere with each other), however
    !AJ I THINK it should be possible to do this in a single array.
    complex, dimension(-ntgrid:,:,:), intent(inout) :: ph,ap,bp,gf_ph,gf_ap,gf_bp
    integer :: is, ic, iex, ig, bnd
    integer :: ik, it, ikit, sendreqnum, recvreqnum, sender, receiver, tag
    integer, dimension (:), allocatable :: recv_handles, send_handles
    integer :: tempfield, head, sub, local_iproc
    complex, dimension(:,:,:,:), allocatable :: tempdata,tempdata3
    complex, dimension(:,:,:), allocatable :: tempdata2
    logical :: use_broadcast = .false.
    allocate(tempdata(-ntgrid:ntgrid,nfield,ntheta0,naky))
    allocate(tempdata3(-ntgrid:ntgrid,nfield,ntheta0,naky))
    
    tempdata = 0.0
   
!CMRQ: cost of Supercell Head Packing?
!CMR: BEGIN  Supercell Head Packing
    !AJ Unpack tmp_sum on the supercell head, tmp_sum has the calculated field update on the supercell head
    do ik=1,self%naky
       if(.not.self%kyb(ik)%is_local) cycle
       do is=1,self%kyb(ik)%nsupercell
          !Skip if we're not the head supercell
          if(.not.self%kyb(ik)%supercells(is)%is_head) cycle
          
          iex=0
          tempfield = 1
          if(fphi>epsilon(0.0)) then
             do ic=1,self%kyb(ik)%supercells(is)%ncell
                it=self%kyb(ik)%supercells(is)%cells(ic)%it_ind
                bnd=0
                if(ic.ne.self%kyb(ik)%supercells(is)%ncell) bnd=1
                do ig=-ntgrid,ntgrid-bnd
                   iex=iex+1
                   tempdata(ig,tempfield,it,ik)=self%kyb(ik)%supercells(is)%tmp_sum(iex)
                enddo
             enddo
             tempfield = tempfield + 1
          end if
          
          if(fapar>epsilon(0.0)) then
             do ic=1,self%kyb(ik)%supercells(is)%ncell
                it=self%kyb(ik)%supercells(is)%cells(ic)%it_ind
                bnd=0
                if(ic.ne.self%kyb(ik)%supercells(is)%ncell) bnd=1
                do ig=-ntgrid,ntgrid-bnd
                   iex=iex+1
                   tempdata(ig,tempfield,it,ik)=self%kyb(ik)%supercells(is)%tmp_sum(iex)
                enddo
             enddo
             tempfield = tempfield + 1
          end if
          
          if(fbpar>epsilon(0.0)) then
             do ic=1,self%kyb(ik)%supercells(is)%ncell
                it=self%kyb(ik)%supercells(is)%cells(ic)%it_ind
                bnd=0
                if(ic.ne.self%kyb(ik)%supercells(is)%ncell) bnd=1
                do ig=-ntgrid,ntgrid-bnd
                   iex=iex+1
                   tempdata(ig,tempfield,it,ik)=self%kyb(ik)%supercells(is)%tmp_sum(iex)
                enddo
             enddo
          end if
          
          
          tempfield  = 1
          
          !<DD>TAGGED:Worth looking at improving this bad memory pattern
          !/Now fix boundary points
          if(self%kyb(ik)%supercells(is)%ncell.gt.1) then
             if(fphi>epsilon(0.0)) then
                do ic=1,self%kyb(ik)%supercells(is)%ncell-1
                   tempdata(ntgrid,tempfield,self%kyb(ik)%supercells(is)%cells(ic)%it_ind,ik)=&
                        tempdata(-ntgrid,tempfield,self%kyb(ik)%supercells(is)%cells(ic+1)%it_ind,ik)
                enddo
                tempfield = tempfield + 1
             end if
             if(fapar>epsilon(0.0)) then
                do ic=1,self%kyb(ik)%supercells(is)%ncell-1
                   tempdata(ntgrid,tempfield,self%kyb(ik)%supercells(is)%cells(ic)%it_ind,ik)=&
                        tempdata(-ntgrid,tempfield,self%kyb(ik)%supercells(is)%cells(ic+1)%it_ind,ik)
                enddo
                tempfield = tempfield + 1
             end if
             if(fbpar>epsilon(0.0)) then
                do ic=1,self%kyb(ik)%supercells(is)%ncell-1
                   tempdata(ntgrid,tempfield,self%kyb(ik)%supercells(is)%cells(ic)%it_ind,ik)=&
                        tempdata(-ntgrid,tempfield,self%kyb(ik)%supercells(is)%cells(ic+1)%it_ind,ik)
                enddo
             end if
             
          endif
       end do
    enddo
    
    !AJ At this point tempdata has the calculated field update for a supercell on the supercell head.
!CMR: END   Supercell Head Packing

!CMRQ: cost of Supercell Head => proc0?
!CMR: BEGIN All Supercell Heads => proc0
    allocate(recv_handles(naky*ntheta0))
    allocate(send_handles(naky*ntheta0))

    call initialise_requests(send_handles)
    call initialise_requests(recv_handles)

    recvreqnum = 1

    !AJ Now we collect all the data from the supercell heads to proc0 for diagnostics.  If
    !AJ we implemented distributed diagnostics then this would not be necessary.
    
    if(self%fm_sub_headsc_p0%nproc.gt.0)then
       call sum_reduce_sub(tempdata,0,self%fm_sub_headsc_p0)
       
       !AJ proc0 unpacks the received field update data and updates it's fields from it.
       !AJ Because proc0 now has all the field update data it can update both the gf_lo
       !AJ fields and the g_lo fields.
       if(proc0) then
          
          tempfield = 1
          if(fphi>epsilon(0.0)) then
             ph = tempdata(:,tempfield,:,:)
             gf_ph = ph
             tempfield = tempfield + 1
          end if
             
          if(fapar>epsilon(0.0)) then
             ap = tempdata(:,tempfield,:,:)
             gf_ap = ap
             tempfield = tempfield + 1
          end if
          
          if(fbpar>epsilon(0.0)) then
             bp = tempdata(:,tempfield,:,:)
             gf_bp = bp
          end if
          
       end if
       
    end if
!CMR: END All Supercell Heads => proc0
       
    !AJ Now send the data back to all processes in the supercell
    !AJ This is less efficient than it needs to be.
    !AJ We could also make use of the fact that in gf_lo a cell is only
    !AJ ever owned by a single owner so we could just send specific 
    !AJ it,ik points to the processes that own them with point-to-point 
    !AJ communications (see the else part of this)

!CMR: BEGIN Supercell Heads => cells
!
!CMR: Note AJ Comment above that USE_BROADCAST=T approach (default) is inefficient!
!
!CMR: Observation, probably need MPI rank to participate in SINGLE supercell
!     MPI rank holding multiple cells might be OK (and avoid comms contention) sometimes:
!        -- if all cells in SAME supercell
!        -- if cells in different supercells that need no COMMS (e.g. supercells with 1 cell)
    if(use_broadcast) then
       do ik=1,self%naky
!CMR: Loop over each supercell (labelled ik, self%kyb(ik)%nsupercell), 
          if(self%kyb(ik)%is_empty) cycle
          do is=1,self%kyb(ik)%nsupercell
             if(self%kyb(ik)%supercells(is)%is_empty) cycle
!CMR: if get here, iproc participates 
             head = self%kyb(ik)%supercells(is)%head_iproc
             sub = self%kyb(ik)%supercells(is)%sc_sub_all%id
             allocate(tempdata2(-ntgrid:ntgrid,nfield,self%kyb(ik)%supercells(is)%ncell))
!CMR: head does all packing ready for broadcast to all supercell
             if(self%kyb(ik)%supercells(is)%is_head) then
                do ic=1,self%kyb(ik)%supercells(is)%ncell                
                   it=self%kyb(ik)%supercells(is)%cells(ic)%it_ind
                   tempfield = 1
                   if(fphi>epsilon(0.0)) then
                      tempdata2(:,tempfield,ic) = tempdata(:,tempfield,it,ik)
                      tempfield = tempfield + 1
                   end if
                   if(fapar>epsilon(0.0)) then
                      tempdata2(:,tempfield,ic) = tempdata(:,tempfield,it,ik)
                      tempfield = tempfield + 1
                   end if
                   if(fbpar>epsilon(0.0)) then
                      tempdata2(:,tempfield,ic) = tempdata(:,tempfield,it,ik)
                   end if
                end do
             end if
             call broadcast_sub(tempdata2,head,sub)
             if(.not.self%kyb(ik)%supercells(is)%is_head) then
!CMR: everyone else unpacks result of broadcast
                do ic=1,self%kyb(ik)%supercells(is)%ncell
                   it=self%kyb(ik)%supercells(is)%cells(ic)%it_ind
                   tempfield = 1
                   if(fphi>epsilon(0.0)) then
                      gf_ph(:,it,ik) = tempdata2(:,tempfield,ic)
                      tempfield = tempfield + 1
                   end if
                   if(fapar>epsilon(0.0)) then
                      gf_ap(:,it,ik) = tempdata2(:,tempfield,ic)
                      tempfield = tempfield + 1
                   end if
                   if(fbpar>epsilon(0.0)) then
                      gf_bp(:,it,ik) = tempdata2(:,tempfield,ic)
                   end if
                end do
             end if
             deallocate(tempdata2)
          end do
       end do

       !AJ Don't use the broadcast to do the sending back of data in the supercell,
       !AJ use point to point communications instead.
    else
!CMR: NB this (more efficient?) approach is OFF by default!
!CMRQ:  Is alternative to USE_BROADCAST=F working?  Is it more efficient? Why not default?
       
       call initialise_requests(send_handles)
       call initialise_requests(recv_handles)
       
       recvreqnum = 1
       sendreqnum = 1
       
       do ik=1,self%naky
          if(self%kyb(ik)%is_empty) cycle
          do is=1,self%kyb(ik)%nsupercell
             if(self%kyb(ik)%supercells(is)%is_empty) cycle
             sender = self%kyb(ik)%supercells(is)%head_iproc_world
             do ic=1,self%kyb(ik)%supercells(is)%ncell                
                it=self%kyb(ik)%supercells(is)%cells(ic)%it_ind
                receiver = proc_id(gf_lo,idx(gf_lo,ik,it))
                !AJ For each supercell that I am head of send the data to the process that owns it in gf_lo
                if(self%kyb(ik)%supercells(is)%is_head) then                
                   !AJ If I am the supercell head and I'm not sending this to proc0 (as that already has all the data)
                   !AJ or to myself then send the data.
                   if(receiver .ne. iproc .and. receiver .ne. 0) then
                      tag = ik*self%ntheta0 + it
                      call nbsend(tempdata(:,:,it,ik),receiver,tag,send_handles(sendreqnum))
                      sendreqnum = sendreqnum + 1
                   !AJ Otherwise, if I'm sending to myself and I'm not proc0 then unpack the 
                   !AJ the data directly from my buffer into the gf fields array.
                   else if(receiver .eq. iproc .and. receiver .ne. 0) then
                      tempfield = 1
                      if(fphi>epsilon(0.0)) then
                         gf_ph(:,it,ik) = tempdata(:,tempfield,it,ik)
                         tempfield = tempfield + 1
                      end if
                      if(fapar>epsilon(0.0)) then
                         gf_ap(:,it,ik) = tempdata(:,tempfield,it,ik)
                         tempfield = tempfield + 1
                      end if
                      if(fbpar>epsilon(0.0)) then
                         gf_bp(:,it,ik) = tempdata(:,tempfield,it,ik)
                      end if
                   end if
                else if(iproc .eq. receiver .and. iproc .ne. 0) then
                   tag = ik*self%ntheta0 + it
                   call nbrecv(tempdata3(:,:,it,ik),sender,tag,recv_handles(recvreqnum))
                   recvreqnum = recvreqnum + 1                   
                end if
             end do
          end do
       end do
       
       call waitall(recvreqnum-1,recv_handles)
       
       do ik=1,self%naky
          if(self%kyb(ik)%is_empty) cycle
          do is=1,self%kyb(ik)%nsupercell
             if(self%kyb(ik)%supercells(is)%is_empty) cycle
             if(self%kyb(ik)%supercells(is)%is_head) cycle
             do ic=1,self%kyb(ik)%supercells(is)%ncell                
                it=self%kyb(ik)%supercells(is)%cells(ic)%it_ind
                receiver = proc_id(gf_lo,idx(gf_lo,ik,it))
                if(receiver .eq. iproc .and. receiver .ne. 0) then                
                   tempfield = 1
                   if(fphi>epsilon(0.0)) then
                      gf_ph(:,it,ik) = tempdata3(:,tempfield,it,ik)
                      tempfield = tempfield + 1
                   end if
                   if(fapar>epsilon(0.0)) then
                      gf_ap(:,it,ik) = tempdata3(:,tempfield,it,ik)
                      tempfield = tempfield + 1
                   end if
                   if(fbpar>epsilon(0.0)) then
                      gf_bp(:,it,ik) = tempdata3(:,tempfield,it,ik)
                   end if
                end if
             end do
          end do
       end do
       
       call waitall(sendreqnum-1,send_handles)
       
    end if

!CMR: END Supercell Heads => cells (gf_fields)

!CMR: BEGIN Supercell Heads => FIRST g_lo rank(ik,it)
       
    call initialise_requests(send_handles)
    call initialise_requests(recv_handles)
    
    recvreqnum = 1
    do ikit=1, g_lo%ikitrange
       ik=g_lo%local_ikit_points(ikit)%ik
       it=g_lo%local_ikit_points(ikit)%it
       if(kwork_filter(it,ik)) cycle
       sender = self%heads(it,ik)
       if(sender .ne. iproc .and. iproc .ne. 0) then
          if(g_lo%ikit_procs_assignment(it,ik)%proc_list(1) .eq. iproc) then
!CMR: if iproc is FIRST processor on list that must receive ik,it in g_lo, set up nbrecv
             tag = ik*self%ntheta0 + it          
             call nbrecv(tempdata3(:,:,it,ik),sender,tag,recv_handles(recvreqnum))             
             recvreqnum = recvreqnum + 1
          end if
       end if
    end do
    
    sendreqnum = 1
    
    do ik=1,self%naky
       if(self%kyb(ik)%is_empty) cycle
       do is=1,self%kyb(ik)%nsupercell
          if(self%kyb(ik)%supercells(is)%is_empty) cycle
          if(self%kyb(ik)%supercells(is)%is_head) then
!CMR: if iproc is supercell head for it,ik, proceed to set up nbsends
             do ic=1,self%kyb(ik)%supercells(is)%ncell                
                it=self%kyb(ik)%supercells(is)%cells(ic)%it_ind
                if(kwork_filter(it,ik)) cycle
                receiver = g_lo%ikit_procs_assignment(it,ik)%proc_list(1) 
                if(receiver .ne. iproc .and. receiver .ne. 0) then            
                   tag = ik*self%ntheta0 + it    
!CMR: nbsend from supercell head to where needed in g_lo
                   call nbsend(tempdata(:,:,it,ik),receiver,tag,send_handles(sendreqnum))
                   sendreqnum = sendreqnum + 1
                else if(receiver .eq. iproc .and. receiver .ne. 0) then
!CMR: unpack data supercell head sends to ITSELF in g_lo
                   tempfield = 1
                   if(fphi>epsilon(0.0)) then
                      ph(:,it,ik) = tempdata(:,tempfield,it,ik)
                      tempfield = tempfield + 1
                   end if
                   if(fapar>epsilon(0.0)) then
                      ap(:,it,ik) = tempdata(:,tempfield,it,ik)
                      tempfield = tempfield + 1
                   end if
                   if(fbpar>epsilon(0.0)) then
                      bp(:,it,ik) = tempdata(:,tempfield,it,ik)
                   end if
                end if
             end do
          end if
       end do
    end do
    
    call waitall(recvreqnum-1,recv_handles)
    
    do ikit=1, g_lo%ikitrange
!CMR: all other cores unpack received data in g_lo
       ik=g_lo%local_ikit_points(ikit)%ik
       it=g_lo%local_ikit_points(ikit)%it
       if(kwork_filter(it,ik)) cycle
       sender = self%heads(it,ik)
       receiver = g_lo%ikit_procs_assignment(it,ik)%proc_list(1) 
       if(sender .ne. iproc .and. receiver .eq. iproc .and. .not. proc0) then
          tempfield = 1
          if(fphi>epsilon(0.0)) then
             ph(:,it,ik) = tempdata3(:,tempfield,it,ik)
             tempfield = tempfield + 1
          end if
          if(fapar>epsilon(0.0)) then
             ap(:,it,ik) = tempdata3(:,tempfield,it,ik)
             tempfield = tempfield + 1
          end if
          if(fbpar>epsilon(0.0)) then
             bp(:,it,ik) = tempdata3(:,tempfield,it,ik)
          end if
       end if
    end do
    
    call waitall(sendreqnum-1,send_handles)

    deallocate(send_handles, recv_handles)
    deallocate(tempdata3)
!CMR: END Supercell Heads => FIRST g_lo rank(ik,it)
    
!CMR: BEGIN   FIRST g_lo rank(ik,it) => OTHER g_lo ranks(ik,it)

    !AJ Now get the data from the process in g_lo that has received a given ik,it point to the rest of the 
    !AJ processes that own it.
    !AJ This could be optimised if the xyblock_comm is available, but as a first start the code below works
    !AJ and doesn't involve global collectives.
    do ik = 1,self%naky
       do it = 1,self%ntheta0       
          
          if(kwork_filter(it,ik)) cycle

          if(g_lo%ikit_procs_assignment(it,ik)%mine) then
             
             sender = g_lo%ikit_procs_assignment(it,ik)%sub_proc_list(1)
             local_iproc = g_lo%ikit_procs_assignment(it,ik)%comm%iproc

             if(sender .eq. local_iproc) then
                tempfield = 1
                if(fphi>epsilon(0.0)) then
                   tempdata(:,tempfield,it,ik) = ph(:,it,ik)
                   tempfield = tempfield + 1
                end if
                
                if(fapar>epsilon(0.0)) then
                   tempdata(:,tempfield,it,ik) = ap(:,it,ik)
                   tempfield = tempfield + 1
                end if
                
                if(fbpar>epsilon(0.0)) then
                   tempdata(:,tempfield,it,ik) = bp(:,it,ik)
                end if
                
             end if

             call broadcast_sub(tempdata(:,:,it,ik),sender,g_lo%ikit_procs_assignment(it,ik)%comm%id)
             
             if(sender .ne. local_iproc .and. .not. proc0) then
                tempfield = 1
                if(fphi>epsilon(0.0)) then
                   ph(:,it,ik) = tempdata(:,tempfield,it,ik)
                   tempfield = tempfield + 1
                end if
                
                if(fapar>epsilon(0.0)) then
                   ap(:,it,ik) = tempdata(:,tempfield,it,ik)
                   tempfield = tempfield + 1
                end if
                
                if(fbpar>epsilon(0.0)) then
                   bp(:,it,ik) = tempdata(:,tempfield,it,ik)
                end if
             end if
          end if
       end do
    end do
    
    deallocate(tempdata)
!CMR: END   FIRST g_lo rank(ik,it) => OTHER g_lo ranks(ik,it)
  end subroutine fm_gather_fields

  !>Update the fields using calculated update
  !AJ We update both the g_lo (phi,apar,bpar) and gf_lo (gf_phi,gf_apar,gf_bpar) fields 
  !AJ as these are needed to make sure we can progress the fields in g_lo and gf_lo both
  !AJ times and thereby not have to transfer the full fields for all processes every iteration, 
  !AJ only the required updates.
  !AJ With some refactoring it should be possible to do all the work on a single set of field 
  !AJ arrays as they are the same size regardless of the layout being used.
  subroutine fm_update_fields(self,phi,apar,bpar,gf_phi,gf_apar,gf_bpar)
    use fields_arrays, only: orig_phi=>phi, orig_apar=>apar, orig_bpar=>bpar
    use fields_arrays, only: orig_gf_phi=>gf_phi, orig_gf_apar=>gf_apar, orig_gf_bpar=>gf_bpar
    use run_parameters, only: fphi, fapar, fbpar
    use kt_grids, only: kwork_filter
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, gf_lo, ik_idx, it_idx
    use mp, only: proc0
    implicit none
    class(fieldmat_type), intent(in) ::  self
    complex,dimension(-ntgrid:,:,:),intent(inout) :: phi,apar,bpar
    complex,dimension(-ntgrid:,:,:),intent(inout) :: gf_phi,gf_apar,gf_bpar
    integer :: ik,it,ikit,igf

    !If we're proc0 then we need to do full array (for diagnostics)
    if(proc0) then
       if(fphi>epsilon(0.0)) then
          phi=phi+orig_phi
          gf_phi=gf_phi+orig_gf_phi
       end if
       if(fapar>epsilon(0.0)) then
          apar=apar+orig_apar
          gf_apar=gf_apar+orig_gf_apar
       end if
       if(fbpar>epsilon(0.0)) then
          bpar=bpar+orig_bpar
          gf_bpar=gf_bpar+orig_gf_bpar
       end if
       return
    else

    !AJ Iterate through each distinct ik,it point I have in g_lo and update the associated fields
    do ikit=1,g_lo%ikitrange
       ik = g_lo%local_ikit_points(ikit)%ik 
       it = g_lo%local_ikit_points(ikit)%it
       if(kwork_filter(it,ik)) cycle
       if(fphi>epsilon(0.0)) phi(:,it,ik)=phi(:,it,ik)+orig_phi(:,it,ik)
       if(fapar>epsilon(0.0)) apar(:,it,ik)=apar(:,it,ik)+orig_apar(:,it,ik)
       if(fbpar>epsilon(0.0)) bpar(:,it,ik)=bpar(:,it,ik)+orig_bpar(:,it,ik)       
    end do

    !AJ Iterate through each gf_lo point I have and update the associated fields
    do igf=gf_lo%llim_proc,gf_lo%ulim_proc
       ik = ik_idx(gf_lo,igf)
       it = it_idx(gf_lo,igf)
       if(kwork_filter(it,ik)) cycle       
       if(fphi>epsilon(0.0)) gf_phi(:,it,ik)=gf_phi(:,it,ik)+orig_gf_phi(:,it,ik)
       if(fapar>epsilon(0.0)) gf_apar(:,it,ik)=gf_apar(:,it,ik)+orig_gf_apar(:,it,ik)
       if(fbpar>epsilon(0.0)) gf_bpar(:,it,ik)=gf_bpar(:,it,ik)+orig_gf_bpar(:,it,ik)       
    end do

    end if

  end subroutine fm_update_fields

  !>Update the fields using the new fields
  subroutine fm_update_fields_newstep(self)
    use fields_arrays, only: phi,apar,bpar,phinew,aparnew,bparnew
    use fields_arrays, only: gf_phi,gf_apar,gf_bpar,gf_phinew,gf_aparnew,gf_bparnew
    use run_parameters, only: fphi, fapar, fbpar
    use kt_grids, only: kwork_filter
    use gs2_layouts, only: g_lo, ik_idx, it_idx, gf_lo
    use mp, only: proc0
    implicit none
    class(fieldmat_type), intent(in) ::  self
    integer :: ik,it,ikit,igf

    !If we're proc0 then we need to do full array (for diagnostics)
    if(proc0) then
       if(fphi>epsilon(0.0)) then
          phi=phinew
          gf_phi=gf_phinew
       end if
       if(fapar>epsilon(0.0)) then
          apar=aparnew
          gf_apar=gf_aparnew
       end if
       if(fbpar>epsilon(0.0)) then 
          bpar=bparnew
          gf_bpar=gf_bparnew
       end if
       return
    else

    !AJ Iterate through each distinct ik,it point I have in g_lo and update the associated fields
    do ikit=1,g_lo%ikitrange
       ik = g_lo%local_ikit_points(ikit)%ik 
       it = g_lo%local_ikit_points(ikit)%it
       if(kwork_filter(it,ik)) cycle
       if(fphi>epsilon(0.0)) phi(:,it,ik)=phinew(:,it,ik)
       if(fapar>epsilon(0.0)) apar(:,it,ik)=aparnew(:,it,ik)
       if(fbpar>epsilon(0.0)) bpar(:,it,ik)=bparnew(:,it,ik)
    end do

    !AJ Iterate through each gf_lo point I have and update the associated fields
    do igf=gf_lo%llim_proc,gf_lo%ulim_proc
       ik = ik_idx(gf_lo,igf)
       it = it_idx(gf_lo,igf)
       if(kwork_filter(it,ik)) cycle       
       if(fphi>epsilon(0.0)) gf_phi(:,it,ik)=gf_phinew(:,it,ik)
       if(fapar>epsilon(0.0)) gf_apar(:,it,ik)=gf_aparnew(:,it,ik)
       if(fbpar>epsilon(0.0)) gf_bpar(:,it,ik)=gf_bparnew(:,it,ik)
    end do
    end if

  end subroutine fm_update_fields_newstep

  !>Write some debug data to file
  subroutine fm_write_debug_data(self)
    use file_utils, only: open_output_file, close_output_file
    use mp, only: iproc, nproc_comm, rank_comm
    implicit none
    class(fieldmat_type), intent(in) :: self
    character(len=20) :: ext
    integer :: unit
    integer :: ik,is,ic,np,ip

    !Make the file extension
    write(ext,'(".debug.",I0,".md")') iproc

    !Open output file
    call open_output_file(unit,ext)

    !Now write data
    write(unit,'("Field matrix debug data for proc : ",I0)') iproc
    write(unit,'(80("="))')
    write(unit,'()')

    !/Here we have details of what cells this proc can see and how
    !they're marked
    write(unit,'("Details of object locality and state")')
    write(unit,'(80("-"))')
    write(unit,'()')
    write(unit,'(7(" ",A13," "))') "ik", "is", "ic", "is_local", "is_empty", "is_all_local", "is_head"
    write(unit,'(7(" ",14("-")))')
    do ik=1,self%naky
       do is=1,self%kyb(ik)%nsupercell
          do ic=1,self%kyb(ik)%supercells(is)%ncell
             write(unit,'(3(I14," "),4(L14," "))') ik,is,ic,&
                  self%kyb(ik)%supercells(is)%cells(ic)%is_local,&
                  self%kyb(ik)%supercells(is)%cells(ic)%is_empty,&
                  self%kyb(ik)%supercells(is)%cells(ic)%is_all_local,&
                  self%kyb(ik)%supercells(is)%is_head
          enddo
       enddo
    enddo
    write(unit,'(" ",7(14("-")," "))')
    write(unit,'()')

    !/Here we see what limits the row blocks have for the cells
    !we actually have
    write(unit,'()')
    write(unit,'("Row limits")')
    write(unit,'(80("-"))')
    write(unit,'()')
    do ik=1,self%naky
       write(unit,'("* ik : ",I0)') ik
!       if(.not.self%kyb(ik)%is_local) cycle
       do is=1,self%kyb(ik)%nsupercell
          write(unit,'("    * is : ",I0)') is
!          if(.not.self%kyb(ik)%supercells(is)%is_local) cycle
          write(unit,'("           ","is_head      - ",L1)') self%kyb(ik)%supercells(is)%is_head
          write(unit,'("           ","is_empty     - ",L1)') self%kyb(ik)%supercells(is)%is_empty
          write(unit,'("           ","is_local     - ",L1)') self%kyb(ik)%supercells(is)%is_local
          write(unit,'("           ","is_all_local - ",L1)') self%kyb(ik)%supercells(is)%is_all_local
          write(unit,'()')
          write(unit,'("        ",4(" ",A8," "))') "ic", "rllim", "rulim", "nrow"
          write(unit,'("        ",4(" ",9("-")))')
          do ic=1,self%kyb(ik)%supercells(is)%ncell
             if(size(self%kyb(ik)%supercells(is)%cells(ic)%rb).gt.0)then
                write(unit,'("        ",4(I9," "))') ic,&
                     self%kyb(ik)%supercells(is)%cells(ic)%rb(1)%row_llim,&
                     self%kyb(ik)%supercells(is)%cells(ic)%rb(1)%row_ulim,&
                     self%kyb(ik)%supercells(is)%cells(ic)%rb(1)%nrow
             endif
          enddo
          write(unit,'("        ",4(" ",9("-")))')
       enddo
    enddo
    write(unit,'()')

    !/Subcommunicators
    write(unit,'()')
    write(unit,'("Subcommunicators")')
    write(unit,'(80("-"))')
    write(unit,'()')
    do ik=1,self%naky
       write(unit,'("* ik : ",I0)') ik
       do is=1,self%kyb(ik)%nsupercell
          write(unit,'("    * is : ",I0)') is
          write(unit,'()')
          write(unit,'("        ",6(" ",A7," "))') "name", "nproc", "iproc", "proc0", "REnproc", "REiproc"
          write(unit,'("        ",6(" ",8("-")))')

          if(self%kyb(ik)%supercells(is)%sc_sub_all%nproc.gt.0)then
             call nproc_comm(self%kyb(ik)%supercells(is)%sc_sub_all%id,np)
             call rank_comm(self%kyb(ik)%supercells(is)%sc_sub_all%id,ip)
          else
             np=-1
             ip=-1
          endif
          write(unit,'("        ",(" ",A8),2(I8," "),1(L8," "),2(I8," "))') "all",&
               self%kyb(ik)%supercells(is)%sc_sub_all%nproc,&
               self%kyb(ik)%supercells(is)%sc_sub_all%iproc,&
               self%kyb(ik)%supercells(is)%sc_sub_all%proc0,&
               np,ip

          if(self%kyb(ik)%supercells(is)%sc_sub_pd%nproc.gt.0)then
             call nproc_comm(self%kyb(ik)%supercells(is)%sc_sub_pd%id,np)
             call rank_comm(self%kyb(ik)%supercells(is)%sc_sub_pd%id,ip)
          else
             np=-1
             ip=-1
          endif
          write(unit,'("        ",(" ",A8),2(I8," "),1(L8," "),2(I8," "))') "pd",&
               self%kyb(ik)%supercells(is)%sc_sub_pd%nproc,&
               self%kyb(ik)%supercells(is)%sc_sub_pd%iproc,&
               self%kyb(ik)%supercells(is)%sc_sub_pd%proc0,&
               np,ip

          if(self%kyb(ik)%supercells(is)%sc_sub_not_full%nproc.gt.0)then
             call nproc_comm(self%kyb(ik)%supercells(is)%sc_sub_not_full%id,np)
             call rank_comm(self%kyb(ik)%supercells(is)%sc_sub_not_full%id,ip)
          else
             np=-1
             ip=-1
          endif
          write(unit,'("        ",(" ",A8),2(I8," "),1(L8," "),2(I8," "))') "nfull",&
               self%kyb(ik)%supercells(is)%sc_sub_not_full%nproc,&
               self%kyb(ik)%supercells(is)%sc_sub_not_full%iproc,&
               self%kyb(ik)%supercells(is)%sc_sub_not_full%proc0,&
               np,ip
          write(unit,'("        ",6(" ",8("-")))')
          write(unit,'()')
       enddo
    enddo
    write(unit,'()')

    !Now close output file
    call close_output_file(unit)
  end subroutine fm_write_debug_data

  subroutine fm_getfieldeq_nogath (self,phi, apar, bpar, fieldeq, fieldeqa, fieldeqp)
    use theta_grid, only: ntgrid
    use dist_fn, only: getan_nogath
    use kt_grids, only: naky, ntheta0
    implicit none
    class(fieldmat_type), intent(in) :: self
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (out) ::fieldeq,fieldeqa,fieldeqp
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: antot, antota, antotp

    call getan_nogath(antot, antota, antotp)
    call self%getfieldeq1_nogath (phi, apar, bpar, antot, antota, antotp, &
         fieldeq, fieldeqa, fieldeqp)

  end subroutine fm_getfieldeq_nogath  

  subroutine fm_getfieldeq1_nogath (self,phi, apar, bpar, antot, antota, antotp, &
       fieldeq, fieldeqa, fieldeqp)
    use gs2_layouts, only: gf_lo,  ik_idx, it_idx
    use theta_grid, only: ntgrid, bmag, delthet, jacob
    use kt_grids, only: naky, ntheta0, aky, kperp2, kwork_filter
    use run_parameters, only: fphi, fapar, fbpar
    use run_parameters, only: beta, tite
    use species, only: spec, has_electron_species
    use dist_fn, only: gamtot,gamtot1, gamtot2, gamtot3, fl_avg, gridfac1, apfac, awgt
    use dist_fn, only: adiabatic_option_switch, adiabatic_option_fieldlineavg, calculate_flux_surface_average
    implicit none
    class(fieldmat_type), intent(in) :: self
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: antot, antota, antotp
    complex, dimension (-ntgrid:,:,:), intent (out) ::fieldeq,fieldeqa,fieldeqp

    integer :: igf, ik, it, is, ic

    if (.not. allocated(fl_avg)) allocate (fl_avg(ntheta0, naky))
    if (.not. has_electron_species(spec)) fl_avg = 0.

    if (.not. has_electron_species(spec)) then
       if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
          call calculate_flux_surface_average(fl_avg,antot)
       end if
    end if

    !AJ TODO: I thinkt his can be replaced by a loop over gf_lo points.
    !Now loop over cells and calculate field equation as required
    do ik=1,self%naky
       !Skip empty cells, note this is slightly different to skipping
       !.not.is_local. Skipping empty is probably faster but may be more dangerous
       if(self%kyb(ik)%is_empty) cycle
       do is=1,self%kyb(ik)%nsupercell
          if(self%kyb(ik)%supercells(is)%is_empty) cycle
          do ic=1,self%kyb(ik)%supercells(is)%ncell
             if(self%kyb(ik)%supercells(is)%cells(ic)%is_empty) cycle

             it=self%kyb(ik)%supercells(is)%cells(ic)%it_ind

             if(kwork_filter(it,ik)) cycle

             !Now fill in data
             if(fphi>epsilon(0.0)) then
                fieldeq(:,it,ik)=antot(:,it,ik)-gamtot(:,it,ik)*gridfac1(:,it,ik)*phi(:,it,ik)
                if(fbpar.gt.epsilon(0.0)) fieldeq(:,it,ik)=fieldeq(:,it,ik)+bpar(:,it,ik)*gamtot1(:,it,ik)
                if(.not.has_electron_species(spec)) fieldeq(:,it,ik)= fieldeq(:,it,ik)+fl_avg(it,ik)
             endif

             if(fapar>epsilon(0.0)) then
                fieldeqa(:,it,ik)=antota(:,it,ik)-kperp2(:,it,ik)*gridfac1(:,it,ik)*apar(:,it,ik)
             endif

             if(fbpar>epsilon(0.0))then
                fieldeqp(:,it,ik)=bpar(:,it,ik)*gridfac1(:,it,ik)+&
                     (antotp(:,it,ik)+bpar(:,it,ik)*gamtot2(:,it,ik)+&
                     0.5*phi(:,it,ik)*gamtot1(:,it,ik))*beta*apfac/bmag(:)**2
             endif
          enddo
       enddo
    enddo
  end subroutine fm_getfieldeq1_nogath

!------------------------------------------------------------------------ 

!================
!PROC_DECOMP
!================
  !>Allocate storage space
  subroutine pc_allocate(self)
    use kt_grids, only: naky, ntheta0
    implicit none
    class(pc_type), intent(inout) :: self
    allocate(self%is_local(ntheta0,naky))
    allocate(self%itik_subcom(ntheta0,naky))
    allocate(self%nproc_per_cell(ntheta0,naky))
    allocate(self%nresp_per_cell(ntheta0,naky))
    allocate(self%navail_per_cell(ntheta0,naky))
  end subroutine pc_allocate

  !>Deallocate storage space
  subroutine pc_deallocate(self)
    implicit none
    class(pc_type), intent(inout) :: self
    if(allocated(self%is_local)) deallocate(self%is_local)
    if(allocated(self%itik_subcom)) deallocate(self%itik_subcom)
    if(allocated(self%nproc_per_cell)) deallocate(self%nproc_per_cell)
    if(allocated(self%nresp_per_cell)) deallocate(self%nresp_per_cell)
    if(allocated(self%navail_per_cell)) deallocate(self%navail_per_cell)
  end subroutine pc_deallocate

  !>Debug printing
  subroutine pc_debug_print(self)
    implicit none
    class(pc_type), intent(in) :: self
    write(dlun,'("Proc type debug print. Index nrow=",I0)') &
         self%current_nresp()
  end subroutine pc_debug_print

  !>Return the current number of rows we're responsible for
  function pc_current_nresp(self)
    implicit none
    class(pc_type), intent(in) :: self
    integer :: pc_current_nresp
    pc_current_nresp=sum(self%nresp_per_cell)
  end function pc_current_nresp

  !>Work out which cells are local
  subroutine pc_find_locality(self)
    use gs2_layouts, only: gf_lo, ik_idx, it_idx
    use mp, only: sum_allreduce
    implicit none
    class(pc_type),intent(inout) :: self
    integer :: it,ik,igf

    !Initialise
    self%is_local=0

    do igf=gf_lo%llim_proc,gf_lo%ulim_proc
       ik=ik_idx(gf_lo,igf)
       it=it_idx(gf_lo,igf)
       self%is_local(it,ik)=1
    enddo

    !Now count how many procs have each cell
    self%nproc_per_cell=self%is_local
    call sum_allreduce(self%nproc_per_cell)
  end subroutine pc_find_locality

  !>Work out how many rows are available in each cell
  subroutine pc_count_avail(self)
    use kt_grids, only: ntheta0, naky
    implicit none
    class(pc_type),intent(inout) :: self
    integer :: it,ik,is

    !Initialise
    self%navail_per_cell=0

    !Loop over cells
    do ik=1,naky
       do it=1,ntheta0
          !If cell isn't local then cycle
          if(self%is_local(it,ik).ne.1) cycle

          !Now find the supercell object with this it,ik
          do is=1,fieldmat%kyb(ik)%nsupercell
             !If this supercell has the it then store the 
             !size of the supercell and exit loop
             if(fieldmat%kyb(ik)%supercells(is)%has_it(it)) then
                self%navail_per_cell(it,ik)=fieldmat%kyb(ik)%supercells(is)%nrow
                exit
             endif
          enddo
       enddo
    enddo

    !Set the total number of responsible rows
    self%navail_tot=sum(self%navail_per_cell)
  end subroutine pc_count_avail

  !>Do we have any cells with passed ik?
  !NOTE: This just says if we can see index
  !not that we have any data there
  function pc_has_ik(self,ik)
    implicit none
    class(pc_type), intent(in) :: self
    integer, intent(in) :: ik
    logical :: pc_has_ik
    pc_has_ik=any(self%is_local(:,ik).eq.1)
  end function pc_has_ik

  !>Do we have any cells with passed it?
  !NOTE: This just says if we can see index
  !not that we have any data there
  function pc_has_it(self,it)
    implicit none
    class(pc_type), intent(in) :: self
    integer, intent(in) :: it
    logical :: pc_has_it
    pc_has_it=any(self%is_local(it,:).eq.1)
  end function pc_has_it
 
  !>A wrapper routine to do the decomposition
  subroutine pc_make_decomp(self)
    use mp, only: mp_abort
    implicit none
    class(pc_type), intent(inout) :: self

    call self%decomp
  end subroutine pc_make_decomp

  !>A routine to reset the pc object
  subroutine pc_reset(self)
    implicit none
    class(pc_type), intent(inout) :: self
    
    !deallocate
    call self%deallocate

    !Could zero out variables but no real need
  end subroutine pc_reset

  !A routine to initialise the object
  subroutine pc_init(self)
    implicit none
    class(pc_type), intent(inout) :: self

    !First allocate arrays
    call self%allocate

    !Now find cell locality
    call self%find_locality

    !Count avail data
    call self%count_avail
  end subroutine pc_init

  !---------------------------
  !Decomposition routines
  !---------------------------

  !AJ The process decomposition works on the gf_lo layout, so if you have a 
  !AJ point in gf_lo that is in a supercell you are in that supercell and you 
  !AJ own the whole cell associated with that point, otherwise you don't have anything 
  !AJ in that supercell, i.e. each cell has a single owner.
  subroutine pc_decomp(self)
    use mp, only: split, comm_type, free_comm
    implicit none
    class(pc_type), intent(inout) :: self
    integer :: ik, is, ic, ifq, it, llim, ulim
    integer :: nrow_tmp


    !/Don't need to change pc%is_local as this is already correct

    !Make sure we do the invert locally as mpi based version not implemented!
    pc%force_local_invert=.true.

    !Now we (re)calculate how much data is available
    call self%count_avail

    !Initialise the number of rows responsible
    self%nresp_per_cell=0

    !Now loop over all cells, setting the row_llim of
    !their row blocks
    do ik=1,fieldmat%naky
       do is=1,fieldmat%kyb(ik)%nsupercell
          nrow_tmp=fieldmat%kyb(ik)%supercells(is)%nrow
          do ic=1,fieldmat%kyb(ik)%supercells(is)%ncell
             !Store it index
             it=fieldmat%kyb(ik)%supercells(is)%cells(ic)%it_ind

             !AJ cell locality is restricted to a single process in gf_lo, meaning for each cell there is only 
             !AJ one process that will match this if statement.
             if(fieldmat%kyb(ik)%supercells(is)%cells(ic)%is_local)then
                
                !Now work out limits
                llim=1
                ulim=llim+nrow_tmp-1


                !Now loop over row blocks to set index
                do ifq=1,fieldmat%kyb(ik)%supercells(is)%cells(ic)%nrb
                   fieldmat%kyb(ik)%supercells(is)%cells(ic)%rb(ifq)%row_llim=llim
                   fieldmat%kyb(ik)%supercells(is)%cells(ic)%rb(ifq)%row_ulim=ulim
                   call fieldmat%kyb(ik)%supercells(is)%cells(ic)%rb(ifq)%set_nrow
                enddo
             else
                do ifq=1,fieldmat%kyb(ik)%supercells(is)%cells(ic)%nrb
                   fieldmat%kyb(ik)%supercells(is)%cells(ic)%rb(ifq)%row_llim=0
                   fieldmat%kyb(ik)%supercells(is)%cells(ic)%rb(ifq)%row_ulim=0
                   call fieldmat%kyb(ik)%supercells(is)%cells(ic)%rb(ifq)%set_nrow
                enddo
             endif

             !Record the number of responsible rows, note we assume all rb have same size
             if(size(fieldmat%kyb(ik)%supercells(is)%cells(ic)%rb).gt.0)then
                self%nresp_per_cell(it,ik)=fieldmat%kyb(ik)%supercells(is)%cells(ic)%rb(1)%nrow
             else
                self%nresp_per_cell(it,ik)=0
             endif
          enddo
       enddo
    enddo
  end subroutine pc_decomp

!------------------------------------------------------------------------ 

!//////////////////////////////
!// MODULE LEVEL ROUTINES
!//////////////////////////////

  !>Routine to initialise
  subroutine init_fields_gf_local
    use antenna, only: init_antenna
    use theta_grid, only: init_theta_grid
    use kt_grids, only: init_kt_grids
    use gs2_layouts, only: init_gs2_layouts
    use mp, only: proc0, mp_abort, barrier
    use file_utils, only: error_unit
    implicit none

    !Early exit if possible
    if (initialised) return

    !Exit if local fields not supported (e.g. compiled with pgi)
    if(.not.fields_gf_local_functional()) call mp_abort("field_option='gf local' not supported with this build")

    !Make sure any dependencies are already initialised and
    !initialise the field matrix object.
    if (debug.and.proc0) write(6,*) "init_fields_gf_local: gs2_layouts"
    call init_gs2_layouts
    if (debug.and.proc0) write(6,*) "init_fields_gf_local: theta_grid"
    call init_theta_grid
    if (debug.and.proc0) write(6,*) "init_fields_gf_local: kt_grids"
    call init_kt_grids
    if (debug.and.proc0) write(6,*) "init_fields_gf_local: init_fields_matrixlocal"
    call init_fields_matrixlocal
    if (debug.and.proc0) write(6,*) "init_fields_gf_local: antenna"
    call init_antenna
   

    !Set the initialised state
    initialised = .true.
    reinit = .false.
  end subroutine init_fields_gf_local

  function fields_gf_local_unit_test_init_fields_matrixlocal()
    logical :: fields_gf_local_unit_test_init_fields_matrixlocal

    call init_fields_gf_local

    ! check that we can finish and then init
    call finish_fields_gf_local
 
    call init_fields_gf_local

    fields_gf_local_unit_test_init_fields_matrixlocal = .true.

  end function fields_gf_local_unit_test_init_fields_matrixlocal

  !>Routine use to initialise field matrix data
  subroutine init_fields_matrixlocal
    use mp, only: proc0, barrier, iproc
    implicit none
    real :: ts,te

    !Exit if we've already initialised
    if(initialised) return
    if(.not.reinit)then
       !Use the fieldmat init routine to create fieldmat
       !structure and fill some of the basic variables
       !if(debug)call barrier
       if(proc0.and.debug) then
          write(dlun,'("Initialising fieldmat.")')
          call cpu_time(ts)
       endif
       call fieldmat%init
       !if(debug)call barrier
       if(proc0.and.debug) then 
          call cpu_time(te)
          write(dlun,'("--Done in ",F12.4," units")') te-ts
       endif

       !Now initialise the parallel decomposition object
       !Note: We don't actually decompose until we make some
       !of the primary subcom so we know how many procs are 
       !available to each object etc.
       !if(debug)call barrier
       if(proc0.and.debug) then
          write(dlun,'("Initialising pc.")')
          call cpu_time(ts)
       endif
       call pc%init
       !if(debug)call barrier
       if(proc0.and.debug) then 
          call cpu_time(te)
          write(dlun,'("--Done in ",F12.4," units")') te-ts
       endif

       !Now use the pc object to set field object locality
       !and make the primary subcommunicators
       !if(debug)call barrier
       if(proc0.and.debug) then
          write(dlun,'("Setting initial locality and making primary subcommunicators.")')
          call cpu_time(ts)
       endif
       call fieldmat%set_is_local
       call fieldmat%make_subcom_1
       !if(debug)call barrier
       if(proc0.and.debug) then 
          call cpu_time(te)
          write(dlun,'("--Done in ",F12.4," units")') te-ts
       endif

       !Now we can setup the decomposition, which basically
       !consists of setting the row limit of row block objects
       !if(debug)call barrier
       if(proc0.and.debug) then
          write(dlun,'("Initialising decomposition.")')
          call cpu_time(ts)
       endif
       call pc%make_decomp
       !if(debug)call barrier
       if(proc0.and.debug) then 
          call cpu_time(te)
          write(dlun,'("--Done in ",F12.4," units")') te-ts
       endif

       !Now set fieldmat locality and allocate data arrays
       !if(debug)call barrier
       if(proc0.and.debug) then
          write(dlun,'("Setting locality and allocating.")')
          call cpu_time(ts)
       endif
       call fieldmat%set_is_local
       call fieldmat%allocate
       !if(debug)call barrier
       if(proc0.and.debug) then 
          call cpu_time(te)
          write(dlun,'("--Done in ",F12.4," units")') te-ts
       endif

       !Need to make cell/supercell level sub-communicators
       !at some point before we try to invert.
       !Can make them at any point after locality defined

       !/Note that to split a sub-comm everybody in the original
       !communicator needs to make the call. We can either split
       !all groups off from WORLD or we can take "nested subsets"
       !Note sure if there is any perfomance benefit in the actual
       !comms (i'd guess not) but it should be faster to split smaller
       !groups and can do it in parallel so the nested subset approach
       !may be better.
       !if(debug)call barrier
       if(proc0.and.debug)then
          write(dlun,'("Creating the secondary subcoms.")')
          call cpu_time(ts)
       endif
       call fieldmat%make_subcom_2

       call fieldmat%calc_sc_heads

       !if(debug)call barrier
       if(proc0.and.debug) then 
          call cpu_time(te)
          write(dlun,'("--Done in ",F12.4," units")') te-ts
       endif
    endif
!!!!!!!!!!
!!NOTE: All of the above should be a one off setup cost unless something
!!fundamental changes with the parallel layout (a change in nproc/layout etc.)
!!The remaining two tasks (populate+prepare) must be repeated each time the
!!relevant physics and numerical parameters change. In practice this just means
!!when the time step changes.
!!!!!!!!!!

    !Now fill the fieldmat with data
    if(read_response)then
       if(proc0.and.debug)then
          write(dlun,'("Populating from dump.")')
          call cpu_time(ts)
       endif
       call read_response_from_file_gf_local
       if(proc0.and.debug) then 
          call cpu_time(te)
          write(dlun,'("--Done in ",F12.4," units")') te-ts
       endif
    else       
       !if(debug)call barrier
       if(proc0.and.debug)then
          write(dlun,'("Populating.")')
          call cpu_time(ts)
       endif
       call fieldmat%populate
       !if(debug)call barrier
       if(proc0.and.debug) then 
          call cpu_time(te)
          write(dlun,'("--Done in ",F12.4," units")') te-ts
       endif
       
       !Now prepare the response matrices
       !Note: Currently prepare just means invert
       !but at some point we may wish to LU decompose
       !or other approaches etc.
       !if(debug)call barrier
       if(proc0.and.debug) then
          write(dlun,'("Preparing.")')
          call cpu_time(ts)
       endif
       call fieldmat%prepare
       !if(debug)call barrier
       if(proc0.and.debug) then 
          call cpu_time(te)
          write(dlun,'("--Done in ",F12.4," units")') te-ts
       endif
    endif
        
    !Dump response to file
    if(dump_response)then
       if(proc0.and.debug)then
          write(dlun,'("Dumping to file.")')
          call cpu_time(ts)
       endif
       call dump_response_to_file_gf_local
       if(proc0.and.debug) then 
          call cpu_time(te)
          write(dlun,'("--Done in ",F12.4," units")') te-ts
       endif
    endif
    
    !Now write debug data
!    call fieldmat%write_debug_data
  end subroutine init_fields_matrixlocal

  !>Reset the fields_gf_local module
  subroutine reset_fields_gf_local
    implicit none
    initialised=.false.
    reinit=.true.
  end subroutine reset_fields_gf_local

  !>Finish the fields_gf_local module
  subroutine finish_fields_gf_local
    implicit none
    call pc%reset
    call fieldmat%reset
    reinit = .false.
    initialised=.false.
  end subroutine finish_fields_gf_local

  !>Calculate the update to the fields
  subroutine getfield_gf_local(phi,apar,bpar,gf_phi,gf_apar,gf_bpar,do_update_in)
    use fields_arrays, only: time_field
    use theta_grid, only: ntgrid
    use mp, only: proc0
    use job_manage, only: time_message
    implicit none
    complex, dimension(-ntgrid:,:,:), intent(out) :: phi,apar,bpar,gf_phi,gf_apar,gf_bpar !Note, these are actually phinew,... in typical usage
    logical, optional, intent(in) :: do_update_in
    logical :: do_update

    do_update=.false.
    if(present(do_update_in)) do_update=do_update_in

    if (proc0) call time_message(.false.,time_field,' Field Solver')

    !Use fieldmat routine to calculate the field update
    !AJ phi(new),apar(new),bpar(new) are in gf_lo here
    call fieldmat%get_field_update(gf_phi,gf_apar,gf_bpar)
   
    !NOTE: We currently calculate omega at every time step so we
    !actually need to gather everytime, which is a pain!
    !We also fill in the empties here.
    !AJ before this  phi(new),apar(new),bpar(new) are in gf_lo
    !AJ before this  phi(new),apar(new),bpar(new) are in g_lo
    call fieldmat%gather_fields(gf_phi,gf_apar,gf_bpar,phi,apar,bpar)


    !AJ phi(new),apar(new),bpar(new) are in g_lo
    !Thi s routine updates *new fields using gathered update
    if(do_update) call fieldmat%update_fields(phi,apar,bpar,gf_phi,gf_apar,gf_bpar)

    if (proc0) call time_message(.false.,time_field,' Field Solver')
  end subroutine getfield_gf_local

  !>Initialise the fields from the initial g, just uses the
  !fields_implicit routine
  subroutine init_allfields_gf_local
    use fields_implicit, only: init_allfields_implicit
    use fields_arrays, only: phinew, aparnew, bparnew, gf_phinew, gf_aparnew, gf_bparnew
    use fields_arrays, only: phi, apar, bpar, gf_phi, gf_apar, gf_bpar
    implicit none
    call init_allfields_implicit(.true.)
    !AJ After the call above phinew, aparnew, bparnew have the initial fields, but in gf_lo layout
    !AJ The call below redistributed them so that phinew,aparnew,bparnew,phi,apar,bpar have the correct 
    !AJ fields in g_lo layout and gf_phinew,gf_aparnew,gf_bparnew,gf_phi,gf_apar,gf_bpar have the correct 
    !AJ fields in gf_lo layouts.  Note that the new and original arrays (i.e. phinew and phi) will have the same
    !AJ values at this point
    call fieldmat%gather_init_fields(gf_phinew,gf_aparnew,gf_bparnew,gf_phi,gf_apar,gf_bpar,phinew,aparnew,bparnew,phi,apar,bpar)
  end subroutine init_allfields_gf_local

  !>This routine advances the solution by one full time step
  subroutine advance_gf_local(istep, remove_zonal_flows_switch)
    use run_parameters, only: fphi, fapar, fbpar, reset
    use fields_implicit, only: remove_zonal_flows
    use fields_arrays, only: phi, apar, bpar, phinew
    use fields_arrays, only: aparnew, bparnew, apar_ext
    use fields_arrays, only: gf_phi, gf_apar, gf_bpar
    use fields_arrays, only: gf_phinew, gf_aparnew, gf_bparnew
    use dist_fn, only: timeadv, exb_shear, g_exb, g_exbfac
    use dist_fn, only: collisions_advance
    use dist_fn_arrays, only: g, gnew
    use antenna, only: antenna_amplitudes, no_driver
    use mp, only: iproc
    implicit none
    integer, intent(in) :: istep
    logical, intent(in) :: remove_zonal_flows_switch
    integer :: diagnostics = 1
    logical :: do_update

    do_update=.true.

    !GGH NOTE: apar_ext is initialized in this call
    if(.not.no_driver) then
       !AJ THIS MAY NOT WORK FOR GF_LO
       write(*,*) 'antenna has not been configured for gf_lo, may produce incorrect results'
       call antenna_amplitudes (apar_ext)
    end if

    !Apply flow shear if active
    if(abs(g_exb*g_exbfac).gt.epsilon(0.)) then 
       !AJ THIS MAY NOT WORK FOR GF_LO
       write(*,*) 'Flow shear has not been configured for gf_lo, may produce incorrect results'
       call exb_shear(gnew,phinew,aparnew,bparnew,istep,field_local_allreduce_sub)
    end if

    !Update g and fields
    !AJ IS THIS FULL COPY NECESSARY?
    g=gnew
    if(do_update)then !Here we tie the "smart" update to do_update
       call fieldmat%update_fields_newstep
    else
       if(fphi.gt.0) then
          phi=phinew
          gf_phi=gf_phinew
       end if
       if(fapar.gt.0) then
          apar=aparnew
          gf_apar=gf_aparnew
       end if
       if(fbpar.gt.0) then
          bpar=bparnew
          gf_bpar=gf_bparnew
       end if
    endif

    !Find gnew given fields at time step midpoint
    call timeadv (phi, apar, bpar, phinew, &
         aparnew, bparnew, istep)
    if(reset) return !Return is resetting

    !Add in antenna driving if present
    !<DD>TAGGED: Should we only this is fapar>0 as well?
    if(.not.no_driver) then
       !AJ THIS MAY NOT WORK FOR GF_LO
       write(*,*) 'antenna has not been configured for gf_lo, may produce incorrect results'
       aparnew=aparnew+apar_ext
    end if

    !Calculate the fields at the next time point
    call getfield_gf_local(phinew,aparnew,bparnew,gf_phinew,gf_aparnew,gf_bparnew,do_update)

    !If we do the update in getfield_gf_local don't do it here
    if(.not.do_update)then
       if(fphi.gt.0) then 
          phinew=phinew+phi
          gf_phinew=gf_phinew+gf_phi
       end if
       if(fapar.gt.0) then
          aparnew=aparnew+apar
          gf_aparnew=gf_aparnew+gf_apar
       end if
       if(fbpar.gt.0) then 
          bparnew=bparnew+bpar
          gf_bparnew=gf_bparnew+gf_bpar
       end if
    endif    

    !Remove zonal component if requested
    if(remove_zonal_flows_switch) then
       !AJ THIS MAY NOT WORK IN GF_LO
       write(*,*) 'zonal flows have not been configured for gf_lo, may produce incorrect results'
       call remove_zonal_flows
    end if

    !Evolve gnew to next half time point
    call timeadv (phi, apar, bpar, phinew, &
         aparnew, bparnew, istep, diagnostics) 

    ! Advance collisions, if separate from timeadv
    call collisions_advance (phi, bpar, phinew, aparnew, bparnew, istep, diagnostics)

  end subroutine advance_gf_local

  !> Routine to dump the current response matrix data
  !! to file. One file per connected domain. Each written
  !! by the head of the supercell.
  subroutine dump_response_to_file_gf_local(suffix)
    use gs2_save, only: gs2_save_response
    use fields_arrays, only: response_file
    implicit none
    character(len=*), optional, intent(in) :: suffix !If passed then use as part of file suffix
    character(len=64) :: suffix_local, suffix_default='.response'
    character(len=256) :: file_name
    complex, dimension(:,:), allocatable :: tmp_arr
    integer :: ik, is

    !Set file suffix
    suffix_local=suffix_default
    if(present(suffix)) suffix_local=suffix

    !First we have to pull data so that head of each supercell has full data
    !Currently we are lazy and just ensure all procs in supercell has full data
    do ik=1,fieldmat%naky
       !If ik isn't local than cycle
       if(.not.fieldmat%kyb(ik)%is_local) cycle

       !Loop over supercells
       do is=1,fieldmat%kyb(ik)%nsupercell
          !If is isn't local than cycle
          if(.not.fieldmat%kyb(ik)%supercells(is)%is_local) cycle

          !Now allocate array
          allocate(tmp_arr(fieldmat%kyb(ik)%supercells(is)%nrow,fieldmat%kyb(ik)%supercells(is)%ncol))

          !First have to pull all row level data
          call fieldmat%kyb(ik)%supercells(is)%pull_rows_to_arr(tmp_arr)
          
          !Now save if on head proc
          if(fieldmat%kyb(ik)%supercells(is)%is_head)then
             !Really we should do this with netcdf but for now we will
             !simply dump a binary file.

             !First make file name
             write(file_name,'(A,"_ik_",I0,"_is_",I0,A)') trim(response_file),ik,is,trim(suffix_local)
             call gs2_save_response(tmp_arr,file_name)
          endif

          !Deallocate
          deallocate(tmp_arr)
       enddo
    enddo

  end subroutine dump_response_to_file_gf_local

  !> Routine to read the current response matrix data
  !! from file dump. One file per connected domain. Each written
  !! by the head of the supercell.
  !! NOTE: Have to have setup communicators etc.
  subroutine read_response_from_file_gf_local(suffix)
    use gs2_save, only: gs2_restore_response
    use mp, only: sum_allreduce_sub
    use fields_arrays, only: response_file
    implicit none
    character(len=*), optional, intent(in) :: suffix !If passed then use as part of file suffix
    character(len=64) :: suffix_local, suffix_default='.response'
    character(len=256) :: file_name
    complex, dimension(:,:), allocatable :: tmp_arr
    integer :: ik, is

    !Set file suffix
    suffix_local=suffix_default
    if(present(suffix)) suffix_local=suffix

    !First we have to pull data so that head of each supercell has full data
    !Currently we are lazy and just ensure all procs in supercell has full data
    do ik=1,fieldmat%naky
       !If ik isn't local than cycle
       if(.not.fieldmat%kyb(ik)%is_local) cycle

       !Loop over supercells
       do is=1,fieldmat%kyb(ik)%nsupercell
          !If is isn't local than cycle
          if(.not.fieldmat%kyb(ik)%supercells(is)%is_local) cycle

          !Now allocate array
          allocate(tmp_arr(fieldmat%kyb(ik)%supercells(is)%nrow,fieldmat%kyb(ik)%supercells(is)%ncol))
          tmp_arr=0.

          !Now if on head proc read
          if(fieldmat%kyb(ik)%supercells(is)%is_head)then
             !Really we should do this with netcdf but for now we will
             !simply dump a binary file.

             !First make file name
             write(file_name,'(A,"_ik_",I0,"_is_",I0,A)') trim(response_file),ik,is,trim(suffix_local)
             call gs2_restore_response(tmp_arr,file_name)
          endif

          !Now we need to broadcast the data to everyone with supercell.
          !NOTE: Really we would want a supercell bound method which will do
          !this so that this module level routine is insensitive to communicators etc.
          !Furthermore we are lazy and use a sum_allreduce_sub to send data to everyone
          call sum_allreduce_sub(tmp_arr,fieldmat%kyb(ik)%supercells(is)%sc_sub_all%id)

          !Now push array to local storage
          call fieldmat%kyb(ik)%supercells(is)%push_arr_to_rows(tmp_arr)

          !Deallocate
          deallocate(tmp_arr)
       enddo
    enddo

  end subroutine read_response_from_file_gf_local

  !> Returns true if GS2 was built in such a way
  !! as to allow this module to work.
  !! Currently does not work with PGI compilers. 
  !! See online discussions.
  function fields_gf_local_functional()
    use runtime_tests, only: compiler_pgi
    logical :: fields_gf_local_functional
    fields_gf_local_functional = (.not. compiler_pgi())
  end function fields_gf_local_functional

end module fields_gf_local
