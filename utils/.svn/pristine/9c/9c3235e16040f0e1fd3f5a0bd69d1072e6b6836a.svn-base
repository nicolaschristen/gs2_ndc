!>This module provides various routines for sorting arrays
!This isn't something we commonly do but these routines are
!useful in a few places.
!
!As with most algorithms different options can be made depending
!on if we want to optimise memory usage, speed etc.
module sorting
  implicit none

  private

  public :: quicksort, insertsort

!Interfaces
  interface quicksort
     module procedure i_quicksort
     module procedure i_1_quicksort
     module procedure i_2_quicksort
     module procedure i_3_quicksort
  end interface

  interface insertsort
     module procedure i_insertsort
     module procedure i_1_insertsort
     module procedure i_2_insertsort
     module procedure i_3_insertsort
  end interface insertsort

  interface swap_elem
     module procedure i_swap_elem
  end interface swap_elem

contains
!////////////////////
!// QUICKSORT
!////////////////////
  !NOTE: Currently we provide several routines for
  !sorting different numbers of arrays based on the same key.
  !It is likely to be faster to sort an "index" array (e.g.
  !an array of integers initially (/1..n/)) and use this to
  !address as many arrays (of any type) that we want to sort.
  !This will however have some additional memory overhead so
  !for now do things this way.

  !NOTE: Simple tests suggest that passing arr_len explicitly
  !can lead to a substantial speedup of the routine compared to
  !deferred shape.

  !Sort an integer array based on integer KEY
  subroutine i_quicksort (arr_len,key)
    implicit none
    !Size of arrays
    integer, intent(in) :: arr_len

    !Arrays to be sorted
    integer, intent(inout), dimension(1:arr_len) :: key
    integer, parameter :: insertion_limit=16 !Lists smaller than this use insertion sort

    !Internals
    integer, dimension(:), allocatable :: stack_inds
    integer :: i, j, right_lim, left_lim, stack_count
    integer :: ak !Used to hold selected elements
    integer :: k !Temp used for middle of range
    logical :: partitioned!, done_exit

    !Intialise vars
    allocate(stack_inds(max(arr_len/2,200)))
    stack_count=0 !How many sections/partitions
    left_lim=1 !Left index of current section
    right_lim=arr_len !Right index of current section
    
    !Infinite loop --> explicit return forced when list sorted
    do while(.true.)
       if(right_lim-left_lim.lt.insertion_limit)then
          !Insertion sort : Loop over current section
          call insertsort(1+right_lim-left_lim,key(left_lim:right_lim))
          !OR can do insert sort internally
          ! do j=left_lim+1,right_lim
          !    !Get current value
          !    ak=key(j)

          !    done_exit=.false.
          !    do i=j-1,left_lim,-1
          !       !Escape loop
          !       if(key(i).le.ak)then
          !          done_exit=.true.
          !          exit
          !       endif
                
          !       !Shift values along one
          !       key(i+1)=key(i)
          !    enddo
          !    if(.not.done_exit) i=left_lim-1
             
          !    !Insert value in correct place in section
          !    key(i+1)=ak
          ! enddo

          !Check if work remaining

          !If stack is empty then done sorting so return
          if(stack_count.eq.0)then
             deallocate(stack_inds)
             return
          endif
          
          !If not get next section limits from stack
          right_lim=stack_inds(stack_count)
          left_lim=stack_inds(stack_count-1)

          !Reduce stack count
          stack_count=stack_count-2
       else
          !What index are we looking at
          k=(left_lim+right_lim)/2

          !Swap values | KEY
          call swap_elem(key,k,left_lim+1)

          !Order so key(left_lim)<=key(left_lim+1)<=key(right_lim)
          if(key(left_lim).gt.key(right_lim)) then
             call swap_elem(key,left_lim,right_lim)
          endif
          if(key(left_lim+1).gt.key(right_lim)) then
             call swap_elem(key,left_lim+1,right_lim)
          endif
          if(key(left_lim).gt.key(left_lim+1)) then
             call swap_elem(key,left_lim,left_lim+1)
          endif

          !Index of left limit of partition
          i=left_lim+1
          !Index of right limit of partition
          j=right_lim

          !Get partition values
          ak=key(left_lim+1)

          !Partition elements
          partitioned=.false.
          do while(.not.partitioned)
             !Search from left for first element bigger than ak
             i=i+1
             do while(key(i).lt.ak)
                i=i+1
             end do

             !Search from right for first element smaller than ak
             j=j-1
             do while(key(j).gt.ak)
                j=j-1
             end do

             !If first smaller element before first larger element then
             !partioning has finished --> Can now sort
             !Otherwise swap elements and try again
             if(j.lt.i)then
                !Indicate partioning complete
                partitioned=.true.

                !Insert partioning element in correct place
                key(left_lim+1)=key(j) ; key(j)=ak
             else
                call swap_elem(key,i,j)
             endif
          end do

          !Update stack counter
          stack_count=stack_count+2

          !At this point i is the start of the left section, which goes up to right_lim
          !and the right section goes from left_lim up to j
          !Update stack_inds values to record section done
          !If left section is larger push to stack and do right
          !otherwise push right to stack and do left
          if(right_lim-i+1.ge.j-left_lim) then
             !Store detail of left section
             stack_inds(stack_count)=right_lim
             stack_inds(stack_count-1)=i

             !Set right lim to one before start of left section
             right_lim=j-1
          else
             !Store details of right section
             stack_inds(stack_count)=j-1
             stack_inds(stack_count-1)=left_lim

             !Set left lim to one after end of right stack
             left_lim=i
          endif
       end if
    enddo
  end subroutine i_quicksort

  !Sort an integer array based on integer KEY
  subroutine i_1_quicksort (arr_len,key,arr1)
    implicit none
    !Size of arrays
    integer, intent(in) :: arr_len

    !Arrays to be sorted
    integer, intent(inout), dimension(1:arr_len) :: key, arr1
    integer, parameter :: insertion_limit=16 !Lists smaller than this use insertion sort

    !Internals
    integer, dimension(:), allocatable :: stack_inds
    integer :: i, j, right_lim, left_lim, stack_count
    integer :: ak, a1 !Used to hold selected elements
    integer :: k !Temp used for middle of range
    logical :: partitioned!, done_exit

    !Intialise vars
    allocate(stack_inds(max(arr_len/2,200)))
    stack_count=0 !How many sections/partitions
    left_lim=1 !Left index of current section
    right_lim=arr_len !Right index of current section
    
    !Infinite loop --> explicit return forced when list sorted
    do while(.true.)
       if(right_lim-left_lim.lt.insertion_limit)then
          !Insertion sort : Loop over current section
          call insertsort(1+right_lim-left_lim,key(left_lim:right_lim),&
               arr1(left_lim:right_lim))
          !OR can do insert sort internally
          ! do j=left_lim+1,right_lim
          !    !Get current value
          !    ak=key(j) ; a1=arr1(j)

          !    done_exit=.false.
          !    do i=j-1,left_lim,-1
          !       !Escape loop
          !       if(key(i).le.ak)then
          !          done_exit=.true.
          !          exit
          !       endif
                
          !       !Shift values along one
          !       key(i+1)=key(i) ; arr1(i+1)=arr1(i)
          !    enddo
          !    if(.not.done_exit) i=left_lim-1
             
          !    !Insert value in correct place in section
          !    key(i+1)=ak ; arr1(i+1)=a1
          ! enddo

          !Check if work remaining

          !If stack is empty then done sorting so return
          if(stack_count.eq.0)then
             deallocate(stack_inds)
             return
          endif
          
          !If not get next section limits from stack
          right_lim=stack_inds(stack_count)
          left_lim=stack_inds(stack_count-1)

          !Reduce stack count
          stack_count=stack_count-2
       else
          !What index are we looking at
          k=(left_lim+right_lim)/2

          !Swap values | KEY
          call swap_elem(key,k,left_lim+1)
          !Swap values | ARR1
          call swap_elem(arr1,k,left_lim+1)

          !Order so key(left_lim)<=key(left_lim+1)<=key(right_lim)
          if(key(left_lim).gt.key(right_lim)) then
             call swap_elem(key,left_lim,right_lim)
             call swap_elem(arr1,left_lim,right_lim)
          endif
          if(key(left_lim+1).gt.key(right_lim)) then
             call swap_elem(key,left_lim+1,right_lim)
             call swap_elem(arr1,left_lim+1,right_lim)
          endif
          if(key(left_lim).gt.key(left_lim+1)) then
             call swap_elem(key,left_lim,left_lim+1)
             call swap_elem(arr1,left_lim,left_lim+1)
          endif

          !Index of left limit of partition
          i=left_lim+1
          !Index of right limit of partition
          j=right_lim

          !Get partition values
          ak=key(left_lim+1)
          a1=arr1(left_lim+1)

          !Partition elements
          partitioned=.false.
          do while(.not.partitioned)
             !Search from left for first element bigger than ak
             i=i+1
             do while(key(i).lt.ak)
                i=i+1
             end do

             !Search from right for first element smaller than ak
             j=j-1
             do while(key(j).gt.ak)
                j=j-1
             end do

             !If first smaller element before first larger element then
             !partioning has finished --> Can now sort
             !Otherwise swap elements and try again
             if(j.lt.i)then
                !Indicate partioning complete
                partitioned=.true.

                !Insert partioning element in correct place
                key(left_lim+1)=key(j) ; key(j)=ak
                arr1(left_lim+1)=arr1(j) ; arr1(j)=a1
             else
                call swap_elem(key,i,j)
                call swap_elem(arr1,i,j)
             endif
          end do

          !Update stack counter
          stack_count=stack_count+2

          !At this point i is the start of the left section, which goes up to right_lim
          !and the right section goes from left_lim up to j
          !Update stack_inds values to record section done
          !If left section is larger push to stack and do right
          !otherwise push right to stack and do left
          if(right_lim-i+1.ge.j-left_lim) then
             !Store detail of left section
             stack_inds(stack_count)=right_lim
             stack_inds(stack_count-1)=i

             !Set right lim to one before start of left section
             right_lim=j-1
          else
             !Store details of right section
             stack_inds(stack_count)=j-1
             stack_inds(stack_count-1)=left_lim

             !Set left lim to one after end of right stack
             left_lim=i
          endif
       end if
    enddo
  end subroutine i_1_quicksort

  !Sort two integer arrays based on integer KEY
  subroutine i_2_quicksort (arr_len,key,arr1,arr2)
    implicit none
    !Size of arrays
    integer, intent(in) :: arr_len

    !Arrays to be sorted
    integer, intent(inout), dimension(1:arr_len) :: key, arr1, arr2 
    integer, parameter :: insertion_limit=16 !Lists smaller than this use insertion sort

    !Internals
    integer, dimension(:), allocatable :: stack_inds
    integer :: i, j, right_lim, left_lim, stack_count
    integer :: ak, a1, a2 !Used to hold selected elements
    integer :: k !Temp used for middle of range
    logical :: partitioned!, done_exit

    !Intialise vars
    allocate(stack_inds(max(arr_len/2,200)))
    stack_count=0 !How many sections/partitions
    left_lim=1 !Left index of current section
    right_lim=arr_len !Right index of current section
    
    !Infinite loop --> explicit return forced when list sorted
    do while(.true.)
       if(right_lim-left_lim.lt.insertion_limit)then
          !Insertion sort : Loop over current section
          
          !Can call insertsort as follows but could be slower (due to creation
          !of array temporaries).
          call insertsort(1+right_lim-left_lim,key(left_lim:right_lim),&
               arr1(left_lim:right_lim),arr2(left_lim:right_lim))
          !OR can do insert sort internally
          ! do j=left_lim+1,right_lim
          !    !Get current value
          !    ak=key(j) ; a1=arr1(j) ; a2=arr2(j)

          !    done_exit=.false.
          !    do i=j-1,left_lim,-1
          !       !Escape loop
          !       if(key(i).le.ak)then
          !          done_exit=.true.
          !          exit
          !       endif
                
          !       !Shift values along one
          !       key(i+1)=key(i) ; arr1(i+1)=arr1(i) ; arr2(i+1)=arr2(i)
          !    enddo
          !    if(.not.done_exit) i=left_lim-1
             
          !    !Insert value in correct place in section
          !    key(i+1)=ak ; arr1(i+1)=a1 ; arr2(i+1)=a2
          ! enddo

          !Check if work remaining

          !If stack is empty then done sorting so return
          if(stack_count.eq.0)then
             deallocate(stack_inds)
             return
          endif
          
          !If not get next section limits from stack
          right_lim=stack_inds(stack_count)
          left_lim=stack_inds(stack_count-1)

          !Reduce stack count
          stack_count=stack_count-2
       else
          !What index are we looking at
          k=(left_lim+right_lim)/2

          !Swap values | KEY
          call swap_elem(key,k,left_lim+1)
          !Swap values | ARR1
          call swap_elem(arr1,k,left_lim+1)
          !Swap values | ARR2
          call swap_elem(arr2,k,left_lim+1)

          !Order so key(left_lim)<=key(left_lim+1)<=key(right_lim)
          if(key(left_lim).gt.key(right_lim)) then
             call swap_elem(key,left_lim,right_lim)
             call swap_elem(arr1,left_lim,right_lim)
             call swap_elem(arr2,left_lim,right_lim)
          endif
          if(key(left_lim+1).gt.key(right_lim)) then
             call swap_elem(key,left_lim+1,right_lim)
             call swap_elem(arr1,left_lim+1,right_lim)
             call swap_elem(arr2,left_lim+1,right_lim)
          endif
          if(key(left_lim).gt.key(left_lim+1)) then
             call swap_elem(key,left_lim,left_lim+1)
             call swap_elem(arr1,left_lim,left_lim+1)
             call swap_elem(arr2,left_lim,left_lim+1)
          endif

          !Index of left limit of partition
          i=left_lim+1
          !Index of right limit of partition
          j=right_lim

          !Get partition values
          ak=key(left_lim+1)
          a1=arr1(left_lim+1)
          a2=arr2(left_lim+1)

          !Partition elements
          partitioned=.false.
          do while(.not.partitioned)
             !Search from left for first element bigger than ak
             i=i+1
             do while(key(i).lt.ak)
                i=i+1
             end do

             !Search from right for first element smaller than ak
             j=j-1
             do while(key(j).gt.ak)
                j=j-1
             end do

             !If first smaller element before first larger element then
             !partioning has finished --> Can now sort
             !Otherwise swap elements and try again
             if(j.lt.i)then
                !Indicate partioning complete
                partitioned=.true.

                !Insert partioning element in correct place
                key(left_lim+1)=key(j) ; key(j)=ak
                arr1(left_lim+1)=arr1(j) ; arr1(j)=a1
                arr2(left_lim+1)=arr2(j) ; arr2(j)=a2
             else
                call swap_elem(key,i,j)
                call swap_elem(arr1,i,j)
                call swap_elem(arr2,i,j)
             endif
          end do

          !Update stack counter
          stack_count=stack_count+2

          !At this point i is the start of the left section, which goes up to right_lim
          !and the right section goes from left_lim up to j
          !Update stack_inds values to record section done
          !If left section is larger push to stack and do right
          !otherwise push right to stack and do left
          if(right_lim-i+1.ge.j-left_lim) then
             !Store detail of left section
             stack_inds(stack_count)=right_lim
             stack_inds(stack_count-1)=i

             !Set right lim to one before start of left section
             right_lim=j-1
          else
             !Store details of right section
             stack_inds(stack_count)=j-1
             stack_inds(stack_count-1)=left_lim

             !Set left lim to one after end of right stack
             left_lim=i
          endif
       end if
    enddo
  end subroutine i_2_quicksort


  !Sort three integer arrays based on integer KEY
  subroutine i_3_quicksort (arr_len,key,arr1,arr2,arr3)
    implicit none
    !Size of arrays
    integer, intent(in) :: arr_len

    !Arrays to be sorted
    integer, intent(inout), dimension(1:arr_len) :: key, arr1, arr2, arr3
    integer, parameter :: insertion_limit=16 !Lists smaller than this use insertion sort

    !Internals
    integer, dimension(:), allocatable :: stack_inds
    integer :: i, j, right_lim, left_lim, stack_count
    integer :: ak, a1, a2, a3 !Used to hold selected elements
    integer :: k !Temp used for middle of range
    logical :: partitioned!, done_exit

    !Intialise vars
    allocate(stack_inds(max(arr_len/2,200)))
    stack_count=0 !How many sections/partitions
    left_lim=1 !Left index of current section
    right_lim=arr_len !Right index of current section
    
    !Infinite loop --> explicit return forced when list sorted
    do while(.true.)
       if(right_lim-left_lim.lt.insertion_limit)then
          !Insertion sort : Loop over current section
          call insertsort(1+right_lim-left_lim,key(left_lim:right_lim),&
               arr1(left_lim:right_lim),arr2(left_lim:right_lim),arr3(left_lim:right_lim))
          !OR can do insert sort internally
          ! do j=left_lim+1,right_lim
          !    !Get current value
          !    ak=key(j) ; a1=arr1(j) ; a2=arr2(j) ; a3=arr3(j)

          !    done_exit=.false.
          !    do i=j-1,left_lim,-1
          !       !Escape loop
          !       if(key(i).le.ak)then
          !          done_exit=.true.
          !          exit
          !       endif
                
          !       !Shift values along one
          !       key(i+1)=key(i) ; arr1(i+1)=arr1(i) ; arr2(i+1)=arr2(i) ; arr3(i+1)=arr3(i)
          !    enddo
          !    if(.not.done_exit) i=left_lim-1
             
          !    !Insert value in correct place in section
          !    key(i+1)=ak ; arr1(i+1)=a1 ; arr2(i+1)=a2 ; arr3(i+1)=a3
          ! enddo

          !Check if work remaining

          !If stack is empty then done sorting so return
          if(stack_count.eq.0)then
             deallocate(stack_inds)
             return
          endif
          
          !If not get next section limits from stack
          right_lim=stack_inds(stack_count)
          left_lim=stack_inds(stack_count-1)

          !Reduce stack count
          stack_count=stack_count-2
       else
          !What index are we looking at
          k=(left_lim+right_lim)/2

          !Swap values | KEY
          call swap_elem(key,k,left_lim+1)
          !Swap values | ARR1
          call swap_elem(arr1,k,left_lim+1)
          !Swap values | ARR2
          call swap_elem(arr2,k,left_lim+1)
          !Swap values | ARR3
          call swap_elem(arr3,k,left_lim+1)

          !Order so key(left_lim)<=key(left_lim+1)<=key(right_lim)
          if(key(left_lim).gt.key(right_lim)) then
             call swap_elem(key,left_lim,right_lim)
             call swap_elem(arr1,left_lim,right_lim)
             call swap_elem(arr2,left_lim,right_lim)
             call swap_elem(arr3,left_lim,right_lim)
          endif
          if(key(left_lim+1).gt.key(right_lim)) then
             call swap_elem(key,left_lim+1,right_lim)
             call swap_elem(arr1,left_lim+1,right_lim)
             call swap_elem(arr2,left_lim+1,right_lim)
             call swap_elem(arr3,left_lim+1,right_lim)
          endif
          if(key(left_lim).gt.key(left_lim+1)) then
             call swap_elem(key,left_lim,left_lim+1)
             call swap_elem(arr1,left_lim,left_lim+1)
             call swap_elem(arr2,left_lim,left_lim+1)
             call swap_elem(arr3,left_lim,left_lim+1)
          endif

          !Index of left limit of partition
          i=left_lim+1
          !Index of right limit of partition
          j=right_lim

          !Get partition values
          ak=key(left_lim+1)
          a1=arr1(left_lim+1)
          a2=arr2(left_lim+1)
          a3=arr3(left_lim+1)

          !Partition elements
          partitioned=.false.
          do while(.not.partitioned)
             !Search from left for first element bigger than ak
             i=i+1
             do while(key(i).lt.ak)
                i=i+1
             end do

             !Search from right for first element smaller than ak
             j=j-1
             do while(key(j).gt.ak)
                j=j-1
             end do

             !If first smaller element before first larger element then
             !partioning has finished --> Can now sort
             !Otherwise swap elements and try again
             if(j.lt.i)then
                !Indicate partioning complete
                partitioned=.true.

                !Insert partioning element in correct place
                key(left_lim+1)=key(j) ; key(j)=ak
                arr1(left_lim+1)=arr1(j) ; arr1(j)=a1
                arr2(left_lim+1)=arr2(j) ; arr2(j)=a2
                arr3(left_lim+1)=arr3(j) ; arr3(j)=a3
             else
                call swap_elem(key,i,j)
                call swap_elem(arr1,i,j)
                call swap_elem(arr2,i,j)
                call swap_elem(arr3,i,j)
             endif
          end do

          !Update stack counter
          stack_count=stack_count+2

          !At this point i is the start of the left section, which goes up to right_lim
          !and the right section goes from left_lim up to j
          !Update stack_inds values to record section done
          !If left section is larger push to stack and do right
          !otherwise push right to stack and do left
          if(right_lim-i+1.ge.j-left_lim) then
             !Store detail of left section
             stack_inds(stack_count)=right_lim
             stack_inds(stack_count-1)=i

             !Set right lim to one before start of left section
             right_lim=j-1
          else
             !Store details of right section
             stack_inds(stack_count)=j-1
             stack_inds(stack_count-1)=left_lim

             !Set left lim to one after end of right stack
             left_lim=i
          endif
       end if
    enddo
  end subroutine i_3_quicksort

  subroutine i_insertsort(arr_len,key)
    implicit none
    !Size of arrays
    integer, intent(in) :: arr_len

    !Arrays to be sorted
    integer, intent(inout), dimension(1:arr_len) :: key

    !Internals
    integer :: i, j
    integer :: ak !Used to hold selected elements
    logical :: done_exit

    do j=2,arr_len
       !Get current value
       ak=key(j)
       
       done_exit=.false.
       do i=j-1,1,-1
          !Escape loop
          if(key(i).le.ak)then
             done_exit=.true.
             exit
          endif
          
          !Shift values along one
          key(i+1)=key(i)
       enddo
       if(.not.done_exit) i=0
       
       !Insert value in correct place in section
       key(i+1)=ak
    enddo
  end subroutine i_insertsort

  subroutine i_1_insertsort(arr_len,key,arr1)
    implicit none
    !Size of arrays
    integer, intent(in) :: arr_len

    !Arrays to be sorted
    integer, intent(inout), dimension(1:arr_len) :: key, arr1

    !Internals
    integer :: i, j
    integer :: ak, a1 !Used to hold selected elements
    logical :: done_exit

    do j=2,arr_len
       !Get current value
       ak=key(j) ; a1=arr1(j)
       
       done_exit=.false.
       do i=j-1,1,-1
          !Escape loop
          if(key(i).le.ak)then
             done_exit=.true.
             exit
          endif
          
          !Shift values along one
          key(i+1)=key(i) ; arr1(i+1)=arr1(i)
       enddo
       if(.not.done_exit) i=0
       
       !Insert value in correct place in section
       key(i+1)=ak ; arr1(i+1)=a1
    enddo
  end subroutine i_1_insertsort

  subroutine i_2_insertsort(arr_len,key,arr1,arr2)
    implicit none
    !Size of arrays
    integer, intent(in) :: arr_len

    !Arrays to be sorted
    integer, intent(inout), dimension(1:arr_len) :: key, arr1, arr2

    !Internals
    integer :: i, j
    integer :: ak, a1, a2 !Used to hold selected elements
    logical :: done_exit

    do j=2,arr_len
       !Get current value
       ak=key(j) ; a1=arr1(j) ; a2=arr2(j)
       
       done_exit=.false.
       do i=j-1,1,-1
          !Escape loop
          if(key(i).le.ak)then
             done_exit=.true.
             exit
          endif
          
          !Shift values along one
          key(i+1)=key(i) ; arr1(i+1)=arr1(i) ; arr2(i+1)=arr2(i)
       enddo
       if(.not.done_exit) i=0
       
       !Insert value in correct place in section
       key(i+1)=ak ; arr1(i+1)=a1 ; arr2(i+1)=a2
    enddo
  end subroutine i_2_insertsort

  subroutine i_3_insertsort(arr_len,key,arr1,arr2,arr3)
    implicit none
    !Size of arrays
    integer, intent(in) :: arr_len

    !Arrays to be sorted
    integer, intent(inout), dimension(1:arr_len) :: key, arr1, arr2, arr3

    !Internals
    integer :: i, j
    integer :: ak, a1, a2, a3 !Used to hold selected elements
    logical :: done_exit

    do j=2,arr_len
       !Get current value
       ak=key(j) ; a1=arr1(j) ; a2=arr2(j) ; a3=arr3(j)
       
       done_exit=.false.
       do i=j-1,1,-1
          !Escape loop
          if(key(i).le.ak)then
             done_exit=.true.
             exit
          endif
          
          !Shift values along one
          key(i+1)=key(i) ; arr1(i+1)=arr1(i) ; arr2(i+1)=arr2(i) ; arr3(i+1)=arr3(i)
       enddo
       if(.not.done_exit) i=0
       
       !Insert value in correct place in section
       key(i+1)=ak ; arr1(i+1)=a1 ; arr2(i+1)=a2 ; arr3(i+1)=a3
    enddo
  end subroutine i_3_insertsort

  subroutine i_swap_elem(arr,ind1,ind2)
    implicit none
    integer, dimension(:), intent(inout) :: arr
    integer, intent(in) :: ind1, ind2
    integer :: tmp
    tmp=arr(ind1)
    arr(ind1)=arr(ind2)
    arr(ind2)=tmp
  end subroutine i_swap_elem
end module sorting
