!------------------------------------------------------------------------------
! MODULE: hdf5_skyax
!------------------------------------------------------------------------------
! DESCRIPTION:
!> @brief
!!This module contains generic HDF5 routines, and specific wrapper routines
!!for IO in SkyAx
!> @details
!!Will add once there's details to use
!------------------------------------------------------------------------------

module hdf5_generic
   use hdf5
   implicit none

   ! later will set programatically
   character(15), parameter :: hdf5File = "skyax_000005.h5"

   interface write_array
      module procedure write_array_1d
      module procedure write_array_2d
      module procedure write_array_3d
      module procedure write_array_4d
   end interface write_array

contains

   subroutine write_array_1d(a)
      implicit none
      real(8), intent(in) :: a(:)

      write(6,*) 'Received array of shape ',shape(a)

      ! print *, "Processing 1D array of size ", size(a)

      ! ... user code here ...
   end subroutine write_array_1d

   subroutine write_array_2d(a)
      implicit none
      real(8), intent(in) :: a(:,:)

      write(6,*) 'Received array of shape ',shape(a)

      ! print *, "Processing 1D array of size ", size(a)

      ! ... user code here ...
   end subroutine write_array_2d

   subroutine write_array_3d(a)
      implicit none
      real(8), intent(in) :: a(:,:,:)

      write(6,*) 'Received array of shape ',shape(a)

      ! print *, "Processing 1D array of size ", size(a)

      ! ... user code here ...
   end subroutine write_array_3d

   subroutine write_array_4d(a)
      implicit none
      real(8), intent(in) :: a(:,:,:,:)

      write(6,*) 'Received array of shape ',shape(a)

      ! print *, "Processing 1D array of size ", size(a)

      ! ... user code here ...
   end subroutine write_array_4d

   ! subroutine write_float_attr(nm,dat,grp_id)
   !    implicit none
   !    character(len=*), intent(in) :: nm
   !    real(8), intent(in)         :: dat
   !    integer(hid_t), intent(in)   :: grp_id

   !    INTEGER(HID_T) :: attr_id, aspace_id
   !    INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
   !    integer :: error

   !    CALL h5screate_simple_f(1, adims, aspace_id, error)
   !    ! call error_exit(error,'h5screate_simple_f')

   !    CALL h5acreate_f(grp_id, nm, H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
   !    ! call error_exit(error,'h5acreate_f')

   !    CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, dat, adims, error)
   !    ! call error_exit(error,'h5awrite_f')

   !    CALL h5aclose_f(attr_id, error)
   !    CALL h5sclose_f(aspace_id, error)

   ! end subroutine write_float_attr
end module hdf5_generic

module hdf5_skyax
   use hdf5_generic

contains
   subroutine write_fields
      implicit none

      real(8), dimension(4,2) :: testArr

      call write_array(testArr)

   end subroutine write_fields
end module hdf5_skyax
