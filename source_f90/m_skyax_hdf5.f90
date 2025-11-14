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

contains

   subroutine write_float_attr(nm,dat,grp_id)
      implicit none
      character(len=*), intent(in) :: nm
      real(8), intent(in)         :: dat
      integer(hid_t), intent(in)   :: grp_id

      INTEGER(HID_T) :: attr_id, aspace_id
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
      integer :: error

      CALL h5screate_simple_f(1, adims, aspace_id, error)
      call error_exit(error,'h5screate_simple_f')

      CALL h5acreate_f(grp_id, nm, H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
      call error_exit(error,'h5acreate_f')

      CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, dat, adims, error)
      call error_exit(error,'h5awrite_f')

      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(aspace_id, error)

   end subroutine write_float_attr
end module hdf5_generic
