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
      module procedure write_r8_array_1d
      module procedure write_r8_array_2d
      module procedure write_r8_array_3d
      module procedure write_r8_array_4d
   end interface write_array

contains

   subroutine write_r8_array_1d(arr,groupID,nm)
      implicit none
      real(8), intent(in) :: arr(:)
      integer(hid_t), intent(in) :: groupID
      character(*), intent(in) :: nm

      INTEGER(HID_T) :: dsetID, spaceID
      INTEGER(HSIZE_T), DIMENSION(1) :: adims
      integer :: error

      adims = shape(arr)

      CALL h5screate_simple_f(1, adims, spaceID, error)
      CALL h5dcreate_f(groupID, nm, H5T_NATIVE_DOUBLE, spaceID, dsetID, error)
      CALL h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, arr, adims, error)
      CALL h5dclose_f(dsetID, error)
   end subroutine write_r8_array_1d

   subroutine write_r8_array_2d(arr,groupID,nm)
      implicit none
      real(8), intent(in) :: arr(:,:)
      integer(hid_t), intent(in) :: groupID
      character(*), intent(in) :: nm

      INTEGER(HID_T) :: dsetID, spaceID
      INTEGER(HSIZE_T), DIMENSION(2) :: adims
      integer :: error

      adims = shape(arr)

      CALL h5screate_simple_f(2, adims, spaceID, error)
      CALL h5dcreate_f(groupID, nm, H5T_NATIVE_DOUBLE, spaceID, dsetID, error)
      CALL h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, arr, adims, error)
      CALL h5dclose_f(dsetID, error)
   end subroutine write_r8_array_2d

   subroutine write_r8_array_3d(arr,groupID,nm)
      implicit none
      real(8), intent(in) :: arr(:,:,:)
      integer(hid_t), intent(in) :: groupID
      character(*), intent(in) :: nm

      INTEGER(HID_T) :: dsetID, spaceID
      INTEGER(HSIZE_T), DIMENSION(3) :: adims
      integer :: error

      adims = shape(arr)

      CALL h5screate_simple_f(3, adims, spaceID, error)
      CALL h5dcreate_f(groupID, nm, H5T_NATIVE_DOUBLE, spaceID, dsetID, error)
      CALL h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, arr, adims, error)
      CALL h5dclose_f(dsetID, error)
   end subroutine write_r8_array_3d

   subroutine write_r8_array_4d(arr,groupID,nm)
      implicit none
      real(8), intent(in) :: arr(:,:,:,:)
      integer(hid_t), intent(in) :: groupID
      character(*), intent(in) :: nm

      INTEGER(HID_T) :: dsetID, spaceID
      INTEGER(HSIZE_T), DIMENSION(4) :: adims
      integer :: error

      adims = shape(arr)

      CALL h5screate_simple_f(4, adims, spaceID, error)
      CALL h5dcreate_f(groupID, nm, H5T_NATIVE_DOUBLE, spaceID, dsetID, error)
      CALL h5dwrite_f(dsetID, H5T_NATIVE_DOUBLE, arr, adims, error)
      CALL h5dclose_f(dsetID, error)
   end subroutine write_r8_array_4d

end module hdf5_generic

module hdf5_skyax
   use hdf5
   use hdf5_generic
   use fields
   implicit none

contains
!---------------------------------------------------------------------------  
! DESCRIPTION: write_fields
!> @brief
!!Saves static densities and fields to HDF5 file; see m_fields.f90 for 
!!definition of fields
!--------------------------------------------------------------------------- 
   subroutine write_fields
      implicit none

      integer :: error
      integer(HID_T) :: fileID, groupID

      ! TODO: move h5fcreate_f to separate routine
      call h5open_f(error)

      CALL h5fcreate_f(hdf5File, H5F_ACC_TRUNC_F, fileID, error)

      !-----------------
      ! Densities
      !-----------------
      call h5gcreate_f(fileID,'densities',groupID,error)

      call write_array(rho, groupID, 'rho')
      call write_array(crho, groupID, 'crho')
      call write_array(tau, groupID, 'tau')
      call write_array(drhols, groupID, 'drhols')
      call write_array(rholsr, groupID, 'rholsr')
      call write_array(rholsz, groupID, 'rholsz')
      call write_array(d2rho, groupID, 'd2rho')
      call write_array(rhoti, groupID, 'rhoti')
      call write_array(paircr, groupID, 'paircr')

      call h5gclose_f(groupID,error)

      !-----------------
      ! Fields
      !-----------------
      call h5gcreate_f(fileID,'fields',groupID,error)

      call write_array(pot, groupID, 'pot')
      call write_array(coulpt, groupID, 'coulpt')
      call write_array(bmass, groupID, 'bmass')
      call write_array(potls, groupID, 'potls')
      call write_array(drbmas, groupID, 'drbmas')
      call write_array(dzbmas, groupID, 'dzbmas')
      call write_array(potadd, groupID, 'potadd')
      call write_array(rhosur, groupID, 'rhosur')
      call write_array(pairpt, groupID, 'pairpt')
      call write_array(potdamp, groupID, 'potdamp')

      call h5gclose_f(groupID,error)

      !-----------------
      ! Wrapping up
      !-----------------

      CALL h5fclose_f(fileID, error)

      call h5close_f(error)

   end subroutine write_fields
end module hdf5_skyax
