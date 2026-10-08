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

   interface write_attribute
      module procedure write_ieee64_attribute
      module procedure write_int32_attribute
   end interface write_attribute

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

   subroutine write_ieee64_attribute(groupID, attrName, attrValue)
      use, intrinsic :: ieee_arithmetic, only : ieee_selected_real_kind

      integer, PARAMETER :: dp = ieee_selected_real_kind(p=15, r=307)

      integer(hid_t), intent(in) :: groupID
      character(*),   intent(in) :: attrName
      real(8),        intent(in) :: attrValue

      real(8)          :: test_r8
      real(dp), target :: attrValue_H5

      integer(hid_t) :: dsetID
      integer(hid_t) :: spaceID
      integer(hid_t) :: attrID
      integer        :: error

      ! Real(8) is not standard Fortran and, when a compiler accepts this, there
      ! is no guarantee that these variables will be IEEE double precision
      ! floating point variables.
      !
      ! Check that real(8) matches the data type that we use to write to file.
      if (KIND(attrValue_h5) /= KIND(test_r8)) then
        write(*,*) "HDF5 real kind does not match SkyAx real kind" 
        STOP
      end if

      ! Explicit cast to hardcoded data attribute data types just to be clear.
      ! Due to the above check, there is no need to check for casting errors
      ! or precision loss.
      attrValue_H5 = REAL(attrValue, kind=dp)

      CALL h5screate_f(H5S_SCALAR_F, spaceID, error)
      CALL h5acreate_f(groupID, TRIM(attrName), H5T_IEEE_F64LE, &
                       spaceID, attrID, error)
      CALL h5awrite_f(attrID, H5T_IEEE_F64LE, C_LOC(attrValue_H5), error)
      CALL h5aclose_f(attrID, error)
      CALL h5sclose_f(spaceID, error)
   end subroutine write_ieee64_attribute

   subroutine write_int32_attribute(groupID, attrName, attrValue)
      integer(hid_t), intent(in) :: groupID
      character(*),   intent(in) :: attrName
      integer,        intent(in) :: attrValue

      integer, target :: attrValue_H5

      integer(hid_t) :: dsetID
      integer(hid_t) :: spaceID
      integer(hid_t) :: attrID
      integer        :: error

      attrValue_H5 = attrValue

      CALL h5screate_f(H5S_SCALAR_F, spaceID, error)
      CALL h5acreate_f(groupID, TRIM(attrName), H5T_STD_I32LE, &
                       spaceID, attrID, error)
      CALL h5awrite_f(attrID, H5T_STD_I32LE, C_LOC(attrValue_H5), error)
      CALL h5aclose_f(attrID, error)
      CALL h5sclose_f(spaceID, error)
   end subroutine write_int32_attribute

end module hdf5_generic

module hdf5_skyax
   use hdf5
   use hdf5_generic
   use fields
   implicit none

   ! Real(8) is not standard Fortran and, when a compiler accepts this, there is
   ! no guarantee that these variables will be IEEE double precision floating
   ! point variables.
   !
   ! For now, we follow the type convention of this code (See m_observ.f90).
   type, public :: observables_t
       real(8) :: binding_energy_MeV
       real(8) :: diffraction_r_fm
       real(8) :: charge_r_rms_fm
       real(8) :: r_N_fm
       real(8) :: skin_thickness_fm
       real(8) :: neutron_gap_MeV
       real(8) :: proton_gap_MeV
       real(8) :: e_ls_p_MeV
       real(8) :: e_ls_n_MeV
       real(8) :: e_ls_n2_208Pb_MeV
       real(8) :: beta2_deformation
   end type observables_t

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

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Write the given SkyAx results to file as attributes.
    !
    !> @todo Write other results to file?
    !>       - Is there a values used to check convergence?
    !>       - Is there a status code?
    !> @todo Confirm units.
    ! 
    !> @param observables physical observable values to write to HDF5 file
    !> @param n_iters     total number of SkyAx iterations taken
    !---------------------------------------------------------------------------
    subroutine write_observables(observables, n_iters)
        character(len=*), PARAMETER :: FILENAME_H5 = "SkyAxObservables.h5"

        character(len=*), PARAMETER :: E_B_ATTR           = "BindingEnergy_MeV"
        character(len=*), PARAMETER :: DIFFRACTION_R_ATTR = "DiffractionRadius_fm"
        character(len=*), PARAMETER :: CHARGE_R_RMS_ATTR  = "ChargeRMSRadius_fm"
        character(len=*), PARAMETER :: R_N_ATTR           = "r_N_fm"
        character(len=*), PARAMETER :: SURF_THICK_ATTR    = "SkinThickness_fm"
        character(len=*), PARAMETER :: N_GAP_ATTR         = "NeutronPairingGap_MeV"
        character(len=*), PARAMETER :: P_GAP_ATTR         = "ProtonPairingGap_MeV"
        character(len=*), PARAMETER :: E_LS_P_ATTR        = "e_ls_p_MeV"
        character(len=*), PARAMETER :: E_LS_N_ATTR        = "e_ls_n_MeV"
        character(len=*), PARAMETER :: E_LS_N2_208PB_ATTR = "e_ls_n2_208Pb_MeV"
        character(len=*), PARAMETER :: BETA2_FINAL_ATTR   = "Beta2QuadrupoleDeformation"
        character(len=*), PARAMETER :: N_ITERS_ATTR       = "NumberIterations"

        type(observables_t), intent(IN) :: observables
        integer,             intent(IN) :: n_iters

        integer(hid_t) :: fileID
        integer        :: error

        error = 0

        CALL h5open_f(error)
        CALL h5fcreate_f(FILENAME_H5, H5F_ACC_EXCL_F, fileID, error)

        CALL write_attribute(fileID, E_B_ATTR,           observables%binding_energy_MeV)
        CALL write_attribute(fileID, DIFFRACTION_R_ATTR, observables%diffraction_r_fm)
        CALL write_attribute(fileID, CHARGE_R_RMS_ATTR,  observables%charge_r_rms_fm)
        CALL write_attribute(fileID, R_N_ATTR,           observables%r_N_fm)
        CALL write_attribute(fileID, SURF_THICK_ATTR,    observables%skin_thickness_fm)
        CALL write_attribute(fileID, N_GAP_ATTR,         observables%neutron_gap_MeV)
        CALL write_attribute(fileID, P_GAP_ATTR,         observables%proton_gap_MeV)
        CALL write_attribute(fileID, E_LS_P_ATTR,        observables%e_ls_p_MeV)
        CALL write_attribute(fileID, E_LS_N_ATTR,        observables%e_ls_n_MeV)
        CALL write_attribute(fileID, E_LS_N2_208PB_ATTR, observables%e_ls_n2_208Pb_MeV)
        CALL write_attribute(fileID, BETA2_FINAL_ATTR,   observables%beta2_deformation)
        CALL write_attribute(fileID, N_ITERS_ATTR,       n_iters)

        CALL h5close_f(error)
    end subroutine write_observables

end module hdf5_skyax
