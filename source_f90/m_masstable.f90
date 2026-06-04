!------------------------------------------------------------------------------
! MODULE: masstable
!------------------------------------------------------------------------------
! DESCRIPTION:
!> @brief
!! Module for mass table and stable-line (dripline) large-scale calculations.
!! Adapted from HFBTHO_v3.00 hfbtho_large_scale.f90 for the SkyAx code.
!!
!! run_mode=0 : single nucleus (default, existing SkyAx behaviour)
!! run_mode=1 : mass table — loop over (Z,N,beta) entries in skyax_MASSTABLE.dat
!! run_mode=2 : stable-line — scan isotopic/isotonic chains from skyax_STABLELINE.dat
!!
!! For parallel (SLURM array) runs set mt_start / mt_end in the &main namelist
!! so that each task processes a disjoint subset of the table rows.
!------------------------------------------------------------------------------
MODULE masstable

  IMPLICIT NONE

  ! ----- run control ---------------------------------------------------------
  INTEGER :: run_mode = 0  !< 0=single, 1=masstable, 2=stableline

  ! ----- mass-table input (skyax_MASSTABLE.dat) ------------------------------
  !  Format:
  !   <header line>
  !   <nRows>          ← N: total number of runs
  !   <column header>
  !   Z  N  beta0       (one row per nucleus; beta0 = initial deformation)
  !
  INTEGER :: mt_nRows = 0    !< Total rows in mass-table file
  INTEGER :: mt_start = 1    !< First row to process (for SLURM array jobs)
  INTEGER :: mt_end   = -1   !< Last  row to process (-1 means mt_nRows)
  LOGICAL :: mt_do_oddA   = .TRUE.  !< .TRUE. = include odd-A nuclei   (Z+N odd)
  LOGICAL :: mt_do_oddodd = .TRUE.  !< .TRUE. = include odd-odd nuclei (Z odd and N odd)
  INTEGER, ALLOCATABLE :: mt_Z(:)    !< Proton  numbers
  INTEGER, ALLOCATABLE :: mt_N(:)    !< Neutron numbers
  REAL(8), ALLOCATABLE :: mt_beta(:) !< Initial deformation (betain)

  ! ----- stable-line input (skyax_STABLELINE.dat) ----------------------------
  !  Format:
  !   <header line>
  !   <nRows>
  !   <column header>
  !   Z  N  direction   (direction=1 → isotopic, direction=-1 → isotonic)
  !
  INTEGER :: sl_nRows = 0
  INTEGER, ALLOCATABLE :: sl_Z(:), sl_N(:), sl_dir(:)

  ! ----- mass-table output arrays (allocated for mt_end-mt_start+1) ----------
  INTEGER :: mt_nout = 0   !< Number of entries stored so far
  INTEGER,  ALLOCATABLE :: mt_Z_out(:), mt_N_out(:)
  REAL(8),  ALLOCATABLE :: mt_beta_in(:)
  ! Deformation
  REAL(8),  ALLOCATABLE :: mt_beta2_out(:), mt_beta2p_out(:), mt_beta2n_out(:)
  REAL(8),  ALLOCATABLE :: mt_beta4_out(:), mt_beta4p_out(:), mt_beta4n_out(:)
  ! Energy and pairing
  REAL(8),  ALLOCATABLE :: mt_E_out(:)
  REAL(8),  ALLOCATABLE :: mt_epairp_out(:), mt_epairn_out(:)
  REAL(8),  ALLOCATABLE :: mt_ferp_out(:),   mt_fern_out(:)
  REAL(8),  ALLOCATABLE :: mt_gapp_out(:),   mt_gapn_out(:)
  ! Iteration and angular momentum
  INTEGER,  ALLOCATABLE :: mt_iter_out(:)
  REAL(8),  ALLOCATABLE :: mt_Kp_out(:), mt_Kn_out(:)
  REAL(8),  ALLOCATABLE :: mt_Jp_out(:), mt_Jn_out(:)
  ! Radii
  REAL(8),  ALLOCATABLE :: mt_rmsc_out(:), mt_rms_out(:)
  REAL(8),  ALLOCATABLE :: mt_rmsp_out(:), mt_rmsn_out(:)
  ! Bound nucleons, variance, convergence
  REAL(8),  ALLOCATABLE :: mt_boundp_out(:), mt_boundn_out(:)
  REAL(8),  ALLOCATABLE :: mt_evar_out(:)   !< E_variance (devh2)
  REAL(8),  ALLOCATABLE :: mt_econ_out(:)   !< dE/diter   (econ)

CONTAINS

  !==========================================================================
  !> Reads skyax_MASSTABLE.dat.  Sets mt_start/mt_end to the valid range.
  !==========================================================================
  SUBROUTINE read_skyax_MassTable
    IMPLICIT NONE
    INTEGER :: lu, i
    lu = 15
    OPEN(lu, FILE='skyax_MASSTABLE.dat', STATUS='old', ACTION='read')
    READ(lu, *)           ! comment/header
    READ(lu, *) mt_nRows  ! N: total number of runs
    READ(lu, *)           ! column header
    ALLOCATE(mt_Z(mt_nRows), mt_N(mt_nRows), mt_beta(mt_nRows))
    DO i = 1, mt_nRows
      READ(lu, *) mt_Z(i), mt_N(i), mt_beta(i)
    END DO
    CLOSE(lu)
    ! Clamp the processing range
    IF (mt_end < 1 .OR. mt_end > mt_nRows) mt_end = mt_nRows
    mt_start = MAX(mt_start, 1)
    IF (mt_start > mt_end) THEN
      WRITE(*, '(a,2i6)') ' WARNING: mt_start > mt_end, resetting to 1:', mt_start, mt_end
      mt_start = 1
    END IF
    WRITE(*, '(a,i6,a,i6,a,i6)') &
      ' MassTable: ', mt_nRows, ' entries read; processing rows ', mt_start, ' to ', mt_end
  END SUBROUTINE read_skyax_MassTable

  !==========================================================================
  !> Reads skyax_STABLELINE.dat.
  !==========================================================================
  SUBROUTINE read_skyax_StableLine
    IMPLICIT NONE
    INTEGER :: lu, i
    lu = 15
    OPEN(lu, FILE='skyax_STABLELINE.dat', STATUS='old', ACTION='read')
    READ(lu, *)
    READ(lu, *) sl_nRows
    READ(lu, *)
    ALLOCATE(sl_Z(sl_nRows), sl_N(sl_nRows), sl_dir(sl_nRows))
    DO i = 1, sl_nRows
      READ(lu, *) sl_Z(i), sl_N(i), sl_dir(i)
    END DO
    CLOSE(lu)
    WRITE(*, '(a,i6,a)') ' StableLine: ', sl_nRows, ' chains to scan'
  END SUBROUTINE read_skyax_StableLine

  !==========================================================================
  !> Allocates output arrays for n entries.
  !==========================================================================
  SUBROUTINE alloc_mt_output(n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    IF (ALLOCATED(mt_Z_out)) THEN
      DEALLOCATE(mt_Z_out, mt_N_out, mt_beta_in)
      DEALLOCATE(mt_beta2_out, mt_beta2p_out, mt_beta2n_out)
      DEALLOCATE(mt_beta4_out, mt_beta4p_out, mt_beta4n_out)
      DEALLOCATE(mt_E_out, mt_epairp_out, mt_epairn_out)
      DEALLOCATE(mt_ferp_out, mt_fern_out, mt_gapp_out, mt_gapn_out)
      DEALLOCATE(mt_iter_out, mt_Kp_out, mt_Kn_out, mt_Jp_out, mt_Jn_out)
      DEALLOCATE(mt_rmsc_out, mt_rms_out, mt_rmsp_out, mt_rmsn_out)
      DEALLOCATE(mt_boundp_out, mt_boundn_out, mt_evar_out, mt_econ_out)
    END IF
    ALLOCATE(mt_Z_out(n), mt_N_out(n), mt_beta_in(n))
    ALLOCATE(mt_beta2_out(n), mt_beta2p_out(n), mt_beta2n_out(n))
    ALLOCATE(mt_beta4_out(n), mt_beta4p_out(n), mt_beta4n_out(n))
    ALLOCATE(mt_E_out(n), mt_epairp_out(n), mt_epairn_out(n))
    ALLOCATE(mt_ferp_out(n), mt_fern_out(n), mt_gapp_out(n), mt_gapn_out(n))
    ALLOCATE(mt_iter_out(n), mt_Kp_out(n), mt_Kn_out(n), mt_Jp_out(n), mt_Jn_out(n))
    ALLOCATE(mt_rmsc_out(n), mt_rms_out(n), mt_rmsp_out(n), mt_rmsn_out(n))
    ALLOCATE(mt_boundp_out(n), mt_boundn_out(n), mt_evar_out(n), mt_econ_out(n))
    mt_nout = 0
  END SUBROUTINE alloc_mt_output

  !==========================================================================
  !> Appends one completed nucleus to the output arrays.
  !==========================================================================
  SUBROUTINE fill_mt_output(Z, N, beta_in, &
      b2, b2p, b2n, b4, b4p, b4n, &
      E, epairp, epairn, ferp, fern, gapp, gapn, &
      niter, Kp, Kn, Jp, Jn, &
      rms_c, rms_t, rms_p, rms_n, &
      Nboundp, Nboundn, Evar, dEcon)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Z, N, niter
    REAL(8), INTENT(IN) :: beta_in
    REAL(8), INTENT(IN) :: b2, b2p, b2n, b4, b4p, b4n
    REAL(8), INTENT(IN) :: E, epairp, epairn, ferp, fern, gapp, gapn
    REAL(8), INTENT(IN) :: Kp, Kn, Jp, Jn
    REAL(8), INTENT(IN) :: rms_c, rms_t, rms_p, rms_n
    REAL(8), INTENT(IN) :: Nboundp, Nboundn, Evar, dEcon

    mt_nout = mt_nout + 1
    mt_Z_out(mt_nout)      = Z
    mt_N_out(mt_nout)      = N
    mt_beta_in(mt_nout)    = beta_in
    mt_beta2_out(mt_nout)  = b2
    mt_beta2p_out(mt_nout) = b2p
    mt_beta2n_out(mt_nout) = b2n
    mt_beta4_out(mt_nout)  = b4
    mt_beta4p_out(mt_nout) = b4p
    mt_beta4n_out(mt_nout) = b4n
    mt_E_out(mt_nout)      = E
    mt_epairp_out(mt_nout) = epairp
    mt_epairn_out(mt_nout) = epairn
    mt_ferp_out(mt_nout)   = ferp
    mt_fern_out(mt_nout)   = fern
    mt_gapp_out(mt_nout)   = gapp
    mt_gapn_out(mt_nout)   = gapn
    mt_iter_out(mt_nout)   = niter
    mt_Kp_out(mt_nout)     = Kp
    mt_Kn_out(mt_nout)     = Kn
    mt_Jp_out(mt_nout)     = Jp
    mt_Jn_out(mt_nout)     = Jn
    mt_rmsc_out(mt_nout)   = rms_c
    mt_rms_out(mt_nout)    = rms_t
    mt_rmsp_out(mt_nout)   = rms_p
    mt_rmsn_out(mt_nout)   = rms_n
    mt_boundp_out(mt_nout) = Nboundp
    mt_boundn_out(mt_nout) = Nboundn
    mt_evar_out(mt_nout)   = Evar
    mt_econ_out(mt_nout)   = dEcon
  END SUBROUTINE fill_mt_output

  !==========================================================================
  !> Writes the collected mass table to outfile.
  !! The first line contains N (the total number of rows written).
  !==========================================================================
  SUBROUTINE print_mass_table(outfile)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: outfile
    INTEGER :: i, lu
    lu = 100
    OPEN(lu, FILE=outfile, STATUS='unknown')
    ! N value at top: number of entries in this file
    WRITE(lu, '(a,i6)') 'N=', mt_nout
    WRITE(lu, '(a)') &
      '#   i   Z    N   beta_in' // &
      '   bet2    bet2p   bet2n   bet4    bet4p   bet4n' // &
      '      energy       epairp   epairn' // &
      '    ferp     fern    gapp    gapn' // &
      '  iter   K_p   K_n   J_p   J_n' // &
      '    rmsc     rms     rmsp    rmsn' // &
      '  N_boundp N_boundn  E_variance   dE/diter'
    DO i = 1, mt_nout
      WRITE(lu, '(3i5,f9.4, 6f8.4, f18.6, 2f9.3, 2f8.3, 2f8.4, i6, 4f6.2, 4f8.4, 2f9.3, 2(1pg12.4))') &
          i, mt_Z_out(i), mt_N_out(i), mt_beta_in(i), &
          mt_beta2_out(i), mt_beta2p_out(i), mt_beta2n_out(i), &
          mt_beta4_out(i), mt_beta4p_out(i), mt_beta4n_out(i), &
          mt_E_out(i), mt_epairp_out(i), mt_epairn_out(i), &
          mt_ferp_out(i), mt_fern_out(i), mt_gapp_out(i), mt_gapn_out(i), &
          mt_iter_out(i), mt_Kp_out(i), mt_Kn_out(i), mt_Jp_out(i), mt_Jn_out(i), &
          mt_rmsc_out(i), mt_rms_out(i), mt_rmsp_out(i), mt_rmsn_out(i), &
          mt_boundp_out(i), mt_boundn_out(i), mt_evar_out(i), mt_econ_out(i)
    END DO
    CLOSE(lu)
    WRITE(*, '(a,a)') ' MassTable output written to ', TRIM(outfile)
  END SUBROUTINE print_mass_table

END MODULE masstable
