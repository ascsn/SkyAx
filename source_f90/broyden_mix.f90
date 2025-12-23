MODULE broyden_mix_mod
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: broyden_vsqf_mix

  ! Persistent history and workspace to avoid frequent allocations
  REAL(8), ALLOCATABLE, SAVE :: x_hist(:,:), r_hist(:,:)
  REAL(8), ALLOCATABLE, SAVE :: df(:,:), dx(:,:)
  REAL(8), ALLOCATABLE, SAVE :: a(:,:), b(:), c(:)
  INTEGER, SAVE :: nvar_save = 0, msave = 0, nhist = 0

CONTAINS

  SUBROUTINE broyden_vsqf_mix(iter, x_old, x_raw, alpha, m_hist)
    INTEGER, INTENT(IN) :: iter, m_hist
    REAL(8), INTENT(IN) :: x_old(:)     ! Input Master Vector (rho, rhot, lambda)
    REAL(8), INTENT(IN OUT) :: x_raw(:) ! Output from HFB step -> optimized x_{n+1}
    REAL(8), INTENT(IN) :: alpha

    INTEGER :: nvar, k, j, i, info
    REAL(8), ALLOCATABLE :: r(:)
    REAL(8) :: regularization

    ! Increased regularization to handle nearly linear dependence near shell closures
    regularization = 1D-8 
    nvar = SIZE(x_old)

    ! 1. Allocation Logic (Optimized to occur only once)
    IF (nvar_save /= nvar .OR. msave /= m_hist) THEN
      nvar_save = nvar
      msave = m_hist
      nhist = 0
      IF (ALLOCATED(x_hist)) DEALLOCATE(x_hist, r_hist, df, dx, a, b, c)
      ALLOCATE(x_hist(nvar, m_hist+1), r_hist(nvar, m_hist+1))
      ALLOCATE(df(nvar, m_hist), dx(nvar, m_hist))
      ALLOCATE(a(m_hist, m_hist), b(m_hist), c(m_hist))
    END IF

    ALLOCATE(r(nvar))
    r = x_raw - x_old   ! residual r_n

    ! 2. Update history buffers
    IF (nhist < m_hist + 1) THEN
      nhist = nhist + 1
    ELSE
      ! Shift history left to make room for newest entry
      x_hist(:, 1:m_hist) = x_hist(:, 2:m_hist+1)
      r_hist(:, 1:m_hist) = r_hist(:, 2:m_hist+1)
    END IF
    x_hist(:, nhist) = x_old
    r_hist(:, nhist) = r

    ! 3. Fallback to Linear Mixing for the first few steps
    ! (Wait until we have enough history for curvature calculation)
    IF (nhist < 3) THEN
      x_raw = x_old + alpha * r
      DEALLOCATE(r)
      RETURN
    END IF

    ! 4. Construct Difference Matrices
    k = MIN(nhist - 1, m_hist)
    DO j = 1, k
      df(:, j) = r_hist(:, nhist-k+j) - r_hist(:, nhist-k+j-1)
      dx(:, j) = x_hist(:, nhist-k+j) - x_hist(:, nhist-k+j-1)
    END DO

    ! 5. Solve Normal Equations: (df^T df) c = df^T r_n
    a(1:k, 1:k) = 0D0
    b(1:k) = 0D0
    DO i = 1, k
      b(i) = DOT_PRODUCT(df(:, i), r)
      DO j = 1, k
        a(i, j) = DOT_PRODUCT(df(:, i), df(:, j))
      END DO
      ! Regularization prevents singularity when system hits a shell-closure stall
      a(i, i) = a(i, i) + regularization 
    END DO

    CALL solve_sym_posdef(k, a(1:k, 1:k), b(1:k), c(1:k), info)

    IF (info /= 0) THEN
      ! Fallback to linear mixing if matrix is ill-conditioned
      x_raw = x_old + alpha * r
    ELSE
      ! 6. Anderson Update: x_{n+1} = x_old + alpha*r - sum[ c_j * (dx_j + alpha*df_j) ]
      x_raw = x_old + alpha * r
      DO j = 1, k
        x_raw = x_raw - (dx(:, j) + alpha * df(:, j)) * c(j)
      END DO
    END IF

    DEALLOCATE(r)

  END SUBROUTINE broyden_vsqf_mix

  SUBROUTINE solve_sym_posdef(n, a, b, x, info)
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN OUT) :: a(n,n)
    REAL(8), INTENT(IN) :: b(n)
    REAL(8), INTENT(OUT) :: x(n)
    INTEGER, INTENT(OUT) :: info
    INTEGER :: i, j, k, piv
    REAL(8) :: maxv, tmp, factor
    REAL(8), ALLOCATABLE :: rhs(:)

    ALLOCATE(rhs(n))
    rhs = b
    info = 0
    DO k = 1, n
      piv = k
      maxv = ABS(a(k, k))
      DO i = k+1, n
        IF (ABS(a(i, k)) > maxv) THEN
          maxv = ABS(a(i, k)); piv = i
        END IF
      END DO
      IF (maxv < 1D-18) THEN
        info = 1
        DEALLOCATE(rhs)
        RETURN
      END IF
      IF (piv /= k) THEN
        DO j = k, n
          tmp = a(k, j); a(k, j) = a(piv, j); a(piv, j) = tmp
        END DO
        tmp = rhs(k); rhs(k) = rhs(piv); rhs(piv) = tmp
      END IF
      DO i = k+1, n
        factor = a(i, k) / a(k, k)
        a(i, k) = 0D0
        DO j = k+1, n
          a(i, j) = a(i, j) - factor * a(k, j)
        END DO
        rhs(i) = rhs(i) - factor * rhs(k)
      END DO
    END DO
    DO i = n, 1, -1
      tmp = rhs(i)
      DO j = i+1, n
        tmp = tmp - a(i, j) * x(j)
      END DO
      x(i) = tmp / a(i, i)
    END DO
    DEALLOCATE(rhs)
  END SUBROUTINE solve_sym_posdef

END MODULE broyden_mix_mod