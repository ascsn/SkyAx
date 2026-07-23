REAL :: dum(11),eqmin(0:3)
CHARACTER(80) :: ain

DO i=1,99999
  READ(5,'(a)') ain
  IF(INDEX(ain,'y (')>0) WRITE(6,'(a)') ain
  IF(INDEX(ain,'beta_3')>0) THEN
    READ(5,'(1x/1x/a)') ain
    WRITE(6,'(a)') ain
  END IF
  IF(INDEX(ain,'tron  lev')>0) THEN
    READ(5,'(1x/1x/1x/1x/1x)')
    eqmin = 1E3
    DO j=1,200
      READ(5,*,ERR=19,END=99) nr,iso,xjz,i1,i2,i3, &
           dum(1:4),spenrg,dum(5:7),equasi
      jact = NINT(xjz-0.5)
      IF(equasi<eqmin(jact)) eqmin(jact)=equasi
    END DO
    19 CONTINUE
    WRITE(6,'(4f8.3)') eqmin
  END IF
END DO

99 STOP
END
