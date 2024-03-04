program main
use modd

implicit none
integer, parameter :: N = 2, LWORK = 11 * N + 20, LIWORK = 20
real(8), parameter :: RTOL = 1e-14, ATOL = 1e-14, & ! локальная ошибка Y(I) не превосходит RTOL*ABS(Y(I))+ATOL
                                       ITOL = 0 ! RTOL и ATOL скаляры (если 1, то векторы и RTOL(I)*ABS(Y(I))+ATOL(I)
real(8) :: WORK(LWORK), RPAR(2), y(N)
integer :: IDID, IWORK(LIWORK), IPAR, IOUT
real(8) :: time, t0
integer :: k, N_surfaces
character(10) :: filename

call InitializeParameters()
N_surfaces = (T_fin - T_ini) / T_step
write(*,*)
write(*,*) 'Число шагов по T:', N_surfaces
do k = 0, N_surfaces
	!=================================================
	! Параметры, влияющие на работу интегратора:
	IOUT = 1 ! 1 - SOLOUT осуществляет вывод; 0 - нет
	WORK = 0.0d0
	!WORK(6) = 0.001d0
	IWORK = 0
	IWORK(1) = 2147483647 ! максимальное допустимое кол-во шагов
	! Поставлено максимальное возможное число, при значении по умолчанию
	! интегратор выдает ошибку "Недостаточно шагов".
	! WORK(5) = 0.04d0 ! beta for dense output
	! WORK(3) = 0.333d0 ! ограничение снизу на отношение нового шага к старому
	! WORK(4) = 6.d0 ! ограничение сверху на отношение нового шага к старому
	! WORK(6) максимальный размер шага, по умолчанию XEND-X.
	! WORK(7) = 0.d0 ! начальный размер шага
	! В методе DOPRI шаг подбирается автоматически. Если результат выглядит 
	! неадекватно, желательно сначала поискать ошибку в шапке main.f90.
	! Например, проверить RTOL и ATOL.
	!=================================================
	t0 = 0.d0
	time = N_revol * period
	IPAR = k * 10
	if (k < 10) then
		write(filename, '(a6,i1,i1,i1)') 'RESULT', 0, 0, k
	else if (k < 100) then
		write(filename, '(a6,i1,i2)') 'RESULT', 0, k
	else
		write(filename, '(a6,i3)') 'RESULT', k
	endif
	open(unit = IPAR, file = trim(filename))
	y(1) = T_ini + k * T_step
	y(2) = dT
	call DOP853(N, func_right_part, t0, y, time, RTOL, ATOL, ITOL, solout, IOUT, WORK, LWORK, IWORK, LIWORK, RPAR, IPAR, IDID)
	write(*,*) 'Шаг:', k, 'Exit status:', IDID
enddo
end

