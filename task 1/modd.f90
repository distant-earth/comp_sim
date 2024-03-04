module modd

implicit none

real(8), parameter :: GM = 1.d0
real(8), parameter :: pi = 4 * atan(1.0d0)
real(8), parameter :: period = 2 * pi
real(8), parameter :: n_mean = 1.d0
real(8) :: ecc, T_ini, T_fin, T_step, dT, N_revol

contains

subroutine InitializeParameters()
open(unit = 1, file = 'INPUT')
read(1,*) ecc
read(1,*) T_ini
read(1,*) T_fin
read(1,*) T_step
read(1,*) dT
read(1,*) N_revol
end subroutine InitializeParameters

function KeplerEquation(t)
real(8) :: t, M, E_old, E, delta, KeplerEquation
! Решаем уравнение Кеплера и находим E:
M = n_mean * t
E_old = M
delta = 1.d0
do while (delta >= 1e-10)
	E = M + ecc * sin(E_old)
	delta = abs(E - E_old)
	E_old = E
enddo
KeplerEquation = E
return
end function KeplerEquation

! Правые части уравнений:
subroutine func_right_part(N, t, y, res, RPAR, IPAR)
integer :: N, IPAR
real(8) :: res(N), y(N), t, RPAR(2)
real(8) :: B_current, T_current
real(8) :: sin_nu, cos_nu
real(8) :: E
T_current = y(1)
B_current = y(2)
E = KeplerEquation(t)
! Косинус и синус истинной аномалии через E:
sin_nu = (sqrt(1 - ecc**2) * sin(E)) / (1 - ecc * cos(E))
cos_nu = (cos(E) - ecc) / (1 - ecc * cos(E))
res(1) = B_current
res(2) = - ((0.25d0 + T_current**2)**(-1.5d0) + ecc * cos_nu) * T_current / (1 + ecc * cos_nu)
end subroutine func_right_part

! Подпрограмма для вывода промежуточных значений в файл:
subroutine solout(NR, t_old, t, y, N, con, icomp, ND, RPAR, IPAR, IRTRN)
integer :: N, NR, ND, IPAR, IRTRN
real(8) :: y(N), t_old, t, con, icomp, RPAR(2)
real(8) :: T_current, B_current, E, sin_nu, cos_nu, z, dz
T_current = y(1)
B_current = y(2)
E = KeplerEquation(t)
! Косинус и синус истинной аномалии через E:
sin_nu = (sqrt(1 - ecc**2) * sin(E)) / (1 - ecc * cos(E))
cos_nu = (cos(E) - ecc) / (1 - ecc * cos(E))
z = T_current * (1 - ecc**2) / (1 + ecc * cos_nu)
dz = (T_current * ecc * sin_nu + dT * (1 + ecc * cos_nu)) / sqrt(1 - ecc**2)
write(IPAR, *) t, y(1), y(2), z, dz
end subroutine solout


end module modd
