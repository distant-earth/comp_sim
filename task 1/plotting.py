import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib
#plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 18})
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Times New Roman'
matplotlib.rcParams['mathtext.it'] = 'Times New Roman:italic'
matplotlib.rcParams['mathtext.bf'] = 'Times New Roman:bold'

def InitializeParameters():
	with open('INPUT') as data:
		ecc = float(data.readline().split()[0])
		T_ini = float(data.readline().split()[0])
		T_fin = float(data.readline().split()[0])
		T_step = float(data.readline().split()[0])
		dT0 = float(data.readline().split()[0])
		time = float(data.readline().split()[0])
	return ecc, T_ini, T_fin, T_step, dT0, time

def MakePlot(T_ini, T_fin, T_step):
	plt.figure(figsize = (14,7))
	plt.gca().xaxis.set_major_locator(MultipleLocator(0.5))
	plt.gca().xaxis.set_major_formatter('{x:.1f}')
	plt.gca().xaxis.set_minor_locator(MultipleLocator(0.1))
	plt.gca().yaxis.set_major_locator(MultipleLocator(0.5))
	plt.gca().yaxis.set_major_formatter('{x:.1f}')
	plt.gca().yaxis.set_minor_locator(MultipleLocator(0.1))
	plt.xlabel(r'$T$')
	plt.ylabel(r"$T$ '")
	plt.xlim([-2,2])
	plt.ylim([-2,2])
	print("\n Подготовка графика (T',T)...")
	for k in range(int((T_fin - T_ini) / T_step) + 1):
		t, T, dT = ReadResult(k, 'T')
		T_points, dT_points = Interpolate(t, T, dT)
		plt.plot(T_points, dT_points, '.', markersize = 1)
	plt.tight_layout()
	plt.savefig('pictures/poincare_T.png', format='png', dpi = 300)
	plt.figure(figsize = (14,7))
	plt.gca().xaxis.set_major_locator(MultipleLocator(0.5))
	plt.gca().xaxis.set_major_formatter('{x:.1f}')
	plt.gca().xaxis.set_minor_locator(MultipleLocator(0.1))
	plt.gca().yaxis.set_major_locator(MultipleLocator(0.5))
	plt.gca().yaxis.set_major_formatter('{x:.1f}')
	plt.gca().yaxis.set_minor_locator(MultipleLocator(0.1))
	plt.xlabel(r'$z$')
	plt.ylabel(r"$z$ '")
	plt.xlim([-2,2])
	plt.ylim([-2,2])
	print("\n Подготовка графика (z',z)...")
	for k in range(int((T_fin - T_ini) / T_step) + 1):
		t, z, dz = ReadResult(k, 'z')
		z_points, dz_points = Interpolate(t, z, dz)
		plt.plot(z_points, dz_points, '.', markersize = 1)
	plt.tight_layout()
	plt.savefig('pictures/poincare_z.png', format='png', dpi = 300)
	plt.show()

def ReadResult(k, par):
	t = []; arr = []; d_arr = []
	if (k < 10):
		filename = 'results/RESULT00' + str(k)
	elif (k < 100):
		filename = 'results/RESULT0' + str(k)
	else:
		filename = 'results/RESULT' + str(k)
	with open(filename) as table:
		print('Обрабатывается:', filename[8:])
		for line in table:
			t.append(float(line.split()[0]))
			arr.append(float(line.split()[1 + 2 * int(par == 'z')]))
			d_arr.append(float(line.split()[2 + 2 * int(par == 'z')]))
	return t, arr, d_arr
					
def Interpolate(t, T, dT):
	T_points = []; dT_points = []
	N = 1
	for i in range(len(t)-10):
		if (t[i] <= 2 * np.pi * N <= t[i+1]):
			T_func = interp1d(t[i-10:i+10], T[i-10:i+10])
			dT_func = interp1d(t[i-10:i+10], dT[i-10:i+10])
			T_points.append(T_func(2 * np.pi * N))
			dT_points.append(dT_func(2 * np.pi * N))
			N += 1
	return T_points, dT_points

ecc, T_ini, T_fin, T_step, dT0, time = InitializeParameters()
MakePlot(T_ini, T_fin, T_step)
