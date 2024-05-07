import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import json
from mpl_toolkits.basemap import Basemap
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import average, fcluster

def InitializeParams():
	print('\nInitializing parameters...')
	with open('input.json', 'r') as input_file:
		params = json.load(input_file)
	month = str(params['observation file (yyyymm)'])
	n_meteors = params['number of meteors to extract (integer/all)']
	proj = params['map projection']
	center = params['central longtitude of the map']
	metric = params['metric (default: cosine)']
	if str(n_meteors).strip().lower() == 'all':
		n_meteors = None
	threshold = params['threshold distance']
	print('done.\n')
	return month, n_meteors, proj, center, metric, threshold

def GetData(month, n_meteors):
	print('Extracting data...')
	filename = 'data/traj_summary_monthly_' + month + '.txt'
	with open(filename, 'r') as datafile:
		if n_meteors is not None:
			row_count = -4 #headers
			for line in datafile:
				if line.strip(): #exclude empty lines
					row_count += 1
			if n_meteors > row_count:
				n_meteors = None
				print(f'using max {row_count} lines available...')
	data = pd.read_table(filename, nrows=n_meteors, delimiter=';', skiprows=4, header=None, usecols = [7,9], names=['ra', 'dec'])
	ra = data['ra'].to_list()
	dec = data['dec'].to_list()
	print('done.\n')
	return ra, dec

def PlotMapGrid(proj="moll", center=270):
	plt.rcParams.update({"text.usetex": True})
	plt.figure(figsize=(8, 6), dpi=300)
	m = Basemap(projection=proj,lon_0=center)
	m.drawmapboundary(fill_color=[0.82, 0.878, 0.961])
	m.drawparallels(np.arange(-90.,90.,10.))
	m.drawmeridians(np.arange(0.,360.,30.))
	for i in np.arange(0, 360, 30):
		if i != abs(180 - center):
			xpt,ypt = m(i,0)
			m.plot(xpt,ypt, marker=f'${i}$', markersize = 4 * len(str(i)), color='black')
	for i in np.arange(-80, 81, 10):
		xpt,ypt = m(0,i)
		m.plot(xpt,ypt, marker=f'${i}$', markersize = 4 * len(str(i)), color='black')
	return m
		
def PlotDataPoints(proj, center, ra, dec, clustered=False, labels=None):
	print('Plotting...')
	m = PlotMapGrid(proj, center)
	x, y = m(ra,dec)
	if clustered == False:
		prefix = 'raw'
		m.scatter(x, y, marker='.', s = 1, color='red', alpha=0.7)
	else:
		prefix = 'clustered'
		n_clusters = len(set(labels))
		for i in range(n_clusters):
			indices = np.where(np.array(labels) == i + 1)
			x, y = m(np.array(ra)[indices], np.array(dec)[indices])
			m.scatter(x, y, marker='.', s = 1)
	if not os.path.exists('./pictures'):
		os.makedirs('./pictures')
	output_filename = 'pictures/' + prefix + '_map.png'
	plt.gca().xaxis.set_inverted(True)
	plt.tight_layout()
	plt.savefig(output_filename)
	print("Plot saved to './pictures'.\n")

def Equatorial2Cartesian(ra, dec):
	ra = np.deg2rad(ra)
	dec = np.deg2rad(dec)
	print('Computing cartesian coords...')
	cartesian = [[np.cos(ra[i]) * np.cos(dec[i]), np.sin(ra[i]) * np.cos(dec[i]), np.sin(dec[i])] for i in range(len(ra))]
	print('done.\n')
	return cartesian

def ComputeDistances(cartesian, metric='cosine'):
	print('Computing distances...')
	dist = pdist(cartesian, metric=metric)
	print('done.\n')
	return dist

def Clustering(ra, dec, threshold=0.5):
	cartes_coords = Equatorial2Cartesian(ra, dec)
	dists = ComputeDistances(cartes_coords)
	print('Clustering...')
	Z = average(dists)
	print('done.\n')
	print('Assigning cluster labels to data points...')
	labels = fcluster(Z, t=threshold, criterion='distance')
	print('done.\n')
	return labels

month, n_meteors, proj, center, metric, threshold = InitializeParams()
RA, Dec = GetData(month, n_meteors)
PlotDataPoints(proj, center, RA, Dec, clustered=False)
cluster_labels = Clustering(RA, Dec, threshold)
PlotDataPoints(proj, center, RA, Dec, clustered=True, labels=cluster_labels)

