#! /usr/bin/env python
from __future__ import division
import os
import matplotlib.pylab as plt
from matplotlib.mlab import griddata
import numpy as np

from base import plot as pf


def big_analysis_plot():
	'''
	'''
	handle = os.path.join(os.path.dirname(os.path.abspath(__file__)),
		'../../analysis.txt')

	data = np.genfromtxt(handle, delimiter='\t', skip_header=5, names=True)
	#print data.dtype.names

	fig = plt.figure(figsize=(14,5))
	pf.AxisFormat(TickDirection='in')
	for i, diop in enumerate(np.unique(data['focus_D'])):
		ind = np.where(data['focus_D'] == diop)
		if diop <= 3: #diop - round(diop, 0) == 0:

			x = data['wavelen_nm'][ind]
			y = np.log10(data['obj_dist_mm'][ind])
			z = data['proportion'][ind]

			Nx = len(np.unique(x))
			Ny = len(np.unique(y))
			extent = (np.min(x), np.max(x), np.min(y), np.max(y))

			#xs, ys = np.mgrid[extent[0]:extent[1], extent[2]:extent[3]]
			#print xs
			
			#resampled = griddata(x, y, z, xs, ys)
			ax = fig.add_subplot(1, 7, i + 1)
			im = ax.imshow(z.reshape(Ny, Nx), extent = extent, aspect='auto',
				interpolation='nearest', origin='lower', vmin=0, vmax=1)
			ax.set_xticks([400, 700])
			ax.set_title('lens ' + str(diop) + 'D', fontsize=20)
			
			if i == 0:
				ax.set_ylabel('object distance log10(mm)')
				
			if i > 0:
				ax.set_yticks([])

	plt.colorbar(im)
	plt.tight_layout()
	plt.show()