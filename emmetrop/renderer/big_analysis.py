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

	fig = plt.figure(figsize=(14,8))
	pf.AxisFormat(TickDirection='in')
	for i, diop in enumerate(np.unique(data['focus_D'])):
		ind = np.where(data['focus_D'] == diop)
		if diop <= 3.5: #diop - round(diop, 0) == 0:

			x = data['wavelen_nm'][ind]
			y = np.log10(data['obj_dist_mm'][ind] / 1000.0)
			z = data['proportion'][ind]

			Nx = len(np.unique(x))
			Ny = len(np.unique(y))
			extent = (np.min(x), np.max(x), np.min(y), np.max(y))

			ax = fig.add_subplot(2, 4, i + 1)
			im = ax.imshow(z.reshape(Ny, Nx), extent = extent, aspect='auto',
				interpolation='nearest', origin='lower', vmin=0, vmax=1)

			ax.set_title('lens ' + str(diop) + 'D', fontsize=20)

			if i == 0 or i == 4:
				ax.set_ylabel('obj. distance log10(m)')
			else:
				ax.set_yticklabels(([]))

			if i > 3:
				ax.set_xlabel('wavelength (nm)')
				ax.set_xticks([400, 600, 800])
			else:
				ax.set_xticks([400, 600, 800])
				ax.set_xticklabels(([]))

	#plt.colorbar(im)
	plt.subplots_adjust(bottom=0.01, hspace=0.01, wspace=0.01)
	plt.tight_layout()
	plt.show()