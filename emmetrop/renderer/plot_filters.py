import matplotlib.pylab as plt
import numpy as np

from base import plot as pf

from emmetrop.eye.movement import brownian_motion
from emmetrop.eye.filter import gauss

def plotMovement():
	'''
	'''
	temp = np.arange(1, 80)
	spat = np.arange(0.1, 100, 0.1)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	pf.AxisFormat()
	pf.TufteAxis(ax, ['bottom', 'left'])
	ax.loglog(spat, brownian_motion(spat, temp), 'k')

	ax.set_ylim([10 ** -5, 1.05])
	ax.set_ylabel('transfer')
	ax.set_xlabel('spatial frequency (cycles / deg)')
	plt.tight_layout()
	plt.show()

def plotGlasses():
	'''
	'''
	spat = np.arange(0.1, 100, 0.1)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	pf.AxisFormat()
	pf.TufteAxis(ax, ['bottom', 'left'])
	ax.loglog(spat, gauss(spat, 15), 'k')

	ax.set_ylim([10 ** -5, 1.05])
	ax.set_ylabel('transfer')
	ax.set_xlabel('spatial frequency (cycles / deg)')
	plt.tight_layout()
	plt.show()