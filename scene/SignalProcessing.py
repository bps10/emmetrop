import numpy as np
import DataManip as dm


## To work on:
## 1. welch3D.

def FindPowerExponents(spectrum):
	"""Find the exponent,m, that provides a best fit for to data using 1/f^m for
	a 2D spectrum.
	
	:param spectrum: power spectrum.
	
	
	:returns:  * M_X: values that are best fits to data in X columns
                 * M_Y: values that are best fits to data in Y columns
	
	.. note:: 
         not working in python yet. See Cone Activity
 
	"""

	xs, ys, = spectrum.shape

	M_Y = np.zeros((ys))
	M_X = np.zeros((xs))

	for i in range(0,xs):
		M_X[i] = FitPowerLaw(spectrum[i,:],2)


	for i in range(0,yx):
		M_Y[i] = FitPowerLaw(np.flipud(spectrum[:,i]),2)
		
	return M_Y

		

def RemoveMean(image):
	"""
	RemoveMean(image)
	
	subtract the mean of an image.
	
	:param image: input image.
	
	:returns: image with mean subtracted

	"""

	x, y = image.shape
	mu = np.sum(image) / (x * y)
	out = image - mu
	
	return out
	
def SpectrumDensity(rawSpectrum):
	"""find the spectrum density (sums to 1), then find log10 to convert
	into decibels.
	
	:param rawSpectrum: Power Spectrum

	:returns: log 10 spectrum normalized to sum to 1. (spectral density in 
                decibels.
	
	
	.. note::
         Called by Welch2d.m and Welch3d.m.
	
	"""
	
	normSpec = rawSpectrum / np.sum(rawSpectrum)
	decibel = np.log10(normSpec)
	
	return decibel

def PowerLaw(freq,power,max_power):
	"""Create an array according to a power law for plotting.
	
	:param freq: frequencies to be plotted on x axis
	:param power: power exponent to use
	:param max_power: maximum power (y axis).

	:returns: array of values according to 1/freq^power . 

	
	"""
	norm_freq = (1.0 / freq**power) / np.max(1.0 / freq**power)
	out = norm_freq * max_power
	return out

def PowerSpectrum2(im,win=2,n1=1,n2=0):
	"""2D spectrum estimation using the modified periodogram.
	This one includes a window function to decrease variance in the estimate.
	
	:param x: input sequence
	:param n1: starting index, x(n1)
	:param n2: ending index, x(n2)
	:param win: The window type
			1 = Rectangular
			2 = Hamming
			3 = Hanning
			4 = Bartlett
			5 = Blackman
	
	:return: spectrum estimate.
	
	.. note::
         If n1 and n2 are not specified the periodogram of the entire
         sequence is computed.
     
	"""
	if n2 == 0:
		n2 = len(im[:,1])
	
	N  = n2 - n1 + 1
	w  = np.ones((N))

	if (win == 2):
		w = np.hamming(N)
	elif (win == 3):
	   w = np.hanning(N)
	elif (win == 4):
	   w = np.bartlett(N)
	elif (win == 5): 
	   w = np.blackman(N);
	
	   

	xs, ys = im.shape
	if xs/ys != 1:
		raise ValueError('Dimensions must be equal')
	

	m = w[:] * w[:][np.newaxis,:]

	U  = np.linalg.norm(w)**2.0 / N**2.0

	fftim = np.abs(np.fft.fftshift(np.fft.fft2(((im) * m))))**2.0 / ( (N**2.0) * U)
	
	return fftim
	
def PowerSpectrum3(im,win=2,n1=1,n2=0):
	"""2D spectrum estimation using the modified periodogram.
	This one includes a window function to decrease variance in the estimate.
	
	:param x: input sequence
	:param n1: starting index, x(n1)
	:param n2: ending index, x(n2)
	:param win: The window type
			1 = Rectangular
			2 = Hamming
			3 = Hanning
			4 = Bartlett
			5 = Blackman
	
	:return: spectrum estimate.
	
	.. note::
         If n1 and n2 are not specified the periodogram of the entire
         sequence is computed.
         Not finished.  See PowerSpectrum3 to develop.
         
	"""
	if n2 == 0:
		n2 = len(im[:,1])
		
	N  = n2 - n1 + 1
	w  = np.ones((N))

	if (win == 2):
		w = np.hamming(N)
	elif (win == 3):
	   w = np.hanning(N)
	elif (win == 4):
	   w = np.bartlett(N)
	elif (win == 5): 
	   w = np.blackman(N)
	
	   

	xs, ys, zs = im.shape
	if xs/ys != 1.0 or xs/zs != 1.0:
		raise ValueError('Dimensions must be equal')
	
	m = np.ones((len(w),len(w),len(w)))
	foo = w[:] * w[:][np.newaxis,:]
	for i in range(0,len(w)):
		m[:,:,i] = foo
	
	U  = np.linalg.norm(w)**3.0 / N**3.0

	fftim = np.abs(np.fft.fftshift(np.fft.fftn(float(im) * m)))**2.0 / ( (N**3.0) * U)
	
	return fftim



def welch2d(x,L = None, over = 0.5, win = 2.0, Density = 1):
	"""2D spectrum estimation using Welch's method.
	The spectrum of a process x is estimated using Welch's method of 
      averaging modified periodograms.

	:param x: input sequence
	:param L: section length 
	:param over: amount of overlap, where 0<over<1, 
	:param win : The window type \n
			1 = Rectangular \n
			2 = Hamming \n
			3 = Hanning \n
			4 = Bartlett \n
			5 = Blackman \n
	
	:param Density: default returns estimate in spectral density measured in
	decibels. Passing a 1 here will turn this off and return the raw
	estimate.
	
	:returns: Welch's estimate of the power spectrum, returned in decibels. 

	.. note::
         Modified from: M.H. Hayes. "Statistical Digital Signal Processing \
         and Modeling" (John Wiley & Sons, 1996).
    
	"""
	xs, ys = x.shape
	
	if L == None:
		L = xs
		
	
	if xs / ys != 1.0:
		raise ValueError('This is a stupid program. Dimensions need to be equal (len(x)=len(y))')
	


	if L < len(x[:,0]) / 2.0:
		raise ValueError('Length must be longer than 1/2 length of x')
	

	if (over >= 1) or (over < 0):
		raise ValueError('Overlap is invalid')


	n0 = (1.0 - over) * L
	n1 = np.array([1.0, 1.0]) - n0
	n2 = np.array([L, L]) - n0
	nsect = int(1.0 + np.floor((len(x) - L) /( n0)))

	Px = 0
	for ix in range(0,nsect):
		n1[0] = n1[0] + n0
		n2[0] = n2[0] + n0
		for iy in range(0,nsect):
			n1[1] = n1[1] + n0
			n2[1] = n2[1] + n0
			Px += PowerSpectrum2(x[ n1[0]:n2[0],n1[1]:n2[1] ],win) / (nsect**2)

	xs, ys = Px.shape

	f2 = np.arange(-xs / 2.0, xs / 2.0)
	f1 = np.arange(-ys / 2.0, ys / 2.0)
	XX, YY = np.meshgrid(f1,f2)
	foo, r = dm.cart2pol(XX,YY)
	if np.mod(xs,2)==1 or np.mod(ys,2)==1:
		r = np.around(r)-1
	else:
		r = np.around(r)

	r = np.array(r,dtype=int)
	### need to possibly use a for loop.
	avg = dm.accum(r.flatten() + 1, Px.flatten())[1:] / dm.accum(r.flatten() + 1)[1:]
	avg = avg[:(xs / 2 + 1)]

	if Density == 1:
		avg = SpectrumDensity(avg)
		
	return avg




### Not yet finished, look at welch2d

def welch3d(x, L_xyz = 0, over = 0.5, win = 2.0, Density = 0):
	"""3D spectrum estimation using Welch's method.
	The spectrum of a process x is estimated using Welch's method
      of averaging modified periodograms.
	
	:param x: input sequence
	:param L: section length 
	:param over: amount of overlap, where 0<over<1, 
	:param win : The window type \n
			1 = Rectangular \n
			2 = Hamming \n
			3 = Hanning \n 
			4 = Bartlett \n 
			5 = Blackman \n
	
	:param Density: default returns estimate in spectral density measured in
	decibels. Passing a 1 here will turn this off and return the raw
	estimate.
	
	:returns: Welch's estimate of the power spectrum, returned in decibels. 

	.. note::
         Modified from: M.H. Hayes. "Statistical Digital Signal Processing \
         and Modeling" (John Wiley & Sons, 1996).
    
	"""
	if L_xyz == 0:
		L_xyz = len(x[:,1])


	if L_xyz < len(x[:,0]) / 2.0:
		raise ValueError('Length must be longer than 1/2 length of x')
	

	if (over >= 1) or (over < 0):
		raise ValueError('Overlap is invalid')


	n0 = (1.0 - over) * L
	n1 = [1.0, 1.0, 1.0] - n0
	n2 = [L, L, L] - n0
	nsect = 1.0 + np.floor((len(x) - L) /( n0))

	Px = 0
	for ix in range(0,nsect):
		n1[0] = n1[0] + n0
		n2[0] = n2[0] + n0
		for iy in range(0,nsect):
			n1[1] = n1[1] + n0
			n2[1] = n2[1] + n0
			for iz in range(0,nsect):
				n1[2] = n1[2] + n0
				n2[2] = n2[2] + n0
				Px = Px + PowerSpectrum2(x[ n1[0]:n2[0],n1[1]:n2[1], n1[2]:n2[2] ],win) / (nsect**2)

	xs, ys, zs, = Px.size
	
	W_s = np.zeros(np.floor(zs), np.floor(xs / 2))
	for i in range(0, np.floor(zs)):
		fft = Px[:,:,i]
		f2 = np.arange(-xs / 2.0, xs / 2.0 - 1.0)
		f1 = np.arange(-ys / 2.0, ys / 2.0 - 1.0)
		XX, YY = np.meshgrid(f1,f2)
		foo, r = dm.cart2pol(XX,YY)
		if mod(xs,2)==1 or mod(ys,2)==1:
			r = round(r)-1
		else:
			r = round(r)


		avg = dm.accum(r.flatten() + 1, Px.flatten())[1:] / dm.accum(r.flatten() + 1)[1:]
		avg = avg[:(xs / 2 + 1)]

	if Density == 0:
		avg = SpectrumDensity(avg)
		
	return avg
	
	
		
