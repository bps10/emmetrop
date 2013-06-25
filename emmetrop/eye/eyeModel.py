import os, sys
import shlex
import numpy as np

if os.path.basename(os.getcwd()) == 'emmetrop':
    sys.path.append("./SchematicEye")
else:
    sys.path.append("./emmetrop/SchematicEye")
import eye as eye


def rad2deg(radians):
    """convert radians to degrees
    """
    return radians * 180.0 / np.pi


class SchematicEye(object):
    """
    This class controls the loading of OSLO data.
    Eventually it will interface with a C++ ray tracer dynamically.
    """
    
    def __init__(self, Use_OSLO=False, OSLO_directory=None):
        """
        """
        if Use_OSLO:
            if not OSLO_directory:
                if os.path.basename(os.getcwd()) == 'emmetrop':
                    p = './OSLO_MTF_DATA/'
                else:
                    p = './emmetrop/OSLO_MTF_DATA/'
            else:
                p = OSLO_directory
            self.loadOSLOData(p)


    def traceEye(self, object_distance=100000, off_axis=0, pupil_size=4, diopters=0):
        '''
        '''

        intensity = eye.py_eye(object_distance, off_axis, pupil_size, diopters)

        return intensity

    def _genPSF(self, intensity, samples):
        '''

        .. todo::
           Need to unhard code xvals and samples (work with C++ - return stuct?)
        '''
        xvals = np.arange(0.0005, 0.2001, 0.0005)
        PSF = np.zeros(samples)
        PSFtotal = np.zeros(samples * 2)

        # we have an integral, therefore take the deriv to get rays / bin
        deriv = np.zeros((samples))
        deriv[0] = intensity[0]
        deriv[1:] = intensity[1:] - intensity[0:-1]

        for i in range(0, samples - 1):

            # account for increasing size of area
            radius0 = xvals[i]
            radius1 = xvals[i + 1]

            # subtract inner and outer circle area to get sliver of interest
            area = (np.pi * radius1 ** 2.0) - (np.pi * radius0 ** 2.0)

            # deriv = amount in each circle; then divide by area
            PSF[i] = deriv[i]  / area 

        # normalize so that each PSF has same integral of 1.

        PSF = PSF / np.sum(PSF)

        PSFtotal[1:samples + 1] = PSF[::-1]
        PSFtotal[samples + 1:] = PSF[1:]


        return PSFtotal

    def genMTF(self, intensity):
        '''
        '''
        samples = len(intensity)
        PSF = self._genPSF(intensity, samples)

        MTF = np.zeros(samples)
        # do the FFT, take only right half
        temp = np.abs(np.fft.fftshift(np.fft.fft(PSF)))
        temp = temp[samples:]

        # make sure we only get real part
        temp = np.real(temp)

        # normalize MTF
        MTF = temp / np.max(temp)
        # temporary for now
        freqs = np.arange(0, len(MTF))
        return MTF, freqs

    def getAxialLength(self):
        """Find the axial length in mm of the optical system used to generate 
        outputed curves.
        
        .. warning::
           Currently hard coded to return 24mm
           
        .. todo::
           Find axial length dynamically. Call getAxialLength function.
        """
        return 24.0

    ### DEPRECIATED OSLO FUNCTIONS ###

    def loadOSLOData(self,p):
        """Load data from OSLO
        
        This is currently a very static function.
        
        :param p: directory path. Passed from init function.
        :type p: path
        
        """
        
        self.INF = self.importOSLOfile(p + 
                                'ONaxisMTFinfFocusNavarrow1999.txt')
        self.TwentyFt = self.importOSLOfile(p +
                                'ONaxisMTF20ftFocusNavarrow1999.txt')
        self.Onemeter = self.importOSLOfile(p + 
                                'ONaxisMTF1mFocusNavarrow1999.txt')
        self.SixteenIn = self.importOSLOfile(p + 
                                'ONaxisMTF16inFocusNavarrow1999.txt')
        
        self.INF_offaxis = self.importOSLOfile(p + 
                                'OFFaxisMTFinfFocusNavarrow1999.txt')
        self.TwentyFt_offaxis = self.importOSLOfile(p + 
                                'OFFaxisMTF20ftFocusNavarrow1999.txt')
        self.Onemeter_offaxis = self.importOSLOfile(p +
                                'OFFaxisMTF1mFocusNavarrow1999.txt')
        self.SixteenIn_offaxis = self.importOSLOfile(p + 
                                'OFFaxisMTF16inFocusNavarrow1999.txt')
        
        self.Sixteen_UnderAccomm = self.importOSLOfile(p +
                            'OFFaxisMTF16inUnderAccom20ftObjNavarrow1999.txt')
        self.Sixteen_SixteenObj_Offaxis = self.importOSLOfile(p +
                                'OFFaxisMTF16inFocus16inObjNavarrow1999.txt')        
        self.Sixteen_TwentyObj_Offaxis = self.importOSLOfile(p + 
                                'OFFaxisMTF16inFocus20ftObjNavarrow1999.txt')
    
        self.TwentyDegOffAxis_InfFoc = self.importOSLOfile(p + 
                                '20degOFFaxisMTFinfFocusNavarrow1999.txt')
        self.FortyDegOffAxis_InfFoc = self.importOSLOfile(p + 
                                '40degOFFaxisMTFinfFocusNavarrow1999.txt')       
 
        ## convert mm to rad (1mm image/24mm axial length)
        self.freqs = self.INF[:,1] / rad2deg(1.0/self.getAxialLength())
        
        
        self.dataPackage = {
                            'onAxis': {'diffract': self.INF[:,4],
                                       'inf': self.INF[:,2], 
                                       '20ft': self.TwentyFt[:,2],
                                       '1m': self.Onemeter[:,2],
                                       '16in': self.SixteenIn[:,2]
                                       },
                                   
                            'offAxis': {'diffract': self.INF_offaxis[:,4],
                                        'inf':self.INF_offaxis[:,2],
                                        '20ft': self.TwentyFt_offaxis[:,2],
                                        '1m': self.TwentyFt_offaxis[:,2],
                                        '16in': self.SixteenIn_offaxis[:,2]
                                        },
                                        
                            'object': {
                            '16inunder': self.Sixteen_UnderAccomm[:,2],
                            '16in16in': self.Sixteen_SixteenObj_Offaxis[:,2],
                            '16in20ft': self.Sixteen_TwentyObj_Offaxis[:,2]
                            },
                                       
                            'farPeriph': {
                                    '20deg': self.TwentyDegOffAxis_InfFoc[:,2],
                                    '40deg': self.FortyDegOffAxis_InfFoc[:,2]
                                    },
                                       
                            'freqs': self.freqs
                            }


    def returnOSLOdata(self):
        """
        Return a dictionary of imported transfer functions from OSLO.
        """
        return self.dataPackage
        

    
    def importOSLOfile(self, OSLOfile):
        """Import a text file with MTF output from OSLO.
        
        :param OSLOfile: name of OSLO output text file to import. \
        Should have 5 columns.
        :type OSLOfile: string
        :returns: an array containing OSLO data: Frequencies, \
        Eye MTF, Diffraction limit MTF.
        :rtype: numpy.array
        
        .. todo::
           import OSLO files, import whole directory into a dictionary.
        
        """
        
        fil = open(OSLOfile)
        fil = fil.read()
        
        MTF = []
        foo = fil
        partitioning = True
        row = 0
        while partitioning:
            
            f = foo.partition('\n')
    
            if f[0]:
                
                parse = shlex.split(f[0])
                floatnum = []
                for num in parse:
                    if row == 0:
                        if num == '--':
                            floatnum.append( 0.0 )
                    if num != '--':
                        floatnum.append( float(num) )
    
                        
                MTF.append(floatnum)
                foo = f[2]
                partitioning = True
            else:
                partitioning = False
        
        MTF = np.array(MTF)
        
        return MTF        

def diffraction(deg, samples, pupil_size_mm, focal_len, ref_index=1.55, wavelength=550.0):
    '''See Appendix B of "Light, the Retinal Image and Photoreceptors"
    Packer & Williams.

    '''
    lam = wavelength / 1000 # convert mm into meters.
    NA = NumericalAperature(ref_index, D=pupil_size_mm, focal_len=focal_len)

    s_0 = NA / lam # convert to radians
    s =  np.linspace(0, s_0, samples)
    print "NA: ", NA, "s_0", s_0

    dif = (2.0 / np.pi) * (np.arccos(s / s_0) - 
            (s / s_0) * np.sqrt(1.0 - (s / s_0) ** 2.0))

    return dif, s

def NumericalAperature(n, theta=None, D=None, focal_len=None):
    '''
    Find the numerical aperature of a system

    :param n: refractive index
    :param theta: angle of marginal rays

    According to the formula

    $$NA = n \\sin(\\theta)$$

    or 

    $$NA = n \\sin(\\arctan(\\frac{D}{2f}))$$
       
    '''
    if D is None and focal_len is None and theta is not None:
        out = n * np.sin(theta)
    elif theta is None and D is not None and focal_len is not None:
        out = n * np.sin(np.arctan(D / (2 * focal_len)))
    else:
        raise IOError("check parameters.")

    return out
