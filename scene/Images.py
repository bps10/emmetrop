import numpy as np
from matplotlib.pylab import imread
from scipy import optimize
from scipy import io as sio
import os

from scene import DataManip as dm
from database import Database as db
from analysis import Information as info

class Images(object):
    
    """Compute the amplitude spectrum of an image.
    
    Currently working on using a database of images.
    
    .. todo::
       * create a cython wrapper around c++ implementation. \
       or write own ray tracer
       * use webkits to download images from web.
       * add help to Qt
       
    """
    
    def __init__(self):
        """
        
        Creates a database instance
        
        """
        
        self.Dbase = db.Database()
        
        if os.path.isdir('./eschaton'):
            p = './eschaton/'
        else:
            p = '.'
            
        try:
            self.Dbase.OpenDatabase(p + '/ImageDatabase.h5')
        except db.DatabaseError:
            self.Dbase.CreateDatabase('ImageDatabase')
            
        self.amp = None
        self.getData()
        
    def getData(self, Directory = ['C:/Data/UPenn_Images/Images/cd01a',
                                   'C:/Data/UPenn_Images/Images/cd02a',
                                   'C:/Data/UPenn_Images/Images/cd32a',
                                   'C:/Data/UPenn_Images/Images/cd38a'],
                                   GroupName = None):
        """Find data in a database or import it if it does not exist.
         
        :param Directory: list of directories of images to analyze.
        :type Directory: list
        :param GroupName: name of groups in directory list. If None a name 
        is approximated based on directory paths
        :type GroupName: str
        
        :returns: Amplitude spectrum
        :return type: numpy.array

        This function will find all jpg files in the list of directories
        check to see if they are in the HDF5 ImageDatabase, import them if
        they are not, and return the amplitude spectrum for the group.
        
        .. todo::
           * ``high prior`` Compute Information for images as well.
           * ``low prior`` make a gui - use the database gui. Introduce images.


        """
        amp_Specs = []        
        for group in Directory:
            GroupName = group[-5:]
            print 'GroupName: ', GroupName
    
            files = dm.getAllFiles(group, suffix='.JPG')
            if self.Dbase.Exists(GroupName) == False:
                self.Dbase.CreateGroup(GroupName)
                
            
            for path in files: 
                
                img = None
                imgBW = None
                amplitude = None
                
                
                name = os.path.basename(path[:-4]) # remove .jpg

                if self.Dbase.Exists(name, GroupName) == False:
    
                    self.Dbase.CreateGroup(name, GroupName)
                      
                # subnode == Image
                if self.Dbase.Exists('Image', 
                                        GroupName + '.' + name) == False:
                                            
                    self.Dbase.CreateGroup('Image', GroupName + '/' + name)
                
                if self.Dbase.Exists('raw_image', 
                                     GroupName + '.' + name + '.Image') == False:
                        
                    img = imread(path)
                    self.Dbase.AddData2Database('raw_image', img, 
                                                GroupName + '.' + name + '.Image')
                
                if self.Dbase.Exists('grayscale', 
                                     GroupName + '.' + name + '.Image') == False:
                                            
                    if img == None:
                        img = self.Dbase.QueryDatabase(GroupName, 
                                                    name + '.' + 'Image', 
                                                    'raw_image')   
                                         
                    imgBW = self.rgb2gray(img)
                    self.Dbase.AddData2Database('grayscale', img, 
                                                   GroupName + '.' + name 
                                                   + '.Image')
                
                # metadata stuff:
                if self.Dbase.Exists('path', GroupName + '.' + name) == False:
    
                    self.Dbase.AddData2Database('path', 
                                                   np.array([path], dtype=str),
                                                   GroupName + '.' + name)
    
                if self.Dbase.Exists('goodFile', GroupName + '.' + name) == False:
                    
                    goodFile = self.getGoodFiles(path)
                    self.Dbase.AddData2Database('goodFile', np.array([goodFile],
                                                                     dtype=bool),
                                                GroupName + '.' + name)
                else:
                    goodFile = self.Dbase.QueryDatabase(GroupName, name, 
                                                        'goodFile')
                                
                                      
                # subnode == Amplitude spectrum
                if self.Dbase.Exists('amplitude', GroupName + '.' + name) == False:
                    
                    self.Dbase.CreateGroup('amplitude', GroupName + '/' + name)
                                  
                if self.Dbase.Exists('raw_amplitude', 
                                        GroupName + '.' + name
                                        + '.' + 'amplitude') == False:
                    
                    if imgBW == None:
                        imgBW = self.Dbase.QueryDatabase(GroupName, 
                                                            name + '.' 'Image', 
                                                            'grayscale')
                        
                    amplitude = self.welch2d(imgBW[500:2000, 500:2000]) 
                    self.Dbase.AddData2Database('raw_amplitude', amplitude, 
                                                   GroupName + '.' + name 
                                                   + '.' + 'amplitude')
                
                if self.Dbase.Exists('amplitude_density', 
                                        GroupName + '.' + name + 
                                        '.' + 'amplitude') == False:
                    
                    if amplitude == None:
                        amplitude = self.Dbase.QueryDatabase(GroupName, 
                                                             name + '.' + 'amplitude',
                                                             'raw_amplitude')
                    amplitude = self.Density(amplitude)                        
                    self.Dbase.AddData2Database('amplitude_density', amplitude, 
                                                   GroupName + '.' + name + '.' 
                                                   + 'amplitude')
                                                   
                    if goodFile != False or goodFile != 'F': 
                        amp_Specs.append(amplitude)
                    else:
                        pass
                # if already exists, query database to get amplitude spectrums here:
                else:
                    if goodFile:
                        amp_Specs.append(self.Dbase.QueryDatabase(GroupName, 
                                                              name + '.' + 'amplitude', 
                                                              'amplitude_density'))
                    else:
                        pass
                    
            self.Dbase.file.flush()
        self.Dbase.CloseDatabase()
    
        self.amp = np.zeros((amp_Specs[0].shape))
        total_images = len(amp_Specs)
        print 'number of images: ', total_images
        for amp in amp_Specs:
            self.amp += amp / total_images


    def computeImageHist(self,img):
        """Compute the histogram of a grayscale (for now) image
        
        :param img: input image.
        :type img: np.array
        
        this function is not complete.
        """
        #image_shape = img.shape
        #BINS = np.linspace(0,4 * 10**4, BIN_RESOLUTION + 1)
        
        pass
        #preallocate memory
        

    def getGoodFiles(self, path):
        """Find files that have warnings and skip them. Marks good files (no 
        warning) with a 1 and bad files with a 0. Usually a file has a warning
        if greater that 5% of pixels are saturated.
        
        :param path: path to .jpg file. Will remove .jpg and append '_AUX.mat' 
        :type path: dir
        
        :returns: bool indicating whether the current image has a warning.
        
        
        This function is called by getData(). The results will be stored in 
        the image database.  Only images that return True will be used in 
        analysis.
        
        """
        # First find an index of images that contain pixel saturations:
        try:
            AUX_file = sio.loadmat(path[:-4] + '_AUX.mat')
            
        except ValueError:
            print "No Auxilary files found. Check directory path or download from: \
            ftp://tofu.psych.upenn.edu/"

        if AUX_file['Image']['warning'][0][0].shape[0] > 0:
            goodFile = False
        else :
            goodFile = True

        return goodFile
    
    
    def estimateInfo(self):
        """Estimate the information in a simple linear cone receptive field.
        
        This function is under development. It will be called by getData 
        function.

        """ 
        information = {}
        information['cone'] = info.ConeInfoFunc()
        information['array'] = info.ArrayInfoFunc()
        information['optPercent'] = info.OptPercentX()
        print information           
        pass
    
    
    def Density(self, dat):
        """Compute the density of a power spectrum.
        
        :param dat: data to normalize.
        :type dat: np.array
        
        :returns: normalized data.
        
        .. note::
           Currently very simply data / sum(data)
        """
    
        return dat / np.sum(dat)
    
    
    def PowerSpectrum2(self, im,win=2,n1=1,n2=0):
        """2D spectrum estimation using the modified periodogram.
        This one includes a window function to decrease variance in \
        the estimate.
        
        :param x: input sequence
        :type x: numpy.array
        :param n1: starting index, x(n1)
        :type n1: int
        :param n2: ending index, x(n2)
        :type n2: int
        :param win: The window type \n
                    1 = Rectangular \n
                    2 = Hamming \n
                    3 = Hanning \n
                    4 = Bartlett \n
                    5 = Blackman \n
        :type win: int
        
        :returns: spectrum estimate.
        :rtype: numpy.array
                
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
        
        fftim = np.abs(np.fft.fftshift(np.fft.fft2(((im) * m)))) / ( (N**2.0) * U)
        
        return fftim
    
    
    
    def welch2d(self, x,L = None, over = 0.5, win = 2.0):
        """2D spectrum estimation using Welch's method.
        The spectrum of a process x is estimated using Welch's \
        method of averaging modified periodograms.
        
        :param x: input sequence
        :param L: section length
        :param over: amount of overlap, where 0<over<1,
        :param win: The window type \n
                        1 = Rectangular \n
                        2 = Hamming \n
                        3 = Hanning \n
                        4 = Bartlett \n
                        5 = Blackman \n
        
        
        :returns: Welch's estimate of the power spectrum, returned in decibels. 
        
        .. note:: 
           Modified from: M.H. Hayes. \
           "Statistical Digital Signal Processing and Modeling" \
           (John Wiley & Sons, 1996).
            
        """
        
        xs, ys = x.shape
        
        if L == None:
            L = xs
        
        
        if xs / ys != 1.0:
            raise ValueError('This is a stupid program. Dimensions need to be \
            equal (len(x)=len(y))')
        
        
        
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
                Px += self.PowerSpectrum2(x[ n1[0]:n2[0],n1[1]:n2[1] ],
                                          win) / (nsect**2)
        
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
        avg = dm.accum(r.flatten() + 1, Px.flatten()
                        )[1:] / dm.accum(r.flatten() + 1)[1:]
        avg = avg[:(xs / 2 + 1)]
        
        
        return avg
    
    
    def PowerLaw(self, xdata, ydata):
        """Create an array according to a power law for plotting.
        
        :param xdata: x values at which to fit a power law.
        :type xdata: numpy.array
        :param ydata: y values used to fit power law.
        :type ydata: numpy.array
        
        :returns: creates a handle to the lambda function to compute a \
        power law fit to the inputed data.
        :rtype: function handle
        
        .. note::   * Previously used only for drawing a power law, \
                        no fitting.
                    * Uncomment code below to reintroduce that functionality.
                    * Fitting code taken from `Scipy`_
                    
                    * See above reference to introduce error bar to fit.
                    
        .. _Scipy: http://www.scipy.org/Cookbook/FittingData
        
        """
        
        logx = np.log10(xdata)
        logy = np.log10(ydata)
        #logyerr = yerr / ydata
        
        # define our (line) fitting function
        fitfunc = lambda p, x: p[0] + p[1] * x
        errfunc = lambda p, x, y: (y - fitfunc(p, x))
        
        pinit = [1.0, -1.0]
        out = optimize.leastsq(errfunc, pinit,
                               args=(logx, logy), full_output=1)
        pfinal = out[0]
        
        #print pfinal
        #print covar
        
        index = pfinal[1]
        amp = 10.0**pfinal[0]
        
        self.powerlaw = lambda x,: amp * (x**index)
    
    
    def rgb2gray(self, rgb):
        """Convert an image from rbg into gray scale.
        
        :param rgb: RGB image to be converted.
        
        :returns: Grayscale image of the same size as RGB input image.
        
        *Stole this one from stackoverflow.*
        
        Formula used in conversion:
        .. math::
           0.2999*R_{channel} + 0.587*G_{channel} + 0.114*B_{channel}
            
        """
        
        r, g, b = np.rollaxis(rgb[...,:3], axis = -1)
        return 0.299 * r + 0.587 * g + 0.114 * b


