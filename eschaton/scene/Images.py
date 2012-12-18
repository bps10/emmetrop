import numpy as np
from matplotlib.pylab import imread
from scipy import io as sio
from scipy import optimize
import os, sys

from eschaton.scene import SignalProcessing as sig
from eschaton.scene import DataManip as dm
from eschaton.database import Database as db


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
        
        try:
            self.Dbase.OpenDatabase('./eschaton/ImageDatabase.h5')
        except db.DatabaseError:
            self.Dbase.CreateDatabase('ImageDatabase')
            print 'created new image database'
                
        self.amp_mean = None
        self.getData()
        
        self.imagexval = np.arange(1,self.amp_mean.shape[0] + 1) / 46.0
        self.PowerLaw(self.imagexval, self.amp_mean)

    def returnImageData(self):
        """
        """
        
        imageData =     {
                        'totalImages': len(self.ampSpecs),
                        'ampSpecs': self.ampSpecs,
                        'ampMean': self.amp_mean,
                        'rawAmp': self.rawAmp,
                        'powerlaw': self.powerlaw,
                        'imagexval': self.imagexval
                        }
                        
        return imageData
               
    def getData(self, Directory = None, GroupName = None):
        """Find data in a database or import it if it does not exist.
         
        :param Directory: list of directories of images to analyze.
        :type Directory: list
        :param GroupName: name of groups in directory list. If None a name \
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
        if not Directory:
            if sys.platform == 'darwin':
                Directory = ['/Users/brianschmidt/Documents/cd01A/']
                index = [-6, -1]
            if sys.platform == 'win32':
                Directory = ['C:/Data/UPenn_Images/Images/cd01A',
                             'C:/Data/UPenn_Images/Images/cd02A',
                             'C:/Data/UPenn_Images/Images/cd32A',
                             'C:/Data/UPenn_Images/Images/cd38A',
                             'C:/Data/UPenn_Images/Images/cd41A',
                             'C:/Data/UPenn_Images/Images/cd58A']
                index = [-5, len(Directory[0])]
    
        self.ampSpecs = []
        self.rawAmp= []
        for group in Directory:
            GroupName = group[index[0]:index[1]]
            print 'GroupName: ', GroupName
    
            files = dm.getAllFiles(group, suffix='.JPG', subdirectories=1)
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
                                         
                    imgBW = sig.rgb2gray(img)
                    self.Dbase.AddData2Database('grayscale', img, 
                                                   GroupName + '.' + name 
                                                   + '.Image')
                
                # metadata stuff:
                if self.Dbase.Exists('path', GroupName + '.' + name) == False:
    
                    self.Dbase.AddData2Database('path', 
                                                   np.array([path], dtype=str),
                                                   GroupName + '.' + name)
    
                if self.Dbase.Exists('goodFile', 
                                     GroupName + '.' + name) == False:
                    
                    goodFile = self.getGoodFiles(path)
                    self.Dbase.AddData2Database('goodFile',
                                                np.array([goodFile], 
                                                         dtype=bool),
                                                GroupName + '.' + name)
                else:
                    goodFile = self.Dbase.QueryDatabase(GroupName, name, 
                                                        'goodFile')
                                
                                      
                # subnode == Amplitude spectrum
                if self.Dbase.Exists('amplitude', 
                                     GroupName + '.' + name) == False:
                    
                    self.Dbase.CreateGroup('amplitude', GroupName + '/' + name)
                 
                # raw amplitude spectrum first                 
                if self.Dbase.Exists('raw_amplitude', 
                                        GroupName + '.' + name
                                        + '.' + 'amplitude') == False:
                    
                    if imgBW == None:
                        imgBW = self.Dbase.QueryDatabase(GroupName, 
                                                            name + '.' 'Image', 
                                                            'grayscale')
                        
                    amplitude = sig.welch2d(imgBW[500:2000, 500:2000]) 
                    self.Dbase.AddData2Database('raw_amplitude', amplitude, 
                                                   GroupName + '.' + name 
                                                   + '.' + 'amplitude')
                    if goodFile != False or goodFile != 'F': 
                        self.rawAmp.append(amplitude)
                    else:
                        pass
                else:
                    if goodFile != False or goodFile != 'F': 
                        self.rawAmp.append(self.Dbase.QueryDatabase(GroupName, 
                                                name + '.' + 'amplitude', 
                                                              'raw_amplitude'))
                    else:
                        pass                    
                
                # then amplitude density
                if self.Dbase.Exists('amplitude_density', 
                                        GroupName + '.' + name + 
                                        '.' + 'amplitude') == False:
                    
                    if amplitude == None:
                        amplitude = self.Dbase.QueryDatabase(GroupName, 
                                                    name + '.' + 'amplitude',
                                                             'raw_amplitude')
                    amplitude = sig.Density(amplitude)                        
                    self.Dbase.AddData2Database('amplitude_density', amplitude, 
                                                   GroupName + '.' + name + '.' 
                                                   + 'amplitude')
                                                   
                    if goodFile != False or goodFile != 'F': 
                        self.ampSpecs.append(amplitude)
                    else:
                        pass
                # if already exists, query database to get amplitude spectrums here:
                else:
                    if goodFile:
                        self.ampSpecs.append(
                                            self.Dbase.QueryDatabase(GroupName, 
                                                    name + '.' + 'amplitude', 
                                                    'amplitude_density'))
                    else:
                        pass
                    
            self.Dbase.file.flush()
        self.Dbase.CloseDatabase()
    
        self.amp_mean = np.zeros((self.ampSpecs[0].shape))
        total_images = len(self.ampSpecs)
        print 'number of images: ', total_images
        for amp in self.ampSpecs:
            self.amp_mean += amp / total_images

    

        
        
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
            print "No Auxilary files found. Check directory path or download \
            from: ftp://tofu.psych.upenn.edu/"

        if AUX_file['Image']['warning'][0][0].shape[0] > 0:
            goodFile = False
        else :
            goodFile = True

        return goodFile


    
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
    
 

