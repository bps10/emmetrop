import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
import shlex
import os

from scene.Images import Images
from analysis import PlottingFun as pf
from analysis import Information as info

class SchematicEyeAnalysis(Images):
    """This class is designed to estimate the response of a linear photoreceptor 
    to an amplitude spectrum derived from natural images.
    The modulation transfer functions are computed in OSLO ray tracing software.
    
    **Legend for figures:**
    
    .. figure:: ../../Figures/Legend.png
       :height: 200px
       :width: 300px
       :align: center
       
       **Fig 1: Use this legend for all figures in this section**
       
    .. todo::
       * single graphic with all subfigures?
       * better plotting options.
       * Write OSLO data into database.
       * self.xval needs to be computed during each iteration. 
         (__getAxialLength__)
    """
    
    def __init__(self):
        """
        .. todo::
           more dynamic OSLO import.
        
        """

        if os.path.isdir('./eschaton'):
            p = './eschaton/OSLO_MTF_DATA/'
        else:
            p = './OSLO_MTF_DATA/'
            
        self.INF = self.importOSLOfile(p + 'ONaxisMTFinfFocusNavarrow1999.txt')
        self.TwentyFt = self.importOSLOfile(p + 'ONaxisMTF20ftFocusNavarrow1999.txt')
        self.Onemeter = self.importOSLOfile(p + 'ONaxisMTF1mFocusNavarrow1999.txt')
        self.SixteenIn = self.importOSLOfile(p + 'ONaxisMTF16inFocusNavarrow1999.txt')
        
        self.INF_offaxis = self.importOSLOfile(p + 'OFFaxisMTFinfFocusNavarrow1999.txt')
        self.TwentyFt_offaxis = self.importOSLOfile(p + 'OFFaxisMTF20ftFocusNavarrow1999.txt')
        self.Onemeter_offaxis = self.importOSLOfile(p + 'OFFaxisMTF1mFocusNavarrow1999.txt')
        self.SixteenIn_offaxis = self.importOSLOfile(p + 'OFFaxisMTF16inFocusNavarrow1999.txt')
        
        self.Sixteen_UnderAccomm = self.importOSLOfile(p + 'OFFaxisMTF16inUnderAccom20ftObjNavarrow1999.txt')
        self.SixteenFocus_SixteenObj_Offaxis = self.importOSLOfile(p + 'OFFaxisMTF16inFocus16inObjNavarrow1999.txt')        
        self.SixteenFocus_TwentyObj_Offaxis = self.importOSLOfile(p + 'OFFaxisMTF16inFocus20ftObjNavarrow1999.txt')
    
        self.TwentyDegOffAxis_InfFoc = self.importOSLOfile(p + '20degOFFaxisMTFinfFocusNavarrow1999.txt')
        self.FortyDegOffAxis_InfFoc = self.importOSLOfile(p + '40degOFFaxisMTFinfFocusNavarrow1999.txt')       
        
        ## convert mm to deg (1mm image/24mm axial length)
        self.xval = self.INF[:,1] / 2.38732415 
        
        Images.__init__(self)
        
        self.powerlaw = None
        self.RF_SPLINE = None
        self.RF_DOG = None
        self.Jay_RF = None
        
        self.DiffractionLim = {}
        self.Infinity = {}
        self.NearFocusFarObject = {}
        self.NearFocusNearObject = {}
        self.UnderAccommFarObject = {}


    ###  DATA PROCESSING FUNCTIONS ###
    ##################################
    
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

        

    def Find_PowerLaw_Placement(self):
        """Find where the :math:`\\frac{1}{f}` power law should be placed.
        
        .. note:: not yet working. Should find len(array)/2. Not much else.
        
        .. todo::
           1 / f heuristic
           
        """
        pass



    def ComputeConeActivity(self, Receptive_Field):
        """Compute the estimated activity of a cone photoreceptor.
        
        :param Receptive_Field: a handle to the spline fitted receptive field.
        :type Receptive_Field: function handle.
        
        .. note:: 
           very preliminary. Completely hard coded.
        
        """
        if Receptive_Field == None:
            raise('Run {0} receptive field function'.format(Receptive_Field))
            
        
        imagexval = np.arange(1,self.amp_mean.shape[0] + 1) / 46.0
        if not self.powerlaw:
            self.PowerLaw(imagexval, self.amp_mean)
            
        powerlaw = self.powerlaw(self.xval[1:])

        RField = interpolate.splev(self.xval[1:], Receptive_Field, der = 0)
        RField = RField / sum(RField)
        
                
        self.DiffractionLim['name'] = 'diffraction limit'
        self.DiffractionLim['yval'] = powerlaw * RField
        self.DiffractionLim['xval'] = self.xval[1:]
        
        self.Infinity['name'] = 'inf focus'
        self.Infinity['yval'] = powerlaw * self.INF[1:, 2] * RField
        self.Infinity['xval'] = self.xval[1:]
        

        self.NearFocusNearObject['name'] = 'near focus, near object'
        self.NearFocusNearObject['yval'] = (powerlaw[:59] *
                                    self.SixteenFocus_SixteenObj_Offaxis[1:60,2] 
                                    * RField[:59])
        self.NearFocusNearObject['xval'] = self.xval[1:60]                
        

        self.NearFocusFarObject['name'] = 'near focus, far object'
        self.NearFocusFarObject['yval'] = (powerlaw[:7] * 
                                        self.SixteenFocus_TwentyObj_Offaxis[1:8,2] 
                                        * RField[:7])
        self.NearFocusFarObject['xval'] = self.xval[1:8]

        self.UnderAccommFarObject['name'] = 'underaccomm, far object'
        self.UnderAccommFarObject['yval'] = (powerlaw[:7] * 
                                            self.Sixteen_UnderAccomm[1:8,2]
                                            * RField[:7])
        self.UnderAccommFarObject['xval'] = self.xval[1:8]
        
        self.TotalActivity()
        self.estimateInfo(Receptive_Field)

                  
    def TotalActivity(self):
        """Compute the estimated activity in a photoreceptor.
        
        :param self: no params yet.
        
        :returns: percent of integrated.
        
        .. note:: 
           Completely hard coded right now.  Would like to change that.
        
        """
        
        self.DiffractionLim['total'] = np.sum(self.DiffractionLim['yval'])
        self.Infinity['total'] = np.sum(self.Infinity['yval'])
        self.NearFocusNearObject['total'] = np.sum(self.NearFocusNearObject['yval'])
        self.NearFocusFarObject['total'] = np.sum(self.NearFocusFarObject['yval'])
        self.UnderAccommFarObject['total'] = np.sum(self.UnderAccommFarObject['yval'])


        self.Infinity['percent'] = (self.Infinity['total'] / 
                                        self.DiffractionLim['total'])
        self.NearFocusFarObject['percent'] = (self.NearFocusFarObject['total'] / 
                                                self.DiffractionLim['total'] )
        self.NearFocusNearObject['percent'] = (self.NearFocusNearObject['total'] /
                                                self.DiffractionLim['total'] )
        self.UnderAccommFarObject['percent'] = (self.UnderAccommFarObject['total'] /
                                                self.DiffractionLim['total'] )

        print ' '
        print 'Total activity (proportion of diffraction)'
        print '-------------------------------------------'
        print self.Infinity['name'], ': ', self.Infinity['percent']
        print self.NearFocusFarObject['name'], ': ', self.NearFocusFarObject['percent']
        print self.NearFocusNearObject['name'], ': ', self.NearFocusNearObject['percent']
        print self.UnderAccommFarObject['name'], ': ', self.UnderAccommFarObject['percent']

        
    def estimateInfo(self, Receptive_Field):
        """Estimate the information in a simple linear cone receptive field.
        
        :param Receptive_Field: type of receptive field to use (FFT, Jay)
        
        This function is under development. It will be called by 
        ComputeConeActivity() function.

        """ 
        if Receptive_Field == None:
            raise('Run {0} receptive field function'.format(Receptive_Field))
                    
        #imagexval = np.arange(1,self.amp_mean.shape[0] + 1) / 46.0
        
        #RField = interpolate.splev(imagexval, Receptive_Field, der = 0)
        RField = interpolate.splev(self.xval[1:], Receptive_Field, der = 0)
        RField = RField / sum(RField)

        cones = [1,2,3,4,5,6]
        total_images = len(self.ampSpecs)
        self.DiffractionLim['info'] = np.zeros(len(cones))
        for amp in self.ampSpecs:
            #Diffraction:
            fooInfo = info.SingleConeEntropyFunc((amp[:60] * 10e+5 * 
                                                self.DiffractionLim['yval'] * 
                                                RField), cones)
            self.DiffractionLim['info'] += fooInfo / total_images
            
            
            #Infinity
            fooInfo = info.SingleConeEntropyFunc((amp[:60] * 10e+5 * 
                                                    self.Infinity['yval'] * 
                                                    RField), cones)
            self.Infinity['info'] = fooInfo

            #Near focus, far object                            
            fooInfo = info.SingleConeEntropyFunc((amp[1:8] * 10e+5 * 
                                            self.NearFocusFarObject['yval'] * 
                                            RField[1:8]), cones)
            self.NearFocusFarObject['info'] = fooInfo 

            #Near focus, near object                            
            fooInfo = info.SingleConeEntropyFunc((amp[1:60] * 10e+5 * 
                                            self.NearFocusNearObject['yval'] * 
                                            RField[1:60]), cones)
            self.NearFocusNearObject['info'] = fooInfo 

            #Under accommodated                           
            fooInfo = info.SingleConeEntropyFunc((amp[1:8] * 10e+5 *
                                            self.UnderAccommFarObject['yval'] * 
                                            RField[1:8]), cones)
                                            
            self.UnderAccommFarObject['info'] = fooInfo 

        print ' '
        print 'Information'
        print '------------'
        print self.DiffractionLim['name'], ': ', self.DiffractionLim['info']
        print self.Infinity['name'], ': ', self.Infinity['info']
        print self.NearFocusFarObject['name'], ': ', self.NearFocusFarObject['info']
        print self.NearFocusNearObject['name'], ': ', self.NearFocusNearObject['info']
        print self.UnderAccommFarObject['name'], ': ', self.UnderAccommFarObject['info']        
 
        plt.figure(figsize=(8,6))
        plt.plot(self.DiffractionLim['info'], label=self.DiffractionLim['name'])
        plt.plot(self.Infinity['info'], label=self.Infinity['name'])
        plt.plot(self.NearFocusFarObject['info'], label=self.NearFocusFarObject['name'])
        plt.plot(self.NearFocusNearObject['info'], label=self.NearFocusNearObject['name'])
        plt.plot(self.UnderAccommFarObject['info'], label=self.UnderAccommFarObject['name'])
        plt.legend(loc='upper left')
        plt.tight_layout()
        plt.show()
        
        
    ### ANALYSIS and PLOTTING FUNCTIONS ###
    #######################################

    
    def DiffOfGaussian(self, excite_SD = 0.5, inhibit_SD = 5.0, plot_opt = True,
                       save_plots = False):
        """Create a difference of gaussians receptive field and plot it.

        :param excite_SD: standard deviation parameter of excitatory gaussian.
        :param inhibit_SD: standard deviation parameter of inhibitory \
        surround gaussian. Default = 5.0.   
        :param plot_opt: decide whether to output plots (True) or not (False)
        :type plot_opt: int
        :param save_plots: decide whether to save plots (True) or not (False)
        :returns: Plot1: difference of gaussians receptive field \n
                    Plot2: FFT spectrum of plot 1.
                
        .. note:: This functions needs to become more flexible. Should \
        eventually add more cone models
        
        
        Currently this function outputs the following:

            
        .. figure:: ../../Figures/ConeRF.png 
           :height: 300px
           :width: 400px
           :align: center   
           
           **Fig 1:** A simple Diff of Gaussian receptive field

           
        .. figure:: ../../Figures/ConeRF_FFT.png  
           :height: 300px
           :width: 400px
           :align: center
           
           **Fig 2:** And converted into a probability density and Fourier 
           transformed

        """
        
        N = 400
        Xvals = np.arange(-15, 15, 10./ (2.0 * N) )
        
        gauss1 = lambda x : 1.0*np.exp(-(x)**2 / (2 * excite_SD**2)) 
        self.y_excite = gauss1(Xvals)

        #self.y_excite = ( max(self.y_excite) / self.y_excite) * a
        
        gauss2 = lambda x : 1.0*np.exp(-(x)**2 / (2 * inhibit_SD**2))     
        y_inhibit = gauss2(Xvals) 

        normFact = sum(self.y_excite) / sum(y_inhibit)
        y_inhibit *= normFact
        

        self.RF_DOG = self.y_excite - y_inhibit
        self.RF_DOG = self.RF_DOG / max(self.RF_DOG)
        
        self.FFT_RF = (np.fft.fftshift(np.abs(np.fft.fft(self.RF_DOG))) 
                        / np.sqrt(2 * N))

        length = np.floor(self.FFT_RF.shape[0] / 2.) + 1       
        ## set up for interpolation
        self.RF_SPLINE = interpolate.splrep(Xvals[length:] * 60,
                                            self.FFT_RF[length:], s=0)
        
        RF = interpolate.splev(Xvals[length:] * 60, self.RF_SPLINE, der = 0)
        RF = RF / np.sum(RF)
            
        if plot_opt:
            fig = plt.figure(figsize=(8,6))
            ax = fig.add_subplot(111)

            pf.AxisFormat()
            pf.TufteAxis(ax, ['left', 'bottom'])
            
            ax.plot(Xvals, self.RF_DOG, 'k', linewidth=2.5)
            #ax.plot(Xvals, self.FFT_RF, linewidth=2.5)
            
    
            ax.get_xaxis().tick_bottom()
            ax.get_yaxis().tick_left()
    
            plt.xlabel('distance (arcmin)')
            plt.ylabel('amplitude')
            
            plt.tight_layout()
            
            if save_plots:
                fig.show()
                fig.savefig('../../Figures/ConeRF.png')
                plt.close()
            else:
                plt.show()
            
            fig2 = plt.figure(figsize=(8,6))
            ax = fig2.add_subplot(111)
            pf.AxisFormat()
            pf.TufteAxis(ax, ['left', 'bottom'])
            
            ax.loglog(Xvals[length:] * 60, RF, 'k',
                          linewidth=2.5)
            
    
            ax.get_xaxis().tick_bottom()
            ax.get_yaxis().tick_left()
    
            plt.ylim([10**-4, 10**-1])
            plt.xlim([self.xval[1], 100])
            plt.xlabel('spatial frequency (cycles / deg)')
            plt.ylabel('density')
            
            plt.tight_layout()
            
            if save_plots:
                fig2.show()
                fig2.savefig('../../Figures/ConeRF_FFT.png')
                plt.close()
            else:
                plt.show()
        

    def DeconstructedFFT(self, plot_option = True, save_plots = False):
        """This one is for Jay to demonstrate the logic behind fft.
        
        :param plot_option: decide whether to plot output (True) or not \
        (False).
        :type plot_option: bool
        :param save_plots: decide whether to save the plot (True) or not\
        (False)
        :type save_option: bool
        
        :returns: two plots.
        
        .. figure:: ../../Figures/RecField_JAY.png
           :height: 300px
           :width: 400px
           :align: center
           
           **Fig 1:** Receptive field (black) with three spatial frequencies \
           (blue = 1 cpd, green = 5 cpd and red = 10 cpd)


        .. figure:: ../../Figures/FFT_JAY.png
           :height: 300px
           :width: 400px
           :align: center
           
           **Fig 2:** Response of receptive field to sine waves. 
           
        """
        
        if self.RF_DOG == None:
            self.DiffOfGaussian(plot_opt = 0)
            
        N = 400
        Xvals = np.arange(-15, 15, 10./ (2.0 * N) )
        length = np.floor(self.FFT_RF.shape[0] / 2.) + 1
        
        spatial_frequencies = np.arange(1, length - 1)

        # allocate memory        
        sine_wave = np.zeros((Xvals.shape[0], spatial_frequencies.shape[0]))
        cone_response = np.zeros(spatial_frequencies.shape[0])        
        for i, thisFrequency in enumerate(spatial_frequencies):
            
            # convert from cpd     (radians    / arcmin)
            Converted_freq = thisFrequency * (2 * np.pi) / 60.0 
            
            sine_wave[:,i] = ( 1.0 + np.sin(Xvals * Converted_freq 
                                        + (np.pi / 2.0) ))/ 2.0
            cone_response[i] = np.sum(sine_wave[:,i] * self.RF_DOG)

        self.Jay_RF = interpolate.splrep(Xvals[length:] * 60,
                                            cone_response, s=0)
        RF = interpolate.splev(Xvals[length:], self.Jay_RF, der = 0)
        RF = RF / np.sum(RF)
        
        
        if plot_option:
            
            fig = plt.figure(figsize=(8,6))
            ax = fig.add_subplot(111)
            pf.AxisFormat()
            pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[7,5])
            
            ax.plot(Xvals, sine_wave[:,1], 'b-', linewidth =2, label='1 cpd')
            ax.plot(Xvals, sine_wave[:,5], 'g-', linewidth = 2, label='5 cpd')
            ax.plot(Xvals, sine_wave[:,10], 'r-', linewidth = 2, 
                    label='10 cpd')
            ax.plot(Xvals, self.RF_DOG, 'k-', linewidth = 3)
            
            plt.xlim([min(Xvals), max(Xvals)])
            plt.ylim([-0.15, 1.05])
    
            plt.xlabel('distance (arcmin)')
            plt.ylabel('amplitude')      
            
            plt.tight_layout()
    
            if save_plots:
                fig.show()
                fig.savefig('../../Figures/RecField_JAY.png')
                plt.close()
            else:
                plt.show()
    
            
            if save_plots:
                fig2 = plt.figure(figsize=(8,6))
                ax = fig2.add_subplot(111)
            else:
                fig = plt.figure(figsize=(8,6))
                ax = fig.add_subplot(111)
                
            pf.AxisFormat()
            pf.TufteAxis(ax, ['left', 'bottom'])
            
            ax.loglog(Xvals[length:] * 60, RF, 'k-', linewidth = 2.5)
            
            plt.ylim([10**-4, 10**-1])
            plt.xlim([self.xval[1], 100])
            
            plt.xlabel('spatial frequency (cycles / deg)')
            plt.ylabel('density')      
            
            plt.tight_layout()
            
    
            if save_plots:
                fig2.show()
                fig2.savefig('../../Figures/FFT_JAY.png')
                plt.close()
            else:
                plt.show()
            
            
    def PlotPowerSpec(self, save_plots = False):
        """Plot the amplitude spec of image.
    
        :param save_plots: decide whether to save plots (True) or not (False)
        :type save_plots: bool        
        
        **This function produces:**
        
        .. figure:: ../../Figures/ampSpec.png
           :height: 300px
           :width: 400px
           :align: center  

           **Fig 1:** An amplitude spectrum of a series of images.
           
        """
        
        imagexval = np.arange(1,self.amp_mean.shape[0] + 1) / 46.0
        if not self.powerlaw:
            self.PowerLaw(imagexval, self.amp_mean)
            
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'])

        ax.loglog(imagexval , self.amp_mean,
                    'r.', markersize = 10, markeredgewidth = 0)
     
        ax.loglog(imagexval, self.powerlaw(imagexval), 'k', linewidth = 2.5)
        
        ax.text(1, 10**-2.4, r'$\frac{1}{\mathit{f}}$', size = 35)

        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('density')
        
        plt.tight_layout()
        
        if save_plots:
            fig.show()
            fig.savefig('../../Figures/ampSpec.png')
            plt.close()
        else:
            plt.show()

    
    def MTFfamilyPlots(self, save_plots = False, plot_option = 1, 
                       legend = False):    
        """Plot a family of MTF curves derived from a schematic eye.
        
        :param save_plots: decide whether to save plots (True) or not (False).
        :type save_plots: bool
        :param plot_option: decide what to plot,
                        1 = on axis
                        2 = off axis
        :param legend: turn legend on (True) or off (False). Default = True. \n
        :type legend: bool
        
        :returns: Modulation transfer function plots for a family of curves.
        
        .. note:: 
           * Need to add options for what to plot. 
           * Legend option is not currently working properly.
        
        **This produces:**
        
        .. figure:: ../../Figures/MTFfamily.png
           :height: 300px
           :width: 400px
           :align: center        
           
           **Fig 1:** A family of MTF curves from OSLO data.
           
        """
        
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], [5,5])

        ''' on axis plots '''
        ax.plot(self.xval, self.INF[:,4], 'k', linewidth = 2.5, 
                label='diffraction ')
        ax.plot(self.xval, self.INF[:,2], 'r', linewidth=2.5, 
                label='infinity') 
        
        if plot_option == 1:
            ax.plot(self.xval, self.TwentyFt[:,2], 'g', linewidth=2.5, 
                    label='20ft')
            ax.plot(self.xval, self.Onemeter[:,2], 'b', linewidth=2.5, 
                    label='1m')
            ax.plot(self.xval[:20], self.SixteenIn[:20,2], 'm', 
                    linewidth=2.5, label='16in')
        
        
        ''' plot2 : off axis plots '''
        
        
        if plot_option ==3:
            ax.plot(self.xval, self.INF_offaxis[:,4], 'k--', linewidth=2.5)
            
            ax.plot(self.xval, self.TwentyFt_offaxis[:,2], 'g--', 
                    linewidth=2.5)
            ax.plot(self.xval, self.Onemeter_offaxis[:,2], 'b--', 
                    linewidth=2.5)
            ax.plot(self.xval[:20], self.SixteenIn_offaxis[:20,2], 'm--', 
                    linewidth=2.5)
        
        
        if plot_option == 2:
            ax.plot(self.xval, self.INF_offaxis[:,2], 'r--', linewidth=2.5)
            
            ax.plot(self.xval[:60], 
                    self.SixteenFocus_SixteenObj_Offaxis[:60,2],
                    'g--', linewidth=2.5, label = 'near focus, near obj')
            ax.plot(self.xval[:20], self.Sixteen_UnderAccomm[:20,2], 
                    'c--', linewidth=2.5, label = 'underacc, far obj')
            ax.plot(self.xval[:20], self.SixteenFocus_TwentyObj_Offaxis[:20,2], 
                    'b--',linewidth=2.5, label = 'near focus, far obj')
        
        if legend: 
            ax.legend(loc='upper right')#title='object, retina')
        
        
        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('modulation')
        
        plt.tight_layout()
        
        if save_plots:
            fig.show()
            if plot_option == 1:
                fig.savefig('../../Figures/MTFfamilyOnAxis.png')
            if plot_option == 2:
                fig.savefig('../../Figures/MTFfamily.png')
            plt.close()
        else:
            plt.show()
        
        """ Legend figure """
        
        if save_plots:
            fig2 = plt.figure(figsize=(8,6))
            ax = fig2.add_subplot(111)
        else:
            fig = plt.figure(figsize=(8,6))
            ax = fig.add_subplot(111)
            
        font = {'weight' : 'norm', 'size'  : 32}
        legend = {'frameon' : False}
        plt.rc('font', **font)
        
        
        ax.plot(0,0, 'k', linewidth=2.5, label='diffraction')
        ax.plot(0,0, 'r', linewidth=2.5, label='infinity') 
        ax.plot(0,0, 'g--', linewidth=2.5, label = 'near focus, near object')
        ax.plot(0,0, 'c--', linewidth=2.5, label = 'underacc, far object')
        ax.plot(0,0, 'b--', linewidth=2.5, label = 'near focus, far object')
        ax.legend(loc = 'center')
        plt.axis('off')
        
        if save_plots:
            fig2.show()
            fig2.savefig('../../Figures/Legend.png')
            plt.close()
        else:
            plt.show()




    def PeripheralPlot(self, save_plots=False, legend=False):
        """
        Plot peripheral MTF with a comparison to Navarro et al 1993 or 
        Williams et al. 1996
        
        :param save_plots: decide whether to save plots (True) or not (False).
        :type save_plots: bool
        :param legend: turn legend on (True) or off (False). Default = True. \n
        :type legend: bool
        
        Currently supports 0, 10, 20, 40 degrees eccentricity.
        
        **This produces:**
        
        .. figure:: ../../Figures/MTFperiphery.png
           :height: 300px
           :width: 400px
           :align: center        
           
           **Fig 1:** A family of MTF curves from experimental data (dotted) 
           and schematic eye.        
        """
        from eye.Optics import MTF
        
        Fovea = MTF(self.xval[:], 0)
        TenDeg = MTF(self.xval[:], 10)
        TwentyDeg = MTF(self.xval[:],20)
        FourtyDeg = MTF(self.xval[:],40)
        
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)

        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], [5,5])
        
        #Navarro et al 1993 analytical func:
        ax.plot(self.xval[:], Fovea, 'm--', label='fovea', linewidth=2.5)
        ax.plot(self.xval[:], TenDeg, 'r--', label='10 deg', linewidth=2.5)
        ax.plot(self.xval[:], TwentyDeg, 'g--', label='20 deg', linewidth=2.5)
        ax.plot(self.xval[:], FourtyDeg, 'b--',label='40 deg', linewidth=2.5)
        
        #OSLO ray trace data:
        ax.plot(self.xval[:], self.INF[:,2], 'm-', label='fovea', 
                linewidth=2.5)
        ax.plot(self.xval[:], self.INF_offaxis[:,2], 'r-', label='10 deg',
                linewidth=2.5)        
        ax.plot(self.xval[:], self.TwentyDegOffAxis_InfFoc[:,2], 'g-', 
                label='20 deg', linewidth=2.5)
        ax.plot(self.xval[:30], self.FortyDegOffAxis_InfFoc[:30,2], 'b-',
                label='40 deg', linewidth=2.5)
                
        if legend: 
            ax.legend(loc='lower left')#,title='object dist, retinal location')
        
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        
        mi, ma = plt.ylim()
        
        #plt.ylim([10**-2.5, 10**0])
        #plt.xlim([self.xval[1], 100])
        
        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('modulation')
        
        plt.tight_layout()
        
        if save_plots:
            fig.show()
            fig.savefig('./Figures/MTFperiphery.png')
            plt.close()
        else:
            plt.show()
                
                    
    def AmpModPlot(self, Receptive_Field = 'FFT', save_plots = False, 
                   legend = False):
        """Plot powerlaw amplitude spectrum after accounting for MTF and MTF \
        + receptive field
        
        :param Receptive_Field: Indicate which receptive field to use. \
                                Option are 'FFT' or 'Jay'
        :type Receptive_Field: str        
        :param legend: turn legend on (True) or off (False). Default = True.
        :type legend: bool
        :param save_plots: decide whether to save a plot (True) or not (False)
        :type save_plots: bool
        :returns: Plot 1: power law times Modulation Transfer Function for 
                  three conditions \n
                  Plot 2: same plot with difference of gaussians receptive 
                  field included as well.  This is the final plot in this \
                  series.
        :rtype: plt.plot

        .. note:: ../../Figures/
           Eventually would like to introduce a heuristic to determine the 
           location of the power law text that is plotted. Currently hard \
           coded.
           
        
        
        **This function plots estimates of the photoreceptor activity**
        
        .. figure:: ../../Figures/MTFfamilyMod.png
           :height: 300px
           :width: 400px
           :align: center    
           
           **Fig 1:** Amplitude * MTF
        
        .. figure:: ../../Figures/MTFfamilyModAmp.png
           :height: 300px
           :width: 400px
           :align: center  
           
           **Fig2:** Amplitude * MTF * Receptive field
           
        """

        if self.RF_SPLINE == None:
            self.DiffOfGaussian()

        if Receptive_Field.lower() == 'jay':
            Rec_Field = self.Jay_RF
            print 'Used Jay receptive field'
            
        elif Receptive_Field.lower() == 'fft':
            Rec_Field = self.RF_SPLINE
            print 'Used FFT receptive field'
            
        else :
            Rec_Field = self.RF_SPLINE
            print 'receptive field not understood, see options. Using FFT.'
            
        self.ComputeConeActivity(Rec_Field)   
        
            
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)

        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'])
        
        powerlaw = self.powerlaw(self.xval[1:])
        
        ax.loglog(self.xval[1:], self.powerlaw(self.xval[1:]), 'k',
                  linewidth = 2.5)
        
        ax.loglog(self.xval[1:], self.powerlaw(self.xval[1:]) * self.INF[1:,2], 
                  'r', linewidth=2.5, label='inf focus')
                  
        ax.loglog(self.xval[1:60], 
                  powerlaw[:59] * self.SixteenFocus_SixteenObj_Offaxis[1:60,2],
                  'g--', linewidth=2.5, label='near focus, near object')          
        ax.loglog(self.xval[1:8], 
                  powerlaw[:7] * self.SixteenFocus_TwentyObj_Offaxis[1:8,2],
                  'b--', linewidth=2.5, label='near focus, far object')
        ax.loglog(self.xval[1:8], 
                  powerlaw[:7] * self.Sixteen_UnderAccomm[1:8,2],
                  'c--', linewidth=2.5, label='underaccom, far object')                  
        
        if legend: 
            ax.legend(loc='lower left')#,title='object dist, retinal location')
        
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        
        mi, ma = plt.ylim()
        
        ax.text(15, 10**-4.0, r'$\frac{1}{\mathit{f}}$', size = 35)
        
        plt.ylim([10**-8, 10**-3])
        plt.xlim([self.xval[1], 100])
        
        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('density')
        
        plt.tight_layout()
        
        if save_plots:
            fig.show()
            fig.savefig('./Figures/MTFfamilyMod.png')
            plt.close()
        else:
            plt.show()
        
        """ 
        plot 2 
        """
       
        if save_plots:
            fig2 = plt.figure(figsize=(8,6))
            ax = fig2.add_subplot(111)
        else:
            fig = plt.figure(figsize=(8,6))
            ax = fig.add_subplot(111)
            
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'])
     

        ax.loglog(self.DiffractionLim['xval'], 
                  self.DiffractionLim['yval'],
                  'k', linewidth = 2.5)

        ax.loglog(self.Infinity['xval'], 
                  self.Infinity['yval'] ,
                  'r', linewidth=2.5, 
                  label= self.Infinity['name'])    
        
        ax.loglog(self.NearFocusNearObject['xval'], 
                  self.NearFocusNearObject['yval'],
                  'g--', linewidth=2.5, 
                  label = self.NearFocusNearObject['name'])   

        ax.loglog(self.NearFocusFarObject['xval'], 
                  self.NearFocusFarObject['yval'],
                  'b--', linewidth=2.5, 
                  label = self.NearFocusFarObject['name'])
                  
        ax.loglog(self.UnderAccommFarObject['xval'], 
                  self.UnderAccommFarObject['yval'] ,
                  'c--', linewidth=2.5, 
                  label =self.UnderAccommFarObject['name'])                  
        
        if legend: 
            ax.legend(loc='upper right')
        

        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        
        mi, ma = plt.ylim()

        ax.text(12, 10**-5.0, r'$\frac{1}{\mathit{f}}$', size = 35)
        
        plt.ylim([10**-8, 10**-3.75])
        plt.xlim([self.xval[1], 100])
        
        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('density')
        
        plt.tight_layout()
       
        if save_plots:
            if Receptive_Field.lower() == 'jay':
                save_name = '../../Figures/MTFfamilyModAmpJAY.png'
            else:
                save_name = '../../Figures/MTFfamilyModAmp.png'
            
            fig2.show()
            fig2.savefig(save_name)     
            plt.close()
            
        else:
            plt.show()