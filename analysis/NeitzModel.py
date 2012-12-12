import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate


from eye.eyeModel import SchematicEye
from scene.Images import Images
from analysis import PlottingFun as pf
from analysis import Information as info

class SchematicEyeAnalysis(Images,SchematicEye):
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
    """
    
    def __init__(self):
        """
        
        """
        
        SchematicEye.__init__(self)
        Images.__init__(self)
        
        self.powerlaw = None
        self.RF_DOG, self.RF_SPLINE, self.FFT_RF = DiffOfGaussian(
                                        excite_SD = 0.5, inhibit_SD = 5.0)
        self.Jay_RF = None
        
        self.Analysis = {}
        self.Analysis['diffract periph'] = {}
        self.Analysis['diffract fovea'] = {}
        self.Analysis['inf perip'] = {}
        self.Analysis['inf fovea'] = {}
        self.Analysis['near focus, far object'] = {}
        self.Analysis['near focus, near object'] = {}
        self.Analysis['underaccomm, far object'] = {}


        

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
            
        powerlaw = self.powerlaw(self.freqs[1:])

        #peripheral RF
        RField_P = interpolate.splev(self.freqs[1:], Receptive_Field, der = 0)
        RField_P = RField_P / sum(RField_P)
        
        #foveal RF
        a, FOVEA_FFT_RF, foo = DiffOfGaussian(excite_SD=0.1, inhibit_SD=0.5)
        RField_F = interpolate.splev(self.freqs[1:], FOVEA_FFT_RF, der = 0)
        RField_F = RField_F / sum(RField_F)
        
        for key in self.Analysis:
            if key == 'diffract periph':
                self.Analysis[key]['retina'] = powerlaw * RField_P
                self.Analysis[key]['freqs'] = self.freqs[1:]

            if key == 'diffract fovea':
                self.Analysis[key]['retina'] = powerlaw * RField_F
                self.Analysis[key]['freqs'] = self.freqs[1:]
                
            if key == 'inf perip':
                self.Analysis[key]['retina'] = (powerlaw * self.INF[1:, 2] * 
                                                RField_P)
                self.Analysis[key]['freqs'] = self.freqs[1:]        

            if key == 'inf fovea':
                self.Analysis[key]['retina'] = (powerlaw * self.INF[1:, 2] * 
                                                RField_F)
                self.Analysis[key]['freqs'] = self.freqs[1:]  
                 
            if key == 'near focus, near object':
                self.Analysis[key]['retina'] = (powerlaw[:59] * 
                                 self.SixteenFocus_SixteenObj_Offaxis[1:60,2] * 
                                                RField_P[:59])
                self.Analysis[key]['freqs'] = self.freqs[1:60]
            
            if key == 'near focus, far object':
                self.Analysis[key]['retina'] = (powerlaw[:7] * 
                                 self.SixteenFocus_TwentyObj_Offaxis[1:8,2] * 
                                                RField_P[:7])
                self.Analysis[key]['freqs'] = self.freqs[1:8]        
            
            if key == 'underaccomm, far object':
                self.Analysis[key]['retina'] = (powerlaw[:7] * 
                                             self.Sixteen_UnderAccomm[1:8,2] * 
                                                RField_P[:7])
                self.Analysis[key]['freqs'] = self.freqs[1:8]        

        
        self.TotalActivity()
        self.estimateInfo(RField_P, RField_F)

                  
    def TotalActivity(self):
        """Compute the estimated activity in a photoreceptor.
        
        :param self: no params yet.
        
        :returns: percent of integrated.
        
        .. note:: 
           Completely hard coded right now.  Would like to change that.
        
        """
        
        for key in self.Analysis:
            self.Analysis[key]['total'] = np.sum(self.Analysis[key]['retina'])
            self.Analysis[key]['percent'] = (self.Analysis[key]['total'] /
                                    self.Analysis['diffract periph']['total'])


        print ' '
        print 'Total activity (proportion of diffraction)'
        print '-------------------------------------------'
        for key in self.Analysis:
            print key, ': ', self.Analysis[key]['percent']


        
    def estimateInfo(self, RField_P, RField_F, print_option=False, 
                     plot_opt=True, save_plots=False, legend=False):
        """Estimate the information in a simple linear cone receptive field.
        
        Information is estimated with Garrigan et al.'s Gaussian approximation
        method.
        
        .. math::
           I_1(S,E) = \\frac{1}{2}\log_2(1+SNR) , 
        
        where :math:`SNR` is the ratio of signal power to noise power.
        
        Information in an array of cones is then represented:
        
        .. math::
           I_N(S,E) = I_1(S,E)N^\delta ,
        
        with :math:`\\delta` representing a derived scale factor.
        
        :param Receptive_Field: type of receptive field to use (FFT, Jay)
        :param print_opt: decide whether to print results (True) or not (False)
        :type print_opt: bool
        :param plot_opt: decide whether to output plots (True) or not (False)
        :type plot_opt: bool
        :param save_plots: decide whether to save plots (True) or not (False)
        :type save_plots: bool
        :param legend: decide whether to add a legend (True) or not (False)
        :type legend: bool
        
        :returns: plot of information estimate
        
        .. warning::
           This function is under development. 
        
        :usage: estimateInfo is called by ComputeConeActivity() function.

        **This function produces:**
        
        .. figure:: ../../Figures/InfoPlot.png 
           :height: 300px
           :width: 400px
           :align: center   
           
        """ 


        
        cones = [1,2,3,4,5,6]
        total_images = len(self.ampSpecs)
        
        for key in self.Analysis:
            self.Analysis[key]['info'] = np.zeros(len(cones))

        for amp in self.rawAmp:
            for key in self.Analysis:
                
                #Diffraction:
                if key == 'diffract periph':
                    fooInfo = info.SingleConeEntropyFunc((amp[:60]**2 *
                                                self.Analysis[key]['retina'] * 
                                                RField_P), cones)
                    self.Analysis[key]['info'] += fooInfo / total_images

                if key == 'diffract fovea':
                    fooInfo = info.SingleConeEntropyFunc((amp[:60]**2 *
                                                self.Analysis[key]['retina'] * 
                                                RField_F), cones)
                    self.Analysis[key]['info'] += fooInfo / total_images
           
                #Infinity
                if key == 'inf perip':
                    fooInfo = info.SingleConeEntropyFunc((amp[:60]**2 * 
                                                    self.Analysis[key]['retina'] * 
                                                    RField_P), cones)
                    self.Analysis[key]['info'] += fooInfo / total_images

                if key == 'inf fovea':
                    fooInfo = info.SingleConeEntropyFunc((amp[:60]**2 * 
                                                    self.Analysis[key]['retina'] * 
                                                    RField_F), cones)
                    self.Analysis[key]['info'] += fooInfo / total_images
    
                if key == 'near focus, far object':                           
                    fooInfo = info.SingleConeEntropyFunc((amp[1:8]**2 *
                                            self.Analysis[key]['retina'] * 
                                            RField_P[1:8]), cones)
                    self.Analysis[key]['info'] += fooInfo / total_images

                if key == 'near focus, near object':            
                    fooInfo = info.SingleConeEntropyFunc((amp[1:60]**2 * 
                                            self.Analysis[key]['retina'] * 
                                            RField_P[1:60]), cones)
                    self.Analysis[key]['info'] += fooInfo / total_images

                if key == 'underaccomm, far object':                          
                    fooInfo = info.SingleConeEntropyFunc((amp[1:8]**2 *
                                            self.Analysis[key]['retina'] * 
                                            RField_P[1:8]), cones)
                                            
                    self.Analysis[key]['info'] += fooInfo / total_images 

        if print_option == True:
            print ' '
            print 'Information'
            print '------------'
            for key in self.Analysis:
                print key, ': ', self.Analysis[key]['info']        
 
        if plot_opt == True:
            
            fig = plt.figure(figsize=(8,6))
            ax = fig.add_subplot(111)

            pf.AxisFormat()
            pf.TufteAxis(ax, ['left', 'bottom'], [5,5])
                    
            plt.plot(cones, self.Analysis['diffract periph']['info'], 
                     'ko-',
                     label='diffract periph',
                     linewidth=2.5, markersize=10)

            plt.plot(cones, self.Analysis['diffract periph']['info'], 
                     'ko-.',
                     label='diffract fova',
                     linewidth=2.5, markersize=10)
                     
            plt.plot(cones, 
                     self.Analysis['inf perip']['info'], 
                     'ro-',
                     label='inf perip',
                     linewidth=2.5, markersize=10)

            plt.plot(cones, 
                     self.Analysis['inf fovea']['info'], 
                     'ro-.',
                     label='inf perip',
                     linewidth=2.5, markersize=10)
                     
                     
            plt.plot(cones, 
                     self.Analysis['near focus, far object']['info'], 
                     'bo--',
                     label='near focus, far object',
                     linewidth=2.5, markersize=10)
                     
            plt.plot(cones, 
                     self.Analysis['near focus, near object']['info'], 
                     'go--',
                     label='near focus, near object',
                     linewidth=2.5, markersize=10)
                     
            plt.plot(cones, 
                     self.Analysis['underaccomm, far object']['info'],
                     'co--',
                     label='underaccomm, far object',
                     linewidth=2.5, markersize=10)

            plt.xlim([min(cones)-0.1, max(cones)+0.1])
            if legend:
                plt.legend(loc='upper left')
            plt.xlabel('cones in mosaic')
            plt.ylabel('entropy (bits)')
            
            plt.tight_layout()
            
            if save_plots:
                fig.show()
                fig.savefig('../../Figures/InfoPlot.png')
                plt.close()
            else:
                plt.show()
            
        
    ### ANALYSIS and PLOTTING FUNCTIONS ###
    #######################################

    

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
            plt.xlim([self.freqs[1], 100])
            
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
        ax.plot(self.freqs, self.INF[:,4], 'k', linewidth = 2.5, 
                label='diffraction ')
        ax.plot(self.freqs, self.INF[:,2], 'r', linewidth=2.5, 
                label='infinity') 
        
        if plot_option == 1:
            ax.plot(self.freqs, self.TwentyFt[:,2], 'g', linewidth=2.5, 
                    label='20ft')
            ax.plot(self.freqs, self.Onemeter[:,2], 'b', linewidth=2.5, 
                    label='1m')
            ax.plot(self.freqs[:20], self.SixteenIn[:20,2], 'm', 
                    linewidth=2.5, label='16in')
        
        
        ''' plot2 : off axis plots '''
        
        
        if plot_option ==3:
            ax.plot(self.freqs, self.INF_offaxis[:,4], 'k--', linewidth=2.5)
            
            ax.plot(self.freqs, self.TwentyFt_offaxis[:,2], 'g--', 
                    linewidth=2.5)
            ax.plot(self.freqs, self.Onemeter_offaxis[:,2], 'b--', 
                    linewidth=2.5)
            ax.plot(self.freqs[:20], self.SixteenIn_offaxis[:20,2], 'm--', 
                    linewidth=2.5)
        
        
        if plot_option == 2:
            ax.plot(self.freqs, self.INF_offaxis[:,2], 'r--', linewidth=2.5)
            
            ax.plot(self.freqs[:60], 
                    self.SixteenFocus_SixteenObj_Offaxis[:60,2],
                    'g--', linewidth=2.5, label = 'near focus, near obj')
            ax.plot(self.freqs[:20], self.Sixteen_UnderAccomm[:20,2], 
                    'c--', linewidth=2.5, label = 'underacc, far obj')
            ax.plot(self.freqs[:20], self.SixteenFocus_TwentyObj_Offaxis[:20,2], 
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
        
        Fovea = MTF(self.freqs[:], 0)
        TenDeg = MTF(self.freqs[:], 10)
        TwentyDeg = MTF(self.freqs[:],20)
        FourtyDeg = MTF(self.freqs[:],40)
        
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)

        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], [5,5])
        
        #Navarro et al 1993 analytical func:
        ax.plot(self.freqs[:], Fovea, 'm--', label='fovea', linewidth=2.5)
        ax.plot(self.freqs[:], TenDeg, 'r--', label='10 deg', linewidth=2.5)
        ax.plot(self.freqs[:], TwentyDeg, 'g--', label='20 deg', linewidth=2.5)
        ax.plot(self.freqs[:], FourtyDeg, 'b--',label='40 deg', linewidth=2.5)
        
        #OSLO ray trace data:
        ax.plot(self.freqs[:], self.INF[:,2], 'm-', label='fovea', 
                linewidth=2.5)
        ax.plot(self.freqs[:], self.INF_offaxis[:,2], 'r-', label='10 deg',
                linewidth=2.5)        
        ax.plot(self.freqs[:], self.TwentyDegOffAxis_InfFoc[:,2], 'g-', 
                label='20 deg', linewidth=2.5)
        ax.plot(self.freqs[:30], self.FortyDegOffAxis_InfFoc[:30,2], 'b-',
                label='40 deg', linewidth=2.5)
                
        if legend: 
            ax.legend(loc='lower left')#,title='object dist, retinal location')
        
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        
        mi, ma = plt.ylim()
        
        #plt.ylim([10**-2.5, 10**0])
        #plt.xlim([self.freqs[1], 100])
        
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

        .. note::
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
        
        powerlaw = self.powerlaw(self.freqs[1:])
        
        ax.loglog(self.freqs[1:], self.powerlaw(self.freqs[1:]), 'k',
                  linewidth = 2.5)
        
        ax.loglog(self.freqs[1:], self.powerlaw(self.freqs[1:]) * self.INF[1:,2], 
                  'r', linewidth=2.5, label='inf focus')
                  
        ax.loglog(self.freqs[1:60], 
                  powerlaw[:59] * self.SixteenFocus_SixteenObj_Offaxis[1:60,2],
                  'g--', linewidth=2.5, label='near focus, near object')          
        ax.loglog(self.freqs[1:8], 
                  powerlaw[:7] * self.SixteenFocus_TwentyObj_Offaxis[1:8,2],
                  'b--', linewidth=2.5, label='near focus, far object')
        ax.loglog(self.freqs[1:8], 
                  powerlaw[:7] * self.Sixteen_UnderAccomm[1:8,2],
                  'c--', linewidth=2.5, label='underaccom, far object')                  
        
        if legend: 
            ax.legend(loc='lower left')#,title='object dist, retinal location')
        
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        
        mi, ma = plt.ylim()
        
        ax.text(15, 10**-4.0, r'$\frac{1}{\mathit{f}}$', size = 35)
        
        plt.ylim([10**-8, 10**-3])
        plt.xlim([self.freqs[1], 100])
        
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

        ax.loglog(self.Analysis['diffract periph']['freqs'], 
                  self.Analysis['diffract periph']['retina'],
                  'k', linewidth = 2.5)

        ax.loglog(self.Analysis['inf perip']['freqs'], 
                  self.Analysis['inf perip']['retina'] ,
                  'r', linewidth=2.5, 
                  label= 'inf perip')    
        
        ax.loglog(self.Analysis['near focus, far object']['freqs'], 
                  self.Analysis['near focus, far object']['retina'],
                  'g--', linewidth=2.5, 
                  label = self.NearFocusNearObject['name'])   

        ax.loglog(self.Analysis['near focus, near object']['freqs'], 
                  self.Analysis['near focus, near object']['retina'],
                  'b--', linewidth=2.5, 
                  label = self.NearFocusFarObject['name'])
                  
        ax.loglog(self.Analysis['underaccomm, far object']['freqs'], 
                  self.Analysis['underaccomm, far object']['retina'] ,
                  'c--', linewidth=2.5, 
                  label =self.UnderAccommFarObject['name'])                  
        
        if legend: 
            ax.legend(loc='upper right')
        

        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        
        mi, ma = plt.ylim()

        ax.text(12, 10**-5.0, r'$\frac{1}{\mathit{f}}$', size = 35)
        
        plt.ylim([10**-8, 10**-3.75])
        plt.xlim([self.freqs[1], 100])
        
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
            
def DiffOfGaussian(excite_SD = 0.5, inhibit_SD = 5.0, plot_opt = False,
                   save_plots = False):
    """Create a difference of gaussians receptive field and plot it.

    :param excite_SD: standard deviation parameter of excitatory gaussian.
    :param inhibit_SD: standard deviation parameter of inhibitory \
    surround gaussian. Default = 5.0.   
    :param plot_opt: decide whether to output plots (True) or not (False)
    :type plot_opt: bool
    :param save_plots: decide whether to save plots (True) or not (False)
    :type save_plots: bool
    
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
    y_excite = gauss1(Xvals)

    
    gauss2 = lambda x : 1.0*np.exp(-(x)**2 / (2 * inhibit_SD**2))     
    y_inhibit = gauss2(Xvals) 

    normFact = sum(y_excite) / sum(y_inhibit)
    y_inhibit *= normFact
    

    RF_DOG = y_excite - y_inhibit
    RF_DOG = RF_DOG / max(RF_DOG)
    
    FFT_RF = (np.fft.fftshift(np.abs(np.fft.fft(RF_DOG))) 
                    / np.sqrt(2 * N))

    length = np.floor(FFT_RF.shape[0] / 2.) + 1       
    ## set up for interpolation
    RF_SPLINE = interpolate.splrep(Xvals[length:] * 60,
                                        FFT_RF[length:], s=0)
    
    RF = interpolate.splev(Xvals[length:] * 60, RF_SPLINE, der = 0)
    RF = RF / np.sum(RF)
        
    if plot_opt:
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)

        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'])
        
        ax.plot(Xvals, RF_DOG, 'k', linewidth=2.5)
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
        plt.xlim([1, 100])
        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('density')
        
        plt.tight_layout()
        
        if save_plots:
            fig2.show()
            fig2.savefig('../../Figures/ConeRF_FFT.png')
            plt.close()
        else:
            plt.show()

    return RF_DOG, RF_SPLINE, FFT_RF


if __name__ == "__main__":
    print 'nothing doing'