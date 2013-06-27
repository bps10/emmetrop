from __future__ import division
import matplotlib.pylab as plt
import numpy as np

from emmetrop.renderer import PlottingFun as pf
from emmetrop.scene import SignalProcessing as sig

class Plotter(object):
    """A plotting repo

    .. todo::
       Add self.ext to allow user to save as png or svg.
    """
    def __init__(self, analysis_args, diff, Analysis, recepitive_field, 
                imageData, plots, save_plots, legend):

        # import our dictionaries of data:
        self.Analysis = Analysis
        self.rec_field = recepitive_field
        self.imageData = imageData
        self.diffract = diff
        
        self.freqs = diff['freqs']
        foo = '../../../bps10.github.com/presentations/static/'
        self.figPath = foo + 'figures/myopiaModel/'
        #options:      
        self.location = ['periph']
        if 'diffract fovea' in self.Analysis:
            self.location.append('fovea')
        
        self.min_dB = -20
        
        # plot the appropriate plots, with options:
        if 'amp' in plots:
            self.plotAmpSpec(_brownian=True, save_plots=save_plots)
        if 'comp' in plots:
            self.plotComparisonMTF(save_plots)
        if 'mtf' in plots:
            self.plotMTFfamily(save_plots)
        if 'plotDoG' in plots:
            self.plotDoG(save_plots)       
        if 'activity' in plots:
            self.plotActivity(save_plots)
        if 'info' in plots:
            self.plotInformation(save_plots, legend)
        if 'seriesPlot' in plots:
            self.plotSeries(analysis_args, save_plots)

    def plotSeries(self, analysis_args, save_plots=False):
        '''
        '''
        from mpl_toolkits.mplot3d import Axes3D

        #for analysis in analysis_args:

        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111, projection='3d')     
        pf.AxisFormat()

        keys = 0
        for key in self.Analysis: keys += 1
        samples = len(self.Analysis[0]['mtf'])
        X = np.zeros((50, keys))
        Y = np.zeros((50, keys))
        Z = np.zeros((50, keys))
        i = 0
        for key in self.Analysis:
            X[:, key] = self.freqs[:50]
            Y[:, key] = np.ones(50) * i 
            Z[:, key] = self.Analysis[key]['mtf'][:50]
            i += 2

        ax.plot_wireframe(X, Y, Z)
        ax.set_xlabel('cycles/degree')
        ax.set_ylabel('degrees')
        plt.show()

        ## Retina plot ##
        #################

        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111, projection='3d')     
        pf.AxisFormat()

        keys = 0
        for key in self.Analysis: keys += 1
        samples = len(self.Analysis[0]['mtf'])
        X = np.zeros((50, keys))
        Y = np.zeros((50, keys))
        Z = np.zeros((50, keys))
        i = 0
        for key in self.Analysis:
            X[:, key] = self.freqs[:50]
            Y[:, key] = np.ones(50) * i 
            contrast = (self.Analysis[key]['retina'][:50] / 
                    np.max(self.Analysis[key]['retina'][:50]))
            Z[:, key] = sig.decibels(contrast)
            i += 2

        ax.plot_wireframe(X, Y, Z)
        ax.set_xlabel('cycles/degree')
        ax.set_ylabel('degrees')
        plt.show()

        
    def plotInformation(self, save_plots=False, legend=False):
        """
        :param cones: list of cones in analysis
        :type cones: list
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
            
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)
    
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], [5,5])
        
        for key in self.Analysis:
            plt.plot(self.Analysis[key]['cones'], self.Analysis[key]['info'],
                     self.Analysis[key]['line']['style'],
                      c = self.Analysis[key]['line']['color'], label=key,
                     linewidth=2.5, markersize=9)
    
    
        plt.xlim([min(self.Analysis[key]['cones'])-0.1, 
                  max(self.Analysis[key]['cones'])+0.1])
        if legend:
            plt.legend(loc='upper left')
        plt.xlabel('cones in mosaic')
        plt.ylabel('entropy (bits)')
        
        plt.tight_layout()
        
        if save_plots:
            fig.show()
            fig.savefig(self.figPath + 'InfoPlot.png')
            plt.close()
        else:
            plt.show()
    
    def plotDoG(self, rec_field, save_plots=False):
    
        """
        Currently this function outputs the following:
    
        :param save_plots: decide whether to save plots (True) or not (False)
        :type save_plots: bool
    
        :returns: * Plot1: difference of gaussians receptive field \n
                  * Plot2: FFT spectrum of plot 1.
                
        .. todo::
           Get an analytical function from Curcio to find cone spacing
           as a function of eccentricity.
            
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
        
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)
    
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], [5,5])
        
        ax.plot(self.rec_field['xvals'], self.rec_field['dog'], 
                'k', linewidth=2.5)
    
        plt.xlabel('distance (arcmin)')
        plt.ylabel('amplitude')
        
        plt.tight_layout()
        
        if save_plots:
            fig.show()
            fig.savefig(self.figPath + 'ConeRF' + loc + '.png')
            plt.close()
        else:
            plt.show()
        
        fig2 = plt.figure(figsize=(8,6))
        ax = fig2.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], [5, 5])
        
        length = self.rec_field['length']
        
        normFFT = (self.rec_field['coneResponse']['fft'][length:] / 
            np.max(self.rec_field['coneResponse']['fft'][length:]))
        ax.semilogx(self.rec_field['xvals'][length:] * 60, 
                  sig.decibels(normFFT), 
                    'k', linewidth=2.5)
    
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
    
        plt.ylim([self.min_dB, 0])
        plt.xlim([self.freqs[1], 100])
        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('contrast sensitivity (dB)')
        
        plt.tight_layout()
        
        if save_plots:
            fig2.show()
            fig2.savefig(self.figPath + 'ConeRF_FFT' + loc + '.png')
            plt.close()
        else:
            plt.show()

            
    def plotAmpSpec(self, _brownian=True, save_plots=False):
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

        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], [5, 3])
        
        plaw = self.imageData['powerlaw'](self.imageData['imagexval'][10:300])
        if _brownian:
            from emmetrop.eye.movement import brownian_motion
            
            temp = np.arange(1, 80, 1)
            spat = self.imageData['imagexval'][10:300]
            movement_filter = brownian_motion(spat, temp)
            #ax.semilogx(self.imageData['imagexval'][10:300], movement_filter, 
            #          'g', linewidth = 2.5)           
            whiteStim = plaw * movement_filter
            whiteStim = sig.decibels(whiteStim)
            ax.semilogx(self.imageData['imagexval'][10:300], whiteStim, 
                      'k', linewidth = 2.5)     
 
        plaw = sig.decibels(plaw)
        

    
        ax.semilogx(self.imageData['imagexval'][10:300], 
                    self.imageData['decibels'][10:300],
                    'r.', markersize = 10, markeredgewidth = 0)
     
        ax.semilogx(self.imageData['imagexval'][10:300],
                    plaw, 
                  'k', linewidth = 2.5)
                  
        ax.set_xlim([0.45, 15])
                      
        #ax.text(1, 10**-2.4, r'$\frac{1}{\mathit{f}}$', size = 35)

        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('spectral density (dB)')
        
        plt.tight_layout()
        
        if save_plots:
            fig.show()
            fig.savefig(self.figPath + 'ampSpec.png')
            plt.close()
        else:
            plt.show()
    
    
    def plotMTFfamily(self, save_plots = False, legend = False):    
        """Plot a family of MTF curves derived from a schematic eye.
        
        :param save_plots: decide whether to save plots (True) or not (False).
        :type save_plots: bool
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
        
        for key in self.Analysis:
            mtf = self.Analysis[key]['mtf']
            ax.plot(self.freqs, mtf,
                      self.Analysis[key]['line']['style'],
                      c = self.Analysis[key]['line']['color'],)
        ax.plot(self.freqs, self.diffract['mtf'], 'k-', linewidth=2)

        if legend: 
            ax.legend(loc='upper right')#title='object, retina')

        #plt.ylim([self.min_dB, 1])
        plt.xlim([0, 60])
            
        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('modulation transfer')
        
        plt.tight_layout()
        
        if save_plots:
            fig.show()
            fig.savefig(self.figPath + 'MTFfamily.png')
            plt.close()
        else:
            plt.show()
            
    def plotComparisonMTF(self, save_plots=False, legend=False):
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
        from emmetrop.eye.Optics import MTF
        from emmetrop.eye.eyeModel import SchematicEye as e
        
        Fovea = MTF(self.freqs, 0)
        TenDeg = MTF(self.freqs, 10)
        TwentyDeg = MTF(self.freqs,20)
        FourtyDeg = MTF(self.freqs,40)
        
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)

        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], [5,5])
        
        #Navarro et al 1993 analytical func:
        ax.plot(self.freqs, Fovea, 'm--', label='fovea', linewidth=2.5)
        ax.plot(self.freqs, TenDeg, 'r--', label='10 deg', linewidth=2.5)
        ax.plot(self.freqs, TwentyDeg, 'g--', label='20 deg', linewidth=2.5)
        ax.plot(self.freqs, FourtyDeg, 'b--',label='40 deg', linewidth=2.5)
        
        #OSLO ray trace data:
        eye = e()
        intensity = eye.traceEye(1e8, 0, 4, 0)
        mtf = eye.genMTF(intensity)
        ax.plot(self.freqs, mtf, 'm-', label='fovea')

        intensity = eye.traceEye(1e8, 10, 4, 0)
        mtf = eye.genMTF(intensity)
        ax.plot(self.freqs, mtf, 'r-', label='10 deg')  

        intensity = eye.traceEye(1e8, 20, 4, 0)
        mtf = eye.genMTF(intensity)      
        ax.plot(self.freqs, mtf, 'g-', label='20 deg')

        intensity = eye.traceEye(1e8, 40, 4, 0)
        mtf = eye.genMTF(intensity)
        ax.plot(self.freqs, mtf, 'b-', label='40 deg')
                
        if legend: 
            ax.legend(loc='lower left')#,title='object dist, retinal location')
        
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        
        mi, ma = plt.ylim()
        
        #plt.ylim([10**-2.5, 10**0])
        plt.xlim([0, 60])
        
        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('modulation transfer')
        
        plt.tight_layout()
        
        if save_plots:
            fig.show()
            fig.savefig(self.figPath + 'MTFperiphery.png')
            plt.close()
        else:
            plt.show()
            
    def plotActivity(self, save_plots = False, legend = False):
        """Plot powerlaw amplitude spectrum after accounting for MTF and MTF \
        + receptive field
               
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
            
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)

        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], [5, 5])

        for key in self.Analysis:
            contrast = (self.Analysis[key]['preCone'] / 
                    np.max(self.Analysis[key]['preCone']))

            ax.semilogx(self.freqs[0:100], sig.decibels(contrast),
                      self.Analysis[key]['line']['style'],
                      c = self.Analysis[key]['line']['color'], label=key)
        contrast = self.diffract['preCone'] / np.max(self.diffract['preCone'])
        ax.plot(self.freqs[0:100], sig.decibels(contrast),
            'k-', linewidth=2)

        if legend: 
            ax.legend(loc='lower left')
        
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
       
        plt.ylim([self.min_dB, 0])
        plt.xlim([0, 60])
        
        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('contrast sensitivity (dB)')
        
        plt.tight_layout()
        
        if save_plots:
            fig.show()
            fig.savefig(self.figPath + 'MTFfamilyMod.png')
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
        pf.TufteAxis(ax, ['left', 'bottom'], [5, 5])

        for key in self.Analysis:
            contrast = (self.Analysis[key]['retina'] / 
                    np.max(self.Analysis[key]['retina']))
            ax.semilogx(self.freqs[0:100], sig.decibels(contrast),
                      self.Analysis[key]['line']['style'],
                      c = self.Analysis[key]['line']['color'], label=key)
        contrast = self.diffract['retina'] / np.max(self.diffract['retina'])                
        ax.plot(self.freqs[0:100], sig.decibels(contrast),
                    'k-', linewidth=2)

        if legend: 
            ax.legend(loc='upper right')
        
        mi, ma = plt.ylim()
   
        plt.ylim([self.min_dB, 0])
        plt.xlim([0, 60])
        
        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('contrast sensitivity (dB)')
        
        plt.tight_layout()
       
        if save_plots:
            save_name = self.figPath + 'MTFfamilyModAmp.png'
            
            fig2.show()
            fig2.savefig(save_name)     
            plt.close()
            
        else:
            plt.show()