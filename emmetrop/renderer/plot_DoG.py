from __future__ import division
import matplotlib.pylab as plt
import numpy as np

from base import plot as pf
from emmetrop.scene import SignalProcessing as sig


def plotDoG(rec_field, min_dB=20, figPath='Figures', save_plots=False):

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
    
    ax.plot(rec_field['xvals'], rec_field['dog'], 
            'k', linewidth=2.5)

    plt.xlabel('distance (arcmin)')
    plt.ylabel('amplitude')
    
    plt.tight_layout()
    
    if save_plots:
        fig.show()
        fig.savefig(figPath + 'ConeRF' + loc + '.png')
        plt.close()
    else:
        plt.show()
    
    fig2 = plt.figure(figsize=(8,6))
    ax = fig2.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], [5, 5])
    
    length = rec_field['length']
    
    normFFT = (rec_field['coneResponse']['fft'][length:] / 
        np.max(rec_field['coneResponse']['fft'][length:]))
    ax.semilogx(rec_field['xvals'][length:] * 60, 
              sig.decibels(normFFT), 
                'k', linewidth=2.5)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    plt.ylim([min_dB, 0])
    plt.xlim([0, 100])
    plt.xlabel('spatial frequency (cycles / deg)')
    plt.ylabel('contrast sensitivity (dB)')
    
    plt.tight_layout()
    
    if save_plots:
        fig2.show()
        fig2.savefig(figPath + 'ConeRF_FFT' + loc + '.png')
        plt.close()
    else:
        plt.show()