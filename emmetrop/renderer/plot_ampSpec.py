from __future__ import division
import matplotlib.pylab as plt
import numpy as np

from base import plot as pf
from emmetrop.scene import SignalProcessing as sig

def plotAmpSpec(imageData, figPath='Figures/', _brownian=True, save_plots=False):
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
    
    plaw = imageData['powerlaw'](imageData['imagexval'][10:300])
    if _brownian:
        from emmetrop.eye.movement import brownian_motion
        
        temp = np.arange(1, 80, 1)
        spat = imageData['imagexval'][10:300]
        movement_filter = brownian_motion(spat, temp)
        #ax.semilogx(imageData['imagexval'][10:300], movement_filter, 
        #          'g', linewidth = 2.5)           
        whiteStim = plaw * movement_filter
        whiteStim = sig.decibels(whiteStim)
        ax.semilogx(imageData['imagexval'][10:300], whiteStim, 
                  'k', linewidth = 2.5)     

    plaw = sig.decibels(plaw)
    


    ax.semilogx(imageData['imagexval'][10:300], 
                imageData['decibels'][10:300],
                'r.', markersize = 10, markeredgewidth = 0)
 
    ax.semilogx(imageData['imagexval'][10:300],
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