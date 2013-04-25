# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pylab as plt
import PlottingFun as pf

def plotAccommodation():
    """
    """    
    
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax, ['left', 'bottom'], [5,4])    
    x = np.arange(5, 21, 1)
    
    y300 = lambda x: 1.07 - (0.02 * x)
    y392 = lambda x: 1.49 - (0.05 * x)
    y480 = lambda x: 1.93 - (0.07 * x)
    y567 = lambda x: 1.83 - (0.05 * x)
    y652 = lambda x: 1.69 - (0.03 * x)
    
    #ax.plot(x, y300(x), 'k', label='3.00', linewidth=2.5)  
    ax.plot(x, y392(x), 'k-', label='3.92', linewidth=2.5) 
    ax.plot(x, y480(x), 'k:', label='4.80', linewidth=2.5) 
    ax.plot(x, y567(x), 'k--', label='5.67', linewidth=2.5) 
    #ax.plot(x, y652(x), 'k', label='6.52', linewidth=2.5) 
    
    ax.legend(title="demand (D)")
    ax.set_ylim([0.45, 1.75])
    ax.set_ylabel('accommodative lag (D)')
    ax.set_xlabel('age')
    '''    
    ax.set_title('4.80D demand')    
    ax.text(0.85, 0.95, 
        ('y = 1.93 - 0.07x'), 
        fontsize=18, 
        horizontalalignment='center',
        verticalalignment='top',
        transform=ax.transAxes)
    '''
    plt.tight_layout()
    plt.show()
    
if __name__ == '__main__':
    plotAccommodation()
    