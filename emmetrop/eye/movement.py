# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np


def brownian_motion(spat_freq, temp_freq):
    '''
    '''
    D = 40 # archmin ** 2 / min
    sqr_spat_freq = spat_freq ** 2
    Q = np.zeros(len(spat_freq))
    for t_freq in temp_freq:
        Q += sqr_spat_freq * D / ((sqr_spat_freq) ** 2 * (D ** 2 / 4) + 
                t_freq ** 2)
    return Q
    

if __name__ == "__main__":
    import matplotlib.pylab as plt
    
    temp = np.arange(1, 751)
    spat = np.arange(0.1, 50, 0.1)
    plt.figure()
    plt.loglog(spat, brownian_motion(spat, temp))
    plt.show()