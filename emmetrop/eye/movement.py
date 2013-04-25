# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np


def brownian_motion(spat_freq, temp_freq):
    '''
    '''
    D = 40 / 60 ** 2 # arcmin ** 2 / min
    # convert freq into angular freq
    #temp_freq *= 2.0 * np.pi

    sqr_spat_freq = (spat_freq) ** 2
    Q = np.zeros(len(spat_freq))
    for t_freq in temp_freq:
        Q += sqr_spat_freq * D / ((sqr_spat_freq ** 2) * (D ** 2.0 / 4.0) + 
                t_freq ** 2.0)
    Q = Q / np.max(Q)
    return Q
    

if __name__ == "__main__":
    import matplotlib.pylab as plt
    
    temp = np.arange(1, 80)
    spat = np.arange(0.1, 50, 0.1)
    plt.figure()
    plt.loglog(spat, brownian_motion(spat, temp))
    plt.show()
    
    '''
    import emmetrop.scene.DataManip as dm
    # need to weight spat freq by occurances in image:    
    xs = 1500
    f = np.arange(-xs / 2.0, xs / 2.0)

    XX, YY = np.meshgrid(f,f)
    foo, r = dm.cart2pol(XX,YY)
    if np.mod(xs,2)==1:
        r = np.around(r)-1
    else:
        r = np.around(r)
    r = np.array(r,dtype=int)
    tot = dm.accum(r.flatten() + 1)[1:]
    tot[:(xs / 2 + 1)]
    '''