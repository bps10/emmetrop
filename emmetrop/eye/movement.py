#! /usr/bin/env python
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
    
    