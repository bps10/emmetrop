from __future__ import division
import numpy as np
import pyximport
pyximport.install(setup_args={"script_args":["--compiler=mingw32"],
                              "include_dirs":np.get_include()},
                  reload_support=True)

"""                  
cdef double computeImageHist(double img):
    Compute the histogram of a grayscale (for now) image
    
    :param img: input image.
    :type img: np.array

    DTYPE = np.float
    image_shape = img.shape
    BINS = np.linspace(0,4 * 10**4, BIN_RESOLUTION + 1)
    
    #preallocate memory
"""