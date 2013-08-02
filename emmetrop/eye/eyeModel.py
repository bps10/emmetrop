import sys, os
import numpy as np

if os.path.basename(os.getcwd()) == 'emmetrop':
    sys.path.append("./SchematicEye")
else:
    sys.path.append("./emmetrop/SchematicEye")
import eye as eye
     

def traceEye(object_distance=1e8, off_axis=0, pupil_size=3, diopters=0, wavelength=550):
    '''
    '''

    intensity = eye.py_eye(object_distance, off_axis, pupil_size, diopters, wavelength)

    return intensity
