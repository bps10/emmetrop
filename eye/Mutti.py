# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 13:46:43 2012

@author: Brian
"""

#import numpy as np

DTE = float(raw_input('Enter DTE : '))
LP = float(raw_input('Enter Lens power : '))
REcornea = float(raw_input('Enter RE cornea : '))

out = ( 1.0 / ((1.0 / ((1.0 / (0.013 - DTE)) + LP )) - 0.13)) - REcornea


print out