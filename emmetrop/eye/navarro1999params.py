from __future__ import division
import numpy as np


def NavarroParams():
    """
    This function implements the equations from Navarro's 1985 paper and the
    follow up wide angle model from 1999.
    
    It takes as input a command line entered float indicating accommodation
    in diopters.  
    
    :returns: prints parameters.
    
    .. note:: 
       using this right now to change OSLO parameters.
    """
    print ' '
    accom = float(raw_input("enter accommodation in diopters:  "))
    print ' '
    
    ant_lens_radius     = 10.2 - ( 1.75 * np.log(accom + 1.0) )
    post_lens_radius    = -6 + ( 0.2294 * np.log(accom + 1.0) )
    aqueous_thickness   = 3.05 - ( 0.05 * np.log(accom + 1.0) )
    lens_thickness      = 4.0 + ( 0.1 * np.log(accom + 1.0) )
    lens_refractive_index = 1.42 + ( 9.0 * 10**-5 * (10.0 * accom + accom**2) )
    anterior_lens_asphericity = -3.1316 - ( 0.34 * np.log(accom + 1.0) )
    posterior_lens_asphericity = -1.0 - ( 0.125 * np.log(accom + 1.0) )
    
    print 'anterior lens radius : ', ant_lens_radius
    print 'posterior lens radius: ', post_lens_radius
    print 'aqueous thickness    : ', aqueous_thickness
    print 'lens thickness       : ', lens_thickness
    print 'lens refractvie index: ', lens_refractive_index
    print 'anterior lens aspher : ', anterior_lens_asphericity
    print 'posterior lens aspher: ', posterior_lens_asphericity

