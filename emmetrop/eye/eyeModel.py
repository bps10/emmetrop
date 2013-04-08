import os
import shlex
import numpy as np

from emmetrop.scene.DataManip import rad2deg

class SchematicEye(object):
    """
    This class controls the loading of OSLO data.
    Eventually it will interface with a C++ ray tracer dynamically.
    """
    
    def __init__(self, OSLO_directory=None):
        """
        
        """
        if not OSLO_directory:
            if os.path.basename(os.getcwd()) == 'emmetrop':
                p = './OSLO_MTF_DATA/'
            else:
                p = './emmetrop/OSLO_MTF_DATA/'
        else:
            p = OSLO_directory
        self.loadOSLOData(p)

        
    def loadOSLOData(self,p):
        """Load data from OSLO
        
        This is currently a very static function.
        
        :param p: directory path. Passed from init function.
        :type p: path
        
        """
        
        self.INF = self.importOSLOfile(p + 
                                'ONaxisMTFinfFocusNavarrow1999.txt')
        self.TwentyFt = self.importOSLOfile(p +
                                'ONaxisMTF20ftFocusNavarrow1999.txt')
        self.Onemeter = self.importOSLOfile(p + 
                                'ONaxisMTF1mFocusNavarrow1999.txt')
        self.SixteenIn = self.importOSLOfile(p + 
                                'ONaxisMTF16inFocusNavarrow1999.txt')
        
        self.INF_offaxis = self.importOSLOfile(p + 
                                'OFFaxisMTFinfFocusNavarrow1999.txt')
        self.TwentyFt_offaxis = self.importOSLOfile(p + 
                                'OFFaxisMTF20ftFocusNavarrow1999.txt')
        self.Onemeter_offaxis = self.importOSLOfile(p +
                                'OFFaxisMTF1mFocusNavarrow1999.txt')
        self.SixteenIn_offaxis = self.importOSLOfile(p + 
                                'OFFaxisMTF16inFocusNavarrow1999.txt')
        
        self.Sixteen_UnderAccomm = self.importOSLOfile(p +
                            'OFFaxisMTF16inUnderAccom20ftObjNavarrow1999.txt')
        self.Sixteen_SixteenObj_Offaxis = self.importOSLOfile(p +
                                'OFFaxisMTF16inFocus16inObjNavarrow1999.txt')        
        self.Sixteen_TwentyObj_Offaxis = self.importOSLOfile(p + 
                                'OFFaxisMTF16inFocus20ftObjNavarrow1999.txt')
    
        self.TwentyDegOffAxis_InfFoc = self.importOSLOfile(p + 
                                '20degOFFaxisMTFinfFocusNavarrow1999.txt')
        self.FortyDegOffAxis_InfFoc = self.importOSLOfile(p + 
                                '40degOFFaxisMTFinfFocusNavarrow1999.txt')       
 
        ## convert mm to rad (1mm image/24mm axial length)
        self.freqs = self.INF[:,1] / rad2deg(1.0/self.getAxialLength())
        
        
        self.dataPackage = {
                            'onAxis': {'diffract': self.INF[:,4],
                                       'inf': self.INF[:,2], 
                                       '20ft': self.TwentyFt[:,2],
                                       '1m': self.Onemeter[:,2],
                                       '16in': self.SixteenIn[:,2]},
                                   
                            'offAxis': {'diffract': self.INF_offaxis[:,4],
                                        'inf':self.INF_offaxis[:,2],
                                        '20ft': self.TwentyFt_offaxis[:,2],
                                        '1m': self.TwentyFt_offaxis[:,2],
                                        '16in': self.SixteenIn_offaxis[:,2]},
                                        
                            'object': {
                            '16inunder': self.Sixteen_UnderAccomm[:,2],
                            '16in16in': self.Sixteen_SixteenObj_Offaxis[:,2],
                            '16in20ft': self.Sixteen_TwentyObj_Offaxis[:,2]},
                                       
                            'farPeriph': {
                                    '20deg': self.TwentyDegOffAxis_InfFoc[:,2],
                                    '40deg': self.FortyDegOffAxis_InfFoc[:,2]},
                                       
                            'freqs': self.freqs
                            }


    def returnOSLOdata(self):
        """
        Return a dictionary of imported transfer functions from OSLO.
        """
        return self.dataPackage
        

    
    def importOSLOfile(self, OSLOfile):
        """Import a text file with MTF output from OSLO.
        
        :param OSLOfile: name of OSLO output text file to import. \
        Should have 5 columns.
        :type OSLOfile: string
        :returns: an array containing OSLO data: Frequencies, \
        Eye MTF, Diffraction limit MTF.
        :rtype: numpy.array
        
        .. todo::
           import OSLO files, import whole directory into a dictionary.
        
        """
        
        fil = open(OSLOfile)
        fil = fil.read()
        
        MTF = []
        foo = fil
        partitioning = True
        row = 0
        while partitioning:
            
            f = foo.partition('\n')
    
            if f[0]:
                
                parse = shlex.split(f[0])
                floatnum = []
                for num in parse:
                    if row == 0:
                        if num == '--':
                            floatnum.append( 0.0 )
                    if num != '--':
                        floatnum.append( float(num) )
    
                        
                MTF.append(floatnum)
                foo = f[2]
                partitioning = True
            else:
                partitioning = False
        
        MTF = np.array(MTF)
        
        return MTF        

    def getAxialLength(self):
        """Find the axial length in mm of the optical system used to generate 
        outputed curves.
        
        .. warning::
           Currently hard coded to return 24mm
           
        .. todo::
           Find axial length dynamically. With C++ can simply call getAxialLength 
           function.
        """
        return 24.0