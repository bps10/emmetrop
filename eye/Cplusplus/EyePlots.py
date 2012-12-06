from __future__ import division
import numpy as np
import matplotlib.pylab as plt
import sys

sys.path.append('C:\\Users\\Brian\\Documents\\Neitz-Lab\\Python')
import PlottingFun
            
            
class EyePlot(PlottingFun.PlottingFun):

    def __init__(self):
        self.data = {}
        self.Intensity = {}

    def fitGaussian(self, data):

        #gaussian = lambda x: 3*np.exp(-(30-x)**2/20.)
        
        #data = gaussian(np.arange(100))
        
        X = np.arange(data.size)
        x = sum(X*data)/sum(data)
        width = np.sqrt(abs(sum((X-x)**2*data)/sum(data)))
        
        MAX = data.max()
        
        fit = lambda t : MAX*np.exp(-(t-x)**2/(2*width**2))
        
        return fit(X)

    def smooth(self, x,window_len=11,window='hanning'):
        """smooth the data using a window with requested size.
        from : http://www.scipy.org/Cookbook/SignalSmooth
        
        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal 
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.
        
        input:
            x: the input signal 
            window_len: the dimension of the smoothing window; should be an odd integer
            window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
                flat window will produce a moving average smoothing.
    
        output:
            the smoothed signal
            
        example:
    
        t=linspace(-2,2,0.1)
        x=sin(t)+randn(len(t))*0.1
        y=smooth(x)
        
        see also: 
        
        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
        scipy.signal.lfilter
     
        TODO: the window parameter could be the window itself if an array instead of a string
        NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
        """
    
        if x.ndim != 1:
            raise ValueError, "smooth only accepts 1 dimension arrays."
    
        if x.size < window_len:
            raise ValueError, "Input vector needs to be bigger than window size."
    
    
        if window_len<3:
            return x
    
    
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
            raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    
    
        s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
        #print(len(s))
        if window == 'flat': #moving average
            w=np.ones(window_len,'d')
        else:
            w=eval('np.'+window+'(window_len)')
    
        y=np.convolve(w/w.sum(),s,mode='valid')
        return y



        
    ## Data Manipulation Methods ##
    ###############################
    
    def FormatData(self, savedata = 'Yes'):

        try:
            self.data = np.load('EyeDataProcessed.npz')
            print 'Successfully loaded previously processed data.'
        except IOError:
            self.data = []
            print 'No processed data found.  Now trying to generate fresh.'
        if self.data == []:
        
            self.data = np.genfromtxt('EyeLSA.csv', delimiter = ',')
            
            self.data = {
                        'lens_accomm':self.data[:,0], 
                        'pupil_size': self.data[:,1], 
                        'age': self.data[:,2], 
                        'accomPower': self.data[:,3], 
                        'relaxPower': self.data[:,4],
                        'defocus': self.data[:,5]
                        }
                    
            if savedata == 'Yes':
                np.savez('EyeDataProcessed.npz',
                        lens_accomm     = self.data['lens_accomm'],
                        pupil_size      = self.data['pupil_size'],
                        age             = self.data['age'],
                        accomPower      = self.data['accomPower'],
                        relaxPower      = self.data['relaxPower'],
                        defocus         = self.data['defocus'])
            
        try:
            self.Intensity = np.load('IntensityDataProcessed.npz')
            print 'Successfully loaded previously processed intensity data.'
        except IOError:
            pass
        if self.Intensity == {}:
                 
                 dat = np.genfromtxt('IntensityBestFocus3mm.csv', delimiter = ',')
                 
                 self.Intensity['downsample'] = 1.0
                 dat = dat[::int(self.Intensity['downsample'])]

                 samples = int( dat.shape[0] / 4.0 )

                 # Preallocate memory:
                 self.Intensity['rawintensity']     = np.zeros((4,samples)) 
                 self.Intensity['intensity']        = np.zeros((4,samples))
                 self.Intensity['pupil_size']       = np.zeros((4,samples))
                 self.Intensity['lens_accom']       = np.zeros((4,samples))
                 self.Intensity['PSF']              = np.zeros((4,samples))
                 self.Intensity['PSFtotal']         = np.zeros((4,(samples * 2)))
                 self.Intensity['PSFkernel']        = np.zeros((4,(samples * 2)))
                 self.Intensity['MTF']              = np.zeros((4,int(samples/2.0) ))
                 self.Intensity['smoothed']         = np.zeros((4,samples))

                 for i in range(0,4):

                     totals = np.max(dat[i*samples:(i+1)*samples,0])
                     self.Intensity['rawintensity'][i,:] = dat[i*samples:(i+1)*samples,0]
                     self.Intensity['intensity'][i,:]    = dat[i*samples:(i+1)*samples,0] / totals
                     self.Intensity['pupil_size'][i,:]   = dat[i*samples:(i+1)*samples,1]
                     self.Intensity['lens_accom'][i,:]   = dat[i*samples:(i+1)*samples,2]
                       
                 deriv = np.zeros((4,samples))
                 deriv[:,0]     = self.Intensity['intensity'][:,0]
                 deriv[:,1:] = self.Intensity['intensity'][:,1:] - self.Intensity['intensity'][:, 0:-1]
                                                
                                                
                 self.Intensity['PSF'][:,0] = 1.0
                 for i in range(1,samples-1):
                     self.Intensity['PSF'][:,i] = self.Intensity['PSF'][:,i-1] - deriv[:,i]
                                 
                 self.Intensity['PSFtotal'][:,1:samples+1] = self.Intensity['PSF'][:,::-1]
                 self.Intensity['PSFtotal'][:,samples+1:]  = self.Intensity['PSF'][:,1:]

                 self.Intensity['trial'] = ([],[],[],[])
                 for i in range(0,samples):
                     for j in range(0,4):
                         
                         self.Intensity['trial'][j].append(self.Intensity['PSF'][j,i])                    
                 self.Intensity['trial'] = np.array(self.Intensity['trial'])   
                 for i in range(0,4):
                     foo = self.Intensity['PSFtotal'][i,:]
                     self.Intensity['PSFkernel'][i,:] = self.fitGaussian(foo)
                     

                     #self.Intensity['PSFkernel'][i,:] = (self.Intensity['PSFkernel'][i,:] / 
                                                    #np.max(self.Intensity['PSFkernel'][i,:]))                     

                     self.Intensity['MTF'][i,:] = np.abs(np.fft.fftshift(
                     np.fft.fft(self.Intensity['PSF'][i,:])))[int(samples/2.0):]
                     
                     self.Intensity['MTF'][i,:] = (self.Intensity['MTF'][i] / 
                                                    np.max(self.Intensity['MTF'][i]))

                     """
                     self.Intensity['MTF'][i,:] = self.smooth(self.Intensity['MTF'][i,:], 
                                                            window_len = 11)[:-10]
                     """
                    
                    
    def FindPlottingData(self, arg1, type1, arg2, type2, arg3=None, type3=None):
        
        ind1 = self.data[type1] == arg1
        ind2 = self.data[type2] == arg2
        
        if arg3 == None:
            self.loc = ind1 & ind2 == True
        else:
            ind3 = self.data[type3] == arg3
            self.loc = ind1 & ind2 & ind3 == True
    
    ## Plotting Methods ##
    ######################
    
    def PowerPlots(self, acc = [0,2,4,6,8], age = [10, 20]):
        
        fig = plt.figure()
        self.AxisFormat(20)
        
        ax = fig.add_subplot(111)
        linestyle = ['ko-','ko--','ko-.','koD']
        for j in range(0, len(age)):
            self.x1dat, self.y1dat, self.topaxis = [],[],[]
            for i in range(0, len(acc)):
                
                self.FindPlottingData(acc[i],'lens_accomm',age[j], 'age', 5, 'pupil_size')
                self.x1dat = np.append(self.x1dat, self.data['lens_accomm'][self.loc])
                self.y1dat = np.append(self.y1dat, self.data['accomPower'][self.loc])
                self.topaxis = np.append(self.topaxis, self.data['defocus'][self.loc])
            
            
            ax.plot(self.x1dat,self.y1dat, linestyle[j], linewidth=2.5, ms=10, 
                    label = '{0} years'.format(age[j]))

        ax.legend(loc = 'upper left')
        
        #ax2 = ax.twiny()
        #ax2.set_xticks(self.x1dat)
        #ax2.set_xticklocations(self.topaxis)
        #ax.set_xticks(self.topaxis)
        ax = self.TufteAxis(ax, ['left', 'bottom'], [5,5], integer='on')
        
        
        plt.xlim([min(acc)-0.1, max(acc)+0.1])
        plt.ylabel('power (D)')
        plt.xlabel('accommodation (D)')
        plt.tight_layout()
        plt.show()
        
    def EncircledIntensityPlot(self, plotRange = 400):
         fig = plt.figure()
         self.AxisFormat(20)
         ax = fig.add_subplot(111)
         deg = 1 / (24*2*np.pi) * 360.0
         x = np.linspace(0,deg,plotRange) * (plotRange / 
         self.Intensity['intensity'].shape[1]) * self.Intensity['downsample']
         for i in range(0,self.Intensity['intensity'].shape[0]):
             ax.plot(x, self.Intensity['intensity'][i,:plotRange].T, linewidth=2.5,
                     label = '{0} D'.format(self.Intensity['lens_accom'][i,0]))
         ax.legend(loc = 'lower right').set_title('lens accommodation')
         ax = self.TufteAxis(ax, ['left', 'bottom'], [5,5], integer='off')
         plt.ylim([0, 1.05])
         plt.ylabel('encircled intensity')
         plt.xlabel('degrees')
         plt.tight_layout()
         plt.show()


        
    def PSFplot(self, plotRange = 400):
         fig = plt.figure()
         self.AxisFormat(20)
         ax = fig.add_subplot(111)
         size = self.Intensity['PSFtotal'].shape[1] / 2.0
         
         
         deg = 1 / (24*2*np.pi) * 360.0
         x = np.linspace(0,deg, plotRange) * (plotRange / 
         self.Intensity['PSFtotal'].shape[1]) * self.Intensity['downsample']

         for i in range(0,self.Intensity['PSFtotal'].shape[0]):
             ax.plot(x, self.Intensity['PSFtotal'][i,size:size+plotRange].T, linewidth=2.5,
                    label = '{0} D'.format(self.Intensity['lens_accom'][i,0]))
         ax.legend(loc = 'upper right').set_title('lens accommodation')
         ax = self.TufteAxis(ax, ['left', 'bottom'], [5,5], integer='off')
         plt.ylim([-0.05, 1.01])
         plt.ylabel('PSF')
         plt.xlabel('degrees')
         plt.tight_layout()
         plt.show()
        
    
    def MTFplot(self):
         fig = plt.figure()
         self.AxisFormat(20)
         ax = fig.add_subplot(111)

         ax.spines['top'].set_visible(False)
         ax.spines['right'].set_visible(False)
         ax.spines['left'].set_position(('outward',10))
         #ax.spines['left'].set_smart_bounds(True)
         ax.spines['bottom'].set_position(('outward',10))
         ax.spines['bottom'].set_smart_bounds(True)
         ax.get_xaxis().tick_bottom()
         ax.get_yaxis().tick_left()

         cycles = (np.arange(0,self.Intensity['MTF'].shape[1]) )#* self.Intensity['downsample'])
         deg = 1 / (24*2*np.pi) * 360.0
         cpd = cycles / deg
         for i in range(0,self.Intensity['MTF'].shape[0]):
             ax.plot(cpd[:], self.Intensity['MTF'][i,:].T, linewidth=2.5,
                         label = '{0} D'.format(self.Intensity['lens_accom'][i,0]))
                         
         ax.legend(loc = 'upper right').set_title('lens accommodation')

         plt.ylim([0, 1.0])
         plt.ylabel('MTF')
         plt.xlabel('cycles / deg')
         plt.tight_layout()
         plt.show()
        
        
    def LSA(self, acc = [0,6], age = [13, 20] ):
        """
        """
        
        fig = plt.figure()
        self.AxisFormat(20)

        ax1 = fig.add_subplot(111)
        
        for j in range(0,len(age)):
            for i in range(0,len(acc)):

                self.FindPlottingData(acc[j],'lens_accomm', age[i], 'age')
                
                self.x1dat = self.data['defocus'][self.loc]
                self.y1dat = self.data['pupil_size'][self.loc]
                
                ax1.plot(self.x1dat,self.y1dat, linewidth=2.5,
                        label='{0} D, {1} y'.format(acc[i], age[j]));
                
        
        ax1.legend(loc = 'lower right').set_title('Accommodation')
        ax1 = self.TufteAxis(ax1,['left','bottom'], [5,5])

        plt.ylabel('pupil size')
        plt.xlabel('diopters')
        plt.tight_layout()
        plt.show()

        
if __name__ == "__main__":

    """
    Want to make the program run to ensure most up to date files are used.
    stdin,stdout = wp.popen4('C:\\Users\\Brian\\Documents\\SchematicEye\\tessar.exe')
    stdin.close()
    print 'pOpen done..'
    print stdout
    """
    # main function can change options
        
    out = EyePlot()
    out.FormatData()
    #out.LSA(acc = [2,6], age = [10, 15])
    #out.PowerPlots()
    #out.EncircledIntensityPlot()
    out.PSFplot()
    out.MTFplot()

