from __future__ import division
import numpy as np
import matplotlib.pylab as plt
from matplotlib.lines import Line2D
import matplotlib.animation as animation

#emmetrop imports
from emmetrop.renderer import PlottingFun as pf
from emmetrop.scene import DataManip as dm

## 1. Figure out relationship between spring constant and number of nodes.

def getParams():
    """
    calling this will get some user input and then create the movie.
    """
    
    choice = raw_input('Do you wish to change default parameters of the myopic eye (y / n) ?')
    
    if choice.lower() == 'y' or choice.lower() == 'yes':
        print 'Which options would you like to change:'
        print ''
        
        params = ['NumberofNodes','age','INTSTEP','LengthofTime','Pressure','K_change','LocationOfLoosen','LengthofLoosen']
        english = ['number of nodes', 'start age', 'integration step', 'length of simulation', 'pressure',
                    'spring constant', 'location of loosening', 'length of loosened section']
        param_range = [ '4-10', '1 - 20', '0.001 - 1', '1 - 20', '0.0001 - 0.1', '0.001 - 0.99', '0 - 10', '1 - 10']
        
        input_params = {}
        
        for i,name in enumerate(english):
            print '', i,'. ', name
        print ''
        
        input = True
        while input:
            while True:
                try:
                    selection_name = int(raw_input("Select parameter (an integer) : "))
                    variable = params[selection_name]
                    break
                except (ValueError,IndexError):
                    print 'Sorry input not understood. Please enter a valid integer.'
                
            selection_value = float(raw_input("Enter new value from {0}  : ".format(param_range[selection_name])))
            
            input_params[str(variable)] = selection_value
            
            print ''
            input_option = (raw_input("Would you like to change anything else (y / n) ?"))
            print ''
            
            if input_option.lower() == 'y' or input_option.lower() == 'yes':
                input = True
            else:
                print ''
                print 'Alright, than off we go.'
                input = False
        
    elif choice.lower() == 'n' or choice.lower() == 'no':
        print 'Alright, than off we go.' 
        input_params = {}
    
    else:
        print 'Woops.  Answer not understood.  Please try again.'
    
    return input_params
    
        
        
        
def EyeGrowthFunc(NumberofNodes=36,age=6,INTSTEP = 0.01, LengthofTime=12,Pressure=0.002,
    Loosen=0,K_change=0.15,LocationOfLoosen=0,LengthofLoosen=2,Animate=0):
    """
    
     Compute the change in node position given a set of input parameters
    
     :param NumberofNodes: Typical = 36.
     :type NumberofNodes: int 
     :param age: Age in years at start of simulation. Typical = 6.
     :type age: float
     :param LengthofTime: Duration of simulation in years. Typical = 12.
     :type LengthofTime: float
     :param INTSTEP: Integration step in fraction of year. Typical = 0.01.
     :type INTSTEP: float
     :param Pressure: Intraocular eye pressure.
     :type Pressure: float
     :param Loosen: Choose whether to loosen springs, 0 = no, 1 = yes.
     :type Loosen: int     
     :param K_change: If loosening, chose how much to change the spring constant,K by. 
     :type K_change: int     
     :param LocationOfLoosen: If loosening, choose location of loosened springs.
     :type LocationOfLoosen: int
     :param LengthofLoosen: If loosening,choose number of springs to loosen.          
     :type LenghtofLoosen: int
     :param Animate: Choose to plot a video of the nodes. Default = 0 (no).
     :type Animate: int
         
     :returns: Nodes structure containing info about node locations

     .. note:: 
        see Hung et al. 2010.
        
    .. todo:
       Move plotting methods to renderer.
    """

    # Generate Nodes as objects
    Nodes = GenerateNodes(NumberofNodes,age=age,Pressure = Pressure,LengthofTime = LengthofTime,INTSTEP = INTSTEP)
    # Generate Springs as objects
    Springs = GenerateSprings(Nodes,loc = LocationOfLoosen,loosening = Loosen,K_change = K_change,
        length = LengthofLoosen)
        
    ## Start of main loop
    
    for time in range(0,int(LengthofTime/INTSTEP)):

        if Nodes['age'][time] < 10.5:
            K = Springs['young']
        elif Nodes['age'][time] >= 10.5:
            K = Springs['old']
        
        Nodes = ComputeSpringForce(Nodes,Springs['spring_eq'],K,time)
        
        Nodes = ComputePressureForce(Nodes,time)
        
        Nodes['xx'] = Nodes['f_px'] - Nodes['f_sx'] # force
        Nodes['yy'] = Nodes['f_py'] - Nodes['f_sy'] # force

        Nodes['xDOT'] += Nodes['xx'] * INTSTEP # update velocity
        Nodes['yDOT'] += Nodes['yy'] * INTSTEP # update velocity

        Nodes['x'][time+1] = Nodes['x'][time] + Nodes['xDOT'] * INTSTEP 
        Nodes['y'][time+1] = Nodes['y'][time] + Nodes['yDOT'] * INTSTEP 

        Nodes = CheckNodes(Nodes,time)
        
        Nodes['modeleyelength'][time] = (abs(min(Nodes['x'][time+1]) * Nodes['radius']) + 
                                        max(Nodes['x'][time+1]) * Nodes['radius'])
    
    if Animate == 1:
        xx = np.zeros(( (LengthofTime / INTSTEP) + 1, NumberofNodes + 1))
        yy = np.zeros(( (LengthofTime / INTSTEP) + 1, NumberofNodes + 1))
        for i in range(0,len(xx)):
        
            xx[i,:] = np.r_[ Nodes['x'][i],Nodes['x'][i,0] ] * Nodes['radius']
            yy[i,:] = np.r_[ Nodes['y'][i],Nodes['y'][i,0] ] * Nodes['radius']
            
        return Nodes, xx, yy
        
    else:    
        return Nodes



    
### Generate Functions ###
##########################
    
def GenerateSprings(Nodes,loc=0,loosening=0,K_change=0.15,length=2,Young_K=0.19,Old_K=0.025):
    """

    Generate a spring object
    
    :param Nodes: structure containing node information
    :type Nodes: struct
    :param loosening: Select whether to loosen any of the springs. Yes = 1, No = 0.
    :type loosening: int
    :param K_change: If loosening, select how much to change the spring constant, K.
    :type K_change: float
    :param loc: If loosening, select the region to loosen. Currently loosening is symmetric
    :type loc: int
    :param len: If loosening, select the number of springs to loosen.
    :type len: int
    
    :returns: spring structure with 'young' and 'old' spring constants and spring equalibrium.
    
    .. note::
       May eventually want to make the loosening area gaussianly distributed.
    """
    NumberofNodes = Nodes['x'][0].shape[0]
    Spring_Eq = np.ones((NumberofNodes)) * (np.max(CircularDist(Nodes['x'][0],Nodes['y'][0])))
    Springs =     {
                'young':np.ones((NumberofNodes)) * Young_K,
                'old': np.ones((NumberofNodes)) * Old_K,
                'spring_eq':Spring_Eq
                }                

    #set the location of weakened springs.
    #multiply K by specific locations of Spring to tighten or weaken it.
    if loosening == 1:

        Springs['young'][loc:loc+length] = Springs['young'][loc:loc+length]*K_change
        if loc == 0:
            Springs['young'][-loc-length:] = Springs['young'][-loc-length:]*K_change
        else:
            Springs['young'][-1-loc-length:-1-loc] = Springs['young'][-1-loc-length:-1-loc]*K_change

        Springs['old'][loc:loc+length] = Springs['old'][loc:loc+length]*K_change
        if loc == 0:
            Springs['old'][-loc-length:] = Springs['old'][-loc-length:]*K_change
        else:
            Springs['old'][-1-loc-length:-1-loc] = Springs['old'][-1-loc-length:-1-loc]*K_change
            
    return Springs


def GenerateNodes(NumberofNodes=36,age=6,Pressure=0.002,LengthofTime=12,INTSTEP=0.01):
    """

    Creates a structure of nodes based on the input parameters.

    :param NumberofNodes: the number of nodes desired.
    :type NumberofNodes: int
    :param age: age at the start of the simulation, in years.
    :type age: float
    :param Pressure: intraocular pressure.
    :type Pressure: float
    :param LengthofTime: duration of the simulation, in years.
    :type LengthofTime: int
    :param INTSTEP: integration step, in fraction of a year.
    :type INTSTEP: float

    :returns: structure of nodes

    .. note::
       Initial node acceleration is still a free parameter.
    """
    
    AGE = np.linspace(age,age + LengthofTime,LengthofTime/INTSTEP)
    
    # Parameters of emmetropic eye from 2004 paper:
    Child_Eye_Length = EyeLength(age) # in mm
    initial_node_acc = NodeRate(age,INTSTEP)

    Radius = Child_Eye_Length/2.0
    F_Px_initial = np.ones((NumberofNodes)) * Pressure * 0.5
    F_Py_initial = np.ones((NumberofNodes)) * Pressure * 0.5

    # memory for x positions: 
    x = np.zeros(( (LengthofTime / INTSTEP) + 1, NumberofNodes))
    y = np.zeros(( (LengthofTime / INTSTEP) + 1, NumberofNodes))
    
    # set the start location of each of the nodes.
    location = np.linspace(0, 2 * np.pi - (2 * np.pi / NumberofNodes), NumberofNodes)
    x[0,:] = np.cos(location) # 16mm typical size @ birth.
    y[0,:] = np.sin(location)

    Pressure = Pressure/NumberofNodes

    THETA, R = dm.cart2pol(x[0,:],y[0,:])
    I_initial, foo = dm.MatlabSort(THETA)
    I_initial = I_initial[np.newaxis,:]
    Initial_Area = np.sum((R**2) * np.pi / len(R))

    # Preallocate memory
    # Memory for checks:
    r  = np.zeros(((LengthofTime/INTSTEP)+1, NumberofNodes))
    theta = np.zeros(((LengthofTime/INTSTEP)+1, NumberofNodes))
    I = np.zeros(((LengthofTime/INTSTEP)+1, NumberofNodes))

    # Memory for velocities.
    xDOT,yDOT = dm.pol2cart(THETA,initial_node_acc * 4)
    xDOT = np.ones((NumberofNodes)) * xDOT
    yDOT = np.ones((NumberofNodes)) * yDOT

    # Memory for force.
    xx = np.zeros((NumberofNodes))
    yy = np.zeros((NumberofNodes))

    # Memory for eye length.
    ModelEyeLength = np.zeros((LengthofTime / INTSTEP,1))

    Nodes = {'x': x,
             'y': y,
             'radius': Radius,
             'f_px_initial': F_Px_initial,
             'f_py_initial': F_Py_initial,
             'f_px': [],
             'f_py': [],
             'node_accel': initial_node_acc,
             'pressure': Pressure,
             'xDOT': xDOT,
             'yDOT': yDOT,
             'xx': xx,
             'yy': yy,
             'initial_area': Initial_Area,
             'i_initial': I_initial[0],
             'r': r,
             'theta': theta,
             'i': I,
             'modeleyelength': ModelEyeLength,
             'age': AGE
             }

    return Nodes

    
### Compute Functions ###
#########################
    
def CircularDist(X,Y):
    """
    
    computes the distance between points on a circle
    Formula = sqrt((X(i) - X(i+1))^2 + (Y(i) - Y(i+1))^2)
    
    
    :params X: array of x values.
    :type X: np.array
    :params Y: array of y values.
    :type Y: np.array
    
    :returns: array of distances in form: (Point(1) - Point(2) ... (Point(N) - Point(1))

    .. note::
       Called by ComputeSpringForce
    """

    NumberofNodes = len(X) - 1
    distance = np.zeros((NumberofNodes + 1))

    for i in range(0,(NumberofNodes)):
        distance[i] = np.sqrt((X[i] - X[i+1])**2 + (Y[i] - Y[i+1])**2)

    # complete the circle, calcute distance between first and last node.
    distance[NumberofNodes] = np.sqrt((X[0] - X[NumberofNodes])**2 + (Y[0] - Y[NumberofNodes])**2)
    
    return distance
    
    
def ComputeSpringForce(Nodes,Spring_Eq,K,time):
    """
    ComputeSpringForce(x,y,Spring_Eq,K,time)
    
    Finds the inward force exerted by the two springs pulling on each node. 
    Uses Hooke's Law and assumes two springs pull on each node.  

    :param Nodes: structure containing nodes
    :type Nodes: struct
    :param Spring_Eq: array of spring equilibrium values.
    :type Spring_Eq: np.array
    :param K: array of K values.
    :type K: np.array
    
    :returns: updated node structure.

    .. note::    
       K values typically change with age.
    
    """
    
    theta,foo = dm.cart2pol(Nodes['x'][time], Nodes['y'][time])
    distance = CircularDist(Nodes['x'][time], Nodes['y'][time])
        
    F_s1 = -K * (Spring_Eq - distance)
    F_s2 = np.r_[F_s1[-1], F_s1[0:-1]]
    F_S = F_s1 + F_s2
    Nodes['f_sx'],Nodes['f_sy'] = dm.pol2cart(theta, F_S)
    
    return Nodes



def ComputePressureForce(Nodes,time):
    """

    Update the pressure force on the nodes.
    
    :param Nodes: node structure.
    :type Nodes: struct

    :returns Nodes: updated node structure.

    """
    F_Px = Nodes['pressure'] * (Nodes['x'][time] / (abs(Nodes['x'][time]) + abs(Nodes['y'][time])))
    F_Py = Nodes['pressure'] * (Nodes['y'][time] / (abs(Nodes['x'][time]) + abs(Nodes['y'][time])))

    hyp = np.hypot(F_Px, F_Py)
    Nodes['f_px'] = F_Px * (Nodes['pressure'] / hyp)
    Nodes['f_py'] = F_Py * (Nodes['pressure'] / hyp)
    
    return Nodes
            

### Check functions ###
#######################            
            
def CheckNodes(Nodes,time):
    """
    
    Make sure that the nodes do not cross over or move in towards the center.
    

    :param Nodes: node structure
    :param time: time in the integration algorithm
    
    :returns: Nodes updated node structure.

    """

    # Conditionals constraining the motion of nodes
    Nodes['theta'][time+1,:], Nodes['r'][time+1,:] = dm.cart2pol(Nodes['x'][time], Nodes['y'][time])
    
    # Ensure that no nodes move inwards
    a = np.nonzero(Nodes['r'][time+1,:] - Nodes['r'][time,:]<0)[0]    
    
    Nodes['i'][time+1,:],n = dm.MatlabSort(Nodes['theta'][time+1,:])
    
    if time > 1:
        # Ensure that no nodes cross over
        c = np.nonzero([Nodes['i'][time+1,:] - Nodes['i_initial']!=0])[0]
        
    else:
        a = np.array([])
        c = np.array([])
    
    if len(a) > 0 and len(c) > 0:
        
        Nodes['r'][time+1,a] = Nodes['r'][time,a]
        d = Nodes['i'][time+1,c]
        Nodes['theta'][time+1,d] = Nodes['theta'][time,d]
        Nodes['x'],Nodes['y'] = dm.pol2cart(Nodes['theta'][time+1,:], Nodes['r'][time+1,:])
        
    elif np.sum(c)>0:
        
        d = Nodes['i'][time+1,c]
        Nodes['theta'][time+1,d] = Nodes['theta'][time,d]
        Nodes['x'],Nodes['y'] = dm.pol2cart(Nodes['theta'][time+1,:], Nodes['r'][time+1,:])
        
    elif np.sum(a)>0:
        
        Nodes['r'][time+1,a]= Nodes['r'][time,a]
        Nodes['x'][time],Nodes['y'][time] = dm.pol2cart(Nodes['theta'][time+1,:], Nodes['r'][time+1,:])
        
    return Nodes

### Plotting Functions ###
##########################
    
def plot_eye(Nodes,axes = None):
    """
    
    Create a movie of eye growth. To be used with EyeGrowthFunc

    :param Nodes: structure containing nodes
    :type Nodes: struct
    :param INTSTEP: time step used for integration
    :type INTSTEP: int

    :returns: plot handle for Node plot.  Used to update during for loop.

    .. note::
       Called in EyeGrowthFunc
    """
    
    
    #set plotting parameters:
    if axes == None:
        fig = plt.figure(figsize=(10, 8))
        axes = fig.add_subplot(111,aspect='equal')
        plt.xlim([-13, 13])
        plt.ylim([-13, 13])
        
    axes.plot(np.r_[ Nodes['x'][0],Nodes['x'][0,0] ] * Nodes['radius'], 
            np.r_[ Nodes['y'][0], Nodes['y'][0,0] ] * Nodes['radius'], 
            '-ok', markerfacecolor = 'k',linewidth = 4, markersize = 10)
                
    axes = pf.TufteAxis(axes,['left','bottom'])
    #axes.set_axis_bgcolor('w')


    return axes
    
    
    
### Animation functions ###
###########################

class EyeAnimation(animation.TimedAnimation):
    """
    """
    
    def __init__(self):
        """
        """

        fig = plt.figure(facecolor = 'w', figsize = [12, 12])
        fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
        
        ax1 = fig.add_subplot(1, 1, 1, aspect='equal')
        
        self.t = np.arange(0, len(Nodes['age']))
        self.x = xx[self.t]
        self.y = yy[self.t]
        self.age = Nodes['age'][self.t]

        plot_eye(Nodes,ax1)

        ax1.set_xlabel('mm',fontsize=20)
        ax1.set_ylabel('mm',fontsize=20)

        self.line1 = Line2D([], [], color='red', linewidth=4)
        self.line1e = Line2D([], [], color='red', marker='o', markeredgecolor='r', markersize=10)
        
        self.text = ax1.text(0.05, 0.05, 'Age: %s'%self.age[0] , fontsize=18, animated = True,
                            transform=ax1.transAxes)
                            
        ax1.add_line(self.line1)
        ax1.add_line(self.line1e)

        
        ax1.set_xlim(-14, 14)
        ax1.set_ylim(-14, 14)
        ax1.set_title('eye growth simulation',fontsize=20)


        ax1.set_xticks([-10, -5, 0, 5, 10])
        ax1.set_yticks([-10, -5, 0, 5, 10])
        
        
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.tight_layout()
        
        animation.TimedAnimation.__init__(self, fig, interval=10, blit=True)

    def _draw_frame(self, framedata):
        """
        """

        i = framedata
        
        self.line1.set_data(self.x[i], self.y[i])
        self.line1e.set_data(self.x[i], self.y[i])

        self.text.set_text('Age: %d'%self.age[i])

        self._drawn_artists = [self.line1, self.line1e, self.text]

    def new_frame_seq(self):
        """
        """
        
        return iter(range(self.t.size))

    def _init_draw(self):
        """
        """
        
        lines =  [self.line1, self.line1e]
        for l in lines:
            l.set_data([], [])



    
def EyeLengthPlot(Nodes1,Nodes2=0,Nodes3=0,FONTSIZE = 20,LINEWIDTH=2.5):
    """
    Plot eye growth reported in Zadnik et al. against various models
    """

    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pf.AxisFormat()
    pf.TufteAxis(ax,['left','bottom'], [5,5])
    
    PaperEyeLength = np.zeros(len(Nodes1['age']))
    for i in range(0,len(Nodes1['age'])):
        PaperEyeLength[i] = EyeLength(Nodes1['age'][i])
    
    ax.plot(Nodes1['age'], Nodes1['modeleyelength'],'--b',
            linewidth = LINEWIDTH)
   
    plt.xlabel('age',fontsize=FONTSIZE)
    plt.ylabel('axial length (mm)', fontsize=FONTSIZE)
    plt.xlim([Nodes1['age'][0], Nodes1['age'][-1]])
    plt.xticks(fontsize=FONTSIZE)
    plt.yticks(fontsize=FONTSIZE)
    
    if isinstance(Nodes2, dict):
        plt.hold(True)
        ax.plot(Nodes1['age'], Nodes2['modeleyelength'],'-b',
                linewidth = LINEWIDTH)
        plt.hold(True)
        ax.plot(Nodes1['age'],PaperEyeLength,'-k',linewidth=LINEWIDTH)
        ax.legend(['Myopia','Emmetropia','Zadnik et al.'],loc=4,
                  prop={'size':FONTSIZE}).get_frame().set_edgecolor('white')

    if isinstance(Nodes3, dict):
        plt.hold(True)
        ax.plot(Nodes1['age'], Nodes3['modeleyelength'],'-b',
                linewidth = LINEWIDTH)
        plt.hold(True)
        ax.plot(Nodes1['age'],PaperEyeLength,'-k',linewidth=LINEWIDTH)
        ax.legend(['Myopia','Emmetropia','Zadnik et al.'],loc=4,
                  prop={'size':FONTSIZE}).get_frame().set_edgecolor('white')
    
    if Nodes2 == 0:
        plt.hold(True)
        ax.plot(Nodes1['age'],PaperEyeLength,'-k',linewidth=LINEWIDTH)
        ax.legend(['Model','Zadnik et al.'],loc=4,
                  prop={'size':FONTSIZE}).get_frame().set_edgecolor('white')
    
    
    plt.tight_layout()
    plt.show(block=False)

    


### Zadink et al. Functions ###
###############################
def EyeLength(age):
    """
    this function derived from Zadnik et al. (2004) - Emmetropic Eye
    Growth

    :param age: age in years at which to compute axial length.
    :type age: float
    
    :returns: eyeLength axial length of eye in mm.

    .. note::
       This function is called by GenerateNodes.m and Eye_Growth_Main.m for
       plotting.
    """
    if age < 10.5:
        eyeLength = 20.189 + 1.258 * np.log(age)
    elif age >= 10.5:
        eyeLength = 21.353 + 0.759 * np.log(age)
    return eyeLength
    
    

def NodeRate(age,timeStep):
    """
    this function uses the best fit model from Zadnik et al. (2004) for axial
    growth to find growth rate over a given timeStep.


    :param age: age in years 
    :param timeStep: duration of time over which to compute the rate (age-timeStep), in years.

    :returns: rate of eye growth at a give age, during a given duration of time.

    .. note::
       Called by GenerateNodes.
    """
    
    if age < 10.5:
        length1 = 20.189 + 1.258 * np.log(age)
        length2 = 20.189 + 1.258 * np.log(age - timeStep)
    elif age >= 10.5:
        length1 = 21.353 + 0.759 * np.log(age)
        length2 = 21.353 + 0.759 * np.log(age - timeStep)
    
    rate = length1 - length2

    return rate

def main():
    """
    Run basic program
    """
    NodesEmmetr = EyeGrowthFunc(Animate=0)     
    NodesMyopia = EyeGrowthFunc(Animate=0,Loosen=1)
 
    EyeLengthPlot(NodesMyopia,NodesEmmetr)
    plt.show()    

    
    
# if enter 'python Eye_Growth.py' into console, this will run.
if __name__ == "__main__":
    input_params = getParams()
    
    NodesEm,xx,yy = EyeGrowthFunc(Animate=0, **input_params)    
    #ani = EyeAnimation()
    #ani.save('eye.mp4', fps=20, codec='mpeg4', clear_temp=True, frame_prefix='_tmp')
    
    NodesMy,xx,yy = EyeGrowthFunc(Animate=0,Loosen=1, **input_params)
    #ani2 = EyeAnimation()
    
    EyeLengthPlot(NodesEm,NodesMy)
    plt.show(block=True)
    
    