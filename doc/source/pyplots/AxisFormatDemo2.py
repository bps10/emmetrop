import matplotlib.pylab as plt
import numpy as np
from renderer import PlottingFun as pf

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

pf.AxisFormat()
   
ax.plot(np.arange(0,100), np.arange(0,100), 'b-', linewidth =2)  

plt.show()    