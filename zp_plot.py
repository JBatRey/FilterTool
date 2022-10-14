import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import patches
from matplotlib.figure import Figure
from matplotlib import rcParams
    
def zplane(z,p,filename=None):
    """Plot the complex z-plane given a transfer function.
    """

    # get a figure/plot
    ax = plt.subplot(111)

    # create the unit circle
    #uc = patches.Circle((0,0), radius=1, fill=False,
                        #color='black', ls='dashed')
    #ax.add_patch(uc)
 
    
    # Plot the zeros and set marker properties    
    t1 = plt.plot(z.real, z.imag, 'bo', ms=5)
    plt.setp( t1, markersize=10.0, markeredgewidth=1.0,
              markeredgecolor='k', markerfacecolor='g')

    # Plot the poles and set marker properties
    t2 = plt.plot(p.real, p.imag, 'rx', ms=5)
    plt.setp( t2, markersize=12.0, markeredgewidth=3.0,
              markeredgecolor='r', markerfacecolor='r')

    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_position('zero')
    ax.spines['top'].set_visible(False)
    ax.set_aspect('equal')
    ax.set_adjustable('datalim')

    # set the ticks
    #r = 1.5; plt.axis('scaled'); plt.axis([-r, r, -r, r])
    #ticks = [-1, -.5, .5, 1]; plt.xticks(ticks); plt.yticks(ticks)

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)
    
