import numpy as np
import matplotlib.ticker as ticker
from numpy.typing import NDArray

def plot_EV(time:NDArray, E_atom:NDArray, V:NDArray, ax1):
    #ax1.set_aspect(1)

    color = 'tab:red'
    ax1.set_ylabel('total energy \n(eV/atom)', color=color)
    ax1.plot(time, E_atom, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel(r'volume (nm$^3$)', color=color)  # we already handled the x-label with ax1
    ax2.plot(time, V/1000, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    ax1.xaxis.set_major_locator(ticker.MultipleLocator(100))
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(10))

    ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.01))
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.002))
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.025))

def plot_theta(time:NDArray, theta:NDArray, color: str, label,  ax, num_to_plot = 10):
    '''
    theta is a list of tilt angle present in each timestep
    '''
    #time = np.array(time)
    #theta = np.array(theta)
    interval = int(len(time)/num_to_plot)
    parts = ax.violinplot(
        theta[::interval], 
        positions=time[::interval], 
        widths=1.5*(time[interval]- time[0]), 
        bw_method = 0.025,
        side='low',showmeans=True)
    for pc in parts['bodies']:
        pc.set_facecolor(color)
    for lc in [parts['cbars'],parts['cmeans'], parts['cmins'], parts['cmaxes']]:
        lc.set_linewidth(0.5)
        lc.set_color(color)

    #X, Y = np.repeat(time, theta.shape[1]), theta.ravel()
    result = np.array([[time[::interval][i], x] for i in range(len(time[::interval])) for x in theta[::interval][i]])
    ax.plot(result[:,0],result[:,1], 'x', color=color, label=label, markersize=2, alpha=0.5)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(2))
    #ax.set_xlabel('time (ps)')
    ax.set_ylabel(r'tilt angle ($^{\circ}$)')
    # ax.set_ylabel(r'g(r) norm.')
    # ax.set_xlabel(r'r ($\AA{}$)')
    # ax.set_xlim([0,r_grid[-1] - 0.5])
