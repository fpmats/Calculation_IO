import numpy as np
import math as math
import matplotlib.pyplot as plt
from ase.io import read
import sys
from monochalcogenpy.plotter import plot_EV, plot_theta
from monochalcogenpy.tilt_angle import tilt_angle

def dashboard(traj, time_step, image_name = './dashboard.png'):
    time = time_step * np.arange(len(traj)) #ps
    fig, axs = plt.subplots(
        2,1,
        figsize=[4,3], 
        height_ratios= [1,1],
        sharex=True,
        #layout="constrained",
        gridspec_kw={'hspace': 0})

    
    E_atom = np.array([atoms.get_total_energy() for atoms in traj])/len(traj[0])
    V = np.array([atoms.get_volume() for atoms in traj])
    plot_EV(time, E_atom, V, axs[0])
    projections = ['bc', 'ac']
    for p,c in zip(projections, ['tab:orange','tab:green']):
        theta = [np.concatenate(list(tilt_angle(atoms, p, thres=0.95, absolute=False).values())) for atoms in traj]
        plot_theta(time,theta,c,p[0],axs[1])
    axs[1].legend(loc='upper right', fontsize=8)
    fig.supxlabel('time (ps)', y=0.025, fontsize=10)
    axs[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    axs[0].set_box_aspect(0.6)
    axs[1].set_box_aspect(0.6)
    #fig.tight_layout()
    #plt.subplots_adjust(left=0.25, bottom=0.2, right=0.8, top=0.8, wspace=None, hspace=None)
    #fig.title(image_name[7:-3])
    fig.savefig(image_name, dpi=300)
    return fig, axs

if __name__=='__main__':
    traj = read(sys.argv[1], ':')
    time_step = float(sys.argv[2])
    fig, axs = dashboard(traj, time_step, image_name = './dashboard.png')
