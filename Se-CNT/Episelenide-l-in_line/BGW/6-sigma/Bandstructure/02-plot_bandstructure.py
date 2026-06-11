#!/usr/bin/env python

# FHJ: This script generates Si.bands.pdf from the output of inteqp
# using numpy & matplotlib. K-point path is hard-coded into the script.

from __future__ import print_function

import sys
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Number of valence states. This is hard coded to set the Fermi energy at zero.
nv = 12
# Minimum/maximum energy to plot relative to the VBM.
emin, emax = -1.5, 3.0

# Load bandstructure file. Format is:
# spin, band, kx, ky, kz, emf, eqp, eqp-emf
# where kx, ky, kz => Cartesian coordinates
# emf, eqp => eV
data = np.loadtxt('bandstructure.dat')

# Bands in this file
bands = np.unique(data[:,1])
# Number of bands
nb = len(bands)
DB = nv   ## index for defect bant
VBM = nv-1   ## index for valence band max
CBM = nv+1 ## index for conduction band min

# We want the "x" axis of our plot to be |dk| in Cartesian coordinates.
# K-points in Cartesian coordinates. We don't select repeated k-points, so only
# select the k-points that appear for the first band.
cond_first_band = data[:,1]==bands[0]
kpoints = data[cond_first_band][:,2:5]
# Number of k-points
nk = len(kpoints)
# Vector dk, difference between neighboring k-points
dk = kpoints[1:] - kpoints[:-1]
# Scalar |dk|
dk = np.linalg.norm(dk, axis=1)
# Note that len(dk) is nk - 1, where nk is the number of k-points. So, we need
# to add an extra dk corresponding to the origin of the plot.
dk = np.insert(dk, 0, 0)
# Now, the x axis is simply the cumulative sum of dk
x = np.cumsum(dk)

# Reshape LDA and GW QP energies. Energies always appear sequentially in the
# file as [ib,nk], where ik is the fast index, corresponding to the k-point,
# and ib is the slow index, corresponding to the band index.
# LDA mean-field energies
elda = data[:,5].reshape((nb,nk))
# GW QP energies
eqp = data[:,6].reshape((nb,nk))

# Set LDA and GW Fermi energy to 0
nv = 12
elda -= np.amax(elda[nv-1])
eqp -= np.amax(eqp[nv-1])

# Report DFT and GW gaps:
gap_lda = np.amin(elda[nv]) - np.amax(elda[nv-1])
gap_gw = np.amin(eqp[nv]) - np.amax(eqp[nv-1])
print('PBE gap: {:.3f} eV'.format(gap_lda))
print('GW gap: {:.3f} eV'.format(gap_gw))

# Cosmetic updates
rc = matplotlib.rc
rc('figure', figsize=(6.0, 8.0))
rc('axes', linewidth=2)
rc('lines', linewidth=2)
rc('font', size=18.0)
rc('xtick.major', size=0.0, pad=8.0)
rc('xtick.minor', size=0.0, pad=8.0)
rc('ytick.major', size=6.0, pad=8.0)
rc('ytick.minor', size=3.0, pad=8.0)

# Plot MF and QP energies
for ib in range(nb):
    plt_qp, = plt.plot(x, eqp[ib], 'black')
    #plt_lda, = plt.plot(x, elda[ib], 'g-')
    if ib==0:
        plt_qp.set_label('GW')
        #plt_lda.set_label('LDA')
    if ib==VBM: 
        plt_qp, = plt.plot(x, eqp[ib], 'red')
    if ib==DB: 
        plt_qp, = plt.plot(x, eqp[ib], '#9600FF')
    if ib==CBM: 
        plt_qp, = plt.plot(x, eqp[ib], 'blue')

# K-point path, manually specified
k_special_index = np.array([0, 30])
k_special_label = np.array(['$\Gamma$', 'Z'])
k_special_x = x[k_special_index]
plt.xticks(k_special_x, k_special_label)

for ksx in k_special_x[1:-1]:
    plt.axvline(ksx, color='k')

# Put Fermi energy
plt.axhline(0,linestyle=':', color='#999999', zorder=0)

# Cosmetic stuff
ax = plt.gca()
for i in ax.get_xticklines() + ax.get_yticklines():
    i.set_markeredgewidth(1.5)

#plt.xlabel('Wavevector')
plt.ylabel('Energy (eV)')
plt.xlim(x[0], x[-1])
plt.ylim(emin, emax)
#plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#           ncol=1, mode="expand", borderaxespad=0.)

plt.savefig('Se-ether-l-NewPos-CNT_bands.pdf', bbox_inches='tight')
