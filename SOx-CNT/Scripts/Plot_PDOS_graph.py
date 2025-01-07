from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.plotter import DosPlotter
import matplotlib.pyplot as plt
import numpy as np

#index = site, orbital,'densities', spin
def pdos_densities(pdos,site,orbital,spin):
    if spin ==0:
        s='1'
    else:
        s='-1'
    return np.array(pdos[site][orbital]['densities'][s])
def tdos_densities(cdos, spin):
    if spin ==0:
        s='1'
    else:
        s='-1'
    return np.array(cdos.as_dict()['densities'][s])

v = Vasprun('vasprun.xml')
cdos = v.complete_dos
pdos=cdos.as_dict()['pdos']

print(cdos)

#Mn vacancies look like this
# 28 X 25 29 X 26
#    24     27

# up_mag = [24,25,26,27]
# dw_mag = [28,29]

C = np.arange(0,363)
S = [364]
#O = [365,366, 367]

energies = np.array(cdos.as_dict()['energies']) - cdos.as_dict()['efermi']
print(energies)

TotDOS = tdos_densities(cdos, 0) 

print(TotDOS)

C_s_u = [sum([pdos_densities(pdos, site,orbital,0) for site in C]) for orbital in ['s']]
#C_s_d = [sum([pdos_densities(pdos, site,orbital,1) for site in C]) for orbital in ['s']]
C_p_u1 = [sum([pdos_densities(pdos, site,orbital,0) for site in C]) for orbital in ['px', 'py', 'pz']]
#C_p_d1 = [sum([pdos_densities(pdos, site,orbital,1) for site in C]) for orbital in ['px', 'py', 'pz']]

print(C_p_u1)

# Mn2_s_up = [sum([pdos_densities(pdos, site,orbital,0) for site in dw_mag]) for orbital in ['s']]
# Mn2_s_dw = [sum([pdos_densities(pdos, site,orbital,1) for site in dw_mag]) for orbital in ['s']]
S_s_u = [sum([pdos_densities(pdos, site,orbital,0) for site in S]) for orbital in ['s']]
S_p_u1 = [sum([pdos_densities(pdos, site,orbital,0) for site in S]) for orbital in ['px', 'py', 'pz']]
#S_s_d = [sum([pdos_densities(pdos, site,orbital,1) for site in S]) for orbital in ['s']]
#S_p_d1 = [sum([pdos_densities(pdos, site,orbital,1) for site in S]) for orbital in ['px', 'py', 'pz']]

#O_s_u = [sum([pdos_densities(pdos, site,orbital,0) for site in O]) for orbital in ['s']]
#O_p_u1 = [sum([pdos_densities(pdos, site,orbital,0) for site in O]) for orbital in ['px', 'py', 'pz']]
#O_s_d = [sum([pdos_densities(pdos, site,orbital,1) for site in O]) for orbital in ['s']]
#O_p_d1 = [sum([pdos_densities(pdos, site,orbital,1) for site in O]) for orbital in ['px', 'py', 'pz']]

#pdos
x = energies #[::2]
#C_s = np.add(C_s_u, C_s_d)
#C_p = np.add(C_p_u, C_p_d)
#print(C_p)

#S_s = np.add(S_s_u, S_s_d)
#S_p = np.add(S_p_u, S_p_d)

#O_s = np.add(O_s_u, O_s_d) 
#O_p1 = np.add(O_p_u, O_p_d)

#O_p2 = np.add(O_p1[0], O_p1[1])
#O_p = np.add( O_p1[2], O_p2 )
#print(O_p) 

C_p3 = np.add(C_p_u1[0], C_p_u1[1])
C_p_u = np.add( C_p_u1[2], C_p3 )

S_p3 = np.add(S_p_u1[0], S_p_u1[1])
S_p_u = np.add( S_p_u1[2], S_p3 )

#O_p3 = np.add(O_p_u1[0], O_p_u1[1])
#O_p_u = np.add( O_p_u1[2], O_p3 )

# y1 = sum(C_s) #[::2]
# y2 = sum(C_p) #[::2]
# y3 = sum(Te_p_up) #[::2]
# y4 = sum(F_p_up) #[::2]


#save to csv
#header = 'energy (eV), Mn1_d_up, Mn1_d_dw, Mn2_d_up, Mn2_d_dw, Te_p_up, Te_p_dw, F_p_up, F_p_dw'
# np.savetxt('pdos.csv',
#            np.vstack([x,y1,y2,y3,y4,y5,y6,y7,y8]).T,
#            delimiter=', ',
#            header = header)

fig,ax = plt.subplots()
# fig.set_size_inches(14, 3)

ax.plot(x, TotDOS, linestyle = "-", color = "#000000")

ax.plot(x, C_s_u[0], linestyle = "-", color = "#00ff80")
ax.plot(x, C_p_u, linestyle = "-", color = "#0000FF")
ax.plot(x, S_s_u[0], linestyle = "-", color = "#ff8000")
ax.plot(x, S_p_u, linestyle = "-", color = "#8000ff")
#ax.plot(x, O_s_u[0], linestyle = "-", color = "#ff46a3")
#ax.plot(x, O_p_u, linestyle = "-", color = "#e50000")

# plt.tick_params(
#     axis='y',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     left=False,      # ticks along the bottom edge are off
#     right=False,         # ticks along the top edge are off
#     labelleft=False) # labels along the bottom edge are off
# ax.tick_params(axis='both', which='major', labelsize=13)
# #ax.tick_params(axis='both', which='minor', labelsize=8)

plt.ylabel('DOS, (state/eV)', fontsize=16)
plt.xlabel('Energy, (eV)', fontsize=16)

#ax.set_xlim([-8,5])
ax.set_xlim([-6,3])
ax.set_ylim([0,10])
ax.legend(['Total DOS','C(s)','C(p)','S(s)','S(p)','O(s)','O(p)'], fontsize=11, loc='upper left')

fig.savefig('CNT-S_ether-d-ringDistortion_dos-zoom2.png', dpi=500)
#fig.savefig('pdos-zoom.png', dpi=300)
