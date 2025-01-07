import numpy as np
from ase.io import read, write
from ase.io.jsonio import read_json
import matplotlib.pyplot as plt

bs = read_json('bs.json')
print(bs)
print(bs.energies)
print(len(bs.energies[0]))

Displ = 3.61

E_sk = bs.energies[0] 
E_sk_2 = E_sk + Displ

fermi = -3.5102 + Displ

kp = np.arange(0,10)
ratio = 2
plt.rcParams["figure.figsize"] = 3.37*2, (3.37)*ratio*2

plt.figure(1)
plt.axhline(y = fermi , color = 'r', linestyle = '-.')
plt.axhline(y = 0.0 , color = 'black', linestyle = '--')
plt.plot(kp, E_sk_2, linestyle = "-", marker = "o", color = "blue")
plt.xlim(0, 9)
plt.ylim(-1.5, 2.5)
plt.ylabel("Energies [ev]")
plt.xticks([0,9], ["$\Gamma$", "Z"])
plt.savefig("CNT-PBE-S-ether-d-ringdistortion.png", bbox_inches='tight')
#plt.show()
