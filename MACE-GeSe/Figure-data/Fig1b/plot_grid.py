from ase.db import connect
import ase.db
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.cm as cm
from itertools import product
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import griddata
import matplotlib.tri as tri
from ase.io.vasp import write_vasp


##-----------------------------------------------------------
## read in the data from the ASE database
##-----------------------------------------------------------

db = connect('./relax-vasp-and-mace.db')

##-----------------------------------------------------------
## analize and graph data as a multi-subplot contour plot
##-----------------------------------------------------------

natoms = db[1].natoms
print("number of atoms:", natoms)
nlayer = natoms/4
print(nlayer)

## Setting up the subplot and what data to pull
##---------------------------------------------------
## MACE-MP = graphs a 1x3 graph with MACE[PBE]+D3, MACE-MP+D3 and PBE+D3 as seen in the SI 
## No-MACE-MP = graphs a 1x2 graph with MACE[PBE]+D3 and PBE+D3 as seen in Fig 1b of the main text
##---------------------------------------------------


# Graph = "MACE-MP"
Graph = "No-MACE-MP"

if Graph == "No-MACE-MP":
    fig, axs = plt.subplots(1,2,figsize=[5,2.5], sharex=True,sharey=True, width_ratios=[0.42,0.55])
    plt.subplots_adjust(wspace=None, hspace=None)
    ## data to pull from database: vasp = PBE+D3, macecalculator = MACE[PBE]+D3
    selections = ['model=vasp', 'model=macecalculator']

if Graph == "MACE-MP":
    fig, axs = plt.subplots(1,3,figsize=[8,2.5], sharex=True,sharey=True, width_ratios=[0.2,0.2,0.25])
    selections = ['model=macecalculator','model=vasp', 'model=mace-mpcalculator']


vmin = np.min([row.energy for row in db.select()])
print("Min energy from total database for all three maodels:", vmin)

## color map to use on the plot 
cmap='RdYlBu'

for selection,ax in zip(selections,axs):
    ax.set_box_aspect(1)

    b_array = np.array([row.b for row in db.select(selection)])
    a_array = np.array([row.a for row in db.select(selection)])
    energy = np.array([row.energy for row in db.select(selection)])
    eng = [row.energy for row in db.select(selection)]
    print(selection)
    b_sort = np.sort(np.unique(b_array))
    a_sort = np.sort(np.unique(a_array))
    ## prints the unique a,b values for the 2D grid and double checks they are the same length 
    print("b values:")
    print(b_sort)
    print("a values:")
    print(a_sort)
    print("length of b list", len(b_sort), "length of a list", len(a_sort))

    energy_all = np.concatenate((energy, np.flip(energy)))
    B,A = np.meshgrid(b_sort, a_sort)
    E = np.zeros_like(B)

    ## finds the min and max energy for the chosen model (i.e PBE, MACE, MACE-MP)
    vmin = np.min([row.energy for row in db.select(selection)])
    vmax = np.max([row.energy for row in db.select(selection)])

    print("vmin of "+str(selection), vmin)
    print("vmax of "+str(selection), vmax)

    vmax = vmin+0.48
    vmax2 = vmin+0.6

    levels = (np.linspace(vmin, vmax, 11)-vmin)*1000/natoms
    levels2 = (np.linspace(vmin, vmax2, 21)-vmin)*1000/nlayer

    if selection == 'model=mace-mpcalculator':
        print("plotting MACE-MP")

        B,A = np.meshgrid(b_sort, a_sort)
        E = np.zeros_like(B)
        for (i,j) in np.array(np.where(B)).T:
            E[i,j] = db.get(f'b={B[i,j]},a={A[i,j]},'+selection).energy

            # per layer
        cs = ax.contourf(B,A,1000*(E-vmin)/nlayer,
                levels2,
                extend='max',
                cmap=cmap,
                )
    
    if selection == 'model=macecalculator' or selection == 'model=vasp':
        
        print(selection)
        energy = np.array([row.energy for row in db.select(selection)])

        ## Finds min energy for selected model and the python index for said energy
        Emin = min(energy)
        Min_ID = np.where(energy == Emin)

        print("length of b,a,energy arrays:", len(b_array), len(a_array), len(energy))
        print("length of sorted b,a arrays:", len(b_sort), len(a_sort))
        id_array = np.array([row.id for row in db.select(selection)])

        if selection == 'model=vasp':
            print("this is the min energy ID:", id_array[Min_ID[0]][0])
            ## loop set up to pull information like min energy point and diagonal energies 
            print("this energy for VASP:")
            print("b, a, energy (eV), database ID")
            for row in db.select(selection): 
                b_t = row.b
                a_t = row.a
                c_t = row.c
                ID_ind = row.id
                ## uncomment to get the min energy point of double well
                # if ID_ind == id_array[Min_ID[0]][0]: 
                #     print(row.b, row.a, row.energy, row.id)
                if b_t == a_t: 
                    print(row.b, row.a, row.energy, row.id)

        if selection == 'model=macecalculator':
            print("this is the min energy ID:", id_array[Min_ID[0]][0])
            ## loop set up to pull information like min energy point and diagonal energies 
            print("this energy for VASP:")
            print("b, a, energy (eV), database ID")
            for row in db.select(selection): 
                b_t = row.b
                a_t = row.a
                c_t = row.c
                ID_ind = row.id
                ## uncomment to get the min energy point of double well
                # if ID_ind == id_array[Min_ID[0]][0]: 
                #     print(row.b, row.a, row.energy, row.id)
                if b_t == a_t: 
                    print(row.b, row.a, row.energy, row.id)

        ##----------------------------------------------------------
        ## NOTE:
        ## this section corrects for rounding errors when a,b were entered in database (e.g 6.9999999 vs 7.0)
        ## It prevents double counting and allows contourf command to be used to generate a smoother 
        ## contour plot than tricontourf command
        ##----------------------------------------------------------
        B_array = []
        A_array = []
        Eng = []
        index = []
        for i in range(0,len(b_array)): 
            if b_array[i] == a_array[i]: 
                print("This is diagonal:", b_array[i], a_array[i])
                index += [i]
        print("index",index)

        for i in range(0,len(b_array)): 
            if b_array[i] != a_array[i]:
                B_array.append(b_array[i])
                A_array.append(a_array[i])
                Eng.append(energy[i])
            if b_array[i] == a_array[i]: 
                for test in index: 
                    if i == test: 
                        B_array.append(b_array[i])
                        A_array.append(a_array[i])
                        Eng.append(energy[i])
        print("length new:", len(B_array), len(A_array), len(Eng))
        for i in range(0,len(B_array)): 
            if B_array[i] == A_array[i]: 
                print("This is diagonal:", B_array[i], A_array[i])
        B_array = np.array(B_array[:])
        B_array = np.around(B_array, decimals=2)
        A_array = np.array(A_array[:])
        A_array = np.around(A_array, decimals=2)
        Eng = np.array(Eng[:])
        print("length rounded:", len(B_array), len(A_array), len(Eng))
        B_sort = np.sort(np.unique(B_array))
        A_sort = np.sort(np.unique(A_array))
        print("length sort:", len(B_sort), len(A_sort))
      
        ##---------------------------------------------

        if selection == 'model=macecalculator': 
            print("Graping MACE[PBE]")
            B,A = np.meshgrid(b_sort, a_sort)
            E = np.zeros_like(B)
            for (i,j) in np.array(np.where(B)).T:
                E[i,j] = db.get(f'b={B[i,j]},a={A[i,j]},'+selection).energy

            # per layer
            cs = ax.contourf(B,A,1000*(E-vmin)/nlayer,
                     levels2,
                     extend='max',
                     cmap=cmap,
                    )

        if selection == 'model=vasp':
            print("Graping PBE")
            B,A = np.meshgrid(b_sort, a_sort)
            E = np.zeros_like(B)
            for (i,j) in np.array(np.where(B)).T:
                E[i,j] = db.get(f'b={B[i,j]},a={A[i,j]},'+selection).energy

            # per layer
            cs = ax.contourf(B,A,1000*(E-vmin)/nlayer,
                     levels2,
                     extend='max',
                     cmap=cmap,
                     
                    )


## Setting up the axies, color bar, legend, and plot name
##---------------------------------------------------------------

if Graph == "No-MACE-MP":
    axs[0].xaxis.set_major_locator(MultipleLocator(0.2))
    axs[0].yaxis.set_major_locator(MultipleLocator(0.2))
    padding = 0.05
    axs[0].set_xlim([a_array.min()-padding,a_array.max()+padding])
    axs[0].set_ylim([b_array.min()-padding,b_array.max()+padding])
    axs[0].set_xlabel(r'b ($\AA{}$)') 
    axs[0].set_ylabel(r'a ($\AA{}$)') 
    axs[1].set_xlabel(r'b ($\AA{}$)')
    
    ## Labels for each plot 
    axs[1].text(3.8,4.9,'MACE[PBE]+D3') 
    axs[0].text(4.0,4.9,'PBE+D3') 

    cbar = fig.colorbar(cs,label="Energy (meV/layer)", pad=0.09, shrink=0.8) ## color bar

    axs[0].tick_params(axis='both', which='major', labelsize=8) 
    axs[1].tick_params(axis='both', which='major', labelsize=8) 

    fig.tight_layout()
    fig.savefig('GeSe_PES-layers.png',dpi=300)


if Graph == "MACE-MP":
    axs[0].xaxis.set_major_locator(MultipleLocator(0.2))
    axs[0].yaxis.set_major_locator(MultipleLocator(0.2))
    padding = 0.05
    axs[0].set_xlim([a_array.min()-padding,a_array.max()+padding])
    axs[0].set_ylim([b_array.min()-padding,b_array.max()+padding])
    axs[0].set_xlabel(r'b ($\AA{}$)') 
    axs[0].set_ylabel(r'a ($\AA{}$)') 
    axs[1].set_xlabel(r'b ($\AA{}$)')
    axs[2].set_xlabel(r'b ($\AA{}$)')

    ## Labels for each plot 
    axs[0].text(3.8,4.87,'MACE+D3, this work') 
    axs[1].text(4.1,4.87,'PBE+D3') 
    axs[2].text(3.95,4.87,'MACE-MP+D3') 

    fig.colorbar(cs,label="Energy (meV/layer)") ## color bar
    

    fig.tight_layout()

    fig.savefig('GeSe_PES-layers-MACE-MP.png',dpi=300)