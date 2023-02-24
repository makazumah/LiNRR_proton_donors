from ase.io import write,read
from ase.io.bader import attach_charges
from ase.units import Bohr
from gpaw import GPAW, restart
from gpaw.analyse.hirshfeld import HirshfeldPartitioning
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from ase.io.cube import read_cube_data
import time

from_scratch = False
c_path = os.path.abspath(os.getcwd())
in_traj= c_path+ '/input.traj'
if os.path.isfile("relax.gpw"):
    
    molecule, _ = restart("relax.gpw", txt="relax.txt")
    molecule.calc = GPAW(
                    h=0.16,
                    txt="bader.txt"
                    )
    molecule.get_potential_energy()
    hf = HirshfeldPartitioning(molecule.calc)
    for atom, charge in zip(molecule, hf.get_charges()):
        atom.charge = charge
    molecule.copy().write('Hirshfeld.traj')
    rho = molecule.calc.get_all_electron_density(gridrefinement=4)
    write('density.cube', molecule, data=rho * Bohr**3)
    
    while not os.path.exists("Hirshfeld.traj"):
        time.sleep(1)
        
    if os.path.isfile("Hirshfeld.traj") and os.path.isfile("density.cube"):
        hirsh = read('Hirshfeld.traj')
        bader_out = hirsh.copy()
                
        os.system("bader -p all_atom -p atom_index -p all_bader -p bader_index density.cube")
        while not os.path.exists("ACF.dat"):
            time.sleep(1)
        
        while not os.path.exists("BCF.dat"):
            time.sleep(1)
        while not os.path.exists("BvIndex.cube"):
            time.sleep(1)
            
        while not os.path.exists("AtIndex.cube"):
            time.sleep(1)
            
        while not os.path.exists("AVF.dat"):
            time.sleep(1) 
            
        
    if os.path.isfile("AtIndex.cube") and os.path.isfile("ACF.dat"):
        attach_charges(bader_out)
        sourceFile = open(c_path + '/Hirsch.txt', 'w')
        

        print('atom Hirshfeld Bader', file =sourceFile)
        for ah, ab in zip(hirsh, bader_out):
            assert ah.symbol == ab.symbol
            print(f'{ah.symbol:4s} {ah.charge:9.2f} {ab.charge:5.2f}', file =sourceFile)

        if os.path.isfile('molecule.pckl'):
            with open('molecule.pckl', 'rb') as fd:
                dens, bader1, atoms = pickle.load(fd)
        else:
            dens, atoms = read_cube_data('density.cube')
            bader1, atoms = read_cube_data('AtIndex.cube')
            x = len(dens) // 2
            dens = dens[x]
            bader1 = bader1[x]
            with open('molecule.pckl', 'wb') as fd:
                pickle.dump((dens, bader1, atoms), fd)

        x0, y0, z0 = atoms.positions[0]
        y = np.linspace(0, atoms.cell[1, 1], len(dens), endpoint=False) - y0
        z = np.linspace(0, atoms.cell[2, 2], len(dens[0]), endpoint=False) - z0
        print(y.shape, z.shape, dens.shape, bader1.shape)
        print(atoms.positions)
        print(dens.min(), dens.mean(), dens.max())
        plt.figure(figsize=(5, 5))
        plt.contourf(z, y, dens, np.linspace(0.01, 0.9, 15))
        # plt.contour(z, y, bader1, [1.5], colors='k')
        # plt.axis(xmin=-2, xmax=2, ymin=-2, ymax=2)
        plt.savefig('molecule-bader.png')
else:
    print("Not ready for Bader Analysis")