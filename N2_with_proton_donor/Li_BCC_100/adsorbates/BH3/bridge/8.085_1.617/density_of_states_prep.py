from gpaw import GPAW, restart
import matplotlib.pyplot as plt
from ase.io import read, write
import os
import numpy as np
from collections import Counter
import pickle
from ase.units import Hartree
from gpaw.utilities.dos import fold

# Initialize slab
from_scratch = False
c_path = os.path.abspath(os.getcwd())
molecule_paths = ["/home/vazumah/Python_Projects/Li_NRR/Li_adsorb_studies/pure_structures_results/neutral_molecules/N2"] # Enter the path for the molecule to calc DOS for


def find_mol_inds(molecule_traj, slab):
    mol_sym = np.array(molecule_traj.get_chemical_symbols())
    slab_sym = np.array(slab.get_chemical_symbols()).astype(str)
    mol_inds = {}
    mol_state_dict= {}

    for mol in set(mol_sym):
        count=1
        result = np.where(slab_sym==mol)[0]
        for pos in result:
            pos_inds = np.arange(pos, pos+len(molecule_traj)).astype(int)
            if pos_inds[-1] < len(slab):
                # print(pos_inds)
                test_state = slab_sym[pos_inds]
                if (test_state==mol_sym).all():
                    mol_key = "".join([tup[0]+str(tup[1]) for tup in sorted(Counter(mol_sym.tolist()).items())])
                    mol_inds[mol_key+'-' +str(count)]= pos_inds
                    # print(test_state)
                    mol_state = range(len(slab))[pos_inds[0]:pos_inds[-1]+1]
                    mol_state_dict[mol_key+'-' +str(count)]= mol_state
                    count+=1
                    
                    
    return mol_inds, mol_state_dict

if os.path.isfile("output.gpw"):
    if os.path.isfile("output.traj") and os.path.getsize("output.traj") > 0:
        slab = read("output.traj")
    else:
        slab = read("history.traj")
    if not os.path.isfile('top.gpw'):
        slab.calc = GPAW(h=0.16,
                        kpts={'size': [4, 4, 1]},
                        occupations={'name': 'fermi-dirac', 'width': 0.05},
                        poissonsolver={'dipolelayer': 'xy'},
                        xc='BEEF-vdW',
                        convergence={'bands': -10},
                        txt='top.txt')
        slab.get_potential_energy()
        slab.calc.write('top.gpw', mode='all')

    for molecule_path in molecule_paths:
        
        if molecule_path[-1]!='/':
            o_molecule_path = molecule_path
            molecule_path = molecule_path + '/'
        if os.path.exists(molecule_path + "relax.traj"):
            molecule_traj = read(molecule_path + "relax.traj")
        elif os.path.exists(molecule_path + "history.traj"):
            molecule_traj = read(molecule_path + "history.traj")
        elif os.path.exists(molecule_path + "input.traj"):
            molecule_traj = read(molecule_path + "input.traj")
        elif os.path.exists(molecule_path + "input.sdf"):
            molecule_traj = read(molecule_path + "input.sdf")
        else:
            molecule_traj = read(o_molecule_path)
        
        mol_inds, mol_state = find_mol_inds(molecule_traj, slab)
        
        for di_key, di_val in enumerate(mol_state):
            molecule = mol_state[di_val]
            mol_ind = mol_inds[di_val]
            mol_calc = slab[mol_ind]
            mol_calc.calc = GPAW(
                                h=0.16,
                                kpts={'size': [4, 4, 1]},
                                occupations={'name': 'fermi-dirac', 'width': 0.05},
                                poissonsolver={'dipolelayer': 'xy'},
                                xc='BEEF-vdW',
                                txt=f"{di_val}.txt"
                            )
            if not os.path.isfile(f"{di_val}.gpw"):
                mol_calc.get_potential_energy()
                mol_calc.calc.write(f"{di_val}.gpw", mode='all')
            c_mol = GPAW(f"{di_val}.gpw")
            
elif os.path.isfile("relax.gpw"):
    mol, calc = restart("relax.gpw", txt="dos.txt")
    e, dos = calc.get_dos(spin=0, npts=2001, width=0.2)
    e_f = calc.get_fermi_level()
    fg1, ax1 = plt.subplots(1,1, figsize=(10,10))
    ax1.plot(e - e_f, dos)
    # ax1.set(xlim=(-20, 10), ylim=( None, 4))
    ax1.set_ylabel('DOS', fontsize=15)
    ax1.set_xlabel('Energy [eV]', fontsize=15)
    ax1.set_title("Molecule Total Density State", fontsize=20)
    ax1.tick_params(axis='both', which='major', labelsize=15)
    ax1.axvline(0, linewidth=2, linestyle='--')
    fg1.savefig('total_pdos.png')
    
else:
    print("Nothing to do DOS on")
