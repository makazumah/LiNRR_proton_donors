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
import pandas as pd
import csv, re
from ase import Atom
import glob
import os

from_scratch = False
c_path = os.path.abspath(os.getcwd())
in_traj= c_path+ '/input.traj'
if os.path.isfile("relax.gpw"):
    
    molecule, calc = restart("relax.gpw", txt="bader.txt")
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
                
        os.system("bader -p all_atom -p atom_index density.cube")
        while not os.path.exists("ACF.dat"):
            time.sleep(1)
        
        # while not os.path.exists("BCF.dat"):
        #     time.sleep(1)
        # while not os.path.exists("BvIndex.cube"):
        #     time.sleep(1)
            
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
elif os.path.isfile("output.gpw"):
    
    slab, calc = restart("output.gpw", txt="bader.txt")
    slab.get_potential_energy()
    hf = HirshfeldPartitioning(slab.calc)
    for atom, charge in zip(slab, hf.get_charges()):
        atom.charge = charge
    slab.copy().write('Hirshfeld.traj')
    rho = slab.calc.get_all_electron_density(gridrefinement=4)
    write('density.cube', slab, data=rho * Bohr**3)
    
    while not os.path.exists("Hirshfeld.traj"):
        time.sleep(1)
        
    if os.path.isfile("Hirshfeld.traj") and os.path.isfile("density.cube"):
        hirsh = read('Hirshfeld.traj')
        bader_out = hirsh.copy()
                
        os.system("bader -p all_atom -p atom_index density.cube")
        while not os.path.exists("ACF.dat"):
            time.sleep(1)
        
        # while not os.path.exists("BCF.dat"):
        #     time.sleep(1)
        # while not os.path.exists("BvIndex.cube"):
        #     time.sleep(1)
            
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

        if os.path.isfile('slab.pckl'):
            with open('slab.pckl', 'rb') as fd:
                dens, bader1, atoms = pickle.load(fd)
        else:
            dens, atoms = read_cube_data('density.cube')
            bader1, atoms = read_cube_data('AtIndex.cube')
            x = len(dens) // 2
            dens = dens[x]
            bader1 = bader1[x]
            with open('slab.pckl', 'wb') as fd:
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
        plt.savefig('slab-bader.png')
else:
    print("Not ready for Bader Analysis")
    
def ACF_Dat_processor(inpath):
    dat_file_loc = inpath+ '/ACF.dat'
    f = open(dat_file_loc, 'r')
    lines = f.read().splitlines()
    f.close()
    clean_lines= []
    for line in lines:
        if ("-"*10) not in line:
            clean_lines.append(line)
    dict_cols = clean_lines[0].split()
    dict_entries = np.array([clean_lines[ind].split() for ind in range(1, len(clean_lines)-3) ]).astype(float)
    dict_cols = [entry for entry in dict_cols if (entry.lower() != 'dist') and (entry.lower() != 'vol')]
    acf_dict = {dict_cols[ind]: dict_entries[:, ind].tolist() for ind in range(len(dict_cols))}
    acf_df = pd.DataFrame.from_dict(acf_dict)
    del acf_df["#"]
    acf_df.rename(columns={'X': 'X (Bohr)',
                            'Y': 'Y (Bohr)',
                            'Z': 'Z (Bohr)',
                            'MIN': 'Min Distance',
                            'ATOMIC': 'Atomic Volume'},
          inplace=True, errors='raise')
    return acf_df

def BaderTxt_processor(inpath):
    bader_loc = inpath +"bader.txt"
    words = ['Positions:', 'Unit cell:']
    bader_lines= []
    with open(bader_loc, 'r') as fp:
        # read all lines in a list
        lines = fp.readlines()
        add_lines = False
        for line in lines:
            # check if string present on a current line
            if line.strip() == words[0]:
                add_lines = True
                # print(words[0], 'string exists in file')
                # print('Line Number:', lines.index(line))
                
            if line.strip() == words[1]:
                add_lines = False
                # print(words[1], 'string exists in file')
                # print('Line Number:', lines.index(line))
                
            if add_lines:
                bader_lines.append(line)
            # if (line.find(words[1]) != -1) or (line.find(words[1]) != -1):
            #     print(words[0], 'string exists in file')
            #     print('Line Number:', lines.index(line))
            #     print('Line:', line)
            # if :
            #     print(words[1], 'string exists in file')
            #     print('Line Number:', lines.index(line))
            #     print('Line:', line)
        del bader_lines[-1]
        del bader_lines[0]
        
    fp.close()
    
    bader_entries = np.array([bader_lines[ind].split() for ind in range(len(bader_lines)) ])
    bader_cols = ['ASE Index','Element', 'X (Ang)','Y (Ang)','Z (Ang)']
    bader_dict = {bader_cols[ind]: bader_entries[:, ind].tolist() for ind in range(len(bader_cols))}
    bader_df = pd.DataFrame.from_dict(bader_dict)
    
    return bader_df

def bader_df_finalization(inpath, write_to_disk=False, overwrite=False, get_nn=4):
    # elements_df, _, = load_elements()
    if inpath[-1]!='/':
        inpath = inpath + '/'
    acf_df= ACF_Dat_processor(inpath)
    bader_df = BaderTxt_processor(inpath)
    full_bader_df = pd.merge(bader_df, acf_df, left_index=True, right_index=True)
    full_bader_df['AN']= 0
    elements_list = list(set(full_bader_df['Element']))
    for element in elements_list:
        An= Atom(element, (0,0,0))
        full_bader_df.loc[full_bader_df['Element']==element, 'AN'] = An.number
    full_bader_df['Charge Gain'] = full_bader_df['AN'] -full_bader_df['CHARGE']
    
    traj_path_1 = inpath + "relax.traj"
    traj_path_2 = inpath + "output.traj"

    if os.path.isfile(traj_path_1):
        # print("Single molecule. Not getting tags")
        slab = read(traj_path_1)
        full_bader_df['Tags'] = slab.get_tags()
    elif os.path.isfile(traj_path_2):
        print('Slab file found')
        slab = read(traj_path_2)
        full_bader_df['Tags'] = slab.get_tags()
    else:
        print('No traj file found')
    
    if get_nn!=0:
        for ii in range(get_nn):
            if ii > len(full_bader_df.index):
                print("Breaking out of loop")
                break
            else:
                the_nearest = slab.get_all_distances().argpartition(ii+1, axis=1)[:, ii+1]
                full_bader_df[f'{ii+1}NN']= (full_bader_df.loc[the_nearest, 'Element']+'-'+full_bader_df.loc[the_nearest, 'ASE Index'].astype(str)).tolist()
                full_bader_df[f'{ii+1}NN charge']= full_bader_df.loc[the_nearest, 'Charge Gain'].tolist()

        
    full_bader_df['Class'] =0
    full_bader_df.loc[full_bader_df['Charge Gain']<=0,'Class'] = 'Nucleophile'
    full_bader_df.loc[full_bader_df['Charge Gain']>0,'Class'] = 'Electrophile'
    
    if write_to_disk:
        outpath = inpath + "full_bader.txt"
        outpath_ex = inpath + "full_bader.xlsx"
        full_bader_df.to_excel(outpath_ex, index = False)
        with open(outpath, 'w') as out_file:
            df_string = full_bader_df.to_string()
            out_file.write(df_string)
            
        out_file.close()
    return full_bader_df


# bader_df_finalization(c_path, write_to_disk=True)
existence_path = c_path + '/full_bader.txt'
if os.path.isfile(existence_path):
    print("File Exists")
else:
    existence_path_2 = c_path + '/bader.txt'
    existence_path_3 = c_path + '/ACF.dat'
    
    if os.path.isfile(existence_path_2) and os.path.isfile(existence_path_3):
        
        print("Requirements Met")
        bader_df_finalization(c_path, write_to_disk=True)
    else:
        print("Couldn't complete file processing")
        