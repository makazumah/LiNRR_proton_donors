from ase.io import read, write
import numpy as np
from ase.visualize import view
import os

c_dir = os.getcwd()
c_path = "/home/vazumah/Python_projects/LiNRR/Li_adsorb_studies/proton_donors/N2_with_proton_donor/Li_BCC_100/adsorbates/HBr/bridge/8.085_1.617/N2_z"
slab = read(c_path + "/output.traj")
hbr_pos = slab.get_positions()[[36, 37],:]
hbr_pos[:,2]= hbr_pos[:,2]-2.
new_slab_pos = slab.get_positions()
new_slab_pos[[36,37],:]=hbr_pos
slab.set_positions(new_slab_pos)
view(slab)
os.chdir(c_path)
write("input.traj", slab)
os.chdir(c_dir)