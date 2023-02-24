from ase.io import read
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase.parallel import barrier
from gpaw import GPAW
import os

import json
import numpy as np

c_path = os.path.abspath(os.getcwd())
out_traj= c_path+ '/output.traj'
slab = read(out_traj)

slab.calc = GPAW(
    h=0.16,
    kpts={'size': [6, 6, 1]},
    occupations={'name': 'fermi-dirac', 'width': 0.05},
    poissonsolver={'dipolelayer': 'xy'},
    xc='RPBE',
    txt="output.txt",
    symmetry={'point_group': False}
)

def calculate_vibrations(atoms):
    indices_to_perturb = np.where(atoms.get_tags() <= 0)[0].tolist()
    vib = Vibrations(atoms, indices = indices_to_perturb, name = "vib.log")
    barrier()
    vib.clean(empty_files=True)
    vib.run()
    vib.summary(log="vib.summary")

    all_eners = vib.get_energies()

    eners = np.real(vib.get_energies())
    eners = eners[eners != 0]
    return eners

eners = calculate_vibrations(slab)
with open("eners.json", "w") as j:
    json.dump(list(eners), j)