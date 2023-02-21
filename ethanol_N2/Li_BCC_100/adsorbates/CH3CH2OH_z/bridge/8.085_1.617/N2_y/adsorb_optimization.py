import os
import numpy as np

from gpaw import GPAW
from gpaw import restart
from ase.io import read
from ase.io import write
from ase.optimize import BFGS
from ase.parallel import barrier
from ase.eos import EquationOfState
from ase.vibrations import Vibrations

from_scratch = False
c_path = os.path.abspath(os.getcwd())
in_traj= c_path+ '/input.traj'
slab = read(in_traj)
# write("input.traj", slab)

if os.path.isfile("output.gpw"):
    try:
        slab, _ = restart("output.gpw", txt="output.txt")
    except:
        from_scratch = True
        if os.path.isfile("output.traj") and os.path.getsize("output.traj") > 0:
            slab = read("output.traj")
        else:
            slab = read("input.traj")
else:
    slab = read("input.traj") # Just an input structure. Just write slab to file
    from_scratch = True

if from_scratch:
    slab.calc = GPAW(
    h=0.16,
    kpts={'size': [4, 4, 1]},
    occupations={'name': 'fermi-dirac', 'width': 0.05},
    poissonsolver={'dipolelayer': 'xy'},
    xc='BEEF-vdW',
    txt="output.txt"
)

def relax(atoms, fmax=0.05, step=0.04):
    atoms.calc.attach(atoms.calc.write, 5, "output.gpw")

    def _check_file_exists(filename):
        #Check if file exists and is not empty
        if os.path.isfile(filename):
            return os.path.getsize(filename) > 0
        else:
            return False

    # check if it is a restart
    barrier()
    if _check_file_exists("output.traj"):
        latest = read("output.traj", index=":")
        # check if already restarted previously and extend history if needed
        if not (_check_file_exists("history.traj")):
            barrier()
            write("history.traj", latest)
        else:
            hist = read("history.traj", index=":")
            hist.extend(latest)
            write("history.traj", hist)

    dyn = BFGS(atoms=atoms, trajectory="output.traj", logfile="qn.log", maxstep=step)
    # if history exists, read in hessian
    if _check_file_exists("history.traj"):
        dyn.replay_trajectory("history.traj")
    # optimize
    dyn.run(fmax=fmax)

relax(slab)
print(f'The energy is given by: {slab.get_potential_energy():.3f} eV')
print('Code run to end')