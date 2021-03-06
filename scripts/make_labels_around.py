import argparse
import numpy as np
import pandas as pd
from mendeleev import element
from pathlib import Path

description = \
    '''
    Code generates ball of labels around a specific label
    '''

'''
Defaults
'''
kbase = Path('/global/scratch/nathan_sandford/kurucz/')
khome = kbase.joinpath('kurucz_home/')
kinput = kbase.joinpath('input')
input_file = kinput.joinpath('labels.h5')
dX_teff = 250
dX_logg = 0.3
dX_v_micro = 0.5
dX_abundances = 0.5
free_abundances = ['Fe', 'Ca', 'Ni', 'Si', 'Ti', 'Co', 'Mg']

parser = argparse.ArgumentParser(description=description)
parser.add_argument("key", help="key corresponding to set of labels in label file")
parser.add_argument("N", type=int, help="Number of labels to generate")
parser.add_argument("teff", metavar="T_eff", type=float, help="Effective Temperature (K)")
parser.add_argument("logg", metavar="log(g)", type=float, help="Surface Gravity")
parser.add_argument("feh", metavar="[Fe/H]", type=float, help="Solar Scaled Iron Abundance")
parser.add_argument("v_micro", type=float, help="Micro-turbulent Velocity (km/s)")
parser.add_argument("--kbase", "-kb",
                    help=f"Base directory (default: {kbase})")
parser.add_argument("--khome", "-kh",
                    help=f"Path to Kurucz Code (default: {khome.name})")
parser.add_argument("--kinput", "-ki", help=f"Input Directory (default: {kinput.name})")
parser.add_argument("--inputfile", "-i", help=f"File containing labels (default: {input_file.name})")
parser.add_argument("--dX_teff", "-dT", type=float,
                    help=f"Radius in T_eff (default: {dX_teff} K)")
parser.add_argument("--dX_logg", "-dg", type=float,
                    help=f"Radius in log(g) (default: {dX_logg})")
parser.add_argument("--dX_v_micro", "-dv", type=float,
                    help=f"Radius in v_micro (default: {dX_v_micro} km/s)")
parser.add_argument("--dX_abundances", "-dX", type=float,
                    help=f"Radius in [X/H] (default: {dX_abundances} dex)")
parser.add_argument("--free_abundances", "-X", nargs='+',
                    help=f"Free abundances (default: {' '.join(free_abundances)})")
args = parser.parse_args()

if args.kbase:
    kbase = Path(args.kbase)
assert kbase.is_dir(), f'Directory {kbase} does not exist'

if args.khome:
    khome = kbase.joinpath(args.khome)
assert khome.is_dir(), f'Directory {khome} does not exist'

if args.kinput:
    kinput = kbase.joinpath(args.kinput)
assert kinput.is_dir(), f'Directory {kinput} does not exist'

if args.inputfile:
    input_file = kinput.joinpath(args.inputfile)
else:
    input_file = kinput.joinpath('labels.h5')
assert input_file.exists(), f'{input_file} does not exist'

if args.dX_teff:
    dX_teff = args.dX_teff
if args.dX_logg:
    dX_logg = args.dX_logg
if args.dX_v_micro:
    dX_v_micro = args.dX_v_micro
if args.dX_abundances:
    dX_abundances = args.dX_abundances
if args.free_abundances:
    free_abundances = args.free_abundances

key = args.key
N = args.N
teff = args.teff
logg = args.logg
feh = args.feh
v_micro = args.v_micro  # km/s

# Initialize Label Dataframe
label_names = ['Teff', 'logg', 'v_micro'] + [x.symbol for x in element(list(range(3,100)))]
labels = np.zeros((N, len(label_names)))
labels = pd.DataFrame(labels, columns=label_names)
labels['Teff'] = teff
labels['logg'] = logg
labels['v_micro'] = v_micro
labels['Fe'] = feh

# Add offsets
labels['Teff'] += np.random.uniform(low=-dX_teff, high=dX_teff, size=N)
labels['logg'] += np.random.uniform(low=-dX_logg, high=dX_logg, size=N)
labels['v_micro'] += np.random.uniform(low=-dX_v_micro, high=dX_v_micro, size=N)
free_abundances = ['Fe', 'Ca', 'Ni', 'Si', 'Ti', 'Co', 'Mg']
labels[free_abundances] += np.random.uniform(low=-dX_abundances, high=dX_abundances,
                                             size=(N, len(free_abundances)))

# Save to labels.h5
labels.to_hdf(input_file, key)
