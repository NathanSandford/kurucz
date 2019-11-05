import argparse
import os
import sys
from pathlib import Path
from shutil import rmtree
import numpy as np
import pandas as pd
from multiprocessing import cpu_count, Pool
from mendeleev import element
from tqdm import tqdm

description = "Script to handle output of ATLAS12+SYNTHE runs"

'''
Defaults
'''
kbase = Path('/global/scratch/nathan_sandford/kurucz')
krun_base = kbase.joinpath('kurucz_run')
kout = kbase.joinpath('kurucz_out')
gen_grid_dir = kout.joinpath('generated_grids')
gen_spec_dir = kout.joinpath('synthetic_spectra')
input_name = 'grids'

'''
Parse Args
'''
parser = argparse.ArgumentParser(description=description)
parser.add_argument("krun", help="kurucz_run sub-directory")
parser.add_argument("output_name", help="Name of output directory")
parser.add_argument("--check_spec", "-check", action="store_true", default=False,
                    help="Checks if spectra all have the same length (default: False)")
parser.add_argument("--sort_output", "-sort", action="store_true", default=False,
                    help="Sorts output into generated_spectra directory (default: False)")
parser.add_argument("--combine_spec", "-combine", action="store_true", default=False,
                    help="Combines all spectra & labels into hdf5 file (default: False)")
parser.add_argument("--kbase", "-kb",
                    help=f"Base directory (default: {kbase})")
parser.add_argument("--krun_base", "-krb",
                    help=f"Path to Directory where code runs {krun_base.name}")
parser.add_argument("--kout", "-ko",
                    help=f"Path to Output Directory {kout.name}")
parser.add_argument("--input_name", "-i", help=f"Name of input directory {input_name}")
parser.add_argument("--gen_grid_dir", "-ggd", help=f"Directory for generated grids {gen_grid_dir}")
parser.add_argument("--gen_spec_dir", "-gsd", help=f"Directory for Synthetic Spectra {gen_spec_dir}")
args = parser.parse_args()

if args.kbase:
    kbase = Path(args.kbase)
assert kbase.is_dir(), f'Directory {kbase} does not exist'

if args.krun_base:
    krun_base = kbase.joinpath(args.krun_base)
assert krun_base.is_dir(), f'Directory {krun_base} does not exist'

if args.kout:
    kout = kbase.joinpath(args.kout)
assert kout.is_dir(), f'Directory {kout} does not exist'

if args.gen_grid_dir:
    gen_grid_dir = kout.joinpath(args.gen_grid_dir)
assert gen_grid_dir.is_dir(), f'Directory {gen_grid_dir} does not exist'

if args.gen_spec_dir:
    gen_spec_dir = kout.joinpath(args.gen_spec_dir)
assert gen_spec_dir.is_dir(), f'Directory {gen_spec_dir} does not exist'

krun = krun_base.joinpath(args.krun)

if args.input_name:
    input_name = args.input_name
input_dir = krun.joinpath(input_name)

output_dir = gen_grid_dir.joinpath(args.output_name)
spec_file = gen_spec_dir.joinpath(args.output_name+'.h5')

'''
Check Spectra
'''
if args.check_spec:
    def get_output_shape(spec):
        try:
            tmp = np.genfromtxt(spec)
            length = tmp.shape[0]
        except:
            length = 0
        return length

    list_spec = list(input_dir.glob('*/spec/*'))
    n_models = len(list_spec)
    print(f'Beginning batch reading of {n_models} spectra w/ {cpu_count()} CPUs')
    with Pool(cpu_count()) as pool:
        lengths = np.array(list(tqdm(pool.imap(get_output_shape, list_spec), total=n_models)))
        n_bad_spec = np.sum(lengths!=lengths.max())
    if n_bad_spec:
        delete = input(f"{n_bad_spec} bad spectra detected. " \
                       + "Remove, 'yes'/'no'? [yes]: ") or "yes"
        if delete == 'yes':
            [rmtree(list_spec[i].parents[1]) for i in range(n_models) if lengths[i] != lengths.max()]
            print('Removed bad spectra.')
            if args.sort_output or args.combine_spec:
                cont = input("Continue sorting output, 'yes'/'no'? [no]: ") or "no"
                if cont == 'no':
                    print("Canceling operation")
                    sys.exit()
    else:
        print('No bad spectra detected!')

'''
Sort Output
'''
if args.sort_output:
    if output_dir.is_dir():
        overwrite = input(f"{output_dir} exists.\n Overwrite, 'yes'/'no'? [no]: ") or "no"
        if overwrite == 'no':
            print("Canceling Operation")
            sys.exit()
        else:
            rmtree(output_dir)
            output_dir.mkdir()

    spec_dir = output_dir.joinpath('spec')
    label_dir = output_dir.joinpath('labels')
    atm_dir = output_dir.joinpath('atm')
    if not output_dir.is_dir():
        output_dir.mkdir()
    if not spec_dir.is_dir():
        spec_dir.mkdir()
    if not label_dir.is_dir():
        label_dir.mkdir()
    if not atm_dir.is_dir():
        atm_dir.mkdir()

    print('Copying Spectra Files...')
    os.system(f"cp {input_dir}/at*/spec/* {spec_dir}")
    print('Copying Atmosphere Files...')
    os.system(f"cp {input_dir}/at*/atm/* {atm_dir}")
    print('Copying Label Files...')
    os.system(f"cp {input_dir}/*/*_labels.npz {label_dir}")
    #print('Zipping Grid Directory')
    #os.chdir(krun)
    #os.system(f"zip -rqm {output_dir}/grids.zip ./{input_name}")
    print("Sorting Completed!")

'''
Combine Spectra
'''
if not args.combine_spec:
    sys.exit()

spec_dir = output_dir.joinpath('spec')
label_dir = output_dir.joinpath('labels')

list_spec = sorted(list(spec_dir.glob('*')))
list_models = [spec_file.name[5:10] for spec_file in list_spec]
list_labels = sorted(list(label_dir.glob('*')))
n_models = len(list_models)

# Array containing Solar elemental abundances from H -> Es
# H and He are fractional abundances while the rest are log10[X/H]
# Added filler 0 to line up indices with atomic number
print('Defining Solar abundance array')
label_names = ['Teff', 'logg', 'v_micro'] + [x.symbol for x in element(list(range(3,100)))]
element_array_solar = np.array([0, 0.92068, 0.07837,
                                -10.99, -10.66, -9.34, -3.61, -4.21,
                                -3.35, -7.48, -4.11, -5.80, -4.44,
                                -5.59, -4.53, -6.63, -4.92, -6.54,
                                -5.64, -7.01, -5.70, -8.89, -7.09,
                                -8.11, -6.40, -6.61, -4.54, -7.05,
                                -5.82, -7.85, -7.48, -9.00, -8.39,
                                -9.74, -8.70, -9.50, -8.79, -9.52,
                                -9.17, -9.83, -9.46, -10.58, -10.16,
                                -20.00, -10.29, -11.13, -10.47, -11.10,
                                -10.33, -11.24, -10.00, -11.03, -9.86,
                                -10.49, -9.80, -10.96, -9.86, -10.94,
                                -10.46, -11.32, -10.62, -20.00, -11.08,
                                -11.52, -10.97, -11.74, -10.94, -11.56,
                                -11.12, -11.94, -11.20, -11.94, -11.19,
                                -12.16, -11.19, -11.78, -10.64, -10.66,
                                -10.42, -11.12, -10.87, -11.14, -10.29,
                                -11.39, -20.00, -20.00, -20.00, -20.00,
                                -20.00, -20.00, -12.02, -20.00, -12.58,
                                -20.00, -20.00, -20.00, -20.00, -20.00,
                                -20.00, -20.00])


def read_spec(i):
    spec = list_spec[i]
    mod = list_models[i]
    label = list_labels[i]
    # Read in spectra file
    tmp = np.genfromtxt(spec_dir.joinpath(f'{spec}'))
    tmp = tmp.transpose()
    # Read in Wavelength, Spectrum, and Theoretical Continuum
    wave = tmp[0] * 10  # AA
    spec = tmp[1]
    cont = tmp[2]
    # Read in labels file output by generate_initial_atm.py
    tmp = np.load(label_dir.joinpath(f'{label}'))
    teff = tmp['teff']
    logg = tmp['logg']
    v_micro = tmp['v_micro']
    # Add Filler 0
    element_array_unscaled = np.concatenate([[0], tmp['elems_array']])
    tmp.close()
    # Scale abundances to solar (except H & He)
    element_array = element_array_unscaled
    element_array[3:] = element_array_unscaled[3:] - element_array_solar[3:]
    # Select which labels to keep
    all_labs = np.concatenate([[teff], [logg], [v_micro], element_array[3:]])
    return wave, spec, cont, all_labs, mod


print('Beginning batch compiling of spectra')
with Pool(cpu_count()) as pool:
    temp = list(tqdm(pool.imap(read_spec, range(n_models)), total=n_models))

wavelength = pd.DataFrame(temp[0][0])
spec = pd.DataFrame.from_dict({_[4]:_[1] for _ in temp})
cont = pd.DataFrame.from_dict({_[4]:_[2] for _ in temp})
labels = pd.DataFrame.from_dict({_[4]:_[3] for _ in temp})
labels.index = label_names

wavelength.to_hdf(spec_file, 'wavelength')
spec.to_hdf(spec_file, 'spectra')
cont.to_hdf(spec_file, 'continuum')
labels.to_hdf(spec_file, 'labels')
print(f'Saved synthetic spectra to {spec_file}')
