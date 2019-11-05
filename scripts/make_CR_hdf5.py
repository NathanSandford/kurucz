#!/usr/bin/env python
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore", message="object name is not a valid Python identifier")

description = "Script to compile CRLB spectra"

'''
defaults
'''
kbase = Path('/global/scratch/nathan_sandford/kurucz')
kout = 'kurucz_out'
spec_dir = 'synthetic_spectra'
crlb_in = 'CRLB_spec' 
crlb_out = 'CRLB_data'
output_specfile = 'reference_spectra_300000.h5'
output_contfile = 'reference_continuum_300000.h5'
output_labelfile = 'reference_labels.h5'

'''
Parse Args
'''
parser = argparse.ArgumentParser(description=description)
parser.add_argument("--kbase", "-kb",
                    help=f"Base directory (default: {kbase})")
parser.add_argument("--kout", "-ko",
                    help=f"Path to Output Directory {kout}")
parser.add_argument("--spec_dir",
                    help=f"Path to Spectra Directory {spec_dir}")
parser.add_argument("--crlb_in", "-i",
                    help=f"Path to Output Directory {crlb_in}")
parser.add_argument("--crlb_out", "-o",
                    help=f"Path to Output Directory {crlb_out}")
parser.add_argument("--specfile",
                    help=f"Name of output spectra file {output_specfile}")
parser.add_argument("--contfile",
                    help=f"Name of output continuum file {output_contfile}")
parser.add_argument("--labelfile",
                    help=f"Name of output label file {output_labelfile}")
args = parser.parse_args()

if args.kbase:
    kbase = Path(args.kbase)
assert kbase.is_dir(), f'Directory {kbase} does not exist'

if args.kout:
    kout = kbase.joinpath(args.kout)
else:
    kout = kbase.joinpath(kout)
assert kout.is_dir(), f'Directory {kout} does not exist'

if args.spec_dir:
    spec_dir = kout.joinpath(args.spec_dir)
else:
    spec_dir = kout.joinpath(spec_dir)
assert spec_dir.is_dir(), f'Directory {spec_dir} does not exist'

if args.crlb_in:
    crlb_in = spec_dir.joinpath(args.crlb_in)
else:
    crlb_in = spec_dir.joinpath(crlb_in)
assert crlb_in.is_dir(), f'Directory {crlb_in} does not exist'

if args.crlb_out:
    crlb_out = kout.joinpath(args.crlb_out)
else:
    crlb_out = kout.joinpath(crlb_out)
assert crlb_out.is_dir(), f'Directory {crlb_out} does not exist'

if args.specfile:
    output_specfile = crlb_out.joinpath(args.specfile)
else:
    print(output_specfile)
    output_specfile = crlb_out.joinpath(output_specfile)
if args.contfile:
    output_contfile = crlb_out.joinpath(args.contfile)
else:
    output_contfile = crlb_out.joinpath(output_contfile)
if args.labelfile:
    output_labelfile = crlb_out.joinpath(args.labelfile)
else:
    output_labelfile = crlb_out.joinpath(output_labelfile)


def load_crlb_spectra_h5(reference, normalize=True):
    """
    Reads in and normalizes Kurucz model spectra from file
    :param reference:
    :return: name of reference spectra
    -high resolution wavelength grid
    -high resolution normalized spectra
    -labels of each spectra
    """
    # Highres Wavelength
    highres_wave = pd.read_hdf(reference, 'wavelength')
    # Unnormalized Spectra
    highres_spectra = pd.read_hdf(reference, 'spectra')
    # Theoretical Continuum
    highres_continuum = pd.read_hdf(reference, 'continuum')
    # Theoretically Normalized Spectra
    highres_norm_spec = highres_spectra / highres_continuum
    # Labels: Teff, logg, v_micro, [Li/H], [Be/H], ..., [Es/H]
    lab = pd.read_hdf(reference, 'labels')
    return highres_wave, highres_norm_spec, highres_continuum, lab


if output_specfile.exists():
    overwrite = input(f"{output_specfile} exists.\n Overwrite, 'yes'/'no'? [no]: ") or "no"
    if overwrite == 'no':
        print("Canceling Operation")
        sys.exit()
    else:
        output_specfile.unlink()
        output_contfile.unlink()
        output_labelfile.unlink()

ref_list = list(crlb_in.glob('*'))
ref_list_name = []
for i, ref in enumerate(ref_list):
    print(f'Reading in {ref.name} ({i+1}/{len(ref_list)})')
    if ref.name[-4:] == '.npz':
        ref_name = ref.name[:-4]
        temp = load_crlb_spectra_npz(ref)
        highres_wavelength = pd.DataFrame(temp[0].T)
        highres_norm_spectra = pd.DataFrame(temp[1].T)
        highres_continuum = pd.DataFrame(temp[2].T)
        labels = pd.DataFrame(temp[3].T)
    elif ref.name[-3:] == '.h5':
        ref_name = ref.name[:-3]
        temp = load_crlb_spectra_h5(ref)
        highres_wavelength = temp[0]
        highres_norm_spectra = temp[1]
        highres_continuum = temp[2]   
        labels = temp[3]
    else:
        print(f'{ref} is not a .npz or a .h5')
        sys.exit()
    highres_wavelength.to_hdf(output_specfile, 'highres_wavelength')
    highres_norm_spectra.to_hdf(output_specfile, ref_name)
    highres_continuum.to_hdf(output_contfile, ref_name)
    labels.to_hdf(output_labelfile, ref_name)
    ref_list_name.append(ref_name)

ref_list = pd.DataFrame(ref_list_name)
ref_list.to_hdf(output_specfile, 'ref_list')
print(f'Saved spectra to {output_specfile}')
ref_list.to_hdf(output_contfile, 'ref_list')
print(f'Saved continuum to {output_contfile}')
ref_list.to_hdf(output_labelfile, 'ref_list')
print(f'Saved labels to {output_labelfile}')

print('Completed')
