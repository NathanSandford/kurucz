import argparse
from pathlib import Path
import pandas as pd
import numpy as np

description = "Combine spectra from several files into one file"

'''
Defaults
'''
kbase = Path('/global/scratch/nathan_sandford/kurucz')
krun_base = 'kurucz_run'
kout = 'kurucz_out'
spec_dir = 'synthetic_spectra'

'''
Parse Args
'''
parser = argparse.ArgumentParser(description=description)
parser.add_argument("input_dir", help="Directory of spectra files to combine")
parser.add_argument("--kbase", "-kb",
                    help=f"Base directory (default: {kbase})")
parser.add_argument("--kout", "-ko",
                    help=f"Path to Output Directory (default: {kout})")
parser.add_argument("--spec_dir",
                    help=f"Path to Spectra Directory (default: {spec_dir})")
parser.add_argument("--output_file", "-o",
                    help="Name of output file (default: {input_dir}.h5)")
args = parser.parse_args()

if args.kbase:
    kbase = Path(args.kbase)
assert kbase.is_dir(), f'Base directory {kbase} does not exist'

if args.kout:
    kout = kbase.joinpath(args.kout)
else:
    kout = kbase.joinpath(kout)
assert kout.is_dir(), f'Output directory {kout} does not exist'

if args.spec_dir:
    spec_dir = kout.joinpath(args.spec_dir)
else:
    spec_dir = kout.joinpath(spec_dir)
assert spec_dir.is_dir(), f'Spectra directory {spec_dir} does not exist'

input_dir = spec_dir.joinpath(args.input_dir)
assert input_dir.is_dir(), f'Input directory {input_dir} does not exist'

if args.output_file:
    if args.output_file[-3:] != '.h5':
        output_file = spec_dir.joinpath(args.output_file+'.h5')
    else:
        output_file = spec_dir.joinpath(args.output_file)
else:
    output_file = spec_dir.joinpath(args.input_dir+f'{input_dir.name}.h5')


print(input_dir.name)
'''
Combine Spectra Files
'''
file_list = list(input_dir.glob('*'))

for i, specfile in enumerate(file_list):
    if i == 0:
        wavelength = pd.read_hdf(specfile, 'wavelength')
        spec = pd.read_hdf(specfile, 'spectra')
        spec.columns = [f'{col}{i}' for col in spec.columns]
        labels = pd.read_hdf(specfile, 'labels')
        labels.columns = [f'{col}{i}' for col in labels.columns]
    else:
        assert np.all(wavelength == pd.read_hdf(specfile, 'wavelength')), \
            f'Wavelength of {specfile.name} does not match wavelength of {file_list[0].name}'
        tmp_spec = pd.read_hdf(specfile, 'spectra')
        tmp_spec.columns = [f'{col}{i}' for col in tmp_spec.columns]
        tmp_labels = pd.read_hdf(specfile, 'labels')
        tmp_labels.columns = [f'{col}{i}' for col in tmp_labels.columns]
        spec = pd.concat([spec, tmp_spec], axis=1)
        labels = pd.concat([labels, tmp_labels], axis=1)

print(f'Saving to {output_file}...')
wavelength.to_hdf(output_file, 'wavelength')
spec.to_hdf(output_file, 'spectra')
labels.to_hdf(output_file, 'labels')
print('Done!')
