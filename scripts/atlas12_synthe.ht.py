#!/usr/bin/env python
import argparse
from pathlib import Path
import os
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--model', help='model name')
parser.add_argument('-kr', '--krun', help='working directory')
args = parser.parse_args()

wrk_dir = args.krun
model_name = args.model

def copy_atlas12(vturb, old_file, new_file):
    replace_vturb = f"'s|VTURB 1.00|VTURB {vturb}|g'"
    replace_krun = f"'s|krun={krun_default}|krun={krun}|g'"
    replace_synthe = f"'s|$khome/synthe/synthe.sh|$krun/synthe/synthe.sh|g'"
    SED_cmd = f"sed -e {replace_krun} -e {replace_vturb} -e {replace_synthe} {old_file} > {new_file}"
    os.system(SED_cmd)
    return


def copy_synthe(old_file, new_file):
    replace_krun = f"'s|krun={krun_default}|krun={krun}|g'"
    SED_cmd = f"sed -e {replace_krun} {old_file} > {new_file}"
    os.system(SED_cmd)
    new_file.chmod(33261)  # add permission to execute
    return


def run_atlas12(atlas12_file, model):
    input_dir = model.name
    output_dir = f'at12_{model.name}'
    atlas12_cmd = f"source {atlas12_file} {output_dir} {input_dir}"
    os.system(atlas12_cmd)
    return


def generate_spectrum(model, label):
    # Prep Files
    print(f'Preparing model {model.name}')
    atlas12_tmp = krun.joinpath(f'atlas12/atlas12_temp_{model.name}.sh')
    with np.load(label) as data:
        vturb = f"{data['v_micro']:1.2f}"
    copy_atlas12(vturb, atlas12, atlas12_tmp)
    # Run Atlas12
    print(f'Generating Spectra for model {model.name} (vturb = {vturb} km/s)')
    run_atlas12(atlas12_tmp, model)
    # Clean Up Files
    atlas12_tmp.unlink()
    print(f'Completed model {model.name}')
    return


def gen_spec_wrap(model):
    a = generate_spectrum(models[model], labels[model])
    return


kbase = Path('/global/scratch/nathan_sandford/kurucz')
khome = kbase.joinpath('kurucz_home')
krun_default = kbase.joinpath('kurucz_run')
krun = Path(wrk_dir)
grid = krun.joinpath('grids')
atlas12 = khome.joinpath('atlas12/atlas12.sh')
synthe = khome.joinpath('synthe/synthe.sh')
synthe_tmp = krun.joinpath(f'synthe/synthe.sh')

model = grid.joinpath(f'{model_name}')
label = model.joinpath(f'{model_name}_labels.npz')

at12_model = model.parent.joinpath(f'at12_{model_name}')

if model_name=='aaaaa':
    print(f'Working Directory: {krun}')

if list(at12_model.glob('spec/*')):
    print(f'Spectra already generated for model {model_name}')
else:
    if not synthe_tmp.is_file():
        copy_synthe(synthe, synthe_tmp)
    generate_spectrum(model, label)
    print(f'Model {args.model} completed')
