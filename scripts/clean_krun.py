#!/usr/bin/env python
import sys
from pathlib import Path
from shutil import rmtree

krun_base = Path('/global/scratch/nathan_sandford/kurucz/kurucz_run')
krun = krun_base.joinpath(sys.argv[1])

atlas12_dir = krun.joinpath('atlas12')
synthe_dir = krun.joinpath('synthe')

print(f'Removing all contents of {atlas12_dir}')
rmtree(atlas12_dir)
atlas12_dir.mkdir()

print(f'Removing all contents of {synthe_dir}')
rmtree(synthe_dir)
synthe_dir.mkdir()

print('Removing all at12_XXXXX directories')
at12_models = list(krun.glob('grids/at12*'))
[rmtree(at12_dir) for i in at12_models]
