#!/usr/bin/env python
from pathlib import Path

krun_base = Path('/global/scratch/nathan_sandford/kurucz/kurucz_run')
krun_list = sorted(list(krun_base.glob('*')))

for krun in krun_list:
    n_spec = len(list(krun.glob('grids/at12_*/spec/*')))
    n_at12 = len(list(krun.glob('grids/at12*')))
    n_total = len(list(krun.glob('grids/?[!t]*')))
    print(f'{krun.name}: {n_spec}/{n_at12}/{n_total}')
