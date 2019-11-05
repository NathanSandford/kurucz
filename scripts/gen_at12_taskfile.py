import sys
from pathlib import Path

knum = sys.argv[1]

kbase = Path('/global/scratch/nathan_sandford/kurucz/')
krun = kbase.joinpath(f'kurucz_run/kurucz_run{knum}')
grid = krun.joinpath('grids')
taskfile = kbase.joinpath(f'jobs/tasks/krun{knum}.task')
scriptfile = kbase.joinpath('scripts/atlas12_synthe.ht.py')

models = sorted([model.name for model in grid.iterdir()])

file = open(taskfile, "w")
for model in models:
    file.write(f"{scriptfile} --krun {krun} --model {model}\n")
