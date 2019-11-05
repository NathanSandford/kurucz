import argparse
from pathlib import Path
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from mendeleev import element

description = \
    '''
    Code generates initial atmospheres for a set of stellar labels in labels.h5
    '''

'''
Defaults
'''
kbase = Path('/global/scratch/nathan_sandford/kurucz/')
khome = kbase.joinpath('kurucz_home/')
krun_base = kbase.joinpath('kurucz_run/')
kinput = kbase.joinpath('input')
input_file = kinput.joinpath('labels.h5')
katm = khome.joinpath('Sync_Spectra_All_Atm/')
name = 'grids'
verbose = 1

parser = argparse.ArgumentParser(description=description)
parser.add_argument("key", help="key corresponding to set of labels in labels.h5")
parser.add_argument("krun", help="kurucz_run sub-directory")
parser.add_argument("--kbase", "-kb",
                    help=f"Base directory (default: {kbase})")
parser.add_argument("--khome", "-kh",
                    help=f"Path to Kurucz Code (default: {khome.name})")
parser.add_argument("--krun_base", "-krb",
                    help=f"Path to Directory where code runs {krun_base.name}")
parser.add_argument("--kinput", "-ki", help=f"Input Directory (default: {kinput.name})"\
)
parser.add_argument("--inputfile", "-i", help=f"File containing labels (default: {input_file.name})")
parser.add_argument("--katm", "-ka", help=f"Directory of initial atm files (default: {katm.name})")
parser.add_argument("--name", "-n",
                    help=f"Directory in krun to save atmospheres (default: {name})")
parser.add_argument("--verbose", "-V", type=int,
                    help=f"Level of verbosity. 0=none, 1=summary, 2=all (default: {verbose})")
args = parser.parse_args()

if args.kbase:
    kbase = Path(args.kbase)
assert kbase.is_dir(), f'Directory {kbase} does not exist'

if args.khome:
    khome = kbase.joinpath(args.khome)
assert khome.is_dir(), f'Directory {khome} does not exist'

if args.krun_base:
    krun_base = kbase.joinpath(args.krun_base)
assert krun_base.is_dir(), f'Directory {krun_base} does not exist'

if args.kinput:
    kinput = kbase.joinpath(args.kinput)
assert kinput.is_dir(), f'Directory {kinput} does not exist'

if args.inputfile:
    input_file = kinput.joinpath(args.inputfile)
else:
    input_file = kinput.joinpath('labels.h5')
assert input_file.exists(), f'{input_file} does not exist'

if args.katm:
    katm = khome.joinpath(args.katm)
assert katm.is_dir(), f'Directory {kinput} does not exist'
    
krun = krun_base.joinpath(args.krun)
assert krun.is_dir(), f'Directory {krun} does not exist'

if args.name:
    name = args.name
output_dir = krun.joinpath(name)
if output_dir.is_dir():
    overwrite = input(f"{output_dir} exists.\n Overwrite, 'yes'/'no'? [no]: ") or "no"
    if overwrite == 'no':
        print("Canceling Operation")
        sys.exit()
    else:
        output_dir.rmdir()
        output_dir.mkdir()
    
if args.verbose:
    verbose = args.verbose

'''Restore Labels'''
labels = pd.read_hdf(input_file, args.key)
num_models = labels.shape[0]
if verbose:
    print('Restored labels for %i models' % num_models)

'''
Initialize Input Arrays w/ Solar Labels:
First two are Teff and logg
Starts with Li (H and He are added separately).
All elements are in the form [X/H].
'''
element_array = np.array([4750., 2.5,
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
element_array = np.repeat(np.array([element_array]), num_models, axis=0).T

'''
Change teff & logg
'''
element_array[0, :] = labels['Teff']
element_array[1, :] = labels['logg']

'''
Scale global metallicity
'''
element_array[2:, :] = element_array[2:, :] + labels['Fe'][np.newaxis, :]

'''
Scale Abundances
'''
feh_temp = labels['Fe']
labels['Fe'] = 0  # avoid doubly adding [Fe/H]
# noinspection PyTypeChecker
elements = [x.symbol for x in element(list(range(3, 100)))]
element_array[2:, :] += labels[elements].T
labels['Fe'] = feh_temp

'''
Renormalize Hydrogen such that X+Y+Z=1
'''
solar_He = 0.07837
solar_He_5s = "{0:.5f}".format(solar_He)
renormed_H = 1. - solar_He - np.sum(10. ** element_array[2:, :], axis=0)

'''
Make Input Kurucz Readable
'''
elems_array_0s = np.copy(element_array).astype("str")
elems_array_2s = np.copy(element_array).astype("str")
elems_array_3s = np.copy(element_array).astype("str")
elems_array_4s = np.copy(element_array).astype("str")
for p1 in range(element_array.shape[0]):
    for p2 in range(element_array.shape[1]):
        elems_array_0s[p1, p2] = '' + "{0:.0f}".format(element_array[p1, p2])
        elems_array_2s[p1, p2] = ' ' + "{0:.2f}".format(element_array[p1, p2])
        elems_array_3s[p1, p2] = ' ' + "{0:.3f}".format(element_array[p1, p2])
        elems_array_4s[p1, p2] = '' + "{0:.4f}".format(element_array[p1, p2])

        if element_array[p1, p2] <= -9.995:
            elems_array_2s[p1, p2] = elems_array_2s[p1, p2][1:]

        if element_array[p1, p2] < -9.9995:
            elems_array_3s[p1, p2] = elems_array_3s[p1, p2][1:]


'''
Restore Grid of Converged Atmospheres
'''
teff_find = []
logg_find = []
atm_find = []
name_find = []

if verbose > 1:
    print('Restoring generated grid from Sync_Spectra_All_Atm')
list_files_find = os.listdir(katm)
for f1 in range(len(list_files_find)):
    name_find.append(list_files_find[f1][5:10])
    teff_find.append(float(list_files_find[f1][12:17]))
    logg_find.append(float(list_files_find[f1][18:22]))
    atm_find.append(list_files_find[f1])

teff_find = np.array(teff_find)
logg_find = np.array(logg_find)
atm_find = np.array(atm_find)
name_find = np.array(name_find)


'''
Find the closest teff and logg as an initial guess.
Add some jitters to overcome convergence problem.
'''
if verbose > 1:
    print('Finding closest teff and logg from grid')
atm_closest = []
for p1 in range(num_models):
    ind_min = np.argmin(np.abs(teff_find - element_array[0, p1])
                        + 1000. * np.abs(logg_find - element_array[1, p1]))
    atm_closest.append(atm_find[ind_min])


'''
Modifying closest atmosphere file
'''
if verbose > 1:
    print('Looping over all models')
for c1 in tqdm(range(num_models), desc='Writing Atms'):
    if verbose > 1:
        print('Modifying closest atmosphere file: {}'.format(atm_closest[c1]))
    lines = open(katm.joinpath(atm_closest[c1])).readlines()
    template = khome.joinpath('template_1.atm')
    open(template, 'w').writelines(lines[44:])

    f = open(template, 'r')
    end_file = f.read()
    f.close()

    renormed_H_5s = f"{renormed_H[c1]:.5f}"

    '''
    Change Abundances
    '''
    start_file = 'TEFF   ' + elems_array_0s[0, c1] \
        + '.  GRAVITY  ' + elems_array_4s[1, c1] + ' LTE \n' \
        + 'TITLE ATLAS12    ' \
        + '                                                               \n' \
        + ' OPACITY IFOP 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 1 0 0 0\n' \
        + ' CONVECTION ON   1.25 TURBULENCE OFF  0.00  0.00  0.00  0.00\n' \
        + 'ABUNDANCE SCALE   1.00000 ABUNDANCE CHANGE 1 ' \
        + renormed_H_5s + ' 2 ' + solar_He_5s + '\n' \
        + ' ABUNDANCE CHANGE  3 ' + elems_array_2s[2, c1] + '  4 ' \
        + elems_array_2s[3, c1] + '  5 ' + elems_array_2s[4, c1] + '  6 ' \
        + elems_array_2s[5, c1] + '  7 ' + elems_array_2s[6, c1] + '  8 ' \
        + elems_array_2s[7, c1] + '\n' \
        + ' ABUNDANCE CHANGE  9 ' + elems_array_2s[8, c1] + ' 10 ' \
        + elems_array_2s[9, c1] + ' 11 ' + elems_array_2s[10, c1] + ' 12 ' \
        + elems_array_2s[11, c1] + ' 13 ' + elems_array_2s[12, c1] + ' 14 ' \
        + elems_array_2s[13, c1] + '\n' \
        + ' ABUNDANCE CHANGE 15 ' + elems_array_2s[14, c1] + ' 16 ' \
        + elems_array_2s[15, c1] + ' 17 ' + elems_array_2s[16, c1] + ' 18 ' \
        + elems_array_2s[17, c1] + ' 19 ' + elems_array_2s[18, c1] + ' 20 ' \
        + elems_array_2s[19, c1] + '\n' \
        + ' ABUNDANCE CHANGE 21 ' + elems_array_2s[20, c1] + ' 22 ' \
        + elems_array_2s[21, c1] + ' 23 ' + elems_array_2s[22, c1] + ' 24 ' \
        + elems_array_2s[23, c1] + ' 25 ' + elems_array_2s[24, c1] + ' 26 ' \
        + elems_array_2s[25, c1] + '\n' \
        + ' ABUNDANCE CHANGE 27 ' + elems_array_2s[26, c1] + ' 28 ' \
        + elems_array_2s[27, c1] + ' 29 ' + elems_array_2s[28, c1] + ' 30 ' \
        + elems_array_2s[29, c1] + ' 31 ' + elems_array_2s[30, c1] + ' 32 ' \
        + elems_array_2s[31, c1] + '\n' \
        + ' ABUNDANCE CHANGE 33 ' + elems_array_2s[32, c1] + ' 34 ' \
        + elems_array_2s[33, c1] + ' 35 ' + elems_array_2s[34, c1] + ' 36 ' \
        + elems_array_2s[35, c1] + ' 37 ' + elems_array_2s[36, c1] + ' 38 ' \
        + elems_array_2s[37, c1] + '\n' \
        + ' ABUNDANCE CHANGE 39 ' + elems_array_2s[38, c1] + ' 40 ' \
        + elems_array_2s[39, c1] + ' 41 ' + elems_array_2s[40, c1] + ' 42 ' \
        + elems_array_2s[41, c1] + ' 43 ' + elems_array_2s[42, c1] + ' 44 ' \
        + elems_array_2s[43, c1] + '\n' \
        + ' ABUNDANCE CHANGE 45 ' + elems_array_2s[44, c1] + ' 46 ' \
        + elems_array_2s[45, c1] + ' 47 ' + elems_array_2s[46, c1] + ' 48 ' \
        + elems_array_2s[47, c1] + ' 49 ' + elems_array_2s[48, c1] + ' 50 ' \
        + elems_array_2s[49, c1] + '\n' \
        + ' ABUNDANCE CHANGE 51 ' + elems_array_2s[50, c1] + ' 52 ' \
        + elems_array_2s[51, c1] + ' 53 ' + elems_array_2s[52, c1] + ' 54 ' \
        + elems_array_2s[53, c1] + ' 55 ' + elems_array_2s[54, c1] + ' 56 ' \
        + elems_array_2s[55, c1] + '\n' \
        + ' ABUNDANCE CHANGE 57 ' + elems_array_2s[56, c1] + ' 58 ' \
        + elems_array_2s[57, c1] + ' 59 ' + elems_array_2s[58, c1] + ' 60 ' \
        + elems_array_2s[59, c1] + ' 61 ' + elems_array_2s[60, c1] + ' 62 ' \
        + elems_array_2s[61, c1] + '\n' \
        + ' ABUNDANCE CHANGE 63 ' + elems_array_2s[62, c1] + ' 64 ' \
        + elems_array_2s[63, c1] + ' 65 ' + elems_array_2s[64, c1] + ' 66 ' \
        + elems_array_2s[65, c1] + ' 67 ' + elems_array_2s[66, c1] + ' 68 ' \
        + elems_array_2s[67, c1] + '\n' \
        + ' ABUNDANCE CHANGE 69 ' + elems_array_2s[68, c1] + ' 70 ' \
        + elems_array_2s[69, c1] + ' 71 ' + elems_array_2s[70, c1] + ' 72 ' \
        + elems_array_2s[71, c1] + ' 73 ' + elems_array_2s[72, c1] + ' 74 ' \
        + elems_array_2s[73, c1] + '\n' \
        + ' ABUNDANCE CHANGE 75 ' + elems_array_2s[74, c1] + ' 76 ' \
        + elems_array_2s[75, c1] + ' 77 ' + elems_array_2s[76, c1] + ' 78 ' \
        + elems_array_2s[77, c1] + ' 79 ' + elems_array_2s[78, c1] + ' 80 ' \
        + elems_array_2s[79, c1] + '\n' \
        + ' ABUNDANCE CHANGE 81 ' + elems_array_2s[80, c1] + ' 82 ' \
        + elems_array_2s[81, c1] + ' 83 ' + elems_array_2s[82, c1] + ' 84 ' \
        + elems_array_2s[83, c1] + ' 85 ' + elems_array_2s[84, c1] + ' 86 ' \
        + elems_array_2s[85, c1] + '\n' \
        + ' ABUNDANCE CHANGE 87 ' + elems_array_2s[86, c1] + ' 88 ' \
        + elems_array_2s[87, c1] + ' 89 ' + elems_array_2s[88, c1] + ' 90 ' \
        + elems_array_2s[89, c1] + ' 91 ' + elems_array_2s[90, c1] + ' 92 ' \
        + elems_array_2s[91, c1] + '\n' \
        + ' ABUNDANCE CHANGE 93 ' + elems_array_2s[92, c1] + ' 94 ' \
        + elems_array_2s[93, c1] + ' 95 ' + elems_array_2s[94, c1] + ' 96 ' \
        + elems_array_2s[95, c1] + ' 97 ' + elems_array_2s[96, c1] + ' 98 ' \
        + elems_array_2s[97, c1] + '\n' \
        + ' ABUNDANCE CHANGE 99 ' + elems_array_2s[98, c1] + '\n' \
        + ' ABUNDANCE TABLE\n' \
        + '    1H   ' + renormed_H_5s \
        + '0       2He  ' + solar_He_5s + '0\n' \
        + '    3Li' + elems_array_3s[2, c1] \
        + ' 0.000    4Be' + elems_array_3s[3, c1] \
        + ' 0.000    5B ' + elems_array_3s[4, c1] \
        + ' 0.000    6C ' + elems_array_3s[5, c1] \
        + ' 0.000    7N ' + elems_array_3s[6, c1] + ' 0.000\n' \
        + '    8O ' + elems_array_3s[7, c1] \
        + ' 0.000    9F ' + elems_array_3s[8, c1] \
        + ' 0.000   10Ne' + elems_array_3s[9, c1] \
        + ' 0.000   11Na' + elems_array_3s[10, c1] \
        + ' 0.000   12Mg' + elems_array_3s[11, c1] + ' 0.000\n' \
        + '   13Al' + elems_array_3s[12, c1] \
        + ' 0.000   14Si' + elems_array_3s[13, c1] \
        + ' 0.000   15P ' + elems_array_3s[14, c1] \
        + ' 0.000   16S ' + elems_array_3s[15, c1] \
        + ' 0.000   17Cl' + elems_array_3s[16, c1] + ' 0.000\n' \
        + '   18Ar' + elems_array_3s[17, c1] \
        + ' 0.000   19K ' + elems_array_3s[18, c1] \
        + ' 0.000   20Ca' + elems_array_3s[19, c1] \
        + ' 0.000   21Sc' + elems_array_3s[20, c1] \
        + ' 0.000   22Ti' + elems_array_3s[21, c1] + ' 0.000\n' \
        + '   23V ' + elems_array_3s[22, c1] \
        + ' 0.000   24Cr' + elems_array_3s[23, c1] \
        + ' 0.000   25Mn' + elems_array_3s[24, c1] \
        + ' 0.000   26Fe' + elems_array_3s[25, c1] \
        + ' 0.000   27Co' + elems_array_3s[26, c1] + ' 0.000\n' \
        + '   28Ni' + elems_array_3s[27, c1] \
        + ' 0.000   29Cu' + elems_array_3s[28, c1] \
        + ' 0.000   30Zn' + elems_array_3s[29, c1] \
        + ' 0.000   31Ga' + elems_array_3s[30, c1] \
        + ' 0.000   32Ge' + elems_array_3s[31, c1] + ' 0.000\n' \
        + '   33As' + elems_array_3s[32, c1] \
        + ' 0.000   34Se' + elems_array_3s[33, c1] \
        + ' 0.000   35Br' + elems_array_3s[34, c1] \
        + ' 0.000   36Kr' + elems_array_3s[35, c1] \
        + ' 0.000   37Rb' + elems_array_3s[36, c1] + ' 0.000\n' \
        + '   38Sr' + elems_array_3s[37, c1] \
        + ' 0.000   39Y ' + elems_array_3s[38, c1] \
        + ' 0.000   40Zr' + elems_array_3s[39, c1] \
        + ' 0.000   41Nb' + elems_array_3s[40, c1] \
        + ' 0.000   42Mo' + elems_array_3s[41, c1] + ' 0.000\n' \
        + '   43Tc' + elems_array_3s[42, c1] \
        + ' 0.000   44Ru' + elems_array_3s[43, c1] \
        + ' 0.000   45Rh' + elems_array_3s[44, c1] \
        + ' 0.000   46Pd' + elems_array_3s[45, c1] \
        + ' 0.000   47Ag' + elems_array_3s[46, c1] + ' 0.000\n' \
        + '   48Cd' + elems_array_3s[47, c1] \
        + ' 0.000   49In' + elems_array_3s[48, c1] \
        + ' 0.000   50Sn' + elems_array_3s[49, c1] \
        + ' 0.000   51Sb' + elems_array_3s[50, c1] \
        + ' 0.000   52Te' + elems_array_3s[51, c1] + ' 0.000\n' \
        + '   53I ' + elems_array_3s[52, c1] \
        + ' 0.000   54Xe' + elems_array_3s[53, c1] \
        + ' 0.000   55Cs' + elems_array_3s[54, c1] \
        + ' 0.000   56Ba' + elems_array_3s[55, c1] \
        + ' 0.000   57La' + elems_array_3s[56, c1] + ' 0.000\n' \
        + '   58Ce' + elems_array_3s[57, c1] \
        + ' 0.000   59Pr' + elems_array_3s[58, c1] \
        + ' 0.000   60Nd' + elems_array_3s[59, c1] \
        + ' 0.000   61Pm' + elems_array_3s[60, c1] \
        + ' 0.000   62Sm' + elems_array_3s[61, c1] + ' 0.000\n' \
        + '   63Eu' + elems_array_3s[62, c1] \
        + ' 0.000   64Gd' + elems_array_3s[63, c1] \
        + ' 0.000   65Tb' + elems_array_3s[64, c1] \
        + ' 0.000   66Dy' + elems_array_3s[65, c1] \
        + ' 0.000   67Ho' + elems_array_3s[66, c1] + ' 0.000\n' \
        + '   68Er' + elems_array_3s[67, c1] \
        + ' 0.000   69Tm' + elems_array_3s[68, c1] \
        + ' 0.000   70Yb' + elems_array_3s[69, c1] \
        + ' 0.000   71Lu' + elems_array_3s[70, c1] \
        + ' 0.000   72Hf' + elems_array_3s[71, c1] + ' 0.000\n' \
        + '   73Ta' + elems_array_3s[72, c1] \
        + ' 0.000   74W ' + elems_array_3s[73, c1] \
        + ' 0.000   75Re' + elems_array_3s[74, c1] \
        + ' 0.000   76Os' + elems_array_3s[75, c1] \
        + ' 0.000   77Ir' + elems_array_3s[76, c1] + ' 0.000\n' \
        + '   78Pt' + elems_array_3s[77, c1] \
        + ' 0.000   79Au' + elems_array_3s[78, c1] \
        + ' 0.000   80Hg' + elems_array_3s[79, c1] \
        + ' 0.000   81Tl' + elems_array_3s[80, c1] \
        + ' 0.000   82Pb' + elems_array_3s[81, c1] + ' 0.000\n' \
        + '   83Bi' + elems_array_3s[82, c1] \
        + ' 0.000   84Po' + elems_array_3s[83, c1] \
        + ' 0.000   85At' + elems_array_3s[84, c1] \
        + ' 0.000   86Rn' + elems_array_3s[85, c1] \
        + ' 0.000   87Fr' + elems_array_3s[86, c1] + ' 0.000\n' \
        + '   88Ra' + elems_array_3s[87, c1] \
        + ' 0.000   89Ac' + elems_array_3s[88, c1] \
        + ' 0.000   90Th' + elems_array_3s[89, c1] \
        + ' 0.000   91Pa' + elems_array_3s[90, c1] \
        + ' 0.000   92U ' + elems_array_3s[91, c1] + ' 0.000\n' \
        + '   93NP' + elems_array_3s[92, c1] \
        + ' 0.000   94Pu' + elems_array_3s[93, c1] \
        + ' 0.000   95Am' + elems_array_3s[94, c1] \
        + ' 0.000   96Cm' + elems_array_3s[95, c1] \
        + ' 0.000   97Bk' + elems_array_3s[96, c1] + ' 0.000\n' \
        + '   98Cf' + elems_array_3s[97, c1] \
        + ' 0.000   99Es' + elems_array_3s[98, c1] + ' 0.000\n'

    file_text = start_file + end_file

    '''
    Name Model
    '''
    naming_index = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"]
    if c1 < 10:
        file_index = 'aaaa' + naming_index[int(str(c1)[0])]

    if 100 > c1 >= 10:
        file_index = 'aaa' + naming_index[int(str(c1)[0])] \
            + naming_index[int(str(c1)[1])]

    if 1000 > c1 >= 100:
        file_index = 'aa' + naming_index[int(str(c1)[0])] \
            + naming_index[int(str(c1)[1])] + naming_index[int(str(c1)[2])]

    if 10000 > c1 >= 1000:
        file_index = 'a' + naming_index[int(str(c1)[0])] \
            + naming_index[int(str(c1)[1])] + naming_index[int(str(c1)[2])] \
            + naming_index[int(str(c1)[3])]

    if 100000 > c1 >= 10000:
        file_index = naming_index[int(str(c1)[0])] \
            + naming_index[int(str(c1)[1])] + naming_index[int(str(c1)[2])] \
            + naming_index[int(str(c1)[3])] + naming_index[int(str(c1)[4])]

    '''
    Check if this Model is Already Generated
    '''
    list_files = os.listdir(khome.joinpath('Sync_Spectra_All'))
    for i in range(len(list_files)):
        list_files[i] = list_files[i][5:10]
    if file_index in list_files and verbose > 1:
        print(f'{file_index} previously generated')
    if not (file_index in list_files):
        if verbose > 1:
            print(f'{file_index} not previously generated')
        # Make Directory
        model_dir = output_dir.joinpath(file_index)
        atm_dir = model_dir.joinpath('atm')
        atm_dir.mkdir(parents=True)
        # Write Text to File
        model_name = f'{file_index}_t0{element_array[0, c1]:.0f}g{element_array[1, c1]:.2f}.atm'
        model_file = atm_dir.joinpath(model_name)
        f = open(model_file, 'w')
        f.write(file_text)
        f.close()
        if verbose > 1:
            print(f"Model {model_name} has been written")
        # Save Labels
        label_file = model_dir.joinpath(f'{file_index}_labels.npz')
        np.savez(label_file,
                 teff=element_array[0, c1],
                 logg=element_array[1, c1],
                 v_micro=labels['v_micro'][c1],
                 elems_array=np.concatenate([[renormed_H[c1]], [solar_He],
                                             element_array[2:, c1]]))

    '''
    Clean up Temporary Files
    '''
    template.unlink()

if verbose:
    print('All done!')
