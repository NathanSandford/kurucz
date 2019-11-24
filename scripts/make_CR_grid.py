import argparse
import os
import sys
from pathlib import Path
import shutil
import numpy as np
from mendeleev import element

description = \
    '''
    Code generates initial atmospheres around a set of stellar labels
    '''

'''
Defaults
'''
kbase = Path('/global/scratch/nathan_sandford/kurucz/')
khome = kbase.joinpath('kurucz_home/')
krun_base = kbase.joinpath('kurucz_run/')
name = 'grids'
alphafe = 0.0
alpha_el = ['O', 'Ne', 'Mg', 'Si', 'S', 'Ar', 'Ca', 'Ti']
include_alpha = True
dX_teff = 50
dX_logg = 0.1
dX_v_micro = 0.1
dX_abundances = 0.05
symmetric = True
verbose = 1

parser = argparse.ArgumentParser(description=description)
parser.add_argument("krun", help=f"kurucz sub-directory where code runs")
parser.add_argument("teff", metavar="T_eff", type=float, help="Effective Temperature (K)")
parser.add_argument("logg", metavar="log(g)", type=float, help="Surface Gravity")
parser.add_argument("feh", metavar="[Fe/H]", type=float, help="Solar Scaled Iron Abundance")
parser.add_argument("v_micro", type=float, help="Micro-turbulent Velocity (km/s)")
parser.add_argument("--kbase", "-kb",
                    help=f"Base directory (default: {kbase})")
parser.add_argument("--khome", "-kh",
                    help=f"Path to Kurucz Code (default: {khome.name})")
parser.add_argument("--krun_base", "-krb",
                    help=f"Path to Directory where code runs {krun_base.name}")
parser.add_argument("--name", "-n",
                    help=f"Directory in krun to save atmospheres (default: {name})")
parser.add_argument("--alphafe", "-afe", type=float,
                    help=f"[alpha/Fe] (default: {alphafe})")
parser.add_argument("--alpha_el", nargs='+',
                    help=f"Free abundances (default: {' '.join(alpha_el)})")
parser.add_argument("--ignore_alpha", action="store_true",
                    help="Ignore [alpha/Fe] perturbations")
parser.add_argument("--dX_teff", "-dT", type=float,
                    help=f"Step size in T_eff (default: {dX_teff} K)")
parser.add_argument("--dX_logg", "-dg", type=float,
                    help=f"Step size in log(g) (default: {dX_logg})")
parser.add_argument("--dX_v_micro", "-dv", type=float,
                    help=f"Step size in v_micro (default: {dX_v_micro} km/s)")
parser.add_argument("--dX_abundances", "-dX", type=float,
                    help=f"Step size in [X/H] (default: {dX_abundances} dex)")
parser.add_argument("--asymmetric", "-A", action="store_true",
                    help="Generate atmospheres asymmetrically about reference point")
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
        shutil.rmtree(output_dir)
        output_dir.mkdir()

if args.alphafe:
    alphafe = alphafe
if args.alpha_el:
    alpha_el = args.alpha_el
alpha_ind = [el.atomic_number for el in element(alpha_el)]
if args.ignore_alpha:
    include_alpha = False

if args.dX_teff:
    dX_teff = args.dX_teff
if args.dX_logg:
    dX_logg = args.dX_logg
if args.dX_v_micro:
    dX_v_micro = args.dX_v_micro
if args.dX_abundances:
    dX_abundances = args.dX_abundances
if args.asymmetric:
    symmetric = False

if args.verbose:
    verbose = args.verbose

teff = args.teff
logg = args.logg
feh = args.feh
v_micro = args.v_micro  # km/s


'''Summary'''
if verbose:
    print(f'Reference Atmosphere: Teff={teff:.0f} K, Log(g)={logg:.1f}, '
          + f'[Fe/H]={feh:.2f}, [alpha/Fe]={alphafe:.2f}, '
          + f'v_micro={v_micro:.2} km/s')

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

if symmetric:
    num_models = 2 * len(element_array) + 3  # 2/(el,logg,teff) + 2/vmicro + 1 ref
    if include_alpha:
        num_models += 2
else:
    num_models = len(element_array) + 2  # 1/(el,logg,teff) + 1/vmicro + 1 ref
    if include_alpha:
        num_models += 1
element_array = np.repeat(np.array([element_array]), num_models, axis=0).T
if verbose:
    print(f'Generating {num_models} models')

'''
Change teff & logg
'''
element_array[0, :] = teff
element_array[1, :] = logg
v_micro_array = v_micro * np.ones(num_models)

'''
Scale global metallicity
'''
element_array[2:, :] += feh
element_array[alpha_ind, :] += alphafe

'''
Apply offsets to all labels
'''
if symmetric:
    element_array[0, 1] += dX_teff
    element_array[0, 2] -= dX_teff
    element_array[1, 3] += dX_logg
    element_array[1, 4] -= dX_logg
    v_micro_array[5] += dX_v_micro
    v_micro_array[6] -= dX_v_micro
    if include_alpha:
        element_array[alpha_ind, 7] += dX_abundances
        element_array[alpha_ind, 8] -= dX_abundances
        a = 9
    else:
        a = 7
    diag_ind = np.diag_indices_from(element_array[2:, a::2])
    element_array[2:, a::2][diag_ind] += dX_abundances
    element_array[2:, (a+1)::2][diag_ind] -= dX_abundances
else:
    element_array[0, 1] += dX_teff
    element_array[1, 2] += dX_logg
    v_micro_array[3] += dX_v_micro
    if include_alpha:
        element_array[alpha_ind, 4] += dX_abundances
        element_array[alpha_ind, 5] -= dX_abundances
        a = 6
    else:
        a = 4
    diag_ind = np.diag_indices_from(element_array[2:, a::2])
    element_array[2:, a::2][diag_ind] += dX_abundances

'''
Renormalize Hydrogen such that X+Y+Z=1
'''
solar_He = 0.07837
solar_He_5s = f"{solar_He:.5f}"
renormed_H = 1. - solar_He - np.sum(10. ** element_array[2:, :], axis=0)


'''
Make Input Kurucz Readable
'''
element_array_0s = np.copy(element_array).astype("str")
element_array_2s = np.copy(element_array).astype("str")
element_array_3s = np.copy(element_array).astype("str")
element_array_4s = np.copy(element_array).astype("str")
for p1 in range(element_array.shape[0]):
    for p2 in range(element_array.shape[1]):
        element_array_0s[p1, p2] = '' + f"{element_array[p1,p2]:.0f}"
        element_array_2s[p1, p2] = ' ' + f"{element_array[p1, p2]:.2f}"
        element_array_3s[p1, p2] = ' ' + f"{element_array[p1, p2]:.3f}"
        element_array_4s[p1, p2] = '' + f"{element_array[p1, p2]:.4f}"

        if element_array[p1, p2] <= -9.995:
            element_array_2s[p1, p2] = element_array_2s[p1, p2][1:]

        if element_array[p1, p2] < -9.9995:
            element_array_3s[p1, p2] = element_array_3s[p1, p2][1:]


'''
Restore Grid of Converged Atmospheres
'''
teff_find = []
logg_find = []
atm_find = []
name_find = []

if verbose > 1:
    print('Restoring generated grid from Sync_Spectra_All_Atm')
Sync_Spectra_All_Atm = khome.joinpath('Sync_Spectra_All_Atm')
list_files_find = os.listdir(Sync_Spectra_All_Atm)
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
for c1 in range(num_models):
    if verbose > 1:
        print('Modifying closest atmosphere file: {}'.format(atm_closest[c1]))
    lines = open(Sync_Spectra_All_Atm.joinpath(atm_closest[c1])).readlines()
    template = khome.joinpath('template_1.atm')
    open(template, 'w').writelines(lines[44:])

    f = open(template, 'r')
    end_file = f.read()
    f.close()

    renormed_H_5s = f"{renormed_H[c1]:.5f}"

    '''
    Change Abundances
    '''
    start_file = 'TEFF   ' + element_array_0s[0, c1] \
                 + '.  GRAVITY  ' + element_array_4s[1, c1] + ' LTE \n' \
                 + 'TITLE ATLAS12    ' \
                 + '                                                               \n' \
                 + ' OPACITY IFOP 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 1 0 0 0\n' \
                 + ' CONVECTION ON   1.25 TURBULENCE OFF  0.00  0.00  0.00  0.00\n' \
                 + 'ABUNDANCE SCALE   1.00000 ABUNDANCE CHANGE 1 ' \
                 + renormed_H_5s + ' 2 ' + solar_He_5s + '\n' \
                 + ' ABUNDANCE CHANGE  3 ' + element_array_2s[2, c1] + '  4 ' \
                 + element_array_2s[3, c1] + '  5 ' + element_array_2s[4, c1] + '  6 ' \
                 + element_array_2s[5, c1] + '  7 ' + element_array_2s[6, c1] + '  8 ' \
                 + element_array_2s[7, c1] + '\n' \
                 + ' ABUNDANCE CHANGE  9 ' + element_array_2s[8, c1] + ' 10 ' \
                 + element_array_2s[9, c1] + ' 11 ' + element_array_2s[10, c1] + ' 12 ' \
                 + element_array_2s[11, c1] + ' 13 ' + element_array_2s[12, c1] + ' 14 ' \
                 + element_array_2s[13, c1] + '\n' \
                 + ' ABUNDANCE CHANGE 15 ' + element_array_2s[14, c1] + ' 16 ' \
                 + element_array_2s[15, c1] + ' 17 ' + element_array_2s[16, c1] + ' 18 ' \
                 + element_array_2s[17, c1] + ' 19 ' + element_array_2s[18, c1] + ' 20 ' \
                 + element_array_2s[19, c1] + '\n' \
                 + ' ABUNDANCE CHANGE 21 ' + element_array_2s[20, c1] + ' 22 ' \
                 + element_array_2s[21, c1] + ' 23 ' + element_array_2s[22, c1] + ' 24 ' \
                 + element_array_2s[23, c1] + ' 25 ' + element_array_2s[24, c1] + ' 26 ' \
                 + element_array_2s[25, c1] + '\n' \
                 + ' ABUNDANCE CHANGE 27 ' + element_array_2s[26, c1] + ' 28 ' \
                 + element_array_2s[27, c1] + ' 29 ' + element_array_2s[28, c1] + ' 30 ' \
                 + element_array_2s[29, c1] + ' 31 ' + element_array_2s[30, c1] + ' 32 ' \
                 + element_array_2s[31, c1] + '\n' \
                 + ' ABUNDANCE CHANGE 33 ' + element_array_2s[32, c1] + ' 34 ' \
                 + element_array_2s[33, c1] + ' 35 ' + element_array_2s[34, c1] + ' 36 ' \
                 + element_array_2s[35, c1] + ' 37 ' + element_array_2s[36, c1] + ' 38 ' \
                 + element_array_2s[37, c1] + '\n' \
                 + ' ABUNDANCE CHANGE 39 ' + element_array_2s[38, c1] + ' 40 ' \
                 + element_array_2s[39, c1] + ' 41 ' + element_array_2s[40, c1] + ' 42 ' \
                 + element_array_2s[41, c1] + ' 43 ' + element_array_2s[42, c1] + ' 44 ' \
                 + element_array_2s[43, c1] + '\n' \
                 + ' ABUNDANCE CHANGE 45 ' + element_array_2s[44, c1] + ' 46 ' \
                 + element_array_2s[45, c1] + ' 47 ' + element_array_2s[46, c1] + ' 48 ' \
                 + element_array_2s[47, c1] + ' 49 ' + element_array_2s[48, c1] + ' 50 ' \
                 + element_array_2s[49, c1] + '\n' \
                 + ' ABUNDANCE CHANGE 51 ' + element_array_2s[50, c1] + ' 52 ' \
                 + element_array_2s[51, c1] + ' 53 ' + element_array_2s[52, c1] + ' 54 ' \
                 + element_array_2s[53, c1] + ' 55 ' + element_array_2s[54, c1] + ' 56 ' \
                 + element_array_2s[55, c1] + '\n' \
                 + ' ABUNDANCE CHANGE 57 ' + element_array_2s[56, c1] + ' 58 ' \
                 + element_array_2s[57, c1] + ' 59 ' + element_array_2s[58, c1] + ' 60 ' \
                 + element_array_2s[59, c1] + ' 61 ' + element_array_2s[60, c1] + ' 62 ' \
                 + element_array_2s[61, c1] + '\n' \
                 + ' ABUNDANCE CHANGE 63 ' + element_array_2s[62, c1] + ' 64 ' \
                 + element_array_2s[63, c1] + ' 65 ' + element_array_2s[64, c1] + ' 66 ' \
                 + element_array_2s[65, c1] + ' 67 ' + element_array_2s[66, c1] + ' 68 ' \
                 + element_array_2s[67, c1] + '\n' \
                 + ' ABUNDANCE CHANGE 69 ' + element_array_2s[68, c1] + ' 70 ' \
                 + element_array_2s[69, c1] + ' 71 ' + element_array_2s[70, c1] + ' 72 ' \
                 + element_array_2s[71, c1] + ' 73 ' + element_array_2s[72, c1] + ' 74 ' \
                 + element_array_2s[73, c1] + '\n' \
                 + ' ABUNDANCE CHANGE 75 ' + element_array_2s[74, c1] + ' 76 ' \
                 + element_array_2s[75, c1] + ' 77 ' + element_array_2s[76, c1] + ' 78 ' \
                 + element_array_2s[77, c1] + ' 79 ' + element_array_2s[78, c1] + ' 80 ' \
                 + element_array_2s[79, c1] + '\n' \
                 + ' ABUNDANCE CHANGE 81 ' + element_array_2s[80, c1] + ' 82 ' \
                 + element_array_2s[81, c1] + ' 83 ' + element_array_2s[82, c1] + ' 84 ' \
                 + element_array_2s[83, c1] + ' 85 ' + element_array_2s[84, c1] + ' 86 ' \
                 + element_array_2s[85, c1] + '\n' \
                 + ' ABUNDANCE CHANGE 87 ' + element_array_2s[86, c1] + ' 88 ' \
                 + element_array_2s[87, c1] + ' 89 ' + element_array_2s[88, c1] + ' 90 ' \
                 + element_array_2s[89, c1] + ' 91 ' + element_array_2s[90, c1] + ' 92 ' \
                 + element_array_2s[91, c1] + '\n' \
                 + ' ABUNDANCE CHANGE 93 ' + element_array_2s[92, c1] + ' 94 ' \
                 + element_array_2s[93, c1] + ' 95 ' + element_array_2s[94, c1] + ' 96 ' \
                 + element_array_2s[95, c1] + ' 97 ' + element_array_2s[96, c1] + ' 98 ' \
                 + element_array_2s[97, c1] + '\n' \
                 + ' ABUNDANCE CHANGE 99 ' + element_array_2s[98, c1] + '\n' \
                 + ' ABUNDANCE TABLE\n' \
                 + '    1H   ' + renormed_H_5s \
                 + '0       2He  ' + solar_He_5s + '0\n' \
                 + '    3Li' + element_array_3s[2, c1] \
                 + ' 0.000    4Be' + element_array_3s[3, c1] \
                 + ' 0.000    5B ' + element_array_3s[4, c1] \
                 + ' 0.000    6C ' + element_array_3s[5, c1] \
                 + ' 0.000    7N ' + element_array_3s[6, c1] + ' 0.000\n' \
                 + '    8O ' + element_array_3s[7, c1] \
                 + ' 0.000    9F ' + element_array_3s[8, c1] \
                 + ' 0.000   10Ne' + element_array_3s[9, c1] \
                 + ' 0.000   11Na' + element_array_3s[10, c1] \
                 + ' 0.000   12Mg' + element_array_3s[11, c1] + ' 0.000\n' \
                 + '   13Al' + element_array_3s[12, c1] \
                 + ' 0.000   14Si' + element_array_3s[13, c1] \
                 + ' 0.000   15P ' + element_array_3s[14, c1] \
                 + ' 0.000   16S ' + element_array_3s[15, c1] \
                 + ' 0.000   17Cl' + element_array_3s[16, c1] + ' 0.000\n' \
                 + '   18Ar' + element_array_3s[17, c1] \
                 + ' 0.000   19K ' + element_array_3s[18, c1] \
                 + ' 0.000   20Ca' + element_array_3s[19, c1] \
                 + ' 0.000   21Sc' + element_array_3s[20, c1] \
                 + ' 0.000   22Ti' + element_array_3s[21, c1] + ' 0.000\n' \
                 + '   23V ' + element_array_3s[22, c1] \
                 + ' 0.000   24Cr' + element_array_3s[23, c1] \
                 + ' 0.000   25Mn' + element_array_3s[24, c1] \
                 + ' 0.000   26Fe' + element_array_3s[25, c1] \
                 + ' 0.000   27Co' + element_array_3s[26, c1] + ' 0.000\n' \
                 + '   28Ni' + element_array_3s[27, c1] \
                 + ' 0.000   29Cu' + element_array_3s[28, c1] \
                 + ' 0.000   30Zn' + element_array_3s[29, c1] \
                 + ' 0.000   31Ga' + element_array_3s[30, c1] \
                 + ' 0.000   32Ge' + element_array_3s[31, c1] + ' 0.000\n' \
                 + '   33As' + element_array_3s[32, c1] \
                 + ' 0.000   34Se' + element_array_3s[33, c1] \
                 + ' 0.000   35Br' + element_array_3s[34, c1] \
                 + ' 0.000   36Kr' + element_array_3s[35, c1] \
                 + ' 0.000   37Rb' + element_array_3s[36, c1] + ' 0.000\n' \
                 + '   38Sr' + element_array_3s[37, c1] \
                 + ' 0.000   39Y ' + element_array_3s[38, c1] \
                 + ' 0.000   40Zr' + element_array_3s[39, c1] \
                 + ' 0.000   41Nb' + element_array_3s[40, c1] \
                 + ' 0.000   42Mo' + element_array_3s[41, c1] + ' 0.000\n' \
                 + '   43Tc' + element_array_3s[42, c1] \
                 + ' 0.000   44Ru' + element_array_3s[43, c1] \
                 + ' 0.000   45Rh' + element_array_3s[44, c1] \
                 + ' 0.000   46Pd' + element_array_3s[45, c1] \
                 + ' 0.000   47Ag' + element_array_3s[46, c1] + ' 0.000\n' \
                 + '   48Cd' + element_array_3s[47, c1] \
                 + ' 0.000   49In' + element_array_3s[48, c1] \
                 + ' 0.000   50Sn' + element_array_3s[49, c1] \
                 + ' 0.000   51Sb' + element_array_3s[50, c1] \
                 + ' 0.000   52Te' + element_array_3s[51, c1] + ' 0.000\n' \
                 + '   53I ' + element_array_3s[52, c1] \
                 + ' 0.000   54Xe' + element_array_3s[53, c1] \
                 + ' 0.000   55Cs' + element_array_3s[54, c1] \
                 + ' 0.000   56Ba' + element_array_3s[55, c1] \
                 + ' 0.000   57La' + element_array_3s[56, c1] + ' 0.000\n' \
                 + '   58Ce' + element_array_3s[57, c1] \
                 + ' 0.000   59Pr' + element_array_3s[58, c1] \
                 + ' 0.000   60Nd' + element_array_3s[59, c1] \
                 + ' 0.000   61Pm' + element_array_3s[60, c1] \
                 + ' 0.000   62Sm' + element_array_3s[61, c1] + ' 0.000\n' \
                 + '   63Eu' + element_array_3s[62, c1] \
                 + ' 0.000   64Gd' + element_array_3s[63, c1] \
                 + ' 0.000   65Tb' + element_array_3s[64, c1] \
                 + ' 0.000   66Dy' + element_array_3s[65, c1] \
                 + ' 0.000   67Ho' + element_array_3s[66, c1] + ' 0.000\n' \
                 + '   68Er' + element_array_3s[67, c1] \
                 + ' 0.000   69Tm' + element_array_3s[68, c1] \
                 + ' 0.000   70Yb' + element_array_3s[69, c1] \
                 + ' 0.000   71Lu' + element_array_3s[70, c1] \
                 + ' 0.000   72Hf' + element_array_3s[71, c1] + ' 0.000\n' \
                 + '   73Ta' + element_array_3s[72, c1] \
                 + ' 0.000   74W ' + element_array_3s[73, c1] \
                 + ' 0.000   75Re' + element_array_3s[74, c1] \
                 + ' 0.000   76Os' + element_array_3s[75, c1] \
                 + ' 0.000   77Ir' + element_array_3s[76, c1] + ' 0.000\n' \
                 + '   78Pt' + element_array_3s[77, c1] \
                 + ' 0.000   79Au' + element_array_3s[78, c1] \
                 + ' 0.000   80Hg' + element_array_3s[79, c1] \
                 + ' 0.000   81Tl' + element_array_3s[80, c1] \
                 + ' 0.000   82Pb' + element_array_3s[81, c1] + ' 0.000\n' \
                 + '   83Bi' + element_array_3s[82, c1] \
                 + ' 0.000   84Po' + element_array_3s[83, c1] \
                 + ' 0.000   85At' + element_array_3s[84, c1] \
                 + ' 0.000   86Rn' + element_array_3s[85, c1] \
                 + ' 0.000   87Fr' + element_array_3s[86, c1] + ' 0.000\n' \
                 + '   88Ra' + element_array_3s[87, c1] \
                 + ' 0.000   89Ac' + element_array_3s[88, c1] \
                 + ' 0.000   90Th' + element_array_3s[89, c1] \
                 + ' 0.000   91Pa' + element_array_3s[90, c1] \
                 + ' 0.000   92U ' + element_array_3s[91, c1] + ' 0.000\n' \
                 + '   93NP' + element_array_3s[92, c1] \
                 + ' 0.000   94Pu' + element_array_3s[93, c1] \
                 + ' 0.000   95Am' + element_array_3s[94, c1] \
                 + ' 0.000   96Cm' + element_array_3s[95, c1] \
                 + ' 0.000   97Bk' + element_array_3s[96, c1] + ' 0.000\n' \
                 + '   98Cf' + element_array_3s[97, c1] \
                 + ' 0.000   99Es' + element_array_3s[98, c1] + ' 0.000\n'

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
    # noinspection PyUnboundLocalVariable
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
        model_name = f'{file_index}_t{element_array[0, c1]:05.0f}g{element_array[1, c1]:.2f}.atm'
        model_file = atm_dir.joinpath(model_name)
        f = open(model_file, 'w')
        f.write(file_text)
        f.close()
        if verbose > 1:
            print(f"Model {model_name} has been written")
        # Save Labels
        label_file = model_dir.joinpath(f'{file_index}_labels.npz')
        np.savez(label_file,
                 alpha_included=include_alpha,
                 teff=element_array[0, c1],
                 logg=element_array[1, c1],
                 v_micro=v_micro_array[c1],
                 elems_array=np.concatenate([[renormed_H[c1]], [solar_He],
                                             element_array[2:, c1]]))

    '''
    Clean up Temporary Files
    '''
    template.unlink()

if verbose:
    print('All done!')
