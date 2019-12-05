#!/usr/bin/env python
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from numpy.fft import rfftfreq
from scipy.interpolate import interp1d

description = \
    '''
    Convolves batch of high-res spectra down to desired resolution and wavelength range.
    '''

'''
Defaults
'''
kbase = Path('/global/scratch/nathan_sandford/kurucz')
kout = 'kurucz_out'
spec_dir = 'synthetic_spectra'
output_dir = 'convolved'
output_tag = '_conv'
R_Samp = 3

'''
Parse Args
'''
parser = argparse.ArgumentParser(description=description)
parser.add_argument("spec_file", metavar="spectra", help="File containing high-res spectra")
parser.add_argument("R_Res", metavar="R", type=float, help="Resolving Power")
parser.add_argument("start_wavelength", metavar="Start", type=float,
                    help="Lower wavelength bound")
parser.add_argument("end_wavelength", metavar="End", type=float,
                    help="Upper wavelength bound")
parser.add_argument("--kbase", "-kb",
                    help=f"Base directory (default: {kbase})")
parser.add_argument("--kout", "-ko",
                    help=f"Relative path to Output Directory from kbase (default: {kout})")
parser.add_argument("--spec_dir", "-sd",
                    help=f"Relative path to Spectra Directory from kout (default: {spec_dir})")
parser.add_argument("--output_dir", "-od",
                    help=f"Relative Path to Output Directory from spec_dir (default: {output_dir})")
parser.add_argument("--abs_output_dir", "-aod",
                    help=f"Absolute Path to Output Directory - overrides output_dir (default: N/A)")
parser.add_argument("--output_tag", "-o",
                    help="tag to append to convolved spectra files (default: {output_tag})")
parser.add_argument("--normalize", '-N', action="store_true", default=False,
                    help=f"Return Normalized Spectra (default: False)")
parser.add_argument("--R_Samp", "-rs", type=float,
                    help=f"Wavelength Sampling in pixels/FWHM (default: {R_Samp})")
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

if args.abs_output_dir:
    output_dir = Path(args.abs_output_dir)
elif args.output_dir:
    output_dir = spec_dir.joinpath(args.output_dir)
else:
    output_dir = spec_dir.joinpath(output_dir)
if not output_dir.is_dir():
    output_dir.mkdir(parents=True, exist_ok=True)
assert output_dir.is_dir(), f'Output directory {output_dir} does not exist'

if args.spec_file[-3:] != '.h5':
    spec_file = spec_dir.joinpath(args.spec_file+'.h5')
else:
    spec_file = spec_dir.joinpath(args.spec_file)

if args.output_tag:
    output_file = output_dir.joinpath(spec_file.name[:-3]+f'{args.output_tag}.h5')
else:
    output_file = output_dir.joinpath(spec_file.name[:-3]+f'{output_tag}.h5')

R_Res = args.R_Res
start_wavelength = args.start_wavelength
end_wavelength = args.end_wavelength
normalize = args.normalize

if args.R_Samp:
    R_Samp = args.R_Samp


def generate_wavelength_template(start_wavelength, end_wavelength,
                                 resolution, truncate=False):
    """

    :param start_wavelength:
    :param end_wavelength:
    :param resolution:
    :param truncate:
    :return:
    """
    '''
    Generate wavelength array with fixed resolution.

    args:
    start_wavelength = minimum wavelength of spectra to include
    end_wavelength = maximum wavelength of spectra to include
    resolution = resolving power of instrument (R = lambda / delta lambda)
    truncate = Boolean, if true, drop final pixel for which lambda > end_wavelength

    returns = wavelength grid of given resolution between start and end wavelengths
    '''
    # TODO: generate_wavelength_template doc-string
    wavelength_template = [start_wavelength]
    wavelength_now = start_wavelength

    while wavelength_now < end_wavelength:
        wavelength_now += wavelength_now / resolution
        wavelength_template.append(wavelength_now)
    wavelength_template = np.array(wavelength_template)

    if truncate:
        wavelength_template = wavelength_template[:-1]
    return wavelength_template


def convolve_spec(wave, spec, resolution, outwave, res_in=None):
    """

    :param wave:
    :param spec:
    :param resolution:
    :param outwave:
    :param res_in:
    :return:
    """
    # TODO: convolve_spec doc-string
    sigma_to_fwhm = 2.355

    width = resolution * sigma_to_fwhm
    sigma_out = (resolution * sigma_to_fwhm) ** -1
    if res_in is None:
        sigma_in = 0
    else:
        sigma_in = (res_in * sigma_to_fwhm) ** -1

    # Trim Wavelength Range
    nsigma_pad = 20.0
    wlim = np.array([outwave.min(), outwave.max()])
    wlim *= (1 + nsigma_pad / width * np.array([-1, 1]))
    mask = (wave > wlim[0]) & (wave < wlim[1])
    wave = wave[mask]
    if spec.ndim == 1:
        spec = spec[mask]
    else:
        spec = spec[:, mask]

    # Make Convolution Grid
    wmin, wmax = wave.min(), wave.max()
    nwave = wave.shape[0]
    nwave_new = int(2 ** (np.ceil(np.log2(nwave))))
    lnwave_new = np.linspace(np.log(wmin), np.log(wmax), nwave_new)
    wave_new = np.exp(lnwave_new)
    fspec = interp1d(wave, spec,
                     bounds_error=False, fill_value='extrapolate')
    spec = fspec(wave_new)
    wave = wave_new

    # Convolve via FFT
    sigma = np.sqrt(sigma_out ** 2 - sigma_in ** 2)
    invres_grid = np.diff(np.log(wave))
    dx = np.median(invres_grid)
    ss = rfftfreq(nwave_new, d=dx)
    taper = np.exp(-2 * (np.pi ** 2) * (sigma ** 2) * (ss ** 2))
    spec_ff = np.fft.rfft(spec)
    ff_tapered = spec_ff * taper
    spec_conv = np.fft.irfft(ff_tapered)

    # Interpolate onto outwave
    fspec = interp1d(wave, spec_conv,
                     bounds_error=False, fill_value='extrapolate')
    return fspec(outwave)


# Generate Wavelength Grid
wave_template \
        = generate_wavelength_template(
                             start_wavelength=start_wavelength,
                             end_wavelength=end_wavelength,
                             resolution=R_Res*R_Samp)
wave_template = pd.DataFrame(wave_template)

# Restore High-Res Spectra
wavelength = pd.read_hdf(spec_file, 'wavelength')
labels = pd.read_hdf(spec_file, 'labels')
spectra = pd.read_hdf(spec_file, 'spectra')
continuum = pd.read_hdf(spec_file, 'continuum')

if normalize:
    spectra /= continuum

print(f'Beginning convolution of {spec_file.name}')
conv_spec = convolve_spec(wave=wavelength[0].values,
                          spec=spectra.values.T,
                          resolution=R_Res,
                          outwave=wave_template[0].values,
                          res_in=300000)
conv_spec = pd.DataFrame(conv_spec.T, columns=spectra.columns)

wave_template.to_hdf(output_file, 'wavelength')
conv_spec.to_hdf(output_file, 'spectra')
labels.to_hdf(output_file, 'labels')
print(f"Saved convolved spectra to {output_file}")
