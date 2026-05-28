"""
Equatorial scattering model and parameter extraction

Author: Óscar Fernández Blanco

Description:
This module implements the model and fitting procedure used to extract
structural parameters from equatorial scattering profiles.

The model is based on a weighted mixture of protofilament populations,
using Bessel functions to describe diffraction intensities.

Notes:
- The implementation is tailored to the experimental setup used in this study.
- Parameter bounds and initialization are empirically chosen.
- Not intended as a general-purpose fitting library.

Part of the fiber diffraction analysis pipeline developed for a PhD thesis.
"""

__author__ = "Óscar Fernández Blanco"
__license__ = "MIT"


import numpy as np
from scipy.optimize import curve_fit
from scipy.special import jn


def calculate_equatorial_parameters(q, i):
    """
    Fit equatorial scattering profile and extract structural parameters.

    Parameters
    ----------
    q : numpy.ndarray
        Scattering vector.
    i : numpy.ndarray
        Intensity profile.

    Returns
    -------
    dict
        Extracted structural parameters.
    """

    # Remove zero-padding regions
    i_start, i_end = np.nonzero(i)[0][0], np.nonzero(i)[0][-1]
    q, i = q[i_start:i_end + 1], i[i_start:i_end + 1]

    # Select fitting window
    idx = np.where((q >= 0.75) & (q <= 1.25))[0]
    box = slice(idx[0], idx[-1] + 1)
    q_fit, i_fit = q[box], i[box]

    # Initial parameters and bounds
    p0 = [5.65, 1., 1., .2, .2, .2, .2, .2]
    bounds = (
        [5.6, 0., 0., 0., 0., 0., 0., 0.],
        [5.7, np.inf, np.inf, 1., 1., 1., 1., 1.]
    )

    # Non-linear least squares fitting
    popt = curve_fit(
        mix_population_intensity,
        q_fit,
        i_fit,
        p0=p0,
        bounds=bounds,
        check_finite=True,
        method='trf',
        max_nfev=1_000_000
    )[0]

    return parse_fit_results(popt)


def parse_fit_results(popt):
    """
    Convert fitted parameters into physically meaningful quantities.

    Parameters
    ----------
    popt : array-like
        Optimized parameters from curve fitting.

    Returns
    -------
    dict
        Structural parameters (radius, average protofilaments, weights, etc.).
    """

    # popt format: [ipfd, s0, sn, w11..w15]
    ipfd = popt[0]

    n_weights = normalize_weights(popt[3:9])
    n_list = np.arange(11, 16, dtype=np.int32)

    avg_n = np.dot(n_list, n_weights)

    # Radius estimation
    r = ipfd / (2.0 * np.sin(np.pi / avg_n))

    return {
        'r': r,
        'avg_n': avg_n,
        'ipfd': ipfd,
        **{f'%{ni}': wi for ni, wi in zip(n_list, n_weights)},
        's0': popt[1],
        'sn': popt[2],
    }


def mix_population_intensity(q, ipfd, s0, sn, *weights):
    """
    Model for mixed protofilament populations.
    """

    n_weights = normalize_weights(weights)
    n_list = np.arange(11, 16, dtype=np.int32)

    avg_n = np.dot(n_list, n_weights)

    den = 2.0 * np.sin(np.pi / avg_n)
    if not np.isfinite(den) or np.abs(den) < 1e-12:
        den = 1e-12

    r = ipfd / den

    # Recalculate ipfd from model constraint
    ipfd = 2 * np.sin((2 * np.pi / avg_n) / 2) * r

    # Enforce physical bounds
    if not (5.6 <= ipfd <= 5.7):
        return np.inf

    return np.sum([
        w * single_population_intensity(q, s0, sn, r / avg_n * n, n)
        for w, n in zip(n_weights, n_list)
    ], axis=0)


def single_population_intensity(q, s0, sn, r, n):
    """
    Intensity contribution of a single protofilament population.
    """

    j_0 = jn(0, r * q)
    j_n = jn(n, r * q)

    return s0 * ((fu(q) * j_0) ** 2 + sn * ((fu(q) * j_n) ** 2))


def fu(q, rt=2.48):
    """
    Form factor of a sphere (used as approximation).
    """

    u = q * rt
    v = (4 / 3) * np.pi * rt ** 3

    with np.errstate(divide='ignore', invalid='ignore'):
        f_u = v * (3 * (np.sin(u) - u * np.cos(u)) / u ** 3)
        f_u = np.where(np.abs(u) < 1e-10, v, f_u)

    return f_u


def normalize_weights(weights):
    """
    Normalize population weights to sum to 1.
    """

    weights = np.array(weights, dtype=np.float64)
    return weights / weights.sum()