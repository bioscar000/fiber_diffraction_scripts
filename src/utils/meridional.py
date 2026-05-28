"""
Meridional peak fitting and parameter extraction

Author: Óscar Fernández Blanco

Description:
This module implements the fitting procedure used to extract structural
parameters from meridional scattering profiles.

The main peak is modeled using a Lorentzian function, and the fitting region
is optimized by scanning different window sizes to maximize the goodness of fit.

Notes:
- The fitting range is determined dynamically around the peak position.
- The procedure is tailored to the experimental datasets used in this study.
- Not intended as a general-purpose peak fitting tool.

Part of the fiber diffraction analysis pipeline developed for a PhD thesis.
"""

__author__ = "Óscar Fernández Blanco"
__license__ = "MIT"


import numpy as np
from scipy.optimize import curve_fit


def calculate_meridional_parameters(q, i):
    """
    Fit the meridional peak and extract structural parameters.

    Parameters
    ----------
    q : numpy.ndarray
        Scattering vector.
    i : numpy.ndarray
        Intensity profile.

    Returns
    -------
    dict
        Extracted parameters (peak position, fit quality, monomer length).
    """

    # Remove zero-padding regions
    i_start, i_end = np.nonzero(i)[0][0], np.nonzero(i)[0][-1]
    q, i = q[i_start:i_end + 1], i[i_start:i_end + 1]

    # Select region where the meridional peak is expected
    idx = np.where((q >= 5.5) & (q <= 6.7))[0]

    # Locate peak center (maximum intensity)
    center = idx[np.argmax(i[idx])]

    def lorentzian(x, y0, a, x0, b):
        """Lorentzian peak model."""
        return y0 + a / (1 + ((x - x0) / b) ** 2)

    best_r2 = -np.inf
    best_popt = None

    # Scan different fitting windows around the peak
    for left in range(2, 15):
        for right in range(3, 15):

            if left + right > 10:
                start = max(center - left, 0)
                end = min(center + right, len(q))

                q_fit, i_fit = q[start:end], i[start:end]

                try:
                    popt = curve_fit(
                        lorentzian,
                        q_fit,
                        i_fit,
                        p0=[
                            np.min(i_fit),
                            np.max(i_fit),
                            np.mean(q_fit),
                            (q_fit[-1] - q_fit[0]) / 4
                        ],
                        method='lm'
                    )[0]

                except RuntimeError:
                    continue

                # Compute R² as goodness of fit metric
                residuals = i_fit - lorentzian(q_fit, *popt)
                ss_res = np.sum(residuals ** 2)
                ss_tot = np.sum((i_fit - np.mean(i_fit)) ** 2)

                r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

                if r2 > best_r2:
                    best_r2 = r2
                    best_popt = popt

    # Handle fitting failure
    if best_popt is None:
        return {
            'q_1nm': np.nan,
            'r2': np.nan,
            'avg_mon_len': np.nan
        }

    # Extract parameters
    return {
        'q_1nm': best_popt[2],
        'r2': best_r2,
        'avg_mon_len': 4 * 2 * np.pi / best_popt[2]
    }