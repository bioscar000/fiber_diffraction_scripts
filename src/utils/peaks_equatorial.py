"""
Equatorial peak detection and parameter extraction

Author: Óscar Fernández Blanco

Description:
This module provides utilities for detecting and refining peak positions
in equatorial scattering profiles, and computing derived structural parameters.

Peak positions are identified using scipy.signal.find_peaks and refined
via parabolic interpolation.

Notes:
- Peak detection parameters are tuned for the datasets used in this study.
- Includes optional plotting for manual inspection of peak detection.
- Not intended as a general-purpose peak detection library.

Part of the fiber diffraction analysis pipeline developed for a PhD thesis.
"""

__author__ = "Óscar Fernández Blanco"
__license__ = "MIT"


import numpy as np
from scipy.signal import find_peaks
from scipy.special import jn_zeros


def parabolic_refinement(q_roi, i_roi, peak_idx):
    """
    Refine peak position using parabolic interpolation.

    Parameters
    ----------
    q_roi : numpy.ndarray
        Scattering vector region.
    i_roi : numpy.ndarray
        Intensity values.
    peak_idx : int
        Index of the peak.

    Returns
    -------
    float
        Refined peak position.
    """

    x0, x1, x2 = q_roi[peak_idx - 1:peak_idx + 2]
    y0, y1, y2 = i_roi[peak_idx - 1:peak_idx + 2]

    h = x1 - x0
    xp = x1 + (h / 2) * (y0 - y2) / (y0 - 2 * y1 + y2)

    return np.float64(xp)


def qmin_to_distance(q, q_min_sep=0.25):
    """
    Convert minimum q separation into index distance.

    Parameters
    ----------
    q : numpy.ndarray
        Scattering vector.
    q_min_sep : float
        Minimum separation in q units.

    Returns
    -------
    int
        Minimum peak distance in index units.
    """

    dq = np.median(np.diff(q))
    return max(1, int(np.ceil(q_min_sep / dq)))


def calculate_peaks_parameters(q, i, title=None, save_plot=False, plot_path=None):
    """
    Detect and refine equatorial peaks, and compute structural parameters.

    Parameters
    ----------
    q : numpy.ndarray
        Scattering vector.
    i : numpy.ndarray
        Intensity profile.
    title : str, optional
        Label for plotting.
    save_plot : bool, default False
        Whether to save diagnostic plot.
    plot_path : str or Path, optional
        Output directory for plot.

    Returns
    -------
    dict
        Peak positions and derived radius values.
    """

    # Remove zero-padding regions
    i_start, i_end = np.nonzero(i)[0][0], np.nonzero(i)[0][-1]
    q, i = q[i_start:i_end + 1], i[i_start:i_end + 1]

    # Peak detection parameters
    peaks_distance = qmin_to_distance(q, q_min_sep=0.18)
    q0 = 0.2
    start = np.searchsorted(q, q0)

    # Detect peaks in region of interest
    j_peaks = find_peaks(i[start:], distance=peaks_distance)[0] + start

    # Estimate spacing between peaks
    diff3 = j_peaks[2] - j_peaks[1]

    # -------- Optional plotting (for debugging/manual validation) --------
    if save_plot and plot_path is not None:

        import matplotlib.pyplot as plt

        plt.plot(q, i, c='k')
        plt.scatter(q[j_peaks], i[j_peaks], c='r')

        plt.axvline(x=q[j_peaks[1]], color='r')
        plt.text(
            q[j_peaks[1]] + 0.1,
            i[j_peaks[1]] + 0.1 * i[j_peaks[1]],
            s=str(q[j_peaks[1]]),
            fontsize=10
        )

        plt.title(title or "Peak detection")

        plot_path.mkdir(parents=True, exist_ok=True)
        plt.savefig(plot_path / f"{title}.png")

        plt.close()

    # -------- Peak refinement --------

    try:
        j01_value = parabolic_refinement(q, i, j_peaks[0])
    except Exception:
        j01_value = np.nan

    try:
        j02_value = parabolic_refinement(q, i, j_peaks[1])
    except Exception:
        j02_value = np.nan

    try:
        j03_value = parabolic_refinement(q, i, j_peaks[2])
    except Exception:
        j03_value = np.nan

    try:
        j04_value = parabolic_refinement(q, i, j_peaks[2] + diff3)
    except Exception:
        j04_value = np.nan

    # -------- Derived parameters --------

    zeros = jn_zeros(1, 4)

    return {
        'j01': j01_value,
        'j02': j02_value,
        'j03': j03_value,
        'j04': j04_value,
        'r_j01': np.divide(zeros[0], j01_value),
        'r_j02': np.divide(zeros[1], j02_value),
        'r_j03': np.divide(zeros[2], j03_value),
        'r_j04': np.divide(zeros[3], j04_value)
    }