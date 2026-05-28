"""
Integration utilities for fiber diffraction analysis

Author: Óscar Fernández Blanco

Description:
This module provides functions for performing 1D azimuthal integration
of diffraction images using pyFAI, including equatorial and meridional cuts.

Notes:
- The azimuthal selection depends on specific experimental beamtime configurations.
- The current implementation is tailored to datasets generated in our lab.
- Users working with different experimental setups may need to adapt these parameters.

Part of the fiber diffraction analysis pipeline developed for a PhD thesis.
"""

__author__ = "Óscar Fernández Blanco"
__license__ = "MIT"


import os

import numpy as np
import pyFAI


def load_calibration_file(poni_file_path):
    """
    Load a pyFAI calibration (.poni) file.

    Parameters
    ----------
    poni_file_path : str or Path
        Path to the calibration file.

    Returns
    -------
    pyFAI.azimuthalIntegrator.AzimuthalIntegrator
        Loaded calibration object.

    Raises
    ------
    FileNotFoundError
        If the calibration file does not exist.
    """

    if os.path.exists(poni_file_path):
        return pyFAI.load(poni_file_path)
    else:
        raise FileNotFoundError(
            f"Calibration file not found at path: {poni_file_path}"
        )


def d1_int(data, cut, ai):
    """
    Perform 1D azimuthal integration over a given angular range.

    Parameters
    ----------
    data : numpy.ndarray
        Diffraction image.
    cut : tuple of float
        Azimuthal angle range (degrees).
    ai : pyFAI AzimuthalIntegrator
        Calibration object.

    Returns
    -------
    tuple
        (q, I) arrays resulting from integration.
    """

    return ai.integrate1d_ng(
        data,
        npt=1300,
        azimuth_range=cut,
        polarization_factor=0.99,
        correctSolidAngle=True,
        method='bbox',
        unit='q_nm^-1',
        safe=True
    )


def calculate_one_d_integration(img, ai, beam_time):
    """
    Compute equatorial and meridional 1D profiles from a diffraction image.

    Parameters
    ----------
    img : numpy.ndarray
        Input diffraction image.
    ai : pyFAI AzimuthalIntegrator
        Calibration object.
    beam_time : str
        Identifier of the experimental beamtime.

    Returns
    -------
    tuple
        q_eq, I_eq, q_me, I_me arrays (float64)
    """

    step = 5
    angle = 45

    # Beamtime-dependent orientation correction
    # NOTE:
    # These values are specific to the experimental setup used in this study
    # and may not be valid for other datasets.
    special_beamtimes = {
    '07_18', '2018_JUL', '2019_MAY',
    '2019_SEP', '11_21', '2021_NOV'
    }

    if beam_time in special_beamtimes:
        res_eq = d1_int(img, (180 + angle - step, 180 + angle + step), ai)
    else:
        res_eq = d1_int(img, (angle - step, angle + step), ai)

    # Meridional profile (perpendicular direction)
    res_me = d1_int(img, (90 + angle - step, 90 + angle + step), ai)

    # Convert outputs to float64 arrays
    (q_eq, I_eq), (q_me, I_me) = [
        np.array(x, dtype=np.float64)
        for x in (res_eq, res_me)
    ]

    return q_eq, I_eq, q_me, I_me