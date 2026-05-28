"""
Preprocessing utilities for fiber diffraction images

Author: Óscar Fernández Blanco

Description:
This module provides basic preprocessing functions for fiber diffraction data,
including loading EDF images and computing average images from a dataset.

Notes:
Part of the fiber diffraction analysis pipeline developed for a PhD thesis.
"""

__author__ = "Óscar Fernández Blanco"
__license__ = "MIT"


import numpy as np
import fabio


def load_edf(edf_path):
    """
    Load an EDF image file.

    Parameters
    ----------
    edf_path : Path
        Path to the EDF file.

    Returns
    -------
    numpy.ndarray
        Image data as a float64 array.

    Raises
    ------
    AssertionError
        If the file is not in EDF format.
    """

    assert edf_path.suffix == '.edf', f"Only .edf files are allowed: {edf_path}"
    return fabio.open(edf_path).data.astype(np.float64)


def calculate_avg_image(folder_path):
    """
    Compute the average image from a folder of EDF files.

    Parameters
    ----------
    folder_path : Path
        Directory containing EDF images.

    Returns
    -------
    numpy.ndarray
        Averaged image.
    """

    # Load all images in folder
    imgs = [load_edf(f) for f in folder_path.iterdir()]

    # Compute pixel-wise average
    avg_image = np.mean(imgs, axis=0)

    return avg_image