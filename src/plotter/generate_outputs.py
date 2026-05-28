"""
Plotting and output utilities for fiber diffraction analysis

Author: Óscar Fernández Blanco

Description:
This module provides utility functions for saving diffraction images
and integrated scattering profiles to disk.

Notes:
- Output paths are relative and assume a standard project structure.
- Directories are created automatically if they do not exist.

Part of the fiber diffraction analysis pipeline developed for a PhD thesis.
"""

__author__ = "Óscar Fernández Blanco"
__license__ = "MIT"


from pathlib import Path

import tifffile
import pandas as pd


def save_difractogram(array, name, output_dir=Path("../outputs/difractograms")):
    """
    Save a diffraction image to disk.

    Parameters
    ----------
    array : numpy.ndarray
        Image data.
    name : str
        Output file name.
    output_dir : Path, optional
        Directory where the image will be saved.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tifffile.imwrite(output_dir / name, array)


def save_profiles(q_eq, i_eq, q_me, i_me, name, output_dir=Path("../outputs/profiles")):
    """
    Save equatorial and meridional profiles to a CSV file.

    Parameters
    ----------
    q_eq : numpy.ndarray
        Equatorial q values.
    i_eq : numpy.ndarray
        Equatorial intensities.
    q_me : numpy.ndarray
        Meridional q values.
    i_me : numpy.ndarray
        Meridional intensities.
    name : str
        Base name for the output file.
    output_dir : Path, optional
        Directory where the file will be saved.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame({
        'q EQ': q_eq,
        'i EQ': i_eq,
        'q ME': q_me,
        'i ME': i_me,
    })

    output_path = output_dir / f"{name}_PROFILES.dat"

    df.to_csv(
        output_path,
        index=False,
        float_format="%.6f",
        sep=","
    )