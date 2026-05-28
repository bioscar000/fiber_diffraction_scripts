"""
Fiber Diffraction Analysis Pipeline

Author: Óscar Fernández Blanco
Contributors: Rebeca París Ogayar

Description:
This is the main execution entry point of the pipeline

Pipeline for processing fiber diffraction datasets:
- Image averaging
- Buffer subtraction
- 1D integration (equatorial and meridional)
- Structural parameter extraction

Notes:
This code was developed for a specific dataset and is not intended as a general-purpose tool.
"""

__author__ = "Óscar Fernández Blanco"
__contact__ = "oscar.fernandez@cib.csic.es"
__credits__ = ["Rebeca París Ogayar"]
__license__ = "MIT"
__status__ = "Research code (final version for thesis)"


from pathlib import Path
import pandas as pd
# Preprocessing utilities
from utils.preprocessing import calculate_avg_image
# Integration utilities
from utils.integration import load_calibration_file, calculate_one_d_integration
# Structural parameters utilities
from utils.equatorial_rebe import calculate_equatorial_parameters
from utils.meridional import calculate_meridional_parameters
from utils.peaks_equatorial import calculate_peaks_parameters
# Plotting utilities
from plotter.plot_utils import save_difractogram, save_profiles


def process_dataset(root_path: Path):
    """
    Process all datasets contained in the given root directory.

    Parameters
    ----------
    root_path : Path
        Root directory containing beamtime folders.

    Returns
    -------
    list of dict
        Extracted structural parameters with metadata.
    """

    results = []

    for bt in root_path.iterdir():
        print(f"Beamtime: {bt}")

        for cal in bt.iterdir():
            print(f"Calibration: {cal}")

            calibration_path = cal / 'calibration.poni'
            poni = load_calibration_file(calibration_path)

            for buffer in cal.iterdir():

                if buffer.name == 'calibration.poni':
                    continue

                print(f"Buffer: {buffer}")

                buffer_path = buffer / 'raw_buffer_files'
                avg_buffer = calculate_avg_image(buffer_path)

                for sample_batch in buffer.iterdir():

                    if sample_batch.name == 'raw_buffer_files':
                        continue

                    print(f"Sample Batch: {sample_batch}")

                    avg_sample = calculate_avg_image(sample_batch)
                    avg_processed_sample = avg_sample: avg_buffer

                    save_difractogram(
                        avg_processed_sample,
                        f'{bt.name}_{sample_batch.name}_AVGimg.tif'
                    )

                    avg_q_eq, avg_i_eq, avg_q_me, avg_i_me = calculate_one_d_integration(
                        avg_processed_sample,
                        poni,
                        bt
                    )

                    save_profiles(
                        avg_q_eq,
                        avg_i_eq,
                        avg_q_me,
                        avg_i_me,
                        f'{bt.name}_{sample_batch.name}'
                    )

                    avg_eq_params_dic = calculate_equatorial_parameters(avg_q_eq, avg_i_eq)
                    avg_me_params_dic = calculate_meridional_parameters(avg_q_me, avg_i_me)
                    avg_j0_params_dic = calculate_peaks_parameters(
                        avg_q_eq,
                        avg_i_eq,
                        f'{bt.name}_{sample_batch.name}'
                    )

                    base_meta = {
                        'beamtime': bt.name,
                        'calibration': cal.name,
                        'buffer': buffer.name,
                        'sample_batch': sample_batch.name
                    }

                    meta_avg = {
                        **base_meta,
                        'sample_file': f'AVG diffractogram {sample_batch.name}',
                        'rep': 'AVG'
                    }

                    results.append({
                        **meta_avg,
                        **avg_eq_params_dic,
                        **avg_me_params_dic,
                        **avg_j0_params_dic
                    })

    return results


def main():
    """
    Main execution function.
    """

    # TODO: Replace with actual dataset path (see README for structure)
    root_path = Path("../data_example")

    results = process_dataset(root_path)

    df = pd.DataFrame(results)

    output_path = Path("../outputs/structural_parameters.csv")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(output_path, index=False)

    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()