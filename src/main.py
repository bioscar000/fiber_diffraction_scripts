from pathlib import Path
import pandas as pd

from plotter.plot_utils import save_difractogram, save_profiles
# from utils.equatorial import calculate_equatorial
# from utils.equatorial_atlas import calculate_equatorial_parameters as cep_atlas
from utils.equatorial_rebe import calculate_equatorial_parameters
from utils.integration import load_calibration_file, calculate_one_d_integration
from utils.meridional import calculate_meridional_parameters
from utils.peaks_equatorial import calculate_peaks_parameters
from utils.preprocessing import load_edf, calculate_avg_image

disk = Path(r'C:\Users\Usuario\Nextcloud2\FD\FD_DATABASES\REBE_THESIS_PREPOSTPTX')
results = []

for bt in disk.iterdir():
    print(f"Beamtime - {bt}")

    for cal in bt.iterdir():
        print(f"Calibration - {cal}")
        calibration_path = cal / 'calibration.poni'
        poni = load_calibration_file(calibration_path)

        for buffer in cal.iterdir():

            if buffer.name == 'calibration.poni': continue
            print(f"Buffer - {buffer}")
            buffer_path = buffer / 'raw_buffer_files'
            avg_buffer = calculate_avg_image(buffer_path)

            for sample_batch in buffer.iterdir():

                if sample_batch.name == 'raw_buffer_files': continue
                print(f"Sample Batch - {sample_batch}")
                avg_sample = calculate_avg_image(sample_batch)
                avg_processed_sample = avg_sample - avg_buffer
                save_difractogram(avg_processed_sample, f'{bt.name}_{sample_batch.name}_AVGimg.tif')

                avg_q_eq, avg_i_eq, avg_q_me, avg_i_me = calculate_one_d_integration(avg_processed_sample, poni, bt)
                save_profiles(avg_q_eq, avg_i_eq, avg_q_me, avg_i_me, f'{bt.name}_{sample_batch.name}')

                avg_eq_params_dic = calculate_equatorial_parameters(avg_q_eq, avg_i_eq)
                avg_me_params_dic = calculate_meridional_parameters(avg_q_me, avg_i_me)
                avg_j0_params_dic = calculate_peaks_parameters(avg_q_eq, avg_i_eq, f'{bt.name}_{sample_batch.name}')

                base_meta = {
                             'beamtime': bt.name,
                             'calibration': cal.name,
                             'buffer': buffer.name,
                             'sample_batch': sample_batch.name
                            }

                meta_avg = {**base_meta, 'sample_file': f'AVG diffractogram {sample_batch.name}', 'rep':'AVG'}
                results.append({**meta_avg, **avg_eq_params_dic, **avg_me_params_dic, **avg_j0_params_dic})

                # rep = 1
                # for sample_file in sample_batch.iterdir():
                #     print(f"Sample - {sample_file}")
                #
                #     sample = load_edf(sample_file)
                #     processed_sample = sample - avg_buffer
                #     save_difractogram(processed_sample, f'{bt.name}_{sample_batch.name}_{rep}.tif')
                #
                #     q_eq, i_eq, q_me, i_me = calculate_one_d_integration(processed_sample, poni, bt)
                #     eq_params_dic = calculate_equatorial_parameters(q_eq, i_eq)
                #     me_params_dic = calculate_meridional_parameters(q_me, i_me)
                #     j0_params_dic = calculate_peaks_parameters(q_eq, i_eq)

                    # meta_single = {**base_meta, 'sample_file': sample_file.name, 'rep':rep}
                    # results.append({**meta_single, **eq_params_dic, **me_params_dic, **j0_params_dic})
                    # rep += 1

df = pd.DataFrame(results)
# df.to_csv("equatorial_results.csv", index=False)
df.to_csv(r"C:\Users\Usuario\Nextcloud2\rebe_thesis_oscar\gs244gtp_structural_parameters_ipfd_constrains.csv", index=False)