import fabio
import numpy as np

def calculate_avg_image(folder_path):
    imgs = [load_edf(f) for f in folder_path.iterdir()]  # load_edf hace el chequeo
    avg_image = np.mean(imgs, axis=0)
    return avg_image

def load_edf(edf_path):
    assert edf_path.suffix == '.edf', f"Only .edf files are allowed: {edf_path}"
    return fabio.open(edf_path).data.astype(np.float64)

# def calculate_avg_buffer_welford(buffer_path):
#     # Welford algorithm
#     buffer_files = [f for f in os.listdir(buffer_path) if f.endswith('.edf')]
#     avg_buffer = None
#     n = 0
#     for buffer_file in sorted(buffer_files):
#         # print(buffer_file)
#
#         img = fabio.open(os.path.join(buffer_path, buffer_file)).data
#         img_float = img.astype(np.float64)
#
#         if avg_buffer is None:
#             avg_buffer = np.zeros_like(img_float, dtype=np.float64)
#
#         n += 1
#         avg_buffer += (img_float - avg_buffer) / n
#
#     return avg_buffer