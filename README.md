# Fiber Diffraction Analysis Pipeline

## Description

This repository contains a Python-based pipeline developed for the analysis of fiber diffraction data in the context of a PhD thesis.

The code performs a complete workflow including:
- Averaging of diffraction images
- Buffer subtraction
- 1D azimuthal integration (equatorial and meridional)
- Extraction of structural parameters
- Peak detection and refinement

The implementation is tailored to specific experimental datasets and beamtime configurations.

---

## Authors

- **Óscar Fernández Blanco** – Development of the analysis pipeline, data processing, and physical modeling.
- **Rebeca París Ogayar** – Design and development of the experimental probes, sample optimization, and synchrotron data acquisition (PhD thesis).

---

## Installation

### Requirements

Python 3.8+ is recommended.

Install dependencies using:

```bash
pip install -r requirements.txt
```

### Main dependencies

- numpy  
- scipy  
- pandas  
- pyFAI  
- fabio  
- tifffile  

---

## Project Structure

```
fiber_diffraction_scripts/
│
├── main.py
├── utils/
│   ├── preprocessing.py
│   ├── integration.py
│   ├── equatorial_model.py
│   ├── meridional.py
│   ├── peaks_equatorial.py
│
├── plotter/
│   └── plot_utils.py
│
├── outputs/
│   ├── difractograms/
│   ├── profiles/
│   └── structural_parameters.csv
│
└── README.md
```

---

## Usage

Edit the input data path in `main.py`:

```python
root_path = Path("path_to_data")
```

Then run:

```bash
python main.py
```

Results will be saved in:

```
outputs/
```

---

## Output

The pipeline generates:

- Diffractograms (`.tif`)  
- 1D profiles (`.dat`)  
- Structural parameters (`structural_parameters.csv`)  

---

## Important Notes

- This code was developed for a **specific experimental setup and dataset**.  
- Several parameters (e.g., beamtime-dependent corrections and fitting ranges) are **not generalizable**.  
- Users working with different data will likely need to adapt:
  - Angular ranges in integration  
  - Peak detection parameters  
  - Model fitting constraints  

---

## Reproducibility

This repository contains the **frozen version of the code used in the PhD thesis**.  

No further development or maintenance is intended.

---

## License

This project is licensed under the MIT License.

---

## Citation

If you use this code, please cite:

```
Fernández, Ó. (2026)  
Fiber diffraction analysis pipeline  
DOI: [to be added after Zenodo deposition]
```

---

## Contact

For questions related to the code:

Óscar Fernández Blanco  
oscar.fernandez@cib.csic.es
