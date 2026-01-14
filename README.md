# spGen — Synthetic Generator of Single Channel sp-ICP-MS Time Scans

**Author**: Pierre-Emmanuel Peyneau  
**Version**: 1.0  
**Date**: 2025-06-05  
**License**: CC-BY 4.0

---

## 🧪 Overview

`spGen` is a Python-based stochastic simulation tool for generating **realistic synthetic time scans** for a single mass-to-charge ratio. It models time series typically encountered in single-particle inductively coupled plasma-mass spectrometry (sp-ICP-MS).

---

## 📁 Inputs

Inputs are specified via an `input.yaml` file in each simulation folder. Below are the required keys and their descriptions:

| Key | Description | Type | Units |
|-----|-------------|------|-------|
| `id` | Identifier for the simulation | `str` | — |
| `n_reads` | Number of readings | `int` | — |
| `dwell_time` | Dwell time | `float` | s |
| `f_transmission` | Ion transmission probability | `float` | — |
| `flux_dissolved` | Analyte ion flux from dissolved species | `float` | s⁻¹ |
| `flux_particles` | Nanoparticle flux | `float` | s⁻¹ |
| `detector` | Type of detector (`ideal`, `non-ideal`) | `str` | — |
| `sir_mean` | Mean of log(single-ion response) | `float` | — |
| `sir_std` | Std dev of log(single-ion response) | `float` | — |
| `size_distribution` | NP size distribution (`dirac`, `lognormal`, `truncated-normal`) | `str` | — |
| `size_particle_mean` | Mean NP diameter | `float` | cm |
| `size_particle_std` | Std dev of NP diameter | `float` | cm |
| `shape_distribution` | Time dispersion distribution (`inverse-gaussian`, `uniform`) | `str` | — |
| `mean` | Mean transport time | `float` | s |
| `scale` | Scale parameter (if inverse-gaussian) | `float` | s |
| `seed` | RNG seed | `int` | — |
| `molar_mass` | Molar mass of isotope | `float` | g/mol |
| `isotopic_abundance` | Isotope abundance | `float` | — |
| `mass_fraction` | Element mass fraction in NP | `float` | — |
| `mass_density` | NP mass density | `float` | g/cm³ |
| `timescan_csv` | Save time scan as CSV | `bool` | — |
| `timescan_csv_header` | Include header in CSV | `bool` | — |
| `timescan_csv_start` | Start time for CSV | `float` | s |

---

## 📝 Example input.yaml

```yaml
id: example01
n_reads: 600000
dwell_time: 0.0001
f_transmission: 0.0001
flux_dissolved: 5.0e+8
flux_particles: 10
detector: non-ideal
sir_mean: -0.11045
sir_std: 0.47
size_distribution: lognormal
size_particle_mean: 60.0e-7
size_particle_std: 10.0e-7
shape_distribution: inverse-gaussian
mean: 0.0010
scale: 0.0400
seed: 1
molar_mass: 197.0
isotopic_abundance: 1.0
mass_fraction: 1.0
mass_density: 19.3
timescan_csv: true
timescan_csv_header: true
timescan_csv_start: 0.0
```
---

## 🚀 How to Run the Program

### 🖥️ Prerequisites

* Python 3.11.5+
* Required Python packages:

  * `numpy`
  * `scipy`
  * `pyyaml`

Install them with pip if needed, using the following command:

```bash
pip install numpy scipy pyyaml
```

---

### 🗂️ Prepare Your Input

1. Create a directory (e.g., `example_run/`).
2. Place your `input.yaml` file in that directory.

   * You can use the [example template](#-example-inputyaml) from the README.

Example structure:

```
example_run/
├── input.yaml
```

---

### ▶️ Run the Program

Use the command-line interface to point to one or more input folders:

```bash
python spGen.py --folder example_run
```

> ✅ You can also run multiple simulations at once:

```bash
python spGen.py --folder example_run example_run2 example_run3
```

---

### 📂 Output Files

Each input folder will receive:

| File                | Description                              |
| ------------------- | ---------------------------------------- |
| `output.yaml`       | Summary of the simulation run            |
| `timescan+<id>.csv` | Time scan data (if `timescan_csv: True`) |
| `events+<id>.csv`   | Particle event time windows              |

----

## 🔍 Citation

If you use this tool in scientific work, please cite the publication mentioned below and this repository.

Peyneau, P. E., Seydoux, L., & Tharaud, M. (2026). Synthetic generation of single-channel single particle ICP-MS time scans.  _J. Anal. At. Spectrom._, 2026, 41, 320-332. DOI: 10.1039/D5JA00232J
