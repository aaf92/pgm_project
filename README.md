# Using Probabilistic Graphical Models to Deconvolute Network Signaling

## How to Run Each Method

### 1. Gene State
- **Script**: `./TLR_data/format_TLR_states.ipynb`
- **Input**: Differential P-site expression data (pre-filtered)
- **Output**: P-site up or down regulation: `TLR_states.csv`

---

### 2. Loopy BP
- **Script**: `./loopy_bp/loopy_bp_final.py`
- **Input**:
  1. Prior network file path  
  2. Gene state file (derived from phosphosite data)  
  3. Output file path
- **Output**: `<output_file_path>/edge_data.csv`

---

### 3. Optimize Loopy BP
- **Script**: `./loopy_bp/optimize_loopy_bp.slurm`
- **Input**:
  1. See **Loopy BP** section  
  2. Dictionary of genes of interest and their first neighborhood (module of interest): `moi_dict.pkl`
- **Output**: `mu_score.csv`  
  Use `optimize_mu_plot.py` to create a performance plot.

---

### 4. Bionic
- **Script**: `./Bionic/submit.bionic.job.sbatch`
- **Input**: `bionic.config.file.json`  
  See the [Bionic paper](https://www.nature.com/articles/s41592-022-01616-x#code-availability) for implementation details.
- **Output**: `bionic_features.tsv`

---

### 5. Evaluation
- **Script**: `./20241203/control_module_eval.ipynb`
- **Input**: 
  - Each network you would like to evaluate  
  - The control pathway of interest network
- **Output**: First neighborhood subnetworks for each gene of interest you query
