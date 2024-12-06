# Using Probabilistic Graphical Models to Deconvolute Network Signaling

## How to run each method:

1. Loopy BP
   - Script: ./loopy_bp/loopy_bp_final.py
   - Input: 1. prior network file path; 2. gene state file (derived from phosphosite data); 3. output file path
   - Ouput: output file path/edge_data.csv
2. Optimize Loopy BP
   - Script: ./loopy_bp/optimize_loopy_bp.slurm
   - Input: 1. see Loopy BP section 2. dictionary of genes of interest and list of correspond first neighborhood (module of interest): moi_dict.pkl
   - Output: mu_score.csv; use optimize_mu_plot.py to create performance plot.
2. Bionic
   - Script: ./Bionic/submit.bionic.job.sbatch
   - Input: bionic.config.file.json; See Bionic paper for implementation details: https://www.nature.com/articles/s41592-022-01616-x#code-availability
   - Output: bionic_features.tsv
3. Gene State
   - Script: ./TLR_data/format_TLR_states.ipynb
   - Input: Differential P-site expression data (pre-filtered)
   - Output: P-site up or down regulation: TLR_states.csv
4. Evaluation
   - Script: ./20241203/control_module_eval.ipynb
   - Input: each network you would like to evaluate and the control pathway of interest network
   - Output: first neighborhood subnetworks for each gene of interest you query
