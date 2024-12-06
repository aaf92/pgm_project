#!/bin/bash

#SBATCH --job loopy_bp
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH --time=6-00:00
#SBATCH --mem=8g

echo Running on `hostname`
echo workdir $SLURM_SUBMIT_DIR


module purge
module load python/ondemand-jupyter-python3.10

cd /ix/djishnu/Aaron_F/PGM_project/loopy_bp

python3 /ix/djishnu/Aaron_F/PGM_project/loopy_bp/optimize_mu_loopy_bp.py \
    /ix/djishnu/Aaron_F/PGM_project/PPI/MusMusculus_lcb_hq.txt \
    /ix/djishnu/Aaron_F/PGM_project/TLR_data/Results/TLR_states.csv \
    /ix/djishnu/Aaron_F/PGM_project/loopy_bp

echo "done!"
