#!/bin/bash

#SBATCH --job bionic
#SBATCH --nodes=1
#SBATCH --cluster=gpu
##SBATCH --partition=smp  
#SBATCH --gres=gpu:1
#SBATCH --partition=a100
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH --time=6-00:00
#SBATCH --mem=238419

echo Running on `hostname`
echo workdir $SLURM_SUBMIT_DIR


module purge
module load cuda/11.1.1
#module load python/3.7.0 venv/wrap

#workon bionic 

module load python/ondemand-jupyter-python3.8

#module load bionic/0.2.5_a100

#source activate /ihome/djishnu/prc44/anaconda3/envs/bionic_gpu_3_8_2023

# source activate /ix/djishnu/Aaron_F/envs/bionic
source activate /ihome/djishnu/prc44/.conda/envs/bionic_gpu

#conda install pytorch==1.9.0 torchvision==0.10.0 torchaudio==0.9.0 cudatoolkit=11.1 -c pytorch -c conda-forge


#export LD_LIBRARY_PATH="/ihome/djishnu/prc44/.conda/envs/bionic_gpu/lib:$LD_LIBRARY_PATH"

cd /ix/djishnu/Aaron_F/PGM_project/Bionic

bionic bionic.config.file.json

echo "done!"


