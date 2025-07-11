#!/bin/bash

#SBATCH --job-name=ruche_keras_train
#SBATCH --output=%x.o%j 
#SBATCH --time=24:00:00 
#SBATCH --ntasks=10
#SBATCH --gres=gpu:1 
#SBATCH --partition=gpu

# [ ! -d output ] && mkdir output

# Module load 
module load anaconda3/2024.06/gcc-13.2.0
module load cuda/12.8.0/gcc-13.2.0

# Activate anaconda environment code
source ~/conda_rc
conda activate MLenv

# Train the network
# python -u /gpfs/users/xiax/SoftwareCompensation/pytorch/project_2fem/train.py
python -u /gpfs/users/xiax/SoftwareCompensation/pytorch/project_2fem/analysis.py