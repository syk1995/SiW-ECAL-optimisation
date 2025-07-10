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

# Activate anaconda environment code
source ~/conda_rc
conda activate tf.ac24

# Train the network
# python /gpfs/users/xiax/SoftwareCompensation/Test.py
# python /gpfs/users/xiax/SoftwareCompensation/Train_step.py
# python /gpfs/users/xiax/SoftwareCompensation/Train.py
py-spy record \
  --output /gpfs/users/xiax/SoftwareCompensation/result/pyspy_train_10.svg \
  --native \
  -- python /gpfs/users/xiax/SoftwareCompensation/Train.py
