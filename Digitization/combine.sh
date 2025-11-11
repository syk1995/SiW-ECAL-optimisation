#!/bin/bash
COMBINE_X=5
COMBINE_Y=5
COMBINE_SI=3
COMBINE_LAYER=2
ABSORBER_LAYER=60
SCRIPT_PATH="/home/llr/ilc/shi/code/Energy-Reco/Digitization"
MY_SCRIPT="${SCRIPT_PATH}/combine_cells.py"

PARTICLE=gamma
Energy_val=(0.05 0.1 0.15 0.2 0.25 0.5 0.75 1.0 2.0 5.0 10.0 20.0 30.0 40.0 50.0 60.0)
#Energy_train=(0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 15.0 25.0 35.0 45.0 55.0 60.0)
Energy_train=(45.0)
ENERGY_LIST=("${Energy_train[@]}")
#PARTICLE=mu-
#ENERGY_LIST=(100.0)

INPUT_PATH="/data_ilc/flc/shi/SiWECAL-Prototype/Simu2025-06/CONF0/${PARTICLE}/Train"
OUTPUT_PATH="/data_ilc/flc/shi/SiWECAL-Prototype/Simu2025-06/CONF0/${PARTICLE}/Train"
CONDA_PATH=/home/llr/ilc/shi/miniconda3
. ${CONDA_PATH}/etc/profile.d/conda.sh
conda activate root_env
export PATH=/data_ilc/flc/shi/miniconda3/envs/root_env/bin:$PATH
which root
which python
root --version
python --version
for ENERGY in "${ENERGY_LIST[@]}"
do
    echo "Running for energy = ${ENERGY} GeV"
    python "${MY_SCRIPT}" \
        --CombineFactor_X $COMBINE_X \
        --CombineFactor_Y $COMBINE_Y \
        --CombineFactor_Si $COMBINE_SI \
        --CombineFactor_layer $COMBINE_LAYER \
        --Energy $ENERGY \
        --Absorber_layer $ABSORBER_LAYER \
        --input_path $INPUT_PATH \
        --output_path $OUTPUT_PATH
done
