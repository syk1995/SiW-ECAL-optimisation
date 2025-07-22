#!/bin/bash
COMBINE_X=2
COMBINE_Y=2
COMBINE_SI=5
COMBINE_LAYER=2
SCRIPT_PATH="/home/llr/ilc/shi/code/Energy-Reco/Digitization"
MY_SCRIPT="${SCRIPT_PATH}/combine_cells.py"
ENV="${SCRIPT_PATH}/env.sh"
#ENV=/home/llr/ilc/shi/code/SiWECAL-Prototype/generation/run_scripts/Simu2025-06/env.sh
INPUT_PATH="/data_ilc/flc/shi/SiWECAL-Prototype/Simu2025-06/CONF0/mu-"
OUTPUT_PATH="/data_ilc/flc/shi/SiWECAL-Prototype/Simu2025-06/CONF0/mu-"
ENERGY_LIST=(200.0)
. "${ENV}"
#export PYTHONPATH=/data_ilc/flc/shi/miniconda3/envs/my_notebook_env/lib/python3.10/site-packages/typing_extensions.py:${PYTHONPATH}
which root
which python
for ENERGY in "${ENERGY_LIST[@]}"
do
    echo "Running for energy = ${ENERGY} GeV"
    python "${MY_SCRIPT}" \
        --CombineFactor_X $COMBINE_X \
        --CombineFactor_Y $COMBINE_Y \
        --CombineFactor_Si $COMBINE_SI \
        --CombineFactor_layer $COMBINE_LAYER \
        --Energy $ENERGY \
        --input_path $INPUT_PATH \
        --output_path $OUTPUT_PATH
done
