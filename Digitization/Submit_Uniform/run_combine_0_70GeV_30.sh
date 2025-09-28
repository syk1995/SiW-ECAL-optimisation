#!/bin/bash
CONDA_PATH=/home/llr/ilc/shi/miniconda3
. /home/llr/ilc/shi/miniconda3/etc/profile.d/conda.sh
conda activate root_torch
export PATH=/data_ilc/flc/shi/miniconda3/envs/root_torch/bin:/data_ilc/flc/shi/miniconda3/envs/root_torch/bin:/data_ilc/flc/shi/miniconda3/bin:/data_ilc/flc/shi/miniconda3/condabin:/home/llr/ilc/shi/.local/bin:/sbin:/usr/sbin:/usr/local/sbin:/usr/share/Modules/bin:/usr/bin:/usr/externals/bin:/usr/sbin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/local/texlive/2019/bin/x86_64-linux:/usr/java/latest/bin:/home/llr/ilc/shi/bin
which root
which python
root --version
python --version
python /home/llr/ilc/shi/code/Energy-Reco/Digitization/combine_cells.py \
    --CombineFactor_X 1 \
    --CombineFactor_Y 1 \
    --CombineFactor_Si 3 \
    --CombineFactor_layer 2 \
    --Absorber_layer 60 \
    --input_file "/home/llr/ilc/shi/data/SiWECAL-Prototype/Simu2025-06/CONF0/gamma/Train/MC/Uniform/0_70GeV_30.root" \
    --output_file "/home/llr/ilc/shi/data/SiWECAL-Prototype/Simu2025-06/CONF0/gamma/Train/Merged_X5.0mm_Y5.0mm_Si0.45mm_layer30_in60/Uniform/0_70GeV_30.root"
