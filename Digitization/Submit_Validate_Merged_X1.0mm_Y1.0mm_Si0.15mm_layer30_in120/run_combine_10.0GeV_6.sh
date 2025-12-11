#!/bin/bash
. /home/llr/ilc/shi/env/root_torch_condor.sh
which root
which python
root --version
python --version
python /home/llr/ilc/shi/code/Energy-Reco/Digitization/combine_cells.py \
    --CombineFactor_X 1 \
    --CombineFactor_Y 1 \
    --CombineFactor_Si 1 \
    --CombineFactor_layer 4 \
    --Absorber_layer 120 \
    --input_file "/home/llr/ilc/shi/data/SiWECAL-Prototype/Simu2025-06/CONF4/gamma/Validate/MC/10.0GeV_6.root" \
    --output_file "/home/llr/ilc/shi/data/SiWECAL-Prototype/Simu2025-06/CONF4/gamma/Validate/Merged_X1.0mm_Y1.0mm_Si0.15mm_layer30_in120/10.0GeV_6.root"
