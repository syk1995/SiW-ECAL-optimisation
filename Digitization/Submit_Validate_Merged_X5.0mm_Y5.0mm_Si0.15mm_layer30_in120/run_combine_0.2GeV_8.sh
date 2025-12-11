#!/bin/bash
. /home/llr/ilc/shi/env/root_torch_condor.sh
which root
which python
root --version
python --version
python /home/llr/ilc/shi/code/Energy-Reco/Digitization/combine_cells.py \
    --CombineFactor_X 5 \
    --CombineFactor_Y 5 \
    --CombineFactor_Si 1 \
    --CombineFactor_layer 4 \
    --Absorber_layer 120 \
    --input_file "/home/llr/ilc/shi/data/SiWECAL-Prototype/Simu2025-06/CONF4/gamma/Validate/MC/0.2GeV_8.root" \
    --output_file "/home/llr/ilc/shi/data/SiWECAL-Prototype/Simu2025-06/CONF4/gamma/Validate/Merged_X5.0mm_Y5.0mm_Si0.15mm_layer30_in120/0.2GeV_8.root"
