#!/bin/bash
. /home/llr/ilc/shi/env/root_torch_condor.sh
which root
which python
root --version
python --version
python /home/llr/ilc/shi/code/Energy-Reco/Digitization/combine_cells.py \
    --CombineFactor_X 3 \
    --CombineFactor_Y 3 \
    --CombineFactor_Si 1 \
    --CombineFactor_layer 4 \
    --Absorber_layer 120 \
    --input_file "/home/llr/ilc/shi/data/SiWECAL-Prototype/Simu2025-06/CONF4/gamma/Validate/MC/1.0GeV_7.root" \
    --output_file "/home/llr/ilc/shi/data/SiWECAL-Prototype/Simu2025-06/CONF4/gamma/Validate/Merged_X3.0mm_Y3.0mm_Si0.15mm_layer30_in120/1.0GeV_7.root"
