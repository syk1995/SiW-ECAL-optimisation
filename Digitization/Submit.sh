#!/bin/bash
pwd=$PWD
SCRIPT_PATH=/home/llr/ilc/shi/code/Energy-Reco/Digitization
MY_SCRIPT=${SCRIPT_PATH}/combine_cells.py
PARTICLE=gamma
COMBINE_X=1
COMBINE_Y=1
COMBINE_SI=3
SI_THICKNESS=0.15
COMBINE_LAYER=2
ABSORBER_LAYER=60
MERGE_NAME=$(printf "Merged_X%.1fmm_Y%.1fmm_Si%.2fmm_layer%.0f_in%d" \
    "$(echo "5*$COMBINE_X" | bc -l)" \
    "$(echo "5*$COMBINE_Y" | bc -l)" \
    "$(echo "$COMBINE_SI * $SI_THICKNESS" | bc -l)" \
    "$(echo "$ABSORBER_LAYER / $COMBINE_LAYER" | bc -l)" \
    "$ABSORBER_LAYER")

echo $MERGE_NAME
DATA_PATH="/home/llr/ilc/shi/data/SiWECAL-Prototype/Simu2025-06/CONF0/${PARTICLE}"
DATA_TYPE=Uniformq  # Train, Validate, Uniform
case $DATA_TYPE in
    Train)
        ENERGY_LIST=("${Energy_train[@]}")
        INPUT_PATH="${DATA_PATH}/${DATA_TYPE}/MC"
        OUTPUT_PATH="${DATA_PATH}/${DATA_TYPE}/${MERGE_NAME}"
        ;;
    Validate)
        ENERGY_LIST=("${Energy_val[@]}")
        INPUT_PATH="${DATA_PATH}/${DATA_TYPE}/MC"
        OUTPUT_PATH="${DATA_PATH}/${DATA_TYPE}/${MERGE_NAME}"
        ;;
    Uniform)
        ENERGY_LIST=("${Energy_uniform[@]}")
        INPUT_PATH="${DATA_PATH}/Train/MC/${DATA_TYPE}"
        OUTPUT_PATH="${DATA_PATH}/Train/${MERGE_NAME}/${DATA_TYPE}"
        ;;
    *)
        echo "[ERROR] Unknown DATA_TYPE: $DATA_TYPE"
        exit 1
        ;;
esac
mkdir -p "$OUTPUT_PATH"
echo "DATA TYPE: $DATA_TYPE"
echo "INPUT PATH: $INPUT_PATH"
echo "OUTPUT PATH: $OUTPUT_PATH"
# ===================== Submit Directory =====================
Submit_DIR="${SCRIPT_PATH}/Submit_${DATA_TYPE}"
rm -rf "$Submit_DIR"
mkdir -p "$Submit_DIR"
cd $Submit_DIR
# ===================== Generate Jobs =====================
for input_file in "${INPUT_PATH}"/*.root; do
    file_name=$(basename "$input_file" .root)    
    JOB_SCRIPT="${Submit_DIR}/run_combine_${file_name}.sh"
    output_file="${OUTPUT_PATH}/${file_name}.root"
    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
CONDA_PATH=/home/llr/ilc/shi/miniconda3
. ${CONDA_PATH}/etc/profile.d/conda.sh
conda activate root_env
export PATH=/data_ilc/flc/shi/miniconda3/envs/root_env/bin:$PATH
which root
which python
root --version
python --version
python ${MY_SCRIPT} \\
    --CombineFactor_X $COMBINE_X \\
    --CombineFactor_Y $COMBINE_Y \\
    --CombineFactor_Si $COMBINE_SI \\
    --CombineFactor_layer $COMBINE_LAYER \\
    --Absorber_layer $ABSORBER_LAYER \\
    --input_file "$input_file" \\
    --output_file "$output_file"
EOF

    chmod +x "$JOB_SCRIPT"
done


# ===================== Submit All =====================
for job in ${Submit_DIR}/*.sh; do
    echo "Submitting job: $job"
    t3submit -singleout "$job"
done
cd $pwd
