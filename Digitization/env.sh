
# conda deactivate
# if [[ $GCC_VERSION -lt 13 ]]; then
#    module reset
#    source /usr/share/Modules/init/sh
#    module use /opt/exp_soft/vo.llr.in2p3.fr/modulefiles_el9
#    module load python/3.12.4
#    module load compilers/gcc/13.x.x
#    #source /opt/exp_soft/llr/root/v6.32-el9-gcc13xx-py3124/etc/init.sh
# else
#    echo "GCC version is already $GCC_VERSION, no need to load."
# fi


CONDA_PATH=/home/llr/ilc/shi/miniconda3
. ${CONDA_PATH}/etc/profile.d/conda.sh
conda activate root_env
export PATH=/data_ilc/flc/shi/miniconda3/envs/root_env/bin:$PATH
which root
which python
root --version
python --version
