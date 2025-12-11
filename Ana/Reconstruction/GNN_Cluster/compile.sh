conda activate GNN
g++ -O3 -Wall -shared  -fPIC     `python3 -m pybind11 --includes`     -I/data_ilc/flc/shi/miniconda3/envs/GNN/lib/python3.10/site-packages/torch/include     -I/data_ilc/flc/shi/miniconda3/envs/GNN/lib/python3.10/site-packages/torch/include/torch/csrc/api/include     BuildCluster.cpp -o ReadRootClusters`python3-config --extension-suffix`
