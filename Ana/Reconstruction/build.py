from torch.utils.cpp_extension import load
import os

# C++ 源文件
src_file = "create_graph.cpp"

# nanoflann 的 include 路径（conda 安装路径）
nanoflann_include = "/data_ilc/flc/shi/miniconda3/envs/GNN/include"

# PyTorch include 路径
torch_include = [
    os.path.join(os.path.dirname(__import__("torch").__file__), "include"),
    os.path.join(os.path.dirname(__import__("torch").__file__), "include/torch/csrc/api/include")
]
root_include = "/data_ilc/flc/shi/miniconda3/envs/GNN/include"
root_lib = "/data_ilc/flc/shi/miniconda3/envs/GNN/lib"
root_libs = ["Core", "RIO", "Net", "Hist", "Tree", "Graf", "Gpad", "TreePlayer"]

# 构建 C++ 扩展
CreateGraph = load(
    name="create_graph",
    sources=[src_file],
    extra_cflags=["-std=c++17"] + [f"-I{nanoflann_include}", f"-I{root_include}"] + [f"-I{p}" for p in torch_include],
    extra_ldflags=[f"-L{root_lib}"] + [f"-l{lib}" for lib in root_libs],
    build_directory=os.getcwd(),
    verbose=True
)

