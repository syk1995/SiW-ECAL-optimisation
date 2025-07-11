import os
import h5py
import numpy as np
import torch
import zarr
import matplotlib.pyplot as plt
from torch.utils.data import Dataset, DataLoader
import torch.nn as nn
from train import DualBranchFEMNet
from train import CaloDataset

def plot_scatter(x, y, xlabel, ylabel, title, filename):
    plt.figure(figsize=(6, 6))
    plt.hist2d(x, y, bins=100, range=[[0, 1], [0, 1]], cmap='Blues', density=True)
    plt.plot([0, 1], [0, 1], 'r--', linewidth=1)  # 参考线
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.colorbar(label='Density')
    plt.savefig(filename)
    plt.close()

def plot_histogram(data, xlabel, ylabel, title, filename):
    plt.figure(figsize=(8, 6))
    plt.hist(data, bins=100, alpha=0.75)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def run_inference(h5_path, model_path, public_path, max_samples=10000, batch_size=512):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"📦 Using device: {device}")

    zarr_data = zarr.open(zarr_path, mode='r')
    total = len(zarr_data['ecal'])

    start = int(total * 0.8)
    end = min(total, start + max_samples)
    test_indices = np.arange(start, end)

    dataset = CaloDataset(h5_path, batch_size=batch_size, indices=test_indices)
    dataloader = DataLoader(dataset, batch_size=1, shuffle=False, num_workers=4)

    model = DualBranchFEMNet()
    model.load_state_dict(torch.load(model_path, map_location=device))
    model.to(device)
    model.eval()

    ecal_fem_true, ecal_fem_pred = [], []
    hcal_fem_true, hcal_fem_pred = [], []

    with torch.no_grad():
        for batch in dataloader:
            batch = {k: v.squeeze(0).to(device) for k, v in batch.items()}
            ecal_fem_pred_batch, hcal_fem_pred_batch = model(
                batch['ecal'], batch['hcal'],
                batch['energy_ecal'], batch['energy_hcal']
            )
            ecal_fem_true.append(batch['fem_ecal'].cpu().numpy())
            ecal_fem_pred.append(ecal_fem_pred_batch.cpu().numpy())
            hcal_fem_true.append(batch['fem_hcal'].cpu().numpy())
            hcal_fem_pred.append(hcal_fem_pred_batch.cpu().numpy())

    ecal_fem_true = np.concatenate(ecal_fem_true)
    ecal_fem_pred = np.concatenate(ecal_fem_pred)
    hcal_fem_true = np.concatenate(hcal_fem_true)
    hcal_fem_pred = np.concatenate(hcal_fem_pred)
    ecal_errors = ecal_fem_pred - ecal_fem_true
    hcal_errors = hcal_fem_pred - hcal_fem_true

    os.makedirs(public_path, exist_ok=True)

    # 绘制散点图
    plot_scatter(
        ecal_fem_true, ecal_fem_pred,
        "ECAL True FEM", "ECAL Predicted FEM",
        "ECAL FEM Prediction vs True",
        os.path.join(public_path, "ecal_fem_scatter.png")
    )
    plot_scatter(
        hcal_fem_true, hcal_fem_pred,
        "HCAL True FEM", "HCAL Predicted FEM",
        "HCAL FEM Prediction vs True",
        os.path.join(public_path, "hcal_fem_scatter.png")
    )
    # 绘制直方图
    plot_histogram(
        ecal_errors, "ECAL FEM Prediction Error", "Frequency",
        "ECAL FEM Prediction Error Distribution",
        os.path.join(public_path, "ecal_fem_error_hist.png")
    )
    plot_histogram(
        hcal_errors, "HCAL FEM Prediction Error", "Frequency",
        "HCAL FEM Prediction Error Distribution",
        os.path.join(public_path, "hcal_fem_error_hist.png")
    )

    # 保存预测结果
    np.save(os.path.join(public_path,"fem_ecal_results.npz"), {
        'true': ecal_fem_true,
        'pred': ecal_fem_pred,
        'errors': ecal_errors
    })
    np.save(os.path.join(public_path,"fem_hcal_results.npz"), {
        'true': hcal_fem_true,
        'pred': hcal_fem_pred,
        'errors': hcal_errors
    })


    print("✅ Inference completed!")

# === 启动 ===
if __name__ == "__main__":
    public_path = "/gpfs/users/xiax/SoftwareCompensation/pytorch/project_2fem/output/"
    zarr_path = "/gpfs/workdir/xiax/dualdata/fem2_05MIP_data/voxelized_data.zarr" 
    model_path = os.path.join(public_path, "best_model.pth")
    run_inference(zarr_path, model_path, public_path, max_samples=200000, batch_size=512)



