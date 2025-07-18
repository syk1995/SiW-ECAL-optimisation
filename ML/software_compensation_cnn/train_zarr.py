import os
import time
import zarr
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, Subset
import matplotlib.pyplot as plt
import json
import time

class CaloDataset(Dataset):
    def __init__(self, zarr_path, batch_size, indices=None):
        self.zarr = zarr.open(zarr_path, mode='r')
        self.batch_size = batch_size
        self.indices = indices if indices is not None else np.arange(len(self.zarr['ecal']))
        self.n_batches = len(self.indices) // self.batch_size

    def __len__(self):
        return self.n_batches
        # return 8

    def __getitem__(self, batch_idx):

        start = batch_idx * self.batch_size
        end = start + self.batch_size
        batch_indices = self.indices[start:end]

        ecal_vox = self.zarr['ecal'][batch_indices]
        hcal_vox = self.zarr['hcal'][batch_indices]
        energy_ecal = self.zarr['energy_ecal'][batch_indices]
        energy_hcal = self.zarr['energy_hcal'][batch_indices]
        fem_ecal = self.zarr['fem_ecal'][batch_indices]
        fem_hcal = self.zarr['fem_hcal'][batch_indices]

        return {
            'ecal': torch.tensor(ecal_vox, dtype=torch.float32).unsqueeze(1),
            'hcal': torch.tensor(hcal_vox, dtype=torch.float32).unsqueeze(1),
            'energy_ecal': torch.tensor(energy_ecal, dtype=torch.float32).unsqueeze(1),
            'energy_hcal': torch.tensor(energy_hcal, dtype=torch.float32).unsqueeze(1),
            'fem_ecal': torch.tensor(fem_ecal, dtype=torch.float32),
            'fem_hcal': torch.tensor(fem_hcal, dtype=torch.float32),
        }


# === 网络模型 ===
class DualBranchFEMNet(nn.Module):
    def __init__(self):
        super().__init__()
        # ECAL branch
        self.ecal_conv = nn.Sequential(
            nn.Conv3d(1, 32, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool3d(2),
            nn.Conv3d(32, 64, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool3d(2),
            nn.Flatten(),
        )
        # HCAL branch
        self.hcal_conv = nn.Sequential(
            nn.Conv3d(1, 32, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool3d(2),
            nn.Conv3d(32, 64, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool3d(2),
            nn.Flatten(),
        )

        dummy_ecal = torch.zeros(1, 1, 40, 30, 40)
        dummy_hcal = torch.zeros(1, 1, 40, 48, 40)
        ecal_feat_size = self.ecal_conv(dummy_ecal).shape[1]
        hcal_feat_size = self.hcal_conv(dummy_hcal).shape[1]

        self.energy_fc = nn.Sequential(
            nn.Linear(2, 32),
            nn.ReLU(),
        )

        fusion_size = ecal_feat_size + hcal_feat_size + 32
        self.fusion_fc = nn.Sequential(
            nn.Linear(fusion_size, 128),
            nn.ReLU(),
            # nn.Dropout(0.3),
            nn.Linear(128, 64),
            nn.ReLU(),
        )
        
        self.fem_ecal_head = nn.Linear(64, 1)
        self.fem_hcal_head = nn.Linear(64, 1)

    def forward(self, ecal_voxel, hcal_voxel, energy_ecal, energy_hcal):
        ecal_feat = self.ecal_conv(ecal_voxel)
        hcal_feat = self.hcal_conv(hcal_voxel)
        energy_feat = self.energy_fc(torch.cat([energy_ecal, energy_hcal], dim=1))

        fused = torch.cat([ecal_feat, hcal_feat, energy_feat], dim=1)
        x = self.fusion_fc(fused)

        fem_ecal_pred = self.fem_ecal_head(x).squeeze(1)
        fem_hcal_pred = self.fem_hcal_head(x).squeeze(1)
        return fem_ecal_pred, fem_hcal_pred


# === 训练和验证主程序 ===
def train_one_epoch(model, dataloader, optimizer, criterion, device, epoch=0, max_batches=None):
    model.train()
    running_loss = 0.0
    num_samples = 0
    batch_logs = []

    for batch_idx, batch in enumerate(dataloader):
        # t0 = time.time()

        ecal = batch['ecal'].squeeze(0).to(device)
        hcal = batch['hcal'].squeeze(0).to(device)
        energy_ecal = batch['energy_ecal'].squeeze(0).to(device)
        energy_hcal = batch['energy_hcal'].squeeze(0).to(device)
        fem_ecal_true = batch['fem_ecal'].squeeze(0).to(device)
        fem_hcal_true = batch['fem_hcal'].squeeze(0).to(device)

        # print(f"model device={next(model.parameters()).device}")
        # print(f"ecal device={ecal.device}")
        # t1 = time.time()

        optimizer.zero_grad()
        fem_ecal_pred, fem_hcal_pred = model(ecal, hcal, energy_ecal, energy_hcal)

        loss_ecal = criterion(fem_ecal_pred, fem_ecal_true)
        loss_hcal = criterion(fem_hcal_pred, fem_hcal_true)
        loss = loss_ecal + loss_hcal

        loss.backward()
        optimizer.step()

        # t2 = time.time()

        num_samples += ecal.size(0)
        running_loss += loss.item() * ecal.size(0)

        batch_logs.append({
            "epoch": epoch,
            "batch_idx": batch_idx,
            "loss": loss.item(),
            # "data_to_gpu_time": round(t1 - t0, 4),
            # "model_step_time": round(t2 - t1, 4),
            # "total_time": round(t2 - t0, 4),
        })

        if max_batches is not None and batch_idx >= max_batches - 1:
            break
    return running_loss / num_samples, batch_logs



def validate(model, dataloader, criterion, device, epoch):
    model.eval()
    val_loss = 0.0
    num_samples = 0
    batch_logs = []

    with torch.no_grad():
        for batch_idx, batch in enumerate(dataloader):
            # t0 = time.time()

            ecal = batch['ecal'].squeeze(0).to(device)
            hcal = batch['hcal'].squeeze(0).to(device)
            energy_ecal = batch['energy_ecal'].squeeze(0).to(device)
            energy_hcal = batch['energy_hcal'].squeeze(0).to(device)
            fem_ecal_true = batch['fem_ecal'].squeeze(0).to(device)
            fem_hcal_true = batch['fem_hcal'].squeeze(0).to(device)

            # t1 = time.time()

            fem_ecal_pred, fem_hcal_pred = model(ecal, hcal, energy_ecal, energy_hcal)
            loss_ecal = criterion(fem_ecal_pred, fem_ecal_true)
            loss_hcal = criterion(fem_hcal_pred, fem_hcal_true)
            loss = loss_ecal + loss_hcal

            # t2 = time.time()

            batch_size = ecal.size(0)
            val_loss += loss.item() * batch_size
            num_samples += batch_size

            batch_logs.append({
                "epoch": epoch,
                "batch_idx": batch_idx,
                "loss": loss.item(),
                # "data_to_gpu_time": round(t1 - t0, 4),
                # "model_step_time": round(t2 - t1, 4),
                # "total_time": round(t2 - t0, 4),
            })

    return val_loss / num_samples, batch_logs



def save_loss_plot(train_losses, val_losses, out_dir):
    plt.figure()
    plt.plot(train_losses, label='Train Loss')
    plt.plot(val_losses, label='Val Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(out_dir, "loss_vs_epoch.png"))
    plt.close()

def save_json_log(data, filename):
    with open(filename, "w") as f:
        json.dump(data, f, indent=2)


if __name__ == "__main__":
    start_time_overall = time.time()
    # 训练配置
    zarr_path = "/gpfs/workdir/xiax/dualdata/fem2_05MIP_data/voxelized_data.zarr" 
    # device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if torch.cuda.is_available():
        device = torch.device("cuda")
        print(f"✅ CUDA is available")
        print(f"🖥️ Device index: {torch.cuda.current_device()}")
        print(f"📟 Device name: {torch.cuda.get_device_name(torch.cuda.current_device())}")
        # print(f"🧠 Memory allocated: {torch.cuda.memory_allocated() / 1024**2:.2f} MB")
    else:
        print("❌ CUDA is not available")

    output_dir = "./output"
    os.makedirs(output_dir, exist_ok=True)

    zarr_data = zarr.open(zarr_path, mode='r')
    n_total = len(zarr_data['ecal'])

    n_train = int(n_total * 0.6)
    n_val = int(n_total * 0.2)

    train_indices = np.arange(0, n_train)
    val_indices = np.arange(n_train, n_train + n_val)
    batch_size = 512

    train_ds = CaloDataset(zarr_path, batch_size=batch_size, indices=train_indices)
    val_ds = CaloDataset(zarr_path, batch_size=batch_size, indices=val_indices)
    
    train_loader = DataLoader(train_ds, batch_size=1, shuffle=False, num_workers=8, pin_memory=True, persistent_workers=True)
    val_loader = DataLoader(val_ds, batch_size=1, shuffle=False, num_workers=8, pin_memory=True, persistent_workers=True)

    model = DualBranchFEMNet().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    criterion = nn.MSELoss()

    best_val_loss = float('inf')
    train_losses = []
    val_losses = []

    epochs = 30
    for epoch in range(epochs):
        start_time = time.time()
        train_loss, train_log = train_one_epoch(model, train_loader, optimizer, criterion, device, epoch=epoch, max_batches=None)
        val_loss, val_log = validate(model, val_loader, criterion, device, epoch=epoch)
        train_losses.append(train_loss)
        val_losses.append(val_loss)

        print(f"Epoch {epoch+1}/{epochs} - Train Loss: {train_loss:.5f}, Val Loss: {val_loss:.5f}, Time: {time.time()-start_time:.1f}s")
        save_json_log(train_log, os.path.join(output_dir, f"epoch_{epoch+1}_train_batch_logs.json"))
        save_json_log(val_log, os.path.join(output_dir, f"epoch_{epoch+1}_val_batch_logs.json"))

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save(model.state_dict(), os.path.join(output_dir, "best_model.pth"))
            print("Best model saved")

    # 保存最终模型和loss曲线
    torch.save(model.state_dict(), os.path.join(output_dir, "final_model.pth"))
    save_loss_plot(train_losses, val_losses, output_dir)
    save_json_log({"train_losses": train_losses, "val_losses": val_losses}, os.path.join(output_dir, "loss_log.json"))


    print(f"Total training time: {(time.time() - start_time_overall) / 60:.2f} minutes")
    print("Training complete")
