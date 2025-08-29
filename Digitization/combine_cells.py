# %%

import os
import sys
import ROOT
import uproot
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import argparse
import time
plt.style.use('/home/llr/ilc/shi/code/Plot_style/cepc.mplstyle')
fig_size = (6, 6)
Cell_X_No=40
Cell_Y_No=40
Layer_No=80
Cell_X=5.0
Cell_Y=5.0
Si_Z=0.15
Cell_Z=2.25
Start_Z=1.575
MIP=(0.0410,0.0861,0.1328,0.1803,0.2282)#MeV
def decode_volid_and_indices(volid_array):
    volid_array = volid_array.astype(np.int64)
    calolayer = volid_array & 0x7F
    abslayer  = (volid_array >> 7) & 0x1
    cellid    = (volid_array >> 8) & 0x1FFF
    index_z = cellid // 1600
    index_y = (cellid % 1600) // 40
    index_x = cellid % 40
    return calolayer, abslayer, index_x, index_y, index_z
def encode_volid_array(calolayer, abslayer, index_x, index_y, index_z):
    cellid = index_x + 40 * index_y + 1600 * index_z
    volid = (calolayer & 0x7F) | ((abslayer & 0x1) << 7) | (cellid << 8)
    return volid.astype(np.int64)

# %%
vectors_out = {
    "simplecaloRO.cellID": ROOT.std.vector('long')(),
    "simplecaloRO.energy": ROOT.std.vector('float')(),
    "simplecaloRO.position.x": ROOT.std.vector('float')(),
    "simplecaloRO.position.y": ROOT.std.vector('float')(),
    "simplecaloRO.position.z": ROOT.std.vector('float')(),
    #"MCParticles.momentum.x": ROOT.std.vector('float')(),
    #"MCParticles.momentum.y": ROOT.std.vector('float')(),
    #"MCParticles.momentum.z": ROOT.std.vector('float')(),
    "MCParticles.p0": ROOT.std.vector('float')()
}

def clear_vectors():
    for v in vectors_out.values():
        v.clear()

def combine_cells(input_file, output_file, CombineFactor_X, CombineFactor_Y, CombineFactor_Si, CombineFactor_layer, Absorber_layer):
    #Set up input
    try:
        with uproot.open(input_file) as f_in:
            tree_in = f_in["events"]
            print(f"[Info] Loaded tree 'events' from {input_file}")
            cellID = tree_in["simplecaloRO.cellID"].array(library="np")#int64
            energy = tree_in["simplecaloRO.energy"].array(library="np")#float32
            pos_x  = tree_in["simplecaloRO.position.x"].array(library="np")#float32
            pos_y  = tree_in["simplecaloRO.position.y"].array(library="np")#float32
            pos_z  = tree_in["simplecaloRO.position.z"].array(library="np")#float32
            MCP_px = tree_in["MCParticles.momentum.x"].array(library="np")
            MCP_py = tree_in["MCParticles.momentum.y"].array(library="np")
            MCP_pz = tree_in["MCParticles.momentum.z"].array(library="np")
    except Exception:
        sys.exit(1)
    #Set up output
    try:
        print(f"[INFO] RECREATE {output_file}")
        f_out = ROOT.TFile(output_file, "RECREATE")
        tree_out = ROOT.TTree("events", "events")
        for name, vec in vectors_out.items():
            tree_out.Branch(name, vec)
    except Exception:
        sys.exit(1)
    #Calo hits
    Combined_Cell_X = Cell_X*CombineFactor_X#mm
    Combined_Cell_Y = Cell_Y*CombineFactor_Y#mm
    Combined_Si = round(Si_Z*CombineFactor_Si, 2)#mm
    Combined_Z = Cell_Z*CombineFactor_layer#mm
    Combined_layerNo = int(Absorber_layer / CombineFactor_layer)
    Combined_X_No = int(Cell_X_No / CombineFactor_X)
    Combined_Y_No = int(Cell_Y_No / CombineFactor_Y)
    print("Combined Parameters")
    print("  Combined_X:", Combined_Cell_X, "mm", "Combined_X_No:", Combined_X_No)
    print("  Combined_Y:", Combined_Cell_Y, "mm", "Combined_Y_No:", Combined_Y_No)
    print("  Combined_Si:", Combined_Si, "mm", "Combined_Factor_Si:", CombineFactor_Si, "Si_Z:", Si_Z)
    print("  Combined_Z:", Combined_Z, "mm", "Selected_LayerNo:", Absorber_layer, "Combined_layerNo:", Combined_layerNo)
    #Loop input and combine cells
    for i_event in range(len(cellID)):
    #for i_event in range(2):    
        clear_vectors()
        if i_event % 1000 == 0:
            print("Processing event:", i_event)
        calolayer, abslayer, index_x, index_y, index_z = decode_volid_and_indices(cellID[i_event])
        energy_i= energy[i_event]
        mask = pos_x[i_event] >10
        # print("  before combining (pos_x == 12.5):")
        # print("    calolayer:", calolayer[mask])
        # print("    abslayer:", abslayer[mask])
        # print("    index_x:", index_x[mask])
        # print("    index_y:", index_y[mask])
        # print("    index_z:", index_z[mask])
        # print("    energy_i:", energy_i[mask])
        # print("    pos_x:", pos_x[i_event][mask])
        # print("    pos_y:", pos_y[i_event][mask])
        # print("    pos_z:", pos_z[i_event][mask])
        #Select Si
        #attention: input calolayer starts from 1, but we output from 0
        mask = (index_z < CombineFactor_Si) & ((calolayer-1)<Absorber_layer) & ((calolayer-1) % CombineFactor_layer == 0)
        calolayer, abslayer, index_x, index_y, index_z, energy_i = calolayer[mask], abslayer[mask], index_x[mask], index_y[mask], index_z[mask], energy_i[mask]
        #Combine X,Y,layer
        Combined_index_x = index_x // CombineFactor_X
        Combined_index_y = index_y // CombineFactor_Y
        Combined_index_z = np.full_like(Combined_index_x, CombineFactor_Si-1, dtype=np.int64)
        Combined_calolayer = (calolayer-1) // CombineFactor_layer
        Combined_abslayer = abslayer
        Combined_cellID = encode_volid_array(Combined_calolayer, Combined_abslayer, Combined_index_x, Combined_index_y, Combined_index_z)
        Combined_calolayer, Combined_abslayer, Combined_index_x, Combined_index_y, Combined_index_z = decode_volid_and_indices(Combined_cellID)
        Combined_cellID, unique_indices = np.unique(Combined_cellID, return_inverse=True)
        Combined_energy = np.zeros_like(Combined_cellID, dtype=np.float32)
        np.add.at(Combined_energy, unique_indices, energy_i)
        Combined_calolayer, Combined_abslayer, Combined_index_x, Combined_index_y, Combined_index_z = decode_volid_and_indices(Combined_cellID)
        Combined_pos_x = (Combined_index_x - Combined_X_No / 2.0 + 0.5) * Combined_Cell_X
        Combined_pos_y = (Combined_index_y - Combined_Y_No / 2.0 + 0.5) * Combined_Cell_Y
        Combined_pos_z = Start_Z + Combined_calolayer * Combined_Z
        mask = Combined_pos_x >10
        # print("  after combining :")
        # print("    calolayer:", Combined_calolayer[mask])
        # print("    abslayer:", Combined_abslayer[mask])
        # print("    index_x:", Combined_index_x[mask])
        # print("    index_y:", Combined_index_y[mask])
        # print("    index_z:", Combined_index_z[mask])
        # print("    energy_i:", Combined_energy[mask])
        # print("    pos_x:", Combined_pos_x[mask])
        # print("    pos_y:", Combined_pos_y[mask])
        # print("    pos_z:", Combined_pos_z[mask])
        #Store output
        #MCP p0
        px0 = MCP_px[i_event][0]
        py0 = MCP_py[i_event][0]
        pz0 = MCP_pz[i_event][0]
        p0 = np.sqrt(px0**2 + py0**2 + pz0**2)

        vectors_out["simplecaloRO.cellID"].assign(Combined_cellID.astype(np.int64))
        vectors_out["simplecaloRO.energy"].assign(Combined_energy.astype(np.float32))
        vectors_out["simplecaloRO.position.x"].assign(Combined_pos_x.astype(np.float32))
        vectors_out["simplecaloRO.position.y"].assign(Combined_pos_y.astype(np.float32))
        vectors_out["simplecaloRO.position.z"].assign(Combined_pos_z.astype(np.float32))
        #vectors_out["MCParticles.momentum.x"].assign(MCP_px[i_event])
        #vectors_out["MCParticles.momentum.y"].assign(MCP_py[i_event])
        #vectors_out["MCParticles.momentum.z"].assign(MCP_pz[i_event])
        vectors_out["MCParticles.p0"].push_back(p0)
        tree_out.Fill()
    #loop event end

    f_in.close()
    f_out.cd()
    tree_out.Write()
    f_out.Close()
    print(f"[Info] Written combined file to {output_file}")
# %%
def main():
    print("[Info] Starting cell combination...")
    start_time = time.time()
    parser = argparse.ArgumentParser(description="Combine cells and write merged ROOT file.")
    parser.add_argument("--CombineFactor_X", type=int, default=1)
    parser.add_argument("--CombineFactor_Y", type=int, default=1)
    parser.add_argument("--CombineFactor_Si", type=int, default=5)
    parser.add_argument("--CombineFactor_layer", type=int, default=1)
    parser.add_argument("--Absorber_layer", type=int, default=80)
    parser.add_argument("--Energy", type=float, default=-1)#-1 for all files in the input path
    parser.add_argument("--input_file", type=str, default="/data_ilc/flc/shi/SiWECAL-Prototype/Simu2025-06/CONF0/mu-")
    parser.add_argument("--output_file", type=str, default="/data_ilc/flc/shi/SiWECAL-Prototype/Simu2025-06/CONF0/mu-")
    args = parser.parse_args()

    Combined_X = Cell_X * args.CombineFactor_X  # mm
    Combined_Y = Cell_Y * args.CombineFactor_Y  # mm
    Combined_Si = round(Si_Z * args.CombineFactor_Si,2)  # mm
    Combined_Z = Cell_Z * args.CombineFactor_layer  # mm
    Absorber_layer = args.Absorber_layer
    Combined_layerNo = int(Absorber_layer / args.CombineFactor_layer)
    
    combine_cells(
        args.input_file,
        args.output_file,
        args.CombineFactor_X,
        args.CombineFactor_Y,
        args.CombineFactor_Si,
        args.CombineFactor_layer,
        Absorber_layer,
    )
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"[Info] Cell combination completed in {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()



