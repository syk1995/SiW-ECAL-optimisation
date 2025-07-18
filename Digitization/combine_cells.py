# %%

import os
import ROOT
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
plt.style.use('/home/llr/ilc/shi/code/Plot_style/cepc.mplstyle')
fig_size = (6, 6)
Cell_X_No=40
Cell_Y_No=40
Layer_No=80
Cell_X=5.0
Cell_Y=5.0
Si_Z=0.15
Cell_Z=2.25


def decode_volid(volid):
    volid = int(volid)
    calolayer = volid & 0x7F                 # bits 0–6
    abslayer  = (volid >> 7) & 0x1           # bit 7
    cellid    = (volid >> 8) & 0x1FFF        # bits 8–20
    return calolayer, abslayer, cellid

def decode_indices(cellid):
    index_z = cellid // 1600
    index_y = (cellid % 1600) // 40
    index_x = cellid % 40
    return index_x, index_y, index_z
def encode_volid(calolayer, abslayer, index_x, index_y, index_z):
    cellid = index_x + 40 * index_y + 1600 * index_z
    volid = (calolayer & 0x7F) | ((abslayer & 0x1) << 7) | (cellid << 8)
    return volid

# %%
def combine_cells(input_file, output_file, CombineFactor_X, CombineFactor_Y, CombineFactor_Si, CombineFactor_layer):
    #Set up input
    f_in = uproot.open(input_file)
    tree_in = f_in["events"]
    cellID = tree_in["simplecaloRO.cellID"].array(library="np")#int64
    energy = tree_in["simplecaloRO.energy"].array(library="np")#float32
    pos_x  = tree_in["simplecaloRO.position.x"].array(library="np")#float32
    pos_y  = tree_in["simplecaloRO.position.y"].array(library="np")#float32
    pos_z  = tree_in["simplecaloRO.position.z"].array(library="np")#float32
    #Set up output
    f_out = ROOT.TFile(output_file, "RECREATE")
    tree_out = ROOT.TTree("events", "events")
    cellID_vec = ROOT.std.vector('int')()
    energy_vec = ROOT.std.vector('float')()
    pos_x_vec  = ROOT.std.vector('float')()
    pos_y_vec  = ROOT.std.vector('float')()
    pos_z_vec  = ROOT.std.vector('float')()
    tree_out.Branch("simplecaloRO.cellID", cellID_vec)
    tree_out.Branch("simplecaloRO.energy", energy_vec)
    tree_out.Branch("simplecaloRO.position.x", pos_x_vec)
    tree_out.Branch("simplecaloRO.position.y", pos_y_vec)
    tree_out.Branch("simplecaloRO.position.z", pos_z_vec)
    Combined_X = Cell_X*CombineFactor_X#mm
    Combined_Y = Cell_Y*CombineFactor_Y#mm
    Combined_Si = Si_Z*CombineFactor_Si#mm
    Combined_Z = Cell_Z*CombineFactor_layer#mm
    Combined_layerNo = int(Layer_No / CombineFactor_layer)
    Combined_X_No = int(Cell_X_No / CombineFactor_X)
    Combined_Y_No = int(Cell_Y_No / CombineFactor_Y)
    Combined_Hits = np.zeros((Combined_X_No, Combined_Y_No, Combined_layerNo), dtype=np.float32)
    print("Combined Parameters")
    print("  Combined_X:", Combined_X, "mm", "Combined_X_No:", Combined_X_No)
    print("  Combined_Y:", Combined_Y, "mm", "Combined_Y_No:", Combined_Y_No)
    print("  Combined_Si:", Combined_Si, "mm")
    print("  Combined_Z:", Combined_Z, "mm", "Combined_layerNo:", Combined_layerNo)
    #Loop input and combine cells
    for i_event in range(len(cellID)):
    #for i_event in range(1):    
        if i_event % 10 == 0:
            print(f"Processing event {i_event}/{len(cellID)}")
        cellID_vec.clear()
        energy_vec.clear()
        pos_x_vec.clear()
        pos_y_vec.clear()
        pos_z_vec.clear()
        Combined_Hits.fill(0)
        for i_hit in range(len(cellID[i_event])):
            calolayer, abslayer, cellid = decode_volid(cellID[i_event][i_hit])
            index_x, index_y, index_z = decode_indices(cellid)
            #Be careful the layer is 1-80, we keep this rule for compatibility 
            if index_z < CombineFactor_Si:
                Combined_Hits[index_x // CombineFactor_X, index_y // CombineFactor_Y, calolayer // CombineFactor_layer - 1] += energy[i_event][i_hit]
            #if pos_x[i_event][i_hit] ==12.5:
                #print("Tag hit at (", index_x, index_y, calolayer, index_z,") with energy", energy[i_event][i_hit], "GeV")
        for i_layer in range(Combined_layerNo):
            for i_x in range(Combined_X_No):
                for i_y in range(Combined_Y_No):
                    if Combined_Hits[i_x, i_y, i_layer] > 0:
                        cellID_vec.push_back(encode_volid(i_layer+1, 1, i_x,i_y, CombineFactor_Si-1))
                        energy_vec.push_back(Combined_Hits[i_x, i_y, i_layer])
                        pos_x_vec.push_back((i_x - Combined_X_No / 2 + 0.5) * Combined_X)
                        pos_y_vec.push_back((i_y - Combined_Y_No / 2 + 0.5) * Combined_Y)
                        pos_z_vec.push_back( 1.575 + i_layer * Combined_Z)# 1.575 mm is position Z for layer 0 sublayer 0
                        #if ((i_x - Combined_X_No / 2 + 0.5) * Combined_X) == 12.5:
                            #print("Tag hit at (", i_x, i_y, i_layer, ") with energy", Combined_Hits[i_x, i_y, i_layer], "GeV")
        tree_out.Fill()
    f_in.close()    
    tree_out.Write()
    f_out.Close()

# %%
def main():

    CombineFactor_X = 1
    CombineFactor_Y = 1
    CombineFactor_Si = 5
    CombineFactor_layer=1
    Combined_X = Cell_X*CombineFactor_X#mm
    Combined_Y = Cell_Y*CombineFactor_Y#mm
    Combined_Si = Si_Z*CombineFactor_Si#mm
    Combined_Z = Cell_Z*CombineFactor_layer#mm
    Combined_layerNo = int(Layer_No / CombineFactor_layer)
    Energy = 100.0#GeV
    data_path = "/data_ilc/flc/shi/SiWECAL-Prototype/Simu2025-06/CONF0/mu-"
    input_file = f"{data_path}/MC/{Energy}GeV.root"
    output_file = f"{data_path}/Merged_X{Combined_X}mm_Y{Combined_Y}mm_Si{Combined_Si}mm_layer{Combined_layerNo}/{Energy}GeV.root"
    output_dir = os.path.dirname(output_file)
    os.makedirs(output_dir, exist_ok=True)
    combine_cells(input_file, output_file, CombineFactor_X, CombineFactor_Y, CombineFactor_Si, CombineFactor_layer)


if __name__ == "__main__":
    main()



