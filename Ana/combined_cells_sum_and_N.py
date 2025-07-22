import numpy as np
import ROOT
import matplotlib.pyplot as plt
import math
from scipy.stats import gamma as gamma_dist
from scipy.stats import lognorm as log_normal_dist
from scipy.stats import norm as gaussian_dist
import os
import awkward as ak
import pandas as pd
import uproot
from collections import defaultdict

def fit_gamma(x, norm, shape, loc, scale):
    """Gamma distribution function with location parameter."""
    return norm * gamma_dist.pdf(x, shape, loc=loc, scale=scale)

def fit_log_normal(x, p0,p1,p2):
    """Log-normal distribution function."""
    par= [p0, p1, p2]
    return par[0] * log_normal_dist.pdf(x, s=par[1], scale=par[2])
# parameters: par[0] = normalization, par[1] = sigma, par[2] = mu
def fit_gaussian(x, p0,p1,p2):
    """Gaussian distribution function."""
    par = [p0, p1, p2]
# parameters: par[0] = normalization, par[1] = sigma, par[2] = mean
    return par[0] * gaussian_dist.pdf(x, loc=par[2], scale=par[1])

def error_correction_factor(layer):
    X0W = 3.504
    X0Si = 93.70 * 1000  # microns!
    Ecal_WThickness = 0.7 / X0W

    d_SiX0 = np.array([650.0/X0Si]*4 + [500.0/X0Si]*6 + [320.0/X0Si]*5)
    d_WX0  = np.array([6]*8 + [8]*7) * Ecal_WThickness
    f0 = d_SiX0[layer] / (d_SiX0[layer] + d_WX0[layer])
    print(f"Layer {layer}: d_SiX0 = {d_SiX0[layer]}, d_WX0 = {d_WX0[layer]}, f0 = {f0}")
    return 1.0 / f0

def get_file(input_energy):
    import uproot
    path = "/home/llr/ilc/shi/data/SiWECAL-Prototype/Simu2025-06/CONF0/mu-/MC/"
    filename = f"{path}{float(input_energy)}GeV.root"
    file = uproot.open(filename)
    print(f"Opened file: {filename}")
    return file

def get_file_gamma(input_energy):
    import uproot
    path = "/home/llr/ilc/shi/data/SiWECAL-Prototype/Simu2025-06/CONF0/gamma/MC/"
    filename = f"{path}{float(input_energy)}GeV.root"
    file = uproot.open(filename)
    print(f"Opened file: {filename}")
    return file

def decode_volid(volid_array):
    volid_array = int(volid_array)
    calo_layer = volid_array & 0x7F                 # bits 0–6
    abslayer   = (volid_array >> 7) & 0x1           # bit 7
    cellid     = (volid_array >> 8) & 0x1FFF        # bits 8–20
    return calo_layer, abslayer, cellid

def decode_indices(cellid):
    index_z = cellid // 1600
    index_y = (cellid % 1600) // 40
    index_x = cellid % 40
    return index_x, index_y, index_z

def rms90(arr):
    """Calculate the RMS and peak of the densest 90% window in an array."""
    arr_sorted = np.sort(arr)
    n = len(arr)
    min_rms = float('inf')
    best_window = None
    window_size = int(0.9 * n)
    if window_size < 2:
        return np.mean(arr), 0.0  # Not enough data for a 90% window

    for i in range(n - window_size + 1):
        window = arr_sorted[i:i + window_size]
        rms = np.std(window)
        if rms < min_rms:
            min_rms = rms
            best_window = window

    # Histogram for peak estimation
    counts, bin_edges = np.histogram(best_window, bins=50, density=False)
    peak_index = np.argmax(counts)
    peak = 0.5 * (bin_edges[peak_index] + bin_edges[peak_index + 1])  # center of peak bin

    rms90_value = np.std(best_window)
    return [peak, rms90_value]


def interval_quantile(x, quant=0.9):
    '''Calculate the shortest interval that contains the desired quantile'''
    x = np.sort(x)
    # the number of possible starting points
    n_low = int(len(x) * (1. - quant))
    # the number of events contained in the quantil
    n_quant = len(x) - n_low

    # Calculate all distances in one go
    distances = x[-n_low:] - x[:n_low]
    i_start = np.argmin(distances)

    return i_start, i_start + n_quant

def fit90(x):
    x = np.sort(x)
    n10percent = int(round(len(x)*0.1))
    n90percent = len(x) - n10percent
    return n10percent, n90percent
#code from /home/llr/ilc/shi/code/Tool

def convert_obj_array_to_ak(obj_array):
    return ak.Array([ak.Array(x) for x in obj_array])

from collections import defaultdict
import numpy as np

def analyse_distributions(energies, number_of_sublayers, merge_factor, file_func):
    raw_data = {}

    for energy in energies:
        file = file_func(energy)  # Use your existing function
        tree = file["events"]

        cellID = tree["simplecaloRO.cellID"].array(library="np")
        energy_arr = tree["simplecaloRO.energy"].array(library="np")
        pos_x  = tree["simplecaloRO.position.x"].array(library="np")
        pos_y  = tree["simplecaloRO.position.y"].array(library="np")
        pos_z  = tree["simplecaloRO.position.z"].array(library="np")

        calolayer_array = []
        index_x_array = []
        index_y_array = []
        index_z_array = []

        all_pos_z = np.concatenate(pos_z)
        print(f"Size of all_pos_z: {len(all_pos_z)}")
        pos_combined_z = np.unique(all_pos_z)
        pos_combined_z_layer = (pos_combined_z * 4) // 9
        unique_layers, first_indices = np.unique(pos_combined_z_layer, return_index=True)
        first_pos_per_layer = pos_combined_z[first_indices]
            
        total_energy_per_event = [np.sum(event) for event in energy_arr]

        for event_cellID in cellID:
            calo_list, x_list, y_list, z_list = [], [], [], []
            for volID in event_cellID:
                calolayer, abslayer, cellid = decode_volid(volID)
                index_x, index_y, index_z = decode_indices(cellid)
                calo_list.append(calolayer)
                x_list.append(index_x)
                y_list.append(index_y)
                z_list.append(index_z)

            calolayer_array.append(np.array(calo_list, dtype=np.int32))
            index_x_array.append(np.array(x_list, dtype=np.int32))
            index_y_array.append(np.array(y_list, dtype=np.int32))
            index_z_array.append(np.array(z_list, dtype=np.int32))


        sum_E_list = []
        N_hits_list = []
        all_hit_energies = []

        for i in range(len(energy_arr)):
            superlayers = []
            raw_layers = calolayer_array[i]
            for layer in raw_layers:
                if layer in superlayers:
                    continue
                superlayers.append(layer)
                superlayers = sorted(superlayers)


        #    current_layer = raw_layers[i]
        #   if current_layer != prev_layer:
        #        prev_layer = current_layer
        #        logical_index += 1
        #    else:
        #        logical_index += 1
        #    logical_layers[i] = logical_index
        #    sub_count += 1

        #   superlayers = logical_layers

            e = energy_arr[i]
            ix = index_x_array[i]
            iy = index_y_array[i]
            iz = index_z_array[i]

            combined_hits = defaultdict(float)
            for sl, x, y, energy_hit in zip(raw_layers, ix, iy, e):
                super_x = x // merge_factor
                super_y = y // merge_factor
                combined_hits[(sl, super_x, super_y)] += energy_hit
            e_list = list(combined_hits.values())
            if all(hit_energy == 0 for hit_energy in e_list):
                # Skip this event — no valid hits
                continue
            all_hit_energies.extend(e_list)      #hit_energies for all events
            sum_E_list.append(np.sum(e_list))    # sum of all energy hits combined in supercells for this event
            N_hits_list.append(len(e_list))      # number of supercells (hits) for this event

            if i == 0:  # Print for the first event only
                print(f"\n=== Event {i} ===")
                print(f"Raw calolayers: {calolayer_array[i]}")
                print(f"Superlayers: {superlayers}")
                print(f"index_x: {index_x_array[i]}")
                print(f"index_y: {index_y_array[i]}")
                print(f"index_z: {index_z_array[i]}")
                print(f"Energy array: {energy_arr[i]}")
                print(f"position_z: {all_pos_z[:20]}")
                print(f"position of superlayers: {pos_combined_z}")
                print(f"position_z layers: {pos_combined_z_layer}")
                print(f"the x-axis: {superlayers}"),
                print(f"the y-axis: {first_pos_per_layer}"),
                print(f"length of the x-axis: {len(superlayers)}"),
                print(f"length of the y-axis: {len(first_pos_per_layer)}")
                print(f"hit_energies: {e_list}")
                print(f"All hit energies: {sorted(all_hit_energies)[:20]}")
                print(f"raw_layers: {raw_layers}")
                #print(f"total energy per event: {total_energy_per_event}")


                
                print(f"\nAfter merging (supercells):")
                for key, val in combined_hits.items():
                    print(f"Layer {key[0]}, x {key[1]}, y {key[2]} -> E = {val:.4f}")

                print(f"\nSum energy (merged): {np.sum(e_list):.4f}")
                print(f"Number of merged hits (supercells): {len(e_list)}")

        raw_data[energy] = {
            "sum_E_arr": np.array(sum_E_list),
            "N_hits_arr": np.array(N_hits_list),
            "all_hit_energies": np.array(all_hit_energies),
            "superlayers": np.array(superlayers),
            "position_z": np.array(first_pos_per_layer),
            "total_energy_per_event": np.array(total_energy_per_event),
            "e_list": np.array(e_list)
        }

    return raw_data

def plot_fits(sum_E_arr, N_hits_arr, energy, calolayers, position_z, all_hit_energies, e_list, e_list_gamma, energy_gamma, N_hits_arr_gamma, sum_E_arr_gamma):
    plt.figure(figsize=(15, 15))
    plt.suptitle(f"Energy: {energy} GeV of mu-, and {energy_gamma} GeV for gamma", fontsize=16)    
    plt.subplot(3,3,1)
    plt.hist(sum_E_arr*1000, bins=1000, color='blue', alpha=0.7)

    #if energy >= 20:
    #    plt.xlim(0.75 * np.mean(sum_E_arr), 1.25 * np.mean(sum_E_arr))
    plt.xlim(0,1000)
    plt.title("Sum Energy Distribution")
    plt.xlabel("Sum Energy (MeV)")
    plt.ylabel("Number of Events")

    plt.subplot(3,3,2)
    plt.hist(N_hits_arr, bins = 1000, color='green', alpha=0.7)
    plt.title("Number of Hits Distribution")
    plt.ylabel("Number of Events")
    plt.xlabel("Number of Hits")
    #plt.xlim(0, max(N_hits_arr) + 1)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.xlim(60,200)
    #plt.savefig(f"/home/llr/ilc/oprea/data/energy_plot_{energy}GeV.png")
    #plt.show()

    plt.subplot(3,3,3)
    plt.hist2d(position_z, calolayers, bins=80)
    plt.colorbar(label='Counts')
    plt.xlabel("Superlayer Position")
    plt.ylabel("Index of Superlayer")
    plt.title("Position Z vs Index Z")

    plt.subplot(3,3,4)
    plt.hist(all_hit_energies*100000, bins=1000, color='purple', alpha=0.7)
    plt.title("Total Energy per Event")
    plt.xlim(0,6000)
    plt.ylim(0,150000)
    plt.xlabel("Total Energy (MeV)")
    plt.ylabel("Number of Events")


    plt.subplot(3,3,6)
    plt.hist(e_list*1000, bins=100, color='orange', alpha=0.7)
    plt.title("Hit Energies Distribution")
    plt.xlabel("Hit Energy (MeV)")
    plt.ylabel("Number of Hits")


    plt.subplot(3,3,5)
    plt.plot(e_list*1000, 'o', markersize=2, color='red', alpha=0.5, label='mu-')
    plt.plot(e_list_gamma*1000, 'o', markersize=2, color='blue', alpha=0.5, label='gamma')
    plt.title("Hit Energies for Mu- and Gammas")
    plt.ylabel("Number of Hits")
    plt.xlabel("Hit Energy (MeV)")

    plt.subplot(3,3,7)
    plt.hist(N_hits_arr_gamma, bins='auto', color='cyan', alpha=0.7)
    plt.title("Number of Hits Distribution for Gamma")
    plt.xlabel("Number of Hits")
    plt.ylabel("Number of Events")

    plt.subplot(3,3,8)
    plt.hist(sum_E_arr_gamma*1000, bins='auto', color='magenta', alpha=0.7)
    plt.title("Sum Energy Distribution for Gamma")
    plt.xlabel("Sum Energy (MeV)")
    plt.ylabel("Number of Events")
    #plt.savefig(f"/mnt/d/plots/combined_cells_sum_and_N_mu-_{energy}GeV_gamma_{energy_gamma}GeV.png") if not os.path.exists("/mnt/d/plots/combined_cells_sum_and_N_mu-_{energy}GeV_gamma_{energy_gamma}GeV.png") else None
    #change the path to the one needed in the server
    #plt.show()

def write_root_file(raw_data, raw_data_gamma, filename):
    """
    ROOT file structure written by this function:

    root_file.root
    ├── mu-/
    │   ├── energy       (Tree: 1 branch, energy labels, one per event)
    │   ├── sum_E        (Tree: 1 branch, sum_E per event)
    │   ├── N_hits       (Tree: 1 branch, N_hits per event)
    │   ├── hit_energy   (Tree: 1 branch, one entry per hit)
    │   ├── position_z   (Tree: 1 branch, one entry per hit)
    │   └── superlayer   (Tree: 1 branch, one entry per hit)
    └── gamma/
        ├── energy
        ├── sum_E
        ├── N_hits
        ├── hit_energy
        ├── position_z
        └── superlayer
    """

    with uproot.recreate(filename) as f:
        for particle_type, dataset in {"mu-": raw_data, "gamma": raw_data_gamma}.items():
            f.mkdir(particle_type)

            # Collect per-event arrays
            all_sum_E = []
            all_N_hits = []
            all_energy_labels_events = []

            for energy, data in dataset.items():
                n_events = len(data["sum_E_arr"])
                all_sum_E.extend(data["sum_E_arr"])
                all_N_hits.extend(data["N_hits_arr"])
                all_energy_labels_events.extend([energy] * n_events)

            # Write separate trees for event-level quantities
            f[f"{particle_type}/energy"] = {"energy": np.array(all_energy_labels_events, dtype="float32")}
            f[f"{particle_type}/sum_E"] = {"sum_E": np.array(all_sum_E, dtype="float32")}
            f[f"{particle_type}/N_hits"] = {"N_hits": np.array(all_N_hits, dtype="int32")}

            # Collect per-hit arrays
            hit_energies = []
            positions_z = []
            superlayers = []

            for energy, data in dataset.items():
                hit_energies.extend(data["all_hit_energies"])
                positions_z.extend(data["position_z"])
                superlayers.extend(data["superlayers"])

            # Write separate trees for hit-level quantities
            f[f"{particle_type}/hit_energy"] = {"hit_energy": np.array(hit_energies, dtype="float32")}
            f[f"{particle_type}/position_z"] = {"position_z": np.array(positions_z, dtype="float32")}
            f[f"{particle_type}/superlayer"] = {"superlayer": np.array(superlayers, dtype="int32")}
            print(f"Written {particle_type} data to {filename}")

energies1 = [0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
energies = [200.0]

# Here you can change the sublayer and merge factor parameters
raw_data = analyse_distributions(energies, number_of_sublayers=5, merge_factor=2, file_func=get_file)
raw_data_gamma = analyse_distributions(energies1, number_of_sublayers=5, merge_factor=2, file_func=get_file_gamma)

for energy_to_plot in energies1:
    plot_fits(
        raw_data[energies[0]]["sum_E_arr"],
        raw_data[energies[0]]["N_hits_arr"],
        energies[0],
        raw_data[energies[0]]["superlayers"],
        raw_data[energies[0]]["position_z"],
        raw_data[energies[0]]["all_hit_energies"],
        raw_data[energies[0]]["e_list"],
        raw_data_gamma[energy_to_plot]["e_list"],
        energy_to_plot,
        raw_data_gamma[energy_to_plot]["N_hits_arr"],
        raw_data_gamma[energy_to_plot]["sum_E_arr"]
    )
    write_root_file(raw_data, raw_data_gamma, f"/home/llr/ilc/shi/data/SiWECAL-Prototype/Simu2025-06/CONF0/combined_cells_sum_and_N_mu-_{energies[0]}GeV_gamma_{energy_to_plot}GeV.root")
    #change the path to the one needed in the server

