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
import numpy as np
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
    path = "/grid_mnt/vol_home/llr/ilc/shi/data/SiWECAL-Prototype/Simu2025-06/CONF0/mu-/MC/"
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

def analyse_distributions(energies, number_of_sublayers, merge_factor):
    raw_data = {}

    for energy in energies:
        file = get_file(energy)  # Use your existing function
        tree = file["events"]

        cellID = tree["simplecaloRO.cellID"].array(library="np")
        energy_arr = tree["simplecaloRO.energy"].array(library="np")

        calolayer_array = []
        index_x_array = []
        index_y_array = []
        index_z_array = []

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


        for i in range(len(energy_arr)):
            raw_layers = calolayer_array[i]
            logical_layers = np.zeros_like(raw_layers)
            logical_index = 0
            sub_count = 0
            prev_layer = raw_layers[0]

            current_layer = raw_layers[i]
            if current_layer != prev_layer:
                prev_layer = current_layer
                sub_count = 0
                logical_index += 1
            elif sub_count == number_of_sublayers:
                sub_count = 0
                logical_index += 1
            logical_layers[i] = logical_index
            sub_count += 1

            superlayers = logical_layers

            e = energy_arr[i]
            ix = index_x_array[i]
            iy = index_y_array[i]
            iz = index_z_array[i]

            combined_hits = defaultdict(float)
            for sl, x, y, energy_hit in zip(superlayers, ix, iy, e):
                super_x = x // merge_factor
                super_y = y // merge_factor
                combined_hits[(sl, super_x, super_y)] += energy_hit
            e_list = list(combined_hits.values())

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
                
                print(f"\nAfter merging (supercells):")
                for key, val in combined_hits.items():
                    print(f"Layer {key[0]}, x {key[1]}, y {key[2]} -> E = {val:.4f}")

                print(f"\nSum energy (merged): {np.sum(e_list):.4f}")
                print(f"Number of merged hits (supercells): {len(e_list)}")

        raw_data[energy] = {
            "sum_E_arr": np.array(sum_E_list),
            "N_hits_arr": np.array(N_hits_list)
        }

    return raw_data

def plot_fits(sum_E_arr, N_hits_arr, energy):
    plt.figure(figsize=(15, 5))
    plt.suptitle(f"Energy: {energy} GeV", fontsize=16)    
    plt.subplot(1,3,1)
    plt.hist(sum_E_arr*1000, bins=1000, color='blue', alpha=0.7)

    #if energy >= 20:
    #    plt.xlim(0.75 * np.mean(sum_E_arr), 1.25 * np.mean(sum_E_arr))
    plt.xlim(0,1000)
    plt.title("Sum Energy Distribution")
    plt.xlabel("Sum Energy")
    plt.ylabel("Number of Events")

    plt.subplot(1,3,2)
    plt.hist(N_hits_arr, bins = 1000, color='green', alpha=0.7)
    plt.title("Number of Hits Distribution")
    plt.ylabel("Number of Events")
    plt.xlabel("Number of Hits")
    #plt.xlim(0, max(N_hits_arr) + 1)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.xlim(60,200)
    #plt.savefig(f"/home/llr/ilc/oprea/data/energy_plot_{energy}GeV.png")
    plt.show()

energies1 = [0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
energies = [200.0]

# Here you can change the sublayer and merge factor parameters
raw_data = analyse_distributions(energies, number_of_sublayers=5, merge_factor=2)
# Plotting the fits for a specific energy
#energy_to_plot = 0.15  # Change this to the desired energy for plotting
for energy_to_plot in energies:
    plot_fits(
        raw_data[energy_to_plot]["sum_E_arr"],
        raw_data[energy_to_plot]["N_hits_arr"],
        energy_to_plot
    )
    

