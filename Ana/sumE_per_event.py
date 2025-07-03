import numpy as np
import ROOT
import matplotlib.pyplot as plt
import math
from scipy.stats import gamma as gamma_dist
from scipy.stats import lognorm as log_normal_dist
from scipy.stats import norm as gaussian_dist
import os
import awkward as ak

def fit_gamma(x, p0,p1,p2):
    """Gamma distribution function."""
# parameters: par[0] = normalization, par[1] = shape, par[2] = scale
    par = [p0, p1, p2]
    return par[0] * gamma_dist.pdf(x, par[1], scale=par[2])
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
    path = "/grid_mnt/vol_home/llr/ilc/shi/data/SiWECAL-Prototype/Simu2025-06/CONF0/gamma/MC/"
    filename = f"{path}{float(input_energy)}GeV.root"
    file = uproot.open(filename)
    print(f"Opened file: {filename}")
    return file

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


def analyse_distributions(energies):
    fit_results = {}
    raw_data = {}
    for energy in energies:
        file = get_file(energy)
        tree = file["events"]
        print(f"Processing energy: {energy} GeV")

        hit_energy = tree["simplecaloRO.energy"].array(library = 'ak')
        hit_cellID = tree["simplecaloRO.cellID"].array(library = 'ak')

        sum_E = ak.sum(hit_energy, axis=1)
        N_hits = ak.num(hit_energy, axis=1)

        print(f"Number of events: {len(sum_E)}")
        print(f"Mean energy: {np.mean(sum_E):.4f} GeV")
        print(f"Mean hits: {np.mean(N_hits):.1f}")

        # If cellID is not available, we can still count the number of hits
        if hit_cellID is not None:
            N_hits = ak.num(hit_cellID, axis=1)
        sum_E = ak.to_numpy(sum_E)
        N_hits = ak.to_numpy(N_hits)

        sum_E_arr = np.array(sum_E)
        N_hits_arr = np.array(N_hits)
        print(f"Energy {energy} GeV: Sum_E mean = {np.mean(sum_E_arr)}, N_hits mean = {np.mean(N_hits_arr)}")
        print(f"Energy {energy} GeV: Sum_E std = {np.std(sum_E_arr)}, N_hits std = {np.std(N_hits_arr)}")
        print(f"Energy {energy} GeV: Sum_E_var = {np.var(sum_E_arr)}, N_hits_var = {np.var(N_hits_arr)}")
        print(f"Energy {energy} GeV: Sum_E min = {np.min(sum_E_arr)}, N_hits min = {np.min(N_hits_arr)}")
        print(f"Energy {energy} GeV: Sum_E max = {np.max(sum_E_arr)}, N_hits max = {np.max(N_hits_arr)}")
        print(f"Energy {energy} GeV: Sum_E median = {np.median(sum_E_arr)}, N_hits median = {np.median(N_hits_arr)}")
        print(f"Energy {energy} GeV: Sum_E 90% RMS = {rms90(sum_E_arr)}, N_hits 90% RMS = {rms90(N_hits_arr)}")
        print(f"N_hits_arr = {N_hits_arr}")
        
    #Fitting the distributions
        from scipy.optimize import curve_fit
        number_of_events, bin_edges= np.histogram(sum_E_arr, bins=50, density=False)
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        start_E, end_E = interval_quantile(sum_E_arr, quant=0.9)
        mean_initial_E = np.mean(sum_E_arr[start_E:end_E])
        std_initial_E = np.std(sum_E_arr[start_E:end_E])
        start_N, end_N = interval_quantile(N_hits_arr, quant=0.9)
        mean_N_initial_N = np.mean(N_hits_arr[start_N:end_N])
        std_N_initial_N = np.std(N_hits_arr[start_N:end_N])

        try:
            popt_gamma_sum_E, _ = curve_fit(fit_gamma, bin_centers, number_of_events, p0=[np.max(number_of_events), 2.0, mean_initial_E], maxfev=10000)
            print(f"Gamma fit parameters for sum_E: {popt_gamma_sum_E}")
        except RuntimeError as e:
            print(f"Gamma fit failed for sum_E with error: {e}")
        norm_gamma = popt_gamma_sum_E[0]
    
        # Define a 20% range around gamma normalization
        norm_lower = 0.8 * norm_gamma
        norm_upper = 1.2 * norm_gamma

        # Bounds: norm, sigma, mu
        bounds = (
            [norm_lower, 0, 0],
            [norm_upper, np.inf, np.inf]
        )
        try:
            popt_log_normal_sum_E, _ = curve_fit(fit_log_normal, bin_centers, number_of_events, p0=[norm_gamma, std_initial_E, std_initial_E], bounds = bounds, maxfev=100000)
            print(f"Log-normal fit parameters for sum_E: {popt_log_normal_sum_E}")
        except RuntimeError as e:
            print(f"Log-normal fit failed for sum_E with error: {e}")
        try:
            popt_gaussian_sum_E, _ = curve_fit(fit_gaussian, bin_centers, number_of_events, p0=[np.max(number_of_events), 2.0 * std_initial_E, mean_initial_E], maxfev=10000)
            print(f"Gaussian fit parameters for sum_E: {popt_gaussian_sum_E}")
        except RuntimeError as e:
            print(f"Gaussian fit failed for sum_E with error: {e}")
        sigma_E = popt_gaussian_sum_E[1]
        mu_E = popt_gaussian_sum_E[2]
        mask_2_sigma_E = (bin_centers >= mu_E - 2 * sigma_E) & (bin_centers <= mu_E + 2 * sigma_E)
        x_fit_2_sigma_E = bin_centers[mask_2_sigma_E]
        y_fit_2_sigma_E = number_of_events[mask_2_sigma_E]
        p0_2_sigma_E = [np.max(y_fit_2_sigma_E), sigma_E, mu_E]
        try:
            popt_gaussian_2_sigma_sum_E, _ = curve_fit(fit_gaussian, x_fit_2_sigma_E, y_fit_2_sigma_E, p0=p0_2_sigma_E, maxfev = 10000)
            print(f"Gaussian 2\u03C3 fit parameters for sum_E: {popt_gaussian_2_sigma_sum_E}")
        except RuntimeError as e:
            print(f"Gaussian 2\u03C3 fit failed for sum_E with error: {e}")

        bins_for_low_E=np.arange(min(N_hits_arr), max(N_hits_arr) + 1.5) - 0.5
        number_of_events, bin_edges = np.histogram(N_hits_arr, bins=bins_for_low_E if len(bins_for_low_E)<50 else 50, density=False)
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        try:
            popt_gamma_N_hits, _ = curve_fit(fit_gamma, bin_centers, number_of_events, p0=[np.max(number_of_events), 5.0, mean_N_initial_N], maxfev=10000)
            print(f"Gamma fit parameters for N_hits: {popt_gamma_N_hits}")
        except RuntimeError as e:
            print(f"Gamma fit failed for N_hits with error: {e}")
        try:
            popt_log_normal_N_hits, _ = curve_fit(fit_log_normal, bin_centers, number_of_events, p0=[np.max(number_of_events), std_N_initial_N, mean_N_initial_N], maxfev=10000)
            print(f"Log-normal fit parameters for N_hits: {popt_log_normal_N_hits}")
        except RuntimeError as e:   
            print(f"Log-normal fit failed for N_hits with error: {e}")
        try:
            popt_gaussian_N_hits, _ = curve_fit(fit_gaussian, bin_centers, number_of_events, p0=[np.max(number_of_events), std_N_initial_N, mean_N_initial_N], maxfev=10000)
            print(f"Gaussian fit parameters for N_hits: {popt_gaussian_N_hits}")
        except RuntimeError as e:
            print(f"Gaussian fit failed for N_hits with error: {e}")
        sigma_N = popt_gaussian_N_hits[1]
        mu_N = popt_gaussian_N_hits[2]
        mask_2_sigma_N = (bin_centers >= mu_N - 2 * sigma_N) & (bin_centers <= mu_N + 2 * sigma_N)
        x_fit_2_sigma_N = bin_centers[mask_2_sigma_N]
        y_fit_2_sigma_N = number_of_events[mask_2_sigma_N]
        p0_2_sigma_N = [np.max(y_fit_2_sigma_N), sigma_N, mu_N]
        try:
            popt_gaussian_2_sigma_N_hits, _ = curve_fit(fit_gaussian, x_fit_2_sigma_N, y_fit_2_sigma_N, p0=p0_2_sigma_N, maxfev=10000)
            print(f"Gaussian 2\u03C3 fit parameters for N_hits: {popt_gaussian_2_sigma_N_hits}")
        except RuntimeError as e:
            print(f"Gaussian 2\u03C3 fit failed for N_hits with error: {e}")

        rms90_sum_E = np.array(rms90(sum_E_arr))
        rms90_N_hits = np.array(rms90(N_hits_arr))


        fit_results[energy] = {"gamma_sum_E" : popt_gamma_sum_E,
                               "log_normal_sum_E" : popt_log_normal_sum_E,
                               "gaussian_sum_E" : popt_gaussian_sum_E,
                               "gaussian_2_sigma_sum_E": popt_gaussian_2_sigma_sum_E,
                               "rms90_sum_E": rms90_sum_E,
                               "gamma_N_hits" : popt_gamma_N_hits,
                               "log_normal_N_hits" : popt_log_normal_N_hits,
                               "gaussian_N_hits" : popt_gaussian_N_hits,
                               "gaussian_2_sigma_N_hits": popt_gaussian_2_sigma_N_hits,
                               "rms90_N_hits": rms90_N_hits}
                               
        raw_data[energy] = {
            "sum_E_arr": sum_E_arr,
            "N_hits_arr": N_hits_arr,
        }
    print (f"raw_data energies: {raw_data.keys()}")



  #      print("Number of hits per event:")
  #      Seeing why the number of events is 0 in certain number of hits
  #      from collections import Counter
  #      hit_counts = Counter(N_hits_arr)
  #
  #      min_hit = int(np.min(N_hits_arr))
  #      max_hit = int(np.max(N_hits_arr))
  #      missing_hits = [h for h in range(min_hit, max_hit+1) if hit_counts[h] == 0]
  #      print(f"Missing hits for energy {energy} GeV: {missing_hits}")
  #      print("Hit counts with few events:")
  #      for h in sorted(hit_counts):
  #          print(f"Hits: {h}, Events: {hit_counts[h]}")
  #      
  #      print(f"Just checking something: Event 18: Hits: {hit_counts[18]}")

    return fit_results, raw_data
energies = [0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0]


#lol

def plot_fits(sum_E_arr, N_hits_arr, fit_results, energy):
    fit_results_energy = fit_results[energy]
    popt_gamma_sum_E = fit_results_energy["gamma_sum_E"]
    popt_log_normal_sum_E = fit_results_energy["log_normal_sum_E"]
    popt_gaussian_sum_E = fit_results_energy["gaussian_sum_E"]
    popt_gaussian_2_sigma_sum_E = fit_results_energy["gaussian_2_sigma_sum_E"]
    popt_gamma_N_hits = fit_results_energy["gamma_N_hits"]
    popt_log_normal_N_hits = fit_results_energy["log_normal_N_hits"]
    popt_gaussian_N_hits = fit_results_energy["gaussian_N_hits"]
    popt_gaussian_2_sigma_N_hits = fit_results_energy["gaussian_2_sigma_N_hits"]
    x_gauss_2_sigma_E = np.linspace(popt_gaussian_2_sigma_sum_E[2] - popt_gaussian_2_sigma_sum_E[1] * 2,
                          popt_gaussian_2_sigma_sum_E[2] + popt_gaussian_2_sigma_sum_E[1] * 2, 500)
    y_gauss_2_sigma_E = fit_gaussian(x_gauss_2_sigma_E, *popt_gaussian_2_sigma_sum_E)
    x_gauss_2_sigma_N = np.linspace(popt_gaussian_2_sigma_N_hits[2] - popt_gaussian_2_sigma_N_hits[1] * 2,
                          popt_gaussian_2_sigma_N_hits[2] + popt_gaussian_2_sigma_N_hits[1] * 2, 500)
    y_gauss_2_sigma_N = fit_gaussian(x_gauss_2_sigma_N, *popt_gaussian_2_sigma_N_hits)

    plt.figure(figsize=(15,4))
    plt.suptitle(f"Energy = {energy} GeV")
    
    plt.subplot(1,3,1)
    plt.hist(sum_E_arr, bins=50, color='blue', alpha=0.7)
    plt.plot(x_gauss_2_sigma_E, y_gauss_2_sigma_E, color='blue', label='Gaussian 2\u03C3 Fit')
    plt.plot(np.linspace(0, np.max(sum_E_arr), 100), fit_gamma(np.linspace(0, np.max(sum_E_arr), 100), *popt_gamma_sum_E), color='red', label='Gamma Fit')
    plt.plot(np.linspace(0, np.max(sum_E_arr), 100), fit_log_normal(np.linspace(0, np.max(sum_E_arr), 100), *popt_log_normal_sum_E), color='green', label='Log-normal Fit')
    plt.plot(np.linspace(0, np.max(sum_E_arr), 100), fit_gaussian(np.linspace(0, np.max(sum_E_arr), 100), *popt_gaussian_sum_E), color='orange', label='Gaussian Fit')
    plt.legend()
    plt.title("Sum Energy Distribution")
    plt.xlabel("Sum Energy")
    plt.ylabel("Number of Events")
    
    bins_for_low_E=np.arange(min(N_hits_arr), max(N_hits_arr) + 1.5) - 0.5
    plt.subplot(1,3,2)
    plt.hist(N_hits_arr, bins = bins_for_low_E if len(bins_for_low_E) < 50 else 50, color='green', alpha=0.7)
    plt.plot(x_gauss_2_sigma_N, y_gauss_2_sigma_N, color='blue', label='Gaussian 2\u03C3 Fit')
    plt.plot(np.linspace(0, np.max(N_hits_arr), 100), fit_gamma(np.linspace(0, np.max(N_hits_arr), 100), *popt_gamma_N_hits), color='red', label='Gamma Fit')
    plt.plot(np.linspace(0, np.max(N_hits_arr), 100), fit_log_normal(np.linspace(0, np.max(N_hits_arr), 100), *popt_log_normal_N_hits), color='green', label='Log-normal Fit')
    plt.plot(np.linspace(0, np.max(N_hits_arr), 100), fit_gaussian(np.linspace(0, np.max(N_hits_arr), 100), *popt_gaussian_N_hits), color='orange', label='Gaussian Fit')
    plt.legend()
    plt.title("Number of Hits Distribution")
    plt.ylabel("Number of Events")
    plt.xlabel("Number of Hits")
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])

   # plt.savefig(f"/home/llr/ilc/oprea/data/energy_plot_{energy}GeV.png")
    plt.show()

def find_peak(fit_results):
    peaks = {}
    for energy in fit_results.keys():
        params = fit_results[energy]
        peaks[energy] = {}

        # For sum_E fits
        x_min = 0
        x_max = max(params["gaussian_sum_E"][2], params["log_normal_sum_E"][2], params["gamma_sum_E"][2]) * 2
        x_sum_E = np.linspace(x_min, x_max, 1000)
        y_gamma = fit_gamma(x_sum_E, *params["gamma_sum_E"])
        y_logn = fit_log_normal(x_sum_E, *params["log_normal_sum_E"])
        y_gauss = fit_gaussian(x_sum_E, *params["gaussian_sum_E"])
        y_gauss2 = fit_gaussian(x_sum_E, *params["gaussian_2_sigma_sum_E"])

        peaks[energy]["gamma_sum_E_peak"] = x_sum_E[np.argmax(y_gamma)]
        peaks[energy]["log_normal_sum_E_peak"] = x_sum_E[np.argmax(y_logn)]
        peaks[energy]["gaussian_sum_E_peak"] = x_sum_E[np.argmax(y_gauss)]
        peaks[energy]["gaussian_2_sigma_sum_E_peak"] = x_sum_E[np.argmax(y_gauss2)]

        # For N_hits fits
        x_N_hits = np.linspace(0, 2 * params["gaussian_N_hits"][2], 1000)
        y_gamma_N = fit_gamma(x_N_hits, *params["gamma_N_hits"])
        y_logn_N = fit_log_normal(x_N_hits, *params["log_normal_N_hits"])
        y_gauss_N = fit_gaussian(x_N_hits, *params["gaussian_N_hits"])
        y_gauss2_N = fit_gaussian(x_N_hits, *params["gaussian_2_sigma_N_hits"])

        peaks[energy]["gamma_N_hits_peak"] = x_N_hits[np.argmax(y_gamma_N)]
        peaks[energy]["log_normal_N_hits_peak"] = x_N_hits[np.argmax(y_logn_N)]
        peaks[energy]["gaussian_N_hits_peak"] = x_N_hits[np.argmax(y_gauss_N)]
        peaks[energy]["gaussian_2_sigma_N_hits_peak"] = x_N_hits[np.argmax(y_gauss2_N)]
    
    return peaks


def add_parameters_to_txt(fit_results, raw_data, filename="fit_parameters.txt"):
    peaks = find_peak(fit_results)
    existing_energies = set()
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith("Energy:"):
                    try:
                        energy_val = float(line.split()[1])
                        existing_energies.add(energy_val)
                    except:
                        continue

    with open(filename, 'a') as f:
        for energy in fit_results.keys():
            if energy in existing_energies:
                print(f"Energy {energy} GeV already exists in file. Skipping.")
                continue
            print(f"Appending data for Energy {energy} GeV.")
            f.write(f"Energy: {energy} GeV\n")
            f.write("Sum_E Fit Parameters:\n")
            f.write(f"Gamma_E_{energy}: {', '.join(f'{x}' for x in fit_results[energy]['gamma_sum_E'])}\n")
            f.write(f"Peak of gamma_E dist at {energy} GeV: {peaks[energy]['gamma_sum_E_peak']}\n")
            f.write(f"Log-normal_E_{energy}: {', '.join(f'{x}' for x in fit_results[energy]['log_normal_sum_E'])}\n")
            f.write(f"Peak of ln_E dist at {energy} GeV: {peaks[energy]['log_normal_sum_E_peak']}\n")
            f.write(f"Gaussian_E_{energy}: {', '.join(f'{x}' for x in fit_results[energy]['gaussian_sum_E'])}\n")
            f.write(f"Peak of gaussian_E dist at {energy} GeV: {peaks[energy]['gaussian_sum_E_peak']}\n")
            f.write(f"Gaussian_2_sigma_E_{energy}: {', '.join(f'{x}' for x in fit_results[energy]['gaussian_2_sigma_sum_E'])}\n")
            f.write(f"Peak of gaussian 2 sigma_E dist at {energy} GeV: {peaks[energy]['gaussian_2_sigma_sum_E_peak']}\n")
            f.write(f"RMS90_sum_E_{energy}: {', '.join(f'{x}' for x in fit_results[energy]['rms90_sum_E'])}\n")
            f.write("N_hits Fit Parameters:\n")
            f.write(f"Gamma_N_{energy}: {', '.join(f'{x}' for x in fit_results[energy]['gamma_N_hits'])}\n")
            f.write(f"Peak of gamma_N dist at {energy} GeV: {peaks[energy]['gamma_N_hits_peak']}\n")
            f.write(f"Log-normal_N_{energy}: {', '.join(f'{x}' for x in fit_results[energy]['log_normal_N_hits'])}\n")
            f.write(f"Peak of ln_N dist at {energy} GeV: {peaks[energy]['log_normal_N_hits_peak']}\n")
            f.write(f"Gaussian_N_{energy}: {', '.join(f'{x}' for x in fit_results[energy]['gaussian_N_hits'])}\n")
            f.write(f"Peak of gaussian_N dist at {energy} GeV: {peaks[energy]['gaussian_N_hits_peak']}\n")
            f.write(f"Gaussian_2_sigma_N_{energy}: {', '.join(f'{x}' for x in fit_results[energy]['gaussian_2_sigma_N_hits'])}\n")
            f.write(f"Peak of gaussian 2 sigma_N dist at {energy} GeV: {peaks[energy]['gaussian_2_sigma_N_hits_peak']}\n")
            f.write(f"RMS90_N_hits_{energy}: {', '.join(f'{x}' for x in fit_results[energy]['rms90_N_hits'])}\n")
            f.write("\n")
            f.write(f"Raw Data for Energy {energy} GeV:\n")
            f.write(f"Sum_E_{energy} Array: {raw_data[energy]['sum_E_arr']}\n")
            f.write(f"N_hits_{energy} Array: {raw_data[energy]['N_hits_arr']}\n")
            f.write("\n")

    print("Done updating fit_parameters.txt.")

    #add the parameters into a txt, then read from the txt into the resolution_and_linearity.py

energies = [0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0]

# Analyze all distributions at given energies
fit_results, raw_data = analyse_distributions(energies)

add_parameters_to_txt(
    fit_results,
    raw_data,
    filename="fit_parameters.txt"
)

# Save all fit parameters and relevant data to a text file
print("Saving to:", os.path.abspath("fit_parameters.txt")) #Seeing where the file is saved

# Plot fits for a specific energy (e.g., 0.05 GeV)
energy_to_plot = 0.05  # Change this to the desired energy for plotting
plot_fits(
    raw_data[energy_to_plot]["sum_E_arr"],
    raw_data[energy_to_plot]["N_hits_arr"],
    fit_results,
    energy_to_plot
)




# This code will analyze the energy distributions, fit them, plot the results, and save the parameters to a text file.
# You can adjust the energies list to analyze different sets of energies as needed.
# The plotting function will visualize the fits for a specific energy, and the parameters will be saved to a text file for later use.
# Make sure to have the necessary libraries installed (ROOT, numpy, matplotlib, scipy) and the ROOT file structure as expected.
# The code assumes that the ROOT files are structured correctly and contain the expected branches.
# If you encounter any issues, check the ROOT file structure and the branch names.
# The code also includes error handling for the fitting process,
# so if a fit fails, it will print an error message but continue processing other fits.




#cOMBINE THE HITS (FOR SUM E AND NUMEBER OF HITS)
#save the comine hits method to a root file
#threshold