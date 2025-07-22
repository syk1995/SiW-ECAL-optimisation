import numpy as np
#import ROOT
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
from scipy.optimize import curve_fit
from scipy.stats import landau as landau_dist

def landau_fit(x, norm, mu, sigma):
    return norm * landau_dist.pdf(x, loc=mu, scale=sigma)

def read_file(input_energy, particle_type):
    path = "/mnt/d/data_roots/version3_combined_cells_sum_and_N_mu-_200.0GeV_gamma_" 
    filename = f"{path}{float(input_energy)}GeV.root"
    with uproot.open(filename) as file:
        # Navigate to the appropriate group (mu- or gamma)
        group = file[particle_type]

        # Read sum_E and N_hits
        sum_E = group["sum_E"].arrays(library="np")["sum_E"]
        N_hits = group["N_hits"].arrays(library="np")["N_hits"]
        #print(f"sum_E array: {sum_E}")
        #print(f"N_hits array: {N_hits}")
    return sum_E, N_hits

def get_hit_energy(input_energy, particle_type):
    path = "/mnt/d/data_roots/version3_combined_cells_sum_and_N_mu-_200.0GeV_gamma_" 
    filename = f"{path}{float(input_energy)}GeV.root"
    with uproot.open(filename) as file:
        # Navigate to the appropriate group (mu- or gamma)
        group = file[particle_type]

        # Read hit energy
        hit_energies_for_one_event = group["hit_energies_for_one_event"].arrays(library="np")["hit_energies_for_one_event"]
        print(f"hit_energy array: {hit_energies_for_one_event}")

    return hit_energies_for_one_event * 1000 # Convert to MeV

def get_hit_energy_for_all_events(input_energy, particle_type):
    path = "/mnt/d/data_roots/version3_combined_cells_sum_and_N_mu-_200.0GeV_gamma_" 
    filename = f"{path}{float(input_energy)}GeV.root"
    with uproot.open(filename) as file:
        # Navigate to the appropriate group (mu- or gamma)
        group = file[particle_type]

        # Read hit energy for all events
        hit_energies_for_all_events = group["hit_energy"].arrays(library="np")["hit_energy"]
        print(f"hit_energy array: {hit_energies_for_all_events}")

    return hit_energies_for_all_events * 1000 # Convert to MeV

def get_mip(input_energy):
    hit_energies_for_one_event = get_hit_energy(input_energy, "mu-")
    number_of_events, bin_edges= np.histogram(np.array(hit_energies_for_one_event), bins=100, density=False)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    hist_peak_E = bin_centers[np.argmax(number_of_events)]
    norm_guess = np.max(number_of_events)
    mu_guess = hist_peak_E if hist_peak_E > 0.1 else 0.2
    sigma_guess = (np.max(hit_energies_for_one_event) - np.min(hit_energies_for_one_event)) / 5

    # Apply mask to ignore 0 bins if needed
    mask = bin_centers > 0.05
    bin_centers = bin_centers[mask]
    number_of_events = number_of_events[mask]

    # Bounds: force mu > 0.05
    bounds = ([0, 0.05, 0], [np.inf, np.max(hit_energies_for_one_event), np.inf])

    fit_params_mu, _ = curve_fit(
        landau_fit,
        bin_centers,
        number_of_events,
        p0=[norm_guess, mu_guess, sigma_guess],
        bounds=bounds,
        maxfev=100000
    )
    print(f"Fit parameters for mu- at {input_energy} GeV: norm={fit_params_mu[0]:.2f}, mu={fit_params_mu[1]:.2f}, sigma={fit_params_mu[2]:.2f}")


    return fit_params_mu

def threshold(input_energy):
    mip = get_mip(input_energy)[1]
    threshold_value = 0.5 * mip
    print(f"Threshold for {input_energy} GeV: {threshold_value:.2f} GeV")
    return threshold_value

def analyse_distribution(input_energy_mu, input_energy_gamma):
    for energy_mu in input_energy_mu:
        print(f"Processing mu- energy: {energy_mu} GeV")
        for energy_gamma in input_energy_gamma:
            print(f"Processing gamma energy: {energy_gamma} GeV")
            sum_E_gamma, N_hits_gamma = read_file(energy_gamma, "gamma")
            hit_energy_mu = get_hit_energy(energy_gamma, "mu-")
            hit_energy_gamma_for_all_events = get_hit_energy_for_all_events(energy_gamma, "gamma")
            threshold_value = threshold(energy_gamma)
            landau_params_mu = get_mip(energy_gamma)
            hit_energy_gamma_thresholded = []
            sum_E_gamma_thresholded = []
            N_hits_gamma_thresholded = []

            for event_hits in hit_energy_gamma_for_all_events:
                event_hits = np.array(event_hits)
                hits_above_thresh = event_hits[event_hits >= threshold_value]

                # Collect flattened hit energies for plotting the energy distribution
                hit_energy_gamma_thresholded.extend(hits_above_thresh)

                # Collect sum E and N hits per event
                sum_E_gamma_thresholded.append(np.sum(hits_above_thresh))
                N_hits_gamma_thresholded.append(len(hits_above_thresh))

            # Convert to numpy arrays for consistency
            hit_energy_gamma_thresholded = np.array(hit_energy_gamma_thresholded)
            sum_E_gamma_thresholded = np.array(sum_E_gamma_thresholded)
            N_hits_gamma_thresholded = np.array(N_hits_gamma_thresholded)

            number_of_events, bin_edges= np.histogram(np.array(hit_energy_gamma_thresholded), bins="auto", density=False)
            bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

            plt.figure(figsize=(10, 6))
            plt.subplot(2,2,3)
            plt.hist(hit_energy_gamma_thresholded, bins='auto', density=False, alpha=0.5, label='Thresholded gamma energy deposition')
            plt.title(f"Energy Deposition Distribution for {energy_gamma} GeV gamma")
            plt.xlabel("Energy Deposition (GeV)")
            plt.ylabel("Number of Events")
            plt.axvline(threshold_value, color='red', linestyle='dashed', linewidth=1, label='Threshold')
            plt.legend()

            plt.subplot(2,2,4)
            plt.hist(N_hits_gamma_thresholded, bins='auto', density=False, alpha=0.5, label='Thresholded gamma hits')
            plt.title(f"Number of Hits Distribution for {energy_gamma} GeV gamma")
            plt.xlabel("Number of Hits")
            plt.ylabel("Number of Events")
            plt.legend()


            plt.subplot(2,2,1)
            x = np.linspace(0, np.max(hit_energy_mu) * 1.1, 500)
            plt.hist(hit_energy_mu, bins=100, density=False, alpha=0.6, label="Hit Energy Distribution")
            plt.plot(x, landau_fit(x, *landau_params_mu), 'r-', label="Landau Fit")
            plt.xlabel("Energy [GeV]")
            plt.ylabel("Counts")
            plt.legend()

            plt.subplot(2, 2, 2)
            plt.hist(sum_E_gamma_thresholded, bins='auto', density=False, alpha=0.5, label='Thresholded gamma sum E')
            plt.title(f"Sum E Distribution for {energy_gamma} GeV gamma")
            plt.xlabel("Sum Energy per Event (GeV)")
            plt.ylabel("Number of Events")
            plt.legend()
            plt.show()



energies_mu = [200.0]
energies_gamma = [0.15, 0.2, 0.25, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
analyse_distribution(energies_mu, energies_gamma)
# Notes for myself to clarify the MIP concept:
# The MIP is used to define the threshold for energy deposition in the calorimeter because it 
# represents the typical energy deposited by a high-energy muon (mu-), which behaves as a minimum ionizing particle.
# Such particles are highly penetrating and deposit a stable, predictable amount of energy per unit length.
# 
# The threshold is used to separate real signals from noise in the calorimeter.
# The MIP is a universal reference for any calorimeter because it depends on the detector's material, thickness,
# and geometry, and thus can also serve as a reference for thresholds, which are similarly detector-dependent.
# 
# Typically, energy depositions below ~0.5 MIP are considered noise, while depositions above this value are 
# considered likely to originate from real charged particles (signal).
#
# The Landau distribution describes the energy deposition fluctuations for MIPs.
# The Most Probable Value (MPV) of this distribution is typically used to define 1 MIP units.


#fie gasesti un mod in care vezi hit_energies_array cum trebe aici, fie il scrii direct in combined hits si-ti
#bagi picioarele in root-ul ma-sii