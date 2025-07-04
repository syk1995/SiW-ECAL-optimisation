import numpy as np
import matplotlib.pyplot as plt
import ast
import math


energies = [0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 60.0]
# Implementing a function to retrieve the params of the fits, for the resolution and linearity plots

def fixing_commas(line):
    parts = line.strip().split(':', 1)
    if len(parts) < 2:
        return None
    values = parts[1].strip()
    if values.startswith('[') and values.endswith(']'):
        inner = values[1:-1].strip()
        fixed = ', '.join(inner.split())
#        print(f"Fixed commas: {fixed}")
        return f"{parts[0]}: [{fixed}]"
    else:
#        print(f"Line without brackets: {line}")
        return line  

def get_fit_params_from_file(filename='/grid_mnt/vol_home/llr/ilc/oprea/code/Energy-Reco/Ana/fit_parameters.txt'):
    gamma_params_sum_E = {}
    log_normal_params_sum_E = {}
    gaussian_params_sum_E = {}
    gaussian_2_sigma_params_sum_E = {}
    rms90_params_sum_E = {}
    gamma_N_hits_params = {}
    log_normal_N_hits_params = {}
    gaussian_N_hits_params = {}
    gaussian_2_sigma_N_hits_params = {}
    rms90_N_hits_params = {}
    # Reading the file and extracting the parameters
    if not filename.endswith('.txt'):
        raise ValueError("Filename must end with '.txt'")
    current_energy = None
    with open(filename, 'r') as f:
        for line in f:
            line = fixing_commas(line.strip())
            if line is None or not line:
                continue
            if line.startswith("Energy"):
                try:
                    current_energy = float(line.split()[1])
                    if current_energy < 1.0:
                        current_energy = current_energy
                    elif current_energy >= 100.0:
                        current_energy = int(current_energy)
                except:
                    continue
            elif current_energy is not None:
                if line.startswith(f"Gamma_E_{current_energy}:"):
                    gamma_params_sum_E[current_energy] = ast.literal_eval(line.split(":", 1)[1].strip())
                elif line.startswith(f"Log-normal_E_{current_energy}:"):
                    log_normal_params_sum_E[current_energy] = ast.literal_eval(line.split(":", 1)[1].strip())
                elif line.startswith(f"Gaussian_E_{current_energy}:"):
                    gaussian_params_sum_E[current_energy] = ast.literal_eval(line.split(":", 1)[1].strip())
                elif line.startswith(f"Gaussian_2_sigma_E_{current_energy}:"):
                    gaussian_2_sigma_params_sum_E[current_energy] = ast.literal_eval(line.split(":", 1)[1].strip())
                elif line.startswith(f"RMS90_sum_E_{current_energy}:"):
                    rms90_params_sum_E[current_energy] = ast.literal_eval(line.split(":", 1)[1].strip())
                elif line.startswith(f"Gamma_N_{current_energy}:"):
                    gamma_N_hits_params[current_energy] = ast.literal_eval(line.split(":", 1)[1].strip())
                elif line.startswith(f"Log-normal_N_{current_energy}:"):
                    log_normal_N_hits_params[current_energy] = ast.literal_eval(line.split(":", 1)[1].strip())
                elif line.startswith(f"Gaussian_N_{current_energy}:"):
                    gaussian_N_hits_params[current_energy] = ast.literal_eval(line.split(":", 1)[1].strip())
                elif line.startswith(f"Gaussian_2_sigma_N_{current_energy}:"):
                    gaussian_2_sigma_N_hits_params[current_energy] = ast.literal_eval(line.split(":", 1)[1].strip())
                elif line.startswith(f"RMS90_N_hits_{current_energy}:"):
                    rms90_N_hits_params[current_energy] = ast.literal_eval(line.split(":", 1)[1].strip())
        print("Extracted parameters:")
        print("Gamma (Sum E):", gamma_params_sum_E)
#        print("Log-normal (Sum E):", log_normal_params_sum_E)
#        print("Gaussian (Sum E):", gaussian_params_sum_E)
#        print("Gaussian 2 sigma (Sum E):", gaussian_params_sum_E)
        print("Gamma (N Hits):", gamma_N_hits_params)
#        print("Log-normal (N Hits):", log_normal_N_hits_params)
#        print("Gaussian (N Hits):", gaussian_N_hits_params)
#        print("Gaussian 2 sigma (N Hits):", gaussian_N_hits_params)
        return {
        'gamma_sum_E': gamma_params_sum_E,
        'log_normal_sum_E': log_normal_params_sum_E,
        'gaussian_sum_E': gaussian_params_sum_E,
        'gaussian_2_sigma_sum_E': gaussian_2_sigma_params_sum_E,
        'rms90_sum_E': rms90_params_sum_E,
        'gamma_N_hits': gamma_N_hits_params,
        'log_normal_N_hits': log_normal_N_hits_params,
        'gaussian_N_hits': gaussian_N_hits_params,
        'gaussian_2_sigma_N_hits': gaussian_2_sigma_N_hits_params,
        'rms90_N_hits': rms90_N_hits_params
    }


def get_peaks_from_file(filename='/grid_mnt/vol_home/llr/ilc/oprea/code/Energy-Reco/Ana/fit_parameters.txt'):
    peaks = {
        'sum_E': {
            'gamma': {},
            'log_normal': {},
            'gaussian': {},
            'gaussian_2_sigma': {}
        },
        'N_hits': {
            'gamma': {},
            'log_normal': {},
            'gaussian': {},
            'gaussian_2_sigma': {}
        }
    }

    if not filename.endswith('.txt'):
        raise ValueError("Filename must end with '.txt'")

    current_energy = None
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("Energy"):
                try:
                    current_energy = float(line.split()[1])
                    if current_energy >= 1.0:
                        current_energy = int(current_energy)
                except Exception:
                    continue

            elif current_energy is not None and line.startswith("Peak of"):
                try:
                    parts = line.split(":")
                    if len(parts) != 2:
                        continue

                    peak_type = parts[0].lower()  # "peak of gamma_E dist at 0.1 gev"
                    value = float(parts[1].strip())

                    # Determine which peak this is
                    if "gamma_e" in peak_type:
                        peaks['sum_E']['gamma'][current_energy] = value
                    elif "log_normal_e" in peak_type or "ln_e" in peak_type:
                        peaks['sum_E']['log_normal'][current_energy] = value
                    elif "gaussian 2 sigma_e" in peak_type:
                        peaks['sum_E']['gaussian_2_sigma'][current_energy] = value
                    elif "gaussian_e" in peak_type:
                        peaks['sum_E']['gaussian'][current_energy] = value

                    elif "gamma_n" in peak_type:
                        peaks['N_hits']['gamma'][current_energy] = value
                    elif "log_normal_n" in peak_type or "ln_n" in peak_type:
                        peaks['N_hits']['log_normal'][current_energy] = value
                    elif "gaussian 2 sigma_n" in peak_type:
                        peaks['N_hits']['gaussian_2_sigma'][current_energy] = value
                    elif "gaussian_n" in peak_type:
                        peaks['N_hits']['gaussian'][current_energy] = value

                except Exception:
                    continue

    return peaks
            
def get_fit_params_for_sum_E(fit_type='gamma'):
    fit_params = get_fit_params_from_file()
    key_map = {
        'gamma': 'gamma_sum_E',
        'log_normal': 'log_normal_sum_E',
        'gaussian': 'gaussian_sum_E',
        'gaussian_2_sigma': 'gaussian_2_sigma_sum_E',
        'rms90': 'rms90_sum_E'
    }
    return fit_params[key_map[fit_type]]

def get_fit_params_for_N_hits(fit_type='gamma'):
    fit_params = get_fit_params_from_file()
    key_map = {
        'gamma': 'gamma_N_hits',
        'log_normal': 'log_normal_N_hits',
        'gaussian': 'gaussian_N_hits',
        'gaussian_2_sigma': 'gaussian_2_sigma_N_hits',
        'rms90': 'rms90_N_hits'
    }
    return fit_params[key_map[fit_type]]


def calculate_resolution_and_linearity_for_sum_E(sum_E_params, peaks_dict, fit_type='gamma'):
    """
    Calculates means, stds, resolutions (sigma/mu), and resolutions_peak (peak/mu) for sum_E.
    peaks_dict should be the output of get_peaks_from_file()['sum_E'].
    """
    means = {}
    stds = {}
    resolutions = {}
    resolutions_peak = {}
    if fit_type not in ['gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', 'rms90']:
        raise ValueError("fit_type must be one of 'gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma, or 'rms90'")
    means[fit_type] = {}
    stds[fit_type] = {}
    resolutions[fit_type] = {}
    resolutions_peak[fit_type] = {}

    # peaks_dict is expected to be e.g. get_peaks_from_file()['sum_E']['gamma'] for fit_type='gamma'
    peaks_for_type = peaks_dict.get(fit_type, {})

    for energy in sorted(sum_E_params.keys()):
        p = sum_E_params[energy]
        peak = peaks_for_type.get(energy, None)
        if fit_type == 'gamma':
            mean = p[1] * p[3]
            std = np.sqrt(p[1]) * p[3]
        elif fit_type == 'log_normal':
            mean = p[2] * np.exp(0.5 * p[1]**2)
            std = p[2] * np.sqrt(np.exp(p[1]**2) - 1) * np.exp(0.5 * p[1]**2)
        elif fit_type == 'gaussian':
            mean = p[2]
            std = p[1]
        elif fit_type == 'gaussian_2_sigma':
            mean = p[2]
            std = p[1]
        elif fit_type == 'rms90':
            mean = p[0]
            std = p[1]
        resolution = std / mean
        # Only compute resolution_peak if peak is available
        resolution_peak = peak / mean if peak is not None else None
        means[fit_type][energy] = mean
        stds[fit_type][energy] = std
        resolutions[fit_type][energy] = resolution
        resolutions_peak[fit_type][energy] = resolution_peak
    return means, stds, resolutions, resolutions_peak

# dictionaries will be something like:
# means = {
#       'gamma' = {
#           energy1 : val1
#           energy2 : val2
#           ...
#           }
#       'log_normal' = {
#           energy1 : val1
#           energy2 : val2
#           ...
#           }
#       'gaussian' = {
#           energy1 : val1
#           energy2 : val2
#           ...
#           }
#   }


def calculate_resolution_and_linearity_for_N_hits(N_hits_params, fit_type='gamma'):
    """
    Calculates means, stds, and resolutions (sigma/mu) for N_hits.
    """
    means = {}
    stds = {}
    resolutions = {}
    resolutions_peak = {}
    if fit_type not in ['gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', 'rms90']:
        raise ValueError("fit_type must be one of 'gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma, or 'rms90'")
    means[fit_type] = {}
    stds[fit_type] = {}
    resolutions[fit_type] = {}
    resolutions_peak[fit_type] = {}

    # peaks_dict is expected to be e.g. get_peaks_from_file()['N_hits']['gamma'] for fit_type='gamma'
    peaks_dict = get_peaks_from_file()['N_hits']
    peaks_for_type = peaks_dict.get(fit_type, {})

    for energy in sorted(N_hits_params.keys()):
        p = N_hits_params[energy]
        peak = peaks_for_type.get(energy, None)
        if fit_type == 'gamma':
            mean = p[1] * p[2]
            std = np.sqrt(p[1]) * p[2]
        elif fit_type == 'log_normal':
            mean = p[2] * np.exp(0.5 * p[1]**2)
            std = p[2] * np.sqrt(np.exp(p[1]**2) - 1) * np.exp(0.5 * p[1]**2)
        elif fit_type == 'gaussian':
            mean = p[2]
            std = p[1]
        elif fit_type == 'gaussian_2_sigma':
            mean = p[2]
            std = p[1]
        elif fit_type == 'rms90':
            mean = p[0]
            std = p[1]
        resolution = std / mean
        resolution_peak = peak / mean if peak is not None else None
        means[fit_type][energy] = mean
        stds[fit_type][energy] = std
        resolutions[fit_type][energy] = resolution
        resolutions_peak[fit_type][energy] = resolution_peak
    return means, stds, resolutions, resolutions_peak
#The dictionary will be the same as for the sum_E function


def plot_resolution_and_linearity_vs_E_for_one_fit_type(energies, means, stds, resolutions, fit_type="gamma"):
    if fit_type not in ['gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', 'rms90']:
        raise ValueError("fit_type must be one of 'gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma' or 'rms90")
    plt.figure(figsize = (15, 4))
    plt.suptitle(f'Resolution and Linearity for {fit_type} fit type (sum E)', fontsize=16)

    plt.subplot(1, 3, 1)
    plt.plot(energies, resolutions, 'o', ms = 3, label=f'Resolution ({fit_type})')
    plt.xlabel('Energy (GeV)')
    plt.ylabel('Energy Resolution (\u03C3 / \u03BC)')
    plt.title(f'Resolution vs Energy ({fit_type})')
    plt.grid()
    plt.legend()

    plt.subplot(1, 3, 2)
    plt.plot(energies, means, 'o', ms = 3, label=f'Mean Energy ({fit_type})')
    plt.xlabel('Primary Energy (GeV)')
    plt.ylabel('Mean Fitted Energy (GeV)')
    plt.title(f'Linearity: Mean vs Beam Energy ({fit_type})')
    plt.grid()
    plt.legend()

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

def plot_resolution_and_linearity_vs_N_for_one_fit_type(energies, means, stds, resolutions, fit_type='gamma'):
    if fit_type not in ['gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', 'rms90']:
        raise ValueError("fit_type must be one of 'gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', or 'rms90")
    plt.figure(figsize=(15, 4))
    plt.suptitle(f'Resolution and Linearity for {fit_type} fit type (N Hits)', fontsize=16)

    plt.subplot(1, 3, 1)
    plt.plot(energies, resolutions, 'o', label=f'Resolution ({fit_type})')
    plt.xlabel('Primary Energy (GeV)')
    plt.ylabel('Energy Resolution (\u03C3 / \u03BC)')
    plt.title(f'Resolution vs Energy ({fit_type})')
    plt.grid()
    plt.legend()

    plt.subplot(1, 3, 2)
    plt.plot(energies, means, 'o', label=f'Mean Energy ({fit_type})', capsize=5)
    plt.xlabel('Primary Energy (GeV)')
    plt.ylabel('Mean Fitted Energy (GeV)')
    plt.title(f'Linearity: Mean vs Beam Energy ({fit_type})')
    plt.grid()
    plt.legend()

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()



def plot_resolution_and_linearity_vs_E_for_all_fit_types(energies, means, stds, resolutions, resolutions_peak):
    peaks = get_peaks_from_file()
    plt.figure(figsize=(15, 9))
    plt.suptitle('Resolution and Linearity for All Fit Types', fontsize=16)

    plt.subplot(2, 3, 4)
    plt.plot(energies, resolutions['gamma'], 'o', ms=3, label='Gamma Resolution', color='blue')
    plt.plot(energies, resolutions['log_normal'], 'o', ms=3, label='Log-normal Resolution', color='yellow')
   # plt.plot(energies, resolutions['gaussian'], 'o', ms=3, label='Gaussian Resolution', color='red')
    plt.plot(energies, resolutions['gaussian_2_sigma'], 'o', ms=3, label='Gaussian 2 Sigma Resolution', color='green')
    plt.plot(energies, resolutions['rms90'], 'o', ms=3, label='RMS90 Resolution', color='magenta')
    plt.xlabel('Primary Energy (GeV)')
    plt.ylabel('Energy Resolution (\u03C3 / \u03BC)')
    plt.title('Resolution vs Energy')
    plt.grid()
    plt.legend()

    plt.subplot(2, 3, 2)
    plt.plot(energies, means['gamma'], 'o', ms=3, label='Gamma Mean Energy', color='blue')
    plt.plot(energies, means['log_normal'], 'o', ms=3, label='Log-normal Mean Energy', color='yellow')
    #plt.plot(energies, means['gaussian'], 'o', ms=3, label='Gaussian Mean Energy', color='red')
    plt.plot(energies, means['gaussian_2_sigma'], 'o', ms=3, label='Gaussian 2 Sigma Mean Energy', color='green')
    plt.plot(energies, means['rms90'], 'o', ms=3, label='RMS90 Mean Energy', color='magenta')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Primary Energy (GeV)')
    plt.ylabel('Mean Fitted Energy (MeV)')
    plt.title('Linearity of Means: Mean vs Beam Energy')
   # plt.grid()

    # Linear-Fits for all distributions
    #plt.subplot(2, 3, 4)
    energies_np = np.array(energies)
    #plt.xscale('log')
    #plt.yscale('log')

    for key, color in zip(['gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', 'rms90'],
                        ['blue', 'yellow', 'red', 'green', 'magenta']):
        y = np.array([means[key][i] for i in range(len(energies))])
        mask = (energies_np > 0) & (y > 0)

        coeffs_lin = np.polyfit(energies_np[mask], y[mask], 1)
        label = f"{key} linear fit: y = {coeffs_lin[0]:.2f}x + {coeffs_lin[1]:.2f}"

        if plt.gca().get_xscale() == 'log' and plt.gca().get_yscale() == 'log':
            coeffs_log = np.polyfit(np.log10(energies_np[mask]), np.log10(y[mask]), 1)
            fit_y = 10**(coeffs_log[1]) * energies_np**coeffs_log[0]
        else:
            fit_y = coeffs_lin[0] * energies_np + coeffs_lin[1]

        plt.plot(energies_np, fit_y, '--', color=color, alpha=0.7)
        
       # plt.xlabel('Primary Energy (GeV)')
       # plt.ylabel('Mean Fitted Energy (MeV)')
       # plt.title('Linear-Fits of all distributions')
        plt.grid()
        plt.legend(loc = 'lower right', fontsize = 'x-small')

    # (Mean Energy - Fitted Energy)/Fitted Energy for all distributions
    plt.subplot(2, 3, 6)
    for key, color in zip(['gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', 'rms90'],
                          ['blue', 'yellow', 'red', 'green', 'magenta']):
        # (E - E_fitted) / E_fitted for each energy
        rel_diff = [np.abs((energy - means[key][energies.index(energy)])) / means[key][energies.index(energy)] for energy in energies]
        plt.plot(energies, rel_diff, 'o-', color=color, label=f'{key}')
    plt.xlabel('Primary Energy (GeV)')
    plt.ylabel('(E - E_fitted)/E_fitted (%)')
    plt.title('Non-linearity Ratio')
    plt.grid()
    plt.ylim(0.975, 0.995)
    plt.legend(fontsize='x-small')

    gamma_ratio = [means['gamma'][i] / peaks['sum_E']['gamma'][energy] for i, energy in enumerate(energies)]
    log_normal_ratio = [means['log_normal'][i] / peaks['sum_E']['log_normal'][energy] for i, energy in enumerate(energies)]
    gaussian_ratio = [means['gaussian'][i] / peaks['sum_E']['gaussian'][energy] for i, energy in enumerate(energies)]
    gaussian_2_sigma_ratio = [means['gaussian_2_sigma'][i] / peaks['sum_E']['gaussian_2_sigma'][energy] for i, energy in enumerate(energies)]

    plt.subplot(2, 3, 3)
    plt.plot(energies, gamma_ratio, 'o', ms=3, label='Gamma Ratio Mean/Peak', color='blue')
    plt.plot(energies, log_normal_ratio, 'o', ms=3, label='Log-normal Ratio Mean/Peak', color='yellow')
   # plt.plot(energies, gaussian_ratio, 'o', ms=3, label='Gaussian Ratio Mean/Peak', color='red')
    plt.plot(energies, gaussian_2_sigma_ratio, 'o', ms=3, label='Gaussian 2 Sigma Ratio Mean/Peak', color='green')
    plt.xlabel('Primary Energy (GeV)')
    plt.ylabel('Mean/peak Ratio')
    plt.title('Mean/Peak vs Beam Energy')
    plt.grid()
    plt.legend()

    peaks_list_gamma = [peaks['sum_E']['gamma'][energy] for energy in energies]
    peaks_list_log_normal = [peaks['sum_E']['log_normal'][energy] for energy in energies]
    peaks_list_gaussian = [peaks['sum_E']['gaussian'][energy] for energy in energies]
    peaks_list_gaussian_2_sigma = [peaks['sum_E']['gaussian_2_sigma'][energy] for energy in energies]

    plt.subplot(2, 3, 5)
    plt.plot(energies, peaks_list_gamma, 'o', ms=3, label='Gamma Linearity', color='blue')
    plt.plot(energies, peaks_list_log_normal, 'o', ms=3, label='Log-Normal Linearity', color='yellow')
   # plt.plot(energies, peaks_list_gaussian, 'o', ms=3, label='Gaussian Linearity', color='red')
    plt.plot(energies, peaks_list_gaussian_2_sigma, 'o', ms=3, label='Gaussian 2 Sigma Linearity', color='green')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Primary Energy (GeV)')
    plt.ylabel('Peaks')
    plt.title('Linearity of Peaks vs Beam Energy')
    plt.grid()
    plt.legend()

    plt.subplot(2,3,1)
    plt.plot(energies, resolutions_peak['gamma'], 'o', ms=3, label='Gamma Resolution Peak', color='blue')
    plt.plot(energies, resolutions_peak['log_normal'], 'o', ms=3, label='Log-normal Resolution Peak', color='yellow')
    #plt.plot(energies, resolutions_peak['gaussian'], 'o', ms=3, label='Gaussian Resolution Peak', color='red')
    plt.plot(energies, resolutions_peak['gaussian_2_sigma'], 'o', ms=3, label='Gaussian 2 Sigma Resolution Peak', color='green')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig("/home/llr/ilc/oprea/data/resolution_and_linearity_vs_E.png")
    plt.show()

def plot_resolution_and_linearity_vs_N_for_all_fit_types(energies, means, stds, resolutions):
    peaks = get_peaks_from_file()
    plt.figure(figsize=(15, 9))
    plt.suptitle('Resolution and Linearity for All Fit Types (N Hits)', fontsize=16)

    plt.subplot(2, 3, 1)
    plt.plot(energies, resolutions['gamma'], 'o', ms=3, label='Gamma Resolution', color='blue')
    plt.plot(energies, resolutions['log_normal'], 'o', ms=3, label='Log-normal Resolution', color='yellow')
    #plt.plot(energies, resolutions['gaussian'], 'o', ms=3, label='Gaussian Resolution', color='red')
    plt.plot(energies, resolutions['gaussian_2_sigma'], 'o', ms=3, label='Gaussian 2 Sigma Resolution', color='green')
    plt.plot(energies, resolutions['rms90'], 'o', ms=3, label='RMS90 Resolution', color='magenta')
    plt.xlabel('Primary Energy (GeV)')
    plt.ylabel('Number of Hits Resolution (\u03C3 / \u03BC)')
    plt.title('Resolution vs Energy')
    plt.grid()
    plt.legend()

    plt.subplot(2, 3, 2)
    plt.plot(energies, means['gamma'], 'o', ms=3, label='Gamma Mean N Hits', color='blue')
    plt.plot(energies, means['log_normal'], 'o', ms=3, label='Log-normal Mean N Hits', color='yellow')
    #plt.plot(energies, means['gaussian'], 'o', ms=3, label='Gaussian Mean N Hits', color='red')
    plt.plot(energies, means['gaussian_2_sigma'], 'o', ms=3, label='Gaussian 2 Sigma Mean N Hits', color='green')
    plt.plot(energies, means['rms90'], 'o', ms=3, label='RMS90 Mean N Hits', color='magenta')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Primary Energy (GeV)')
    plt.ylabel('Mean Fitted N Hits')
    plt.title('Linearity of Means: Mean vs Beam Energy')
    plt.grid()
    plt.legend()

    # Linear-Fits for all distributions
    energies_np = np.array(energies)
    for key, color in zip(['gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', 'rms90'],
                          ['blue', 'yellow', 'red', 'green', 'magenta']):
        y = np.array([means[key][i] for i in range(len(energies))])
        mask = (energies_np > 0) & (y > 0)
        coeffs_lin = np.polyfit(energies_np[mask], y[mask], 1)
        label = f"{key} linear fit: y = {coeffs_lin[0]:.2f}x + {coeffs_lin[1]:.2f}"
        if plt.gca().get_xscale() == 'log' and plt.gca().get_yscale() == 'log':
            coeffs_log = np.polyfit(np.log10(energies_np[mask]), np.log10(y[mask]), 1)
            fit_y = 10**(coeffs_log[1]) * energies_np**coeffs_log[0]
        else:
            fit_y = coeffs_lin[0] * energies_np + coeffs_lin[1]
        plt.plot(energies_np, fit_y, '--', color=color, alpha=0.7)
        plt.grid()
        plt.legend(loc='lower right', fontsize='x-small')

    # (Mean N_hits - Fitted N_hits)/Fitted N_hits for all distributions
    plt.subplot(2, 3, 6)
    for key, color in zip(['gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', 'rms90'],
                          ['blue', 'yellow', 'red', 'green', 'magenta']):
        rel_diff = [np.abs((energy - means[key][energies.index(energy)])) / means[key][energies.index(energy)] for energy in energies]
        plt.plot(energies, rel_diff, 'o-', color=color, label=f'{key}')
    plt.xlabel('Primary Energy (GeV)')
    plt.ylabel('(N - N_fitted)/N_fitted (%)')
    plt.title('Non-linearity Ratio')
    plt.grid()

    plt.legend(fontsize='x-small')

    gamma_ratio = [means['gamma'][i] / peaks['N_hits']['gamma'][energy] for i, energy in enumerate(energies)]
    log_normal_ratio = [means['log_normal'][i] / peaks['N_hits']['log_normal'][energy] for i, energy in enumerate(energies)]
    gaussian_ratio = [means['gaussian'][i] / peaks['N_hits']['gaussian'][energy] for i, energy in enumerate(energies)]
    gaussian_2_sigma_ratio = [means['gaussian_2_sigma'][i] / peaks['N_hits']['gaussian_2_sigma'][energy] for i, energy in enumerate(energies)]

    plt.subplot(2, 3, 3)
    plt.plot(energies, gamma_ratio, 'o', ms=3, label='Gamma Ratio Mean/Peak', color='blue')
    plt.plot(energies, log_normal_ratio, 'o', ms=3, label='Log-normal Ratio Mean/Peak', color='yellow')
    #plt.plot(energies, gaussian_ratio, 'o', ms=3, label='Gaussian Ratio Mean/Peak', color='red')
    plt.plot(energies, gaussian_2_sigma_ratio, 'o', ms=3, label='Gaussian 2 Sigma Ratio Mean/Peak', color='green')
    plt.xlabel('Primary Energy (GeV)')
    plt.ylabel('Mean/Peak Ratio')
    plt.title('Mean/Peak vs Beam Energy')
    plt.grid()
    plt.legend()

    peaks_list_gamma = [peaks['N_hits']['gamma'][energy] for energy in energies]
    peaks_list_log_normal = [peaks['N_hits']['log_normal'][energy] for energy in energies]
    #peaks_list_gaussian = [peaks['N_hits']['gaussian'][energy] for energy in energies]
    peaks_list_gaussian_2_sigma = [peaks['N_hits']['gaussian_2_sigma'][energy] for energy in energies]

    plt.subplot(2, 3, 5)
    plt.plot(energies, peaks_list_gamma, 'o', ms=3, label='Gamma Linearity', color='blue')
    #plt.plot(energies, peaks_list_log_normal, 'o', ms=3, label='Log-Normal Linearity', color='yellow')
    #plt.plot(energies, peaks_list_gaussian, 'o', ms=3, label='Gaussian Linearity', color='red')
    plt.plot(energies, peaks_list_gaussian_2_sigma, 'o', ms=3, label='Gaussian 2 Sigma Linearity', color='green')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Primary Energy (GeV)')
    plt.ylabel('Peaks')
    plt.title('Linearity of Peaks vs Beam Energy')
    plt.grid()
    plt.legend()

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig("/home/llr/ilc/oprea/data/resolution_and_linearity_vs_N.png")
    plt.show()

def analyze_N_hits_resolution_and_linearity():
    # Retrieve fit parameters for each fit type for N_hits, including gaussian_2_sigma and rms90
    fit_params_N = {
        'gamma': get_fit_params_for_N_hits('gamma'),
        'log_normal': get_fit_params_for_N_hits('log_normal'),
        'gaussian': get_fit_params_for_N_hits('gaussian'),
        'gaussian_2_sigma': get_fit_params_for_N_hits('gaussian_2_sigma'),
        'rms90': get_fit_params_for_N_hits('rms90')
    }

    # Use the same set of energies for all fit types (assuming keys are the same)
    energies_for_plot_N = sorted(fit_params_N['gamma'].keys())

    # Calculate means, stds, and resolutions for each fit type
    means_N = {}
    stds_N = {}
    resolutions_N = {}
    for fit_type in ['gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', 'rms90']:
        means, stds, resolutions = calculate_resolution_and_linearity_for_N_hits(
            fit_params_N[fit_type], fit_type
        )
        means_N[fit_type] = means[fit_type]
        stds_N[fit_type] = stds[fit_type]
        resolutions_N[fit_type] = resolutions

    # Build lists for plotting
    means_N_list = {ft: [means_N[ft][e] for e in energies_for_plot_N] for ft in ['gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', 'rms90']}
    stds_N_list = {ft: [stds_N[ft][e] for e in energies_for_plot_N] for ft in ['gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', 'rms90']}
    resolutions_N_list = {ft: [resolutions_N[ft][e] for e in energies_for_plot_N] for ft in ['gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', 'rms90']}
    print(f"means N_hits list: {means_N_list}\n")
    print(f"stds N_hits list: {stds_N_list}\n")
    print(f"resolutions N_hits list: {resolutions_N_list}")

    # Plot all fit types together for N_hits
    plot_resolution_and_linearity_vs_N_for_all_fit_types(
        energies_for_plot_N,
        means_N_list,
        stds_N_list,
        resolutions_N_list
    )

    plot_resolution_and_linearity_vs_N_for_one_fit_type(
        energies_for_plot_N,
        means_N['log_normal'],
        stds_N['log_normal'],
        resolutions_N['log_normal']
    )

# --- sum_E code (already written, do not delete) ---

def analyze_sum_E_resolution_and_linearity():
    # Retrieve fit parameters for each fit type for sum_E, including gaussian_2_sigma and rms90
    fit_params_E = {
        'gamma': get_fit_params_for_sum_E('gamma'),
        'log_normal': get_fit_params_for_sum_E('log_normal'),
        'gaussian': get_fit_params_for_sum_E('gaussian'),
        'gaussian_2_sigma': get_fit_params_for_sum_E('gaussian_2_sigma'),
        'rms90': get_fit_params_for_sum_E('rms90')
    }

    # Use the same set of energies for all fit types (assuming keys are the same)
    energies_for_plot_E = sorted(fit_params_E['gamma'].keys())

    # Calculate means, stds, and resolutions for each fit type
    means_E = {}
    stds_E = {}
    resolutions_E = {}
    for fit_type in ['gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', 'rms90']:
        means, stds, resolutions = calculate_resolution_and_linearity_for_sum_E(
            fit_params_E[fit_type], fit_type
        )
        means_E[fit_type] = means[fit_type]
        stds_E[fit_type] = stds[fit_type]
        resolutions_E[fit_type] = resolutions

    # Build lists for plotting
    means_E_list = {ft: [means_E[ft][e] for e in energies_for_plot_E] for ft in ['gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', 'rms90']}
    stds_E_list = {ft: [stds_E[ft][e] for e in energies_for_plot_E] for ft in ['gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', 'rms90']}
    resolutions_E_list = {ft: [resolutions_E[ft][e] for e in energies_for_plot_E] for ft in ['gamma', 'log_normal', 'gaussian', 'gaussian_2_sigma', 'rms90']}
    print(f"means sum_E list: {means_E_list}\n")
    print(f"stds sum_E list: {stds_E_list}\n")
    print(f"resolutions sum_E list: {resolutions_E_list}")

    # Plot all fit types together for sum_E
    plot_resolution_and_linearity_vs_E_for_all_fit_types(
        energies_for_plot_E,
        means_E_list,
        stds_E_list,
        resolutions_E_list
    )

    plot_resolution_and_linearity_vs_E_for_one_fit_type(
        energies_for_plot_E,
        means_E['log_normal'],
        stds_E['log_normal'],
        resolutions_E['log_normal']
    )

analyze_N_hits_resolution_and_linearity()





