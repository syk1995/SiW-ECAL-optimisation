import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

def mean_hits_from_file(file_path):
    with uproot.open(file_path) as file:
        tree = file["events"]
        # List all branches
        branches = tree.keys()
        if "simplecaloRO" not in branches:
            raise ValueError(f"'simplecaloRO' branch not found in {file_path}")

        try:
            hits = tree.arrays(filter_name="simplecaloRO/*", library="ak")

        except Exception:
            hits = tree.array("simplecaloRO", library="ak")

        # Count hits: usually count by 'cellID' if available
        if "simplecaloRO.cellID" in hits.fields:
            n_hits = ak.count(hits["simplecaloRO.cellID"], axis=1)
        elif "cellID" in hits.fields:
            n_hits = ak.count(hits["cellID"], axis=1)
        else:
            # If no cellID, count by any field length (e.g. energy)
            first_field = hits.fields[0]
            n_hits = ak.count(hits[first_field], axis=1)

        mean_hits = ak.mean(n_hits)
        return mean_hits

def get_mean_hits_vs_energy(energies, file_pattern):
    mean_hits = []
    for E in energies:
        file_path = file_pattern.format(E)
        try:
            mean_hits_val = mean_hits_from_file(file_path)
            mean_hits.append(mean_hits_val)
            print(f"Energy {E} GeV: Mean hits = {mean_hits_val:.2f}")
        except Exception as e:
            print(f"Failed to process {file_path}: {e}")
            mean_hits.append(np.nan)
    return np.array(energies), np.array(mean_hits)

def add_threshold(hits, threshold):
    if 'energy' not in hits.fields:
       raise ValueError('Hits array does not have an energy field')
    mask = hits['energy'] > threshold
    filtered_hits = hits[mask]
    return filtered_hits

def load_hits(file_path):
    with uproot.open(file_path) as file:
        tree = file["events"]
        try:
            hits = tree.arrays(filter_name="simplecaloRO/*", library="ak")
        except Exception:
            hits = tree.array("simplecaloRO", library="ak")
    return hits


if __name__ == "__main__":
    energies = [0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 60.0]
    file_pattern = "/grid_mnt/vol_home/llr/ilc/shi/data/SiWECAL-Prototype/Simu2025-06/CONF0/e-/MC/{}GeV.root"

    energies, mean_hits = get_mean_hits_vs_energy(energies, file_pattern)
    
    valid = ~np.isnan(mean_hits)
    if np.count_nonzero(valid) < 2:
        print("Not enough valid points for fitting.")
    else:
        spline = UnivariateSpline(energies[valid], mean_hits[valid], s=0)

        # Plot calibration points and spline curve
        energies_range = np.linspace(np.min(energies[valid]), np.max(energies[valid]), 500)
        spline_hits = spline(energies_range)

        plt.figure(figsize=(8,5))
        plt.plot(energies[valid], mean_hits[valid], 'o', label='Calibration points')
        plt.plot(energies_range, spline_hits, '-', label='Spline fit E(N)')
        plt.xlabel('Energy (GeV)')
        plt.ylabel('Number of hits')
        plt.title('Energy reconstruction from number of hits')
        plt.grid(True)
        plt.legend()
        plt.show()





        file_path = "/grid_mnt/vol_home/llr/ilc/shi/data/SiWECAL-Prototype/Simu2025-06/CONF0/e-/MC/10.0GeV.root"
        with uproot.open(file_path) as file:
            tree = file["events"]
            # Read the first event's simplecaloRO hits
            hits = tree.arrays(filter_name="simplecaloRO", entry_start=0, entry_stop=1, library="ak")
            simplecaloRO_hits = hits["simplecaloRO"][0]
            
            print("Fields of hits in first event:", simplecaloRO_hits.fields)
            # Print a few hits to inspect energy values
            print(simplecaloRO_hits[:5])




        # Define energy reconstruction function using spline
        def reco_energy(n_hits):
            return spline(n_hits)

        print(f"Reconstructed energy for 100 hits: {reco_energy(100):.2f} GeV")
        hits_10GeV = load_hits(file_pattern.format(10.0))
        print(f"Number of events at 10 GeV: {len(hits_10GeV['simplecaloRO.cellID'])}")

        # Filter hits in the first event only
        hits_first_event = hits_10GeV[0]
        print(f"Number of hits in first event: {len(hits_first_event)}")
        print(f"First 5 hits in first event: {hits_first_event[:5]}")
        filtered_hits = add_threshold(hits_first_event, threshold=0.0005)
        print(f"Filtered hits in first event (energy > 0.0005): {filtered_hits}")

#ask about the dotted splines
#ask about the mip
#ask why matplotlib.pyplot does not work
