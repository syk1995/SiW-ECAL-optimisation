import os
import glob
import math
import uproot
import numpy as np
import tensorflow as tf
from tensorflow.keras.layers import Conv3D, MaxPooling3D, BatchNormalization, Flatten, Dense, Input, Concatenate, Dropout
import re
from tensorflow.keras.models import load_model
import matplotlib.pyplot as plt

# file_pattern = "/gpfs/workdir/xiax/old_data/*_ready.root"
file_pattern = "/gpfs/workdir/xiax/ready/*_ready.root"
x_index = []
y_cnn_output_fem = []
y_ccn_output_energy = []
y_actual_fem = []
y_delta_fem = []
energies_meas = []
energies_true = []
y_delta_energy = []
input_shape_conv = (50, 50, 50, 1)
# print(glob.glob(file_pattern))
# pattern = r'\d+\.\d+'
# for it, file_path in enumerate(glob.glob(file_pattern)):
#     match = float(re.search(pattern, file_path).group())
#     if 9.5 < match < 10.5:
#         print(it, file_path)

file_index = 937 #454 10GeV #14 145gev 

print(glob.glob(file_pattern)[file_index:file_index+1][0])
for index in range(1000):
    file_path = glob.glob(file_pattern)[file_index:file_index+1][0]
    print(file_path)
    with uproot.open(file_path) as f:
        hist = f[f'pos/pos_{index}']
        fem_value = f["fem_evt"].values()[index]
        energy = f['energy_evt'].values()[index]

    energies_meas.append(energy)
    pattern = r'\d+\.\d+'
    match = re.search(pattern, file_path)
    # print(match.group())
    energies_true.append(float(match.group()))
        

    # Preprocess the histogram values by normalizing them
    max_value = np.max(hist.values())
    if max_value != 0:
        normalized_hist = np.array(hist.values() / max_value, dtype=np.float32)
    else:
        normalized_hist = np.zeros(hist.values().shape, dtype=np.float32)

    # Prepare the input data for the model
    input_shower = np.array([normalized_hist])
    input_shower = input_shower.reshape((1,) + input_shape_conv)

    input_energy = energy
    # Use the model's predict() method to obtain the predicted output
    model = load_model('/gpfs/users/xiax/SoftwareCompensation/model/mse_ssh_sigmoid_1epoch_000_450.h5')
    predicted_output = model.predict([input_shower, np.array([input_energy])])

    # Compare the predicted output with the actual number
    # print(f"Index {index}, Actual Number: {fem_value}, Predicted Number: {predicted_output[0][0]}")
    x_index.append(index)
    y_cnn_output_fem.append(predicted_output[0][0])
    # y_ccn_output_energy.append(predicted_output[1][0])
    y_actual_fem.append(fem_value)
    y_delta_fem.append(abs(fem_value - predicted_output[0][0]))
    # y_delta_energy.append(abs(energy - predicted_output[1][0][0]))
    # print(predicted_output)
    # print(predicted_output[0][0], predicted_output[1][0], fem_value, energy)

plt.plot(x_index, y_actual_fem, label='Actual Fem')
plt.plot(x_index, y_cnn_output_fem, label='Predicted Fem')

plt.xlabel('Index')
plt.ylabel('Fem')
plt.legend()
# plt.show()
plt.savefig("/gpfs/users/xiax/SoftwareCompensation/result/Fem_1epoch_000_450.pdf", format="pdf")

plt.plot(x_index, y_delta_fem, label='Difference')
plt.xlabel('Index')
plt.ylabel('Difference')
plt.legend()
# plt.show()
plt.savefig("/gpfs/users/xiax/SoftwareCompensation/result/Difference_1epoch_000_450.pdf", format="pdf")

plt.hist(y_delta_fem, bins=100)
plt.xlabel('Difference')
plt.ylabel('Frequency')
# plt.show()
plt.savefig("/gpfs/users/xiax/SoftwareCompensation/result/Difference_Frequency_1epoch_000_450.pdf", format="pdf")
