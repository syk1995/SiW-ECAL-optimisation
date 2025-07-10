import os
import glob
import math
import uproot
import numpy as np
import tensorflow as tf
from tensorflow.keras.layers import Conv3D, MaxPooling3D, BatchNormalization, Flatten, Dense, Input, Concatenate, Dropout
import re
import time

start_time = time.time()
print(f"Start Time: {start_time:.2f} s")

gpus = tf.config.list_physical_devices('GPU')

if gpus:
    print("Available GPU devices:")
    for gpu in gpus:
        print(gpu)
    print("For train.")
else:
    print("No GPU devices found.")

input_shape_conv = (50, 50, 50, 1)
additional_input_shape = (1,)

model = tf.keras.Sequential([
    Conv3D(64, (5, 5, 5), activation='relu', input_shape=input_shape_conv),
    Conv3D(32, (3, 3, 3), activation='relu'),
    MaxPooling3D(pool_size=(2, 2, 2)),
    
    Conv3D(32, (3, 3, 3), activation='relu'),
    Conv3D(32, (3, 3, 3), activation='relu'),
    BatchNormalization(),
    MaxPooling3D(pool_size=(2, 2, 2)),
    
    Conv3D(32, (3, 3, 3), activation='relu'),
    Conv3D(6, (3, 3, 3), activation='relu'),
    MaxPooling3D(pool_size=(2, 2, 2)),
    
    Flatten(),
])

additional_input = Input(shape=additional_input_shape)
concatenated_output = Concatenate()([model.output, additional_input])
fc_output = Dense(512, activation='relu')(concatenated_output)
fc_output = Dense(128, activation='relu')(fc_output)
fc_output = Dense(32, activation='relu')(fc_output)
fc_output = Dropout(0.2)(fc_output)

fem_prediction = Dense(1, activation='sigmoid', name='fem_output')(fc_output)

# Output layer for energy prediction
# energy_prediction = Dense(1, activation='linear', name='energy_output')(fc_output)


final_model = tf.keras.Model(inputs=[model.input, additional_input], outputs=fem_prediction)

# final_model = tf.keras.Model(inputs=[model.input, additional_input], outputs=[fem_prediction, energy_prediction])

# Print model summary
final_model.summary()

# Compile the model
final_model.compile(optimizer='adam', loss="mean_squared_error", metrics=['mape','mse','mae'])

compile_time = time.time()
print(f"Compile Time: {compile_time - start_time:.2f} s")

# Save the model
# final_model.save('/gpfs/users/xiax/SoftwareCompensation/model/ssh_sigmoid_2epoch_full.h5')

class CustomSequence(tf.keras.utils.Sequence):
    def __init__(self, sample_paths, batch_size):
        self.sample_paths = sample_paths
        self.batch_size = batch_size

    def __len__(self):
        return math.ceil(len(self.sample_paths) / self.batch_size)

    def __getitem__(self, idx):
        low = idx * self.batch_size
        high = min(low + self.batch_size, len(self.sample_paths))
        sample_batch = self.sample_paths[low:high]

        x_shower_batch = []
        x_energy_batch = []
        y_fem_batch = []
        # y_energy_batch = []

        for sample in sample_batch:
            with uproot.open(sample["file_path"]) as f:
                fem_values = f["fem_evt"].values()
                energy_values  = f['energy_evt'].values()
                hist = f[f'pos/pos_{sample["hist_index"]}']
                
            max_value = np.max(hist.values())

            if max_value != 0:
                fem = round(fem_values[sample["hist_index"]], 3)
                energy = round(energy_values[sample["hist_index"]], 3)
                normalized_hist = np.array(hist.values() / max_value, dtype=np.float32)
            else:
                fem = 0
                energy = 0
                normalized_hist = np.zeros(hist.values().shape, dtype=np.float32)
            
            # pattern = r'\d+\.\d+'
            # match = re.search(pattern, sample["file_path"])

            x_shower_batch.append(normalized_hist)
            x_energy_batch.append(energy)
            y_fem_batch.append(fem)
            # y_energy_batch.append(float(match.group()))


        # Convert lists to numpy arrays
        x_shower_batch = np.array(x_shower_batch, dtype=np.float32)
        x_energy_batch = np.array(x_energy_batch, dtype=np.float32)
        y_fem_batch = np.array(y_fem_batch, dtype=np.float32)
        # y_energy_batch = np.array(y_energy_batch, dtype=np.float32)
        
        # Construct x_batch as a tuple
        x_batch = (x_shower_batch, x_energy_batch)
        # y_batch = (y_fem_batch, y_energy_batch)
        y_batch = y_fem_batch
        
        return x_batch, y_batch

# from tensorflow.keras.callbacks import Callback
batch_size = 256
file_pattern = "/gpfs/workdir/xiax/ready/*_ready.root"

file_list = glob.glob(file_pattern)[:10]
# print(file_list)
expanded_file_list = []
for file in file_list:
    for index in range(1000):
        expanded_file_list.append({"file_path": file,"hist_index": index})

# Split your file list into training and validation sets
split_index = int(len(expanded_file_list) * 0.7)
print(f"Split Index: {split_index}")
train_files = expanded_file_list[:split_index]
val_files = expanded_file_list[split_index:]

# Create instances of the CustomSequence class for training and validation
train_sequence = CustomSequence(train_files, batch_size)
val_sequence = CustomSequence(val_files, batch_size)

load_time = time.time()
print(f"Load Time: {load_time - compile_time:.2f}s")

# Train your model using fit method and include the callback
history = final_model.fit(train_sequence,
                epochs=1,
                validation_data=val_sequence)

import matplotlib.pyplot as plt
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.legend(['train','val'])
plt.savefig('/gpfs/users/xiax/SoftwareCompensation/result/loss_function.pdf')

fit_time = time.time()
print(f"Fit Time: {fit_time - load_time:.2f}s")

model_name = '/gpfs/users/xiax/SoftwareCompensation/model/mse_ssh_sigmoid_1epoch_0_900.h5'
final_model.save(model_name)
print(f"Model saved as: {model_name}")

total_time = time.time()
print(f"Save Time: {total_time - fit_time :.2f} s")