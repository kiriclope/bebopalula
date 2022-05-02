import numpy as np
from tensorflow import keras
from keras.callbacks import ModelCheckpoint
from matplotlib import pyplot as plt
from keras import backend as K

import sys
sys.path.insert(0, '../') 
import utils.constants as gv 
from utils.options import *

import utils.get_data as data
from utils.get_days import * 

import utils.preprocessing as pp
import utils.plot_utils as pl 

# create the dataset
data_n = 100000
timepoints = 40
sign_vec = (np.random.randint(2, size=(data_n, 1)) * 2) - 1

input = np.random.rand(data_n, timepoints) - 0.5  # random values centered around 0
input = input + (sign_vec * 0.2)
input = input[:, :, np.newaxis]  # reshape for RNN
output = sign_vec * 3

# import the dataset
data.get_days() # do not delete that !!

X_S1, X_S2 = get_X_S1_X_S2_days_task(day=options['day'],
                                     stimulus=options['stimulus'],
                                     task=options['task'],
                                     trials=options['trials']) 

X_S1, X_S2 = pp.preprocess_X_S1_X_S2(X_S1, X_S2,
                                     scaler=options['scaler_BL'],
                                     center=options['center_BL'],
                                     scale=options['scale_BL'],
                                     avg_mean=options['avg_mean_BL'],
                                     avg_noise=options['avg_noise_BL'],
                                     unit_var=options['unit_var']) 

input = np.vstack(X_S1, X_S2)
input = input[:, :, np.newaxis]  # reshape for RNN

timepoints = input.shape[-1]

# plot the dataset
plt.figure()
for i in range(100):
   c = 'k' if sign_vec[i] == 1 else 'r'
   plt.plot(input[i, :], c)

plt.show()

# build model
units = 50
model = keras.models.Sequential()
# have to return sequences to get the actual unrolled response across time
model.add(keras.layers.SimpleRNN(units, return_sequences=True, input_shape=(timepoints, 1)))
model.add(keras.layers.Dense(1))
model.summary()

# set a history callback to save weights after each epoch
weights_file_prefix = "tmp/model_"
model.save(weights_file_prefix + "untrained.h5")
checkpoint = ModelCheckpoint(weights_file_prefix + "{epoch}.h5", save_freq='epoch')

# train
epochs = 10
model.compile(
   loss="mse",
   optimizer="sgd",
   metrics=["mse"]
)

history = model.fit(
   input, output, batch_size=1000, epochs=epochs, callbacks=[checkpoint]
)

# try a prediction on the trained network
n_examples = 100
neuron = 0

plt.figure()
for i in range(n_examples):
   input_test = input[i][np.newaxis, :, :]
   y_pred = model.predict(input_test)

   # get activations for rnn layer
   out_func = K.function([model.input], [model.layers[0].output])
   out_vals = out_func([input_test])
   out_activation = out_vals[0][0, :, :]

   c = 'k' if sign_vec[i] == 1 else 'r'
   plt.plot(out_activation[:, neuron], c)
plt.ylim(-1, 1)
plt.show()# repeat for the untrained network
model.load_weights(weights_file_prefix + "untrained.h5")
plt.figure()
for i in range(n_examples):
   input_test = input[i][np.newaxis, :, :]
   y_pred = model.predict(input_test)

   # get activations for rnn layer
   out_func = K.function([model.input], [model.layers[0].output])
   out_vals = out_func([input_test])
   out_activation = out_vals[0][0, :, :]

   c = 'k' if sign_vec[i] == 1 else 'r'
   plt.plot(out_activation[:, neuron], c)
plt.ylim(-1, 1)
plt.show()

