#!/usr/bin/env python3
import warnings
warnings.filterwarnings('ignore')

import os 
os.environ['CUDA_VISIBLE_DEVICES'] = '-1' # disable GPU usage
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

import keras
from keras.models import Sequential
from keras.layers import Dense
from keras.regularizers import L1L2
# from scikeras.wrappers import KerasClassifier 
from keras.wrappers.scikit_learn import KerasClassifier
from keras.utils import np_utils

# from tensorflow.config import run_functions_eagerly
# run_functions_eagerly(True)

def build_keras_model(alpha, input_dim, activation=None, optimizer='rmsprop', loss=keras.losses.BinaryCrossentropy(from_logits=True), metrics=['accuracy']):
    
    model = Sequential() 
    
    model.add( Dense(1,
                     activation=activation, 
                     kernel_regularizer=L1L2(l1=alpha, l2=1.0-alpha), 
                     input_dim = input_dim) 
    ) 
    
    # print(model.summary()) 
    
    model.compile(optimizer=optimizer,
                  loss=loss,
                  metrics=metrics
    )
    
    return model

def keras_logreg(alpha=0, activation=None, input_dim=4, optimizer='rmsprop', loss=keras.losses.BinaryCrossentropy(from_logits=True), metrics=['accuracy']):
    
    return KerasClassifier(build_fn=build_keras_model,
                           epochs=10,
                           batch_size=1,
                           verbose=0) 

