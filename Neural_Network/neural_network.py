import pandas as pd
import sys
import keras
import numpy as np
import scipy as sp
import math
import helper

#Import data
df = pd.read_csv('Sample_Data.csv', sep=',')
number_cells = df.shape[0]
print(df.shape)
print(df.head())

# Set size
# train_size = 27848
train_size = 60
train_data = df[0:train_size]
train_data = train_data.fillna(0)

# Training/testing sets
y = train_data['Cell_class']
X = train_data.drop(['Cell_ID','Behavior','Bregma','Centroid_X','Centroid_Y',
    'Cell_class','Neuron_cluster_ID'], axis = 1)

# Define a neuron
class Neuron:

    def __init__(self, activation, bias = None):
        self.b = bias
        self.a = activation
        self.edges = []
        self.neighbors = []

    def update(self):
        helper.sigmoid()

# Define a layer of neurons
class Layer():

    def __init__(self):
        pass

# Define an edge connecting a neuron in one layer to a neuron in the next layer
class Edge():

    def __init__(self, src, target, weight):
        if type(src) != Neuron or type(target) != Neuron:
            print("Edge __init__ parameter type error")
            sys.exit()

        self.src = src
        self.target = target
        self.w = weight

# Define network setup
class Network():

    def __init__(self, data_input):
        self.data = data_input
        self.first_layer_size = data_input.shape[0]
        self.primary_layer = []

    # Train neural network
    def run(self, gtruth):
        pass

###############################################################
# RUN
###############################################################
n = Network(helper.drop_empty_rows(helper.genes_only(df)))