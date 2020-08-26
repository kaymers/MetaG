import keras
from keras import Sequential
from keras.layers import Dense
import numpy as np
import pyVectorizer
from keras import backend as K
import pandas as pd
from os import path

def normalize_over_axis1(arr2d):
    for i in range(arr2d.shape[0]):
        hori_sum = sum(arr2d[i])
        if hori_sum!= 0 :arr2d[i] = (arr2d[i])/hori_sum

def build_model(layers_sizes, activation ):
    model = Sequential()
    print(layers_sizes[1], layers_sizes[0])
    model.add(Dense(layers_sizes[1], input_dim=layers_sizes[0],activation = activation))
    # d = Dense(64, input_dim=136, activation="tanh")
    # model.add(d)
    # model.add(Dense(64, input_dim=136,activation = "tanh"))

    for i in range(2,len(layers_sizes)):
        model.add(Dense(layers_sizes[i], activation = activation))
    model.compile(loss='mse', optimizer=keras.optimizers.Adam())
    return model

def train(model,epochs,vectors):
    model.fit(vectors,vectors,epochs=epochs)

def get_points(model, vectors):
    middle_layer=model.layers[int(len(model.layers) / 2 - 1)]
    input_layer=model.layers[0]
    mapping_function = K.function([input_layer.input], [middle_layer.output])
    points = mapping_function(vectors)[0]
    return points

def process(input_path, k,epochs, activation, layers_sizes, output_path):
    model = build_model(layers_sizes,activation)
    vectors = np.array(pyVectorizer.vectorize_file(input_path, k)).astype(np.float32)
    normalize_over_axis1(vectors)
    train(model,epochs,vectors)
    r=get_points(model,vectors)
    print(r)
    dataset = pd.DataFrame({'x': r[:, 0], 'y': r[:, 1]}) if (r.shape[1] == 2) else pd.DataFrame({'x': r[:, 0], 'y': r[:, 1], 'z': r[:, 2]})

    if output_path: dataset.to_csv(output_path)

    x = layers_sizes[len(layers_sizes) // 2]
    layers_sizes[len(layers_sizes) // 2] = 'x'

    # Saving settings
    settings_path = 'settings/sample_settings.csv'
    input_path = str(input_path)
    settings = {'FASTA File': [input_path.split("/")[-1].split('\\')[-1]],
                'Iterations': [epochs],
                'k in kmer': [k],
                'layers_sizes': str(layers_sizes).replace("'", "").replace(" ", ""),
                '2D': 'Yes' if x == 2 else "",
                '3D': 'Yes' if x == 3 else "",
                'Activation': [activation]
                }

    pd.DataFrame(settings).to_csv(settings_path, mode='a', header=not (path.exists(settings_path)), index=False)
