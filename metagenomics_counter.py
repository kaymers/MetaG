import json

import numpy as np
import pyVectorizer
import pandas as pd


def count_process(input_path, k, output_path):
    vectors = np.array(pyVectorizer.vectorize_file(input_path, k)).astype(np.float32)
    acgt = pyVectorizer.count_acgt_file(input_path)

    with open(output_path/"_acgt.json", 'w') as outfile:
        json.dump(acgt, outfile)


def freq_process(input_path, k, output_path):
    vectors = np.array(pyVectorizer.vectorize_file(input_path, k)).astype(np.float32)

    vectorSum = np.zeros(len(pyVectorizer.all_kmers(k)),dtype=np.int64)

    for i in vectors:
        vectorSum += np.array(i).astype(np.int64)
    file_name = str(k) + "_freqs.json"
    with open(output_path/file_name, 'w') as outfile:
        json.dump((k,vectorSum.tolist()), outfile)
