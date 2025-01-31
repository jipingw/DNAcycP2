import keras
import pandas as pd
import numpy as np
from numpy import array
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed
import tensorflow as tf
import os
from typing import List, Tuple

detrend_int_original = 0.029905550181865692
detrend_slope_original = 0.973293125629425
normal_mean_original = -0.18574825868055558
normal_std_original = 0.4879013326394626

detrend_int_smooth = 0.001641373848542571
detrend_slope_smooth = 1.0158132314682007
# Mean and stdev of smoothed C0 for Tiling library:
# (calculated from/on the scale of normalized Cn values)
normal_mean_smooth = -0.011196041799376931
normal_std_smooth = 0.651684644408004

def dnaOneHot(sequence):
    code = np.array([
        [1, 0, 0, 0],  # A / a
        [0, 1, 0, 0],  # C / c
        [0, 0, 1, 0],  # G / g
        [0, 0, 0, 1],  # T / t
        [0, 0, 0, 0]   # N / n
    ])
    
    mapping = np.zeros(128, dtype=int)
    mapping[ord('A')] = 0
    mapping[ord('C')] = 1
    mapping[ord('G')] = 2
    mapping[ord('T')] = 3
    mapping[ord('N')] = 4
    mapping[ord('a')] = 0
    mapping[ord('c')] = 1
    mapping[ord('g')] = 2
    mapping[ord('t')] = 3
    mapping[ord('n')] = 4

    indices = np.fromiter((mapping[ord(char)] for char in sequence), dtype=int)
    onehot_encoded_seq = code[indices]

    return onehot_encoded_seq

def construct_predict_fn(model):
    @tf.function(reduce_retracing=True)
    def predict_fn(x):
        return model(x, training=False)
    return predict_fn

def cycle_fasta(inputfile, folder_path, chunk_size, num_threads):
    network_final = keras.models.load_model(folder_path)
    smooth = True if "smooth" in folder_path else False
    genome_file = SeqIO.parse(open(inputfile),'fasta')

    predict_fn = construct_predict_fn(network_final)

    def process_chunk(args) -> Tuple[np.ndarray, np.ndarray]:
        """
        Process a chunk of the sequence data with both forward and reverse predictions.
        """
        onehot_sequence, start_ind, end_ind, network_final = args
        
        # Extract local sequence window
        ind_local = np.arange(start_ind, end_ind)
        onehot_sequence_local = onehot_sequence[np.arange(ind_local[0] - 25, ind_local[-1] - 24)[:, None] + np.arange(50)]
        onehot_sequence_local = onehot_sequence_local.reshape((-1, 50, 4, 1))
        
        # Create reverse sequence
        onehot_sequence_local_reverse = np.flip(onehot_sequence_local, [1, 2])
        
        # Convert to TensorFlow tensors with fixed shape
        onehot_sequence_local = tf.cast(onehot_sequence_local, tf.float32)
        onehot_sequence_local_reverse = tf.cast(onehot_sequence_local_reverse, tf.float32)
        
        # Make predictions using the optimized prediction function
        fit_local = predict_fn(onehot_sequence_local).numpy().reshape(-1)
        fit_local_reverse = predict_fn(onehot_sequence_local_reverse).numpy().reshape(-1)
        
        return fit_local, fit_local_reverse
    
    ret = {}

    for fasta in genome_file:
        chrom = fasta.id
        genome_sequence = str(fasta.seq)
        print(f"Sequence length for ID {chrom}: {len(genome_sequence)}", flush=True)
        onehot_sequence = dnaOneHot(genome_sequence)
        onehot_sequence = array(onehot_sequence)
        onehot_sequence = onehot_sequence.reshape((onehot_sequence.shape[0],4,1))
        print("Predicting cyclizability...", flush=True)

        print(f"Chunk size: {chunk_size}, num threads: {num_threads}", flush=True)

        sequence_length = onehot_sequence.shape[0] - 49
        start_indices = range(25, sequence_length + 25, chunk_size)
        end_indices = [min(start + chunk_size, sequence_length + 25) for start in start_indices]

        chunk_args = [
            (onehot_sequence, start, end, predict_fn) 
            for start, end in zip(start_indices, end_indices)
        ]

        fit = []
        fit_reverse = []

        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            for i, (fit_local, fit_local_reverse) in enumerate(executor.map(process_chunk, chunk_args)):
                fit.append(fit_local)
                fit_reverse.append(fit_local_reverse)
                
                # Progress reporting for long sequences
                total_processed = sum(len(chunk) for chunk in fit)
                if onehot_sequence.shape[0] > 10**7:
                    print(f"\t Completed predictions on {total_processed} out of {sequence_length} sequences", flush=True)

        fit = np.concatenate(fit)
        fit_reverse = np.concatenate(fit_reverse)
        if smooth:
            fit = detrend_int_smooth + (fit + fit_reverse) * detrend_slope_smooth / 2
            fit2 = fit * normal_std_smooth + normal_mean_smooth
        else:
            fit = detrend_int_original + (fit + fit_reverse) * detrend_slope_original / 2
            fit2 = fit * normal_std_original + normal_mean_original
        n = fit.shape[0]
        positions = np.arange(25, 25 + n)
        fitall = np.column_stack((positions, fit, fit2))
        if smooth:
            fitall = pd.DataFrame(fitall, columns=["position", "C0S_norm", "C0S_unnorm"])
        else:
            fitall = pd.DataFrame(fitall, columns=["position", "C0_norm", "C0_unnorm"])
        fitall = fitall.astype({"position": int})
        ret[f"cycle_{chrom}"] = fitall
    
    return ret

def cycle(sequences, folder_path):
    if type(sequences) == str:
        sequences = [sequences]
    network_final = keras.models.load_model(folder_path)
    smooth = True if "smooth" in folder_path else False
    X = []
    all50 = True
    print("Reading sequences...")
    for sequence_nt in sequences:
        if len(sequence_nt) != 50:
            all50=False
        X.append(dnaOneHot(sequence_nt))

    if all50:
        print("Predicting cyclizability...")
        X = array(X)
        X = X.reshape((X.shape[0],50,4,1))
        X_reverse = np.flip(X,[1,2])

        model_pred = network_final.predict(X)
        model_pred_reverse = network_final.predict(X_reverse)

        if smooth:
            model_pred = detrend_int_smooth + (model_pred + model_pred_reverse) * detrend_slope_smooth / 2
        else:
            model_pred = detrend_int_original + (model_pred + model_pred_reverse) * detrend_slope_original / 2
        output_cycle = model_pred.flatten()
        if smooth:
            output_cycle2 = [item * normal_std_smooth + normal_mean_smooth for item in output_cycle]
        else:
            output_cycle2 = [item * normal_std_original + normal_mean_original for item in output_cycle]
    else:
        print("Not all sequences are length 50, predicting every subsequence...")
        output_cycle = []
        lenX = len(X)
        for j, onehot_loop in enumerate(X):
            l = len(onehot_loop)
            onehot_loop = array(onehot_loop)
            onehot_loop = onehot_loop.reshape((l,4,1))
            onehot_loops = []
            for i in range(l-49):
                onehot_loops.append(onehot_loop[i:i+50])
            onehot_loops = array(onehot_loops)
            onehot_loops_reverse = np.flip(onehot_loops,[1,2])
            if l > 1000:
                # Provide status bar for long sequences:
                cycle_local = network_final.predict(onehot_loops)
                cycle_local_reverse = network_final.predict(onehot_loops_reverse)
            else:
                # No status bar for short sequences (verbose=0):
                cycle_local = network_final.predict(onehot_loops, verbose=0)
                cycle_local_reverse = network_final.predict(onehot_loops_reverse, verbose=0)
            if smooth:
                cycle_local = detrend_int_smooth + (cycle_local + cycle_local_reverse) * detrend_slope_smooth/2
            else:
                cycle_local = detrend_int_original + (cycle_local + cycle_local_reverse) * detrend_slope_original/2
            cycle_local = cycle_local.reshape(cycle_local.shape[0])
            output_cycle.append(cycle_local)
            if j%10==9:
                print(f"Completed {j+1} out of {lenX} total sequences")
        if smooth:
            output_cycle2 = [item * normal_std_smooth + normal_mean_smooth for item in output_cycle]
        else:
            output_cycle2 = [item * normal_std_original + normal_mean_original for item in output_cycle]

    if smooth:
        ret = []
        for i in range(len(output_cycle)):
            cur_length = 1 if all50 else len(output_cycle[i])
            df = pd.DataFrame({
                "position": np.arange(25, 25 + cur_length),
                "C0S_norm": output_cycle[i],
                "C0S_unnorm": output_cycle2[i]
            })
            ret.append(df)
    else:
        ret = []
        for i in range(len(output_cycle)):
            cur_length = 1 if all50 else len(output_cycle[i])
            df = pd.DataFrame({
                "position": np.arange(25, 25 + cur_length),
                "C0_norm": output_cycle[i],
                "C0_unnorm": output_cycle2[i]
            })
            ret.append(df)
    ret2 = [ret, output_cycle, output_cycle2]
    return ret2
