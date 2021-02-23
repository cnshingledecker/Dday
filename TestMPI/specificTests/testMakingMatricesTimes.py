# This file tests: the time it takes (collectively and on average per matrix) to make 100 matrices of various sizes using the np.random.rand function

import numpy as np
import time, os

array_num_rows_and_num_cols = [32,64,128,256,512,1024,2048] # Each element is the size of a matrix (it is N, where a matrix that is N x N is made)
times_ms = [0,0,0,0,0,0,0]
for i in range(0,len(array_num_rows_and_num_cols)):
    length = array_num_rows_and_num_cols[i]
    for j in range(0, 100):
        ts = time.time()
        array = np.random.rand(length, length)
        te = time.time()
        times_ms[i] += (te-ts)*1000

print(times_ms)
for i in range(0, len(times_ms)):
    print("Matrix size: " + str(array_num_rows_and_num_cols[i]) +" x " + str(array_num_rows_and_num_cols[i]) + " "*(25-len("Matrix size: " + str(array_num_rows_and_num_cols[i]) +" x " + str(array_num_rows_and_num_cols[i]))) + " It took " + str(float(times_ms[i]/100)) + " "*(20-len(str(float(times_ms[i]/100)))) + " ms to make 100 matrices of this size, for an average time of " + str(times_ms[i]/(100*array_num_rows_and_num_cols[i])) + " "*(22-len(str(times_ms[i]/(100*array_num_rows_and_num_cols[i])))) +  " ms per matrix.")
