# len_data = 78

import numpy as np
import itertools, math
# num_chunks = 4
# num_chunks_with_extra = (len_data / num_chunks - int(len_data / num_chunks)) * num_chunks
# print("Num chunks with extra is " + str(num_chunks_with_extra))
# minimum_chunk_size = int(len_data/num_chunks)
# start = 0
# end = 0
# for i in range(0, num_chunks):
#     print("\n\n")
#     length = minimum_chunk_size
#     if(i < num_chunks_with_extra):
#         length += 1
#     end = start + length - 1
#     # print("I is " + str(i) + "Start is " + str(start) + ", End is " + str(end))
#         # print("Inside if; i is " + str(i) + "Start is " + str(start) + ", End is " + str(end))
#     print("Chunk " + str(i) + " gets elements start:end --> " + str(start) + ":" + str(end))
#     start = end+1


# data = [list(np.linspace(1.7,2.7,1)), list(np.linspace(1,2,34)), list(np.linspace(0.3,0.35,1))]
# j = [list(np.linspace(1.7,2.7,1)), list(np.linspace(1,2,34)), list(np.linspace(0.3,0.35,1))]
# data = itertools.product(*data)

# k = -1
# diff = 0
# print(type(j[1][3]))
# for i in data:
#     k += 1
#     diff += (j[1][k] - i[1])
#     print(type(i[1]))
#     print(type(i))
#     # print("Difference is " + str(j[0][0] - i[0]) + " for 0, " + str(j[1][k] - i[1]) + " for 1, and " + str(j[2][0] - i[2]) + " for 2")

# print("Total difference is " + str(diff))

i = [
        [
            [10,1,-7],
            [3,2,5],
            [11,1,1],
            [3,5,8]
        ],
        [
            [2,6,1],
            [4,5,1],
            [4,2,4],
            [3,7,0]
        ]
    ]

fitting_factors_and_least_rmsd = [1e80,0,0,0]
for mini_chunk in i:
    for combination in mini_chunk:
        # Do work of inserting new parameters; running the run.sh file for this processor (includes monaco), and calculation of a performance measure
        #
        #
            
        fake_performance_measure = math.sqrt(sum([fitting_factor**2 for fitting_factor in combination]))
        if (fake_performance_measure < fitting_factors_and_least_rmsd[0]):
            fitting_factors_and_least_rmsd[0] = fake_performance_measure
            fitting_factors_and_least_rmsd[1] = combination[0]
            fitting_factors_and_least_rmsd[2] = combination[1]
            fitting_factors_and_least_rmsd[3] = combination[2]

print(fitting_factors_and_least_rmsd)

j = [
    [math.sqrt(9+25+64), 3,5,8],
    [math.sqrt(41), 2,6,1],
    [6.0,4,2,4],
    [math.sqrt(9+49), 3,7,0]
]

print("\n")
print("Before sorting: ")
print(j)
print("\n")

j.sort()

print("After sorting: ")
print(j)
print("\n")

print("Smallest performance metric value and associated fitting factors: ")
print(j[0])