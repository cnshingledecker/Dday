# This file tests: Finding the fitting factor combination with the least value of the fake performance metric,
#                  the accuracy of calculating a fake performance measure (test is printing it out to console), 
#                  and sorting a list of n-element lists (n > 1)
import math

i = [   # I is essentially a list of 8 lists that is as if it had been split up into 2 chunks of 4 elements each
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

fitting_factors_and_least_rmsd = [1e80,0,0,0] # 1e80 is a flag value (the fake_performance measure will be less than it)
for mini_chunk in i: 
    for combination in mini_chunk: # Calculate fake_performance measure and save it in the fitting_factors_and_least_rmsd list       
        fake_performance_measure = math.sqrt(sum([fitting_factor**2 for fitting_factor in combination]))
        if (fake_performance_measure < fitting_factors_and_least_rmsd[0]):
            fitting_factors_and_least_rmsd[0] = fake_performance_measure
            fitting_factors_and_least_rmsd[1] = combination[0]
            fitting_factors_and_least_rmsd[2] = combination[1]
            fitting_factors_and_least_rmsd[3] = combination[2]

print("Fitting factors and associated (is smallest) performance metric: " + str(fitting_factors_and_least_rmsd))

j = [ # An array to test how the sort function works with sorting lists with more than 1 element
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