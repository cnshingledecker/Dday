import math

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