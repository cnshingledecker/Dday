import os, time

num_fitting_factor_combinations = [1,2,4,8,12,16,20,40,60,100,200,300,400]
startTimesOriginal = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
startTimesParallel = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

index = -1
for number in num_fitting_factor_combinations:
    index += 1

    file = open("experimental_data/num_fitting_factor_combinations.csv", "w")
    file.write(str(number))
    file.close()

    print("Running original script with input size " + str(number))
    startTime = int(round(time.time() * 1000))
    os.system("python3 baragiola_optimization_generalization_rmsd_generalization.py")
    endTime = int(round(time.time() * 1000))
    startTimesOriginal[index] = (endTime - startTime)/1000

    print("Running original script with input size " + str(number))
    startTime = int(round(time.time() * 1000))
    os.system("mpiexec -n 4 python3 baragiola_generalized_parallel.py")
    endTime = int(round(time.time() * 1000))
    os.system("rm -r baragiola_files_processor*")
    startTimesParallel[index] = (endTime - startTime)/1000

os.system("./clean.sh")

with open("runtime_comparison.csv", "w") as runtime_csv:
    i = -1
    for element in num_fitting_factor_combinations:
        i += 1
        runtime_csv.write(str(element) + "," + str(startTimesOriginal[i]) + "," + str(startTimesParallel[i]) + "\n")
