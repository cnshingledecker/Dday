import os, time

num_fitting_factor_combinations = [1,2,5,10,20,40,60,100,200]
startTimesOriginal = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

startTimesParallel = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
]

num_processors = [1, 2, 3, 4]

index = -1
for number in num_fitting_factor_combinations:
    index += 1
    print("Running original script with input size " + str(number))
    file = open("experimental_data/num_fitting_factor_combinations.csv", "w")
    file.write(str(number))
    file.close()
    startTime = int(round(time.time() * 1000))
    os.system("python3 runtime_test_baragiola_optimization_generalization_rmsd_generalization.py")
    endTime = int(round(time.time() * 1000))
    startTimesOriginal[index] = (endTime - startTime)/60000



for num in num_processors:
    num_processors_file = open("experimental_data/num_processors.csv", "w")
    num_processors_file.write(str(num))
    num_processors_file.close()

    index = -1
    for number in num_fitting_factor_combinations:
        index += 1

        file = open("experimental_data/num_fitting_factor_combinations.csv", "w")
        file.write(str(number))
        file.close()

        print("Running tests with " + str(num) + " processors (this is the parallel version)." + ", with " + str(number) + " fitting factor combinations.")
        startTime = int(round(time.time() * 1000))
        os.system("mpiexec -n " + str(num) + " python3 runtime_test_baragiola_generalized_parallel.py")
        endTime = int(round(time.time() * 1000))
        os.system("rm -r baragiola_files_processor*")
        startTimesParallel[num-1][index] = (endTime - startTime)/60000

os.system("./clean.sh")

with open("runtime_comparison.csv", "w") as runtime_csv:
    i = -1
    for element in num_fitting_factor_combinations:
        i += 1
        string_to_write = str(element) + "," + str(startTimesOriginal[i]) + ","
        for j in range(0, len(num_processors) - 1):
            string_to_write += str(startTimesParallel[num_processors[j] - 1][i]) + ", "
        string_to_write += str(startTimesParallel[num_processors[len(num_processors) - 1] - 1][i])
        string_to_write += "\n"
        runtime_csv.write(string_to_write)
