import time
import os
from exportable_custom_functions import is_float, is_int

# startTime1 = int(round(time.time() * 1000))
# os.system("python3 baragiola_optimization.py")
# endTime1 = int(round(time.time() * 1000))

startTime2 = int(round(time.time() * 1000))
os.system("python3 baragiola_optimization_generalization_rmsd_generalization.py") # Run the generalized baragiola optimization
endTime2 = int(round(time.time() * 1000))

startTime3 = int(round(time.time() * 1000))
os.system("python3 baragiola_optimization_test.py") # Run the original baragiola optimization file (the modified copy used for testing)
endTime3 = int(round(time.time() * 1000))


# print("Original (baragiola_optimization.py) took " + str((endTime1 - startTime1) / 1000) + " seconds")
print("Generalization (baragiola_optimization_generalization_rmsd_generalization.py) took " + str((endTime2 - startTime2) / 1000) + " seconds")
print("Modified original (baragiola_optimization_test.py) took " + str((endTime3 - startTime3) / 1000) + " seconds")


# Compare the results from the original optimization and the generalized optimization
num_errors = 0

# Compare the output in the results files; note: gen is short for generalization
results_original = open('resultsFile_3', 'r')
results_gen_rmsd_gen = open('resultsFile_2', 'r')
results_original_lines = results_original.readlines()
results_gen_rmsd_gen_lines = results_gen_rmsd_gen.readlines()

if(len(results_original_lines) == len(results_gen_rmsd_gen_lines)):
    for i in range (0, len(results_original_lines)): # It doesn't matter which list length we use, because if we made it inside this if statement, the lengths are the same
        if(len(results_original_lines[i].strip()) > 0 and len(results_gen_rmsd_gen_lines[i].strip()) > 0): # Skips over lines that look empty in the files, but are really whitespace (the line consists of a positive number of whitespaces)
            results_original_lines_val = results_original_lines[i][0:23].rstrip()
            results_gen_rmsd_gen_lines_val = results_original_lines[i][0:23].rstrip()
            if(results_original_lines_val != results_gen_rmsd_gen_lines_val): # If the fitting factor or the RMSD is not the same
                num_errors += 1
                print("Error: the values on line " + str(i) + " do not match. They are " + results_original_lines_val + " and " + results_gen_rmsd_gen_lines_val)  
else:
    num_errors += 1
    print("Error: A different number of fitting factor combinations were processed.")

results_original.close()
results_gen_rmsd_gen.close()

# Compare the results in the lines of the photo_processes files produced
photo_processes_original = open('photo_processes_3.dat', 'r')
photo_processes_gen_rmsd_gen = open('photo_processes_2.dat', 'r')
photo_processes_original_lines = photo_processes_original.readlines()
photo_processes_gen_rmsd_gen_lines = photo_processes_gen_rmsd_gen.readlines()

# Each of the following two lines is a boolean array that holds (for each line in the respective photo_processes file) whether or not the line is empty
line_contains_text_photo_processes_original = []        
line_contains_text_photo_processes_gen_rmsd_gen = []     

num_non_empty_lines_photo_processes_original = 0
for i in range(0, len(photo_processes_original_lines)):
    if(len(photo_processes_original_lines[i].strip()) > 1):
        num_non_empty_lines_photo_processes_original += 1
        line_contains_text_photo_processes_original.append(True)
    else:
        line_contains_text_photo_processes_original.append(False)

num_non_empty_lines_photo_processes_gen_rmsd_gen = 0
for i in range(0, len(photo_processes_gen_rmsd_gen_lines)):
    if(len(photo_processes_gen_rmsd_gen_lines[i].strip()) > 1):
        num_non_empty_lines_photo_processes_gen_rmsd_gen += 1
        line_contains_text_photo_processes_gen_rmsd_gen.append(True)
    else:
        line_contains_text_photo_processes_gen_rmsd_gen.append(False)

if(num_non_empty_lines_photo_processes_gen_rmsd_gen == num_non_empty_lines_photo_processes_original):  # If the files contain the same number of non-empty lines
    # Due to the fact that the boolean arrays (created above for each file) contain true and false values representing the non-emptiness or emptiness of a line.
    #     the below two lines (for each file) keep track of the index of the next-to-be-processed line in the list of lines from their respective files
    index_line_to_compare_original = 0  
    index_line_to_compare_gen_rmsd_gen = 0

    for i in range(0, num_non_empty_lines_photo_processes_original):    # It doesn't matter which list length we use, because if we made it inside this if statement, the lengths are the same
        # The below for loops ensure that the correct (non-empty) line from each photo_processes file is being compared
        while(line_contains_text_photo_processes_original[index_line_to_compare_original] == False): 
            index_line_to_compare_original += 1
        while(line_contains_text_photo_processes_gen_rmsd_gen[index_line_to_compare_gen_rmsd_gen] == False):
            index_line_to_compare_gen_rmsd_gen += 1
        
        original_line_as_list = photo_processes_original_lines[index_line_to_compare_original].split()
        gen_line_as_list = photo_processes_gen_rmsd_gen_lines[index_line_to_compare_gen_rmsd_gen].split()

        if(line_contains_text_photo_processes_original[i] and line_contains_text_photo_processes_gen_rmsd_gen[i]): 
            # We need to find the first index in the original line (as a list) where there is a fitting factor value (it is the first number besides the first element of the list, hence why i (below) is set to 1)
            index_first_float_value_in_original_line = 1  # The first entry in the line is a float, but we are looking for the first fitting factor (the 2nd float value in the line); starting at 1 therefore makes the implementation of the below while loop easier
            # The while loop below skips the first value in the list because it is an integer and would thus be counted as a float. 
            while(index_first_float_value_in_original_line < len(original_line_as_list) and not(is_float(original_line_as_list[index_first_float_value_in_original_line]))):  # Increments the index value until the first float value (fitting factor value) is found
                index_first_float_value_in_original_line += 1

            # We need to find the first index in the generalization line (as a list) where there is a fitting factor value (it is the first number besides the first element of the list, hence why i (below) is set to 1)
            index_first_float_value_in_gen_line = 1  # The first entry in the line is a float, but we are looking for the first fitting factor (the 2nd float value in the line); starting at 1 therefore makes the implementation of the below while loop easier
            # The while loop below skips the first value in the list because it is an integer and would thus be counted as a float. 
            while(index_first_float_value_in_gen_line < len(gen_line_as_list) and not(is_float(gen_line_as_list[index_first_float_value_in_gen_line]))):  # Increments the index value until the first float value (fitting factor value) is found
                index_first_float_value_in_gen_line += 1

            for j in range(0, 3): # Compare the 3 float values from each line and see if they are numerically the same (even if they are not the same as strings, which is their data type in the list)
                if(float(original_line_as_list[index_first_float_value_in_original_line + j]) != float(gen_line_as_list[index_first_float_value_in_gen_line + j])):
                    num_errors += 1
                    print("Error: the " + str(j) + " fitting factor of the line (from the photo_processes files (one each) for the original optimization and the generalization) is not the same. Below are the two lines that are not equal, where the first line is from the original optimization and the second line is from the generalization: " + photo_processes_original_lines[i] + photo_processes_gen_rmsd_gen_lines[i])
else:
    num_errors += 1
    print("Error: The photo_processes.dat files for the original optimization and the generalization have different amounts of lines.")

photo_processes_gen_rmsd_gen.close()
photo_processes_original.close()

if(num_errors == 0):
    print("Success! The generalized optimization produced the same results as the original optimization.")