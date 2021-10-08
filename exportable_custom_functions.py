import math, os
import numpy as np

def is_int(val):
    try:
        int(val)
        return True
    except ValueError:
        return False

def is_float(val):
    try:
        float(val)
        return True
    except ValueError:
        return False

def find_nearest_index(value_to_find, pos, in_order_list):
    index = 0
    smallest_index = 0
    smallest_diff = 1e20
    for value in in_order_list:
        diff = abs(float(value_to_find) - float(value[pos]))
        if diff < smallest_diff:
            smallest_index = index
            smallest_diff = diff
        elif diff > smallest_diff:
            return smallest_index
        index += 1
    return smallest_index

# How to test this function? I'm not sure how to write an automated test to see if all the start and end indices are correctly calculated? Figure out good test cases.
def split_list_chunks(list_var, mini_chunk_size): # list_var is a list of n chunks that are each split up into mini_chunk_size evenly-sized (as much as possible) chunks
    new_split_data = []
    for i in range(0, len(list_var)): # Len(list_var) is the number of big chunks (num_chunks, if split_list was called using list_var (before calling this function))
        new_chunk = [] # For each big chunk, split_list_chunks replaces it with a list that contains the mini-chunks (which collectively contain all of the elements that were in the big chunk)
        for j in range(0, math.ceil(len(list_var[i])/mini_chunk_size)):
            start = mini_chunk_size*j
            end = mini_chunk_size*j + mini_chunk_size
            if(end >= len(list_var[i])):
                end = len(list_var[i])
            new_chunk.append(list_var[i][start:end])
        new_split_data.append(new_chunk)
    return new_split_data

# How write automated test for below function to see if all of the start and end indices are correctly calculated? Figure out good test cases.
def split_list(list_var, num_chunks): # list_var is a list of n items that is split up into num_chunks evenly-sized (as much as possible) chunks
    data_to_return = []
    num_chunks_with_extra = round((len(list_var) / num_chunks - int(len(list_var) / num_chunks)) * num_chunks)
    minimum_chunk_size = int(len(list_var) / num_chunks)
    start = 0
    end = 0
    for i in range(0, num_chunks):
        length = minimum_chunk_size
        if(i < num_chunks_with_extra):
            length += 1
        end = start + length - 1
        if(end >= len(list_var)):
            end = len(list_var) - 1
        data_to_return.append([list_var[j] for j in range(start, end+1)])
        start = end+1 # Start is updated for the next iteration of the loop
    return data_to_return

def modifyModelInpValues(linesToModify): # required format of linesToModify is specified below
    # linesToModify is a list.
    # Format of the lists that go into the above list: [lineNumber, variableValue, variableName], where
                                                       #     - lineNumber is an integer that is the line number (0-indexed) of the line that the user wants to replace the variable's value,
                                                       #     - variableValue is the value that the user wants to set the variable to (must be a number and will be converted to scientific notation with precision to 4 digits), and
                                                       #     - variableName is a string that is the name of the variable that the user wants to change the value of
                                                       #           (this should be the variable name (exact upper/lowercase) from the line where the variable is)
                                                       # List of the lines desired to be modified. Each list (1 for each line) contains (in the following order)
                                                       #     The line number (0-indexed) of the file which is desired to be modified,
                                                       #     The value that is desired to be written in for that variable, and
                                                       #     The  name of the variable (exact upper/lower casing) of the row
                                                       #         (this is done so the computer can check whether the variable that is desired to be changed is
                                                       #         at the specified line number)
                                                       # Note: do not put two line arrays within this array that have the same line number. 
                                                       #       If this is done, only the array of these with the lowest variableValue will not be used to change the line value.
    linesToModify.sort() # Done so the interpreter can not have to search through the entire list upon processing each line to see if each line number specified there 
                        #     is the current line number 
    
    lineIndex = -1 # Lines are processed as 0-indexed for this code, and this is incremented upon looking at every new line
    linesToModifyIndex = 0 # Start at the first line number and data that might need to be replaced; the index to be used to access elements of the 'linesToModify' array
    while(linesToModify[linesToModifyIndex][0] < 0):
        linesToModifyIndex += 1 # Skip over negative line numbers (only should happen due to user error)

    outputModelFile = open("modelTemp.inp", "w") # A temp file whose contents will be written into model.inp upon completion of the below 'with' block

    with open("model.inp", "r", newline='') as inpFile:
        for line in inpFile:
            lineIndex += 1
            lineToWrite = line # The value of this variable will be written (may be modified if the user wants to change the value for that line)
            if(linesToModifyIndex < len(linesToModify)): # If there are still lines that may need to be checked
                if(len(line) >= 35):    # len(line) >= 35 is for: If the current line is one with a variable name and a value (based on how model.inp is set up)
                    if(linesToModify[linesToModifyIndex][0] == lineIndex):
                        lineText = line[0:(len(linesToModify[linesToModifyIndex][2].strip()))].strip() # The spots where the variable name would be in the line 
                        userVarName = linesToModify[linesToModifyIndex][2].strip() # What the user said the variable name was

                        if(lineText == userVarName): # If the line is the one with the variable name specified in the second spot in the array and the line number is the one specified during the array
                            lineToWrite = line[0:24] + np.format_float_scientific(float(linesToModify[linesToModifyIndex][1]), precision = 4, unique=False) + line[34:len(line)] # Replace the variable value with the user-specified value
                        linesToModifyIndex += 1
                        while(linesToModifyIndex < len(linesToModify) - 1 and linesToModify[linesToModifyIndex][0] == linesToModify[linesToModifyIndex - 1][0]): # While there are more line arrays with the same line index as the one just processed
                            linesToModifyIndex += 1
            outputModelFile.write(lineToWrite)
        inpFile.close()
    outputModelFile.close()

    os.system("cat modelTemp.inp > model.inp")
    os.system("rm modelTemp.inp")