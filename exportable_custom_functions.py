import math

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