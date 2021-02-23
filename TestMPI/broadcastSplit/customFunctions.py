import math

# How to test this function? I'm not sure how to write an automated test to see if all the start and end indices are correctly calculated? Figure out good test cases.
def split_array_chunks(array, mini_chunk_size): # Array is a list of n chunks that are each split up into mini_chunk_size evenly-sized (as much as possible) chunks
    new_split_data = []
    for i in range(0, len(array)): # Len(array) is the number of big chunks (num_chunks, if split_array was called using array (before calling this function))
        new_chunk = [] # For each big chunk, split_array_chunks replaces it with a list that contains the mini-chunks (which collectively contain all of the elements that were in the big chunk)
        for j in range(0, math.ceil(len(array[i])/mini_chunk_size)):
            start = mini_chunk_size*j
            end = mini_chunk_size*j + mini_chunk_size
            if(end >= len(array[i])):
                end = len(array[i])
            new_chunk.append(array[i][start:end])
        new_split_data.append(new_chunk)
    return new_split_data

# How write automated test for below function to see if all of the start and end indices are correctly calculated? Figure out good test cases.
def split_array(array, num_chunks): # Array is a list of n items that is split up into num_chunks evenly-sized (as much as possible) chunks
    data_to_return = []
    num_chunks_with_extra = round((len(array) / num_chunks - int(len(array) / num_chunks)) * num_chunks)
    minimum_chunk_size = int(len(array) / num_chunks)
    start = 0
    end = 0
    for i in range(0, num_chunks):
        length = minimum_chunk_size
        if(i < num_chunks_with_extra):
            length += 1
        end = start + length - 1
        if(end >= len(array)):
            end = len(array) - 1
        data_to_return.append([array[j] for j in range(start, end+1)])
        start = end+1 # Start is updated for the next iteration of the loop
    return data_to_return