import math

# How to test this function? I'm not sure how to write an automated test to see if all the start and end indices are correctly calculated? Figure out good test cases.
def split_array_chunks(array, mini_chunk_size):
    new_split_data = []
    for i in range(0, len(array)):
        # print("Chunk size is " + str(len(array[i])))
        new_chunk = []
        for j in range(0, math.ceil(len(array[i])/mini_chunk_size)):
            start = mini_chunk_size*j
            end = mini_chunk_size*j + mini_chunk_size
            if(end >= len(array[i])):
                end = len(array[i])
            # print("Start: " + str(start) + " "*(5-len(str(start))) + " End: " + str(end-1))
            new_chunk.append(array[i][start:end])
        new_split_data.append(new_chunk)
    return new_split_data

def split_array_test(array, num_chunks):
    testing = True
    if(testing is True):
        print("\n\n\nSplitting an array of " + str(len(array)) + " elements into " + str(num_chunks) + " chunks.")
    data_to_return = []
    num_chunks_with_extra = round((len(array) / num_chunks - int(len(array) / num_chunks)) * num_chunks)
    minimum_chunk_size = int(len(array) / num_chunks)
    if(testing is True):
        pass
        # print("Num chunks with extra is " + str(num_chunks_with_extra))
        # print("Min chunk size is " + str(minimum_chunk_size))
    total = 0
    start = 0
    end = 0
    for i in range(0, num_chunks):
        length = minimum_chunk_size
        if(i < num_chunks_with_extra):
            length += 1
        end = start + length - 1
        if(end >= len(array)):
            end = len(array) - 1
        if(testing is True):
            print("    Start: " + str(start) + " "*(6-len(str(start))) + " End: " + str(end) + " "*(5-len(str(end))) + "Size: " + str(end-start+1))
        total += (end-start+1)
        data_to_return.append([array[j] for j in range(start, end+1)])
        start = end+1
    # if(testing == True):
    assert total==len(array)
    return data_to_return

# How write automated test for below function to see if all of the start and end indices are correctly calculated? Figure out good test cases.
def split_array(array, num_chunks):
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
        start = end+1
    return data_to_return