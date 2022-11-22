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

bestRMSD = 1e80

bestRMSDLines = ["", "", "", "", "", "", ""]
for i in range(0, 4):
    with open(f"core{i}results", "r") as resultsFile:
        fileLines = resultsFile.readlines()
        fileLinesSplit = split_list(fileLines, int((len(fileLines) + 1) / 8))
        for resultSet in fileLinesSplit:
            rmsd = float(resultSet[6][0:resultSet[6].find(" ")])
            if(rmsd < bestRMSD):
                bestRMSD = rmsd
                for i in range(0, 7):
                    bestRMSDLines[i] = resultSet[i]

with open("serverResultsParallel", "w") as bestResultsFile:
    bestResultsFile.writelines(bestRMSDLines)
