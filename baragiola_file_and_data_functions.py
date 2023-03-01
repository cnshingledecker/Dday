def modelCSVFileName():
    return "./csv/bO3.csv"

# For the parallel scripts
def num_processors_to_use():
    return 4

def min_field_width():
    return 30

def serialAllOutputFileName():
    return "resultsSerialAll.txt"

def serialBestResultsOutputFileName():
    return "resultsSerialBest.txt"

def parallelAllOutputForCoreFileName():
    return "resultsParallelAll.txt"

def parallelBestResultsFileName():
    return "resultsParallelBest.txt"

def parallelTrialNuIonNuSameAllOutputForCoreFileName():
    return "resultsParallelTrialNuIonNuSameAll.txt"

def parallelTrialNuIonNuSameBestResultsFileName():
    return "resultsParallelTrialNuIonNuSameBest.txt"

def process_model_data(data):
    initialO2 = 5.7E22
    return (float(data) / initialO2) * 100