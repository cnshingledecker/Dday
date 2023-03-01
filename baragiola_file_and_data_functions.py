def modelCSVFileName():
    return "./csv/bO3.csv"

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