import time
import os

# startTime1 = int(round(time.time() * 1000))
# os.system("python3 baragiola_optimization.py")
# endTime1 = int(round(time.time() * 1000))

startTime2 = int(round(time.time() * 1000))
os.system("python3 baragiola_optimization_generalization_rmsd_generalization.py")
endTime2 = int(round(time.time() * 1000))

startTime3 = int(round(time.time() * 1000))
os.system("python3 baragiola_optimization_test.py")
endTime3 = int(round(time.time() * 1000))


# print("Original (baragiola_optimization.py) took " + str((endTime1 - startTime1) / 1000) + " seconds")
print("Generalization (baragiola_optimization_generalization_rmsd_generalization.py) took " + str((endTime2 - startTime2) / 1000) + " seconds")
print("Modified original (baragiola_optimization_test.py) took " + str((endTime3 - startTime3) / 1000) + " seconds")
