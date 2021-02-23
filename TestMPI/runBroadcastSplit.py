import os, time

ts = time.time() * 1000
os.system("mpiexec -n 4 python3 broadcastSplit.py")
te = time.time() * 1000
print("Runtime of parallel version was " + str(te-ts) + " milliseconds.")

ts = time.time() * 1000
os.system("python3 nonParallelBroadcastSplit.py")
te = time.time() * 1000
print("Runtime of non-parallel version was " + str(te-ts) + " milliseconds.")