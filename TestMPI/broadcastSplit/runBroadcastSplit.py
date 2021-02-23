# This file tests: the runtime of parallel and non-parallel scripts.
#                  Note: do not run this file, as it will not work as intended with the current versions of broadcastSplit.py and nonParallelBroadcastSplit.py

import os, time

ts = time.time() * 1000
os.system("mpiexec -n 4 python3 broadcastSplit.py")
te = time.time() * 1000
print("Runtime of parallel version was " + str(te-ts) + " milliseconds.")

ts = time.time() * 1000
os.system("python3 nonParallelBroadcastSplit.py")
te = time.time() * 1000
print("Runtime of non-parallel version was " + str(te-ts) + " milliseconds.")