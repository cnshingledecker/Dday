from mpi4py import MPI
import numpy as np
import itertools

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    a = list(np.linspace(1.7,2.7,10))
    b = list(np.linspace(1,2,1))
    c = list(np.linspace(0.3,0.35,1))
    data = [a,b,c]
    data = itertools.product(*data)
elif rank != 0:
    data = None
    
data = comm.bcast(data, root=0)

if rank >= 0:
    print("This came to node " + str(rank) + ": ")
    for j in data:
        print(str(rank) + "             " + str(j))
