from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    data = {'a': 7, 'b': 3.14}
    comm.send(data, dest=1)
    data2 = {'a': 57, 'b': 19}
    comm.send(data2, dest=2)
    data3 = {'a': "Nick", 'b': "Alarm"}
    comm.send(data3, dest=3)
elif rank > 0:
    data = comm.recv(source=0)
    print(str(data) + " came to node of rank " + str(rank))