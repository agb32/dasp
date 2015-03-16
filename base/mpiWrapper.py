"""
comm requires:
send(data,rank,tag)
receive(data,rank,tag)
rank
size
barrier()
broadcast(data,rank)
abort()
"""

mpitype=None
try:
    import Scientific.MPI
    mpitype="Scientific.MPI"
except:
    pass

if mpitype==None:
    try:
        import mpi4py.MPI
        mpitype="mpi4py.MPI"
    except:
        pass

if mpitype==None:
    raise Exception("No suitable MPI modules found")

if mpitype=="Scientific.MPI":
    comm=Scientific.MPI.world.duplicate()
elif mpitype=="mpi4py.MPI":
    class mympi:
        def __init__(self):
            self.comm=mpi4py.MPI.COMM_WORLD
            self.send=self.comm.Send
            self.receive=self.comm.Recv
            self.rank=self.comm.rank
            self.size=self.comm.size
            self.barrier=self.comm.Barrier
            self.broadcast=self.comm.Bcast
            self.abort=self.comm.Abort
    comm=mympi()
