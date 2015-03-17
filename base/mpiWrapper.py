"""
comm requires:
send(data,rank,tag)
receive(data,rank,tag) -> data,rank,tag,nelements
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
            #self.receive=self.comm.Recv
            self.rank=self.comm.rank
            self.size=self.comm.size
            self.barrier=self.comm.Barrier
            self.broadcast=self.comm.Bcast
            self.abort=self.comm.Abort
        def receive(self,data,rank,tag):
            #print "receive data type: %s %s (rank %d, tag %d), rms %g"%(data.dtype.char,str(data.shape),rank,tag,data.std())
            self.comm.Recv(data,rank,tag)
            #print "Received std: %g"%data.std()

            #if data==None:
            #    rt=0
            #else:
            #    rt=1
            rt=1
            return data,rank,tag,rt

        #def send(self,data,rank,tag):
        #    #print "send data type: %s %s (rank %d, tag %d)"%(data.dtype.char,str(data.shape),rank,tag)
        #    self.comm.Send(data,rank,tag)
    comm=mympi()
