#$Id: SockConnTest.py,v 1.4 2008/03/05 15:43:09 ali Exp $
import socket,serialise,os,numpy,time
class SockConnTest:
    """Class for testing purposes"""
    def __init__(self,port,host=os.environ["HOSTNAME"]):
        self.port=port
        self.host=host
        self.sock=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.connect((host,port))
    def test(self):
        """Send test data"""
        a=numpy.array([[1,22,33],[4,55,66]],"i")
        b=numpy.array([[11,122,133],[14,155,166]],"i")
        b[1]=a[0]
        print "Contiguous?",b.iscontiguous(),b
        serialise.Send(["hello there",10,(1,2,3),{5:2,3:4},a,b],self.sock)
        print "sent..."

if __name__=="__main__":
    """Testing purposes"""
    import sys
    if len(sys.argv)==1:
        s=SockConnTest(9000)
        s.test()
    else:
        s=SockConnTest(8000)
        s.test()
        data=serialise.ReadMessage(s.sock.fileno())
        print data
