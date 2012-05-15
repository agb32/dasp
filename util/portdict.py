#!/usr/local/bin/python
"""A simple daemon which can record all open simulation ports, users, files and ranks."""
import socket,select,threading,sys
sys.path.insert(1,"/var/opt/pce/home/ali/cvsstuff/aosim")
import util.serialise as serialise
class portdict:
    def __init__(self,hostList=["129.234.187.10","192.168.3.30","10.128.27.86"],port=8999,startThread=1,debug=0):
        self.lsockList=[]
        self.selIn=[]
        self.debug=debug
        for host in hostList:
            lsock=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            lsock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR,1)
            lsock.bind((host,port))
            lsock.listen(1)
            self.lsockList.append(lsock)
            self.selIn.append(lsock)
        self.selOut=[]
        self.go=1
        self.simList={}
        if startThread:
            self.t=threading.Thread(target=self.loop,args=())
            self.t.setDaemon(0)
            self.t.start()
    def loop(self):
        if self.debug:
            print "portdict entering main loop"
        while self.go:
            selOut=[]
            rtr,rtw,err=select.select(self.selIn,selOut,self.selIn)
            for s in err:
                #remove from dict, and act on it...
                self.close(s)
            for s in rtr:#ready to read the socket
                if s in self.lsockList:#listening socket - new client...
                    if self.debug:
                        print 'Accepting new client...'
                    try:
                        conn,raddr=s.accept()
                        if self.debug:
                            print "from",raddr
                        self.selIn.append(conn)
                    except:
                        if self.debug:
                            print "Failed to accept client"
                else:#read socket
                    if self.readsock(s)==-1:#connection closed - sim ended
                        self.close(s)
    def close(self,sock):
        """Close socket and stop listening on it
        @param sock: Socket
        @type sock: socket.socket instance
        """
        sock.close()
        if self.debug:
            print "Closed",sock
        if sock in self.simList.keys():
            del(self.simList[sock])
        if sock in self.selIn:
            self.selIn.remove(sock)
        #print "Closed ",sock
    def readsock(self,sock):
        data=serialise.ReadMessage(sock.fileno())
        if not data:
            return -1
        if self.debug:
            print "Got data:",data
        if data[0]=="add":#add this socket to the list
            if self.debug:
                print "Adding data to simList"
            host=data[1]
            port=data[2]
            rank=data[3]
            user=data[4]
            simID=data[5]
            file=data[6]
            if len(data)>7:
                batchno=data[7]
            else:
                batchno=""
            
            #print type(port),type(rank)
            self.simList[sock]=(host,port,rank,user,batchno,simID,file)
            if self.debug:
                print self.simList
        elif data[0]=="get":#get the list of sockets and return it...
            data=self.simList.values()
            if self.debug:
                print "Sending",data,"simList is",self.simList
            serialise.Send(data,sock)
        return 0
            
            
if __name__=="__main__":
    try:
        p=portdict(startThread=0)
    except: # port in use
        import util.analyse
        sim=util.analyse.analyse()
        data=sim.queryPortDict()
        print "Connected simulations:"
        print data
    else: # run the main loop
        p.loop()
