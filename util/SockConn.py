#dasp, the Durham Adaptive optics Simulation Platform.
#Copyright (C) 2004-2016 Alastair Basden and Durham University.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as
#published by the Free Software Foundation, either version 3 of the
#License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

import socket,serialise,types,os,select,cPickle,thread,time
import traceback
import util.ConnObj as ConnObj
#import Scientific.MPI
import sys
class SockConn:
    
    """Class to open a listening socket, and start a select loop in a
    new thread.  Each simulation process will start a sockConn object.
    It opens a listening socket (if it cannot bind to the given port,
    it will simply try increasing port numbers until it succeeds), and
    enters a main loop (in a separate thread from the main process)
    waiting for clients to connect.  Once connected, clients (such as
    ConnObj objects) send instances of ConnMsg objects (actually
    ConnMsg.pickle() objects), which are used to instruct the main
    simulation process.
    
    A received command can be executed immediately (not advisable
    because the simulation will be at an unknown state), or at the end
    of an iteration cycle.  Additionally, the same command can be
    executed at the end of every iteration cycle if desired.  These
    commands have access to the entire simulation process, and so can
    be used to alter state, change parameters, or (more normally)
    return data to the client process (user).     

    Class variables (important to simulation programmer):
     - port - int, port number to listen on
     - host - string, hostname
     - globals - global dictionary
    Class variables (not important to simulation programmer):
     - go - int, determines whether to listen or not
     - fwd - tuple or None, whether a forwarding socket
     - lsock - socket object, listening socket
     - printmsg - int, whether to print messages
     - selIn - list, for select()
     - selOut - list, for select()
     - cmdList - list, of commands to be executed
     - rptCmdList - list, of commands to be executed every iteration
     - fsock - forwading socket

    @cvar port: port number to listen on
    @type port: int
    @cvar host: hostname
    @type host: String
    @cvar go: determines whether to listen or not
    @type go: Int
    @cvar fwd: whether a forwarding socket
    @type fwd: Tuple (hostname,port) or None
    @cvar lsock: listening socket
    @type lsock: socket.socket instance
    @cvar printmsg: whether to print messages
    @type printmsg: Int
    @cvar selIn: for select()
    @type selIn: List
    @cvar selOut: for select()
    @type selOut: List
    @cvar cmdList: commands to be executed
    @type cmdList: List
    @cvar rptCmdList: commands to be executed every iteration
    @type rptCmdList: List
    @cvar fsock: forwading socket
    @type fsock: socket.socket instance
    @cvar globals: global dictionary
    @type globals: Dict
    """
    def __init__(self, port, host="", fwd=None,globals=None,startThread=1,listenSTDIN=1,mpiWorldSize=1):
        """Opens a listening port, and either acts on commands send, or if
        fwd is (host,port), forwards commands to this port (little used)
        @param port: Port number to listen on
        @type port: Int
        @param host: hostname
        @type host: String
        @param fwd: Forwarding socket information
        @type fwd: None or Tuple
        @param globals: Global dictionary
        @type globals: Dict
        @param startThread: Whether to start the listening loop
        @type startThread: Int
        """
        self.port=port
        self.host=host
        self.go=1
        self.fwd=fwd
        self.printmsg=0
        if port is None:
            self.lsock=None
            self.selIn=[]
        else:
            self.lsock=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.lsock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR,1)
            bound=0
            cnt=0
            while bound==0 and cnt<1000:
                try:
                    self.lsock.bind((host,port))
                except:
                    print "WARNING Couldn't bind to port %d.  "%port,
                    port+=mpiWorldSize#Scientific.MPI.world.size#inc by number of processes in this mpi run
                    print "INFORMATION Trying port %d"%port
                    cnt+=1
                    self.port=port
                else:
                    bound=1
                    print "INFORMATION Bound to port %d"%port
            self.lsock.listen(1)
            self.selIn=[self.lsock]
        self.listenSTDIN=listenSTDIN
        if mpiWorldSize>1:#Scientific.MPI.world.size>1:
            self.listenSTDIN=0
        if self.listenSTDIN:
            self.selIn.append(sys.stdin)
        self.selOut=[]
        self.cmdList=[]
        self.rptCmdList=[]
        self.fsock=None
        self.globals=globals
        if self.fwd!=None:
            self.fsock=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.fsock.connect(self.fwd)
            self.selIn.append(self.fsock)
        if startThread:
            time.sleep(1)
            thread.start_new_thread(self.loop,())
    def setGlobals(self,globs):
        """sets the globals which are the objects we have access to...
        @param globs: typically globals()
        @type globs: Dict"""
        self.globals=globs
    def endLoop(self):
        """End the loop"""
        self.go=0
        
    def loop(self):
        """probably start this in a new thread... listens on the socket, and
        does stuff depending on what it gets"""
        while self.go:
            if self.fsock==-1:
                break
            selOut=[]
            try:
                rtr,rtw,err=select.select(self.selIn,selOut,self.selIn)
            except:
                traceback.print_exc()
                print "WARNING in select - continuing..."
                time.sleep(1)
                rtr=[]
                rtw=[]
                err=[]
            for s in err:
                #remove from dict, and act on it...
                self.close(s)
            for s in rtr:#ready to read the socket
                if s==self.lsock:#listening socket - new client...
                    print 'INFORMATION Accepting new client...'
                    conn,raddr=s.accept()
                    print "from %s"%str(raddr)
                    self.selIn.append(conn)
                elif s==sys.stdin:
                    #print "Got from stdin"
                    try:
                        data=s.readline()[:-1]
                    except:
                        data="Error reading from stdin"
                    #print self.globals
                    try:
                        ft=self.globals["ctrl"].frametime
                        ti=self.globals["ctrl"].thisiter
                    except:
                        ti=0
                        ft=1.
                    if ft==0:
                        ft=1.
                    print "INFORMATION At iteration %d, taking %gs per frame (%g fps) on port %s"%(ti,ft,1./ft,str(self.port))
                    if len(data)>0:
                        if data[0]=="s":
                            try:
                                txt=""
                                ctrl=self.globals["ctrl"]
                                scienceList=ctrl.compList
                                dlist=[]
                                for s in scienceList:
                                    if hasattr(s,'strParams') and s not in dlist:
                                        dlist.append(s)
                                        for i in range(len(s.thisObjList)):
                                            txt+=s.thisObjList[i].idstr+': '+s.strParams(i)+'\n'
                                txt+='Iter %d frametime %g (mean %g) batch %d\n%s'%(ctrl.thisiter,ctrl.frametime,ctrl.meanTiming.sum()/ctrl.thisiter,ctrl.batchno,ctrl.simID)
                                print txt
                            except:
                                traceback.print_exc()
                        elif data[0]=="h":
                            print "HELP: <ret> for iteration number\ns<ret> for science information"
                        else:
                            print data,len(data)
                        
                elif s==self.fsock:#forward data... (return to client)
                    data=s.recv(1024)
                    length=len(data)
                    if length==0:
                        self.close(s)
                    else:
                        for sk in self.selIn:
                            if sk!=self.lsock and sk!=self.fsock:
                                try:
                                    sk.sendall(data)
                                except:
                                    self.close(sk)
                else:#readsock or forward data...
                    if self.fsock!=None:#forward the data...
                        data=s.recv(1024)
                        if len(data)==0:
                            self.close(s)
                        else:
                            try:
                                self.fsock.sendall(data)
                            except:
                                serialise.Send(["warning",None,"couldn't forward data"],s)
                    else:#readsock.
                        #print 'Reading socket...'
                        if self.readsock(s)==-1:#connection closed
                            try:
                                print 'INFORMATION Closing socket %s'%str(s.getpeername())
                            except:
                                print "INFORMATION Closing socket"
                            self.close(s)
            
    def close(self,sock):
        """Close socket and stop listening on it
        @param sock: Socket
        @type sock: socket.socket instance
        """
        sock.close()
        self.selIn.remove(sock)
        #print "Closed ",sock

    def readsock(self,sock):
        """Reads the socket, obtaining 'data'.
        data.action will be a control word, e.g. 'cmd'
        If data.action=='now':
        data.command will be command to be exec'd.
        data.ret (if present) will be list of things to return.
        If data.action=='rpt':#repeat command
        data.command will be command to be exec'd.
        data.ret (if present) will be list of things to return.
        If data.action=='cmd':#command to execute during break...
        data.command will be command to be exec'd.
        data.ret (if present) will be list of things to return.
        If data.action=='del':#delete command from regular/cmd list
        this will be deleted...
        @param sock: Socket to read
        @type sock: socket.socket instance
        """
        try:
            tmp=serialise.ReadMessage(sock.fileno())
        except:
            print "WARNING Connection reset by peer."
            return -1
        if not tmp:
            return -1 #connection closed.
        #print "socketdata:",tmp
        data=ConnObj.ConnMsg(None,None)
        data.unpickle(tmp)
        #data=cPickle.loads(data[0])#unpickle a ConnMsg object.
        action=data.action
        tag=data.tag
        print "INFORMATION got data: %s %s %s %s"%(str(data.action),str(data.command),str(data.tag),str(data.ret))
        #cmd=data.pop(0)
        if action=="now":
            if self.globals==None:
                self.cmdList.append([data,sock])
            else:
                self.execCmd(data,sock)
        elif action=="cmd":
            self.cmdList.append([data,sock])
        elif action[:3]=="rpt":
            freq=1
            if len(action)>3:
                freq=int(action[3:])
                if freq<1:
                    freq=1
            self.rptCmdList.append([data,sock,freq,0])
        elif action=="del":
            remlist=[]
            #print "Action del:",data,sock
            for cmd in self.rptCmdList:
                #print cmd
                try:
                    if (data.command=="" or data.command==None or cmd[0].command==data.command) and sock==cmd[1] and (tag==None or cmd[0].tag==tag):
                        #if cmd[:2]==[data,sock]:
                        remlist.append(cmd)
                except:
                    print "ERROR deleting",data,sock,cmd
                    print data.command
            for cmd in remlist:
                print "INFORMATION Deleting action",cmd
                self.rptCmdList.remove(cmd)
            while [data,sock] in self.cmdList:
                print "Deleting action:",[data,sock]
                self.cmdList.remove([data,sock])
        elif action=="add":
            #prepend data.command to config.postList...
            print "INFORMATION SockConn - got data action add"
            if type(data.command)==type(()) and len(data.command)==2 and type(data.command[0])==type(""):
                if type(self.globals)!=type(None):
                    print "INFORMATION Adding to config.postList - %s"%data.command[0]
                    self.globals["ctrl"].config.postAdd(data.command)
                    print "Added post variable %s to config.postList"%data.command[0]
                else:
                    print "ERROR Cannot add post variable %s to config, SockConn has no globals"%data.command[0]
            else:
                print "ERROR SockConn - action add not received with non-valid data, should be tuple of (str,data)."
        else:
            print action,data
            serialise.Send(["warning",tag,"data not understood"],sock)

    def execCmd(self,data,sock):
        """Execute a command.
        @param data: The command to be executed
        @type data: ConnMsg instance
        @param sock: A socket
        @type sock: socket.socket instance
        @return:  The executed value
        @rtype: User defined
        """
        #command=data[0]
        command=data.command
        tag=data.tag
        rtval=0
##         if len(data)>1:
##             rtObjects=data[1]
##             if type(rtObjects)!=types.ListType:
##                 rtObjects=[rtObjects]
##         else:
##             rtObjects=None
        rtObjects=data.ret
        d={}#newly created variables go in here...
        if self.printmsg:
            print "INFORMATION Executing",command
        try:
            exec command in self.globals,d
        except Exception,msg:
            data.retryCnt+=1
            txt="Will retry"
            if data.retryCnt>data.maxRetry:
                rtval=1
                txt="Cancelling"
            print "ERROR Command Exec failed1:",command,msg
            #print str(sys.exc_info()),str(sys.exc_info()[1].args)
            #print self.globals
            if sock!=None:
                try:
                    serialise.Send(["error",tag,"Command execution failed (%s): %s %s"%(txt,command,msg)],sock)
                except:
                    print "ERROR Serialise failed in SockConn.execCmd - couldn't send error message"
        except:
            data.retryCnt+=1
            txt="Will retry"
            if data.retryCnt>data.maxRetry:
                rtval=1
                txt="Cancelling"
            print "ERROR Command exec failed2:",command
            if sock!=None:
                try:
                    serialise.Send(["error",tag,"Command execution failed (%s): %s"%(txt,command)],sock)
                except:
                    print "ERROR Serialise failed in SockConn.execCmd - couldn't send error message"
        else:
            rt={}#will return rt to the user.
            l=0
            if rtObjects==None:#return all objects created...
                rtObjects=d.keys()
            elif type(rtObjects) not in [types.ListType,types.TupleType]:
                rtObjects=[rtObjects]
            for key in rtObjects:
                if d.has_key(key):
                    rt[key]=d[key]
                    l+=1
                elif key not in [""," "]:
                    print "ERROR key",key,"not found"
            #sock.send(serialise.Serialise(["data",rt]))
            if l>0:#don't reply if nothing to reply with!
                if sock!=None:
                    if self.printmsg:
                        print "INFORMATION Sending data over socket"
                    try:
                        serialise.Send(["data",tag,rt],sock)
                    except:
                        rtval=1
                        print "ERROR Serialise failed in SockConn.execCmd for tag %s with keys %s (error: %s %s)"%(str(tag),str(rt.keys()),str(sys.exc_info()),str(sys.exc_info()[1].args))

        return rtval
    def doCmdList(self,thisiter):
        """This will be called at the end of every next() iteration, and executes everything waiting to be executed.
        @param thisiter:  Iteration number
        @type thisiter: Int
        """
        rlist=[]
        if self.globals!=None:
            for data,sock in self.cmdList:#only execute this once...
                self.execCmd(data,sock)
                rlist.append([data,sock])
            for r in rlist:#...so now remove from the list.
                try:
                    self.cmdList.remove(r)
                except:
                    pass
            rlist=[]
            for data in self.rptCmdList:#repeat these continually...
                freq=data[2]
                nextiter=data[3]
                if nextiter<=thisiter:
                    if self.execCmd(data[0],data[1])==1:#failed
                        rlist.append(data)
                    if freq>=0:
                        data[3]=thisiter+freq
                    else:#remove it from the list.
                        rlist.append(data)
            for r in rlist:
                try:
                    self.rptCmdList.remove(r)
                except:
                    pass
    def printList(self):
        """Prints the command list"""
        print "INFORMATION rpt command list:"
        for r in self.rptCmdList:
            print r
if __name__=="__main__":
    """Testing function"""
    import sys
    if len(sys.argv)==1:
        s=SockConn(9000)
        s.loop()
    elif sys.argv[1]=="f":
        s=SockConn(8000,fwd=(socket.gethostname(),9000))
        s.loop()

