#!/usr/bin/env python
"""Python code for analysing a simulation - either from the command line or
from a gui.  This allows the simulation state to be analysed and changed"""
import types,socket,select
import util.serialise as serialise
## class Connection:
##     def __init__(self,host,port):
##         self.host=host
##         self.port=port
##         self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
##         try:
##             self.sock.connect((host,port))
##         except:
##             print "Couldn't connect"
##             self.sock=None
##     def close(self):
##         if self.sock!=None:
##             self.sock.close()
##         self.sock=None

class analyse:
    """Class for analysing (and altering) simulation state

    Class variables (important to simulation programmer):
     - connList - list, of current socket connections.
     - savedTag - tag to use with command if one isn't supplied by user
     - recDataList - data received from simulation but not yet dealt with
     - dataProcessDict - dict, default functions to apply when data from a given socket with a given tag is returned
    @cvar connList: Current open socket connections
    @type connList: List
    @cvar savedTag: Default tag to use with data if none supplied by user
    @type savedTag: Int
    @cvar recDataList: Data received from simulation, but not yet dealt with
    @type recDataList: List
    @cvar dataProcessDict: Default responce functions for data from simulation
    @type dataProcessDict: Dict
    """
    def __init__(self):
        """Initialise simulation analyser"""
        self.connList=[]
        self.savedTag=2**30#note tag can be anything - int, string etc.
        self.recDataList=[]
        self.dataProcessDict={}#dictionary with key=(socket,tag), value=method
    def openConnection(self,host,port):
        """Open a connection to simulation.
        @param host: Hostname
        @type host: String
        @param port: Port to connect to
        @type port: Int
        @return: Opened connection
        @rtype: socket.socket instance
        """
        #conn=Connection(host,port)
        #can use conn.getpeername() to get the IP/port.
        conn=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            print "Connecting..."
            conn.connect((host,port))
            print "Connected..."
            self.connList.append(conn)
        except:
            print "Couldn't connect"
            conn=None
        return conn
    def closeConnection(self,conn=None):
        """Close a connection.
        @param conn: The connection to close (or list of) (assume all, if conn is None)
        @type conn: None, List or socket.socket instance
        """
        if conn==None:
            conn=self.connList
        if type(conn)!=types.ListType:
            conn=[conn]
        closedList=[]
        for c in conn:
            c.close()
            closedList.append(c)
        for c in closedList:
            if c in self.connList:
                self.connList.remove(c)
    def parseConnList(self,connList):#get it into a list form.
        """Put None, a socket.socket instance or a list into list form.
        @param connList: The connection(s)
        @type connList: None, socket.socket instance or List
        @return: List of socket.socket instances
        @rtype: List
        """
        if connList==None:
            connList=self.connList
        if type(connList)!=types.ListType and type(connList)!=types.TupleType:
            connList=[connList]
        return connList
    def queryPortDict(self,host="129.234.187.10",port=8999):
        conn=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        conn.connect((host,port))
        serialise.Send(["get"],conn)
        data=serialise.ReadMessage(conn.fileno())
        conn.close()
        return data
        
    def getTag(self):
        """Get a new tag.
        @return: The new tag
        @rtype: String
        """
        self.savedTag+=1
        if self.savedTag==2**31:
            self.savedTag=2**30
        return str(self.savedTag)
    def stop(self,connList=None):
        """Stop the simulation
        @param connList:  List of sockets to stop
        @type connList: None, List or socket.socket instance
        """
        #print "Stopping simulation"
        self.execute("ctrl.go=0",connList=connList)
        #self.execute("ctrl.end()",action="now",connList=connList)
    def kill(self,connList=None):
        """Stop the simulation
        @param connList:  List of sockets to stop
        @type connList: None, List or socket.socket instance
        """
        #print "Stopping simulation"
        self.execute("import os;os.system('killall mpipython')",connList=connList)
        #self.execute("ctrl.end()",action="now",connList=connList)
    def pause(self,connList=None):
        """Pause the simulation
        @param connList:  List of sockets to stop
        @type connList: None, List or socket.socket instance
        """
        #print "Sending pause command"
        self.execute("ctrl.paused=1",connList=connList)
    def run(self,niters=None,connList=None):
        """Run (unpause) the simulation
        @param niters: Number of iterations to run for (or None to run indefinately).
        @type niters: None or Int
        @param connList:  List of sockets to stop
        @type connList: None, List or socket.socket instance
        """
        #print "Unpausing simulation"
        if niters!=None:
            self.execute("ctrl.addnextniters(%d)"%niters,connList=connList)
        self.execute("ctrl.paused=0",connList=connList)
    def getObjects(self,regexp=None,base=None,printit=0,ignoreList=["config","parent"],connList=None,readsock=1,indent="    "):
        """
        Use this method to query a running simulation about existing objects.
        Base can be None (everything) or a string which should be an object
        within the simulation, e.g. 'Cent' or 'wfs.wfs_int' or similar.
        ignoreList is a list of strings to be ignored... e.g. ["config"] would
        ignore cent.config etc.
        If printit is 1, the simulation will print out what it is returning
        (useful for debugging process).
        return a list of tuples of (socket,objectstring)
        @param base: Object to start with (or None for whole simulation)
        @type base: None or String
        @param printit: Whether to print results in simulation terminal
        @type printit: Int
        @param ignoreList: Object to ignore
        @type ignoreList: Int
        @param connList: Connection list
        @type connList: None, List or socket.socket instance
        @param readsock: Whether to read the socket for possible response (blocking)
        @type readsock: Int
        @param indent:  Indentation to use
        @type indent: String
        @return: The data obtained from simulation
        @rtype: List
        """
        connList=self.parseConnList(connList)
        ignoreList=str(ignoreList)
        if base==None:
            base="None"
        else:
            base="'"+base+"'"
        if regexp==None:
            regexp="None"
        else:
            regexp="'"+regexp+"'"
        cmdtxt="txt=ctrl.getTree(regexp=%s,base=%s,printit=%d,ignoreList=%s,indent='%s')"%(regexp,base,printit,ignoreList,indent)
        self.execute(cmdtxt,tag="getObj",rt="txt",connList=connList)
        dataList=[]
        if readsock:
            for conn in connList:#sockets with data ready to read...
                data=self.read(connList=conn,tag="getObj",block=1)
                dataList.append((data[0][0],data[0][3]["txt"]))
        return dataList
        
    def execute(self,command,tag="NOTAG",rt=None,action="cmd",connList=None):
        """Tell the simulation to execute a command
        tag can be anything (string).  Gets sent back via serialise.
        @param command: The command to be exec'd
        @type command: String
        @param tag: data tag
        @type tag: String
        @param rt: Data to return
        @type rt: None or String
        @param action: The action to perform
        @type action: String
        @param connList: socket list
        @type connList: None, List or socket.socket instance
        @return: The tag
        @rtype: String
        """
        if action not in ["cmd","del","now"]:
            if action[0:3]!="rpt":
                print "Action %s not recognised"%action
                if action=="add":
                    print "Note - action add not allowed in util.analyse..."
                return
        connList=self.parseConnList(connList)
        if tag=="NOTAG":
            if action=="del":
                tag=None
            else:
                tag=self.getTag()
        lst=[action,command,rt,tag]
        remlist=[]
        for conn in connList:
            try:
                serialise.Send(lst,conn)
            except:
                remlist.append(conn)
        for c in remlist:
            connList.remove(c)
        return tag
    def process(self,connList=None,readSocks=1):
        """Process any input data which has an entry in the
        self.dataProcessDict
        @param connList: socket list
        @type connList: None, List or socket.socket instance
        @param readSocks: Whether to read sockets or not
        @type readSocks: Int
        @return: List of data handled
        @rtype: List
        """
        connList=self.parseConnList(connList)
        if readSocks:
            self.readSockets(connList)
        remList=[]
        for data in self.recDataList:
            key=(data[0],data[2])
            #print "key:",key
            if self.dataProcessDict.has_key(key):
                #socket and tag match...
                if self.dataProcessDict[key](data)==1:#failed command
                    pass
                remList.append(data)
        for d in remList:
            self.recDataList.remove(d)
        return remList#return list of things handled.
    def addCallback(self,conn,tag,method):
        """Add a callback to the data process dictionary
        @param conn: socket
        @type conn: socket.socket instance
        @param tag: The tag to act on
        @type tag: String
        @param method: The function to call
        @type method: Function
        """
        self.dataProcessDict[(conn,tag)]=method
    def getData(self):
        """Get a list of any data returned from all sockets of the simulation.
        This is provided for convenience, rather than use, as only the data
        is returned, so you will not be able to tell which socket this has come
        from, and what the associated tag is, or what the data was.  It is
        envisaged for use in the simple case where only one connection is
        open, and where you are only expecting one piece of data (or at least
        know what the datas will be).
        @return: A list of requested data
        @rtype: List
        """
        retn=[]
        dl=self.read()
        for d in dl:
            for k in d[3].keys():
                retn.append(d[3][k])
        return retn
    
    def read(self,connList=None,tag=None,maxlen=0,block=0):
        """get a list of the datas that match given socket and tag.
        @param connList: socket list
        @type connList: None, List or socket.socket instance
        @param tag: The tag to look for
        @type tag: String
        @param maxlen: Maximum number of data items to return
        @type maxlen: Int
        @param block: Whether to block while reading socket
        @type block: Int
        @return: The requested data
        @rtype: List
        """
        connList=self.parseConnList(connList)
        dataList=self.readSpecific(connList,tag,["data"],maxlen)
        if len(dataList)==0 and block!=0:#no data found - try to read for more, blocking
            self.readSockets(connList,block)
            dataList=self.readSpecific(connList,tag,["data"],maxlen)
        return dataList

    
    def readError(self,connList=None,tag=None,errorList=["error","warning"],maxlen=0):
        """
        Read any errors from the simulation
        @param connList: socket list
        @type connList: None, List or socket.socket instance
        @param tag: The tag to look for
        @type tag: String
        @param errorList: Items counted as errors
        @type errorList: List of String
        @param maxlen: Maximum number of items to return
        @type maxlen: Int
        @return: Error messages
        @rtype: List
        """
        return self.readSpecific(connList,tag,errorList,maxlen)

    def readSpecific(self,connList=None,tag=None,errorList=["error","warning"],maxlen=0):
        """get a list of warnings/errors/data
        @param connList: socket list
        @type connList: None, List or socket.socket instance
        @param tag: The tag to look for
        @type tag: String
        @param errorList: Items counted as errors
        @type errorList: List of String
        @param maxlen: Maximum number of items to return
        @type maxlen: Int
        @return: The requested data
        @rtype: List
        """
        connList=self.parseConnList(connList)
        self.readSockets(connList)
        errList=[]
        if type(tag)!=types.NoneType and type(tag)!=types.ListType and type(tag)!=types.TupleType:
            tag=[tag]
        for data in self.recDataList:
            if data[0] in connList and data[1] in errorList:
                if tag==None or (data[2] in tag):
                    errList.append(data)
                    if maxlen>0 and len(errList)>=maxlen:
                        break
        for e in errList:
            self.recDataList.remove(e)
        return errList
    def readSockets(self,connList=None,block=0):
        """read all sockets and store results ready for later treatment...
        @param connList: socket list
        @type connList: None, List or socket.socket instance
        @param block: Whether to block when reading the sockets
        @type block: Int
        """
        connList=self.parseConnList(connList)
        if block==0:
            ready=select.select(connList,[],[],0.0)[0]
        else:
            ready=select.select(connList,[],[])[0]
        closedList=[]
        for conn in ready:#sockets with data ready to read...
            try:
                data=serialise.ReadMessage(conn.fileno())
            except:
                print "util.analyse: Error in serialise.ReadMessage"
                data=None
            if type(data)==types.NoneType:
                #print "analyse: closing connection"
                self.closeConnection(conn)
                closedList.append(conn)
            else:
                self.recDataList.append([conn]+data)
        return closedList





if __name__=="__main__":
    import sys
    portlist=[9000]
    hostname=socket.gethostname()
    if len(sys.argv)>1:
        tmp=sys.argv[1]
        try:
            port=int(tmp)
            hostname="localhost"
        except:
            hostname=tmp
    if len(sys.argv)>2:
        portlist=eval(sys.argv[2])
    if type(portlist)==type(0):
        portlist=[portlist]
    txt="data=''\nscienceList=globals().get('scienceList',[])\nfor s in scienceList:\n for i in range(len(s.thisObjList)):\n  data+=s.thisObjList[i].idstr+': '+s.strParams(i)+'\\n'\ndata+='Iter %d frametime %g batch %d\\n%s'%(ctrl.thisiter,ctrl.frametime,ctrl.batchno,ctrl.simID)\nprint data"
    if len(sys.argv)>3:
        txt="data=None\n"+sys.argv[3]
    a=analyse()
    connList=[]
    for port in portlist:
        connList.append((port,a.openConnection(hostname,port)))
    a.execute(txt,rt="data",action="now")
    for port,conn in connList:
        a.readSockets(connList=[conn],block=1)
        for data in a.recDataList:
            if len(data)>3:
                print "port %d:"%port
                print data[3]["data"]
        a.recDataList=[]
        a.closeConnection(conn)
    a.closeConnection()
