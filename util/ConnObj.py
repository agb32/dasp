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

import cPickle,os,socket,serialise,types
class ConnMsg:
    """Class to store commands sent to the simulation from outside (via
    socket).  Simulation programmer does not need to know the details of
    this class

    Class variables (not important to simulation programmer):
     - action - string, the command to be send, one of "now","rpt","rptN" (integer N), "cmd" or "del".
     - command - string, the command to be exec'd at the right time
     - ret - string, the return value created by exec'ing the command
     - tag - string, the tag associated with the data.
     - dataList - list, a list that can be unpickled to fill a ConnMsg object
    @cvar command: The command to be exec'd at the correct time
    @type command: String
    @cvar action: Determine when to exec the command, one of "now" (exec
    immediately),"rpt" (repeat at the end of every iteration), "rptN" (integer
    N, repeat at the end of every N iterations), "cmd" (exec at the end of the
    #current iteration) or "del" (delete a rpt or cmd action).
    @type action: String
    @cvar ret: The name of the value to be returned (created by the exec).
    @type ret: String
    @cvar tag: A tag for the returned data
    @type tag: String
    """
    
    #"""Possible actions are "now", "rpt", "cmd", "del".
    #now will execute the command immediately, rpt will execute every
    #cycle, cmd will execute once during the break, and del will delete
    #a scheduled command.
    #commmand is a text to be exec'd.
    #ret is None (everything) or a list of things to be returned
    #tag can be used for identification purposes.
    #"""
    def __init__(self,command,action,ret=None,tag=None,dataList=None,maxRetry=10):
        """Initialise, and unpickle the dataList if not None
        @param command: The command to be exec'd at the correct time
        @type command: String
        @param action: Determine when to exec the command, one of "now" (exec
        immediately),"rpt" (repeat at the end of every iteration), "rptN" (integer
        N, repeat at the end of every N iterations), "cmd" (exec at the end of the
        #current iteration) or "del" (delete a rpt or cmd action).
        @type action: String
        @param ret: The name of the value to be returned (created by the exec).
        @type ret: String
        @param tag: A tag for the returned data
        @type tag: String
        @param dataList: A list, length 4, which can be unpickled.
        @type dataList: List
        """
        self.command=command
        self.action=action
        self.maxRetry=maxRetry
        self.retryCnt=0
        if ret!=None and type(ret)!=types.ListType:
            self.ret=[ret]
        else:
            self.ret=ret
        self.tag=tag
        if dataList!=None:
            self.unpickle(dataList)
    def pickle(self):
        """Create a list of strings from self.
        @return: A list of [command, action, ret, tag]
        @rtype: List
        """
        return [self.command,self.action,self.ret,self.tag]
    def unpickle(self,data):
        """Populate self from a list of length 4.
        @param data: A list to use to populate self
        @type data: List
        """
        self.action=data[0]
        self.command=data[1]
        if len(data)>2:
            self.ret=data[2]
        else:
            self.ret=None
        if len(data)>3:
            self.tag=data[3]
        else:
            self.tag=None
    
class ConnObj:
    """A connection object, used to connect to a simulation.

    Class variables (important to simulation programmer):
     - port - int, port number to connect to
     - host - string, hostname to connect to
    Class variables (not important to simulation programmer):
     - sock - the socket object
    @cvar port: Port number to connect to
    @type port: Int
    @cvar host: Host name to connect to
    @type host: String
    @cvar sock: Internet socket
    @type sock: socket.socket object
    """
    
    def __init__( self, port, host=socket.gethostname() ):
        """Open a connection with the server
        @param port: Port to connect to
        @type port: Int
        @param host: Hostname to connect to
        @type host: String
        """
        self.port=port
        self.host=host
        self.sock=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.connect((host,port))
    def send(self,connMsg):
        """Send a command to the server.  Server may reply
        (if reply requested) or may not.
        @param connMsg: The command to be sent
        @type connMsg: ConnMsg object"""
        #str=cPickle.dumps(connMsg,protocol=cPickle.HIGHEST_PROTOCOL)
        #serialise.Send([str],self.sock)
        serialise.Send(connMsg.pickle(),self.sock)

    def recv(self):
        """Receive data from a socket.  
        Expecting data=[info,tag,whatever]
        if info=="data" whatever will be a dictionary of returned data.
        if info=="warning", whatever will be a string warning message.
        tag is the tag sent, or None if not obtainable.
        @return: The data received from the socket
        @rtype: List
        """
        data=serialise.ReadMessage(self.sock.fileno())
        #print data
        
        return data

        
if __name__=="__main__":
    """test communications... and then quit the simulation"""
    connMsg=ConnMsg("atmos=Atmos.outputData","cmd")
    connObj=ConnObj(9000)
    connObj.send(connMsg)
    data=connObj.recv()

    if data[0]=="data":
        print "Got data, tag, shape, type:",data[0],data[1],data[2]["atmos"].shape,type(data[2]["atmos"]),data[2]["atmos"].typecode()
    elif data[0]=="warning":
        print "Got warning:",data[0],data[1],data[2]
    else:
        print "Unknown response:",data
    print "Now quitting the simulation:"
    connMsg.command="ctrl.go=0"
    connMsg.ret=[]
    connObj.send(connMsg)
    print "Quit"
