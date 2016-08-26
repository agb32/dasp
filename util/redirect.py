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

#!/usr/local/bin/python
import sys,os,socket,select
#class mysock(socket.socket):
#    def __init__(self,t1,t2):
#        socket.socket.__init__(t1,t2)
#        self.write=self.send
#        self.read=self.recv
def redirect():
    """Will attempt to open a listening port.  If succeeds, will then connect clients, and print any message that they pass.
    If can't open the listening port, it means that something else is the server.  Will then attempt to connect to this server, and redirect any print statements to the server.
    To use this, run one instance of "python redirect.py".  Then, all other modules that wish to redirect output can "import redirect; redirect.redirect()".
    """
    
    stdouthost="aipc52.phyaig.dur.ac.uk"
    stdoutport=45454
    listener=1
    if os.environ["HOSTNAME"]==stdouthost:#open listening socket
        sock=socket.socket(socket.AF_INET,socket.SOCK_STREAM)
        addr=(stdouthost,stdoutport)
        try:
            sock.bind(addr)
        except:
            listener=0
        if listener==1:
            print "I am stdout server process (if this isn't what you expected"
            print "then run python redirect.py first)."
            sock.listen(1)
            selOut=[]
            selIn=[sock]
            sdict={}
            while 1:
                rtr,rtw,err=select.select(selIn,selOut,selIn)
                for s in err:
                    #remove from dict, and act on it...
                    print "Closing",sdict[s]
                    s.close()
                    selIn.remove(s)
                    del(sdict[s])
                for s in rtr:#ready to read the socket
                    if s==sock:#new client
                        conn,raddr=s.accept()
                        selIn.append(conn)
                        sdict[conn]="("+raddr[0]+" "+str(raddr[1])+"):"
                        print "Connected from",raddr
                    else:
                        txt=s.recv(1024)
                        if len(txt)==0:
                            print "closing",sdict[s]
                            s.close()
                            selIn.remove(s)
                            del(sdict[s])
                        else:
                            print sdict[s],txt,


    else:
        listener=0
    if listener==0:
        addr=(stdouthost,stdoutport)
        print "Redirecting stdout and stderr",addr
        stdoutsock = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
        try:
            stdoutsock.connect(addr)
        except:
            print "Error connecting to host.  Could be a firewall problem?"
            raise
        #sys.stdout=stdoutsock.fileno()
        sys.stdout=stdoutsock.makefile("w",1)
        sys.stderr=sys.stdout


if __name__=="__main__":
    redirect()
