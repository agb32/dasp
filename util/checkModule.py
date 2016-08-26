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
import inspect,sys,string
"""Not yet used.  Check a science module to see if it will be compatible with the GUI - ie that the GUI can query it etc in the right way."""
def dothecheck(str):
    """Run the check
    @param str: String to be executed
    @type str: String
    @return: Tuple of txt and object
    @rtype: Tuple
    """
    print "Executing:",str
    d={}
    exec str in {},d
    return d["txt"],d["obj"]
def check(module,object):
    """Prepare string for checking module, and perform the check.
    Does (suspect) checking on the module for various methods etc to be present.
    @param module: The module name for checking
    @type module: String
    @param object: The object name for checking
    @type object: String
    @return: Error, messages
    @rtype: Tuple
    """
    err=0
    msg=""
    d={}
    str="import inspect,"+module+"; obj="+module+"."+object+"; txt=inspect.getsource("+module+"."+object+")"
    try:
        lines,obj=dothecheck(str)
    except:
        err=1
        msg+="Could not create object %s.%s\n"%(module,object)
        msg+="Possibly this may be solved by deleting all pyc objects"
    if err==0:
        args=[]
        arglist=[]
        lines=lines.split("\n")#d["txt"].split("\n")
        #obj=d["obj"]
        for line in lines:
            dpos=string.find(line,"def __init__")
            hpos=string.find(line,'#')
            if hpos>-1:
                line=line[:hpos]
            if dpos>-1 and (hpos==-1 or dpos<hpos):
                line=line[dpos+12:]
                ppos=string.find(line,'(')
                epos=string.find(line,'):')
                if ppos>-1:
                    line=line[ppos+1:]
                if epos>0:
                    line=line[:epos-1]
                print line
                line=string.split(line,',')
                for l in line:
                    t=string.split(l,'=')
                    if len(t)>1:
                        args.append(None)
                    arglist.append(t[0])
                break
        print "Argument list:",arglist
        obj=obj(*args)
        try:
            userargs=obj.getInitArgs()
        except:
            err=1
            msg+="getInitArgs method not found\n"
        if err==0:
            userargkeys=userargs.keys()
            print userargs
            for arg in arglist:
                if arg!="self" and arg not in userargkeys:
                    msg+= "ERROR:  arg %s not specified in getInitArgs\n"%arg
                    err=1
        try:
            paramlist=obj.getParams()
        except:
            err=1
            msg+="getParams method not found\n"
        try:
            intype=obj.getInputType()
        except:
            err=1
            msg+="getInitArgs method not found\n"
        try:
            outtype=obj.getOutputTypes()
        except:
            err=1
            msg+="getInitArgs method not found\n"


    if err:
        print "\nERRORS:"
        print msg
        print "Module has failed the test - please amend\n"
    else:
        print "\nCongratulations - module has passed the test."
        print "However, this doesn't gaurantee that it will work, just increases"
        print "the chances of doing so!\n"
    return err,msg
if __name__=="__main__":
    check(sys.argv[1],sys.argv[2])
