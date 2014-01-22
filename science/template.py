# science/???.py
#
# aosim ???
#
# Author:   ???
# Date:     ???/201?
# Version:  not for production
# Platform: ??? (where it has been tested, either Darwin or Linux---usually)
# Python:   2.?.? (say which version of Python it was tested on)
# Numpy:    1.?.? (say which version of numpy it was tested with)

# Lines marked with '#?' are those that may be important but don't appear
# to be necessary.
#
# This code based on science/WFSspatialFilter.py

import base.aobase
import numpy
import time


class aosimNessuno(base.aobase.aobase):
    """aosim ???
    [todo, epydoc compatible documentation]
    """
    def __init__(self,parent,config,args={},forGUISetup=0,passThrough=False,debug=None,idstr=None):
        """ [todo, epydoc compatible documentation]
        """
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)

        self.active=self.config.getVal( "active", default=0, raiseerror=0 ) 
           # \/ load configuration variables
#(example)        self.yyy=self.config.getval(
#(example)              "nessuno",
#(example)              default=none,
#(example)              raiseerror=0 ) 
           # \/ sparse or not?
        self.sparse=self.config.getVal( "sparse", default=0, raiseerror=0 )
        self.nValid=0

           # \/ do initialisation calculations
        self.obj=self._doInitializationCode()

    def _doInitializationCode(self):
        pass
   
    def generateNext(self,msg=None):
        """[todo, epydoc compatible documentation]
        """
        if self.generate==1:
            if self.newDataWaiting:#this is always 1
                if self.parent.dataValid==1:
                    self.inputInvalid=0
                else:
                    self.inputInvalid=1
                    self.dataValid=0
            if self.inputInvalid==0: # there was an input => do stuff!
                self.dataValid=1
                self.nValid+=1
                if not self.active:
                   self.outputData=self.parent.outputData
                else:
                   # apply the relevant code
                   self.outputData=self._doSomeProcessing(
                          self.parent.outputData )
    
    def _doSomeProcessing(self,ip):
        return ip

    def plottable(self,objname="$OBJ"):
        """[todo] write epydoc compatible description"""
        op=""
        if not self.sparse:
           cmd="%s.outputData"%(objname)
        else:
           cmd="%s.outputData.todense()"%(objname)
        op+=('<plot title="ip or op" '+
             'cmd="data=%s" ret="data" type="pylab" '+
             'when="rpt" palette="jet"/>\n')%(cmd)
             # /\ options other than once include: once

        op+=('<plot title="On/off" when="cmd" ret="feedback" texttype="1" '+
             'wintype="mainwindow">'+
             '<cmd>%s.active=1-self.active</cmd></plot>')%(objname)
        return(op)
