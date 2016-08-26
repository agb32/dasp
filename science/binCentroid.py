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
"""A module to bin the x and y centroid signals.  Parent should be
wfscent, child should be a reconstructor.

Basically bins/interpolates down/up the centroid arrays to a certain
specified size.  Can be useful eg for extracting tip/tilt signals or similar.

No point being resource sharing, since about the only array is the
output, which is passed straight to reconstructor and might overwrite
each other if sharing.

Note, vdmUser does not do the same thing, since it relies entirely on
interpolating, which means that averaging all x and y centroids for
example gives the wrong result.

"""

import cmod.binimg
import numpy
import base.aobase
class BinCentroid(base.aobase.aobase):
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        """Important parameters here are the nsubx for the output.  Other parameters can simply be taken from the input array size.
        """
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.nsubx=self.config.getVal("wfs_nsubx")
        if forGUISetup==1:
            self.outputData=[(self.nsubx,self.nsubx,2),numpy.float32]
        else:
            self.outputData=numpy.zeros((self.nsubx,self.nsubx,2),numpy.float32)
            self.pupil=self.config.getVal("pupil")
            self.binsubflag=numpy.zeros((self.nsubx,self.nsubx),numpy.float32)

    def finalInitialisation(self):
        #Count the number of valid subaps in each bin.
        self.parentSubflag=self.pupil.getSubapFlag(nsubx=self.parentList[0].outputData.shape[0]).astype(numpy.float32)
        cmod.binimg.binimg(self.parentSubflag,self.binsubflag)
    def generateNext(self):
        """Compression main loop - compress along a line of sight to get DM actuator values."""
        if self.generate==1:
            if self.newDataWaiting:
                if self.parent.dataValid==1:
                    self.dataValid=1
                    self.binCentroids(self.parent.outputData)
                else:
                    print "binCentroid waiting for data from parent, but not valid"
                    self.dataValid=0
        else:
            self.dataValid=0

    def binCentroids(self,inputData):
        cmod.binimg.binimg(inputData[:,:,0],self.outputData[:,:,0])
        cmod.binimg.binimg(inputData[:,:,1],self.outputData[:,:,1])
        #now take the average.
        self.outputData[:,:,0]/=self.binsubflag
        self.outputData[:,:,1]/=self.binsubflag
        
    def plottable(self,objname="$OBJ"):
        """Return a XML string which contains the commands to be sent
        over a socket to obtain certain data for plotting.  The $OBJ symbol
        will be replaced by the instance name of the object - e.g.
        if scrn=mkscrns.Mkscrns(...) then $OBJ would be replaced by scrn."""
        if self.idstr[0]==None or self.idstr[0]=="":
            id=""
        else:
            id=" (%s)"%self.idstr[0]
        txt="""<plot title="binCentroid x%s" cmd="data=%s.outputData[:,:,0]" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
        txt+="""<plot title="binCentroid y%s" cmd="data=%s.outputData[:,:,1]" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
        txt+="""<plot title="binCentroid subflag y%s" cmd="data=%s.binsubflag" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
        return txt
	
