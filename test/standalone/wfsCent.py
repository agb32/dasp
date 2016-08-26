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
"""Example of using the wfsCent class outside of the simulation.
The example here uses a aobase instance as a 'null' instance.
This is the parent to the wfsCent instance.
Setting the appropriate parameters in aobase allows for the right
format of data to enter wfsCent, and then be processed as expected.
The output data can then be extracted from wfsCent, and visualized.
If wfs_rowint is defined, this code could check that output is
correct as per a rolling shutter.
"""
from __future__ import print_function

import base.aobase
import base.readConfig
import numpy
import science.wfscent
import sys
import util.Ctrl
try:
   import pylab
except ImportError:
   no_pylab=None

orig_stdout=sys.stdout
sys.stdout=util.Ctrl.myStdout(314519) # stdout = attractive

## start of instance creation and setup

config=base.readConfig.AOXml("wfsCent.xml") 
aobase=base.aobase.aobase(None,config,idstr="a")
aobase.dataValid=1
wfscent=science.wfscent.wfscent(aobase,config,idstr="a")
wfscent.initialise(aobase,"wfs")
wfscent.finalInitialisation()
wfscent.wfscentObj=wfscent.thisObjList[0].wfscentObj

## end of instance creation and setup

sys.stdout=orig_stdout # reset stdout

print("\n"+"-"*79)

orig_stdout=sys.stdout
sys.stdout=util.Ctrl.myStdout(314519) # stdout = attractive
wfs_npix=config.getVal("wfs_n")
wfs_nsub=config.getVal("wfs_nsubx")
sys.stdout=orig_stdout # reset stdout
#
i=0
   # \/ input data
aobase.outputData=numpy.zeros([wfs_npix*wfs_nsub,wfs_npix*wfs_nsub])
while not wfscent.dataValid:
   wfscent.generateNext()
   i+=1
wfs_int,wfs_tstep=[ config.getVal(var) for var in "wfs_int","tstep" ]
assert i==( wfs_int/wfs_tstep ), "Wrong number of integrations?"
print("Using {0:d} integrations per calculation".format(i))
wfs_offsets=wfscent.outputData.copy()
   # /\ output data
#
# Check validity/invalidity is okay
wfscent.generateNext()
for dummy in xrange(1,i):
   assert wfscent.dataValid==0, "Unexpected validity"
   wfscent.generateNext()
assert wfscent.dataValid, "Unexpected invalidity"
#
wfscent_op=wfscent.outputData-wfs_offsets
#assert wfscent_op.var()==0.0, "Expected zero"
print("Zero input:: o/p variance={0:5.3g}".format(
      wfscent_op.var()))
#
aobase.outputData=numpy.ones([wfs_npix*wfs_nsub,wfs_npix*wfs_nsub])
for dummy in xrange(i):
   wfscent.generateNext()
#
assert wfscent.dataValid, "Unexpected invalidity"
wfscent_op=wfscent.outputData-wfs_offsets
#assert wfscent_op.var()<1e-8, "Expected zero" # there will be noise
print("Ones input:: o/p variance={0:5.3g}".format(
      wfscent_op.var()))
#
# make a gradient
aobase.outputData=numpy.add.outer(
   numpy.linspace(-1,1,wfs_npix*wfs_nsub), numpy.zeros(wfs_npix*wfs_nsub) )
for dummy in xrange(i):
   wfscent.generateNext()
#
assert wfscent.dataValid, "Unexpected invalidity"
wfscent_op=wfscent.outputData-wfs_offsets
assert wfscent_op[:,:,0].var()<wfscent_op[:,:,1].var(), "Unexpected signal"
print("Grad input:: o/p variance={0[0]:5.3g}/{0[1]:5.3g}".format(
      wfscent_op.reshape([-1,2]).var(axis=0)))
#
aobase.outputData=aobase.outputData.T
for dummy in xrange(i):
   wfscent.generateNext()
#
assert wfscent.dataValid, "Unexpected invalidity"
wfscent_op=wfscent.outputData-wfs_offsets
assert wfscent_op[:,:,0].var()>wfscent_op[:,:,1].var(), "Unexpected signal"
print("Grad^T input:: o/p variance={0[0]:5.3g}/{0[1]:5.3g}".format(
      wfscent_op.reshape([-1,2]).var(axis=0)))
if 'no_pylab' in dir():
   print("(A) Need analysis code here")
else:
   for k in (0,1):
      pylab.subplot(3,2,5+k)
      pylab.title("j=n/a,k={0:d}".format(k))
      pylab.imshow( wfscent_op[:,:,k] )
   
#
# try rolling shutter::
orig_stdout=sys.stdout
sys.stdout=util.Ctrl.myStdout(314519) # stdout = attractive
#
try:
   wfs_rowint=config.getVal("wfs_rowint")
except Exception:
   wfs_rowint=None
#
sys.stdout=orig_stdout

if wfs_rowint!=None:
   r=wfs_rowint/wfs_tstep
   l=wfs_int/wfs_tstep-r
   assert not r%1 and not l%1, "fraction r/l"
   r,l=int(r),int(l)
   print("r,l,intT={0:d},{1:d},{2:d}".format(r,l,r+l))
   #
   # algorithm to predict start integration is then::
   # if ( integrationNum<int(row*(wfs_nsub-1)**-1.0*l+0.5) or
   #      integrationNum>int(row*(wfs_nsub-1)**-1.0*l+r-0.5) )):
   #   not using this integration step
   # else
   #   use this integration step
   #
   # test by alternating sign of input at halfway point of integrations
   wfscent_ops=[]
   aobase.outputData=numpy.add.outer(
         numpy.linspace(-1,1,wfs_npix*wfs_nsub),
         numpy.linspace(-1,1,wfs_npix*wfs_nsub) )
   print("Looping:",end="") ; sys.stdout.flush()
   assert wfscent.dataValid, "Unexpected invalidity"
   print("|",end="") ; sys.stdout.flush()
   for manequin in xrange((r+l)-1):
      if manequin>=l:
         aobase.outputData*=-1
         print("~",end="") ; sys.stdout.flush()
      else:
         print("=",end="") ; sys.stdout.flush()
      wfscent.generateNext()
      assert not wfscent.dataValid, "Unexpected validity"
   wfscent.generateNext()
   print("|",end="") ; sys.stdout.flush()
   wfscent_op=wfscent.outputData-wfs_offsets
   wfscent_ops.append( wfscent_op-wfs_offsets )
   #
   assert wfscent.dataValid, "Unexpected invalidity"
   for manequin in xrange((r+l)-1):
      if manequin<=l:
         aobase.outputData*=-1
         print("~",end="") ; sys.stdout.flush()
      else:
         print("=",end="") ; sys.stdout.flush()
      wfscent.generateNext()
      assert not wfscent.dataValid, "Unexpected validity"
   wfscent.generateNext()
   assert wfscent.dataValid, "Unexpected invalidity"
   print("|",end="") ; sys.stdout.flush()
   wfscent_op=wfscent.outputData-wfs_offsets
   wfscent_ops.append( wfscent_op-wfs_offsets )
   
   print("(done)") ; sys.stdout.flush()
   if 'no_pylab' in dir():
      print("(B) Need analysis code here")
   else:
      for j in (0,1):
         for k in (0,1):
            pylab.subplot(3,2,1+j+k*2)
            pylab.title("j={0:d},k={1:d}".format(j,k))
            pylab.imshow( wfscent_ops[j][:,:,k] )

