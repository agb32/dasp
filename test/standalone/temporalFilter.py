"""Example of using the temporalFilter class outside of the simulation.
  The example here uses a aobase instance as a 'null' instance.
  This is the parent to the temporalFilter instance.

  There are no necessary parameters in the parameter file: they are added
  here.
"""
from __future__ import print_function

import base.aobase
import base.readConfig
import numpy
import science.temporalFilter
import sys
import util.Ctrl
try:
   import pylab
except ImportError:
   no_pylab=None

orig_stdout=sys.stdout
sys.stdout=util.Ctrl.myStdout(999) # stdout = attractive

## start of instance creation and setup

config=base.readConfig.AOXml("temporalFilter.xml") 
# \/ test with a simple alternating input
#
config.setVal( "filterM", [[0.25,0.25,0.25,0.25],[0,0,0,0]] )
aobase=base.aobase.aobase(None,config,idstr="a")
aobase.initialise(aobase,"aobase1")
aobase.finalInitialisation()
aobase.dataValid=1
tFilt=science.temporalFilter.aosimTemporalFilter(aobase,config,idstr="a")
tFilt.initialise(aobase,"tfilt1")
tFilt.finalInitialisation()

sys.stdout=orig_stdout # reset stdout
print("\n"+"-"*79)
## end of instance creation and setup
#
print("Starting basic test...")
   # \/ input data
filterM = numpy.array( config.getVal("filterM") )
aobase.outputData = numpy.ones( filterM.shape )

for opExpectation in (0.25,0,0.75,0,1,0):
   tFilt.generateNext()
   #
   assert tFilt.dataValid==1, "Unexpected invalidity"
   assert tFilt.outputData.var()==0, "Variance not zero"
   assert tFilt.outputData.mean()==opExpectation,\
         "Differing value from expectation:"+str(opExpectation)
   assert tFilt.outputData.shape==( filterM.shape[1], ), "Wrong shape"

print("...complete")


   # \/ reinstantiate the filter
sys.stdout=util.Ctrl.myStdout(998) # stdout = attractive
tFilt=science.temporalFilter.aosimTemporalFilter(aobase,config,idstr="a")
tFilt.initialise(aobase,"tfilt2")
tFilt.finalInitialisation()
sys.stdout=orig_stdout # reset stdout
print("\n"+"-"*79)
#
print("Starting data invalid in parent test...")
   # \/ input data
filterM = numpy.array( config.getVal("filterM") )
aobase.outputData = numpy.ones( filterM.shape )

for counter, opExpectation in enumerate([0.25,0,0.5,0,0.75,0,1,0,1,0,1]):
   aobase.dataValid = 1-counter%2 # flip-flop the data validity
   tFilt.generateNext()
   #
   assert tFilt.dataValid==1, "Unexpected invalidity"
   assert tFilt.outputData.var()==0, "Variance not zero"
   assert tFilt.outputData.mean()==opExpectation,\
         "Differing value from expectation:"+str(opExpectation)
   
print("...switching...")

   # \/ new input data
for counter, opExpectation in enumerate([0.25,0,0.5,0,0.75,0,1,0,1,0,1]):
   ipHistory= [numpy.random.normal(0,1,size=filterM.shape[1]) for
         dummy in range(filterM.shape[0]) ]
   aobase.outputData = numpy.array(ipHistory)
   aobase.dataValid = 1-counter%2 # flip-flop the data validity
   tFilt.generateNext()
   #
   assert tFilt.dataValid==1, "Unexpected invalidity"
   if counter%2==0:
      assert tFilt.outputData.var()==0, "Variance not zero"
   else:
      assert (tFilt.outputData.var()>0 and
              tFilt.outputData.var()<2), "Variance too big?"

print("...complete")

