# science/recorder.py
#
# aosim recording module
#
# Author:   Nazim Bharmal
# Date:     Jan/2014
# Version:  not for production
# Platform: Darwin
# Python:   2.6.7
# Numpy:    1.5.1

# The basic concept is to open a new file upon reaching the StartNiter,
# and then append until reaching the EndNiter

import template
numpy=template.numpy
import util.FITS

class aosimRecorder(template.aosimNessuno):
    """aosim recorder, saving data to a file
    [todo, epydoc compatible documentation]
    """
    def _doInitializationCode(self):
           # \/ load configuration variables
        self.active=self.config.getVal(
              "active",
              default=1, raiseerror=0 ) # whether active or not
        self.separate=self.config.getVal(
              "recordSeparated",
              default=0, raiseerror=0 ) # write separated or a combined file(s)
        self.startnValid=int(self.config.getVal(
              "recordStartNiter",
              default=1, raiseerror=0 )) # when to start recording (inclusive)
        self.endnValid=int(self.config.getVal(
              "recordEndNiter",
              default=1, raiseerror=0 )) # when to end recording (exclusive)
        if self.startnValid>=self.endnValid:
              errMsg="Require EndNiter>StartNiter"
              if self.active:
                  raise ValueError( "recorder: "+errMsg )
              else:
                  print "WARNING:(**recorder_"+str(self.idstr)+"**):"+errMsg
        self.fnameStem=self.config.getVal(
              "recordFilenameStem",
              default="aosim_recorder_{0:s}-".format(str(self.parent.idstr)), 
              raiseerror=0 ) # stem for the filename to record to
        if self.separate:
            self._formFileName=lambda self,n :\
               "{0:s}{1:d}.spfits".format(self.fnameStem,n)
        else:
            self._formFileName=lambda self,n :\
               "{0:s}.fits".format(self.fnameStem)
        self._checkFilesExistence()
        if self.active:
           print("INFORMATION(**recorder_"+str(self.idstr)+"**):"+
               "from %s upto %s"%( str(self.startnValid),str(self.endnValid)) )
  

    def _checkFilesExistence(self):
        import os
        def __test4path(self,i):
           fn=self._formFileName(self,i)
           if os.path.exists( fn ):
               print("WARNING(**recorder_"+str(self.idstr)+"**): file '"+fn
                        +"' exists, will be appended to.")
        if self.separate:
           for i in range(self.startnValid,self.endnValid):
              __test4path(self,i)
        else:
           __test4path(self,None)
   
    def _doSomeProcessing(self,ip):
        if self.nValid>=self.startnValid and self.nValid<self.endnValid:
               # \/ thisFrac will equal zero when reaching (approximately)
               #  every 10% of the total number
            thisFrac=( (self.nValid-self.startnValid)
                        %int(0.1*(self.endnValid-self.startnValid)) )
            if thisFrac==0: 
               print("INFORMATION:(**recorder_"+str(self.idstr)
                  +"**) recording ("
                  +str(int(
                      (100*(self.nValid-self.startnValid))
                     /(self.endnValid-self.startnValid) ))
                  +"%)")
            if not self.sparse:
               util.FITS.Write( ip, self._formFileName(self,self.nValid), 
                    "RECORDER= "+str(self.nValid), "a" )
            else:
               util.FITS.saveSparse( ip, self._formFileName(self,self.nValid),
                     "a", "RECORDER= "+str(self.nValid) )
