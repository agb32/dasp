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
#from __future__ import print_function
# Find python paths
#
# $Id: pypath.py,v 1.2 2005/11/17 13:40:55 ali Exp $
#
if __name__=="__main__":
   import sys
   
   op="" ; term=sys.argv[1]
   if term=="base":
      op+=sys.prefix
   elif term=="version":
      op+="%d.%d"%sys.version_info[:2]#"{0[0]:d}.{0[1]:d}".format( sys.version_info) - undid the nab obfuscation.
   elif term=="lib":
      for x in sys.path:
         op+=" -I%s" % (x)
   elif term=="site-packages":
      from distutils.sysconfig import get_python_lib;
      op+=get_python_lib()
   elif term=="numpy":
      import numpy
      op=numpy.get_include()
   else:
      print "Unknown option"
      sys.exit(1)
   print op
