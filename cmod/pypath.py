from __future__ import print_function
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
      op+="{0[0]:d}.{0[1]:d}".format( sys.version_info)
   elif term=="lib":
      for x in sys.path:
         op+=" -I%s" % (x)
   elif term=="site-packages":
      from distutils.sysconfig import get_python_lib;
      op+=get_python_lib()
   else:
      print("Unknown option")
      sys.exit(1)
   print(op,end="")
