# Find python paths
#
# $Id: pypath.py,v 1.2 2005/11/17 13:40:55 ali Exp $
# 
if __name__=="__main__":
   import sys

   if sys.argv[1]=="base":
      print sys.prefix,
   elif sys.argv[1]=="version":
      print "%d.%d" % (sys.version_info[0],sys.version_info[1])
   elif sys.argv[1]=="lib":
      pathstr=""
      for x in sys.path:
         pathstr+=" -I%s" % (x)
      print pathstr,
   else:
      print "Unknown option"
      sys.exit(1)
