#mklmodule is a module that lets one use Intel MKL SVD and GEMM routines.

from distutils.core import setup, Extension
import sys,os.path,os,string
idnumpy=[sys.prefix+'/lib/python%d.%d/site-packages/numpy/core/include'%(sys.version_info[0],sys.version_info[1]),sys.prefix+'/include']
cont=0
if os.path.exists("/opt/intel/mkl"):
    versions=os.listdir("/opt/intel/mkl")
    versions=map(lambda x:string.split(x,"."),versions)
    versions.sort()
    if len(versions)>0:
        version=string.join(versions[-1],".")
    
        mklinclude=["/opt/intel/mkl/%s/include"%version]
        ld=[sys.prefix+'/lib']
        mkllib=["/opt/intel/mkl/%s/lib/em64t"%version]
        print "Using MKL /opt/intel/mkl/%s/lib/em64t"%version
        cont=1
if cont==0:
    print "MKL library not found - not making mklmodule"
else:
    mkl=Extension('mklmodule',
                  include_dirs=idnumpy+mklinclude,
                  library_dirs=ld+mkllib,
                  libraries=["mkl_lapack","mkl_intel_ilp64","mkl_intel_thread","mkl_core","guide","pthread"],
                  extra_compile_args=["-DMKL_ILP64"],
                  extra_link_args=["-lmkl_lapack","-lmkl_intel_ilp64","-lmkl_intel_thread","-lmkl_core","-lguide","-lpthread","-lm"],
                  sources=["mklmodule.c"]
                  )
              
    setup (ext_modules = [mkl])
