#mklmodule is a module that lets one use Intel MKL SVD and GEMM routines.
import traceback
from distutils.core import setup, Extension
import sys,os.path,os,string
idnumpy=[sys.prefix+'/lib/python%d.%d/site-packages/numpy/core/include'%(sys.version_info[0],sys.version_info[1]),sys.prefix+'/include']
cont=0
if os.path.exists("/opt/intel/mkl"):
    versions=os.listdir("/opt/intel/mkl")
    versions=map(lambda x:string.split(x,"."),versions)
    v=[]
    for vv in versions:
        if len(vv)==4:
            v.append(vv)
    versions=v
    versions.sort()
    if len(versions)>0:
        version=string.join(versions[-1],".")
    else:
        version=""
    mklinclude=["/opt/intel/mkl/%s/include"%version]
    ld=[sys.prefix+'/lib']
    mkllib=["/opt/intel/mkl/%s/lib/em64t"%version]
    if not os.path.exists(mkllib[0]):
        mkllib=["/opt/intel/mkl/%s/lib/intel64"%version]
    if os.path.exists(mkllib[0]):
        print "Using MKL %s"%mkllib[0]
        cont=1
        if os.path.exists(mkllib[0]+"/libmkl_lapack.so"):
            lapack="mkl_lapack"
        elif os.path.exists(mkllib[0]+"/libmkl_lapack95_ilp64.so"):
            lapack="mkl_lapack95_ilp64.so"
        else:
            cont=0
            lapack=None
        libs=[lapack,"mkl_intel_ilp64","mkl_intel_thread","mkl_core","guide","pthread"]
        linkArgs=["-l%s"%lapack,"-lmkl_intel_ilp64","-lmkl_intel_thread","-lmkl_core","-lguide","-lpthread","-lm"]
    else:
        cont=0
if cont==0:
    print "MKL library not found - not making mklmodule"
else:
    mkl=Extension('mklmodule',
                  include_dirs=idnumpy+mklinclude,
                  library_dirs=ld+mkllib,
                  libraries=libs,
                  extra_compile_args=["-DMKL_ILP64"],
                  extra_link_args=linkArgs,
                  sources=["mklmodule.c"]
                  )
              
    setup (ext_modules = [mkl])
files=os.listdir("/opt")
acmllist=[]
for f in files:
    if f.startswith("acml"):
        acmllist.append(f)
cont=0
acmllist.sort()
if len(acmllist)>0 and os.path.exists("/opt/%s"%acmllist[-1]):
    #versions=os.listdir("/opt/intel/mkl")
    #versions=map(lambda x:string.split(x,"."),versions)
    #versions.sort()
    #if len(versions)>0:
    #    version=string.join(versions[-1],".")
    if 1:
        #mklinclude=["/opt/intel/mkl/%s/include"%version]
        acmlinclude=["/opt/%s/gfortran64_mp/include/"%acmllist[-1]]
        ld=[sys.prefix+'/lib']
        #mkllib=["/opt/intel/mkl/%s/lib/em64t"%version]
        acmllib=["/opt/%s/gfortran64_mp/lib/"%acmllist[-1]]
        print "Using ACML /opt/%s/gfortran64_mp/lib/"%acmllist[-1]
        cont=1
        print "You now should export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/%s/gfortran64_mp/lib if this isn't already on the path"%acmllist[-1]
if cont==0:
    print "ACML library not found - not making acmlmodule"
else:
    acml=Extension('acmlmodule',
                  include_dirs=idnumpy+acmlinclude,
                  library_dirs=ld+acmllib,
                  libraries=["acml_mp","pthread"],
                   #extra_compile_args=["-DMKL_ILP64"],
                  extra_link_args=["-lacml_mp","-lpthread","-lm"],
                  sources=["acmlmodule.c"]
                  )
    try:
        setup (ext_modules = [acml])
    except:
        traceback.print_exc()
        print "UNABLE TO COMPILE ACML MODULE... CONTINUING"

cont=0
atpath=None
if os.path.exists("/usr/lib64/atlas/libatlas.so"):
    atpath="/usr/lib64/atlas/"
elif os.path.exists("/usr/lib/libatlas.so"):
    atpath="/usr/lib/"
if atpath!=None:
    #versions=os.listdir("/opt/intel/mkl")
    #versions=map(lambda x:string.split(x,"."),versions)
    #versions.sort()
    #if len(versions)>0:
    #    version=string.join(versions[-1],".")
    if 1:
        #mklinclude=["/opt/intel/mkl/%s/include"%version]
        atlaslib=[atpath]
        ld=[sys.prefix+'/lib']
        #mkllib=["/opt/intel/mkl/%s/lib/em64t"%version]
        atlasinclude=["/usr/include/","/usr/include/atlas"]
        print "Using atlas/lapack"
        cont=1
if cont==0:
    print "atlas/lapack library not found - not making atlasmodule"
else:
    atlas=Extension('atlasmodule',
                  include_dirs=idnumpy+atlasinclude,
                  library_dirs=ld+atlaslib,
                  libraries=["atlas","pthread"],
                   #extra_compile_args=["-DMKL_ILP64"],
                  extra_link_args=["-latlas","-lptcblas","-lclapack","-lpthread","-lm"],
                  sources=["atlasmodule.c"]
                  )
              
    try:
        setup (ext_modules = [atlas])
    except:
        traceback.print_exc()
        print "UNABLE TO COMPILE ATLAS MODULE... CONTINUING"
