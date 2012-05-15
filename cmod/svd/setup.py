from distutils.core import setup, Extension
import sys,os.path
id=[sys.prefix+'/include/python%d.%d/Numeric'%(sys.version_info[0],sys.version_info[1]),sys.prefix+'/include']
idnumpy=[sys.prefix+'/lib/python%d.%d/site-packages/numpy/core/include'%(sys.version_info[0],sys.version_info[1]),sys.prefix+'/include']
ld=[sys.prefix+'/lib']
svd = Extension('svdmodule',
                 include_dirs = idnumpy+["."],#["SVDLIBC"],
		library_dirs=ld+["."],#+["SVDLIBC"],
#			runtime_library_dirs=['/usr/local/lib'],
#		libraries=["fftw3f"],
#                extra_compile_args=["-pthread"],
		extra_link_args=["-lsvd"],#,'-lgsl','-lgslcblas','-lm','-pthread'],
	        sources = ['svdmodule.c']
		)
 

setup (ext_modules = [svd])
