from distutils.core import setup, Extension
import sys,os.path
id=[sys.prefix+'/include/python%d.%d/Numeric'%(sys.version_info[0],sys.version_info[1]),sys.prefix+'/include']
idnumpy=[sys.prefix+'/lib/python%d.%d/site-packages/numpy/core/include'%(sys.version_info[0],sys.version_info[1]),sys.prefix+'/include']
ld=[sys.prefix+'/lib']
fft=Extension('fftmodule',
              include_dirs=idnumpy,
              library_dirs=ld,
              libraries=["pthread","fftw3f_threads","fftw3f"],
              extra_compile_args=[],
              extra_link_args=["-lpthread","-lfftw3f_threads","-lm"],#,"-lfftw3f"
              sources=["fftmodule.c"]
              )
              
cent = Extension('centmodule',
                 include_dirs = idnumpy,
		library_dirs=ld,
#			runtime_library_dirs=['/usr/local/lib'],
		libraries=["fftw3f"],
                extra_compile_args=["-pthread"],
		extra_link_args=["-lfftw3f",'-lgsl','-lgslcblas','-lm','-pthread'],
	        sources = ['centmodule.c']
		)
binimg=Extension('binimgmodule',
		include_dirs=idnumpy,
		sources=['binimgmodule.c'],
		library_dirs=ld,
		extra_link_args=['-lm'],
		)
imgnoise=Extension('imgnoisemodule',
		include_dirs=idnumpy,
		sources=['imgnoisemodule.c'],
		library_dirs=ld+[os.path.realpath('..'),os.path.realpath('.')],
		extra_link_args=['-lcrecipes','-lm'],
		)


utils=Extension('utilsmodule',
		include_dirs=idnumpy,
		sources=['utils.c'],
		library_dirs=ld,
		extra_link_args=['-lm'],
                extra_objects = ['mvm.o']
		)
sor=Extension('sormodule',
              include_dirs=idnumpy,
              sources=['sormodule.c'],
              library_dirs=ld,
#              extra_link_args=['-lm'],
              )
interp=Extension('interpmodule',
                 include_dirs=idnumpy,
                 sources=['interpmodule.c'],
                 library_dirs=ld+[os.path.realpath('..'),os.path.realpath('.')],
                 extra_link_args=['-lcrecipes','-lgsl','-lgslcblas','-lm'],
                 extra_objects = ['interpolate.o']
                 )
phaseCov=Extension('phaseCovmodule',
                 include_dirs=idnumpy,
                 sources=['phaseCovmodule.c'],
                 library_dirs=ld+[os.path.realpath('..'),os.path.realpath('.')],
                   extra_compile_args=["-pthread"],
                 extra_link_args=['-lcrecipes','-lgsl','-lgslcblas','-lm','-lpthread'],
                 )
zernike=Extension('zernikemodule',
                 include_dirs=idnumpy,
                 sources=['zernikemodule.c','josesubs.c'],
                 library_dirs=ld,
                 extra_link_args=['-lm'],
                 )
xpoke=Extension('xpokemodule',
                 include_dirs=idnumpy,
                 sources=['xpokemodule.c','josesubs.c'],
                 library_dirs=ld,
                 extra_link_args=['-lm'],
                 )
psfparams=Extension('psfparamsmodule',
                 include_dirs=idnumpy,
                 sources=['psfparamsmodule.c','josesubs.c'],
                 library_dirs=ld,
                 extra_link_args=['-lm'],
                 )
                 

setup (ext_modules = [fft,cent,binimg,imgnoise,utils,sor,interp,phaseCov,zernike,xpoke,psfparams])
