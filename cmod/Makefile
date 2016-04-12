# Makefile for aosim
# Ali Bharmal <n.a.bharmal@durham.ac.uk>	13/06/2005
#
# $Id: Makefile,v 1.7 2011/12/06 05:14:39 ali Exp $
#
# Updated 2012Aug09 by Urban Bitenc.
#

install:
	echo "Making rww modules"
	(cd src/ ; make install )
	echo "Making shmem module"
	(cd shmem/ ; make install )
#	echo "Making svd module"
#	(cd svd/ ; make install)
# commented out by UB, 2012 Aug 8th:
# (fails to build)
	echo "Making MKL/AMD module"
	(cd mkl ; python setup.py install --install-lib=../)
clean:
	echo "Cleaning rww modules"
	(cd src/ ; make clean )
	echo "Cleaning shmem module" 
	(cd shmem/ ; make clean )
	echo "Cleaning svd module"
	(cd svd/; make clean)

#	UB, 2012 Aug 09:
#	The next line (rm *.so) removes the .so files that were copied to cmod/ by the
#	"/src/Makefile: python setup.py install --install-lib=../" comand.
#	setup.py can not clean these files; they have to be removed "manually".
#	(Before 2012Aug09 they were removed by "rm ../*.so" in shmem/Makefile.generic.)
	##rm *.so # will fail if no .so files
	for x in *.so ; do if [ -r $$x.o ] ; then rm $$x.o ; fi ; done
