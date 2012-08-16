# Makefile for aosim
# Ali Bharmal <n.a.bharmal@durham.ac.uk>	13/06/2005
#
# $Id: Makefile,v 1.7 2011/12/06 05:14:39 ali Exp $
#

install:
	echo "Making rww modules"
	(cd src/ ; make install )
	echo "Making shmem module"
	(cd shmem/ ; make install )
	echo "Making svd module"
	(cd svd ; make install)
# commented out by UB, 2012 Aug 8th:
# (fails to build)
#	echo "Making MKL module"
#	(cd mkl ; python setup.py install --install-lib=../)
clean:
	echo "Cleaning rww modules"
	(cd src/ ; make clean )
	echo "Cleaning shmem module" 
	(cd shmem/ ; make clean )
	echo "Cleaning svd module"
	(cd svd; make clean)

