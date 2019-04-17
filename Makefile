# Makefile for aosim
# Ali Bharmal <n.a.bharmal@durham.ac.uk>	13/06/2005
#
# $Id: Makefile,v 1.21 2006/11/07 11:17:18 ali Exp $
#

install:
	@printf 'Touching __init__.py files.\n'
	touch base/__init__.py util/__init__.py science/__init__.py cmod/__init__.py gui/__init__.py gui/selectionbox/__init__.py gui/myFileSelection/__init__.py gui/textbox/__init__.py gui/simctrl/__init__.py gui/dialog/__init__.py gui/pylab/__init__.py
	chmod +x gui/simctrl/simctrl.py
	chmod +x gui/simsetup/simsetup.py
	chmod +x gui/paramgui/paramgui.py
	(cd gui/bin && ln -fs ../simctrl/simctrl.py && ln -fs ../paramgui/paramgui.py && ln -fs ../simsetup/simsetup.py && ln -fs ../../util/analyse.py && ln -fs ../../util/daspgrep.py && ln -fs ../simctrl/simctrl.py daspctrl.py && ln -fs ../simsetup/simsetup.py daspsetup.py && ln -fs ../../util/analyse.py daspanalyse.py)
	chmod +x util/portdict.py
	(cd gui/bin && ln -fs ../../util/portdict.py)
#echo "Compiling"
	@printf 'Compiling\n'
	(cd cmod/ ; make install )
#echo Remember to export PYTHONPATH=${PWD} and add gui/bin to your path.
# color for printf
#31=red, 32=green, 33=yellow,34=blue, 35=pink, 36=cyan, 37=white
	@printf '\033[31mRemember to:\nexport PYTHONPATH=$$PYTHONPATH:%s\nexport PATH=$$PATH:%s/gui/bin\033[m\n' ${PWD} ${PWD}
clean:
	@printf 'Cleaning c modules\n'
	(cd cmod/ ; make clean )
	rm -f base/__init__.py util/__init__.py science/__init__.py cmod/__init__.py gui/__init.py__ gui/selectionbox/__init__.py gui/myFileSelection/__init__.py gui/textbox/__init__.py gui/simctrl/__init__.py gui/dialog/__init__.py
	for fn in $( find . -name \*.pyc ) ; do if [ -a ${fn} ] ; then rm ${fn} ; fi ; done

#	for x in base/__init__.py util/__init__.py science/__init__.py cmod/__init__.py gui/__init.py__ gui/selectionbox/__init__.py cmod/Numfftw3/__init__.py gui/myFileSelection/__init__.py gui/textbox/__init__.py gui/simctrl/__init__.py gui/dialog/__init__.py ; do if [ -a $$x ] ; then rm $$x ; fi ; done

ubuntu1404:
	sudo apt-get install python-dev fftw3-dev libatlas-dev gsl-bin libgsl0-dev libatlas-base-dev python-scipy nfs-common screen glade python-glade2 python-matplotlib python-mpi4py
	cat INSTALL.txt
ubuntu1604:
	sudo apt-get install python-dev fftw3-dev libatlas-dev gsl-bin libgsl0-dev libatlas-base-dev python-scipy nfs-common screen glade python-glade2 python-matplotlib python-mpi4py libopenmpi-dev libgsl-dev libatlas3-base
	cat INSTALL.txt
	echo "Please run sudo ln -s /usr/lib/libatlas.so.3 /usr/lib/libatlas.so"

#If don't want to use the git repository for this, can try: wget http://github.com/xianyi/OpenBLAS/tarball/v0.2.13
openblas:
	(cd ~/ && git clone https://github.com/xianyi/OpenBLAS.git)
	(cd ~/OpenBLAS && make)
	sudo mkdir /opt/openblas
	(cd ~/OpenBLAS && sudo make PREFIX=/opt/openblas install)
	echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/openblas/lib
