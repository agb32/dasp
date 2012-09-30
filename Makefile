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
	(cd gui/bin && ln -fs ../simctrl/simctrl.py && ln -fs ../paramgui/paramgui.py && ln -fs ../simsetup/simsetup.py)
	chmod +x util/portdict.py
	(cd gui/bin && ln -fs ../../util/portdict.py)
#echo "Compiling"
	@printf 'Compiling\n'
	(cd cmod/ ; make install )
#echo Remember to export PYTHONPATH=${PWD} and add gui/bin to your path.
# color for printf
#31=red, 32=green, 33=yellow,34=blue, 35=pink, 36=cyan, 37=white
	@printf '\033[31mRemember to export PYTHONPATH=%s and add gui/bin to your path.\033[m\n' ${PWD}
clean:
	@printf 'Cleaning c modules\n'
	(cd cmod/ ; make clean )
	rm -f base/__init__.py util/__init__.py science/__init__.py cmod/__init__.py gui/__init.py__ gui/selectionbox/__init__.py gui/myFileSelection/__init__.py gui/textbox/__init__.py gui/simctrl/__init__.py gui/dialog/__init__.py
	find . -name \*.pyc | xargs rm

#	for x in base/__init__.py util/__init__.py science/__init__.py cmod/__init__.py gui/__init.py__ gui/selectionbox/__init__.py cmod/Numfftw3/__init__.py gui/myFileSelection/__init__.py gui/textbox/__init__.py gui/simctrl/__init__.py gui/dialog/__init__.py ; do if [ -a $$x ] ; then rm $$x ; fi ; done
