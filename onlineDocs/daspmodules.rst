DASP modules
============

DASP modules are the key building blocks of a simulation.  These can be
connected together graphically using daspsetup.py, automatically using
daspbuilder.py, or manually in a text editor.  daspbuilder.py creates
an xml file which can be loaded into daspsetup.py.  Both
daspbuilder.py and daspsetup.py create a python simulation file, which
can then be executed.

Key simulation modules
======================

Phase screens
-------------
Generate translating atmospheric phase screens with Von Karman
statistics.

.. automodule:: science.iscrn
		:members:
		:undoc-members:
		:show-inheritance:

Pupil phase
-----------

Generate the phase at the telescope pupil for a given direction and a
given wavelength.

.. automodule:: science.iatmos
		:members:
		:undoc-members:
		:show-inheritance:


Zonal DM
--------

A DM surface for a given direction (if not ground conjugate) and given
wavelength.  Represented by an actuator map with interpolation
(various options) between the actuators.

.. automodule:: science.xinterp_dm
		:members:
		:undoc-members:
		:show-inheritance:


Modal DM
--------

A DM comprised of Zernike modes.  Often used for a tip-tilt mirror.

.. automodule:: science.zdm
		:members:
		:undoc-members:
		:show-inheritance:


SHS WFS
-------

A SHS WFS.  Input is phase, output are slope measurements.  In
between, full Fourier propagation, noise, centroiding etc.

.. automodule:: science.wfscent
		:members:
		:undoc-members:
		:show-inheritance:


Pyramid WFS
-----------

A Pyramid WFS.  Input is phase, output are slopes.  Modulation and
noise sources are included.

.. automodule:: science.pyramid
		:members:
		:undoc-members:
		:show-inheritance:


Wide field SHS WFS (solar)
--------------------------

A SHS WFS suitable for solar AO modelling.  Anisoplanatic effects are
evident.  Input is the phase screens and DM surfaces, output are the slopes.


.. automodule:: science.widefield
		:members:
		:undoc-members:
		:show-inheritance:


Reconstruction
--------------

A wavefront reconstruction module, suitable for SCAO and tomographic
reconstruction.  Inputs are the wavefront sensor measurements, outputs
are the DM values.

.. automodule:: science.tomoRecon
		:members:
		:undoc-members:
		:show-inheritance:


Science PSF generation
----------------------

A science PSF module.  Input is the phase at the pupil, at a specified
wavelength (and for a specified direction).  A Fourier propagation is
then performed to give a noiseless PSF.  Analysis to get Strehl,
ensquared energy, FWHM, etc are then performed.

This module can also perform lucky imaging.

.. automodule:: science.science
		:members:
		:undoc-members:
		:show-inheritance:

Physical (Fresnel) propagation
------------------------------

A module for physical propagation of light through the atmosphere to
the telescope pupil.

.. automodule:: science.physProp
		:members:
		:undoc-members:
		:show-inheritance:


Other simulation modules
========================

Not widely used, and may or may not work.
