#dasp, the Durham Adaptive optics Simulation Platform.
#Copyright (C) 2004-2016 Alastair Basden and Durham University.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as
#published by the Free Software Foundation, either version 3 of the
#License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#mpirun -np 1 -hostlist n1-c437  /usr/local/bin/mpipython $PWD/thisfile.py
#Python code created using the simulation setup GUI...
#Order of execution may not be quite optimal - you can always change by hand
#for large simulations - typically, the order of sends and gets may not be
#quite right.  Anyway, enjoy...
import util.Ctrl
import science.infScrn
import science.infAtmos
import science.xinterp_dm
import science.wfscent
import science.tomoRecon
import science.science
ctrl=util.Ctrl.Ctrl(globals=globals())
if not "nopoke" in ctrl.userArgList:ctrl.doInitialPokeThenRun()
#ctrl.doInitialPokeThenRun()
scrn=science.infScrn.infScrn(None,ctrl.config,idstr="L0")
atmos=science.infAtmos.infAtmos({"L0":scrn,},ctrl.config,idstr="a")
dm=science.xinterp_dm.dm(None,ctrl.config,idstr="dma")
wfs=science.wfscent.wfscent(dm,ctrl.config,idstr="a")
recon=science.tomoRecon.recon({"a":wfs,},ctrl.config,idstr="ngs")
#sci=science.science.science(dm,ctrl.config,idstr="m")
dm.newParent({"1":atmos,"2":recon,},"dma")
import cProfile, pstats, StringIO
pr = cProfile.Profile()
pr.enable()
ctrl.mainloop([scrn,atmos,dm,wfs,recon])
pr.disable()
s = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print s.getvalue()
print "Simulation finished..."


