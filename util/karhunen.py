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
"""For generating a karhunen-loeve modes...
I think this is correct... but haven't tested to see...
Signs may be different?
"""
import util.FITS,scipy.linalg,numpy,util.zernikeMod,string

def KLmodes(nmodes,zern=None,npup=None,nzern=None,noll=None,dtype=numpy.float32,retvt=0):
    """nmodes is the number of modes to specify...
    noll is None, or a noll matrix, with piston term removed.
    npup is the number of pixels across the pupil if zern is not specified
    zern is a zernike array (with piston removed) or None (in which case it will be created).
    nzern is the number of zernike modes to use - taken from zern if zern specified.
    """
    if zern==None:
        if npup==None:
            raise Exception("util.karhunen.KLmodes: need to specify npup if zern not specified")
        if nmodes==None:
            raise Exception("util.karhunon.KLmodes: need to specify nmodes if zern not specified")
        if nzern==None:
            nzern=nmodes
        zern=util.zernikeMod.Zernike(npup,nzern+1,computeInv=0).zern[1:]
    else:
        npup=zern.shape[1]
        nzern=zern.shape[0]

    kl=numpy.zeros((nmodes,npup,npup),dtype)
    if noll==None:
        noll=util.FITS.Read(string.join(__file__.split("/")[:-2]+["util","noll.fits"],"/"))[1]
    noll=noll[:nzern,:nzern]
    u,w,vt=scipy.linalg.svd(noll)
    #Now, vt[x] gives the zernike coefficients for mode x.
    for i in range(nmodes):
        for j in range(nzern):
            kl[i]+=vt[i,j]*zern[j]
    if retvt:
        return kl,vt
    else:
        return kl
