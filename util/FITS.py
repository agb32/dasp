#!/usr/bin/env python
# $Id: FITS.py,v 1.33 2011/01/24 08:00:40 ali Exp $
#
# Functions to read and write FITS image files
import string
import numpy
##import Numeric ## commented out by UB 2012May31


#from Numeric import *
import os.path,os
error = 'FITS error'
#
# Read a FITS image file
#
# Returns a header, image tuple
# The header is returned in two forms, a raw list of header records, and
# a parsed dictionary for keyword access.
#
# By default, the image is returned as floats, using the BSCALE and BZERO
# keywords where available. By setting the asFloat keyword to zero it is 
# returned unscaled and in the (presumably more compact) numeric format
# that it was stored in
# 
def Read(filename, asFloat = 1,savespace=1,doByteSwap=1,compliant=1,memmap=None) :
    """if savespace is set, the array will maintain its type if asfloat is set.
    If doByteSwap is not set, no byteswap will be done if little endian - if this is the case, the file is not actually fits compliant
    if memmap is used, the file is memmapped rather than read.  If this is used, it is recommended that memmap="r" and that the file is non-byte-swapped.
    """
    file = open(filename, "r")
    done=0
    filelen=os.path.getsize(filename)
    returnVal=[]
    while done==0:
        rawHeader = []
        header = {}
        buffer = file.read(2880)
        if buffer[:6] != 'SIMPLE' and buffer[:8]!="XTENSION":
            if compliant:
                print buffer[:80]
                raise Exception(error+ 'Not a simple fits file')
            else:
                print "WARNING - non compliant FITS file"
                return returnVal
        while(1) :
            for char in range(0,2880,80) :
                line = buffer[char:char+80]
                rawHeader.append(line)
                key = string.strip(line[:8]).strip("\0")
                if key :
                    val = line[9:]
                    val = string.strip(val).strip("\0")
                    if val :
                        if val[0] == "'" :
                            try:
                                pos = string.index(val,"'",1)
                                val = val[1:pos]
                            except:
                                val=val[1:]
                        else :
                            pos = string.find(val, '/')
                            if pos != -1 :
                                val = val[:pos]
                    header[key] = val
            if header.has_key('END') : break
            buffer = file.read(2880)
        naxis = string.atoi(header['NAXIS'])
        shape = []
        for i in range(1,naxis+1) :
            shape.append(string.atoi(header['NAXIS%d' % i]))
        shape.reverse()
        numPix = 1
        for i in shape :
            numPix = numPix * i
        bitpix = string.atoi(header['BITPIX'])
        if bitpix == 8 :
            typ = numpy.uint8
        elif bitpix == 16 :
            typ = numpy.int16
        elif bitpix == 32 :
            typ = numpy.int32
        elif bitpix == -32 :
            typ = numpy.float32
            bitpix = 32
        elif bitpix == -64 :
            typ = numpy.float64
            bitpix = 64
        if memmap==None:
            data=numpy.fromfile(file,typ,count=numPix)
        else:
            nel=reduce(lambda x,y:x*y,shape)
            data=numpy.memmap(filename,dtype=typ,mode=memmap,offset=file.tell())[:nel]
            file.seek(bitpix/8*nel,1)
        numByte = numPix * bitpix/8
        #data = file.read(numByte)
        #data = numpy.fromstring(data, dtype=typ)
        #data.savespace(1)
        data.shape = shape
        if numpy.little_endian and doByteSwap:
            if header.has_key("UNORDERD") and header["UNORDERD"]=='T':
                pass
            else:
                data.byteswap(True)
        if asFloat :
            bscale = string.atof(header.get('BSCALE', '1.0'))
            bzero = string.atof(header.get('BZERO', '0.0'))
            if savespace:
                if bscale!=1:
                    data*=bscale#array(bscale,typecode=typ)
                if bzero!=0:
                    data+=bzero#array(bzero,typecode=typ)
            else:
                data = data*bscale + bzero
        returnVal.append({ 'raw' : rawHeader, 'parsed' : header})
        returnVal.append(data)
        ntoread=2880-numByte%2880
        if ntoread!=0 and ntoread!=2880:
            file.read(ntoread)
        #print "Read 1 hdu at %d/%d"%(file.tell(),filelen)
        if file.tell()==filelen:
            done=1
    return returnVal#( { 'raw' : rawHeader, 'parsed' : header},  data  )

#
# Save a numeric array to a FITS image file
#
# The FITS image dimensions are taken from the array dimensions.
# The data is saved in the numeric format of the array, or as 32-bit
# floats if a non-FITS data type is passed.
#
# A list of extra header records can be passed in extraHeader. The
# SIMPLE, BITPIX, NAXIS* and END records in this list are ignored.
# Header records are padded to 80 characters where necessary.
#
def Write(data, filename, extraHeader = None,writeMode='w',doByteSwap=1,preserveData=1) :
    """Writes data to filename, with extraHeader (string or list of strings).  If writeMode="a" will overwrite existing file, or if "a", will append to it.  If doByteSwap==1, then will do the byteswap on systems that require it, to preserve the FITS standard.  If preserveData==1, will then byteswap back after saving to preserve the data.
    """

## OLD (before 2012May31; UB):
#    if type(data)==Numeric.ArrayType:
#        typ = data.typecode()
#    else:#assume numpy...
#        typ=data.dtype.char
#    if typ == Numeric.UnsignedInt8 : bitpix = 8
#    elif typ==Numeric.Int16 : bitpix = 16
#    elif typ==Numeric.Int32 : bitpix = 32
#    elif typ==Numeric.Float32: bitpix=-32
#    elif typ==Numeric.Float64: bitpix=-64

## NEW (UB, 2012May31, 6 lines):
    typ=data.dtype.char
    if   typ=='b': bitpix =  8 # UnsignedInt8
    elif typ=='s': bitpix = 16 # Int16
    elif typ=='i': bitpix = 32 # Int32
    elif typ=='f': bitpix =-32 # Float32
    elif typ=='d': bitpix =-64 # Float64

    #elif typ=="H": bitpix=16
    elif typ=="h": bitpix=16
    else :
        print "FITS: Converting type %s to float32"%typ
##	data = data.astype(Numeric.Float32) # OLD, UB 2012May31
	data = data.astype('f') # Float32   # NEW, UB 2012May31
	bitpix = -32
    shape = list(data.shape)
    shape.reverse()
    naxis = len(shape)
    if writeMode=='w' or os.path.exists(filename)==False:
        header = [ 'SIMPLE  = T']
        writeMode='w'
    else:
        header=  [ "XTENSION= 'IMAGE'"]
    header+=['BITPIX  = %d' % bitpix, 'NAXIS   = %d' % naxis ]
    for i in range(naxis) :
	keyword = 'NAXIS%d' % (i + 1)
	keyword = string.ljust(keyword, 8)
	header.append('%s= %d' % (keyword, shape[i]))
    header.append('EXTEND  = T')
    if doByteSwap==0 and numpy.little_endian:
        header.append('UNORDERD= T')
    if extraHeader != None :
        if type(extraHeader)==type(""):
            extraHeader=[extraHeader]
	for rec in extraHeader :
	    try :
		key = string.split(rec)[0]
	    except IndexError :
		pass
	    else :
		if key != 'SIMPLE' and \
		   key != 'BITPIX' and \
		   key[:5] != 'NAXIS' and \
		   key != 'END' :
		    header.append(rec)
    header.append('END')
    header = map(lambda x: string.ljust(x,80)[:80], header)
    header = string.join(header,'')
    numBlock = (len(header) + 2880 - 1) / 2880
    header = string.ljust(string.join(header,''), numBlock*2880)
    file = open(filename, writeMode)
    file.write(header)
    if numpy.little_endian and doByteSwap:
##        if type(data)==Numeric.ArrayType: # ?? OLD, UB 2012May31
        if False:                           # ?? NEW, UB 2012May31
            data = data.byteswapped()
        else:
            data.byteswap(True)
    data.tofile(file)
    #data = data.tostring()
    #file.write(data)
    numBlock = (data.itemsize*data.size + 2880 - 1) / 2880
    padding = ' ' * (numBlock*2880 - data.itemsize*data.size)
    file.write(padding)
    if numpy.little_endian and doByteSwap and preserveData==1:
##        if type(data)==Numeric.ArrayType: # ?? OLD, UB 2012May31
        if False:                           # ?? NEW, UB 2012May31
            pass
        else:#preserve the data...
            data.byteswap(True)


def ReadHeader(filename, asFloat = 1) :
    file = open(filename, "r")
    header = {}
    rawHeader = []
    buffer = file.read(2880)
    if buffer[:6] != 'SIMPLE' :
	raise Exception(error+ 'Not a simple fits file')
    while(1) :
	for char in range(0,2880,80) :
	    line = buffer[char:char+80]
	    rawHeader.append(line)
	    key = string.strip(line[:8])
	    if key :
		val = line[9:]
		val = string.strip(val)
		if val :
		    if val[0] == "'" :
			pos = string.index(val,"'",1)
			val = val[1:pos]
		    else :
			pos = string.find(val, '/')
			if pos != -1 :
			    val = val[:pos]
		header[key] = val
	if header.has_key('END') : break
	buffer = file.read(2880)
    return( { 'raw' : rawHeader, 'parsed' : header} )

def MakeHeader(shape,dtype,extraHeader=None,doByteSwap=1,extension=0):
    """Return a text string which can be used as a header"""
    if dtype=="b": bitpix=8
    elif dtype=="B": bitpix=8
    elif dtype=="s": bitpix=16
    elif dtype=="h": bitpix=16
    elif dtype=="i": bitpix=32
    elif dtype=="f": bitpix=-32
    elif dtype=="d": bitpix=-64
    else: raise Exception("Unknown datatype in MakeHeader")
    shape=list(shape)
    shape.reverse()
    naxis=len(shape)
    if extension:
        header=["XTENSION= 'IMAGE'"]
    else:
        header=["SIMPLE  = T"]
    header+=["BITPIX  = %d"%bitpix,"NAXIS   = %d"%naxis]
    for i in range(naxis):
	keyword = 'NAXIS%d' % (i + 1)
	keyword = string.ljust(keyword, 8)
	header.append('%s= %d' % (keyword, shape[i]))
    header.append('EXTEND  = T')
    if doByteSwap==0 and numpy.little_endian:
        header.append('UNORDERD= T')
    if extraHeader != None :
        if type(extraHeader)==type(""):
            extraHeader=[extraHeader]
	for rec in extraHeader :
	    try :
		key = string.split(rec)[0]
	    except IndexError :
		pass
	    else :
		if key != 'SIMPLE' and \
		   key != 'BITPIX' and \
		   key[:5] != 'NAXIS' and \
		   key != 'END' :
		    header.append(rec)
    header.append('END')
    header = map(lambda x: string.ljust(x,80)[:80], header)
    header = string.join(header,'')
    numBlock = (len(header) + 2880 - 1) / 2880
    header = string.ljust(string.join(header,''), numBlock*2880)
    return header


def WriteKey(file,key,value=None,comment=None):
    """Will write to file object file at current position."""
    if file.tell()%80!=0:
        print "ERROR: util.FITS.WriteKey - not at the start of a header line"
        raise Exception("FITS error")
    txt=key+" "*(8-len(key))
    if value!=None:
        txt+="= "+str(value)
    if comment!=None:
        txt+=" /"+comment
    txt=txt+" "*(80-len(txt))
    txt=txt[:80]
    file.write(txt)
def EndHeader(file):
    WriteKey(file,"END")
    pos=file.tell()
    pmod=pos%2880
    if pmod>0:
        file.write(" "*(2880-pmod))
def WriteComment(file,comment):
    if file.tell()%80!=0:
        print "ERROR: util.FITS.WriteComment - not at the start of a header line"
        raise Exception("FITS error")
    while len(comment)>0:
        c="COMMENT "+comment[:72]
        comment=comment[72:]
        c=c+" "*(80-len(c))
        file.write(c)

def savecsc(csc,filename,hdr=None,doByteSwap=1):
    """Save a scipy.sparse.csc sparse matrix"""
    if type(hdr)==type(""):
        hdr=[hdr]
    elif type(hdr)==type(None):
        hdr=[]
    hdr.append("SHAPE   = %s"%str(csc.shape))
    rowind=csc.rowind[:csc.indptr[-1]].view(numpy.int32)
    indptr=csc.indptr.view(numpy.int32)
    
    Write(csc.data[:csc.indptr[-1]],filename,extraHeader=hdr,doByteSwap=doByteSwap)
    Write(rowind,filename,writeMode="a",doByteSwap=doByteSwap)
    Write(indptr,filename,writeMode="a",doByteSwap=doByteSwap)
def loadcsc(filename,doByteSwap=1):
    """Load a scipy.sparse.csc sparse matrix"""
    f=Read(filename,savespace=1,doByteSwap=doByteSwap)
    if len(f)==2:
        print "WARNING - loadcsc - %s is not a sparse matrix"%filename
        mx=f[1].astype("f")
    elif len(f)>2:
        import scipy.sparse
        mx=scipy.sparse.csc_matrix((numpy.array(f[1],copy=0),numpy.array(f[3],copy=0).view(numpy.uint32),numpy.array(f[5],copy=0).view(numpy.uint32)),eval(f[0]["parsed"]["SHAPE"]))
    return mx

def saveSparse(sp,filename,hdr=None,doByteSwap=1):
    """Save a sparse matrix - csc or csr."""
    if type(hdr)==type(""):
        hdr=[hdr]
    elif type(hdr)==type(None):
        hdr=[]
    hdr.append("SHAPE   = %s"%str(sp.shape))
    hdr.append("FORMAT  = '%s'"%sp.format)
    Write(sp.data[:sp.indptr[-1]],filename,extraHeader=hdr,doByteSwap=doByteSwap)
    if sp.format=="csr":
        ind=sp.colind[:sp.indptr[-1]]
    elif sp.format=="csc":
        ind=sp.rowind[:sp.indptr[-1]]
    else:
        raise Exception("Sparse matrix type not yet implemented")
    Write(ind.view(numpy.int32),filename,writeMode="a",doByteSwap=doByteSwap)
    Write(sp.indptr.view(numpy.int32),filename,writeMode="a",doByteSwap=doByteSwap)

def loadSparse(filename,doByteSwap=1):
    """load a scipy.sparse matrix - csc or csr."""
    f=Read(filename,savespace=1,doByteSwap=doByteSwap)
    if len(f)==2:
        print "WARNING - loadSparse - %s is not a sparse matrix"%filename
        mx=f[1]
    elif len(f)==6 and len(f[1].shape)==1:
        import scipy.sparse
        if f[0]["parsed"].has_key("FORMAT"):
            fmt=f[0]["parsed"]["FORMAT"]
        else:
            print "Warning - loadSparse %s - assuming csc"%filename
            fmt="csc"
        print "Format is %s"%fmt
        shape=eval(f[0]["parsed"]["SHAPE"])
        indptr=f[5].view(numpy.uint32)
        if fmt=="csc":
            if indptr.shape[0]==(shape[1]+1)*2:
                indptr=indptr.view(numpy.uint64)
            mx=scipy.sparse.csc_matrix((f[1],f[3].view(numpy.uint32),indptr),eval(f[0]["parsed"]["SHAPE"]))
        elif fmt=="csr":
            if indptr.shape[0]==(shape[0]+1)*2:
                indptr=indptr.view(numpy.uint64)
            mx=scipy.sparse.csr_matrix((f[1],f[3].view(numpy.uint32),indptr),eval(f[0]["parsed"]["SHAPE"]))
        else:
            raise Exception("sparse matrix style not known %s"%filename)
    else:
        mx=f[1]
        print "util.FITS.loadSparse - matrix style not known %s, assuming dense matrix"%filename
    return mx

def loadBlockMatrix(filename,doByteSwap=1):
    f=Read(filename,savespace=1,doByteSwap=doByteSwap)
    #if len(f)==2:
    #    mx=f[1]
    #else:
    import util.blockMatrix
    mx=util.blockMatrix.BlockMatrix()
    mx.assign(f[1::2])
    return mx
def saveBlockMatrix(mx,filename,extraHeader=None,doByteSwap=1):
    wm="w"
    for b in mx.blockList:#append the blocks to the file.
        Write(b,filename,writeMode=wm,extraHeader=extraHeader,doByteSwap=doByteSwap)
        wm="a"
    
