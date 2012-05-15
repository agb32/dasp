#$Id: serialise.py,v 1.12 2009/09/11 08:17:48 ali Exp $

"""The DFB serialise library, modified by AGB for improved
functionality.  This module allows efficient serialising of large
arrays (and other data) so that they can be send efficiently over a
socket (or stored in a file), rather like pickle, but with greater
efficiency (most arrays are sent direct, without converting to a
string). """

# Python implementation of functions to read and write serialised data
try:
    import Numeric
except:
    print "Not importing Numeric."
    class dummy:
        ArrayType="Numeric not imported"
    Numeric=dummy()
import types, os, cPickle,numpy


Tuple = 'TUPLE_TYPE'
Char = 'CHAR_TYPE'
translate = [#currently no more that 128 (was 16) different types allowed (was 8).
    types.StringType, 
    "b",#numpy.int8,
    "h",#numpy.int16,
    "i",#numpy.int32,
    None,
    "f",#numpy.float32,
    "d",#numpy.float64,
    types.TupleType,
    types.DictType,#added by agb
    types.IntType,#added by agb
    types.FloatType,#added by agb
    types.ListType,#added by agb
    "Pickled",#added by agb
    'l',#added by agb
    'D',#added by agb
    'L',#added by agb
    'F',#added by agb
    'I',#added by agb
    'H',#added by agb
]

def Deserialise(string, start = None, end = None):
    """Deserialise data obtained back into its original form
    @param string: The string to be deserialised
    @type string: String
    @param start: Starting point
    @type start: Int or None
    @param end: Ending point
    @type end: Int or None
    @return: The data
    @rtype: List
    """
    if (start == None): start = 0
    if (end == None): end = len(string)
    self = []
    while(start < end):
        typ, endian, length = DecodeHeader(string[start:])
        start = start + 5
        if (typ == types.TupleType):
            data = tuple(Deserialise(string, start, start + length))
        elif typ==types.ListType:
            data=Deserialise(string,start,start+length)
        elif (typ==types.DictType):#added by agb
            d = Deserialise(string,start, start+length)
            #list of form [key,data,key,data,...]
            data={}
            while len(d)>0:
                key=d.pop(0)
                val=d.pop(0)
                if type(key)==numpy.ndarray:#key must be int or string.
                    key=key[0]
                data[key]=val
        elif (typ == types.StringType):
            if (string[start+length-1] == "\0"):
                data = string[start:start + length - 1]
            else:
                data = string[start:start + length]
        elif typ==types.FloatType:
            data=numpy.fromstring(string[start:start+length],numpy.float64)
            if endian!=numpy.little_endian:
                data=float(data.byteswap()[0])
            else:
                data=float(data[0])
        elif typ==types.IntType:
            data=numpy.fromstring(string[start:start+length],numpy.int32)
            if endian!=numpy.little_endian:
                data=int(data.byteswap()[0])
            else:
                data=int(data[0])
        elif typ==None:
            data=None
        elif typ=="Pickled":#agb
            if string[start+length-1]=="\0":
                data=string[start:start+length-1]
            else:
                data=string[start:start+length]
            data=cPickle.loads(data)
        else:
            #print "deserialise numpy"
            shapelen=numpy.fromstring(string[start:start+4],numpy.int32)
            if endian!=numpy.little_endian:
                shapelen=shapelen.byteswap()
            shapelen=shapelen[0]*4
            shape=numpy.fromstring(string[start+4:start+4+shapelen],numpy.int32)
            if endian!=numpy.little_endian:
                shape=tuple(shape.byteswap())
            else:
                shape=tuple(shape)
            data =numpy.fromstring(string[start+4+shapelen:start+length],typ)
            if endian != numpy.little_endian:
                data = data.byteswap()
            data.shape=shape#=Numeric.reshape(data,shape)
        self.append(data)
        start = start + length
    return(self)


def ReadMessage(infile):
    """Read the next serialised entity (most likely a tuple) from a file/socket
    and return it deserialised.
    Returns None if there is no data but will raise an exception if there
    is a finite but insufficient amount of data
    @param infile: File with a read method.
    @type infile: Object
    @return: The data
    @rtype: List
    """
    header = os.read(infile, 5)
    headerSize = len(header)
    if headerSize == 0:
        return None
    if headerSize != 5:
        raise IOError, "Serial packet header only %d bytes long" \
              % headerSize 
    type, endian, length = DecodeHeader(header)
    readlength = 0
    body = ''
    while readlength < length:
        body = body + os.read(infile, length-readlength)
        readlength = len(body)
        #print "serialise readlength:",readlength
    # Strip off extra tupling because we know it's a singleton
    return Deserialise(header+body)[0]

def DecodeHeader(string):
    """Decode the header of the string
    @param string: Serialised string
    @type string: String
    @return: Type, endian and length
    @rtype: Tuple
    """
    byte = ord(string[0])
    length = (  (ord(string[1]) << 24)
              | (ord(string[2]) << 16)
              | (ord(string[3]) << 8)
              | (ord(string[4])))
    endian = byte & 1
    type = (byte >> 1) & 0x7f#agb - this causes problems...
    type = translate[type]
    return(type, endian, length)

def Serialise(value):
    """Serialise a value
    @param value: The data to be serialised
    @type value: User defined
    @return: Serialised value
    @rtype: String
    """
    thisType = type(value)
    if thisType == types.IntType:
        value = numpy.array(value, numpy.int32)#"i"
    elif thisType == types.FloatType:
        value = numpy.array(value, numpy.float64)#"d"
    elif thisType==Numeric.ArrayType:#THIS SHOULD SAY NUMERIC not numpy.
        value=numpy.array(value,copy=0)
        thisType=type(value)
    #thisType = type(value)
    
    if thisType == numpy.ndarray:#Numeric.ArrayType
        #print "serialise numpy1",value.dtype
        try:
            headerByte=translate.index(value.dtype.char)#typecode())
        except:#see at start of this file for allowed data types...
            print "Datatype %s of array not known in serialise - conversion will fail."%str(value.dtype.char)
            raise
        stringValue=numpy.array(len(value.shape),numpy.int32).tostring()
        stringValue+=numpy.array(value.shape,numpy.int32).tostring()
        stringValue+= value.tostring()
    elif thisType==types.IntType:#added by agb
        headerByte=translate.index(types.IntType)
        stringValue=value.tostring()
    elif thisType==types.FloatType:#added by agb
        headerByte=translate.index(types.FloatType)
        stringValue=value.tostring()
    elif thisType == types.StringType:
        headerByte=translate.index(thisType)
        stringValue = value + '\0'
    elif thisType == types.TupleType:
        headerByte=translate.index(types.TupleType)
        stringValue = ''
        for element in value:
            stringValue = stringValue + Serialise(element)
    elif thisType == types.ListType:#added by agb
        headerByte=translate.index(types.ListType)
        stringValue = ''
        for element in value:
            stringValue = stringValue + Serialise(element)
    elif thisType == types.DictType:#added by agb
        headerByte=translate.index(types.DictType)
        stringValue=''
        for key in value:
            stringValue = stringValue + Serialise(key)+Serialise(value[key])
    elif thisType==types.NoneType:#added by agb
        headerByte=translate.index(None)
        stringValue=''
    else: #added by agb:
        print "WARNING: serialise pickling object with type %s (could be inefficient)"%type(value)
        stringValue=cPickle.dumps(value)+'\0'
        headerByte=translate.index("Pickled")
    #else:#added by agb
    #    print "SERIALISE FAILED - unrecognised type:",thisType
    #    raise "SERIALISE ERROR"
    headerByte = (headerByte << 1) | numpy.little_endian#Numeric.LittleEndian
    length = numpy.array(len(stringValue), numpy.int32)
    if numpy.little_endian: length=length.byteswap()
    return chr(headerByte)+length.tostring()+stringValue

def SerialiseToList(value,sendList=[]):
    """Puts serialised data into a list which should then be sent over a
    socket.  This reduces the overhead by not copying
    numeric arrays into a string.  This should reduce overhead for large
    numeric arrays.
    @param value:  Value to be serialised
    @type value: User defined
    @param sendList: List to be sent
    @type sendList: List of arrays and strings
    @return: The length of the list
    @rtype: Int
    """
    thisType = type(value)
    if thisType == types.IntType:
        value = numpy.array(value, numpy.int32)#"i"
    elif thisType == types.FloatType:
        value = numpy.array(value, numpy.float64)#"d"
    elif thisType==Numeric.ArrayType:#THIS SHOULD SAY NUMERIC not numpy.
        value=numpy.array(value,copy=0)
        thisType=type(value)

    #thisType = type(value)
    length=None
    hdr=None
    ending=None
    hdrSize=5
    if thisType == numpy.ndarray:
        #print "serialise numpy",value.dtype,translate.index(value.dtype),value.dtype==numpy.float64
        try:
            headerByte=translate.index(value.dtype.char)#code())
        except:
            print "util.serialise - headerByte unknown '%s'"%value.dtype.char
        hdr=numpy.array(len(value.shape),numpy.int32).tostring()
        hdr+=numpy.array(value.shape,numpy.int32).tostring()
        if len(value.shape)==0:
            length=value.itemsize
        else:
            length=reduce(lambda x,y:x*y,value.shape)*value.itemsize
        length+=len(hdr)
        if not value.flags.contiguous:#make contiguous so can be sent...
            value=numpy.array(value)
        #print headerByte,hdr
    elif thisType==types.IntType:#added by agb
        headerByte=translate.index(types.IntType)
        length=value.itemsize
    elif thisType==types.FloatType:#added by agb
        headerByte=translate.index(types.FloatType)
        length=value.itemsize
    elif thisType == types.StringType:
        headerByte=translate.index(thisType)
        ending=['\0']#need to send value+'\0'
        length = len(value)+1
    elif thisType == types.TupleType:
        headerByte=translate.index(types.TupleType)
        ending=[]
        length=0
        for element in value:
            length+=SerialiseToList(element,ending)+hdrSize
    elif thisType == types.ListType:#added by agb
        headerByte=translate.index(types.ListType)
        ending=[]
        length=0
        for element in value:
            length+=SerialiseToList(element,ending)+hdrSize
        
    elif thisType == types.DictType:#added by agb
        headerByte=translate.index(types.DictType)
        ending=[]
        length=0
        for key in value:
            length+=SerialiseToList(key,ending)+hdrSize
            length+=SerialiseToList(value[key],ending)+hdrSize
    elif thisType==types.NoneType:#added by agb
        headerByte=translate.index(None)
        ending=[]
        length=0
    else:#added by agb:
        print "WARNING: serialiseToList pickling object with type %s (could be inefficient)"%type(value)
        value=cPickle.dumps(value)
        headerByte=translate.index("Pickled")
        length=len(value)+1
        ending=['\0']
    #else:#added by agb
    #    print "SERIALISE FAILED - unrecognised type:",thisType
    #    raise Exception("SERIALISE ERROR")
    headerByte = (headerByte << 1) | numpy.little_endian
    lngth=numpy.array(length,numpy.int32)
    if numpy.little_endian: lngth=lngth.byteswap()
    txt=chr(headerByte)+lngth.tostring()
    if hdr!=None:
        txt+=hdr
    if len(sendList)>0 and type(sendList[-1])==types.StringType:
        sendList[-1]+=txt
    else:
        sendList.append(txt)
    if type(value) not in [types.ListType,types.DictType,types.TupleType]:
        if type(value)==types.StringType and len(sendList)>0 and type(sendList[-1])==types.StringType:
            sendList[-1]+=value
        elif type(value)!=types.NoneType:
            sendList.append(value)
    if ending!=None:
        for t in ending:
            if type(t)==types.StringType and len(sendList)>0 and type(sendList[-1])==types.StringType:
                sendList[-1]+=t
            else:
                sendList.append(t)
        #sendList+=ending
    return length

def Send(value,sock,verbose=0):
    """Send a message/value over a socket.
    @param value: Data to be sent
    @type value: User defined
    @param sock: Socket
    @type sock: socket.socket instance
    @param verbose: Flag, whether to print messages
    @type verbose: Int
    """
    l=[]
    #print "serialisetolist"
    SerialiseToList(value,l)
    #print "sertolist"
    #print l,len(l)
    for v in l:
        if verbose==1:
            if type(v)==types.StringType:
                print v
            elif type(v)==numpy.ndarray:#Numeric.ArrayType:
                print "numpy array, shape",v.shape
            else:
                print "Type unknown",type(v)
        if sock!=None:
            try:
                sock.sendall(v)
            except:
                print "Error - sock.sendall failed - couldn't send serialised data"
                raise
        #print "sent...",v
if __name__ == "__main__":
    """Testing function"""
    import sys
    option = sys.argv[1]
    fileName = sys.argv[2]
    if option == 'w':
        file = open(fileName, 'w')
        file.write(Serialise(["python test first part!",
                              Numeric.array([1.2345678,3.45,4.5], Numeric.Float32),
                              numpy.array([[1,2,3],[4,5,6]], numpy.int16),
                              numpy.array(99,numpy.int32),
                              ("inner string", 0x12345), 12345.6789]))
        file.write(Serialise(["python test second part!", "the end"]))
    else:
        file = open(fileName, 'r')
        while 1:
            a = ReadMessage(file.fileno())
            if a is None: break
            print a
