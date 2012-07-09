#
# This function checks the input and decides, whether the multiplication is done by 
#
#        cmod.utils.dot
# or by 
#        numpy.dot.
#
# Currently (July 2nd, 2012) cmod.utils.dot does only the matrix-vector multiplication. When more
# functionality is added to cmod.utils.dot (e.g. matrix-matrix multiplication), this must be
# changed here to be included in aosim.
#

# Possible usages:
# x = dot(a, b)
# x = dot(a, b, c)
# x = dot(a, b, nthr)
# x = dot(a, b, c, nthr)
#
# a, b ... input arrays
# r    ... 1-D array - result vector (if you provide it on input, the function does 
#                      not need to allocate memory; increase in speed of about 0.3%)
# nthr ... the number of threads:
#           if nthr <= 0 : sysconf(_SC_NPROCESSORS_ONLN)/2 threads for big matrices,
#                          1 thread for small matrices (m*n < 500.000)
#           if nthr >= 1 : nthr number of threads are used

def dot(a,b,c=None,nthr=0):
    """Calls cmod.utils.dot or numpy.dot, depending on the input."""

    ## If it is a matrix-vector multiplication, call the cmod.utils.dot:
    if len(a.shape)==2 and len(b.shape)==1 and a.flags.contiguous and b.flags.contiguous:
        import cmod.utils
        c=cmod.utils.dot(a, b, c, nthr) # Arguments: input matrix, input vector, 
                                        #            output vector, number of threads

    ## otherwise call the numpy.dot:
    else:
        import numpy
        if c==None:
            c=numpy.dot(a,b)
        else:
            c[:]=numpy.dot(a,b)

    return c
