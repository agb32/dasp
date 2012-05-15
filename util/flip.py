#$Id: flip.py,v 1.10 2008/03/05 15:43:09 ali Exp $
"""Utility functions for flipping an array"""

### Flip array quadrants for FFTs ####################################


#import Numeric
import numpy as na


def fliparray(mx):
	"""Flip the data within each quadrant, but keep quatrants in
	same place, e.g.:
	
 	 >>> a
	 array([[   1,    2,    3,    4],
	        [  11,   22,   33,   44],
		[ 111,  222,  333,  444],
		[1111, 2222, 3333, 4444]],'i')
	 >>> flip.fliparray(a).astype("i")
	 array([[  22,   11,   44,   33],
	        [   2,    1,    4,    3],
		[2222, 1111, 4444, 3333],
		[ 222,  111,  444,  333]],'i')
	       
        Bad version of flip, incompatible with FFT coordinate definition for non-WFS applications (ask FA for further details) - however, is still correct for WFS.
	@param mx: The array to be flipped
	@type mx: Array
	@return: Flipped array
	@rtype: Array
       """
	n =mx.shape[0]
	n2=mx.shape[1]
	if(n!=n2):
		print 'Fliparray:  Array must be square.'
		return mx

	mx0=na.zeros((n,n),mx.dtype)

	mx0[:n/2,:n/2]        	 = mx[-n/2-1:-n-1:-1,-n/2-1:-n-1:-1]
	mx0[n/2:n,n/2:n]	 = mx[-1:-n/2-1:-1,-1:-n/2-1:-1]
	mx0[:n/2,n/2:n]	 	 = mx[-n/2-1:-n-1:-1,-1:-n/2-1:-1]
	mx0[n/2:n,:n/2]      	 = mx[-1:-n/2-1:-1,-n/2-1:-n-1:-1]

	return mx0

def fliparray2(mx,output=None):
	"""Flip the quadrants, but keep quatrants in
	same place, e.g.:
	
	 >>> a
	 array([[   1,    2,    3,    4],
	        [  11,   22,   33,   44],
		[ 111,  222,  333,  444],
		[1111, 2222, 3333, 4444]],'i')
	 >>> flip.fliparray2(a).astype("i")
	 array([[ 333,  444,  111,  222],
	        [3333, 4444, 1111, 2222],
		[   3,    4,    1,    2],
		[  33,   44,   11,   22]],'i')
		
	Correct version of flip for non-WFS usage, compatible with FFT coordinate definition
	This is effectively a rotation of 180 degrees from the result of
	fliparray().
	Note, fliparray2(a)[::-1,::-1]===fliparray(a)
	So, not really sure this does what it should!!!  Its just mirrored in x and y...
	@param mx: The array to be flipped
	@type mx: Array
	@return: Flipped array
	@rtype: Array
       """
	n =mx.shape[0]
	n2=mx.shape[1]
	if(n!=n2):
		print 'Fliparray:  Array must be square.'
		return mx
	if type(output)==type(None):
		mx0=na.zeros((n,n),mx.dtype)
	else:
		mx0=output
	n=n/2
	mx0[:n,:n]        	 = mx[n:,n:,]
	mx0[n:,n:,]        	 = mx[:n,:n]
	mx0[n:,:n]        	 = mx[:n,n:,]
	mx0[:n,n:,]        	 = mx[n:,:n]
	
	return mx0


