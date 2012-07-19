try:
    from scipy.sparse import linalg
except:
    print "Unable to import scipy.sparse.linalg in regularisation.py - continuing - but may raise an error if this is used..."
import numpy
#import util.dot as quick


class lsreg:
    def __init__(self):
        pass
    
    def regmatrix1(self,pokemx,reg_value):
        """Add a regularization term reg_value^2 to the least squares formulation
        for iterative solvers. Returns the matrix A for the Ax=b form.
        A = H^TH + a^2I, b= H^Ts, H=poke matrix, a=reg_value, s= centroids"""
        pokemxt=numpy.transpose(pokemx)
        regmx=reg_value**2*numpy.identity(pokemx.shape[0])
        A=quick.dot(pokemx,pokemxt)+regmx
        return A
    
    def solvecg(self,A,b,x0=None,maxiter=40):
        """Solve Ax=b using scipy conjugate gradient routine"""
        x=numpy.zeros(A.shape[0])
        x,info=linalg.cg(A,b,x0=x0,maxiter=maxiter)
        return x

    def wafflemx(self,pokemx,reg_value):
        """This is done directly from theory (D.Gavel 2002) but wasn't tested with simulated waffle so far
        Add a waffle filter (convolution) matrix with a regularization term to create iterative solver operator A=H^tH + a^2F^tF"""
        pokemxt=numpy.transpose(pokemx)
        nact=pokemx.shape[0]
                     
        #First create the subaperture  convolution filter in one dimension
        h0=numpy.concatenate(([1.,-1],numpy.zeros((nact-2))),axis=1)#convolution starts at upper left subaperture, h0 indicates first row waffle points along x axis
        h1=-h0 #like before, the second row
# For each non-zero filter sequence (h0,h1)create 1d convolution matrix
        h0mat=h1mat=numpy.zeros((nact,nact))
        for i in range(nact):
            h0mat[:,i]=numpy.roll(h0,i)
            h1mat[:,i]=numpy.roll(h1,i)

        #Now for the 2D convolution matrix
        #Build a block-column of 2d conv. matrix F and roll as before
        Frow=numpy.concatenate([h0mat,h1mat,numpy.zeros((nact**2-2*nact,nact))],axis=0)
        F=numpy.zeros((nact**2,nact**2))
        for i in range(0,nact**2,nact):
            F[:,i:i+nact]=numpy.roll(Frow,i,axis=0)

        # Waffle matrix V=F^tF
        V=quick.dot(F.transpose(),F)
        #According to the theory this only penalizes local waffle, so piston penalty should be added as in following lines. suspended temporarily to be tested manually
        #pistmx=numpy.ones((pokemx.shape[0],pokemx.shape[0]))
        #regmx=wafmx + reg_value**2*numpy.identity(pokemx.shape[0])
        A=quick.dot(pokemx,pokemxt) +V # reg_value**2*wafmx #+ pistmx
        return A,V
    
def invert(pmx,regval,large=0):
    #assumes pmx.shape=nacts,ncents
    """For large matrices use util.computeRecon"""
    if pmx.shape[0]>pmx.shape[1]:
        pmx=pmx.T
    if large:#do it the larger inversion way
        pmxTpmx=quick.dot(pmx.T,pmx)
        s=pmxTpmx.shape
        pmxTpmx.shape=pmxTpmx.size,
        pmxTpmx[::s[0]+1]+=regval
        pmxTpmx.shape=s
        ipmxTpmx=numpy.linalg.inv(pmxTpmx)
        rmx=quick.dot(ipmxTpmx,pmx.T)
    else:#do it the smaller inversion way
        pmxpmxT=quick.dot(pmx,pmx.T)
        s=pmxpmxT.shape
        pmxpmxT.shape=pmxpmxT.size,
        pmxpmxT[::s[0]+1]+=regval
        pmxpmxT.shape=s
        ipmxpmxT=numpy.linalg.inv(pmxpmxT)
        rmx=quick.dot(pmx.T,ipmxpmxT)
    return rmx
