import numpy,scipy.special,util.zernikeMod
def nollCoeff(i,j):
    """
    Generate coefficients for the Noll matrix.
    
    Seems to do a reasonable job...

    i,j are termed with piston==1

    Factor, which is approx 7.2e-3 given by Roddier was obtained from TB mathemetica code.
    """
    n,m,s=util.zernikeMod.nm(i,1)
    np,mp,sp=util.zernikeMod.nm(j,1)
    if s!=sp:
        return 0.
    if mp!=m:
        return 0.#only terms with same azimuthal freq are correlated.
    #arr=numpy.zeros((n,n),numpy.float64)
    pow=((n+np-2*m)/2.)
    sign=(-1)**pow
    gamma=scipy.special.gamma
    gamfn=gamma(14./3)*gamma((n+np-5./3)/2)/(gamma((n-np+17./3)/2)*gamma((np-n+17./3)/2)*gamma((n+np+23./3)/2))
    factor=(2*numpy.pi)**(11./3)*0.0457654*(1./2)**(5./3)/((2*numpy.pi)**(14./3))*numpy.pi
    arr=factor*numpy.sqrt((n+1)*(np+1))*sign*numpy.pi**(8./3)*gamfn
    return arr

def noll(n,dr0=1.):
    """noll matrix size n (including piston)
    If needed, this could easily be updated to return a sparse matrix format... (ask AGB).
    Generate a noll matrix of size n.  This is equation 3.14 in
    Roddiers book AO in astronomy.
    dr0 is d/r0 for the telescope

    """
    arr=numpy.zeros((n,n),numpy.float64)
    for i in range(1,n):#ignore piston
        for j in range(1,n):#ignore piston
            arr[i,j]=nollCoeff(i+1,j+1)
    return arr*dr0**(5./3)
