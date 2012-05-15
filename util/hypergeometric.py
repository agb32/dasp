def hyp2F1(a,b,c,z,maxiter=1000):
    """Hypergeometric function 2F1
    Compare with scipy.special.hyp2f1
    """
    go=1
    t=1.
    res=1.
    k=0
    while go:
        tnew=t*(a+k)*(b+k)*z/((c+k)*(k+1))
        res+=tnew
        t=tnew
        #print k,tnew,res
        k+=1
        if tnew*1e12<res:
            go=0
        if k==maxiter:
            print "util.hypergeometric - hyp2F1 failed to converge after %d iterations"%maxiter
            go=0
            res=None
    return res

def hyp2F3(a,b,c,d,e,z,maxiter=1000):
    """Hypergeometric function 2F3
    Used by Winter for zernike phase variances for von-Karman phase screens.
    Taken from functions.wolfram.com and (to work out what the definitions
    mean!) www.efunda.com/math/hypergeometric 
    """
    go=1
    t=1
    res=1.
    k=0
    while go:
        tnew=t*(a+k)*(b+k)*(c+k)*z/((d+k)*(e+k)*(k+1))
        res+=tnew
        t=tnew
        #print k,tnew,res
        k+=1
        if tnew*1e12<res:
            go=0
        if k==maxiter:
            print "util.hypergeometric - hyp2f3 failed to converge after %d iterations"%maxiter
            go=0
            res=None
    return res
