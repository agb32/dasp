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
