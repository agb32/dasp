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
import math,random,numpy
#Numeric
def mypoissondist(mean,nsamp=1024):
    """Generates a poisson distribution by considering the probability distribution and then randomising the order.  This ensures that even after the sample has been used many times, the overall probability is still correct.  See py/poisson/poisson.py for more details.  This is used for the FPGA wfscent stuff"""
    if mean>32:
        raise Exception("Mean value too large, use normal distribution")
    y=[]
    x=range(64)
    out=[]
    for i in range(64):
        tmp=int(round(calcpprob(i,mean)*nsamp))
        y.append(tmp)
        for j in range(tmp):
            out.append(i)
    #print len(out)
    while len(out)<nsamp:
        out.append(int(round(mean)))
    out=out[:nsamp]
    #now randomise the order...
    d={}
    for i in out:
        k=random.random()
        while d.has_key(k):
            k=random.random()
        d[k]=i
    l=d.keys()
    l.sort()
    out=[]
    for i in l:
        out.append(d[i])
    return numpy.array(out,numpy.uint8)
def fact(x):
    f=1.
    i=1.
    while i<=x:
        f*=i
        i+=1
    return f
def calcpprob(x,mean):
    if mean<=0 or x<0:
        return 0
    return math.exp(-mean)*(float(mean)**x)/fact(x)
