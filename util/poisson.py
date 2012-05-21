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
