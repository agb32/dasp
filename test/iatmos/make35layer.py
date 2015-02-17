import numpy
def makeESO35layer(profile=0,dmListIndices=[0,5,8,11,14,17,21,23,28]):
    """profile=0 for median, 1-4 for each quartile.
    Email from Phil, on 15/1/2014.

    """
    atmosData=numpy.array([
1,3.00E+01,5.50E+00,2.42E+01,2.26E+01,2.51E+01,2.55E+01,2.36E+01,
2,9.00E+01,5.50E+00,1.20E+01,1.12E+01,1.16E+01,1.19E+01,1.31E+01,
3,1.50E+02,5.10E+00,9.69E+00,1.01E+01,9.56E+00,9.32E+00,9.81E+00,
4,2.00E+02,5.50E+00,5.90E+00,6.40E+00,5.84E+00,5.57E+00,5.77E+00,
5,2.45E+02,5.60E+00,4.73E+00,4.15E+00,3.70E+00,4.50E+00,6.58E+00,
6,3.00E+02,5.70E+00,4.73E+00,4.15E+00,3.70E+00,4.50E+00,6.58E+00,
7,3.90E+02,5.80E+00,4.73E+00,4.15E+00,3.70E+00,4.50E+00,6.58E+00,
8,6.00E+02,6.00E+00,4.73E+00,4.15E+00,3.70E+00,4.50E+00,6.58E+00,
9,1.13E+03,6.50E+00,3.99E+00,3.12E+00,3.25E+00,4.19E+00,5.41E+00,
10,1.88E+03,7.00E+00,3.24E+00,2.26E+00,3.47E+00,4.04E+00,3.20E+00,
11,2.63E+03,7.50E+00,1.62E+00,1.13E+00,1.74E+00,2.02E+00,1.60E+00,
12,3.50E+03,8.50E+00,2.61E+00,2.21E+00,3.00E+00,3.04E+00,2.18E+00,
13,4.50E+03,9.50E+00,1.56E+00,1.33E+00,1.80E+00,1.82E+00,1.31E+00,
14,5.50E+03,1.15E+01,1.04E+00,8.83E-01,1.20E+00,1.21E+00,8.70E-01,
15,6.50E+03,1.75E+01,1.00E+00,1.47E+00,1.30E+00,8.56E-01,3.74E-01,
16,7.50E+03,2.30E+01,1.20E+00,1.77E+00,1.56E+00,1.03E+00,4.49E-01,
17,8.50E+03,2.60E+01,4.01E-01,5.90E-01,5.21E-01,3.42E-01,1.50E-01,
18,9.50E+03,2.90E+01,1.40E+00,2.06E+00,1.82E+00,1.20E+00,5.24E-01,
19,1.05E+04,3.20E+01,1.30E+00,1.92E+00,1.69E+00,1.11E+00,4.87E-01,
20,1.15E+04,2.70E+01,7.01E-01,1.03E+00,9.12E-01,5.99E-01,2.62E-01,
21,1.25E+04,2.20E+01,1.60E+00,2.32E+00,1.87E+00,1.43E+00,7.96E-01,
22,1.35E+04,1.45E+01,2.59E+00,3.75E+00,3.03E+00,2.31E+00,1.29E+00,
23,1.45E+04,9.50E+00,1.91E+00,2.76E+00,2.23E+00,1.70E+00,9.49E-01,
24,1.55E+04,6.30E+00,9.87E-01,1.43E+00,1.15E+00,8.79E-01,4.90E-01,
25,1.65E+04,5.50E+00,6.17E-01,8.92E-01,7.20E-01,5.50E-01,3.06E-01,
26,1.75E+04,6.00E+00,4.01E-01,5.80E-01,4.68E-01,3.57E-01,1.99E-01,
27,1.85E+04,6.50E+00,2.47E-01,3.57E-01,2.88E-01,2.20E-01,1.22E-01,
28,1.95E+04,7.00E+00,2.16E-01,3.12E-01,2.52E-01,1.92E-01,1.07E-01,
29,2.05E+04,7.50E+00,1.85E-01,2.68E-01,2.16E-01,1.65E-01,9.19E-02,
30,2.15E+04,8.00E+00,1.36E-01,1.96E-01,1.58E-01,1.21E-01,6.74E-02,
31,2.25E+04,8.50E+00,1.11E-01,1.61E-01,1.30E-01,9.89E-02,5.51E-02,
32,2.35E+04,9.00E+00,6.17E-02,8.92E-02,7.20E-02,5.50E-02,3.06E-02,
33,2.45E+04,9.50E+00,9.26E-02,1.34E-01,1.08E-01,8.24E-02,4.59E-02,
34,2.55E+04,1.00E+01,4.94E-02,7.13E-02,5.76E-02,4.40E-02,2.45E-02,
35,2.65E+04,1.00E+01,4.32E-02,6.24E-02,5.04E-02,3.85E-02,2.14E-02]).astype(numpy.float32)
    atmosData.shape=35,8
    atmosData[:,3:]/=atmosData[:,3:].sum(0)
    heightScale=[1.0,0.88,0.92,1.02,1.30]
    atmosData[:,1]*=heightScale[profile]
    strList=atmosData[:,3+profile]
    hList=atmosData[:,1]
    vList=atmosData[:,2]
    dirList=(numpy.arange(atmosData.shape[0])%2)*180.
    hListdm=atmosData[dmListIndices,1].copy()
    hListdm[0]=0
    strListdm=atmosData[dmListIndices,3+profile].copy()
    strListdm/=strListdm.sum()
    r0=[0.157,0.234,0.178,0.139,0.097][profile]
    l0=25.
    #r0=0.129/numpy.cos(30*numpy.pi/180.)**0.6
    return hList,strList,vList,hListdm,strListdm,l0,r0,dmListIndices

def makeESO35layer6dm(profile=0,dmListIndices=[0,8,11,14,17,21]):
    return makeESO35layer(profile,dmListIndices)

def makeESO35layer9dm(profile=0,dmListIndices=[0,5,8,11,14,17,21,23,28]):
    return makeESO35layer(profile,dmListIndices)

def makeESO35layer9dmLower(profile=0,dmListIndices=[0,5,8,11,14,17,20,22,26]):
    return makeESO35layer(profile,dmListIndices)

def makeESO35layer11dm(profile=0,dmListIndices=[0,5,8,11,14,17,19,21,23,28,32]):
    return makeESO35layer(profile,dmListIndices)

def makeESO35layer12dm(profile=0,dmListIndices=[0,5,8,11,13,15,17,19,21,23,28,32]):
    return makeESO35layer(profile,dmListIndices)

def makeESO35layer12dmLower(profile=0,dmListIndices=[0,5,8,11,13,15,17,18,20,22,24,28]):
    return makeESO35layer(profile,dmListIndices)

def makeESO35layerNdmEqual(profile=0,heights=numpy.arange(0,20001,2000),doubleGnd=0):
    """Integrates between DM heights to get strengths.
    If doubleGnd, then doubles the ground layer DM strength - possibly this should be used because this one only integrates half the turbulence.
    """
    hList,strList,vList,hListdm,strListdm,l0,r0,dmListIndices=makeESO35layer(profile=profile)
    endHeight=0
    h=0
    cn=numpy.zeros((len(heights),),numpy.float32)
    for i in range(len(heights)):
        startHeight=endHeight
        if i==len(heights)-1:#last one - include the rest of turb
            endHeight=hList[-1]
            if endHeight<startHeight:
                endHeight=startHeight
        else:#the mean of the 2 dm heights
            endHeight=(heights[i]+heights[i+1])/2.
        #Now integrate the cn2 between these heights.
        if h==0:#special first case - incase hList[0]>0
            cn[i]=hList[0]*strList[0]
        elif h>0 and h<len(hList):#integrate from startHeight upto hList[h]
            sh=(startHeight+hList[h])/2.
            s=strList[h-1]+(sh-hList[h-1])/(hList[h]-hList[h-1])*(strList[h]-strList[h-1])
            cn[i]+=s*(hList[h]-startHeight)
            print startHeight,hList[h],s,cn[i]
        while h<len(hList) and hList[h]<endHeight:
            if h+1<len(hList):
                if hList[h]>=startHeight:
                    if hList[h+1]<=endHeight: #Now integrate between the two layers
                        cn[i]+=(hList[h+1]-hList[h])*(strList[h]+strList[h+1])/2
                        print hList[h],hList[h+1],(strList[h]+strList[h+1])/2,cn[i]
                    else:#integrate up to endHeight
                        #Mean pos between h1 and endHeight:
                        eh=(endHeight+hList[h])/2.
                        s=strList[h]+(eh-hList[h])/(hList[h+1]-hList[h])*(strList[h+1]-strList[h])
                        cn[i]+=s*(endHeight-hList[h])
                        print hList[h],endHeight,s,cn[i]
                else:
                    pass
            else:#got to last height... assume no turb after this height.
                break
            h+=1
    if doubleGnd:
        cn[0]*=2
    cn/=cn.sum()
    return hList,strList,vList,heights,cn,l0,r0,dmListIndices

def make9layer():
    strList=numpy.array([0.5224,0.0260,0.0444,0.116,0.0989,0.0295,0.0598,0.0430,0.06])
    hList=numpy.array([47.,140.,281.,562.,1125.,2250.,4500.,9000.,18000.])
    vList=numpy.array([4.55,  12.61,  12.61,   8.73,   8.73,  14.55,  24.25,  38.8 , 20.37])
    r0=0.135
    l0=30
    return hList,strList,vList,hList,strList,l0,r0,numpy.arange(9)

def plotit(flist=[makeESO35layer,
                  ],
           subplot=None,
           vscale=1000.,
           fmt=".",
           outfile=None,
           plotline=1,
           plotdm=1,logx=0,logdata=0,plotvelocity=1,
           labels=None,colours=None,plotbar=0,xlab=None,ylab=None
           ):
    """To plot different options - eg different profile, dm heights etc, can use e.g. flist=[lambda:make35layer.makeESO35layer(profile=2)]
    """
    import pylab
    pylab.clf()
    if xlab!=None:
        pylab.xlabel(xlab)
    if ylab!=None:
        pylab.ylabel(ylab)
        
    for i in range(len(flist)):
        fun=flist[i]
        if logx:
            pylab.xscale("log")
        if subplot!=None:
            pylab.subplot(subplot[0],subplot[1],i+1)
        hnew,snew,vnew,hdm,sdm=fun()[:5]
        vnew2=numpy.zeros((2,vnew.size),vnew.dtype)
        vnew2[1]=vnew
        if logdata:
            snew=numpy.log10(snew)
            vnew2[1]=numpy.log10(vnew2[1])

        if labels==None:
            l=None
        else:
            l=labels[i]
        if colours!=None:
            colour=colours[i]
        else:
            colour="green"
        if plotvelocity:
            pylab.errorbar(snew,hnew,xerr=vnew2/vscale,fmt=fmt)
        else:
            pylab.plot(snew,hnew,fmt,color=colour)
        if plotdm:
            pylab.plot(sdm/sdm[0]*snew[0],hdm,'r+',markersize=20,color="red",ls="None")
        if plotline:
            pylab.plot(snew,hnew,color=colour,label=l)
            l=None
        if plotbar:
            print hnew,snew
            left=0
            if logx:
                left=1e-06
            pylab.barh(numpy.array(hnew),numpy.array(snew),height=200.,label=l,color=colour,edgecolor=colour,log=logx,alpha=0.3,left=left)
            #for j in range(len(snew)):
            #    pylab.bar(left=1e-06,height=hnew[j]+100,width=snew[j],bottom=hnew[j]-100,label=l,color=colour)
        if labels!=None:
            pylab.legend(loc=1)
        #pylab.xlim([1e-06,0.6])
    if outfile==None:
        pylab.show()
    else:
        pylab.savefig(outfile)



#make35layer.plotit([make35layer.makeESO35layer12dmLower,make35layer.make9layer,make40layer.makeESO40layer9dm],labels=["35 layer","9 layer","40 layer"],plotvelocity=0,colours=["red","green","blue"],plotline=1,plotbar=0,logdata=0,xlab="Cn2/fraction",ylab="Height / m",logx=1,outfile="out2.png")
#make35layer.plotit([make35layer.makeESO35layer12dmLower,make35layer.make9layer,make40layer.makeESO40layer9dm],labels=["35 layer","9 layer","40 layer"],plotvelocity=0,colours=["red","green","blue"],plotline=1,plotbar=1,logdata=0,xlab="Cn2/fraction",ylab="Height / m",logx=0,outfile="out1.png")
#make35layer.plotit([make35layer.makeESO35layer12dmLower,make35layer.make9layer,make40layer.makeESO40layer9dm],labels=["35 layer","9 layer","40 layer"],plotvelocity=0,colours=["red","green","blue"],plotline=1,plotbar=1,logdata=0,xlab="Cn2/fraction",ylab="Height / m",logx=1,outfile="logbarline.png")
#make35layer.plotit([make35layer.makeESO35layer12dmLower,make35layer.make9layer,make40layer.makeESO40layer9dm],labels=["35 layer","9 layer","40 layer"],plotvelocity=0,colours=["red","green","blue"],plotline=0,plotbar=1,logdata=0,xlab="Cn2/fraction",ylab="Height / m",logx=1,outfile="logbar.png")
