"""Fit a 3D parabola to an image.

See:
Least Squares Fitting of Data
David Eberly
Geometric Tools, LLC
http://www.geometrictools.com/
Copyright c 1998-2012. All Rights Reserved.
Created: July 15, 1999
Last Modified: February 9, 2008
http://www.geometrictools.com/Documentation/LeastSquaresFitting.pdf
"""
import numpy


def makeFitMatrix(size):
    #First create the matrix
    mx=numpy.zeros((6,6),"f")
    vec=numpy.zeros((6,),"f")
    if type(size) in [type([]),type(())]:
        sizey=size[0]
        sizex=size[1]
    else:
        sizey=sizex=size
    for y in range(sizey):
        for x in range(sizex):
            mx[0,0]+=x**4
            mx[0,1]+=x**3*y
            mx[0,2]+=x**2*y**2
            mx[0,3]+=x**3
            mx[0,4]+=x**2*y
            mx[0,5]+=x**2
            mx[1,0]+=x**3*y
            mx[1,1]+=x**2*y**2
            mx[1,2]+=x*y**3
            mx[1,3]+=x**2*y
            mx[1,4]+=x*y**2
            mx[1,5]+=x*y
            mx[2,0]+=x**2*y**2
            mx[2,1]+=x*y**3
            mx[2,2]+=y**4
            mx[2,3]+=x*y**2
            mx[2,4]+=y**3
            mx[2,5]+=y**2
            mx[3,0]+=x**3
            mx[3,1]+=x**2*y
            mx[3,2]+=x*y**2
            mx[3,3]+=x**2
            mx[3,4]+=x*y
            mx[3,5]+=x
            mx[4,0]+=x**2*y
            mx[4,1]+=x*y**2
            mx[4,2]+=y**3
            mx[4,3]+=x*y
            mx[4,4]+=y**2
            mx[4,5]+=y
            mx[5,0]+=x**2
            mx[5,1]+=x*y
            mx[5,2]+=y**2
            mx[5,3]+=x
            mx[5,4]+=y
            mx[5,5]+=1
    # Invert the matrix
    imx=numpy.linalg.inv(mx)
    return mx,imx
def makeFitVector(data):
    # Now the vector
    vec=numpy.zeros((6,),"f")
    for y in range(data.shape[0]):
        for x in range(data.shape[1]):
            vec[0]+=data[y,x]*x**2
            vec[1]+=data[y,x]*x*y
            vec[2]+=data[y,x]*y**2
            vec[3]+=data[y,x]*x
            vec[4]+=data[y,x]*y
            vec[5]+=data[y,x]
    return vec


def doFit(data):
    imx=makeFitMatrix(data.shape)[1]
    vec=makeFitVector(data)
    # Find the coefficients:
    p=numpy.dot(imx,vec)
    return p#[a,b,c,d,e,f] with a*x^2+b*x*y+c*y^2+d*x+e*y+f

"""
d/dx = 2ax+by+d
d/dy = 2cy+bx+e
d2/dxdy=b

Solve d/dx=0, d/dy=0.
y=(-e-bx)/2c

2ax-be/2c-b^2x/2c+d=0
x(2a-b^2/2c)=be/2c-d
x=(be/2c-d)/(2a-b^2/2c)

"""
def getMaxPos(data):
    p=doFit(data)
    xmax=(p[1]*p[4]/(2*p[2])-p[3])/(2.*p[0]-p[1]*p[1]/(2.*p[2]))
    ymax=-(p[4]-p[1]*xmax)/(2.*p[2])
    return xmax,ymax

def getMaxPosQuadInterp(data):
    if data.shape!=(3,3):
        raise Exception("Should be 3x3, is %s"%str(data.shape))
    a2=(data[2,1]-data[0,1])/2.
    a3=(data[2,1]-2*data[1,1]+data[0,1])/2.
    a4=(data[1,2]-data[1,0])/2.
    a5=(data[1,2]-2*data[1,1]+data[1,0])/2.
    a6=(data[2,2]-data[0,2]-data[2,0]+data[0,0])/4.
    y=(2*a2*a5-a4*a6)/(a6*a6-4*a3*a5)
    x=(2*a3*a4-a2*a6)/(a6*a6-4*a3*a5)
    return x+1,y+1
    
def fitPeak(data,size=3,minimum=1,quadInterp=0):
    """Finds the parabola max/min for a region of size x size around the minimum or maximum of dat."""
    if minimum:
        pos=numpy.argmin(data)
    else:
        pos=numpy.argmax(data)
    x=pos%data.shape[1]
    y=pos//data.shape[1]
    fx=x-size//2
    fy=y-size//2
    if fx<0:
        fx=0
    if fx>data.shape[1]-size:
        fx=data.shape[1]-size
    if fy<0:
        fy=0
    if fy>data.shape[0]-size:
        fy=data.shape[0]-size
    tx=fx+size
    ty=fy+size
    #print fx,tx,fy,ty
    data=data[fy:ty,fx:tx]
    if quadInterp==1 and size==3:
        xp,yp=getMaxPosQuadInterp(data)
    else:
        xp,yp=getMaxPos(data)
    return xp+fx,yp+fy




def getMaxPosOld(data):
    p=doFit(data)
    xmax = ( 2*p[2]*p[3]/p[1]**2 - p[4]/p[1] ) / (1-4*p[2]*p[0]/p[1]**2)
    ymax = (-p[3] - 2*p[0]*xmax)/p[1]
    return xmax,ymax


def fitGaussian(data,retP=0):
    """Fit a 2D gaussian to the data.
    Gaussian form is:
    g . exp( - (x-a)^2/(2sx^2) - (y-b)^2/(2sy^2) )


    """
    ldata=numpy.where(data<=0,-90.,numpy.log(data))
    p=doFit(ldata)
    a=-p[3]/2./p[0]
    b=-p[4]/2./p[2]
    sx2=-1/(2*p[0])
    sy2=-1/(2*p[2])

    g=numpy.exp(p[5]+a*a/(2*sx2) + b*b/(2*sy2))
    if retP:
        return g,a,b,sx2,sy2,p
    else:
        return g,a,b,sx2,sy2

def makeGaussian(g,a,b,sx2,sy2,size=10):
    x=numpy.arange(size)
    y=numpy.arange(size)
    data=g*numpy.exp(-((x-a)**2)/(2*sx2) - ((y[None,:].T-b)**2)/(2*sy2) )
    return data

def makeParabola(cx,cy,size=10):
    x=numpy.arange(size)
    y=numpy.arange(size)
    data=(x-cx)**2+(y[None,:]-cy).T**2
    return data

def test(data=None):
    if data==None:
        size=10
        x=numpy.arange(size)
        y=numpy.arange(size)
        data=(x-3.)**2+(y[None,:]-5.).T**2

# Now make the reconstructed image.
    p=doFit(data)
    print p
    x=numpy.zeros(data.shape,"f")
    y=x.copy()
    x[:]=numpy.arange(data.shape[1])
    y.T[:]=numpy.arange(data.shape[0])

    fitted=p[0]*x**2+p[1]*x*y+p[2]*y**2+p[3]*x+p[4]*y+p[5]
    xmax = ( 2*p[2]*p[3]/p[1]**2 - p[4]/p[1] ) / (1-4*p[2]*p[0]/p[1]**2)

    ymax = (-p[3] - 2*p[0]*xmax)/p[1]

    print xmax,ymax
    return fitted,data
#Now find x, y at the maximum.
#This occurs when d/dx and d/dy ==0.
#d/dx = 2*p[0]*x + p[1]*y + p[3]
#d/dy = 2*p[2]*y + p[1]*x + p[4]

#Solve for x,y:
"""
y = (-p[3] - 2*p[0]*x)/p[1]

x = (-p[4] - 2*p[2]*y)/p[1]

x = (-p[4] - 2*p[2]*(-p[3] - 2*p[0]*x)/p[1])/p[1]

  = -p[4]/p[1] + 2*p[2]*p[3]/p[1]**2 + 4*p[2]*p[0]*x/p[1]**2

x * (1-4*p[2]*p[0]/p[1]**2) = -p[4]/p[1] + 2*p[2]*p[3]/p[1]**2

So, we can now get the position of the max of the parabola:
"""
if __name__=="__main__":
    fitted,data=test()
    import pylab
    pylab.ion()
    pylab.subplot(131)
    pylab.imshow(fitted,interpolation="nearest")
    pylab.subplot(132)
    pylab.imshow(data,interpolation="nearest")
    pylab.subplot(133)
    pylab.imshow(data-fitted,interpolation="nearest")
    raw_input()
