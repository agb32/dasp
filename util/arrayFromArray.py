def arrayFromArray(arr,shape,dtype):
    a=arr.view(dtype)
    a.shape=(reduce(lambda x,y:x*y,a.shape),)
    a=a[:reduce(lambda x,y:x*y,shape)]
    a.shape=shape
    return a
