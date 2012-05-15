/*
 * shmem.c
 *
 * Alistair Basden <a.g.basden@durham.ac.uk>    13/06/2005
 *
 * $Id: shmem.c,v 1.23 2009/10/23 14:18:48 ali Exp $
 */

#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>
#include <semaphore.h>
//#include <linux/sem.h>
//#include <sys/stat.h>
#define OPLEN 1
#if defined(__GNU_LIBRARY__) && !defined(_SEM_SEMUN_UNDEFINED) || defined(__APPLE__) 
     /* union semun is defined by including <sys/sem.h> */
#else
     /* according to X/OPEN we have to define it ourselves */
   union semun {
       int val;                    /* value for SETVAL */
       struct semid_ds *buf;       /* buffer for IPC_STAT, IPC_SET */
       unsigned short int *array;  /* array for GETALL, SETALL */
       struct seminfo *__buf;      /* buffer for IPC_INFO */
       };
#endif
static PyObject *ShmemError;

//SEMMSL - max number of semaphores per semaphore set (250) - in linux/sem.h
//SEMMNI - max number of semaphore sets (128) - in linux/sem.h
//SEMMNS - max number of semaphores (250*128) - in linux/sem.h
#define SEMMNI 128
static int semidArray[SEMMNI];
static int allFinished;//this is used to determine whether EIVAL and EIDRM errors are treated as error or not.

typedef struct filenamelist{
    char *fname;
    void *next;
} fnamelist;

static fnamelist *fnlist=NULL;

int saveSemID(int semid){//returns -1 if error.
    int i,done;
    done=0;
    for(i=0; (i<SEMMNI && done==0); i++){
	if(semidArray[i]==-1){
	    semidArray[i]=semid;
	    done=1;
	}
    }
    return done-1;
}
int deleteSemID(int semid){
    int i,done;
    done=0;
    for(i=0; (i<SEMMNI && done==0); i++){
	if(semidArray[i]==semid){
	    semidArray[i]=-1;
	    semctl(semid,0,IPC_RMID);//remove semaphore set.
	    done=1;
	}
    }
    return done-1;
}
static PyObject* shmem_allFinished(PyObject *self,PyObject *args){
  int i=1;
  if(!PyArg_ParseTuple(args,"|i",&i)){
    printf("Usage: value for allFinished (optional)\n");
    return NULL;
  }
  allFinished=i;
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* shmem_deleteAllSems(PyObject *self,PyObject *args){
  int i;
  for(i=0; i<SEMMNI; i++){
    if(semidArray[i]!=-1){
      semctl(semidArray[i],0,IPC_RMID);//remove semaphore set.
      semidArray[i]=-1;
    }
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* shmem_cleanUp(PyObject *self,PyObject *args){
  int i;
  fnamelist *fnlistptr;
  allFinished=1;
  for(i=0; i<SEMMNI; i++){
    if(semidArray[i]!=-1){
      semctl(semidArray[i],0,IPC_RMID);//remove semaphore set.
      semidArray[i]=-1;
    }
  }
  //now delete temporary files.
  while(fnlist!=NULL){
      fnlistptr=fnlist;
      fnlist=fnlist->next;
      unlink(fnlistptr->fname);
      free(fnlistptr);
  }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* shmem_deleteSemid(PyObject *self,PyObject *args){
    int semid;
    if(!PyArg_ParseTuple(args,"i",&semid)){
	printf("Usage: semid\n");
	return NULL;
    }
    if(deleteSemID(semid)==-1){
	printf("semid %d not deleted\n",semid);
	return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}
static PyObject *shmem_getCurrentSemids(PyObject *self,PyObject *args){
  npy_intp dims[1];
  dims[0]=SEMMNI;
  return PyArray_SimpleNewFromData(1,dims,NPY_INT,(void*)semidArray);
}



static PyObject* shmem_shmemopen(PyObject *self, PyObject *args){//tuple of dimensions,type,name
  char *type,*name;
  int flag;
  npy_intp dims[4];
  int idims[4];
  int fd,i;
  int size,itype;
  char* addr;
  PyObject* dim;
  int printerr=1;
  idims[0]=idims[1]=idims[2]=idims[3]=-1;
  if (!PyArg_ParseTuple(args, "Ossi|i", &dim,&type,&name,&flag,&printerr)){
    printf("Usage1: (sizex,y,z,zz),type,name,createFlag, printerr flag (optional, default 1)\n");
    return NULL;
  }
  if (!PyArg_ParseTuple(dim,"i|iii",&idims[0],&idims[1],&idims[2],&idims[3])){
    printf("Usage2: (sizex,y,z,zz),type,name,createFlag\n");
    return NULL;
  }
  for(i=0; i<4; i++){
    //printf("dim[%d]=%d\n",i,idims[i]);
    dims[i]=idims[i];
  }
  fd=shm_open(name, O_RDWR|(O_CREAT*flag),0777);//O_RDWR | O_CREAT
  if(fd==-1){
    if(printerr)
      printf("shm_open failed for %s:%s\n",name,strerror(errno));
    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
    close(fd);
    return NULL;
  }
  size=dims[0]*(dims[1]>0?dims[1]:1)*(dims[2]>0?dims[2]:1)*(dims[3]>0?dims[3]:1);
  if(*type=='f'){
    size*=sizeof(float);
    itype=NPY_FLOAT;
  }else if(*type=='d'){
    size*=sizeof(double);
    itype=NPY_DOUBLE;
  }else if(*type=='i'){
    size*=sizeof(int);
    itype=NPY_INT;
  }else if(*type=='s'){
    size*=sizeof(short);
    itype=NPY_SHORT;
  }else if(*type=='c'){
    size*=sizeof(char);
    itype=NPY_CHAR;
  }else if(*type=='u'){
    size*=sizeof(char);
    itype=NPY_UBYTE;
  }else if(*type=='b'){
    size*=sizeof(char);
    itype=NPY_BYTE;
  }else{
    printf("Type code not recognised - assuming int\n");
    itype=NPY_INT;
    size*=sizeof(int);
  }
  if(ftruncate(fd,size) == -1 ) {
    printf("ftruncate failed: %s\n",strerror( errno ) );
    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
    close(fd);
    return NULL;
  }
  addr=(char*)mmap(0,size,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
  close(fd);
  if(addr==MAP_FAILED){
    printf("mmap failed: %s\n",strerror( errno ) );
    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
    return NULL;
  }
  //printf( "Map addr is 0x%08x\n", (int)addr );
  //printf("%d %d %d %s %s %d\n",dim1,dim2,dim3,type,name,flag);
  return PyArray_SimpleNewFromData(dims[1]>0?(dims[2]>0?(dims[3]>0?4:3):2):1,dims,itype,addr);
  
}

static PyObject* shmem_shmemunlink(PyObject *self, PyObject *args){//tuple of dimensions,type,name
  char *name;
  if(!PyArg_ParseTuple(args, "s", &name)){
    printf("Usage: name\n");
    return NULL;
  }
  if(shm_unlink(name)){
    printf("unlink failed: %s\n",strerror(errno));
    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;

}
static PyObject* shmem_shmemunmap(PyObject *self, PyObject *args){//tuple of dimensions,type,name
  int size,esize,i;
  PyObject *Narray;
  PyArrayObject *NumericArray;
  if(!PyArg_ParseTuple(args, "O", &Narray)){
    printf("Usage: NumericArray\n");
    return NULL;
  }
  NumericArray=(PyArrayObject*)Narray;
  if(NumericArray->descr->type_num==NPY_CHAR){
    esize=sizeof(char);
  }else if(NumericArray->descr->type_num==NPY_UBYTE){
    esize=sizeof(char);
  }else if(NumericArray->descr->type_num==NPY_BYTE){
    esize=sizeof(char);
  }else if(NumericArray->descr->type_num==NPY_CHAR){
    esize=sizeof(char);
  }else if(NumericArray->descr->type_num==NPY_STRING){
    esize=sizeof(char);
  }else if(NumericArray->descr->type_num==NPY_SHORT){
    esize=sizeof(short);
  }else if(NumericArray->descr->type_num==NPY_INT){
    esize=sizeof(int);
  }else if(NumericArray->descr->type_num==NPY_LONG){
    esize=sizeof(long);
  }else if(NumericArray->descr->type_num==NPY_FLOAT){
    esize=sizeof(float);
  }else if(NumericArray->descr->type_num==NPY_DOUBLE){
    esize=sizeof(double);
  }else if(NumericArray->descr->type_num==NPY_CFLOAT){
    esize=sizeof(float)*2;
  }else if(NumericArray->descr->type_num==NPY_CDOUBLE){
    esize=sizeof(double)*2;
  }else if(NumericArray->descr->type_num==NPY_OBJECT){
    esize=sizeof(PyObject*);
  }else if(NumericArray->descr->type_num==NPY_STRING){
    esize=sizeof(char);
  }else{
    printf("Type not known - assuming int (%d)\n",NumericArray->descr->type_num);
    esize=sizeof(int);
  }
  size=1;
  for(i=0; i<NumericArray->nd; i++){
    size*=NumericArray->dimensions[i];
  }
  if(munmap(NumericArray->data,size*esize)==-1){
    printf("Unmap (%d) failed: %s\n",size*esize,strerror(errno));
    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
  }
  //now create new array... this is a bit dodgy, but not really sure what to do, to make sure that the array is really unmapped...
  //Maybe this is better:
  /*free(NumericArray->strides);
  free(NumericArray->dimensions);
  NumericArray->nd=1;
  NumericArray->dimensions=(npy_intp*)malloc(sizeof(npy_intp));
  NumericArray->strides=(npy_intp*)malloc(sizeof(npy_intp));
  NumericArray->data=(char*)malloc(esize);
  memset(NumericArray->data,0,esize);*/
  for(i=0; i<NumericArray->nd; i++){
    NumericArray->dimensions[0]=0;
    NumericArray->strides[0]=esize;
  }
  NumericArray->data=NULL;
  Py_INCREF(Py_None);
  return Py_None;
  
}

static PyObject* shmem_newsemid(PyObject *self,PyObject *args,PyObject *kwds){
  //create a new semid, which can be used for semaphores etc...
  //creates a temporary file for this to be generated from.
  char *name=NULL;
  int fd;
  char *txt=NULL;
  key_t key;
  int projid=1;
  int semid;
  int nSems=1;//one semaphore in the set by default.
  int existingFile=0;//file doesn't exist.
  fnamelist *fnlistptr=NULL;
  static char *kwlist[]={"name","projid","existingFile","nSems",NULL};
  if(!PyArg_ParseTupleAndKeywords(args,kwds,"|siii",kwlist,&name,&projid,&existingFile,&nSems)){
    printf("Usage: name(optional - temp file will then be /tmp/nameXXXXXX), projid (optional), existingFile(if file already exists), nSems(optional - default 1)\n");
    return NULL;
  }
  if(existingFile==0){
    if(name==NULL){
      asprintf(&txt,"/tmp/aosimXXXXXX");
    }else{
      asprintf(&txt,"/tmp/%sXXXXXX",name);
    }
    if((fd=mkstemp(txt))==-1){
      printf("Failed to create unique tempfilename for semid %s\n",txt);
      free(txt);
      return NULL;
    }//this tmpfile should be automaticall deleted on program exit.
  }else{
    if(name==NULL){
      printf("Failed to create semid from NULL file\n");
      return NULL;
    }else{
      txt=strdup(name);
    }
  }
  //add txt to the file list.  These files can then be deleted when cleanUp is called.  
  fnlistptr=(fnamelist*)malloc(sizeof(fnamelist));
  fnlistptr->fname=txt;
  fnlistptr->next=fnlist;
  fnlist=fnlistptr;
  if((key=ftok(txt,projid))==-1){
    printf("Couldn't get key: %s\n",strerror(errno));
    //free(txt);
    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
    return NULL;
  }
  //free(txt);
  if((semid=semget(key,nSems,0700|IPC_CREAT))==-1){
    printf("semget failed (%d, %d): %s\n",nSems,(int)key,strerror(errno));
    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
    return NULL;
  }else{
      if(saveSemID(semid)==-1)
	  printf("WARNING: call to saveSemID with semid=%d failed\n",semid);
  }
  return Py_BuildValue("i",semid);
}
static PyObject* shmem_initSemaphore(PyObject *self, PyObject *args){
  //use this to put initial values in a semaphore set.
  union semun argument;
  int semid,semNo,semVal;
  if(!PyArg_ParseTuple(args, "iii", &semid,&semNo,&semVal)){
    printf("Usage: semid, semNo, semVal\n");
    return NULL;
  }
  argument.val=semVal;
  if(semctl(semid,semNo,SETVAL,argument)==-1){
    printf("semctl failed: %s\n",strerror(errno));
    PyErr_SetFromErrno(ShmemError);
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* shmem_semop(PyObject *self,PyObject *args){
  //use this to perform a blocking semop...
  int semid,nSemNo,nSemOp,val,i;
  struct sembuf *operations=NULL;
  PyObject *semNoList=NULL,*semOpList=NULL;
  if(!PyArg_ParseTuple(args,"iOO",&semid,&semNoList,&semOpList)){
    printf("Usage: semid, semNoList, semOpList\n");
    return NULL;
  }
  if(semNoList==NULL || semOpList==NULL){
    printf("ERROR: semNoList or semOpList are NULL\n");
    return NULL;
  }
  if(PyInt_Check(semNoList) && PyInt_Check(semOpList)){//is integer not list
    if(!PyInt_Check(semOpList)){
      printf("ERROR: BOTH args must be integers, or lists of same length\n");
      return NULL;
    }
    nSemNo=1;
    operations=malloc(sizeof(struct sembuf)*nSemNo);
    memset(operations,0,sizeof(struct sembuf)*nSemNo);
    operations[0].sem_num=(unsigned short)PyInt_AsLong(semNoList);
    operations[0].sem_op=(short)PyInt_AsLong(semOpList);
  }else if(PyList_Check(semNoList) && PyList_Check(semOpList)){
    nSemNo=PyList_Size(semNoList);
    nSemOp=PyList_Size(semOpList);
    if(nSemNo!=nSemOp || nSemOp==-1){
      printf("Error:  nSemNo!=nSemOp, %d!=%d or both are negative\n",nSemNo,nSemOp);
      return NULL;
    }
    operations=malloc(sizeof(struct sembuf)*nSemNo);
    memset(operations,0,sizeof(struct sembuf)*nSemNo);
    for(i=0; i<nSemNo; i++){
      operations[i].sem_num=(unsigned short)PyInt_AsLong(PyList_GetItem(semNoList,i));
      operations[i].sem_op=(short)PyInt_AsLong(PyList_GetItem(semOpList,i));
    }
  }else{
    printf("ERROR: semNoList is not a list or int\n");
    return NULL;
  }
  Py_BEGIN_ALLOW_THREADS;
  val=semop(semid,operations,nSemNo);
  Py_END_ALLOW_THREADS;
  free(operations);
  if(val==-1){
      if(allFinished==0){
	  if(errno==EAGAIN){
	      printf("timeout received while doing shmem.semop:%s\n",strerror(errno));
	  }else if(errno==EINTR){
	      printf("Signal recved while acquiring sem/waiting:%s\n",strerror(errno));
	  }else if(errno==EINVAL){
	      printf("sem/wait doesn't exist (semid=%d): %s\n",semid,strerror(errno));
	  }else if(errno==EIDRM){
	      printf("sem/wait has been removed from system: %s\n",strerror(errno));
	  }else{
	      printf("Failed to get sem/wait (unexpected error - possibly because errno not thread safe?): %s\n",strerror(errno));
	  }
	  PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
	  return NULL;
      }else{
	  printf("All finished... shmem.semop failed\n");
      }
  }//operations complete...
  Py_INCREF(Py_None);
  return Py_None;
}
 
static PyObject* shmem_semdel(PyObject *self,PyObject *args){
    int semid;
    if(!PyArg_ParseTuple(args,"i",&semid)){
	printf("Usage: semid\n");
	return NULL;
    }
    if(semctl(semid,0,IPC_RMID)==-1){
	printf("Call to semctl IPC_RMID with semid %d failed\n",(int)semid);
	PyErr_SetFromErrno(ShmemError);
	return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}
static PyObject* shmem_getSemValue(PyObject *self, PyObject *args){
  int semid,semNo,semVal;
  if(!PyArg_ParseTuple(args, "ii", &semid,&semNo)){
    printf("Usage: semid, semNo\n");
    return NULL;
  }
  if((semVal=semctl(semid,semNo,GETVAL))==-1){
    printf("semctl GETVAL failed: %s\n",strerror(errno));
    PyErr_SetFromErrno(ShmemError);
    return NULL;
  }
  return Py_BuildValue("i",semVal);
}
// For Darwin: has a bug (feature?) in semop 
// If can get the sem0, wait for sem1 to be 0, then add 1 to sem 2, then wait
// for sem3 to be zero.
static PyObject *shmem_acquireAndWait(PyObject *self, PyObject *args){
  // used to acquire a semaphore with zero timeout, and then block on another
  // semaphore used as an event.  This is done atomically.  This is used e.g.
  // in Splitter.py.  The semaphores should be set up correctly previously.
  int semid,val;
  struct sembuf operations[1];
  int rtval=1;//success...
  if(!PyArg_ParseTuple(args, "i", &semid)){
    printf("Usage: semid.  Warning: behaves differently on Linux and Darwin.\n");
    return NULL;
  }
  operations[0].sem_num=0;
  operations[0].sem_op=-1;//acquire a semaphore (first semaphore)
  operations[0].sem_flg=IPC_NOWAIT;//error if can't acquire immediately.
  //begin/end added 051125 - though not needed because won't block
  //Py_BEGIN_ALLOW_THREADS;//allow multi threading around blocking ops
  val=semop(semid,operations,1);
  //Py_END_ALLOW_THREADS;//allow multi threading around blocking ops
  if(val==-1 && errno!=EAGAIN && allFinished==0){
    //failed to get block (NOT timeout - which is allowed...)
    if(errno==EINTR){
      printf("Signal recved while acquiring sem/waiting:%s\n",strerror(errno));
    }else if(errno==EINVAL){
      printf("sem/wait doesn't exist: %s\n",strerror(errno));
    }else if(errno==EIDRM){
      printf("sem/wait has been removed from system: %s\n",strerror(errno));
    }else{
      printf("Failed to get sem/wait (unexpected error - possibly because errno not thread safe?): %s\n",strerror(errno));
    }
    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
    return NULL;
  }else if(val==-1){//timed out - ie couldn't get first semaphore (if called in Splitter, this means its the last thread to call...
    if(allFinished==0)
      rtval=0;
    else
      rtval=-1;
  }else{//semaphore obtained so now wait on event.
    operations[0].sem_num=1;
    operations[0].sem_op=0;//wait on the second semaphore.
    operations[0].sem_flg=0;
    Py_BEGIN_ALLOW_THREADS;//allow multi threading around blocking ops
    val=semop(semid,operations,1);// if successful, blocks waiting for 
    Py_END_ALLOW_THREADS;         // "event" (semnum 2)
    if(val==-1){
	if(allFinished==0){
	    printf("Failed while waiting on event\n");
	    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
	    return NULL;
	}else{
	    rtval=-1;
	}
    }

#ifdef __APPLE__
    //now add one to sem2.
    operations[0].sem_num=2;
    operations[0].sem_op=1;//add one to the third semaphore.
    operations[0].sem_flg=0;
    //Py_BEGIN_ALLOW_THREADS;//allow multi threading around blocking ops
    val=semop(semid,operations,1);// if successful, blocks waiting for 
    //Py_END_ALLOW_THREADS;         // "event" (semnum 2)
    if(val==-1){
	if(allFinished==0){
	    printf("Failed while increasing sem2\n");
	    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
	    return NULL;
	}else{
	    rtval=-1;
	}
    }
    //now wait for sem3.
    operations[0].sem_num=3;
    operations[0].sem_op=0;//wait on 4th sem.
    operations[0].sem_flg=0;
    Py_BEGIN_ALLOW_THREADS;//allow multi threading around blocking ops
    val=semop(semid,operations,1);// if successful, blocks waiting for 
    Py_END_ALLOW_THREADS;         // "event" (semnum 2)
    if(val==-1){
	if(allFinished==0){
	    printf("Failed while sem3 waiting on event\n");
	    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
	    return NULL;
	}else{
	    rtval=-1;
	}
    }
#endif
  }
  return Py_BuildValue("i",rtval);//return 0 if couldn't get first semaphore.
}
// For Darwin: has a bug (feature?) in semop 
// Set sem3 to 1, Set sem1 to zero, then wait for sem2 to reach n-1.  Then set
// sem1 to 1, and set set3 to zero.
static PyObject* shmem_trigEvent(PyObject *self,PyObject *args){
  // used to set and clear an event on semnum 1 of the set.  Used e.g. in
  // Splitter.py.
  int semid,val,waitno=1;
  struct sembuf operations[2];
  if(!PyArg_ParseTuple(args, "i|i", &semid,&waitno)){
    printf("Usage: semid, waitno (required only for Darwin).  Behaves differently on Linux and Darwin\n");
    return NULL;
  }

#ifdef __APPLE__
  operations[0].sem_num=3; // the "event" semaphore.
  operations[0].sem_op=1;  // set the event.
  operations[0].sem_flg=0; // shouldn't have to block...
  val=semop(semid,operations,1);
  if(val==-1){
      printf("Failed while setting sem3\n");
      PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
      return NULL;
  }

  operations[0].sem_num=1; // the "event1" semaphore.
  operations[0].sem_op=-1; // trigger the event.
  operations[0].sem_flg=IPC_NOWAIT; // shouldn't have to block...
  Py_BEGIN_ALLOW_THREADS;  // allow multi threading around blocking ops
  val=semop(semid,operations,1);
  Py_END_ALLOW_THREADS;    // allow multi threading around blocking ops
  if(val==-1){//failed. 
      printf("set of the event semaphore failed...\n");
      if(errno==EAGAIN){
	  printf("Timeout whilst trying to set/clear event flag: %s\n",strerror(errno));
      }else if(errno==EINTR){
	  printf("Signal received while trying to set/clear event flag: %s\n",strerror(errno));
      }else if(errno==EINVAL){
	  printf("event doesn't exist: %s\n",strerror(errno));
      }else if(errno==EIDRM){
	  printf("event has been removed from system: %s\n",strerror(errno));
      }else{
	  printf("Failed to get/set/clear event: %s\n",strerror(errno));
      }
      PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
      return NULL;
  }
  
  operations[0].sem_num=2;//the "event2" semaphore.
  operations[0].sem_op=-waitno;//wait until all threads have incremented
  operations[0].sem_flg=0;//should have to block...
  Py_BEGIN_ALLOW_THREADS;//allow multi threading around blocking ops
  val=semop(semid,operations,1);
  Py_END_ALLOW_THREADS;//allow multi threading around blocking ops
  if(val==-1){//failed.
      printf("set of the event2 semaphore failed...\n");
      if(errno==EAGAIN){
	  printf("Timeout whilst trying to set/clear event flag: %s\n",strerror(errno));
      }else if(errno==EINTR){
	  printf("Signal received while trying to set/clear event flag: %s\n",strerror(errno));
      }else if(errno==EINVAL){
	  printf("event doesn't exist: %s\n",strerror(errno));
      }else if(errno==EIDRM){
	  printf("event has been removed from system: %s\n",strerror(errno));
      }else{
	  printf("Failed to get/set/clear event: %s\n",strerror(errno));
      }
      PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
      return NULL;
  }
  

  operations[0].sem_num=1;//the "event" semaphore.
  operations[0].sem_op=1;//set the event.
  operations[0].sem_flg=0;//shouldn't have to block...
  val=semop(semid,operations,1);
  if(val==-1){
      printf("Failed while setting sem1 to 1\n");
      PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
      return NULL;
  }
  


  operations[0].sem_num=3;//the "event" semaphore.
  operations[0].sem_op=-1;//set the event.
  operations[0].sem_flg=IPC_NOWAIT;//shouldn't have to block...
  Py_BEGIN_ALLOW_THREADS;//allow multi threading around blocking ops
  val=semop(semid,operations,1);
  Py_END_ALLOW_THREADS;//allow multi threading around blocking ops
  if(val==-1){//failed.
      printf("set of the event3 semaphore failed...\n");
      if(errno==EAGAIN){
	  printf("Timeout whilst trying to set/clear event flag: %s\n",strerror(errno));
      }else if(errno==EINTR){
	  printf("Signal received while trying to set/clear event flag: %s\n",strerror(errno));
      }else if(errno==EINVAL){
	  printf("event doesn't exist: %s\n",strerror(errno));
      }else if(errno==EIDRM){
	  printf("event has been removed from system: %s\n",strerror(errno));
      }else{
	  printf("Failed to get/set/clear event: %s\n",strerror(errno));
      }
      PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
      return NULL;
  }

#else // not Darwin- assume linux...
  int i;

  operations[0].sem_num=1;//the "event" semaphore.
  operations[0].sem_flg=IPC_NOWAIT;//shouldn't have to block...
  for(i=-1; i<2; i+=2){//when i==-1 will set event, when i==+1, will clear it.
    operations[0].sem_op=i;//set then clear the event (things are blocked wait on a zero).
    //operations[1].sem_num=1;
    //operations[1].sem_op=+1;//clear the event.
    //operations[1].sem_flg=0;
    //Begin/end threads probably not needed because we want this to
    //happen atomically.  We know the value will be 1, so okay to decrement.  Non-blocking call anyway...
    //Py_BEGIN_ALLOW_THREADS;//allow multi threading around blocking ops
    val=semop(semid,operations,1);
    //Py_END_ALLOW_THREADS;//allow multi threading around blocking ops
    if(val==-1){//failed.
      printf("%s of the event semaphore failed...\n",i==-1?"Set":"Clear");
      if(errno==EAGAIN){
	printf("Timeout whilst trying to set/clear event flag: %s\n",strerror(errno));
      }else if(errno==EINTR){
	printf("Signal received while trying to set/clear event flag: %s\n",strerror(errno));
      }else if(errno==EINVAL){
	printf("event doesn't exist: %s\n",strerror(errno));
      }else if(errno==EIDRM){
      printf("event has been removed from system: %s\n",strerror(errno));
      }else{
	printf("Failed to get/set/clear event: %s\n",strerror(errno));
      }
      PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
      return NULL;
    }
  }
#endif
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* shmem_newlock(PyObject *self,PyObject *args,PyObject *keywds){
  //Create a new lock which can be used as a simple lock, or as a rw lock (many readers, one writer).  
  char *name;
  int projid=1;
  key_t key;
  int semid,i;
  int realfileflag=0;
  char *txt=NULL;
  int settoone=0;//this should be used for the process creating the lock
  int nReaders=0;
  union semun argument;
  static char *kwlist[] = {"name", "projid", "realfileflag","setone","nReaders",NULL};
  if(!PyArg_ParseTupleAndKeywords(args,keywds, "s|iiii",kwlist, &name,&projid,&realfileflag,&settoone,&nReaders)){
    printf("Usage: name, projid(optional),realfileflag(optional),setone(optional), nReaders (optional)\n");
    return NULL;
  }
  if(realfileflag==0){
    asprintf(&txt,"/dev/shm%s",name);
  }else{
    asprintf(&txt,"%s",name);
  }
  if((key=ftok(txt,projid))==-1){
    printf("Couldn't get key (%s): %s\n",txt,strerror(errno));
    free(txt);
    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
    return NULL;
  }
  free(txt);
  //printf("Key as int: %d\n",key);
  if((semid=semget(key,2+nReaders,0700|IPC_CREAT))==-1){//note - only creates if not already created.
    printf("semget failed (%d, %d): %s\n",2+nReaders,(int)key,strerror(errno));
    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
    return NULL;
  }else{
      if(saveSemID(semid)==-1)
	  printf("WARNING: call to saveSemID with semid=%d failed\n",semid);
  }
  //printf("Semid=%d\n",semid);
  if(settoone==1){
    argument.val = 1;
    if(semctl(semid, 0, SETVAL, argument)==-1){//one writer only.
      printf("semctl failed: %s\n",strerror(errno));//set to 1 owner initially.
      PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
      return NULL;
    }
    argument.val=0;
    if(semctl(semid,1,SETVAL,argument)==-1){//nothing reading yet
      printf("semctl failed: %s\n",strerror(errno));
    }
    for(i=2; i<2+nReaders; i++){//data not yet written, so nReaders
      if(semctl(semid,i,SETVAL,argument)==-1){//bits dont need to be set.
	printf("semctl failed: %s\n",strerror(errno));
	PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
	return NULL;
      }
    }
  }
  //printf("semval: %d %d\n",semctl(semid,0,GETVAL),semctl(semid,1,GETVAL));
  return Py_BuildValue("i", semid);
}

static PyObject* shmem_getlock(PyObject *self, PyObject *args){//tuple of dimensions,type,name
  //grabs a semaphore (write) lock (or waits until it becomes available and nothing has a read lock).
  //This grabs the write lock (the only lock if not bothering with read locks).
  //Process:  block until all readCnBit==0, get wrl, wait for rl==0, set all readCntBits, then return (do the writing).
  float timeout=-1.0;//if negative, no timeout required
  int semid,val;
  int nReaders=0,i;
  struct timespec tout;
  struct sembuf operations[OPLEN];
  if(!PyArg_ParseTuple(args, "i|fi", &semid,&timeout,&nReaders)){
    printf("Usage: semid, timeout(optional),nReaders(optional)\n");
    return NULL;
  }
  tout.tv_sec=(int)timeout;
  tout.tv_nsec=(timeout-(int)timeout)*1e9;
  if(timeout==0.0)
    tout.tv_sec=tout.tv_nsec=0;
  for(i=2; i<nReaders+2; i++){//wait until all data has been read...
    operations[0].sem_num=i;
    operations[0].sem_op=0;//block until data has been read.
    operations[0].sem_flg=0;
    //val=semtimedop(semid, operations,OPLEN, timeout<0.0?NULL:&tout);
    if(timeout>=0.0)
	printf("\aWARNING: timeout=%g, but only <0 supported\n",timeout);
    Py_BEGIN_ALLOW_THREADS;//allow multi threading around blocking ops
    val=semop(semid, operations,OPLEN);
    Py_END_ALLOW_THREADS;
    if(val==-1){
      //failed to get block (probably timeout).
      printf("Couldn't block until all data was read... %d\n",i);
      if(errno==EAGAIN){//timeout
	printf("Timeout aquiring readcntbits for write lock: %s\n",strerror(errno));
      }else if(errno==EINTR){
	printf("Signal received while aquiring lock: %s\n",strerror(errno));
      }else if(errno==EINVAL){
	printf("Lock doesn't exist: %s\n",strerror(errno));
      }else if(errno==EIDRM){
	printf("Lock has been removed from system: %s\n",strerror(errno));
      }else{
	printf("Failed to get lock: %s\n",strerror(errno));
      }
      PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
      return NULL;
    }
  }
  //so, now everything has been read.  So, get the write lock...
  operations[0].sem_num = 0;//first semaphore only
  operations[0].sem_op = -1;// decrement by one (grab write lock)
  operations[0].sem_flg = 0;// do not permit undo
  // perform the operation on one semaphore
  //val=semtimedop(semid, operations,OPLEN, timeout<0.0?NULL:&tout);
  if(timeout>=0.0)
      printf("\aWARNING: timeout=%g, but only <0 supported\n",timeout);
  Py_BEGIN_ALLOW_THREADS;//allow multi threading around blocking ops
  val=semop(semid, operations,OPLEN);
  Py_END_ALLOW_THREADS;
  if(val==-1){
    //failed to get lock (probably timeout).
    printf("Couldn't grab write lock\n");
    printf("Sem values:%d %d\n",semctl(semid,0,GETVAL),semctl(semid,1,GETVAL));
    if(errno==EAGAIN){//timeout
      printf("Timeout aquiring write lock: %s\n",strerror(errno));
    }else if(errno==EINTR){
      printf("Signal received while aquiring lock: %s\n",strerror(errno));
    }else if(errno==EINVAL){
      printf("Lock doesn't exist: %s\n",strerror(errno));
    }else if(errno==EIDRM){
      printf("Lock has been removed from system: %s\n",strerror(errno));
    }else{
      printf("Failed to get lock: %s\n",strerror(errno));
    }
    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
    return NULL;
  }else{//got write lock successfully...
    //Now check that the readlock (sem1) is 0, ie not being read.
    operations[0].sem_num=1;//second (read) semaphore only
    operations[0].sem_op=0;//wait until no more readers
    operations[0].sem_flg=0;//don't permit undo.
    //val=semtimedop(semid, operations,OPLEN, timeout<0.0?NULL:&tout);
    if(timeout>=0.0)
	printf("\aWARNING: timeout=%g, but only <0 supported\n",timeout);
    Py_BEGIN_ALLOW_THREADS;//allow multi threading around blocking ops
    val=semop(semid, operations,OPLEN);
    Py_END_ALLOW_THREADS;
    if(val==-1){
      printf("Failed while waiting for readers to free\n");
      printf("Sem values:%d %d\n",semctl(semid,0,GETVAL),semctl(semid,1,GETVAL));
      //failed to get lock (probably timeout).
      if(errno==EAGAIN){//timeout
	printf("Timeout aquiring readlock for writelock: %s\n",strerror(errno));
      }else if(errno==EINTR){
	printf("Signal received while aquiring lock: %s\n",strerror(errno));
      }else if(errno==EINVAL){
	printf("Lock doesn't exist: %s\n",strerror(errno));
      }else if(errno==EIDRM){
	printf("Lock has been removed from system: %s\n",strerror(errno));
      }else{
	printf("Failed to get lock: %s\n",strerror(errno));
      }
      //should now release write lock since we failed to get the read lock
      operations[0].sem_num = 0;//first semaphore only
      operations[0].sem_op = +1;// increment by one (release write lock)
      operations[0].sem_flg = 0;// do not permit undo
      // perform the operation on one semaphore (this always succeeds).
      //semtimedop(semid, operations,OPLEN, timeout<0.0?NULL:&tout);
      if(timeout>=0.0)
	  printf("\aWARNING: timeout=%g, but only <0 supported\n",timeout);
      //Py_BEGIN_ALLOW_THREADS;//non-blocking - not needed.
      val=semop(semid, operations,OPLEN);
      //Py_END_ALLOW_THREADS;//allow multi threading around blocking ops
      PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
      return NULL;
    }else{//we've got the write lock, and there are no readers...
      for(i=2; i<2+nReaders; i++){//set all the data to unread...
	operations[0].sem_num=i;
	operations[0].sem_op=+1;//always succeeds.
	operations[0].sem_flg=0;
	//semtimedop(semid,operations,OPLEN,timeout<0.0?NULL:&tout);
	if(timeout>=0.0)
	    printf("\aWARNING: timeout=%g, but only <0 supported\n",timeout);
	//Py_BEGIN_ALLOW_THREADS;//non-blocking - not needed.
	val=semop(semid, operations,OPLEN);
	//Py_END_ALLOW_THREADS;//allow multi threading around blocking ops
	if(val==-1){
	  printf("shmemmodule Failed when setting flags to one\n");
	  PyErr_SetFromErrno(ShmemError);
	  return NULL;
	}
      }
    }
  }
  Py_INCREF(Py_None);
  return Py_None;//only if succeeded to get both locks.

}
static PyObject* shmem_freelock(PyObject *self, PyObject *args){//tuple of dimensions,type,name
  //Frees a semaphore lock - so that other processes can now grab it.
  //Free a write lock or non-specific lock (if read locks not used)...
  int semid,val;
  struct sembuf operations[OPLEN];
  if(!PyArg_ParseTuple(args, "i", &semid)){
    printf("Usage: semid\n");
    return NULL;
  }
  operations[0].sem_num = 0;// first semaphore only
  operations[0].sem_op = +1;// increment by one (free)
  operations[0].sem_flg = 0;// do not permit undo
  // perform the operation on one semaphore
  //Py_BEGIN_ALLOW_THREADS;//non-blocking - always succeeds.
  val=semop(semid, operations,OPLEN);
  //Py_END_ALLOW_THREADS;//allow multi threading around blocking ops
  if(val==-1){
    //failed to get lock (probably timeout).
    if(errno==EAGAIN){//timeout
      printf("Timeout freeing lock: %s\n",strerror(errno));
    }else if(errno==EINTR){
      printf("Signal received while freeing lock: %s\n",strerror(errno));
    }else if(errno==EINVAL){
      printf("Lock doesn't exist: %s\n",strerror(errno));
    }else if(errno==EIDRM){
      printf("Lock has been removed from system: %s\n",strerror(errno));
    }else{
      printf("Failed to free lock: %s\n",strerror(errno));
    }
    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;

}

static PyObject* shmem_getreadlock(PyObject *self, PyObject *args){//tuple of dimensions,type,name
  //grabs a semaphore read lock (waits until nothing has the write lock).
  //This increments the read lock count.
  //Process:  unset cntBit (if used), get wrl, inc rl, free wrl, then return (do the reading).
  float timeout=-1.0;//no timeout if <0
  int semid,val;
  int myBit=-1;
  struct timespec tout;
  struct sembuf operations[OPLEN];
  if(!PyArg_ParseTuple(args, "i|fi", &semid,&timeout,&myBit)){
    printf("Usage: semid, timeout(optional),myBit(optional)\n");
    return NULL;
  }
  tout.tv_sec=(int)timeout;
  tout.tv_nsec=(timeout-(int)timeout)*1e9;
  if(timeout==0.0)
    tout.tv_sec=tout.tv_nsec=0;
  if(myBit!=-1){//see if the readCountBit is set - if so, un set it.
    operations[0].sem_num=2+myBit;//myBit semaphore only
    operations[0].sem_op=-1;//decrement (unset or block until new data)
    operations[0].sem_flg=0;//no undo
    //val=semtimedop(semid,operations,OPLEN,timeout<0.0?NULL:&tout);
    if(timeout>=0.0)
	printf("\aWARNING: timeout=%g, but only <0 supported\n",timeout);
    Py_BEGIN_ALLOW_THREADS;//allow multi threading around blocking ops
    val=semop(semid, operations,OPLEN);
    Py_END_ALLOW_THREADS;
    if(val==-1){
      printf("Error while waiting for myBit %d to be set (new data)\n",myBit);
      if(errno==EAGAIN){//timeout
	printf("Timeout aquiring myBit lock to acquire read lock: %s\n",strerror(errno));
      }else if(errno==EINTR){
	printf("Signal received while aquiring lock: %s\n",strerror(errno));
      }else if(errno==EINVAL){
	printf("Lock doesn't exist: %s\n",strerror(errno));
      }else if(errno==EIDRM){
	printf("Lock has been removed from system: %s\n",strerror(errno));
      }else{
	printf("Failed to get lock: %s\n",strerror(errno));
      }
      PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
      return NULL;//couldn't get the write lock (to ensure nothing is writing)
    }else{
      printf("myBit %d obtained\n",myBit);
    }
  }
  //we now know that there is (or will be when the write lock is released) new
  //data.  So, get a readlock for this data (first getting the write lock)...
  operations[0].sem_num = 0;//first semaphore only (write lock)
  operations[0].sem_op = -1;// decrement by one (grab write lock)
  operations[0].sem_flg = 0;// do not permit undo
  // perform the operation on one semaphore
  //get the write lock
  //val=semtimedop(semid, operations,OPLEN, timeout<0.0?NULL:&tout);
  if(timeout>=0.0)
      printf("\aWARNING: timeout=%g, but only <0 supported\n",timeout);
  Py_BEGIN_ALLOW_THREADS;//allow multi threading around blocking ops
  val=semop(semid, operations,OPLEN);
  Py_END_ALLOW_THREADS;
  if(val==-1){
    //failed to get lock (probably timeout).
    if(errno==EAGAIN){//timeout
      printf("Timeout aquiring write lock to acquire read lock: %s\n",strerror(errno));
    }else if(errno==EINTR){
      printf("Signal received while aquiring lock: %s\n",strerror(errno));
    }else if(errno==EINVAL){
      printf("Lock doesn't exist: %s\n",strerror(errno));
    }else if(errno==EIDRM){
      printf("Lock has been removed from system: %s\n",strerror(errno));
    }else{
      printf("Failed to get lock: %s\n",strerror(errno));
    }
    //now re-set the readCntBit, since we haven't read...
    if(myBit!=-1){
      operations[0].sem_num=2+myBit;
      operations[0].sem_op=+1;//always succeeds
      operations[0].sem_flg=0;
      //semtimedop(semid,operations,OPLEN,timeout<0.0?NULL:&tout);//can't fail.
      if(timeout>=0.0)
	  printf("\aWARNING: timeout=%g, but only <0 supported\n",timeout);
      Py_BEGIN_ALLOW_THREADS;//allow multi threading around blocking ops
      val=semop(semid, operations,OPLEN);
      Py_END_ALLOW_THREADS;//allow multi threading around blocking ops
      if(val==-1){
	printf("shmemmodule Failed when resetting readCntBit\n");
	PyErr_SetFromErrno(ShmemError);
	return NULL;
      }
    }
    //didn't get the write lock, so no need to release it.
    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
    return NULL;//couldn't get the write lock (to ensure nothing is writing)
  }else{//got write lock successfully...

    //now increase the readlock count by 1.
    operations[0].sem_num=1;//second (read) semaphore only
    operations[0].sem_op=+1;//inc the readlock cnt by one (a new reader).
    operations[0].sem_flg=0;//don't permit undo.
    //semtimedop(semid, operations,OPLEN, timeout<0.0?NULL:&tout);//can't fail.
    if(timeout>=0.0)
	printf("\aWARNING: timeout=%g, but only <0 supported\n",timeout);
    //Py_BEGIN_ALLOW_THREADS;//non-blocking - not needed.
    val=semop(semid, operations,OPLEN);
    //Py_END_ALLOW_THREADS;//allow multi threading around blocking ops
    if(val==-1){
      printf("shmemmodule Failed when increasing readlock count\n");
      PyErr_SetFromErrno(ShmemError);
      return NULL;
    }
    //now release the write lock.
    operations[0].sem_num=0;//second (read) semaphore only
    operations[0].sem_op=+1;//release the write lock.
    operations[0].sem_flg=0;//don't permit undo.
    //semtimedop(semid, operations,OPLEN, timeout<0.0?NULL:&tout);//cant fail
    if(timeout>=0.0)
	printf("\aWARNING: timeout=%g, but only <0 supported\n",timeout);
    //Py_BEGIN_ALLOW_THREADS;//non-blocking - not needed
    val=semop(semid, operations,OPLEN);
    //Py_END_ALLOW_THREADS;//allow multi threading around blocking ops
    if(val==-1){
      printf("shmemmodule Failed when releasing writelock\n");
      PyErr_SetFromErrno(ShmemError);
      return NULL;
    }
  }
  Py_INCREF(Py_None);
  return Py_None;//only if succeeded to get both locks.

}

static PyObject* shmem_freereadlock(PyObject *self, PyObject *args){//tuple of dimensions,type,name
  //Frees a semaphore lock - so that other processes can now grab it.
  //Free a read lock...
  int semid,val;
  struct sembuf operations[OPLEN];
  if(!PyArg_ParseTuple(args, "i", &semid)){
    printf("Usage: semid\n");
    return NULL;
  }
  operations[0].sem_num = 1;// second (read) semaphore only
  operations[0].sem_op = -1;// decrement read count by one (free a read)
  operations[0].sem_flg = 0;// do not permit undo
  // perform the operation on one semaphore
  Py_BEGIN_ALLOW_THREADS;//allow multi threading around blocking ops
  val=semop(semid, operations,OPLEN);
  Py_END_ALLOW_THREADS;
  if(val==-1){
    //failed to get lock (probably timeout).
    if(errno==EAGAIN){//timeout
      printf("Timeout freeing read lock: %s\n",strerror(errno));
    }else if(errno==EINTR){
      printf("Signal received while freeing read lock: %s\n",strerror(errno));
    }else if(errno==EINVAL){
      printf("Read Lock doesn't exist: %s\n",strerror(errno));
    }else if(errno==EIDRM){
      printf("Read Lock has been removed from system: %s\n",strerror(errno));
    }else{
      printf("Failed to free read lock: %s\n",strerror(errno));
    }
    PyErr_SetFromErrno(ShmemError);//PyExc_IOError);
    return NULL;
  }
  //printf("readlock freed, semval: %d %d\n",semctl(semid,0,GETVAL),semctl(semid,1,GETVAL));
  Py_INCREF(Py_None);
  return Py_None;

}

static PyObject *shmem_queryNumericArray(PyObject *self,PyObject *args){
  PyObject *Narray,*tup;
  PyArrayObject *NumericArray;
  //int i;
  int raiseError=1;
  if(!PyArg_ParseTuple(args, "O|i", &Narray,&raiseError)){
    printf("Usage: NumericArray (optional, raiseError=1 if non-contiguous)\n");
    return NULL;
  }
  NumericArray=(PyArrayObject*)Narray;
  if(raiseError==1 && (NumericArray->flags&1)==0){
    printf("Error - Array is non-contiguous\n");
    return NULL;
  }
  switch(NumericArray->nd){
      case 1:
	  tup=Py_BuildValue("(i)",NumericArray->strides[0]);
	  break;
      case 2:
	  tup=Py_BuildValue("(ii)",NumericArray->strides[0],NumericArray->strides[1]);
	  break;
      case 3:
	  tup=Py_BuildValue("(iii)",NumericArray->strides[0],NumericArray->strides[1],NumericArray->strides[2]);
	  break;
      case 4:
	  tup=Py_BuildValue("(iiii)",NumericArray->strides[0],NumericArray->strides[1],NumericArray->strides[2],NumericArray->strides[3]);
	  break;
      default:
	  tup=Py_BuildValue("()");
  }
  return Py_BuildValue("(lO)",(long)NumericArray->data,tup);
}

static PyMethodDef ShmemMethods[] = {
  {"deleteAllSems",shmem_deleteAllSems,METH_VARARGS,
   "Delete all semaphore sets."},
  {"deleteSemid",shmem_deleteSemid,METH_VARARGS,
   "Delete a semid semaphore set."},
  {"getCurrentSemids", shmem_getCurrentSemids,METH_VARARGS,
   "Get array of current semids."},
  {"allFinished",shmem_allFinished,METH_VARARGS,
   "Set the allFinished flag..."},
  {"cleanUp",shmem_cleanUp,METH_VARARGS,
   "Clean up the simulation (set allFinished flag and remove semaphores)."},
  {"shmemopen",  shmem_shmemopen, METH_VARARGS,
   "Open a shared memory array."},
  {"shmemunlink",  shmem_shmemunlink, METH_VARARGS,
   "Unlink a shared memory array."},
  {"shmemunmap",  shmem_shmemunmap, METH_VARARGS,
   "Unmap a shared memory array."},
  {"open",  shmem_shmemopen, METH_VARARGS,
   "Open a shared memory array."},
  {"unlink",  shmem_shmemunlink, METH_VARARGS,
   "Unlink a shared memory array."},
  {"unmap",  shmem_shmemunmap, METH_VARARGS,
   "Unmap a shared memory array."},
  {"newsemid",(PyCFunction)shmem_newsemid, METH_VARARGS|METH_KEYWORDS,
   "Create a new semid"},
  {"initSemaphore",shmem_initSemaphore,METH_VARARGS,
   "Initialise semaphore values"},
  {"semop",shmem_semop,METH_VARARGS,
   "Perform an semop operation"},
  {"semdel",shmem_semdel,METH_VARARGS,
   "Delete a semaphore set"},
  {"getSemValue",shmem_getSemValue,METH_VARARGS,
   "Get semaphore current value"},
  {"acquireAndWait",shmem_acquireAndWait,METH_VARARGS,
   "Acquire a count semaphore and block on an event semaphore"},
  {"trigEvent",shmem_trigEvent,METH_VARARGS,
   "Trigger an event on a semaphore set"},
  {"newlock",(PyCFunction)shmem_newlock,METH_VARARGS|METH_KEYWORDS,
   "Create a new lock."},
  {"getlock",  shmem_getlock, METH_VARARGS,
   "Grab a (write) lock."},
  {"freelock",  shmem_freelock, METH_VARARGS,
   "Free a (write) lock."},
  {"getreadlock",  shmem_getreadlock, METH_VARARGS,
   "Grab a read lock."},
  {"freereadlock",  shmem_freereadlock, METH_VARARGS,
   "Free a read lock."},
  {"queryNumericArray",  shmem_queryNumericArray, METH_VARARGS,
   "Find out about a numeric array."},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};
//PyMODINIT_FUNC 
void initshmem(void)
{
  int i;
  PyObject *m;
  PyImport_AddModule("shmem");
  m=Py_InitModule("shmem", ShmemMethods);
  import_array();
  ShmemError = PyErr_NewException("shmem.error", NULL, NULL);
  Py_INCREF(ShmemError);
  PyModule_AddObject(m, "error", ShmemError);
  allFinished=0;
  for(i=0; i<SEMMNI; i++){
      semidArray[i]=-1;
  }
}
int
main(int argc, char *argv[])
{
  /* Pass argv[0] to the Python interpreter */
  Py_SetProgramName(argv[0]);
  
  /* Initialize the Python interpreter.  Required. */
  Py_Initialize();
  
  /* Add a static module */
  initshmem();
  return 0;
}
