#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/*
Here are functions and structures to deal with a memory limited sparse
array.  Basically, the user specifies a maximum number of data points
for the array and this is then never exceeded.  If a new value is
inserted into a full array, it will only be inserted if absolutely larger
than the minimum value, and this minimum value will then be removed
from the array.
If an existing value is changed to zero, it will be removed.

*/

#include "sparsemem.h"
#define SAFE_FREE(a) {if (a) {free(a); a = NULL;}}

dSpMem *smNewFromExisting(unsigned int ndata,unsigned int rows, unsigned int cols, double* data,unsigned int *indptr, unsigned int *rowind){
  //create a sparse mem struct from existing data/row/ind memory.
  dSpMem *sp;
  sp=(dSpMem*)malloc(sizeof(dSpMem));
  if(sp){
    sp->typ='d';
    sp->ndata=ndata;
    sp->rows=rows;
    sp->cols=cols;
    sp->alloced=0;
    sp->data=data;//(double*)malloc(sizeof(double)*ndata);
    sp->indptr=indptr;//(int*)malloc(sizeof(int)*(ndata+1));
    sp->rowind=rowind;//(int*)malloc(sizeof(int)*ndata);
    if(sp->indptr)
      memset(sp->indptr,0,sizeof(unsigned int)*(cols+1));
    sp->min=0.;
    sp->rowmin=0;
    sp->indmin=0;
    sp->cnt=0;
    if(!sp->data || !sp->indptr || !sp->rowind){
      SAFE_FREE(sp);
    }
  }
  return sp;
}



dSpMem *smNewSparseMem(unsigned int ndata,unsigned int rows,unsigned int cols){
  dSpMem *sp;
  sp=(dSpMem*)malloc(sizeof(dSpMem));
  if(sp){
    sp->typ='d';
    sp->ndata=ndata;
    sp->rows=rows;
    sp->cols=cols;
    sp->alloced=1;
    sp->data=(double*)malloc(sizeof(double)*ndata);
    sp->indptr=(unsigned int*)malloc(sizeof(unsigned int)*(cols+1));
    sp->rowind=(unsigned int*)malloc(sizeof(unsigned int)*ndata);
    if(sp->indptr)
      memset(sp->indptr,0,sizeof(unsigned int)*(cols+1));
    sp->min=0.;
    sp->rowmin=0;
    sp->indmin=0;
    sp->cnt=0;
    if(!sp->data || !sp->indptr || !sp->rowind){
      SAFE_FREE(sp->data);
      SAFE_FREE(sp->indptr);
      SAFE_FREE(sp->rowind);
      SAFE_FREE(sp);
    }
  }
  return sp;
}

void smFreeSparseMem(dSpMem *sp){
  if(sp!=NULL && sp->alloced){
    SAFE_FREE(sp->data);
    SAFE_FREE(sp->indptr);
    SAFE_FREE(sp->rowind);
  }
  SAFE_FREE(sp);
}

unsigned int getRowPosition(unsigned int s, unsigned int e, unsigned int row, unsigned int *rowind){
  //Find the position in rowind[s->e] where row should go.
  //Eventually this might be replaced by a quicksort algorithm...
  //Assumes that rowind[s->e] is ordered small to large.
  while(s<e && rowind[s]<row)
    s++;
  return s;
}

unsigned int getRowPositionQuick(unsigned int s, unsigned int e, unsigned int row, unsigned int *rowind){
  unsigned int m,so,eo;
  double dm;
  if(e-s<10)//don't bother with quicksort method...
    return getRowPosition(s,e,row,rowind);
  eo=e;//save the original positions.
  so=s;
  dm=(s+e)/2;//compute the middle position.
  m=(unsigned int)(dm+0.5);
  //Finish when either m=so, m=eo or rowind[m-1]<row && rowind[m]>=row.
  while(m>so && m<eo && !(rowind[m-1]<row && rowind[m]>=row)){
    if(rowind[m]<row){
      s=m;
      dm=(dm+e)/2;
    }else if(rowind[m]>row){
      e=m;
      dm=(s+dm)/2;
    }else{//rowind[m]==row - exact match.

    }
    m=(unsigned int)(dm+0.5);
  }
  return m;
}

unsigned int smGetIndxForRow(dSpMem *sp,unsigned int row, unsigned int col, unsigned int suggestion){
  //Finds the index for data[]/rowind[] which corresponds to position
  //row,col, or the nearest next row available, if this entry does not exist.
  //Suggestion can be provided as a suggested starting point for the
  //search, or can be -1 to be ignored.
  unsigned int s=sp->indptr[col];
  unsigned int e=sp->indptr[col+1];
  if(suggestion<s || suggestion>=e){//bad suggestion
    //do a quick search.
    suggestion=getRowPositionQuick(s,e,row,sp->rowind);
  }else{//search around suggestion...
    while(suggestion>=s && sp->rowind[suggestion]>row)
      suggestion--;
    while(suggestion<e && sp->rowind[suggestion]<row)
      suggestion++;
  }
  return suggestion;
}


unsigned int findNewMin(dSpMem *sp){
  double minval;
  unsigned int colmin;
  unsigned int i;
  unsigned int rowmin,j;
  minval=-1;//fabs(sp->data[0]);
  rowmin=0;
  colmin=0;
  for(i=0; i<sp->cols; i++){
    for(j=sp->indptr[i]; j<sp->indptr[i+1]; j++){
      if(minval<0 || fabs(sp->data[j])<minval){
	minval=fabs(sp->data[j]);
	rowmin=j;
	colmin=i;
      }
    }
  }
  sp->rowmin=rowmin;
  sp->indmin=colmin;
  sp->min=minval;
  return rowmin;
}

void removeValue(dSpMem *sp,unsigned int rowind, unsigned int col){
  //assumes the value at this position does exist.
  unsigned int findmin=0,c;
  unsigned int elementsMove;
  if(rowind==sp->rowmin)
    findmin=1;//we're removing the current min, so will need to re-find it.
  elementsMove=sp->cnt-(rowind+1);
  memmove(&sp->data[rowind],&sp->data[rowind+1],elementsMove*sizeof(double));
  memmove(&sp->rowind[rowind],&sp->rowind[rowind+1],elementsMove*sizeof(unsigned int));
  for(c=col+1; c<=sp->cols; c++){
    sp->indptr[c]--;
  }
  sp->cnt--;
  if(findmin==1){
    findNewMin(sp);
  }else{//may need to adjust sp->rowmin if it has been shifted.
    if(sp->rowmin>rowind)
      sp->rowmin--;
  }
}


int smReplaceData(dSpMem *sp,unsigned int rowind,unsigned int col, double data){
  double fdata;
  //This assumes the user knows that the entry for thsi data already exists.
  //If data==0, the entry will be removed.
  if(data==0.){
    removeValue(sp,rowind,col);
    return 0;
  }else{
    fdata=fabs(data);
    sp->data[rowind]=data;
    if(fdata<sp->min || sp->min==0){//we're the new minimum...
      //findNewMin(sp);
      sp->min=fdata;
      sp->rowmin=rowind;
      sp->indmin=col;
    }else if(sp->rowmin==rowind){//we're replacing the min, so find the new one
      findNewMin(sp);
    }
  }
  return 1;
}


double smGet(dSpMem *sp,unsigned int row, unsigned int col){
  unsigned int e,i;
  if(row<0 || col<0 || row>=sp->rows || col>=sp->cols){
    printf("sparsemem.c: Getting from illegal row/col: %d %d\n",row,col);
    return 0.;
  }
  e=(unsigned int)sp->indptr[col+1];
  i=getRowPositionQuick((unsigned int)sp->indptr[col],e,row,sp->rowind);
  if(i<e && sp->rowind[i]==row){
    return sp->data[i];
  }
  return 0.;//not found...
}

int smInsert(dSpMem *sp,unsigned int row, unsigned int col, double data){
  //Attempt to insert a value at row,col.  Returns 1 if inserted, 0 otherwise, -1 if row/col illegal.
  int inserted=1;
  double fdata;
  unsigned int i,e;
  unsigned int c;
  unsigned int elementsMove;
  if(row<0 || col<0 || row>=sp->rows || col>=sp->cols){
    printf("sparsemem.c - inserting at illegal row/col: %d %d\n",row,col);
    return -1;
  }
  fdata=fabs(data);
  if(data==0.){
    //look to see if an entry exists - if so, remove it... if not, do nothing.
    e=(unsigned int)sp->indptr[col+1];
    i=getRowPositionQuick((unsigned int)sp->indptr[col],e,row,sp->rowind);
    if(i<e && sp->rowind[i]==row){
      //remove this value...
      removeValue(sp,i,col);
    }
    return 0;
  }
  
  //take care of a special case - where data<sp->min but row/col is at the position of sp->min - it should therefore be inserted at this position, and be the new min... also take care of the case where entry at row/col already exists - should always be replaced in this case, even if the value is less than min.  Note, if replacing with zero, it should be removed.
  if(sp->cnt==sp->ndata){//no space for the value
    //look to see if an entry exists at this point - if so, replace it, regardless of the new value...
    e=(unsigned int)sp->indptr[col+1];
    i=getRowPositionQuick((unsigned int)sp->indptr[col],e,row,sp->rowind);
    if(i<e && sp->rowind[i]==row){
      smReplaceData(sp,i,col,data);
    }else if(fdata>sp->min){//remove the current minimum value... and insert this one.
      //for efficiency, want to remove and add at the same time - to avoid shifting some data twice.
      //The point we want to remove is in col sp->indmin, row sp->rowind[sp->rowmin]
      //find the point at which we wish to insert.
      //i=sp->indptr[col];
      //e=sp->indptr[col+1];
      //move to the right point...
      //i=getRowPosition(i,e,row,sp->rowind);
      //Note, now, rowind[i]>row.
      //So, we want to insert at this point.
      if(sp->rowmin==i){//the value we're overwriting happens to be at the same point we wish to insert (no memmove needed). 
	sp->rowind[i]=row;
	sp->data[i]=data;
	//may need to increment the indptr, if i==rowmin but i is one to be inserted at the end of a row!
	if(sp->indmin<col){
	  for(c=sp->indmin+1; c<=col; c++){
	    sp->indptr[c]--;
	  }
	}else if(sp->indmin>col){
	  for(c=col+1; c<=sp->indmin; c++){
	    sp->indptr[c]++;
	  }
	}
      }else if(sp->rowmin<i){
	elementsMove=i-sp->rowmin;
	memmove(&sp->data[sp->rowmin],&sp->data[sp->rowmin+1],elementsMove*sizeof(double));
	memmove(&sp->rowind[sp->rowmin],&sp->rowind[sp->rowmin+1],elementsMove*sizeof(unsigned int));
	for(c=sp->indmin+1; c<=col; c++){
	  sp->indptr[c]--;
	}
	i--;
	sp->data[i]=data;
	sp->rowind[i]=row;
      }else{//sp->rowmin>i
	elementsMove=sp->rowmin-i;
	memmove(&sp->data[i+1],&sp->data[i],elementsMove*sizeof(double));
	memmove(&sp->rowind[i+1],&sp->rowind[i],elementsMove*sizeof(unsigned int));
	for(c=col+1; c<=sp->indmin; c++){
	  sp->indptr[c]++;
	}
	sp->data[i]=data;
	sp->rowind[i]=row;
      }
      //now find the new min value and its position...
      findNewMin(sp);
    }else{//data does not warrent entry - not large enough... so do nothing
      inserted=0;
    }
  }else{//space in the array... so insert it and find the new min.
    i=sp->indptr[col];
    e=sp->indptr[col+1];
    //move to the right point...
    i=getRowPositionQuick(i,e,row,sp->rowind);//get the position at which it should be inserted.
    if(i==e){//reached the end for this row, so insert here...
      elementsMove=sp->cnt-i;
      if(elementsMove>0){
	memmove(&sp->data[i+1],&sp->data[i],elementsMove*sizeof(double));
	memmove(&sp->rowind[i+1],&sp->rowind[i],elementsMove*sizeof(unsigned int));
	if(sp->rowmin>=i)
	  sp->rowmin++;//increment the min position.
      }
      sp->data[i]=data;
      sp->rowind[i]=row;
      for(c=col+1; c<sp->cols+1; c++){
	sp->indptr[c]++;
      }
    }else if(sp->rowind[i]==row){//replace this value
      smReplaceData(sp,i,col,data);
      sp->cnt--;//will be re-incremented later...
      //sp->data[i]=data;
      //if(i==sp->rowmin)
      //findNewMin(sp);
    }else{//insert before this value (ie shift this value up).
      elementsMove=sp->cnt-i;
      if(elementsMove>0){
	memmove(&sp->data[i+1],&sp->data[i],elementsMove*sizeof(double));
	memmove(&sp->rowind[i+1],&sp->rowind[i],elementsMove*sizeof(unsigned int));
	if(sp->rowmin>=i)
	  sp->rowmin++;//increment the min position.
      }
      sp->data[i]=data;
      sp->rowind[i]=row;
      for(c=col+1; c<sp->cols+1; c++){
	sp->indptr[c]++;
      }
    }
    if(fdata<sp->min || sp->min==0.){
      //update min to the new value...
      sp->min=fdata;
      sp->rowmin=i;
      sp->indmin=col;
    }
    sp->cnt++;
  }
  return inserted;
}



/****************************************************************************
Now the floating point version...
****************************************************************************/

fSpMem *smNewFromExistingFloat(unsigned int ndata,unsigned int rows, unsigned int cols, float* data,unsigned int *indptr, unsigned int *rowind){
  //create a sparse mem struct from existing data/row/ind memory.
  fSpMem *sp;
  sp=(fSpMem*)malloc(sizeof(fSpMem));
  if(sp){
    sp->typ='f';
    sp->ndata=ndata;
    sp->rows=rows;
    sp->cols=cols;
    sp->alloced=0;
    sp->data=data;//(float*)malloc(sizeof(float)*ndata);
    sp->indptr=indptr;//(int*)malloc(sizeof(int)*(ndata+1));
    sp->rowind=rowind;//(int*)malloc(sizeof(int)*ndata);
    if(sp->indptr)
      memset(sp->indptr,0,sizeof(unsigned int)*(cols+1));
    sp->min=0.;
    sp->rowmin=0;
    sp->indmin=0;
    sp->cnt=0;
    if(!sp->data || !sp->indptr || !sp->rowind){
      SAFE_FREE(sp);
    }
  }
  return sp;
}



fSpMem *smNewSparseMemFloat(unsigned int ndata,unsigned int rows,unsigned int cols){
  fSpMem *sp;
  sp=(fSpMem*)malloc(sizeof(fSpMem));
  if(sp){
    sp->typ='f';
    sp->ndata=ndata;
    sp->rows=rows;
    sp->cols=cols;
    sp->alloced=1;
    sp->data=(float*)malloc(sizeof(float)*ndata);
    sp->indptr=(unsigned int*)malloc(sizeof(unsigned int)*(cols+1));
    sp->rowind=(unsigned int*)malloc(sizeof(unsigned int)*ndata);
    if(sp->indptr)
      memset(sp->indptr,0,sizeof(unsigned int)*(cols+1));
    sp->min=0.;
    sp->rowmin=0;
    sp->indmin=0;
    sp->cnt=0;
    if(!sp->data || !sp->indptr || !sp->rowind){
      SAFE_FREE(sp->data);
      SAFE_FREE(sp->indptr);
      SAFE_FREE(sp->rowind);
      SAFE_FREE(sp);
    }
  }
  return sp;
}

void smFreeSparseMemFloat(fSpMem *sp){
  if(sp!=NULL && sp->alloced){
    SAFE_FREE(sp->data);
    SAFE_FREE(sp->indptr);
    SAFE_FREE(sp->rowind);
  }
  SAFE_FREE(sp);
}

unsigned int getRowPositionFloat(unsigned int s, unsigned int e, unsigned int row, unsigned int *rowind){
  //Find the position in rowind[s->e] where row should go.
  //Eventually this might be replaced by a quicksort algorithm...
  //Assumes that rowind[s->e] is ordered small to large.
  while(s<e && rowind[s]<row)
    s++;
  return s;
}

unsigned int getRowPositionQuickFloat(unsigned int s, unsigned int e, unsigned int row, unsigned int *rowind){
  unsigned int m,so,eo;
  float dm;
  if(e-s<10)//don't bother with quicksort method...
    return getRowPositionFloat(s,e,row,rowind);
  eo=e;//save the original positions.
  so=s;
  dm=(s+e)/2;//compute the middle position.
  m=(unsigned int)(dm+0.5);
  //Finish when either m=so, m=eo or rowind[m-1]<row && rowind[m]>=row.
  while(m>so && m<eo && !(rowind[m-1]<row && rowind[m]>=row)){
    if(rowind[m]<row){
      s=m;
      dm=(dm+e)/2;
    }else if(rowind[m]>row){
      e=m;
      dm=(s+dm)/2;
    }else{//rowind[m]==row - exact match.

    }
    m=(unsigned int)(dm+0.5);
  }
  return m;
}

unsigned int smGetIndxForRowFloat(fSpMem *sp,unsigned int row, unsigned int col, unsigned int suggestion){
  //Finds the index for data[]/rowind[] which corresponds to position
  //row,col, or the nearest next row available, if this entry does not exist.
  //Suggestion can be provided as a suggested starting point for the
  //search, or can be -1 to be ignored.
  unsigned int s=sp->indptr[col];
  unsigned int e=sp->indptr[col+1];
  if(suggestion<s || suggestion>=e){//bad suggestion
    //do a quick search.
    suggestion=getRowPositionQuickFloat(s,e,row,sp->rowind);
  }else{//search around suggestion...
    while(suggestion>=s && sp->rowind[suggestion]>row)
      suggestion--;
    while(suggestion<e && sp->rowind[suggestion]<row)
      suggestion++;
  }
  return suggestion;
}


unsigned int findNewMinFloat(fSpMem *sp){
  float minval;
  unsigned int colmin;
  unsigned int i;
  unsigned int rowmin,j;
  minval=-1;//fabsf(sp->data[0]);
  rowmin=0;
  colmin=0;
  for(i=0; i<sp->cols; i++){
    for(j=sp->indptr[i]; j<sp->indptr[i+1]; j++){
      /*if(j>=sp->ndata){
	printf("Error in sparsemem findNewMinFloat - j=%u cnt=%u ndata=%u i=%u ind=%u %u\n",j,sp->cnt,sp->ndata,i,sp->indptr[i],sp->indptr[i+1]);
	for(i=0; i<=20; i++)
	  printf("%u  \n",sp->indptr[i]);

	exit(0);

	}*/
      if(minval<0 || fabsf(sp->data[j])<minval){
	minval=fabsf(sp->data[j]);
	rowmin=j;
	colmin=i;
      }
    }
  }
  sp->rowmin=rowmin;
  sp->indmin=colmin;
  sp->min=minval;
  if(minval==-1)
    printf("Error - sparsemem -> findNewMinFloat found -1\n");
  //printf("New min %g at %u %u (%u %u %u) %u %u\n",minval,rowmin,colmin,sp->indptr[0],sp->indptr[1],sp->indptr[2],i,j);
  return rowmin;
}

void removeValueFloat(fSpMem *sp,unsigned int rowind, unsigned int col){
  //assumes the value at this position does exist.
  int findmin=0,c;
  unsigned int elementsMove;
  if(rowind==sp->rowmin)
    findmin=1;//we're removing the current min, so will need to re-find it.
  elementsMove=sp->cnt-(rowind+1);
  memmove(&sp->data[rowind],&sp->data[rowind+1],elementsMove*sizeof(float));
  memmove(&sp->rowind[rowind],&sp->rowind[rowind+1],elementsMove*sizeof(unsigned int));
  for(c=col+1; c<=sp->cols; c++){
    sp->indptr[c]--;
  }
  sp->cnt--;
  if(findmin==1){
    findNewMinFloat(sp);
  }else{//may need to adjust sp->rowmin if it has been shifted.
    if(sp->rowmin>rowind)
      sp->rowmin--;
  }
}


int smReplaceDataFloat(fSpMem *sp,unsigned int rowind,unsigned int col, float data){
  float fdata;
  //This assumes the user knows that the entry for thsi data already exists.
  //If data==0, the entry will be removed.
  //printf("smReplaceDataFloat\n");
  if(data==0.){
    removeValueFloat(sp,rowind,col);
    return 0;
  }else{
    fdata=fabsf(data);
    sp->data[rowind]=data;
    if(fdata<sp->min || sp->min==0){//we're the new minimum...
      //findNewMin(sp);
      sp->min=fdata;
      sp->rowmin=rowind;
      sp->indmin=col;
    }else if(sp->rowmin==rowind){//we're replacing the min, so find the new one
      findNewMinFloat(sp);
    }
  }
  return 1;
}


float smGetFloat(fSpMem *sp,unsigned int row, unsigned int col){
  unsigned int e,i;
  if(row<0 || col<0 || row>=sp->rows || col>=sp->cols){
    printf("sparsemem.c: Getting from illegal row/col: %d %d\n",row,col);
    return 0.;
  }
  e=(unsigned int)sp->indptr[col+1];
  i=getRowPositionQuickFloat((unsigned int)sp->indptr[col],e,row,sp->rowind);
  if(i<e && sp->rowind[i]==row){
    return sp->data[i];
  }
  return 0.;//not found...
}

int smInsertFloat(fSpMem *sp,unsigned int row, unsigned int col, float data){
  //Attempt to insert a value at row,col.  Returns 1 if inserted, 0 otherwise, -1 if row/col illegal.
  int inserted=1;
  float fdata;
  unsigned int i,e;
  unsigned int c;
  unsigned int elementsMove;
  //printf("%u %u %u %u %u %g %g\n",sp->indptr[0],sp->indptr[1],sp->cnt,row,col,data,sp->min);
  if(row<0 || col<0 || row>=sp->rows || col>=sp->cols){
    printf("sparsemem.c - inserting at illegal row/col: %d %d\n",row,col);
    return -1;
  }
  fdata=fabsf(data);
  if(data==0.){
    //look to see if an entry exists - if so, remove it... if not, do nothing.
    e=(unsigned int)sp->indptr[col+1];
    i=getRowPositionQuickFloat((unsigned int)sp->indptr[col],e,row,sp->rowind);
    if(i<e && sp->rowind[i]==row){
      //remove this value...
      removeValueFloat(sp,i,col);
    }
    return 0;
  }
  
  //take care of a special case - where data<sp->min but row/col is at the position of sp->min - it should therefore be inserted at this position, and be the new min... also take care of the case where entry at row/col already exists - should always be replaced in this case, even if the value is less than min.  Note, if replacing with zero, it should be removed.
  if(sp->cnt==sp->ndata){//no space for the value
    //look to see if an entry exists at this point - if so, replace it, regardless of the new value...
    e=(unsigned int)sp->indptr[col+1];
    i=getRowPositionQuickFloat((unsigned int)sp->indptr[col],e,row,sp->rowind);
    //printf("%u %u %u %g %u %g %u %u %u\n",i,e,sp->rowmin,sp->data[sp->rowmin],sp->indmin,sp->min,sp->indptr[0],sp->indptr[1],sp->indptr[2]);
    if(i<e && sp->rowind[i]==row){
      smReplaceDataFloat(sp,i,col,data);
    }else if(fdata>sp->min){//remove the current minimum value... and insert this one.
      //for efficiency, want to remove and add at the same time - to avoid shifting some data twice.
      //The point we want to remove is in col sp->indmin, row sp->rowind[sp->rowmin]
      //find the point at which we wish to insert.
      //i=sp->indptr[col];
      //e=sp->indptr[col+1];
      //move to the right point...
      //i=getRowPosition(i,e,row,sp->rowind);
      //Note, now, rowind[i]>row.
      //So, we want to insert at this point.
      if(sp->rowmin==i){//the value we're overwriting happens to be at the same point we wish to insert (no memmove needed). 
	sp->rowind[i]=row;
	sp->data[i]=data;
	//may need to increment the indptr, if i==rowmin but i is one to be inserted at the end of a row!
	if(sp->indmin<col){
	  for(c=sp->indmin+1; c<=col; c++){
	    sp->indptr[c]--;
	  }
	}else if(sp->indmin>col){
	  for(c=col+1; c<=sp->indmin; c++){
	    sp->indptr[c]++;
	  }
	}
      }else if(sp->rowmin<i){
	elementsMove=i-sp->rowmin;
	memmove(&sp->data[sp->rowmin],&sp->data[sp->rowmin+1],elementsMove*sizeof(float));
	memmove(&sp->rowind[sp->rowmin],&sp->rowind[sp->rowmin+1],elementsMove*sizeof(unsigned int));
	for(c=sp->indmin+1; c<=col; c++){
	  sp->indptr[c]--;
	}
	i--;
	sp->data[i]=data;
	sp->rowind[i]=row;
      }else{//sp->rowmin>i
	elementsMove=sp->rowmin-i;
	memmove(&sp->data[i+1],&sp->data[i],elementsMove*sizeof(float));
	memmove(&sp->rowind[i+1],&sp->rowind[i],elementsMove*sizeof(unsigned int));
	for(c=col+1; c<=sp->indmin; c++){
	  sp->indptr[c]++;
	}
	sp->data[i]=data;
	sp->rowind[i]=row;
      }
      //now find the new min value and its position...
      findNewMinFloat(sp);
    }else{//data does not warrent entry - not large enough... so do nothing
      inserted=0;
    }
  }else{//space in the array... so insert it and find the new min.
    i=sp->indptr[col];
    e=sp->indptr[col+1];
    //move to the right point...
    i=getRowPositionQuickFloat(i,e,row,sp->rowind);//get the position at which it should be inserted.
    if(i==e){//reached the end for this row, so insert here...
      elementsMove=sp->cnt-i;
      if(elementsMove>0){
	memmove(&sp->data[i+1],&sp->data[i],elementsMove*sizeof(float));
	memmove(&sp->rowind[i+1],&sp->rowind[i],elementsMove*sizeof(unsigned int));
	if(sp->rowmin>=i)
	  sp->rowmin++;//increment the min position.
      }
      sp->data[i]=data;
      sp->rowind[i]=row;
      for(c=col+1; c<sp->cols+1; c++){
	sp->indptr[c]++;
      }
    }else if(sp->rowind[i]==row){//replace this value
      smReplaceDataFloat(sp,i,col,data);
      sp->cnt--;//will be re-incremented later...
      //sp->data[i]=data;
      //if(i==sp->rowmin)
      //findNewMin(sp);
    }else{//insert before this value (ie shift this value up).
      elementsMove=sp->cnt-i;
      if(elementsMove>0){
	memmove(&sp->data[i+1],&sp->data[i],elementsMove*sizeof(float));
	memmove(&sp->rowind[i+1],&sp->rowind[i],elementsMove*sizeof(unsigned int));
	if(sp->rowmin>=i)
	  sp->rowmin++;//increment the min position.
      }
      sp->data[i]=data;
      sp->rowind[i]=row;
      for(c=col+1; c<sp->cols+1; c++){
	sp->indptr[c]++;
      }
    }
    if(fdata<sp->min || sp->min==0.){
      //update min to the new value...
      sp->min=fdata;
      sp->rowmin=i;
      sp->indmin=col;
    }
    sp->cnt++;
  }
  return inserted;
}


int sparseDotSparse(float *dataa,unsigned int *colind,unsigned int *aindptr,int ninda,float *datab,unsigned int *rowind,unsigned int *bindptr,int nindb,float **resdataptr,unsigned int **rescolindptr,unsigned int* resindptr,unsigned int *resn,int tn){
  //perform dot product of a (csr) and b (csc), returning result as csr format.  Size of resdata must be equal to dataa + datab size.
  int respos=0;
  float *resdata;
  unsigned int *rescolind;
  float val;
  unsigned int i,j,posa,posb,oldsize;
  float *tmpf;
  unsigned int *tmpi;
  resdata=*resdataptr;
  rescolind=*rescolindptr;
  resindptr[0]=0;
  for(i=0; i<ninda-1; i++){
    if(tn==0){
      printf("%d/%d     \r",i,ninda);
      fflush(NULL);
    }
    for(j=0; j<nindb-1; j++){
      posa=aindptr[i];
      posb=bindptr[j];
      val=0.;
      while(posa<aindptr[i+1] && posb<bindptr[j+1]){
	if(colind[posa]<rowind[posb])
	  posa++;
	else if(colind[posa]>rowind[posb])
	  posb++;
	else{
	  val+=dataa[posa]*datab[posb];
	  posa++;
	  posb++;
	}
      }
      if(val!=0.){
	if(respos==(*resn)){
	  oldsize=*resn;
	  (*resn)=(int)((*resn)*1.2);
	  if(*resn<oldsize){
	    printf("\nWarning - using more than 4G values...ERROR\n");
	    return 1;
	  }
	  printf("Resizing to %d (thread %d)\n",*resn,tn);
	  tmpf=malloc(sizeof(float)*(*resn));
	  memcpy(tmpf,resdata,sizeof(float)*oldsize);
	  free(resdata);
	  *resdataptr=tmpf;
	  resdata=*resdataptr;
	  tmpi=malloc(sizeof(int)*(*resn));
	  memcpy(tmpi,rescolind,sizeof(int)*oldsize);
	  free(rescolind);
	  *rescolindptr=tmpi;
	  rescolind=*rescolindptr;
	}
	resdata[respos]=val;
	rescolind[respos]=j;
	respos++;
      }
    }
    resindptr[i+1]=respos;
  }
  printf("Done thread (%d) %d         \n",tn,i);
  return 0;
}

int denseDotCSC(float *a,int arows, int acols, float *datab,unsigned int *browind, unsigned int *bindptr,int nindb,float **resdataptr,unsigned int **resrcindptr, unsigned int *resindptr, unsigned int *resn,float valmin,int restype,unsigned int resnmax,int pr){
  //perform dot product of dense matrix with CSC matrix.  Output is either csc or csr depending on resIsCsc value.
  //a is a arows x acols matrix
  //Return values: 0 = okay, -1 = realloc failed, >0 maxelements reached at row N.
  int respos=0;
  float *resdata;
  unsigned int *resrcind;
  float val;
  long i,j,oldsize;
  unsigned int posb;
  float *tmpf;
  unsigned int *tmpi;
  resdata=*resdataptr;
  //printf("%p\n",resdata);
  if(restype==0){//result is csc
    resrcind=*resrcindptr;
    resindptr[0]=0;
    for(i=0; i<nindb-1; i++){
      if(pr){
	printf("%ld/%d    \r",i,nindb-1);
	fflush(NULL);
      }
      for(j=0; j<arows; j++){
	posb=bindptr[i];
	val=0.;
	while(posb<bindptr[i+1]){
	  val+=datab[posb]*a[j*acols+browind[posb]];
	  posb++;
	}
	if(fabsf(val)>=valmin){
	  if(respos==(*resn)){
	    oldsize=*resn;
	    (*resn)=((unsigned int)(*resn)*1.2);
	    if(*resn>resnmax)
	      *resn=resnmax;
	    if(*resn<=oldsize){
	      printf("\nWARNING - can no longer increase size - please select a larger minimum value - got to %ld/%d (oldsize %ld requested %u\n",i,nindb-1,oldsize,*resn);
	      return (int)i;
	    }
	    printf("Resizing to %u at row %ld\n",*resn,i);
	    //tmpf=malloc(sizeof(float)*(*resn));
	    //memcpy(tmpf,resdata,sizeof(float)*oldsize);
	    //free(resdata);
	    if((tmpf=realloc(resdata,sizeof(float)*(*resn)))==NULL){
	      printf("realloc failed\n");
	      return -1;
	    }
	    *resdataptr=tmpf;
	    resdata=*resdataptr;
	    //tmpi=malloc(sizeof(int)*(*resn));
	    //memcpy(tmpi,resrcind,sizeof(int)*oldsize);
	    //free(resrcind);
	    if((tmpi=realloc(resrcind,sizeof(int)*(*resn)))==NULL){
	      printf("realloc tmpi failed\n");
	      return -1;
	    }
	    *resrcindptr=tmpi;
	    resrcind=*resrcindptr;
	  }
	  resdata[respos]=val;
	  resrcind[respos]=(unsigned int)j;
	  respos++;
        }
      }
      resindptr[i+1]=respos;
    }
  }else if(restype==1){//result is CSR...
    resrcind=*resrcindptr;
    resindptr[0]=0;
    for(i=0; i<arows; i++){
      if(pr){
	printf("%ld/%d     \r",i,arows);
	fflush(NULL);
      }
      for(j=0; j<nindb-1; j++){
	posb=bindptr[j];
	val=0.;
	while(posb<bindptr[j+1]){
	  val+=datab[posb]*a[i*acols+browind[posb]];
	  posb++;
	}
	if(fabsf(val)>=valmin){
	  if(respos==(*resn)){
	    oldsize=*resn;
	    (*resn)=((int)(*resn)*1.2);
	    if(*resn>resnmax)
	      *resn=resnmax;
	    if(*resn<=oldsize){
	      printf("WARNING - can no longer increase size - please select a larger minimum value - got to %ld/%d (oldsize %ld requested %u\n",i,arows,oldsize,*resn);
	      return (int)i;
	    }
	    printf("Resizing to %u at row %ld\n",*resn,i);
	    if((tmpf=realloc(resdata,sizeof(float)*(*resn)))==NULL){
	      printf("realloc failed\n");
	      return -1;
	    }
	    //tmpf=malloc(sizeof(float)*(*resn));
	    //memcpy(tmpf,resdata,sizeof(float)*oldsize);
	    //free(resdata);
	    *resdataptr=tmpf;
	    resdata=*resdataptr;
	    //tmpi=malloc(sizeof(int)*(*resn));
	    //memcpy(tmpi,resrcind,sizeof(int)*oldsize);
	    //free(resrcind);
	    if((tmpi=realloc(resrcind,sizeof(int)*(*resn)))==NULL){
	      printf("realloc tmpi failed\n");
	      return -1;
	    }
	    *resrcindptr=tmpi;
	    resrcind=*resrcindptr;
	  }
	  resdata[respos]=val;
	  resrcind[respos]=(unsigned int)j;
	  respos++;

	}
      }
      resindptr[i+1]=respos;
    }
  }else{//result is dense...
    //printf("arows\n");
    for(i=0; i<arows; i++){
      if(pr){
	printf("%ld/%d\r",i,arows);
	fflush(NULL);
      }
      for(j=0; j<nindb-1; j++){
	posb=bindptr[j];
	val=0.;
	while(posb<bindptr[j+1]){
	  val+=datab[posb]*a[i*acols+browind[posb]];
	  posb++;
	}
	if(fabsf(val)<valmin){
	  val=0.;
	}else if(isnan(val)){
	  val=0.;
	  //printf("Got NAN at %d %d\n",i,j);
	}
	//printf("%d %d %d %d %d %g\n",i*acols+browind[posb],i,j,nindb,i*(nindb-1)+j,val);
	//resdata[0]=val;
	resdata[(long)i*(long)(nindb-1)+(long)j]=val;
      }
    }
  }
  if(pr)
    printf("\nDone %ld                   \n",i);
  return 0;
}

int countInstances(float *a,long size,float *vals,long *cnt, int n){
  //vals is a sorted array containing the histogram steps for which the number of instances in a should be counted.  Results are put in cnt.
  long i;
  int j,k;
  memset(cnt,0,sizeof(long)*n);
  for(i=0; i<size; i++){
    for(j=n-1; j>=0; j--){
      if(fabsf(a[i])>=vals[j]){
	for(k=0; k<=j; k++){
	  cnt[k]++;
	}
	j=-1;
      }
    }
  }
  return 0;
}
int sparsifyCsr(float *a,int nrows,int ncols,int lda,int ldb,float *data,unsigned int *colind, unsigned int *indptr,unsigned int n,float val){
  //sparsify a dense matrix...
  int i,j;
  long apos=0;
  unsigned int pos=0;
  printf("sparsifyCsr nrows %d ncols %d lda %d ldb %d n %d val %g\n",nrows,ncols,lda,ldb,n,val);
  indptr[0]=0;
  for(i=0; i<nrows; i++){
    apos=(long)i*(long)lda;
    for(j=0; j<ncols; j++){
      //printf("%d %d\n",i,j);
      if(fabsf(a[apos])>=val){
	if(pos==n){
	  printf("sparsifyCsr error - out of space (pos %d at %d/%d, %d/%d)\n",pos,i,nrows,j,ncols);
	  return 1;
	}
	//printf("%ld %d\n",apos,pos);
	data[pos]=a[apos];
	colind[pos]=j;
	pos++;
      }
      apos+=(long)ldb;
      
    }
    indptr[i+1]=pos;
  }
  return 0;
}

int densifyCsr(float *a, int nrows, int ncols, int lda, int ldb, float *data, unsigned int *colind,unsigned int *indptr,unsigned int n){
  //put the sparse matrix into a.  For some reason, csr.todense() segments on a large (18GB) matrix.
  //For a c contiguous matrix, ldb==1, lda==ncols.
  long indx;
  int i,j;
  for(i=0; i<nrows; i++){
    for(j=indptr[i]; j<indptr[i+1]; j++){
      indx=(long)colind[j]*(long)ldb+(long)i*lda;
      a[indx]=data[j];
    }
  }
  return 0;
}
int densifyCsrL(float *a, int nrows, int ncols, int lda, int ldb, float *data, unsigned int *colind,unsigned long *indptr,unsigned int n){
  //put the sparse matrix into a.  For some reason, csr.todense() segments on a large (18GB) matrix.
  //For a c contiguous matrix, ldb==1, lda==ncols.
  long indx;
  long i,j;
  for(i=0; i<nrows; i++){
    for(j=indptr[i]; j<indptr[i+1]; j++){
      indx=(long)colind[j]*(long)ldb+(long)i*lda;
      a[indx]=data[j];
    }
  }
  return 0;
}
