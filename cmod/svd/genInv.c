//stuff here is related to things within ritvec.
#include <stdio.h>
#include <math.h>
#include "sparsemem.h"

int convertToSparse(ArrUnion *genInv,double minGIVal, int transpose,ArrUnion *spGenInv){
  //Here, a dense matrix is converted to sparse format.  The size of
  //the dense matrix, genInv is taken from spGenInv and transpose.
  //The best way to do this is probabyl to have transpose set (and have it also set when you compute
  //the genInv in the first place.
  float fval;
  double dval;
  int rout,cout,cin;//rin not used - comment out - UB 2012Aug08
  long colstep,rowstep;
  long ndata,indx,ii;
  int i,j;
  int rtval=0;
  float minGIValf=(float)minGIVal;
  if((genInv->typ=='d' || genInv->typ=='f') && (spGenInv->typ=='D' || spGenInv->typ=='F')){
    if(spGenInv->typ=='F'){
      rout=spGenInv->data.fSp->rows;
      cout=spGenInv->data.fSp->cols;
      ndata=spGenInv->data.fSp->ndata;
    }else{
      rout=spGenInv->data.dSp->rows;
      cout=spGenInv->data.dSp->cols;
      ndata=spGenInv->data.dSp->ndata;
    }
    if(transpose){
      //      rin=cout; not used - comment out - UB 2012Aug08
      cin=rout;
      colstep=1;
      rowstep=(long)cin;
    }else{
      //      rin=rout; not used - comment out - UB 2012Aug08
      cin=cout;
      colstep=(long)cin;
      rowstep=1;
    }
    printf("Converting to sparse - ndata=%ld, rows=%d, cols=%d types %c %c rowstep %ld colstep %ld\n",ndata,rout,cout,genInv->typ,spGenInv->typ,rowstep,colstep);
    //Now, iterate over rows then columns...
    //We do rows first, because this is best for the sparse memory format.
    indx=0;
    if(spGenInv->typ=='D'){
      if(genInv->typ=='d'){//converting from dense double to sparse double...
	spGenInv->data.dSp->indptr[0]=0;
	for(i=0; i<cout; i++){
	  for(j=0; j<rout; j++){
	    ii=j*colstep+i*rowstep;
	    if(fabs(dval=genInv->data.ddata[ii])>minGIVal){//a non-zero value - insert it.
	      if(indx>=ndata){
		printf("ERROR: sparse memory full, ignoring rest of data %d %d %d %d.\n",i,j,cout,rout);
		j=rout;
		i=cout;
		rtval=1;
	      }else{
		spGenInv->data.dSp->data[indx]=dval;
		spGenInv->data.dSp->rowind[indx]=j;
		indx++;
	      }
	    }
	  }
	  spGenInv->data.dSp->indptr[i+1]=indx;
	}
      }else{//converting from dense float to sparse double...
	spGenInv->data.dSp->indptr[0]=0;
	for(i=0; i<cout; i++){
	  for(j=0; j<rout; j++){
	    ii=j*colstep+i*rowstep;
	    if(fabs(dval=(double)genInv->data.fdata[ii])>minGIVal){//a non-zero value - insert it.
	      if(indx>=ndata){
		printf("ERROR: sparse memory full, ignoring rest of data %d %d %d %d.\n",i,j,cout,rout);
		j=rout;
		i=cout;
		rtval=1;
	      }else{
		spGenInv->data.dSp->data[indx]=dval;
		spGenInv->data.dSp->rowind[indx]=j;
		indx++;
	      }
	    }
	  }
	  spGenInv->data.dSp->indptr[i+1]=indx;
	}
      }
    }else{
      if(genInv->typ=='d'){//converting from dense double to sparse float
	spGenInv->data.fSp->indptr[0]=0;
	for(i=0; i<cout; i++){
	  for(j=0; j<rout; j++){
	    ii=j*colstep+i*rowstep;
	    if(fabs(dval=genInv->data.ddata[ii])>minGIVal){//a non-zero value - insert it.
	      if(indx>=ndata){
		printf("ERROR: sparse memory full, ignoring rest of data %d %d %d %d.\n",i,j,cout,rout);
		j=rout;
		i=cout;
		rtval=1;
	      }else{
		spGenInv->data.fSp->data[indx]=(float)dval;
		spGenInv->data.fSp->rowind[indx]=j;
		indx++;
	      }
	    }
	  }
	  spGenInv->data.fSp->indptr[i+1]=indx;
	}
      }else{//converting from dense float to sparse float.
	spGenInv->data.fSp->indptr[0]=0;
	for(i=0; i<cout; i++){
	  //printf("%d\n",i);
	  for(j=0; j<rout; j++){
	    ii=j*colstep+i*rowstep;
	    //if(i==84894){
	    //  printf("%d %ld %ld %ld %ld\n",j,ii,colstep,rowstep,indx);
	    //}
	    if(fabsf(fval=genInv->data.fdata[ii])>minGIValf){//a non-zero value - insert it.
	      if(indx>=ndata){
		printf("ERROR: sparse memory full, ignoring rest of data i=%d j=%d cout=%d rout=%d ii=%ld colstep=%ld rowstep=%ld indx=%ld.\n",i,j,cout,rout,ii,colstep,rowstep,indx);
		j=rout;
		i=cout;
		rtval=1;
	      }else{
		spGenInv->data.fSp->data[indx]=fval;
		spGenInv->data.fSp->rowind[indx]=j;
		indx++;
	      }
	    }
	  }
	  spGenInv->data.fSp->indptr[i+1]=indx;
	}
      }
    }
  }else{
    printf("Error - wrong input/output array types\n");
    rtval=1;
  }
  return rtval;
}
