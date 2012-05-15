
/* JOSE analysis subroutines  */

#include <stdio.h>
#include <stdlib.h>
#include "jose.h"

/* allocate memory for a 2D integer array */

int **alloc2d_int(int m, int n)		/* m rows,  n columns  */
{
	int i;
	int **a;

	a=calloc(m, sizeof(int *));
	for(i=0;i<m;++i){
		a[i] = calloc(n, sizeof(int));
	}
	return a;
}


/* allocate memory for a 2D float array */

float **alloc2d_float(int m, int n)
{
	int i;
	float **a;

	a=calloc(m, sizeof(float *));
	for(i=0;i<m;++i){
		a[i] = calloc(n, sizeof(float));
	}
	return a;
}

int free2d_float(float **a,int m,int n){
    int i;
    for(i=0; i<m; i++){
	free(a[i]);
    }
    free(a[i]);
    return 0;
}

/* (Try to) open a file. Return file pointer */

FILE	*openfile(char *filename, char *mode)
/*      ********		*/
{
	FILE	*ifile;

	ifile=fopen(filename,mode);
	if(ifile==NULL)
	{
		printf("Can't open file:  %s\n",filename);
		exit(2);
	}
	else{
		printf("Opened file:  %s\n",filename);
	}
	return ifile;
}




