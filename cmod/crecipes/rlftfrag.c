#include <stdlib.h>
#define NRANSI
#include "nrutil.h"
#define N2 256
#define N3 256

int main(void)		/* example1 */
{
	void rlft3(float ***data, float **speq, unsigned long nn1,
		unsigned long nn2, unsigned long nn3, int isign);
	float ***data, **speq;

	data=f3tensor(1,1,1,N2,1,N3);
	speq=matrix(1,1,1,2*N2);
	rlft3(data,speq,1,N2,N3,1);
	rlft3(data,speq,1,N2,N3,-1);
	free_matrix(speq,1,1,1,2*N2);
	free_f3tensor(data,1,1,1,N2,1,N3);
	return 0;
}
#undef N3
#undef N2

#define N1 32
#define N2 64
#define N3 16

int main(void)		/* example2 */
{
	void rlft3(float ***data, float **speq, unsigned long nn1,
		unsigned long nn2, unsigned long nn3, int isign);
	int j;
	float ***data,**speq;

	data=f3tensor(1,N1,1,N2,1,N3);
	speq=matrix(1,N1,1,2*N2);
	rlft3(data,speq,N1,N2,N3,1);
	free_matrix(speq,1,N1,1,2*N2);
	free_f3tensor(data,1,N1,1,N2,1,N3);
	return 0;
}
#undef N1
#undef N2
#undef N3

#define N 32

int main(void)		/* example3 */
{
	void rlft3(float ***data, float **speq, unsigned long nn1,
		unsigned long nn2, unsigned long nn3, int isign);
	int j;
	float fac,r,i,***data1,***data2,**speq1,**speq2,*sp1,*sp2;

	data1=f3tensor(1,N,1,N,1,N);
	data2=f3tensor(1,N,1,N,1,N);
	speq1=matrix(1,N,1,2*N);
	speq2=matrix(1,N,1,2*N);

	rlft3(data1,speq1,N,N,N,1);
	rlft3(data2,speq2,N,N,N,1);
	fac=2.0/(N*N*N);
	sp1 = &data1[1][1][1];
	sp2 = &data2[1][1][1];
	for (j=1;j<=N*N*N/2;j++) {
		r = sp1[0]*sp2[0] - sp1[1]*sp2[1];
		i = sp1[0]*sp2[1] + sp1[1]*sp2[0];
		sp1[0] = fac*r;
		sp1[1] = fac*i;
		sp1 += 2;
		sp2 += 2;
	}
	sp1 = &speq1[1][1];
	sp2 = &speq2[1][1];
	for (j=1;j<=N*N;j++) {
		r = sp1[0]*sp2[0] - sp1[1]*sp2[1];
		i = sp1[0]*sp2[1] + sp1[1]*sp2[0];
		sp1[0] = fac*r;
		sp1[1] = fac*i;
		sp1 += 2;
		sp2 += 2;
	}
	rlft3(data1,speq1,N,N,N,-1);
	free_matrix(speq2,1,N,1,2*N);
	free_matrix(speq1,1,N,1,2*N);
	free_f3tensor(data2,1,N,1,N,1,N);
	free_f3tensor(data1,1,N,1,N,1,N);
	return 0;
}
#undef N
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 21%. */
