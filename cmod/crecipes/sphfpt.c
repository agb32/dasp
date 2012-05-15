#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nrutil.h"
#define N1 2
#define N2 1
#define NTOT (N1+N2)
#define DXX 1.0e-4

int m,n;
float c2,dx,gmma;

int nn2,nvar;
float x1,x2,xf;

int main(void)  /* Program sphfpt */
{
	void newt(float x[], int n, int *check,
		void (*vecfunc)(int, float [], float []));
	void shootf(int n, float v[], float f[]);
	int check,i;
	float q1,*v1,*v2,*v;

	v=vector(1,NTOT);
	v1=v;
	v2 = &v[N2];
	nvar=NTOT;
	nn2=N2;
	dx=DXX;
	for (;;) {
		printf("input m,n,c-squared\n");
		if (scanf("%d %d %f",&m,&n,&c2) == EOF) break;
		if (n < m || m < 0) continue;
		gmma=1.0;
		q1=n;
		for (i=1;i<=m;i++) gmma *= -0.5*(n+i)*(q1--/i);
		v1[1]=n*(n+1)-m*(m+1)+c2/2.0;
		v2[2]=v1[1];
		v2[1]=gmma*(1.0-(v2[2]-c2)*dx/(2*(m+1)));
		x1 = -1.0+dx;
		x2=1.0-dx;
		xf=0.0;
		newt(v,NTOT,&check,shootf);
		if (check) {
			printf("shootf failed; bad initial guess\n");
		} else {
			printf("\tmu(m,n)\n");
			printf("%12.6f\n",v[1]);
		}
	}
	free_vector(v,1,NTOT);
	return 0;
}

void load1(float x1, float v1[], float y[])
{
	float y1 = (n-m & 1 ? -gmma : gmma);
	y[3]=v1[1];
	y[2] = -(y[3]-c2)*y1/(2*(m+1));
	y[1]=y1+y[2]*dx;
}

void load2(float x2, float v2[], float y[])
{
	y[3]=v2[2];
	y[1]=v2[1];
	y[2]=(y[3]-c2)*y[1]/(2*(m+1));
}

void score(float xf, float y[], float f[])
{
	int i;

	for (i=1;i<=3;i++) f[i]=y[i];
}
#undef N1
#undef N2
#undef NTOT
#undef DXX
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 21%. */
