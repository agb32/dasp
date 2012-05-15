#include <stdio.h>
#define NRANSI
#include "nrutil.h"
#define N2 1

int m,n;
float c2,dx,gmma;

int nvar;
float x1,x2;

int main(void)  /* Program sphoot */
{
	void newt(float x[], int n, int *check,
		void (*vecfunc)(int, float [], float []));
	void shoot(int n, float v[], float f[]);
	int check,i;
	float q1,*v;

	v=vector(1,N2);
	dx=1.0e-4;
	nvar=3;
	for (;;) {
		printf("input m,n,c-squared\n");
		if (scanf("%d %d %f",&m,&n,&c2) == EOF) break;
		if (n < m || m < 0) continue;
		gmma=1.0;
		q1=n;
		for (i=1;i<=m;i++) gmma *= -0.5*(n+i)*(q1--/i);
		v[1]=n*(n+1)-m*(m+1)+c2/2.0;
		x1 = -1.0+dx;
		x2=0.0;
		newt(v,N2,&check,shoot);
		if (check) {
			printf("shoot failed; bad initial guess\n");
		} else {
			printf("\tmu(m,n)\n");
			printf("%12.6f\n",v[1]);
		}
	}
	free_vector(v,1,N2);
	return 0;
}

void load(float x1, float v[], float y[])
{
	float y1 = (n-m & 1 ? -gmma : gmma);
	y[3]=v[1];
	y[2] = -(y[3]-c2)*y1/(2*(m+1));
	y[1]=y1+y[2]*dx;
}

void score(float xf, float y[], float f[])
{
	f[1]=(n-m & 1 ? y[1] : y[2]);
}

void derivs(float x, float y[], float dydx[])
{
	dydx[1]=y[2];
	dydx[2]=(2.0*x*(m+1.0)*y[2]-(y[3]-c2*x*x)*y[1])/(1.0-x*x);
	dydx[3]=0.0;
}
#undef N2
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 21%. */
