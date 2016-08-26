/*
#dasp, the Durham Adaptive optics Simulation Platform.
#Copyright (C) 2004-2016 Alastair Basden and Durham University.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as
#published by the Free Software Foundation, either version 3 of the
#License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//acmlfft.c: C interface to FFT routines in ACMl
//used to link with g2c library
//written I think by FA...

#include <acml.h>
/*
//ACML info functions
void ACMLVERSION(int *major, int *minor, int *patch)
{
    acmlversion(major,minor,patch);
}

void ACMLINFO(void)
{
    acmlinfo();
}

//ACML FFT functions
//Single complex functions
void CFFT1D(int mode, int n, complex *x, complex *comm, int *info)
{
    cfft1d(mode,n,x,comm,info);
}
*/
extern void cfft1d(int mode, int n, complex *x, complex *comm, int *info);
extern void cfft1dx(int mode, float scale, int inpl, int n, complex *x, int incx, complex *y, int incy, complex *comm, int *info);
extern void cfft1m(int mode, int nseq, int n, complex *x, complex *comm, int *info);
extern void cfft1mx(int mode, float scale, int inpl, int nseq, int n, complex *x, int incx1, int incx2, complex *y, int incy1, int incy2, complex *comm, int *info);
extern void cfft2d(int mode, int m, int n, complex *x, complex *comm, int *info);
extern void cfft2dx(int mode, float scale, int ltrans, int inpl, int m, int n, complex *x, int incx1, int incx2, complex *y, int incy1, int incy2, complex *comm, int *info);
extern void cfft3d(int mode, int l, int m, int n, complex *x, complex *comm, int *info);
extern void cfft3dx(int mode, float scale, int ltrans, int inpl, int l, int m, int n, complex *x, complex *y, complex *comm, int *info);
extern void csfft(int mode, int n, float *x, float *comm, int *info);
extern void csfftm(int nseq, int n, float *x, float *comm, int *info);
extern void dzfft(int mode, int n, double *x, double *comm, int *info);
extern void dzfftm(int nseq, int n, double *x, double *comm, int *info);
extern void scfft(int mode, int n, float *x, float *comm, int *info);
extern void scfftm(int nseq, int n, float *x, float *comm, int *info);
extern void zdfft(int mode, int n, double *x, double *comm, int *info);
extern void zdfftm(int nseq, int n, double *x, double *comm, int *info);
extern void zfft1d(int mode, int n, doublecomplex *x, doublecomplex *comm, int *info);
extern void zfft1dx(int mode, double scale, int inpl, int n, doublecomplex *x, int incx, doublecomplex *y, int incy, doublecomplex *comm, int *info);
extern void zfft1m(int mode, int nseq, int n, doublecomplex *x, doublecomplex *comm, int *info);
extern void zfft1mx(int mode, double scale, int inpl, int nseq, int n, doublecomplex *x, int incx1, int incx2, doublecomplex *y, int incy1, int incy2, doublecomplex *comm, int *info);
extern void zfft2d(int mode, int m, int n, doublecomplex *x, doublecomplex *comm, int *info);
extern void zfft2dx(int mode, double scale, int ltrans, int inpl, int m, int n, doublecomplex *x, int incx1, int incx2, doublecomplex *y, int incy1, int incy2, doublecomplex *comm, int *info);
extern void zfft3d(int mode, int l, int m, int n, doublecomplex *x, doublecomplex *comm, int *info);
extern void zfft3dx(int mode, double scale, int ltrans, int inpl, int l, int m, int n, doublecomplex *x, doublecomplex *y, doublecomplex *comm, int *info);



/*

Instructions for compilations

gcc -m64 acmlfft.c -fPIC -shared -lacml -lg2c -o acmlfft.so


*/

