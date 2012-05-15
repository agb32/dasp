# Linux or Windows:
CC = gcc -Wall -O4# -march=i486
# CC = icc -w1 -O3 -march=i486
CFLAGS = -fPIC
# Macintosh:
ifeq ($(HOSTTYPE),powerpc)
  CC = cc -pipe -O3 -Wall -fno-common -arch ppc
endif

LIBS=-lm
OBJ=svdlib.o svdutil.o las2.o sparsemem.o ritvecfff.o ritvecffd.o ritvecfdf.o ritvecfdd.o ritvecdff.o ritvecdfd.o ritvecddf.o ritvecddd.o genInv.o

#svd: Makefile main.o libsvd.a
#	${CC} ${CFLAGS} -o svd main.o libsvd.a ${LIBS}
#	mv -f $@ ${HOSTTYPE}/$@
#	ln -s ${HOSTTYPE}/$@ $@

#main.o: Makefile main.c svdlib.h
#	${CC} ${CFLAGS} -c main.c
install: svdmodule.so
	echo Made SVD module

svdmodule.so: libsvd.a svdmodule.c setup.py
	rm -f ../svdmodule.so
	rm -rf build
	python setup.py install --install-lib=../

libsvd.a: Makefile ${OBJ}
	rm -f $@
	ar cr $@ ${OBJ}
	ranlib $@
#libsvd.a: ${HOSTTYPE} ${OBJ}
#	rm -f $@ ${HOSTTYPE}/$@
#	ar cr $@ ${OBJ}
#	ranlib $@
#	mv -f $@ ${HOSTTYPE}/$@
#	ln -s ${HOSTTYPE}/$@ $@

svdlib.o: Makefile svdlib.h svdlib.c
	${CC} ${CFLAGS} -c svdlib.c
svdutil.o: Makefile svdutil.c svdutil.h
	${CC} ${CFLAGS} -c svdutil.c
las2.o: Makefile las2.c svdlib.h svdutil.h
	${CC} ${CFLAGS} -c las2.c
sparsemem.o: Makefile sparsemem.c sparsemem.h
	${CC} ${CFLAGS} -c sparsemem.c
ritvec: ritvecfff.o ritvecffd.o ritvecfdf.o ritvecfdd.o ritvecdff.o ritvecdfd.o ritvecddf.o ritvecddd.o
	echo Made ritvec

ritvecfff.o: Makefile ritvec.c ritvec.h svdlib.h svdutil.h sparsemem.h
	${CC} ${CFLAGS} -DSFLOAT -DSTOREFLOAT -DVTFLOAT -c ritvec.c -o ritvecfff.o
ritvecffd.o: Makefile ritvec.c ritvec.h svdlib.h svdutil.h sparsemem.h
	${CC} ${CFLAGS} -DSFLOAT -DSTOREFLOAT -DVTDOUBLE -c ritvec.c -o ritvecffd.o
ritvecfdf.o: Makefile ritvec.c ritvec.h svdlib.h svdutil.h sparsemem.h
	${CC} ${CFLAGS} -DSFLOAT -DSTOREDOUBLE -DVTFLOAT -c ritvec.c -o ritvecfdf.o
ritvecfdd.o: Makefile ritvec.c ritvec.h svdlib.h svdutil.h sparsemem.h
	${CC} ${CFLAGS} -DSFLOAT -DSTOREDOUBLE -DVTDOUBLE -c ritvec.c -o ritvecfdd.o
ritvecdff.o: Makefile ritvec.c ritvec.h svdlib.h svdutil.h sparsemem.h
	${CC} ${CFLAGS} -DSDOUBLE -DSTOREFLOAT -DVTFLOAT -c ritvec.c -o ritvecdff.o
ritvecdfd.o: Makefile ritvec.c ritvec.h svdlib.h svdutil.h sparsemem.h
	${CC} ${CFLAGS} -DSDOUBLE -DSTOREFLOAT -DVTDOUBLE -c ritvec.c -o ritvecdfd.o
ritvecddf.o: Makefile ritvec.c ritvec.h svdlib.h svdutil.h sparsemem.h
	${CC} ${CFLAGS} -DSDOUBLE -DSTOREDOUBLE -DVTFLOAT -c ritvec.c -o ritvecddf.o
ritvecddd.o: Makefile ritvec.c ritvec.h svdlib.h svdutil.h sparsemem.h
	${CC} ${CFLAGS} -DSDOUBLE -DSTOREDOUBLE -DVTDOUBLE -c ritvec.c -o ritvecddd.o
genInv.o: Makefile genInv.h genInv.c sparsemem.h
	${CC} ${CFLAGS} -c genInv.c -o genInv.o


clean: 
	rm *.o
	rm libsvd.a
	rm -rf build
#clean: 
#	rm *.o
#	rm libsvd.a
#	rm svd
#	rm -rf ${HOSTTYPE}

$(HOSTTYPE):
	if test ! -d $(HOSTTYPE); \
	then mkdir $(HOSTTYPE); fi
