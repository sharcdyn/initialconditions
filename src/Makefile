FLAGS           =
FFLAGS           =
CPPFLAGS         =
FPPFLAGS         =
CLEANFILES       =

### This can be defined from outside:
# make COMP=gfortran
UNAME = ${shell uname}

COMP=gfortran
STATIC=no

F77C= $(COMP)
F90C= $(COMP)
F90OPTS=-O2
F77OPTS=-O2
LINKOPTS=-I. -llapack

MKDIR_P=mkdir -p
BIN=../bin
.PHONY: ${BIN}
	

OBJS= diag-lapack.o

all: ${BIN} init init_rv xyz2geom

${BIN}:
	${MKDIR_P} ${BIN}

init: $(OBJS) sharc_init.o
	$(F90C) $(F90OPTS) -o ${BIN}/sharc_init.exe sharc_init.o $(OBJS) \
	$(LINKOPTS)

init_rv:	sharc_init_rv.o
	$(F90C) $(F90OPTS) -o ${BIN}/sharc_init_rv.exe sharc_init_rv.o $(OBJS) \
        $(LINKOPTS)

xyz2geom:	xyz2geom.o
	$(F90C) $(F90OPTS) -o ${BIN}/xyz2geom.exe xyz2geom.o $(OBJS) \
        $(LINKOPTS)

clean:
	rm -f *.o *.mod

.SUFFIXES : .f90 .f

.f90.o:
	$(F90C) $(F90OPTS) -c $<

.f.o:
	$(F77C) $(F77OPTS) -c $<

