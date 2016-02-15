###############################################################################
# Sloan Digital Sky Survey (SDSS)
# IDL support code for products: idlmapper, idlspec2d
#
# Slatec code in Fortran with C wrappers for calling from IDL
#
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
.c.o :
	$(CC) -c $(CCCHK) $(CFLAGS) $(X_CFLAGS) $*.c

.f.o :
	$(F77) -c $(FFLAGS) $(X_CFLAGS) $*.f

#
INC = ../../include
CFLAGS  = -I$(IDL_DIR)/external $(SDSS_CFLAGS) -DCHECK_LEAKS -I$(INC) 

LIB = $(IDLUTILS_DIR)/lib

RECIPES = nrecipes.o nrutil.o
OBJECTS = \
	bint4.o bndacc.o bndsol.o bnfac.o bnslv.o bsplvn.o bspvd.o bspvn.o \
	bvalu.o chkder.o csevl.o d1mach_Linux.o dasum.o daxpy.o dbint4.o dbintk.o dbndac.o \
	dbndsl.o dbnfac.o dbnslv.o dbspvd.o dbspvn.o dbvalu.o dckder.o dcopy.o \
	dcsevl.o ddot.o defc.o defcmn.o denorm.o derf.o derfc.o dfdjc3.o dfspvn.o \
	dgeco.o dgefa.o dgefs.o dgesl.o dh12.o dintrv.o dmpar.o dnls1.o dnls1e.o \
	dp1vlu.o dpcoef.o dpolft.o dpsort.o dqrfac.o dqrslv.o dscal.o dsidi.o \
	dsifa.o dsort.o dswap.o dwupdt.o efc.o efcmn.o enorm.o erf.o erfc.o fdjac3.o \
	fdump.o h12.o i1mach_Linux.o idamax.o idlwrappers.o initds.o inits.o intrv.o \
	isort.o j4save.o lmpar.o pcoef.o polfit.o pvalue.o qrfac.o qrsolv.o  \
	r1mach_Linux.o radb2.o radb3.o radb4.o radb5.o radbg.o radf2.o radf3.o \
	radf4.o radf5.o radfg.o rand.o rfftb.o rfftb1.o rfftf.o rfftf1.o rffti.o \
	rffti1.o rgauss.o runif.o rwupdt.o saxpy.o scopy.o sdot.o snls1.o snls1e.o \
	sscal.o ssort.o sswap.o xercnt.o xerhlt.o xermsg.o xerprn.o xersve.o \
	xgetua.o xsetf.o 

all : $(LIB)/libslatec.$(SO_EXT)

$(LIB)/libslatec.$(SO_EXT) : $(OBJECTS)
	$(LD) $(X_LD_FLAGS) -o $(LIB)/libslatec.$(SO_EXT) $(OBJECTS)

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install :
	- cp Makefile $(IDLUTILS_DIR)/src/slatec
	- cp *.f $(IDLUTILS_DIR)/src/slatec
	- cp *.c $(IDLUTILS_DIR)/src/slatec

clean :
	- /bin/rm -f *~ core *.o so_locations
