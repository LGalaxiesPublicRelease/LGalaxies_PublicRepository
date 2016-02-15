###############################################################################
# Sloan Digital Sky Survey (SDSS)
# IDL support code for products: idlmapper, idlspec2d
#
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
.c.o :
	$(CC) -c $(CCCHK) $(CFLAGS) $(X_CFLAGS) $*.c
#.f.o :
#	$(F77) -c $(CCCHK) $(CFLAGS) $(X_CFLAGS) $*.f
#
INC = ../../include
CFLAGS  = $(SDSS_CFLAGS) -DCHECK_LEAKS -I$(INC)

LIB = $(IDLUTILS_DIR)/lib

OBJECTS = \
	grow_obj.o \
	pop_image.o \
	idl_photfrac.o \
	idl_reject_cr_psf.o \
	photfrac.o \
	reject_cr_psf.o \
	interp_profmean.o \
	p_cisi.o \
	p_midpnt.o \
	p_qromo.o \
	p_polint.o \
	p_utils.o 
#	idl_interp_profmean.o
#	ppvalu.o
#	tautsp.o
#	interv.o

all : $(LIB)/libimage.$(SO_EXT)

$(LIB)/libimage.$(SO_EXT): $(OBJECTS) ph.h
	$(LD) $(X_LD_FLAGS) -o $(LIB)/libimage.$(SO_EXT) $(OBJECTS) -lm
#	$(LD) $(X_LD_FLAGS) -o $(LIB)/libimage.$(SO_EXT) $(OBJECTS) $(MAKE_FTNLIB)
#	nm -s ../lib/libimage.$(SO_EXT)

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install : 
	-@ cp *.c $(IDLUTILS_DIR)/src/image
	-@ cp *.h $(IDLUTILS_DIR)/src/image
	-@ cp Makefile $(IDLUTILS_DIR)/src/image

clean :
	- /bin/rm -f *~ core *.o so_locations
