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
#
INC = ../../include
CFLAGS = $(SDSS_CFLAGS) -DCHECK_LEAKS -I$(INC)

LIB = $(IDLUTILS_DIR)/lib

RECIPES = nrecipes.o nrutil.o
# NNLS = 	nnls.o diff.o g1.o g2.o h12.o
OBJECTS = \
	arravsigclip.o \
	arravsigmask.o \
	arrmedian.o \
	ccorrelate.o \
	idl_mmsparse.o \
	idl_mmeval.o \
	memshift.o \
	$(RECIPES)

all : $(LIB)/libmath.$(SO_EXT)


$(LIB)/libmath.$(SO_EXT): $(OBJECTS)
	$(LD) $(X_LD_FLAGS) -o $(LIB)/libmath.$(SO_EXT) $(OBJECTS)
#	nm -s ../lib/libmath.$(SO_EXT)

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install : 
	-@ cp *.c $(IDLUTILS_DIR)/src/math
	-@ cp *.h $(IDLUTILS_DIR)/src/math
	-@ cp Makefile $(IDLUTILS_DIR)/src/math

clean :
	- /bin/rm -f *~ core *.o so_locations
