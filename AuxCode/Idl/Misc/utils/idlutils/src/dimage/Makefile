###############################################################################
# IDL utilities for image analysis
#
# M.R. Blanton
###############################################################################

SHELL = /bin/sh
#
.c.o :
	$(CC) -c $(CCCHK) $(CFLAGS) $(X_CFLAGS) $*.c
#
INC = ../../include
CFLAGS  = $(SDSS_CFLAGS) -DCHECK_LEAKS -I$(INC)

LIB = $(IDLUTILS_DIR)/lib

RECIPES = 
OBJECTS = \
	idl_dmedsmooth.o \
	dmedsmooth.o \
	dselip.o \
	dfind.o \
	idl_dsmooth.o \
	dsmooth.o \
	dsigma.o \
	idl_dobjects_multi.o \
	dobjects_multi.o 

all : $(LIB)/libdimage.$(SO_EXT)

$(LIB)/libdimage.$(SO_EXT): $(OBJECTS)
	$(LD) $(X_LD_FLAGS) -o $(LIB)/libdimage.$(SO_EXT) $(OBJECTS)
#	nm -s ../lib/libdimage.$(SO_EXT)

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install : 
	-@ cp *.c $(IDLUTILS_DIR)/src/dimage
	-@ cp *.h $(IDLUTILS_DIR)/src/dimage
	-@ cp Makefile $(IDLUTILS_DIR)/src/dimage

clean :
	- /bin/rm -f *~ core *.o so_locations
