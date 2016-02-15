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
CFLAGS  = $(SDSS_CFLAGS) -DCHECK_LEAKS -I$(INC)

LIB = $(IDLUTILS_DIR)/lib

OBJECTS = \
	evilmath.o \
	trace_crude.o \
	trace_fweight.o \
	trace_gweight.o

all : $(LIB)/libtrace.$(SO_EXT)


$(LIB)/libtrace.$(SO_EXT): $(OBJECTS)
	$(LD) $(X_LD_FLAGS) -o $(LIB)/libtrace.$(SO_EXT) $(OBJECTS)
#	nm -s ../lib/libtrace.$(SO_EXT)

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install : 
	-@ cp *.c $(IDLUTILS_DIR)/src/trace
	-@ cp *.h $(IDLUTILS_DIR)/src/trace
	-@ cp Makefile $(IDLUTILS_DIR)/src/trace

clean :
	- /bin/rm -f *~ core *.o so_locations
