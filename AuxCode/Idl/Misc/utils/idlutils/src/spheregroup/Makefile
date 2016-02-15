###############################################################################
# Sloan Digital Sky Survey (SDSS)
# IDL utilities for group finding and matching quickly
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
	spheregroup.o \
	spherematch.o \
	chunks.o \
	rarange.o \
	separation.o \
	friendsoffriends.o \
	chunkfriendsoffriends.o 

all : $(LIB)/libspheregroup.$(SO_EXT)


$(LIB)/libspheregroup.$(SO_EXT): $(OBJECTS)
	$(LD) $(X_LD_FLAGS) -o $(LIB)/libspheregroup.$(SO_EXT) $(OBJECTS)
#	nm -s ../lib/libspheregroup.$(SO_EXT)

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install : 
	-@ cp *.c $(IDLUTILS_DIR)/src/spheregroup
	-@ cp *.h $(IDLUTILS_DIR)/src/spheregroup
	-@ cp Makefile $(IDLUTILS_DIR)/src/spheregroup

clean :
	- /bin/rm -f *~ core *.o so_locations
