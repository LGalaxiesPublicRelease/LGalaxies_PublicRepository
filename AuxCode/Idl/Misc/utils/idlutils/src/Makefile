###############################################################################
# Sloan Digital Sky Survey (SDSS)
# IDL support code for products: idlmapper, idlspec2d
#
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
.c.o :
	$(CC) -c $(CCCHK) $(CFLAGS) $*.c
#
CFLAGS  = $(SDSS_CFLAGS) -DCHECK_LEAKS -I../include

# Note that either mangle1.4/src or mangle1.2 must be built before idlmangle!
SUBDIRS = image mangle1.4/src idlmangle math spheregroup trace dimage

all :
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) all ); \
	done

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install :
	@ for f in $(SUBDIRS); do \
		(mkdir $(IDLUTILS_DIR)/src/$$f; cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) install ); \
	done
	-@ cp Makefile $(IDLUTILS_DIR)/src

clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) clean ); \
	done
