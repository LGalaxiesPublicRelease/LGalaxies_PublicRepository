###############################################################################
# Sloan Digital Sky Survey (SDSS)
# IDL support code for products: idlmapper, idlspec2d
#
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
SUBDIRS = data pro text

all :;

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install :
	- cp Makefile $(IDLUTILS_DIR)/goddard
	- cp README $(IDLUTILS_DIR)/goddard
	- cp *.txt $(IDLUTILS_DIR)/goddard
	@ for f in $(SUBDIRS); do \
		mkdir $(IDLUTILS_DIR)/goddard/$$f; \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) install ); \
	done


clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; /bin/rm -f core *~ ); \
	done
