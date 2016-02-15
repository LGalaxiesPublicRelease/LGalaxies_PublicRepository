###############################################################################
# Sloan Digital Sky Survey (SDSS)
# IDL support code for products: idlmapper, idlspec2d
#
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
SUBDIRS = filters objects sdss

all :;

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install :
	@ /bin/rm -rf $(IDLUTILS_DIR)/data
	@ mkdir $(IDLUTILS_DIR)/data
	@ for f in $(SUBDIRS); do \
		(mkdir $(IDLUTILS_DIR)/data/$$f; \
		echo Copying $$f; \
                cp $$f/* $(IDLUTILS_DIR)/data/$$f ); \
	done
	-@ cp Makefile $(IDLUTILS_DIR)/data

clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; /bin/rm -f core *~ ); \
	done
