###############################################################################
# Sloan Digital Sky Survey (SDSS)
# IDL support code for products: idlmapper, idlspec2d
#
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
SUBDIRS = database fits_bintable image misc sdas_table tv \
	astro disk_io fits_table jhuapl plot structure \
	astrom fits idlphot math sdas robust sockets

all :;

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install :
	@ for f in $(SUBDIRS); do \
		(mkdir $(IDLUTILS_DIR)/goddard/pro/$$f; \
		echo Copying $$f; \
		cp $$f/README $(IDLUTILS_DIR)/goddard/pro/$$f; \
		cp $$f/*.pro $(IDLUTILS_DIR)/goddard/pro/$$f ); \
	done
	-@ cp Makefile $(IDLUTILS_DIR)/goddard/pro
	-@ cp README $(IDLUTILS_DIR)/goddard/pro

clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; /bin/rm -f core *~ ); \
	done
