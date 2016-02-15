###############################################################################
# Sloan Digital Sky Survey (SDSS)
# IDL support code for products: idlmapper, idlspec2d
#
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
all :
	idl mk_doc

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install :
	@ rm -rf $(IDLUTILS_DIR)/doc/www
	@ mkdir $(IDLUTILS_DIR)/doc/www
	- cp mk_doc $(IDLUTILS_DIR)/doc
	- cp www/*.html $(IDLUTILS_DIR)/doc/www
	- cp Makefile $(IDLUTILS_DIR)/doc

clean :
	- /bin/rm -f *~ core www/*~ www/core


