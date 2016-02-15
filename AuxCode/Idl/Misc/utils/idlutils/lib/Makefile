###############################################################################
# Sloan Digital Sky Survey (SDSS)
# IDL support code for products: idlmapper, idlspec2d
#
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
all :;

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install :
	- cp *.$(SO_EXT) $(IDLUTILS_DIR)/lib
	- cp *.a $(IDLUTILS_DIR)/lib
	- cp Makefile $(IDLUTILS_DIR)/lib

clean :
	- /bin/rm -f *~ core *.so *.a *.dylib

