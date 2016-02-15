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
	-@ cp Makefile $(IDLUTILS_DIR)/goddard/data
	-@ cp README $(IDLUTILS_DIR)/goddard/data
	-@ cp ebcasc.dat $(IDLUTILS_DIR)/goddard/data

clean :
	- /bin/rm -f *~ core
