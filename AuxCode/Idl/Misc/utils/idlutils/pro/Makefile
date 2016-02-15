###############################################################################
# Sloan Digital Sky Survey (SDSS)
# IDL support code for products: idlmapper, idlspec2d
#
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
SUBDIRS = 2mass astrom bspline coord cosmography djsphot dust first fits healpix image mangle math misc mcmc mpeg mpfit physics plot psconfig psf rgbcolor sdss slatec spheregroup TeXtoIDL trace yanny

all :;

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install :
	@ /bin/rm -rf $(IDLUTILS_DIR)/pro
	@ mkdir $(IDLUTILS_DIR)/pro
	@ for f in $(SUBDIRS); do \
		(mkdir $(IDLUTILS_DIR)/pro/$$f; \
		echo Copying $$f; \
                cp $$f/*.pro $(IDLUTILS_DIR)/pro/$$f ); \
	done
	-@ cp Makefile $(IDLUTILS_DIR)/pro

clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; /bin/rm -f core *~ ); \
	done
