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
#install :
#	/bin/cp Makefile $(IDLUTILS_DIR)/bin
#	/bin/cp build_links $(IDLUTILS_DIR)/bin
#	/bin/cp evilinstall $(IDLUTILS_DIR)/bin
#	/bin/cp evilmake $(IDLUTILS_DIR)/bin
#	/bin/cp eviltarball $(IDLUTILS_DIR)/bin
#	/bin/cp install_ups_table $(IDLUTILS_DIR)/bin
#	/bin/cp memtrack $(IDLUTILS_DIR)/bin
#	/bin/cp pstogif $(IDLUTILS_DIR)/bin

EXECS1 = Makefile build_links evilmake eviltarball memtrack pstogif install_ups_table
EXECS2 = balkanize drangle harmonize map poly2poly polyid ransack snap unify
EXECS3 = mrb_balkanize mrb_snap mrb_unify weight

EXECS = $(EXECS1) $(EXECS2) $(EXECS3)

install :
	@ for f in $(EXECS); do \
		(/bin/cp $$f $(IDLUTILS_DIR)/bin/ ); \
	done


clean :
	- /bin/rm -f *~ core
