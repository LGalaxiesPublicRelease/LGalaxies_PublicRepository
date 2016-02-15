###############################################################################
# LOCAL CONFIGURATION: EDIT AS DESIRED
###############################################################################

#------
# Linux

#CC = gcc
#CC = cc
#CFLAGS = -g -Wall -DDEBUG -DLINUX
CFLAGS = $(X_CFLAGS) -D$(UNAME)

# f2c
#F77 = f77
#FFLAGS = -u -a -C -g -DF2C
FFLAGS = $(X_FFLAGS) -DF2C
#FFLAGS = -u -O -DSUN

## gnu g77 (which assumes no SAVE)
#F77 = g77
#FFLAGS = -Wimplicit -g -O3 -DG77

all: mangler

##----
## Sun
#
#CC = cc
#CFLAGS = -O -DSUN
#
##CC = gcc
##CFLAGS = -O3 -DGCC
#
#F77 = f77
##FFLAGS = -u -O -pg -DSUN
#FFLAGS = -u -O -DSUN
#
### gnu g77 (which assumes no SAVE)
##F77 = g77
##FFLAGS = -Wimplicit -g -O3 -DG77
#
#all: mangler

#-----------------
# starlink library

# path to search for the starlink library
#ISLIB = -L$(IDLUTILS_DIR)/lib

# starlink library should be called libstar.a
# Edit as desired,
# or make a soft link from the actual starlink library to libstar.a:
#    ln -s <name_of_starlink_library> libstar.a
#SLIB = -lstar

###############################################################################
# END OF LOCAL CONFIGURATION
###############################################################################

ILIB = -L.
LLIB = -lmangle

PROGS = mrb_balkanize balkanize harmonize map polyid poly2poly ransack snap unify weight 
# excised rotate because of failed starlink dependency
#rotate

COBJ = advise_fmt.o convert.o copy_format.o copy_poly.o drandom.o findtop_.o garea.o gcmlim.o gptin.o gspher.o gvert.o gvlim.o gvphi.o ikrand.o msg.o new_poly.o new_vert.o partition_poly.o places.o prune_poly.o rdangle.o rdline.o rdmask.o rdspher.o scale.o snap_poly.o strdict.o vmid.o weight_fn.o wrangle.o wrho.o wrmask.o wrspher.o add_parent.o

FOBJ = azel.s.o azell.s.o felp.s.o fframe.s.o findtop.s.o garea.s.o gcmlim.s.o gphi.s.o gptin.s.o gsphera.s.o gspher.s.o gsubs.s.o gvert.s.o gvlim.s.o gvphi.s.o iylm.s.o twodf.o twoqz.o wlm.s.o wrho.s.o

# objects not already present (will try to extract them from libmangle.a)
OBJ := $(strip $(foreach objfile, $(COBJ) $(FOBJ), $(shell if [ ! -r $(objfile) ]; then echo $(objfile); fi)))

ifeq ($(CC), gcc)
	CFLAGS := $(CFLAGS) -DGCC
endif

mangle:
	@ if [ -r libmangle.a -a "$(OBJ)" ] ; then \
		echo "extracting files from libmangle.a" ; \
		ar xo libmangle.a $(OBJ) &> /dev/null ; \
		rm -f __.SYMDEF ; \
	fi
	@ make mangler
	@ rm -f $(COBJ) $(FOBJ)

mangler:
	make libmangle.a $(PROGS)
	cp $(PROGS) $(IDLUTILS_DIR)/bin
	cp libmangle.a $(IDLUTILS_DIR)/lib

clean:
	rm -f *.o core

cleaner:
	rm -f *.o core libmangle.a

cleanest:
	rm -f *.o core libmangle.a $(PROGS)

libmangle.a: $(COBJ) $(FOBJ) Makefile
	@ if [ -r libmangle.a ] ; then echo "updating libmangle.a" ; else echo "creating libmangle.a" ; fi
	ar ru libmangle.a $(COBJ) $(FOBJ)
	ranlib libmangle.a

balkanize: balkanize.o libmangle.a Makefile
	$(F77) $(FFLAGS) -o balkanize balkanize.o $(ILIB) $(LLIB)

mrb_balkanize: mrb_balkanize.o libmangle.a Makefile
	$(F77) $(FFLAGS) -o mrb_balkanize mrb_balkanize.o $(ILIB) $(LLIB)

harmonize: harmonize.o libmangle.a Makefile
	$(F77) $(FFLAGS) -o harmonize harmonize.o $(ILIB) $(LLIB)

map: map.o libmangle.a Makefile
	$(F77) $(FFLAGS) -o map map.o $(ILIB) $(LLIB)

polyid: polyid.o libmangle.a Makefile
	$(F77) $(FFLAGS) -o polyid polyid.o $(ILIB) $(LLIB)

poly2poly: poly2poly.o libmangle.a Makefile
	$(F77) $(FFLAGS) -o poly2poly poly2poly.o $(ILIB) $(LLIB)

ransack: ransack.o libmangle.a Makefile
	$(F77) $(FFLAGS) -o ransack ransack.o $(ILIB) $(LLIB)

rotate: rotate.o libmangle.a Makefile
	$(F77) $(FFLAGS) -o rotate rotate.o $(ILIB) $(LLIB) $(ISLIB) $(SLIB)

snap: snap.o libmangle.a Makefile
	$(F77) $(FFLAGS) -o snap snap.o $(ILIB) $(LLIB)

unify: unify.o libmangle.a Makefile
	$(F77) $(FFLAGS) -o unify unify.o $(ILIB) $(LLIB)

weight: weight.o libmangle.a Makefile
	$(F77) $(FFLAGS) -o weight weight.o $(ILIB) $(LLIB)

advise_fmt.o: angunit.h format.h scale.h advise_fmt.c
	$(CC) $(CFLAGS) -c advise_fmt.c
balkanize.o: parse_args.c defaults.h format.h polygon.h scale.h usage.h balkanize.c
	$(CC) $(CFLAGS) -c balkanize.c
mrb_balkanize.o: parse_args.c defaults.h format.h polygon.h scale.h usage.h mrb_balkanize.c
	$(CC) $(CFLAGS) -c mrb_balkanize.c
convert.o: pi.h polygon.h vertices.h convert.c
	$(CC) $(CFLAGS) -c convert.c
copy_format.o: format.h copy_format.c
	$(CC) $(CFLAGS) -c copy_format.c
copy_poly.o: polygon.h copy_poly.c
	$(CC) $(CFLAGS) -c copy_poly.c
drandom.o: drandom.c
	$(CC) $(CFLAGS) -c drandom.c
findtop_.o: findtop_.c
	$(CC) $(CFLAGS) -c findtop_.c
garea.o: logical.h polygon.h garea.c
	$(CC) $(CFLAGS) -c garea.c
gcmlim.o: polygon.h gcmlim.c
	$(CC) $(CFLAGS) -c gcmlim.c
gptin.o: logical.h polygon.h gptin.c
	$(CC) $(CFLAGS) -c gptin.c
gspher.o: harmonics.h logical.h polygon.h gspher.c
	$(CC) $(CFLAGS) -c gspher.c
gvert.o: logical.h polygon.h gvert.c
	$(CC) $(CFLAGS) -c gvert.c
add_parent.o: logical.h polygon.h add_parent.c
	$(CC) $(CFLAGS) -c add_parent.c
gvlim.o: logical.h polygon.h gvlim.c
	$(CC) $(CFLAGS) -c gvlim.c
gvphi.o: polygon.h gvphi.c
	$(CC) $(CFLAGS) -c gvphi.c
harmonize.o: parse_args.c defaults.h format.h harmonics.h polygon.h scale.h usage.h harmonize.c
	$(CC) $(CFLAGS) -c harmonize.c
ikrand.o: ikrand.c
	$(CC) $(CFLAGS) -c ikrand.c
map.o: parse_args.c angunit.h defaults.h format.h harmonics.h inputfile.h scale.h usage.h map.c
	$(CC) $(CFLAGS) -c map.c
msg.o: msg.c
	$(CC) $(CFLAGS) -c msg.c
new_poly.o: polygon.h new_poly.c
	$(CC) $(CFLAGS) -c new_poly.c
new_vert.o: vertices.h new_vert.c
	$(CC) $(CFLAGS) -c new_vert.c
partition_poly.o: polygon.h partition_poly.c
	$(CC) $(CFLAGS) -c partition_poly.c
places.o: places.c
	$(CC) $(CFLAGS) -c places.c
poly2poly.o: parse_args.c defaults.h format.h polygon.h scale.h usage.h poly2poly.c
	$(CC) $(CFLAGS) -c poly2poly.c
polyid.o: parse_args.c angunit.h defaults.h format.h inputfile.h polygon.h scale.h usage.h polyid.c
	$(CC) $(CFLAGS) -c polyid.c
prune_poly.o: polygon.h prune_poly.c
	$(CC) $(CFLAGS) -c prune_poly.c
ransack.o: parse_args.c angunit.h defaults.h format.h polygon.h scale.h usage.h vertices.h ransack.c
	$(CC) $(CFLAGS) -c ransack.c
rdangle.o: rdangle.c
	$(CC) $(CFLAGS) -c rdangle.c
rdline.o: inputfile.h rdline.c
	$(CC) $(CFLAGS) -c rdline.c
rdmask.o: format.h inputfile.h polygon.h scale.h vertices.h rdmask.c
	$(CC) $(CFLAGS) -c rdmask.c
rdspher.o: harmonics.h rdspher.c
	$(CC) $(CFLAGS) -c rdspher.c
rotate.o: parse_args.c parse_fopt.c angunit.h defaults.h format.h inputfile.h polygon.h scale.h usage.h rotate.c
	$(CC) $(CFLAGS) -c rotate.c
scale.o: scale.h vertices.h scale.c
	$(CC) $(CFLAGS) -c scale.c
snap.o: parse_args.c defaults.h format.h polygon.h scale.h usage.h snap.c
	$(CC) $(CFLAGS) -c snap.c
snap_poly.o: polygon.h snap_poly.c
	$(CC) $(CFLAGS) -c snap_poly.c
strdict.o: strdict.c
	$(CC) $(CFLAGS) -c strdict.c
unify.o: parse_args.c defaults.h format.h polygon.h scale.h usage.h unify.c
	$(CC) $(CFLAGS) -c unify.c
vmid.o: polygon.h vmid.c
	$(CC) $(CFLAGS) -c vmid.c
weight.o: parse_args.c defaults.h format.h polygon.h scale.h usage.h vertices.h weight.c
	$(CC) $(CFLAGS) -c weight.c
weight_fn.o: inputfile.h weight_fn.c
	$(CC) $(CFLAGS) -c weight_fn.c
wrangle.o: wrangle.c
	$(CC) $(CFLAGS) -c wrangle.c
wrho.o: harmonics.h wrho.c
	$(CC) $(CFLAGS) -c wrho.c
wrmask.o: format.h pi.h polygon.h vertices.h wrmask.c
	$(CC) $(CFLAGS) -c wrmask.c
wrspher.o: harmonics.h wrspher.c
	$(CC) $(CFLAGS) -c wrspher.c

azell.s.o: azell.s.f
	$(F77) $(FFLAGS) -c azell.s.f
azel.s.o: azel.s.f
	$(F77) $(FFLAGS) -c azel.s.f
felp.s.o: frames.par felp.s.f
	$(F77) $(FFLAGS) -c felp.s.f
fframe.s.o: frames.par radian.par fframe.s.f
	$(F77) $(FFLAGS) -c fframe.s.f
findtop.s.o: heapsort.inc findtop.s.f
	$(F77) $(FFLAGS) -c findtop.s.f
garea.s.o: pi.par garea.s.f
	$(F77) $(FFLAGS) -c garea.s.f
gcmlim.s.o: pi.par gcmlim.s.f
	$(F77) $(FFLAGS) -c gcmlim.s.f
gphi.s.o: pi.par gphi.s.f
	$(F77) $(FFLAGS) -c gphi.s.f
gptin.s.o: gptin.s.f
	$(F77) $(FFLAGS) -c gptin.s.f
gsphera.s.o: pi.par gsphera.s.f
	$(F77) $(FFLAGS) -c gsphera.s.f
gspher.s.o: pi.par gspher.s.f
	$(F77) $(FFLAGS) -c gspher.s.f
gsubs.s.o: pi.par gsubs.s.f
	$(F77) $(FFLAGS) -c gsubs.s.f
gvert.s.o: pi.par gvert.s.f
	$(F77) $(FFLAGS) -c gvert.s.f
gvlim.s.o: pi.par gvlim.s.f
	$(F77) $(FFLAGS) -c gvlim.s.f
gvphi.s.o: pi.par gvphi.s.f
	$(F77) $(FFLAGS) -c gvphi.s.f
iylm.s.o: pi.par iylm.s.f
	$(F77) $(FFLAGS) -c iylm.s.f
twodf.o: pi.par twodf.f
	$(F77) $(FFLAGS) -c twodf.f
twoqz.o: mangdir.data mangdir.inc twoqz.f
	$(F77) $(FFLAGS) -c twoqz.f
wlm.s.o: pi.par wlm.s.f
	$(F77) $(FFLAGS) -c wlm.s.f
wrho.s.o: pi.par wrho.s.f
	$(F77) $(FFLAGS) -c wrho.s.f

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install : 
	-@ cp *.c $(IDLUTILS_DIR)/src/mangle1.2
	-@ cp *.f $(IDLUTILS_DIR)/src/mangle1.2
	-@ cp *.h $(IDLUTILS_DIR)/src/mangle1.2
	-@ cp Makefile $(IDLUTILS_DIR)/src/mangle1.2

clean :
	- /bin/rm -f *~ core *.o so_locations
