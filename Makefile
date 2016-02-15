EXEC   = L-Galaxies

OBJS   = ./code/main.o ./code/io_tree.o ./code/init.o ./code/cool_func.o \
     ./code/save.o ./code/save_galtree.o \
     ./code/mymalloc.o ./code/read_parameters.o \
	 ./code/peano.o ./code/allvars.o ./code/age.o ./code/update_type_two.o \
	 ./code/metals.o \
	 ./code/model_infall.o \
	 ./code/model_cooling.o \
	 ./code/model_starformation_and_feedback.o \
	 ./code/model_reincorporation.o \
	 ./code/model_mergers.o \
	 ./code/model_dust.o \
	 ./code/model_misc.o \
	 ./code/model_disrupt.o \
	 ./code/model_stripping.o \
	 ./code/scale_cosmology.o

INCL   = ./code/allvars.h  ./code/proto.h  Makefile


# Either include the default set of Makefile options, or define your own
#include My_Makefile_options
#include My_Makefile_options_MCMC
include My_Makefile_options_MCMC_HaloModel

# Choose your system type (needs to match an entry in Makefile_compilers)
SYSTYPE = "MyMachine"
include Makefile_compilers




LIBS   =   -g $(LDFLAGS) -lm  $(GSL_LIBS)  $(RLIBS) -lgsl -lgslcblas 

CFLAGS =   -g $(OPTIONS) $(OPT) -DCOMPILETIMESETTINGS=\""$(OPT)"\" $(OPTIMIZE) $(GSL_INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

#$(OBJS): $(INCL) My_Makefile_options Makefile_compilers
#$(OBJS): $(INCL) My_Makefile_options_MCMC Makefile_compilers
$(OBJS): $(INCL) My_Makefile_options_MCMC_haloModel Makefile_compilers

clean:
	rm -f $(OBJS)

tidy:
	rm -f $(OBJS) .$(EXEC)

# use next target to generate metadata about the result files
# uses -E compiler option to preprocess the allvars.h file, stores result in allvars.i
# then calls awk scripts from ./awk/ folder to extract cleand-up version of GALAXY_OUTPUT struct
# and generate different representations of use for post-processing the result 	
metadata:
	${CC_MD} ${OPT} ${CFLAGS} -E ./code/allvars.h -o ./code/allvars.i
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_TypeString.awk > ./AuxCode/awk/L-Galaxies_Types.txt
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_DDL.awk > ./AuxCode/awk/L-Galaxies_DDL.sql	
ifeq (NORMALIZEDDB,$(findstring NORMALIZEDDB,$(OPT)))
	awk -f ./AuxCode/awk/extractSFH_BIN.awk ./code/allvars.i |awk -f ./AuxCode/awk/SFH_BIN_2_DDL.awk >> ./AuxCode/awk/L-Galaxies_DDL.sql
else
	awk -f ./AuxCode/awk/extractSFH_Time.awk ./code/allvars.i |awk -f ./AuxCode/awk/SFH_Time_2_DDL.awk >> ./AuxCode/awk/L-Galaxies_DDL.sql
endif	
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/idl/GALAXY_OUTPUT_2_IDL_struct.awk >  ./AuxCode/awk/idl/LGalaxy.pro
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/idl/GALAXY_OUTPUT_2_IDL_hists.awk > ./AuxCode/awk/idl/LGalaxy_plot.pro
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/idl/GALAXY_OUTPUT_2_IDL_testfloats.awk > ./AuxCode/awk/idl/LGalaxy_testfloats.pro
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/idl/GALAXY_OUTPUT_2_IDL_zerofloats.awk > ./AuxCode/awk/idl/LGalaxy_zerofloats.pro
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_LGalaxy.awk > ./AuxCode/awk/L-Galaxies.h
	awk -f ./AuxCode/awk/extractGALAXY_OUTPUT.awk ./code/allvars.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_FileFormat.awk > ./AuxCode/awk/L-Galaxies_FileFormat.csv
	awk -f ./AuxCode/awk/extractSFH_BIN.awk ./code/allvars.i |awk -f ./AuxCode/awk/MOMAF_INPUT_2_MoMaFGalaxy.awk >> ./AuxCode/awk/L-Galaxies.h

metadata_db:
	awk -f ./AuxCode/awk/extract_struct_metals.awk ./code/allvars.i > ./AuxCode/awk/structs.dat
	awk -f ./AuxCode/awk/extract_struct_elements.awk ./code/allvars.i >> ./AuxCode/awk/structs.dat
	awk -f ./AuxCode/awk/extract_struct_GALAXY_OUTPUT.awk ./code/allvars.i >> ./AuxCode/awk/structs.dat
	
