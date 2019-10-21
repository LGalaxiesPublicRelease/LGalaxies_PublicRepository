EXEC  = L-Galaxies

# Default object files (others may be added with -D options)
OBJS  = ./code/main.o \
	./code/allvars.o \
	./code/age.o \
	./code/cool_func.o \
	./code/io_tree.o \
	./code/init.o \
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
	./code/mymalloc.o \
	./code/peano.o \
	./code/read_parameters.o \
	./code/save.o \
	./code/save_galtree.o \
	./code/scale_cosmology.o \
	./code/update_type_two.o

# The following is used only to set dependencies
INCL  = Makefile \
	./code/allvars.h \
	./code/h_funcs.h \
	./code/h_galaxy_output.h \
	./code/h_galaxy_tree_data.h \
	./code/h_galaxy.h \
	./code/h_halo_data.h \
	./code/h_halo_ids_data.h \
	./code/h_halo_aux_data.h \
	./code/h_lightcone.h \
	./code/h_metals.h \
	./code/h_params.h \
	./code/h_variables.h \
	./code/proto.h
ifeq (ALL_SKY_LIGHTCONE,$(findstring ALL_SKY_LIGHTCONE,$(OPT)))
INCL  += ./code/lightcone.h
endif

# Either include the default set of Makefile options, or define your own
# include Makefile_options
include My_Makefile_options
#include My_Makefile_options_MCMC

# Choose your system type (copy an entry from Makefile_compilers)
include My_Makefile_compilers

LIBS   =   -g $(LDFLAGS) -lm  $(GSL_LIBS)  $(RLIBS) -lgsl -lgslcblas $(HDF5_LIBS) -lhdf5_serial -lhdf5_serial_hl
#LIBS   =   -g $(LDFLAGS) -lm  $(GSL_LIBS)  $(RLIBS) -lgsl -lgslcblas $(HDF5_LIBS) -lhdf5 -lhdf5_hl

CFLAGS =   -g $(OPTIONS) $(OPT) -DCOMPILETIMESETTINGS=\""$(OPT)"\" $(OPTIMIZE) $(GSL_INCL) $(HDF5_INCL)

all: metadata $(EXEC)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) Makefile My_Makefile_compilers My_Makefile_options

clean:
	rm -f $(OBJS)

tidy:
	rm -f $(OBJS) .$(EXEC)

# use next target to generate metadata about the result files
# uses -E compiler option to preprocess the allvars.h file, stores result in allvars.i
# uses -CC compiler option to save comments, needed for HDF5 output
# then calls awk scripts from ./awk/ folder to extract cleand-up version of GALAXY_OUTPUT struct
# and generate different representations of use for post-processing the result 	
metadata:
	${CC} ${OPT} ${CFLAGS} -E -CC ./code/h_galaxy_output.h -o ./code/h_galaxy_output.i
	${CC} ${OPT} ${CFLAGS} -E -CC ./code/h_metals.h -o ./code/h_metals.i
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_TypeString.awk > ./AuxCode/awk/output/L-Galaxies_Types.txt
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_DDL.awk > ./AuxCode/awk/output/L-Galaxies_DDL.sql	
ifeq (NORMALIZEDDB,$(findstring NORMALIZEDDB,$(OPT)))
	awk -f ./AuxCode/awk/extract_SFH_BIN.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/SFH_BIN_2_DDL.awk >> ./AuxCode/awk/output/L-Galaxies_DDL.sql
else
	awk -f ./AuxCode/awk/extract_SFH_Time.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/SFH_Time_2_DDL.awk >> ./AuxCode/awk/output/L-Galaxies_DDL.sql
endif	
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_IDL_struct.awk >  ./AuxCode/awk/output/idl/LGalaxy.pro
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_IDL_hists.awk > ./AuxCode/awk/output/idl/LGalaxy_plot.pro
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_IDL_testfloats.awk > ./AuxCode/awk/output/idl/LGalaxy_testfloats.pro
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_IDL_zerofloats.awk > ./AuxCode/awk/output/idl/LGalaxy_zerofloats.pro
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_LGalaxy.awk > ./AuxCode/awk/output/L-Galaxies.h
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_FileFormat.awk > ./AuxCode/awk/output/L-Galaxies_FileFormat.csv
	awk -f ./AuxCode/awk/extract_SFH_BIN.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/MOMAF_INPUT_2_MoMaFGalaxy.awk >> ./AuxCode/awk/output/L-Galaxies.h
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_python_struct.awk >  ./AuxCode/awk/output/python/LGalaxy.py
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_HDF5.awk > ./code/io_hdf5.h
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT_props.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_prop_2_HDF5_proptable.awk > ./input/hdf5_field_props.txt

	awk -f ./AuxCode/awk/extract_struct_metals.awk ./code/h_metals.i > ./AuxCode/awk/output/structs.dat
	awk -f ./AuxCode/awk/extract_struct_elements.awk ./code/h_metals.i >> ./AuxCode/awk/output/structs.dat
	awk -f ./AuxCode/awk/extract_struct_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i >> ./AuxCode/awk/output/structs.dat
