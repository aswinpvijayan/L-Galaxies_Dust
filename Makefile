#  Copyright (C) <2016>  <L-Galaxies>

#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>

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
include My_Makefile_options
#include My_Makefile_options_MCMC
#include My_Makefile_options_MCMC_HaloModel

# Choose your system type (needs to match an entry in Makefile_compilers)
<<<<<<< HEAD
#SYSTYPE = "MyMachine"
include Makefile_compilers
#include My_Makefile_compilers
=======
SYSTYPE = "MyMachine"
#include Makefile_compilers
include My_Makefile_compilers
>>>>>>> 39853857269bc0dcfaeb394b679ec6c80393ee4b

LIBS   =   -g $(LDFLAGS) -lm  $(GSL_LIBS)  $(RLIBS) -lgsl -lgslcblas 
ifeq (HDF5_OUTPUT,$(findstring HDF5_OUTPUT,$(OPT)))
LIBS += $(HDF5_LIBS) -lhdf5 -lhdf5_hl
endif

CFLAGS =   -g $(OPTIONS) $(OPT) -DCOMPILETIMESETTINGS=\""$(OPT)"\" $(OPTIMIZE) $(GSL_INCL) $(HDF5_INCL)
ifeq (HDF5_OUTPUT,$(findstring HDF5_OUTPUT,$(OPT)))
CFLAGS += $(HDF5_INCL)
endif

all: metadata $(EXEC)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) My_Makefile_options Makefile_compilers
#$(OBJS): $(INCL) My_Makefile_options_MCMC Makefile_compilers
#$(OBJS): $(INCL) My_Makefile_options_MCMC_haloModel Makefile_compilers

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
	${CC_MD} ${OPT} ${CFLAGS} -E -CC ./code/h_galaxy_output.h -o ./code/h_galaxy_output.i
	${CC_MD} ${OPT} ${CFLAGS} -E -CC ./code/h_metals.h -o ./code/h_metals.i
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_TypeString.awk          > ./AuxCode/awk/output/L-Galaxies_Types.txt
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_DDL.awk                 > ./AuxCode/awk/output/L-Galaxies_DDL.sql	
ifeq (NORMALIZEDDB,$(findstring NORMALIZEDDB,$(OPT)))
	awk -f ./AuxCode/awk/extract_SFH_BIN.awk       ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/SFH_BIN_2_DDL.awk                       >> ./AuxCode/awk/output/L-Galaxies_DDL.sql
else
	awk -f ./AuxCode/awk/extract_SFH_Time.awk      ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/SFH_Time_2_DDL.awk                      >> ./AuxCode/awk/output/L-Galaxies_DDL.sql
endif	
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_IDL_struct.awk          > ./AuxCode/awk/output/idl/LGalaxy.pro
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_IDL_hists.awk           > ./AuxCode/awk/output/idl/LGalaxy_plot.pro
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_IDL_testfloats.awk      > ./AuxCode/awk/output/idl/LGalaxy_testfloats.pro
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_IDL_zerofloats.awk      > ./AuxCode/awk/output/idl/LGalaxy_zerofloats.pro
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_LGalaxy.awk             > ./AuxCode/awk/output/L-Galaxies.h
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_FileFormat.awk          > ./AuxCode/awk/output/L-Galaxies_FileFormat.csv
	awk -f ./AuxCode/awk/extract_SFH_BIN.awk       ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/MOMAF_INPUT_2_MoMaFGalaxy.awk           >> ./AuxCode/awk/output/L-Galaxies.h
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_2_python_struct.awk       > ./AuxCode/awk/output/python/LGalaxy.py
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_3_HDF5.awk                > ./code/io_hdf5.h
	awk -f ./AuxCode/awk/extract_GALAXY_OUTPUT_props.awk ./code/h_galaxy_output.i |awk -f ./AuxCode/awk/GALAXY_OUTPUT_prop_2_HDF5_proptable.awk > ./input/hdf5_field_props.txt

	awk -f ./AuxCode/awk/extract_struct_metals.awk ./code/h_metals.i                                                                      > ./AuxCode/awk/output/structs.dat
	awk -f ./AuxCode/awk/extract_struct_elements.awk ./code/h_metals.i                                                                    >> ./AuxCode/awk/output/structs.dat
	awk -f ./AuxCode/awk/extract_struct_GALAXY_OUTPUT.awk ./code/h_galaxy_output.i                                                        >> ./AuxCode/awk/output/structs.dat
