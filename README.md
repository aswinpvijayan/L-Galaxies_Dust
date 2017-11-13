# L-Galaxies Dust

This repo holds the implementation of the detailed dust model into the public released version of the Henriques2015 L-Galaxies model. This adds a model of dust production from AGB stars, supernovae as well as grain growth in molecular clouds. It also includes a model of dust destruction from supernovae shocks. 


## Installing and Running the model

A step by step series of instructions to run this version of L-Galaxies on your system. 

* ```git clone https://github.com/scottclay/Lgalaxies_Dust.git```

You will probably need to change the filepaths or obtain the following required treefiles etc. from elsewhere:

* Input file 
	* FirstFile
	* LastFile
	* CoolFunctionsDir
	* SpecPhotDir
	* SimulationDir
	* OutputDir
	
* The following compiler options have been added to ```My_Makefile_options```
	* ```OPT += -DDETAILED_DUST``` - Switch on dust model with default options
	* ```OPT += -DDUST_AGB``` - Switch on AGB dust production (default on)
	* ```OPT += -DDUST_SNII``` - Switch on SNII dust production (default on)
	* ```OPT += -DDUST_SNIA``` - Switch on SNIA dust production (default on)
	* ```OPT += -DDUST_GROWTH``` - Switch on grain growth dust production (default on)
	* ```OPT += -DDUST_DESTRUCTION``` - Switch on dust destruction (default on)
	* ```OPT += -DFULL_DUST_RATES``` - Output dust production/destruction rates (default on)
	* ```OPT += -DREDUCED_OUTPUT``` - Shorten output (mainly cutting SFH bins)
	* ```OBJS += ./code/model_dustyields.o ```
	* ```OBJS += ./code/dustyields_functions.o```
	* ```OBJS += ./code/dustyields_integrals.o```
	* ```OBJS += ./code/dustyields_read_tables.o```

To compile the model:
```
make
```

To run the model locally on box 5 of MR   ```./L-Galaxies input.1```

To run the model locally on box 40 of MRII ```./L-Galaxies input.1```

To run the model on apollo for all MR at Sussex ```qsub batch_MR_apollo.sh``` 

To run the model on apollo for box 1 and 2 of MRII at Sussex ```qsub -l m_mem_free=220G batch_MRII_apollo1.sh``` 

To run the model on apollo for box 3-512 of MRII at Sussex ```qsub batch_MRII_apollo2.sh``` 


 