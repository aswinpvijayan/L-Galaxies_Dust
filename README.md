# L-Galaxies Dust

This repo holds the implementation of the detailed dust model into the public released version of the Henriques2015 L-Galaxies model. This adds a model of dust production from AGB stars, supernovae as well as grain growth in molecular clouds. It also includes a model of dust destruction from supernovae shocks. The dust production arising from AGB and supernovae implemented following dust production tables from Zhukovska et al. 2008. This is implemented thorugh look up tables that injects dust based on the metallicity. All of the recipes for dust injection and destruction by SNe as well as integrating the code to the already present ```L-Galaxies``` framework were done by Dr. Scott Clay. The changes to his model is in the way dustgrowth and destruction in the ISM are implemented. Documentation of the public version of the code can be found at http://galformod.mpa-garching.mpg.de/public/LGalaxies/

## Changes to Henriques 2015 you might expect

* Detailed dust modelling. Each galaxy now outputs the dust mass as a function of constituent elements (i.e. C,O,Si etc.). You should sum these for the total dust mass. 
* The dust model is called in ```main.c``` (```update_dust_mass```) immediately AFTER the call to the chemical enrichment model (```update_yields_and_return_mass```)
* Changes to various models, transfer functions, header files etc. to account for the dust. 
* New file: ```dustyields_read_tables.c``` - Similar to yields_read_tables but for dust. Run before L-Galaxies. Reads in the mass and metallicity dependent dust yield tables. 
* New file: ```dustyields_integrals.c``` - Similar to yields_integrals.c but for dust. Run before L-Galaxies. Pre-calculates the normalised ejecta rates at every timestep, assuming 1 Msun populations. 
* New file: ```model_dustyields.c``` - Holds the new detailed dust modelling code.
* New functions to calculate the molecular gas fraction following three implementations: Blitz & Rosolowsky 2006, Gnedin & Kratsow 2011 and Krumholz et al. 2013. This can be implemented by switching on ```Obreshkow``` or ```GK11``` or ```Krumholz``` in  My_Makefile_options. Based on the output from these functions the cold ISM is split into moelcular clouds and the diffused ISM. The dust is injected into the two using this ratio. Dust growth happens only in molecular clouds, with exchange between the two media taking place on an effective exchange timescale.



## Changes to Henriques 2015 you might not expect

* If dust is switched on, the call to ```SN_feedback``` that was at the end of model_yields is now called AFTER the call to the dust model in ```main.c```. 
* My_Makefile_options has a new option, ```REDUCED_OUTPUT```, which is switched on by default. The dust model requires SFHs to be switched on, but the SFH output when combined with the chemical enrichment model and the full dust model is very large. In order to make it easier to play with, as well as to take up less hard drive space, ```REDUCED_OUTPUT``` cuts the SFH from the output. 


## Installing and Running the model

A step by step series of instructions to run this version of L-Galaxies on your system. 

You will probably need to change the filepaths or obtain treefiles/coolfunctions etc. from elsewhere (check out http://galformod.mpa-garching.mpg.de/public/LGalaxies/running_the_model.php for general instructions on running the Henriques 2015 version of L-Galaxies, as well as where to find the tree files and Stellar Population Synthesis tables):

* Edit the following in the input file:
	* FirstFile
	* LastFile
	* CoolFunctionsDir
	* SpecPhotDir
	* SimulationDir
	* OutputDir
	
* The following compiler options have been added to ```My_Makefile_options```

    * ```OPT += -DObreshkow``` - Calculate the molecular gas fraction using Blitz & Rosolowsky 2006 
    * ```OPT += -DGK11``` - Calculate the molecular gas fraction using Gnedin & Kratsow 2011
    * ```OPT += -DKrumholz``` - Calculate the molecular gas fraction using Fu et al. 2013
	* ```OPT += -DDETAILED_DUST``` - Switch on dust model with default options
	* ```OPT += -DDUST_AGB``` - Switch on AGB dust production (default on)
	* ```OPT += -DDUST_SNII``` - Switch on SNII dust production (default on)
	* ```OPT += -DDUST_SNIA``` - Switch on SNIA dust production (default on)
	* ```OPT += -DDUST_GROWTH``` - Switch on grain growth dust production (default on)
	* ```OPT += -DDUST_TEXCH``` - Set the dust exchange timescale (set to 20.0E6 yrs now)
    * ```OPT += -DDUST_TACC0``` - Set the accretion timescale in molecular clouds (set to 0.5E5 yrs now)
	* ```OPT += -DDUST_DESTRUCTION``` - Switch on dust destruction (default on)
	* ```OPT += -DFULL_DUST_RATES``` - Output dust production/destruction rates (default on)
	* ```OPT += -DREDUCED_OUTPUT``` - Shorten output (mainly cutting SFH bins)
	* ```OBJS += ./code/model_dustyields.o ```
	* ```OBJS += ./code/dustyields_integrals.o```
	* ```OBJS += ./code/dustyields_read_tables.o```

To compile the model ```make```

To run the model locally on box 5 of MR   ```./L-Galaxies ./input/input_MR_W1_PLANCK_local.par```

To run the model locally on box 40 of MRII ```./L-Galaxies ./input/input_MRII_W1_PLANCK_local.par```

To run the model on apollo for all MR at Sussex ```qsub batch_MR_script.sh``` 

To run the model on apollo for box 1 and 2 of MRII at Sussex ```qsub -l m_mem_free=220G batch_MRII_script_1_2.sh``` 

To run the model on apollo for box 3-512 of MRII at Sussex ```qsub batch_MRII_script_3_512.sh``` 


