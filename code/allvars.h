/*  Copyright (C) <2016>  <L-Galaxies>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/> */

/** @file allvars.h
  */
#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#define MIN_ALLOC_NUMBER       1000
#define ALLOC_INCREASE_FACTOR  1.1
#define ALLOC_DECREASE_FACTOR  0.7

#define PRECISION_LIMIT 1.e-7

#define  min(x,y)  ((x)<(y) ?(x):(y))
#define  max(x,y)  ((x)>(y) ?(x):(y))
#define  wrap(x,y) ( (x)>((y)/2.) ? ((x)-(y)) : ((x)<(-(y)/2.)?((x)+(y)):(x)) )
#define  pow2(x)   ((x)*(x))
#define  pow3(x)   ((x)*(x)*(x))

#define  terminate(x) {char termbuf[5000]; sprintf(termbuf, "code termination on task=%d, function %s(), file %s, line %d: %s\n", ThisTask, __FUNCTION__, __FILE__, __LINE__, x); printf("%s", termbuf); fflush(stdout); endrun(1);}


#define  mymalloc(x, y)            mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  mymalloc_movable(x, y, z) mymalloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__)

#define  myrealloc(x, y)           myrealloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  myrealloc_movable(x, y)   myrealloc_movable_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)

#define  myfree(x)                 myfree_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)
#define  myfree_movable(x)         myfree_movable_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)

#define  report_memory_usage(x, y) report_detailed_memory_usage_of_largest_task(x, y, __FUNCTION__, __FILE__, __LINE__)


#ifdef GALAXYTREE
#define  CORRECTDBFLOAT(x)  ((fabs(x)<(1.e-30) || isnan(x)) ?(0.0):(x))
#else
#define  CORRECTDBFLOAT(x) x
#endif


//WATCH OUT! In the case of MCMC running both MR and MRII the larger value is used to "allocate" all the arrays
//inside the code its LastDarkMatterSnapShot+1 that defines the extent of the loops
//(in MCMC MR_plus_MRII mode this are not always identical)
#ifdef MRII
#define  MAXSNAPS  68     /* Number of snapshots in the dark matter simulation */
#else

#ifdef PHOENIX
#define  MAXSNAPS  72
#else

#ifdef CATERPILLAR
#define  MAXSNAPS  256
#else

#define  MAXSNAPS  64  //NORMAL MILLENNIUM

#endif //CATERPILLAR
#endif //PHOENIX
#endif //MRII

#define  MAXGALFAC 2.3 /*1.5/2.3 - maximum fraction of satellite without a halo (for memory allocation)*/

#define  STEPS 20		/* Number of integration intervals between two snapshots */

#define  ALLOCPARAMETER 50.  /* new definition !!! THIS HAS TO BE 50 !!! DONT EVER EVER EVER CHANGE !!! */

//To understand the units in the code read through set_units in init.c!!!
#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  PROTONMASS  1.6726e-24
#define  HUBBLE      3.2407789e-18   /* in h/sec */

//To understand the units in the code read through set_units in init.c!!!

#define UNITLENGTH_IN_CM                   3.08568e+24	// Mpc - WATCH OUT, distances in the code are in Mpc/h
#define UNITMASS_IN_G                      1.989e+43	// 10^10Msun - WATCH OUT, masses in the code are in 10^10Msun/h
#define UNITVELOCITY_IN_CM_PER_S           100000	// Km/s - WATCH OUT, this are the correct units in the code km/s
#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#ifdef GALAXYTREE
#undef  NOUT
#define NOUT MAXSNAPS
#endif


#ifdef STAR_FORMATION_HISTORY
#define SFH_NMERGE 3  //  SFH_NMERGE=Nmax+1 (Nmax used in Shamshiri2014)

#ifdef CATERPILLAR
#define SFH_NBIN 24 //  CATERPILLAR - 256 snapshots
#else
#define SFH_NBIN 20
#endif //CATERPILLAR

#endif //STAR_FORMATION_HISTORY

#include "h_metals.h"
#include "h_galaxy_output.h"

struct galaxy_tree_data
{
  int HaloGalIndex;
  int IndexStored;
  int SnapNum;
  int GalID;
  int FirstProgGal;
  int NextProgGal;
  int LastProgGal;
  int DescendantGal;
  int MainLeaf;
  int TreeRoot;
  int FOFCentralGal;
  int Done;
}
 *GalTree;

/*Structure with all the data associated with galaxies (this is not the same as the output!)*/
struct GALAXY			/* Galaxy data */
{
  int HeapIndex;
  int GalTreeIndex;
  int NextGalaxy;
#ifdef GALAXYTREE
  int FirstProgGal;
#endif
  int Type;
  int HaloNr;
  long long MostBoundID;
  int SnapNum;
  int CentralGal;  //own ID for types 0 and 1, unless 1's have little dark matter and already merging to type 0. For 2's its the merger centre
  float CentralMvir;
  /* properties of subhalo at the last time this galaxy was a central galaxy */
  float Pos[3];
  float MergCentralPos[3];
  float Vel[3];
  float Pos_notupdated[3];
  float Vel_notupdated[3];
#ifdef HALOPROPERTIES
  float HaloM_Mean200, HaloM_Crit200, HaloM_TopHat;
  float HaloPos[3];
  float HaloVel[3];
  float HaloVelDisp;
  float HaloVmax;
#endif
  float HaloSpin[3];
  float GasSpin[3];
  float StellarSpin[3];
  int   Len;   
  float Mvir;
  float Rvir;
  float Vvir;
  float Vmax;
  float InfallVmax;
  float InfallVmaxPeak; // km/s - Max previous Vmax at infall
  int InfallSnap;
  float InfallHotGas;
  float InfallHotGasRadius;
  float HotRadius;
  /* baryonic reservoirs */
  float ColdGas;
  float BulgeMass;
  float DiskMass;
  float HotGas;
  float EjectedMass;
  float BlackHoleMass;
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals MetalsColdGas;
  struct metals MetalsBulgeMass;
  struct metals MetalsDiskMass;
  struct metals MetalsHotGas;
  struct metals MetalsEjectedMass;
#ifdef METALS_SELF
  struct metals MetalsHotGasSelf;
#endif
#else
  float MetalsColdGas;
  float MetalsBulgeMass;
  float MetalsDiskMass;
  float MetalsHotGas;
  float MetalsEjectedMass;
#ifdef METALS_SELF
  float MetalsHotGasSelf;
#endif
#endif //DETAILED_METALS_AND_MASS_RETURN
#ifdef TRACK_BURST
  float BurstMass;
#endif

  /* misc */
  float PrimordialAccretionRate;
  float CoolingRate;
  float CoolingRate_beforeAGN;
  float CoolingRadius;
  float CoolingGas;
  float QuasarAccretionRate;
  float RadioAccretionRate;
  float AGNheatingFromCentral;
  float Sfr;
  float SfrBulge;
  float StarMerge;
  float XrayLum;
  float BulgeSize;
  float StellarDiskRadius;
  float GasDiskRadius;
#ifdef GALAXYTREE
  int   DisruptOn;
#endif
  // float halfradius;
  //float periradius;
  float CosInclination; //angle between galaxy spin and the z-axis
  float OriMergTime;
  float MergTime;
  float OriMvir;
  float OriRvir;
  float MergeSat;
  float DistanceToCentralGal[3];
  int MergeOn;
  float ICM;
#ifdef DETAILED_METALS_AND_MASS_RETURN
   struct metals MetalsICM;
 #else
   float MetalsICM;
 #endif

  /* luminosities in various bands */
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifndef POST_PROCESS_MAGS
#ifdef OUTPUT_REST_MAGS
  float Lum[NMAG][NOUT];
  float YLum[NMAG][NOUT];
  float LumBulge[NMAG][NOUT];
  float YLumBulge[NMAG][NOUT];
  float LumDust[NMAG][NOUT];
#ifdef ICL
  float ICLLum[NMAG][NOUT];
#endif
#endif //OUTPUT_REST_MAGS

#ifdef COMPUTE_OBS_MAGS
  float ObsLum[NMAG][NOUT];
  float ObsYLum[NMAG][NOUT];
  float ObsLumBulge[NMAG][NOUT];
  float ObsYLumBulge[NMAG][NOUT];
  float ObsLumDust[NMAG][NOUT];
#ifdef ICL
  float ObsICL[NMAG][NOUT];
#endif

#ifdef OUTPUT_MOMAF_INPUTS
  float dObsLum[NMAG][NOUT];
  float dObsYLum[NMAG][NOUT];
  float dObsLumBulge[NMAG][NOUT];
  float dObsYLumBulge[NMAG][NOUT];
  float dObsLumDust[NMAG][NOUT];
#ifdef ICL
  float dObsICL[NMAG][NOUT];
#endif
#endif
#endif //COMPUTE_OBS_MAGS

#endif //ndef POST_PROCESS_MAGS
#endif //COMPUTE_SPECPHOT_PROPERTIES

  float MassWeightAge[NOUT];
#ifdef STAR_FORMATION_HISTORY
  int sfh_ibin; //Index of highest bin currently in use
  double sfh_age; //Time in years of last call to sph_update_bins
  int sfh_flag[SFH_NBIN];
  float sfh_dt[SFH_NBIN]; //Size of time interval in units of years
  float sfh_t[SFH_NBIN]; //Time at low-redshift edge of bin in same units
  int sfh_Nbins[SFH_NBIN]; //Number of bins on the time interval
  float sfh_DiskMass[SFH_NBIN]; //Stellar mass in disk, in bin in standard units
  float sfh_BulgeMass[SFH_NBIN]; //Stellar mass in bulge, in bin in standard units
  float sfh_ICM[SFH_NBIN]; //Stellar mass in ICM, in bin in standard units
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals sfh_MetalsDiskMass[SFH_NBIN]; //Metals locked up in stars in disk.
  struct metals sfh_MetalsBulgeMass[SFH_NBIN]; //Metals locked up in stars in bulge.
  struct metals sfh_MetalsICM[SFH_NBIN]; //Metals locked up in stars in ICM.
#else
  float sfh_MetalsDiskMass[SFH_NBIN]; //Metals locked up in stars in disk.
  float sfh_MetalsBulgeMass[SFH_NBIN]; //Metals locked up in stars in bulge.
  float sfh_MetalsICM[SFH_NBIN]; //Metals locked up in stars in ICM.
#endif
#ifdef TRACK_BURST
  float sfh_BurstMass[SFH_NBIN]; //Stellar mass formed in bursts, in standard units.
#endif
#endif //STAR_FORMATION_HISTORY

#ifdef INDIVIDUAL_ELEMENTS
  struct elements sfh_ElementsDiskMass[SFH_NBIN];
  struct elements sfh_ElementsBulgeMass[SFH_NBIN];
  struct elements sfh_ElementsICM[SFH_NBIN];

  struct elements DiskMass_elements;
  struct elements BulgeMass_elements;
  struct elements ColdGas_elements;
  struct elements HotGas_elements;
  struct elements ICM_elements;
  struct elements EjectedMass_elements;
  struct elements ColdGasDiff_elements;
  struct elements ColdGasClouds_elements;
  float mu_gas;
#endif //INDIVIDUAL_ELEMENTS

#ifdef DETAILED_DUST
#ifdef FULL_DUST_RATES
	struct DustRates DustColdGasRates;		
#endif
	struct elements DustColdGasDiff_elements; // ? // Mass of elements in the diffused phase locked up in dust (in ColdGas)
    struct elements DustColdGasClouds_elements; // ? // Mass of elements in the clouds locked up in dust (in ColdGas)
    float t_acc; // ? // dust growth rate
    float f_i[9]; // ? // fraction of dust in diffused media
    float f_c[9]; // ? // fraction of dust in molecular media
    float f_cmax[9]; // ? // max fraction of molecular media
#endif //DETAILED_DUST

} *Gal, *HaloGal;


// Documentation can be found in the database
struct halo_data
{
	/* merger tree pointers */
	int Descendant;
	int FirstProgenitor;
	int NextProgenitor;
	int FirstHaloInFOFgroup;
	int NextHaloInFOFgroup;

  /* properties of halo */
	int Len;
	float M_Mean200, M_Crit200, M_TopHat;
	float Pos[3];
	float Vel[3];
	float VelDisp;
	float Vmax;
	float Spin[3];
	long long MostBoundID;

  /* original position in subfind output */
	int SnapNum;
	int FileNr;
	int SubhaloIndex;
	float SubHalfMass;
}
  *Halo, *Halo_Data;


// Documentation can be found in the database
#ifndef MCMC
extern struct  halo_ids_data
{
 long long HaloID;
 long long FileTreeNr;
 long long FirstProgenitor;
 long long LastProgenitor;
 long long NextProgenitor;
 long long Descendant;
 long long FirstHaloInFOFgroup;
 long long NextHaloInFOFgroup;
#ifdef MRII
 long long MainLeafID; 
#endif
 double    Redshift;
 int       PeanoKey;
 int       dummy;      /* need to use this padding for 64bit alignment */
} *HaloIDs, *HaloIDs_Data;
#else
extern struct  halo_ids_data
{
 long long FirstHaloInFOFgroup;
} *HaloIDs, *HaloIDs_Data;
#endif


// Documentation can be found in the database
struct halo_aux_data  /* auxiliary halo data */
{
	int DoneFlag;
	int HaloFlag;
	int NGalaxies;
	int FirstGalaxy;
	float M_Mean200_Unscaled;
	float M_Crit200_Unscaled;
	float Pos_Unscaled[3];
	float Vel_Unscaled[3];
	float Vmax_Unscaled;
	float Spin_Unscaled[3];
}
 *HaloAux;

// For some unkown reason need NOT to do extern on the following
char *inputFile; // Name of the input file

extern int FirstFile;		/* first and last file for processing */
extern int LastFile;

extern int Ntrees;		/* number of trees in current file */
extern double AllocValue_MaxHaloGal;
extern double AllocValue_MaxGal;
extern double AllocValue_MaxGalTree;

extern int MaxGal;		/* Maximum number of galaxies allowed for Gal[] array */
extern int NHaloGal, MaxHaloGal;
extern int NGalTree, MaxGalTree;
extern int *HaloGalHeap;
extern int IndexStored;

extern int LastSnapShotNr;

extern int LastDarkMatterSnapShot;
#ifdef MR_PLUS_MRII //OPTION for MCMC
extern int LastDarkMatterSnapShot_MR;
extern int LastDarkMatterSnapShot_MRII;
#endif


extern char SpecPhotDir[512];
extern char PhotPrefix[50];
extern char SpecPhotIMF[50];
extern char McFile[512];
extern char FileWithFilterNames[512];
extern char CoolFunctionsDir[512];
extern char OutputDir[512];
/* in case a second parameter is given as argument to the code, this will be taken as a
 * temporary outputdir to allow fast I/O. OutputDir will be replaced by this directory
 * and in the end everything will be moved to the FinalOutputDir (original OutputDir
 * given in input.par )*/
extern char FinalOutputDir[512];
extern char FileNameGalaxies[512];
extern char SimulationDir[512];
extern char FileWithOutputRedshifts[512];

extern char FileWithZList[512];
//variable used to scale to a different cosmology
extern char FileWithZList_OriginalCosm[512];
#ifdef MR_PLUS_MRII  //OPTION for MCMC
extern char FileWithZList_MR[512];
extern char FileWithZList_OriginalCosm_MR[512];
extern char FileWithZList_MRII[512];
extern char FileWithZList_OriginalCosm_MRII[512];
#endif

extern double ScalePos;
extern double ScaleMass;

#ifdef SPECIFYFILENR
extern char   FileNrDir[512];
extern int    ListInputFilrNr[111];
#endif

extern int TotHalos;
extern int TotGalaxies[NOUT];
extern int *TreeNgals[NOUT];

extern int *FirstHaloInSnap;

extern int *TreeNHalos;
extern int *TreeFirstHalo;

extern void *TreeAuxData;


extern double MaxMemSize;

extern size_t AllocatedBytes;
extern size_t HighMarkBytes;
extern size_t FreeBytes;

extern int ThisTask, NTask;

#ifdef GALAXYTREE
extern int GalCount;
extern int TotGalCount;
#endif

/* Cosmological parameters */
extern double BaryonFrac;
extern double Sigma8;
extern double Omega;
extern double OmegaLambda;
extern double Hubble_h;
extern double Omega_OriginalCosm;
extern double OmegaLambda_OriginalCosm;
extern double Hubble_h_OriginalCosm;
//SIMULATION RELATED
extern double PartMass;
extern double BoxSize;
extern double PartMass_OriginalCosm;
extern double BoxSize_OriginalCosm;
#ifdef MR_PLUS_MRII  //OPTION for MCMC
extern double PartMass_MR;
extern double BoxSize_MR;
extern double PartMass_OriginalCosm_MR;
extern double BoxSize_OriginalCosm_MR;
extern double PartMass_MRII;
extern double BoxSize_MRII;
extern double PartMass_OriginalCosm_MRII;
extern double BoxSize_OriginalCosm_MRII;
#endif


/* flags */
extern int ReionizationModel;
extern int DiskRadiusModel;
extern int StarFormationModel;
extern int FeedbackReheatingModel;
extern int FeedbackEjectionModel;
extern int FateOfSatellitesGas;
extern int ReIncorporationModel;
extern int AGNRadioModeModel;
extern int DiskInstabilityModel;
extern int BHGrowthInDiskInstabilityModel;
extern int HotGasStripingModel;
extern int DisruptionModel;
extern int StarBurstModel;
extern int BulgeFormationInMinorMergersOn;
extern int MetallicityOption;

/* parameters */
extern double Reionization_z0;
extern double Reionization_zr;
extern double Yield;
extern double RecycleFraction;
extern double ThreshMajorMerger;
extern double MergerTimeMultiplier;
extern double RamPressureStrip_CutOffMass;
extern double SfrEfficiency;
extern double SfrColdCrit;
extern double SfrBurstEfficiency;
extern double SfrBurstSlope;
extern double AgnEfficiency;
extern double BlackHoleGrowthRate;
extern double BlackHoleSeedMass;
extern double BlackHoleCutoffVelocity;
extern double FeedbackReheatingEpsilon;
extern double ReheatPreVelocity;
extern double ReheatSlope;
extern double FeedbackEjectionEfficiency;
extern double EjectPreVelocity;
extern double EjectSlope;
extern double ReIncorporationFactor;
extern double EnergySNcode, EnergySN;
extern double EtaSNcode, EtaSN;

extern double
	UnitTime_in_s,
	UnitPressure_in_cgs,
	UnitDensity_in_cgs,
	UnitCoolingRate_in_cgs,
	UnitEnergy_in_cgs,
	UnitTime_in_Megayears, //Using time as stored in the code, this gives Myr/h
	UnitTime_in_years,
#ifdef HALOMODEL
	RhoCrit,
#endif
	G,
	Hubble,
	a0, ar;

extern int ListOutputSnaps[NOUT];
extern float ListOutputRedshifts[NOUT];

extern double ZZ[MAXSNAPS];
extern double AA[MAXSNAPS];
//variable used to scale to a different cosmology
extern double AA_OriginalCosm[MAXSNAPS];

extern double Age[MAXSNAPS];

extern int    Zlistlen;

extern gsl_rng *random_generator;


extern int    NumMergers;


/*  tabulated stuff */

#ifdef STAR_FORMATION_HISTORY
/* SFH_ is the reference structure for storing the star formation histories in
 * logarithmic bins. It is computed in init.c generating a binning structure for
 * each snapshot/time step. In the code galaxy structures are adjusted with respect
 * to this structure at each step. */
extern double SFH_t[MAXSNAPS][STEPS][SFH_NBIN]; //Time to present at the lower edge of the bin (code units)
extern double SFH_dt[MAXSNAPS][STEPS][SFH_NBIN]; //Time width of the bin (code units)
extern int SFH_Nbins[MAXSNAPS][STEPS][SFH_NBIN]; //Number of bins merged in each bin (only useful for the merging algorithm)
extern int SFH_ibin[MAXSNAPS][STEPS]; //Last active bin
#ifdef DETAILED_METALS_AND_MASS_RETURN
extern double tau_t[STEPS*MAXSNAPS]; //Time-to-z=0 of every timestep in the code. (Used for SNe rates in yield_integrals.c)
extern double tau_dt[STEPS*MAXSNAPS];//Width of every timestep in the code. (Used for SNe rates in yield_integrals.c)
#endif
#endif //STAR_FORMATION_HISTORY


#ifdef DETAILED_METALS_AND_MASS_RETURN

//Number of interpolated points within the mass ranges for the four types of yield table:
//Edit given by Rob Yates 27/09/2018
//Number of interpolated points within the mass ranges for the four types of yield table:
#define LIFETIME_MASS_NUM 150
#define LIFETIME_Z_NUM 6
#define AGB_MASS_NUM 59 //55 //ROB: 59, when going from 0.85 to 7 Msun
#define AGB_Z_NUM 3
#ifdef PORTINARI
#define SNII_MASS_NUM 85  //ROB: 85, from 6 <= M[Msun] <= 120. Change SNII_MIN_MASS and SNII_MAX_MASS for shorter ranges.
#define SNII_Z_NUM 5
#endif
#ifdef CHIEFFI
#define SNII_MASS_NUM 81 //ROB: 56 if 7 <= M[Msun] <= 50. 81 if 7 <= M[Msun] <= 120. (NB: You can set SNII_MASS_NUM 81, and SNII_MAX_MASS 50. But DON"T put SNII_MASS_NUM > 81 ever!)
#define SNII_Z_NUM 6
#endif
#define SNIA_MASS_NUM 83 //48 //Number increased after extending range to cover M2 masses (07-02-12)

//Mass ranges for the different modes of ejection:
#define AGB_MIN_MASS 0.85
#define AGB_MAX_MASS 7.0 //6.0
#define SNIA_MIN_MASS 3.0
#define SNIA_MAX_MASS 16.0
#define SNIA_MIN_TIME 35.0*1.0e6
#define SNIA_MAX_TIME 21.0*1.0e9
#ifdef PORTINARI
#define SNII_MIN_MASS 7.0 //6.0
#define SNII_MAX_MASS 120.0
#endif
#ifdef CHIEFFI
#define SNII_MIN_MASS 7.0
#define SNII_MAX_MASS 120.0 //50.0
#endif

int ELETOBIGCOUNTA;
int FRACCOUNTA;

//Arrays that yield tables are written to:
float lifetimeMasses[LIFETIME_MASS_NUM];
float lifetimeMetallicities[LIFETIME_Z_NUM];
float lifetimes[LIFETIME_Z_NUM][LIFETIME_MASS_NUM];
float AGBMasses[AGB_MASS_NUM]; //Initial star masses [Msun]
float AGBMetallicities[AGB_Z_NUM]; //Initial star metallicities [Msun]
float AGBEjectedMasses[AGB_Z_NUM][AGB_MASS_NUM]; //Total mass ejected [Msun]
float AGBTotalMetals[AGB_Z_NUM][AGB_MASS_NUM]; //Total metal YIELD ejected [Msun]
float AGBYields[AGB_Z_NUM][11][AGB_MASS_NUM]; //YIELD ejected, for each element [Msun]
float SNIIMasses[SNII_MASS_NUM];
float SNIIMetallicities[SNII_Z_NUM];
float SNIIEjectedMasses[SNII_Z_NUM][SNII_MASS_NUM];
float SNIITotalMetals[SNII_Z_NUM][SNII_MASS_NUM];
float SNIIYields[SNII_Z_NUM][11][SNII_MASS_NUM];
#ifndef DTD
float SNIaMasses[SNIA_MASS_NUM];
float SNIaEjectedMasses[SNIA_MASS_NUM];
float SNIaTotalMetals[SNIA_MASS_NUM];
float SNIaYields[42][SNIA_MASS_NUM];
#else
float SNIaYields[42];
#endif

//Integrated yields arrays:
float NormSNIIMassEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
float NormSNIIMetalEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
#ifdef INDIVIDUAL_ELEMENTS
float NormSNIIYieldRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM][NUM_ELEMENTS];
#endif
float NormAGBMassEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
float NormAGBMetalEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
#ifdef INDIVIDUAL_ELEMENTS
float NormAGBYieldRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM][NUM_ELEMENTS];
#endif
float NormSNIaMassEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
float NormSNIaMetalEjecRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM];
#ifdef INDIVIDUAL_ELEMENTS
float NormSNIaYieldRate[STEPS*MAXSNAPS][SFH_NBIN][LIFETIME_Z_NUM][NUM_ELEMENTS];
#endif

//Arrays used to plot SNe rates from SFH bins (yield_integrals.c):
float TheSFH[SFH_NBIN];
float SNIIRate[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
float SNIaRate[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
float AGBRate[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
//Arrays used to plot SNe rates from SFH-timesteps (calc_SNe_rates.c):
float TheSFH2[STEPS*MAXSNAPS];
float SNIIRate2[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
float SNIaRate2[STEPS*MAXSNAPS][LIFETIME_Z_NUM];
float AGBRate2[STEPS*MAXSNAPS][LIFETIME_Z_NUM];


//IMF parameters (for chemical enrichment):
#define IMF_SLOPE 2.3 //2.15 //High-mass slope of Chabrier IMF. (2.3 = normal Chabrier IMF. <2.3 = top-heavy, >2.3 = bottom-heavy.)

//SNIa parameters:
#define A_FACTOR 0.028 //Fraction of mass from all objects between SNIa_MIN_MASS and SNIA_MAX_MASS that comes from SN-Ia. //0.028 preferred in Yates+13.
#ifdef DTD
#define SNIAEJECMASS 1.2300971 //Total mass (and total metals) ejected by a SNIa explosion in Msun //Value form original yield table (42 elements): 1.3740855. //Value when only considering 11 elements: 1.2300971
#ifdef BIMODALDTD
	//#define DTD_NORM 0.903206 //(26Myrs - 21Gyrs)
	//#define DTD_NORM 0.896668 //For P98 Z=0.02 lifetimes (35Myrs - 17Gyrs)
	#define DTD_NORM 0.900348 //(35Myrs - 21Gyrs)
#endif
#ifdef CUSTOMDTD
	//#define DTD_NORM 0.524836 //(26Myrs - 21Gyrs)
	//#define DTD_NORM 0.606746 //For P98 Z=0.02 lifetimes (35Myrs - 17Gyrs, for ~32% in prompt component)
	#define DTD_NORM 0.610431 //(35Myrs - 21Gyrs, for ~32% in prompt component)
#endif
#ifdef GAUSSIANDTD
	#define DTD_NORM = 1.0 //For P98 Z=0.02 lifetimes (35Myrs - 17Gyrs) and (35Myrs - 21Gyrs)
	#define TAUCHARAC 1.0 //Characteristic delay time for SNe-Ia (i.e. peak of Gaussian distribution) in Gyrs //default: 2.0
	#define SIGMA_TD 0.2*TAUCHARAC //0.2 for narrow-DTD, 0.5 for wide_DTD
#endif
#ifdef POWERLAWDTD
	//#define DTD_NORM 7.21863 // (26Myrs - 21Gyrs)
	#define DTD_NORM 6.72574 //6.72544 // (35Myrs - 21Gyrs)
	//#define DTD_NORM 6.56087 //For P98 Z=0.02 lifetimes (35Myrs - 17Gyrs)
	//#define DTD_NORM 6.22432 // (35Myrs - 11Gyrs)
	//#define DTD_NORM 6.02197 // (50Myrs - 17Gyrs)
	//#define DTD_NORM 6.35503 // (40Myrs - 17Gyrs)
	#define DTD_SLOPE -1.12 //Slope of power law, according to Maoz et al. (2012)
#endif
#ifdef RUITERDTD
	//#define DTD_NORM 1.09545 //(26Myrs - 21Gyrs)
	#define DTD_NORM 1.09545 //(35Myrs - 21Gyrs)
	//#define DTD_NORM 1.08422 //For P98 Z=0.02 lifetimes (35Myrs - 17Gyrs)
	#define TAUCHARAC 0.5 //Peak of Gaussian (prompt) component [in Gyrs]
	#define SIGMA_TD 0.2*TAUCHARAC //Width of Gaussian (prompt) component
	#define DTD_SLOPE -2.0 //Slope of power law (delayed) component (see Ruiter et al. 2012)
#endif
#endif

double TotalMassReturnedToColdDiskGas;
double TotalMassReturnedToHotGas;

#endif //DETAILED_METALS_AND_MASS_RETURN

#ifdef DETAILED_DUST

// Numbers for reading in dust tables in dustyields_read_tables.c
#define AGB_DUST_MASS_NUM 28
#define AGB_DUST_METAL_NUM 4 
#define AGB_DUST_TYPE_NUM 11

float AGBDustMasses[AGB_DUST_MASS_NUM]; //Initial star masses [Msun]
float AGBDustMetallicities[AGB_DUST_METAL_NUM]; //Initial star metallicities [Msun]
float AGBDustCreated[AGB_DUST_METAL_NUM][AGB_DUST_MASS_NUM][AGB_DUST_TYPE_NUM]; //Total mass ejected [Msun]

float NormAGBDustYieldRate[(STEPS*MAXSNAPS)][SFH_NBIN][LIFETIME_Z_NUM][AGB_DUST_TYPE_NUM]; 

// Variables to hold amount of metals created in prev time step in model_yields.c
float SNII_prevstep_Cold_Si[SFH_NBIN];
float SNII_prevstep_Cold_Fe[SFH_NBIN];
float SNII_prevstep_Cold_Cb[SFH_NBIN];
float SNIa_prevstep_Cold_Fe[SFH_NBIN];

// For storing metallicity in model_yields.c
int Zi_saved[SFH_NBIN];
double Zi_disp_saved[SFH_NBIN];

#endif //DETAILED_DUST


#ifdef COMPUTE_SPECPHOT_PROPERTIES
// SSP PHOT_TABLES - magnitues of starburst population as a function of age

#ifdef M05
#define SSP_NAGES 220		// Age grid of the SSP tables
#define SSP_NMETALLICITES 4			// Number of Metalicities used
#ifdef SPEC_PHOTABLES_ON_THE_FLY
#define SSP_NLambda 1221
#endif
#endif

#ifdef BC03
#define SSP_NAGES 221		// Age grid of the SSP tables
#define SSP_NMETALLICITES 6			// Number of Metalicities used
#ifdef SPEC_PHOTABLES_ON_THE_FLY
#define SSP_NLambda 1221
#endif
#endif

#ifdef CB07
#define SSP_NAGES 221		// Age grid of the SSP tables
#define SSP_NMETALLICITES 6			// Number of Metalicities used
#ifdef SPEC_PHOTABLES_ON_THE_FLY
#define SSP_NLambda 1238
#endif
#endif

//table containing the Metallicity grid of the SSP tables (converted to log10)
extern float SSP_logMetalTab[SSP_NMETALLICITES];
//table containing the Age grid of the SSP tables (originally in years, converted to log10(internal time units 1e12 Yrs/h))
extern float SSP_logAgeTab[SSP_NAGES];
//table containing redshift (different from the one in the code when scaling to future times)
extern float RedshiftTab[MAXSNAPS];
extern float LumTables[NMAG][SSP_NMETALLICITES][MAXSNAPS][SSP_NAGES];
extern float FilterLambda[NMAG+1];//wavelength of each filter + 1 for V-band

#ifdef SPEC_PHOTABLES_ON_THE_FLY
#define MAX_NLambdaFilter 1000
extern int NLambdaFilter[NMAG];
//VEGA
#define NLambdaVega 3303
#endif

//DUST EXTINCTION
#define ExpTauBCBulge 0.5	// constant extinction for young stars in bulges.
#define MUWIDTH  0.2
#define MUCENTER 0.3
extern long mu_seed;

#endif //COMPUTE_SPECPHOT_PROPERTIES

/*For H2 formation recipe - Not Supported*/
#define RHO_LEN 101
#define Z_LEN 13

extern size_t HighMark;

#ifdef UPDATETYPETWO
extern int NtotHalos, TotIds, Nids, TotSnaps, OffsetIDs;
extern int *CountIDs_halo, *OffsetIDs_halo, *CountIDs_snaptree, *OffsetIDs_snaptree;
extern long long *IdList;
extern float *PosList, *VelList;
#endif


extern int Hashbits;
extern int NumWrittenInParallel;
extern double ScaleFactor;	// factor by which to multiply a position to get its ph index (after floring)


#ifdef USE_MEMORY_TO_MINIMIZE_IO
extern char *ptr_auxdata, *ptr_treedata, *ptr_dbids, *ptr_galaxydata, *ptr_galsnapdata[NOUT];
extern size_t offset_auxdata, offset_treedata, offset_dbids;
extern size_t offset_galaxydata, maxstorage_galaxydata, filled_galaxydata;
extern size_t offset_galsnapdata[NOUT], maxstorage_galsnapdata[NOUT], filled_galsnapdata[NOUT];
#endif

#endif

extern float Reion_z[46],Reion_Mc[46];

extern FILE *tree_file;
extern FILE *treeaux_file;
extern FILE *treedbids_file;
extern FILE *FdGalTree;
extern FILE *FdGalTreeSFH;
extern FILE *FdGalDumps[NOUT];

#ifdef HDF5_OUTPUT
int b[NOUT];
struct GALAXY_OUTPUT galaxy_output_hdf5[NOUT][NRECORDS_APP];
#endif //HDF5_OUTPUT




