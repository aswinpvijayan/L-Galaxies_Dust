/**
 * Galaxy structure for output
 */
#ifdef LIGHT_OUTPUT
struct GALAXY_OUTPUT
{
  int   Type; // Galaxy type: 0 for central galaxies of a main halo, 1 for central galaxies in sub-halos, 2 for satellites without halo.
  int   SnapNum; // The snapshot number where this galaxy was identified.
  float CentralMvir; // 10^10/h Msun virial mass of background (FOF) halo containing this galaxy
  float CentralRvir; // 10^10/h Msun virial mass of background (FOF) halo containing this galaxy
  float Pos[3]; // 1/h Mpc - Galaxy Positions
  float Mvir; // 10^10/h Msun - Virial mass of the subhalo the galaxy is/was the center of.
  float Rvir; // Mpc/h - Virial radius of the subhalo the galaxy is/was the center of.
  float Vvir; // km/s -	Virial velocity of the subhalo the galaxy is/was the center of.
  float DistanceToCentralGal[3];

  /* baryonic reservoirs */
  float ColdGas; // 10^10/h Msun - Mass in cold gas.
  float BulgeMass; // 10^10/h Msun - Mass in the bulge
  float DiskMass;
  float HotGas; // 10^10/h Msun - Mass in hot gas
  float BlackHoleMass; // 10^10/h Msun - Mass in black hole

  /* magnitudes in various bands */
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef OUTPUT_OBS_MAGS
  float ObsMagDust[NMAG]; // dust corrected, rest-frame absolute mags
#endif
#ifdef OUTPUT_REST_MAGS
  float MagDust[NMAG]; // dust corrected, rest-frame absolute mags
#endif
#endif
};
#else
#pragma pack(1)  //structure alignment for 1 Byte.
struct GALAXY_OUTPUT
{
#ifdef GALAXYTREE
  long long GalID; /** ID of galaxy, unique within simulation and SAM run.*/
  long long HaloID; // Unique ID of MPA halo containing this galaxy
#endif
#ifdef MBPID
  long long MostBoundID; // Most bound particle at centre of subhalo last associated with this galaxy.  Put here as want all 8-byte blocks together at top of output record.
#endif
#ifdef GALAXYTREE
  long long FirstProgGal;	// Main progenitor of this galaxy. Also the first progenitor in a linked list representation of the merger tree.
  long long NextProgGal;	// Next progenitor of this galaxy in linked list representation of merger tree
  long long LastProgGal;	// Galaxies with id between this galaxyId and this lastProgenitorId form the merger tree rooted in this galaxy.
  long long FOFCentralGal;
  long long FileTreeNr;
  long long DescendantGal;	// Pointer to the descendant of this galaxy in its merger tree; -1 if there is no descendant
  long long MainLeafId;
  long long TreeRootId;
  long long SubID;
  long long MMSubID; // fofId, the subhaloid of the subhalo at the center of the fof group
  int   PeanoKey; // Peano-Hilbert key, (bits=8), for position in 500/h Mpc box
  float Redshift; // redshift of the snapshot where this galaxy resides
#endif
  int   Type; // None // Galaxy type: 0 for central galaxies of a main halo; 1 for central galaxies in sub-halos; 2 for satellites without halo.
#ifndef GALAXYTREE
  int   SnapNum; // None // The snapshot number where this galaxy was identified.
  int   HaloIndex; // None // Unique ID of MPA halo containing this galaxy
#endif
#ifdef HALOPROPERTIES
  float HaloM_Mean200, HaloM_Crit200, HaloM_TopHat;
  float HaloPos[3];
  float HaloVel[3];
  float HaloVelDisp;
  float HaloVmax;
  float HaloSpin[3];
#endif
  float LookBackTimeToSnap; // yr // The time from a given snapshot to z=0, in years
  float CentralMvir; // 10^10/h Msun // virial mass of background (FOF) halo containing this galaxy
  float CentralRvir; // Mpc/h // Rvir of background (FOF) halo containing this galaxy
  float DistanceToCentralGal[3]; // Mpc/h? // Separation of galaxy and the central galaxy in the halo
  /* properties of subhalo at the last time this galaxy was a central galaxy */
  float Pos[3]; // 1/h Mpc // Galaxy Positions
  float Vel[3]; // km/s // Galaxy Velocities
  int   Len; // None // Number of particles in the associated subhalo  
  float Mvir; // 10^10/h Msun // Virial mass of the subhalo the galaxy is/was the centre of.
  float Rvir; // Mpc/h // Virial radius of the subhalo the galaxy is/was the centre of.
  float Vvir; // km/s //	Virial velocity of the subhalo the galaxy is/was the centre of.
  float Vmax; // km/s // Maximum rotational velocity of the subhalo, or the last value for type 2's galaxies.
  float GasSpin[3]; // Mpc/h km/s // specific ColdGas Spin
  float StellarSpin[3]; // Mpc/h km/s // specific Stellar Disk? Spin
  float InfallVmax; // km/s // Vmax at infall
  float InfallVmaxPeak; // km/s // Max previous Vmax at infall
  int InfallSnap; // None // Snapnum at infall
  float InfallHotGas; // 10^10/h Msun ? // ?
  float HotRadius; //Mpc/h // Radius of the hot gas
  /*dynamical friction merger time*/
  float OriMergTime; // ? // ?
  float MergTime; // ? // ?
  /* baryonic reservoirs */
  float ColdGas; // 10^10/h Msun // Mass in cold gas.
  float StellarMass; // 10^10/h Msun // Disk+Bulge
  float BulgeMass; // 10^10/h Msun // Mass in the bulge
  float DiskMass; // 10^10/h Msun // Mass in the disk
  float HotGas; // 10^10/h Msun // Mass in hot gas
  float EjectedMass; // 10^10/h Msun // Mass in ejected gas
  float BlackHoleMass; // 10^10/h Msun // Mass in black hole
  /* ICL magnitude and mass*/
  float ICM; // 10^10/h Msun // mass in intra-cluster stars for type 0 & 1
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals MetalsColdGas; // 10^10/h Msun // Mass in metals in cold gas.
  struct metals MetalsBulgeMass; // 10^10/h Msun // Mass in metals in the bulge
  struct metals MetalsDiskMass; // 10^10/h Msun // Mass in metals in the disk
  struct metals MetalsHotGas; // 10^10/h Msun // Mass in metals in the hot gas
  struct metals MetalsEjectedMass; // 10^10/h Msun //Mass in metals in the ejected gas
  struct metals MetalsICM; // 10^10/h Msun // total mass in metals in intra-cluster stars, for type 0,1
#ifdef METALS_SELF
  struct metals MetalsHotGasSelf; // hot gas metals that come from self
#endif
#else
  float MetalsColdGas; // 10^10/h Msun // Mass in metals in cold gas.
  float MetalsStellarMass; // 10^10/h Msun // Mass in metals in the bulge+disk
  float MetalsBulgeMass; // 10^10/h Msun // Mass in metals in the bulge
  float MetalsDiskMass; // 10^10/h Msun // Mass in metals in the disk
  float MetalsHotGas; // 10^10/h Msun // Mass in metals in the hot gas
  float MetalsEjectedMass; // 10^10/h Msun // Mass in metals in the ejected gas
    float MetalsICM; // 10^10/h Msun // total mass in metals in intra-cluster stars, for type 0,1
#ifdef METALS_SELF
  float MetalsHotGasSelf; // 10^10/h Msun // hot gas metals that come from self
#endif
#endif //DETAILED_METALS_AND_MASS_RETURN
#ifdef TRACK_BURST
  float BurstMass; // 10^10/h Msun // Mass formed in starbursts
#endif
  /* misc */
    float PrimordialAccretionRate; // Msun/yr // Accretion rate of primordial gas.
    float CoolingRadius; // Mpc/h // The radius within which the cooling time scale is shorter than the dynamical timescale
    float CoolingRate; // Msun/yr // Cooling rate of the hot gas
    float CoolingRate_beforeAGN; // Msun/yr // What the cooling rate of the hot gas would have been if there was no AGN feedback.
    float QuasarAccretionRate; // Msun/yr // Rate at which cold gas is accreted into the central black hole in the quasar mode.
    float RadioAccretionRate; // Msun/yr // Rate at which hot gas is accreted into the central black hole in the radio mode.
    float Sfr; // Msun/yr // Star formation rate
#ifdef OUTPUT_RINGS
    float SfrRings[RNUM]; // Msun/yr // Star formation rate within each annular ring
#endif     //H2_AND_RINGS
    float SfrBulge; // Msun/yr // Star formation rate in bulge.
    float XrayLum; // log10(erg/sec) // (log_10 of) X-Ray luminosity
    float BulgeSize; // Mpc/h // Half mass radius of bulge
    float StellarDiskRadius; // Mpc/h // Size of the stellar disk, 3x the scale length.
    float GasDiskRadius; // Mpc/h // Size of the gas disk, 3x the scale length.
    float StellarHalfMassRadius; // Mpc/h // stellar Half mass radius
    float StellarHalfLightRadius; // Mpc/h // stellar Half light radius
    float CosInclination; // deg // Inclination of the galaxy. Derived from the angle between the stellar spins of the galaxy and the z-axis
    int   DisruptOn; // None // 0: galaxy merged onto merger center 1: galaxy was disrupted before merging onto its descendant, matter went into ICM of merger center;
    int   MergeOn; // None // 0: standard delucia-like merger behaviour for type 1 galaxy; 1: galaxy mass > halo mass, separate dynamical friction time calculated ....

  /* magnitudes in various bands */
#ifdef COMPUTE_SPECPHOT_PROPERTIES
#ifdef OUTPUT_REST_MAGS
  float MagDust[NMAG]; // None // dust corrected rest-frame absolute mags
  float Mag[NMAG]; // None // rest-frame absolute mags
  float MagBulge[NMAG]; // None // rest-frame absolute mags for the bulge
#ifdef ICL
  float MagICL[NMAG]; // None  // rest-frame absolute mags of ICL
#endif
#endif

#ifdef OUTPUT_OBS_MAGS
  float ObsMagDust[NMAG]; // None // dust-corrected obs-frame absolute mags
  float ObsMag[NMAG]; // None // obs-frame absolute mags
  float ObsMagBulge[NMAG]; // None // obs-frame absolute mags for the bulge
#ifdef ICL
  float ObsMagICL[NMAG];  // None // observer-frame absolute mags for intra-cluster light
#endif
#ifdef OUTPUT_MOMAF_INPUTS
  // define luminosities as if the galaxy were one snapshot earlier, i.e. higher redshift, than its actual snapshot
  float dObsMagDust[NMAG];
  float dObsMag[NMAG];
  float dObsMagBulge[NMAG];
#ifdef ICL
  float dObsMagICL[NMAG];
#endif	//
#ifdef KITZBICHLER
  // define luminosities as if the galaxy were one snapshot later, i.e. lower redshift, than its actual snapshot
  float dObsMagDust_forward[NMAG];
  float dObsMag_forward[NMAG];
  float dObsMagBulge_forward[NMAG];
#ifdef ICL
  float dObsMagICL_forward[NMAG];
#endif	//
#endif //KITZBICHLER
#endif	//OUTPUT_MOMAF_INPUTS
#endif	//OUTPUT_OBS_MAGS

#endif  //COMPUTE_SPECPHOT_PROPERTIES

  float MassWeightAge; // yr // Age weighted by stellar mass
#ifdef  POST_PROCESS_MAGS
  float rbandWeightAge; // yr // Age weighted by rband flux
#endif

#ifdef STAR_FORMATION_HISTORY
#ifndef REDUCED_OUTPUT
  int sfh_ibin; //Index of highest bin currently in use
  int sfh_numbins; // number of non empty bins
  float sfh_DiskMass[SFH_NBIN];
  float sfh_BulgeMass[SFH_NBIN];
  float sfh_ICM[SFH_NBIN];
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals sfh_MetalsDiskMass[SFH_NBIN]; // Metals locked up in stars in disk.
  struct metals sfh_MetalsBulgeMass[SFH_NBIN]; // Metals locked up in stars in bulge.
  struct metals sfh_MetalsICM[SFH_NBIN]; // Metals locked up in stars in ICM.
#else
  float sfh_MetalsDiskMass[SFH_NBIN]; // Metals locked up in stars in disk.
  float sfh_MetalsBulgeMass[SFH_NBIN]; //Metals locked up in stars in bulge.
  float sfh_MetalsICM[SFH_NBIN]; // Metals locked up in stars in ICM.
#endif
#ifdef TRACK_BURST
  float sfh_BurstMass[SFH_NBIN]; // Mass formed in starbursts
#endif
#endif //REDUCED_OUTPUT
#endif //STAR_FORMATION_HISTORY

#ifdef INDIVIDUAL_ELEMENTS
  #ifndef REDUCED_OUTPUT
  struct elements sfh_ElementsDiskMass[SFH_NBIN];
  struct elements sfh_ElementsBulgeMass[SFH_NBIN];
  struct elements sfh_ElementsICM[SFH_NBIN];
  #endif //REDUCED_OUTPUT
  //float DiskMass_elements[ELEMENT_NUM];
  struct elements DiskMass_elements;
  struct elements BulgeMass_elements;
  struct elements ColdGas_elements;
  struct elements HotGas_elements;
  struct elements ICM_elements;
  struct elements EjectedMass_elements;
#endif //INDIVIDUAL_ELEMENTS

#ifdef DETAILED_DUST
	#ifdef FULL_DUST_RATES
	struct DustRates DustISMRates;		
	#endif
	struct elements Dust_elements; //Dust is now stored as an elements array
#endif //DETAILED_DUST


};

// next only of interest to DB output, which generally requires complete tree
#ifdef STAR_FORMATION_HISTORY
struct SFH_BIN {
	long long GalID; // ID of the galaxy
	short snapnum; // snapnum of the galaxy, repeated here for faster lookups of times etc
    short sfh_ibin; //Index of highest bin currently in use
//    float sfh_time; //time to present at the middle of bin in years.
//    float sfh_dt; //time width of bin in years.
  float sfh_DiskMass;
  float sfh_BulgeMass;
  float sfh_ICM;
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals sfh_MetalsDiskMass; // Metals locked up in stars in disk.
  struct metals sfh_MetalsBulgeMass; // Metals locked up in stars in bulge.
  struct metals sfh_MetalsICM; // Metals locked up in stars in ICM.
#else
  float sfh_MetalsDiskMass; // Metals locked up in stars in disk.
  float sfh_MetalsBulgeMass; //Metals locked up in stars in bulge.
  float sfh_MetalsICM; // Metals locked up in stars in ICM.
#endif
#ifdef TRACK_BURST
  float sfh_BurstMass; // Mass formed in starbursts
#endif

#ifdef INDIVIDUAL_ELEMENTS
  struct elements sfh_ElementsDiskMass;
  struct elements sfh_ElementsBulgeMass;
  struct elements sfh_ElementsICM;

  //float DiskMass_elements[ELEMENT_NUM];
  struct elements DiskMass_elements;
  struct elements BulgeMass_elements;
  struct elements ColdGas_elements;
  struct elements HotGas_elements;
  struct elements ICM_elements;
  struct elements EjectedMass_elements;
#endif //INDIVIDUAL_ELEMENTS
};

struct SFH_Time
{
 	int snapnum; // snapnum
 	int bin; // index of current bin
 	double lookbacktime; // lookback time in years (???) to center of current bin
 						// proposal: in output write the start of the bin and its end, rather than center and dt
 	double dt; // width of the current bin in years (???)
 	int nbins; // # of highest resolution bins used to create current bin
};
#endif //STAR_FORMATION_HISTORY
#pragma pack()   //structure alignment ends.
#endif //When LIGHT_OUTPUT is not defined
