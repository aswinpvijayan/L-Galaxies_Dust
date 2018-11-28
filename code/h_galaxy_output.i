# 1 "./code/h_galaxy_output.h"
# 1 "/lustre/scratch/astro/ap629/sc558/Dust_Rob//"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "./code/h_galaxy_output.h"
/**
 * Galaxy structure for output
 */
# 35 "./code/h_galaxy_output.h"
#pragma pack(1)
struct GALAXY_OUTPUT
{
# 59 "./code/h_galaxy_output.h"
  int Type; // None // Galaxy type: 0 for central galaxies of a main halo; 1 for central galaxies in sub-halos; 2 for satellites without halo.

  int SnapNum; // None // The snapshot number where this galaxy was identified.
  int HaloIndex; // None // Unique ID of MPA halo containing this galaxy
# 72 "./code/h_galaxy_output.h"
  float LookBackTimeToSnap; // yr // The time from a given snapshot to z=0, in years
  float CentralMvir; // 10^10/h Msun // virial mass of background (FOF) halo containing this galaxy
  float CentralRvir; // Mpc/h // Rvir of background (FOF) halo containing this galaxy
  float DistanceToCentralGal[3]; // Mpc/h? // Separation of galaxy and the central galaxy in the halo
  /* properties of subhalo at the last time this galaxy was a central galaxy */
  float Pos[3]; // 1/h Mpc // Galaxy Positions
  float Vel[3]; // km/s // Galaxy Velocities
  int Len; // None // Number of particles in the associated subhalo  
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

  struct metals MetalsColdGas; // 10^10/h Msun // Mass in metals in cold gas.
  struct metals MetalsBulgeMass; // 10^10/h Msun // Mass in metals in the bulge
  struct metals MetalsDiskMass; // 10^10/h Msun // Mass in metals in the disk
  struct metals MetalsHotGas; // 10^10/h Msun // Mass in metals in the hot gas
  struct metals MetalsEjectedMass; // 10^10/h Msun //Mass in metals in the ejected gas
  struct metals MetalsICM; // 10^10/h Msun // total mass in metals in intra-cluster stars, for type 0,1
# 129 "./code/h_galaxy_output.h"
  /* misc */
    float PrimordialAccretionRate; // Msun/yr // Accretion rate of primordial gas.
    float CoolingRadius; // Mpc/h // The radius within which the cooling time scale is shorter than the dynamical timescale
    float CoolingRate; // Msun/yr // Cooling rate of the hot gas
    float CoolingRate_beforeAGN; // Msun/yr // What the cooling rate of the hot gas would have been if there was no AGN feedback.
    float QuasarAccretionRate; // Msun/yr // Rate at which cold gas is accreted into the central black hole in the quasar mode.
    float RadioAccretionRate; // Msun/yr // Rate at which hot gas is accreted into the central black hole in the radio mode.
    float Sfr; // Msun/yr // Star formation rate



    float SfrBulge; // Msun/yr // Star formation rate in bulge.
    float XrayLum; // log10(erg/sec) // (log_10 of) X-Ray luminosity
    float BulgeSize; // Mpc/h // Half mass radius of bulge
    float StellarDiskRadius; // Mpc/h // Size of the stellar disk, 3x the scale length.
    float GasDiskRadius; // Mpc/h // Size of the gas disk, 3x the scale length.
    float StellarHalfMassRadius; // Mpc/h // stellar Half mass radius
    float StellarHalfLightRadius; // Mpc/h // stellar Half light radius
    float CosInclination; // deg // Inclination of the galaxy. Derived from the angle between the stellar spins of the galaxy and the z-axis
    int DisruptOn; // None // 0: galaxy merged onto merger center 1: galaxy was disrupted before merging onto its descendant, matter went into ICM of merger center;
    int MergeOn; // None // 0: standard delucia-like merger behaviour for type 1 galaxy; 1: galaxy mass > halo mass, separate dynamical friction time calculated ....

  /* magnitudes in various bands */


  float MagDust[2]; // None // dust corrected rest-frame absolute mags
  float Mag[2]; // None // rest-frame absolute mags
  float MagBulge[2]; // None // rest-frame absolute mags for the bulge
# 191 "./code/h_galaxy_output.h"
  float MassWeightAge; // yr // Age weighted by stellar mass
# 224 "./code/h_galaxy_output.h"
  struct elements DiskMass_elements; // ? // Mass of elements in disk
  struct elements BulgeMass_elements; // ? // Mass of elements in bulge
  struct elements ColdGas_elements; // ? // Mass of elements in cold gas
  struct elements HotGas_elements; // ? // Mass of elements in hot gas
  struct elements ICM_elements; // ? // Mass of elements in ICM
  struct elements EjectedMass_elements; // ? // Mass of elements in ejected gas
  struct elements ColdGasDiff_elements; // Msol // Mass of elements in the diffused phase (in ColdGas)
  struct elements ColdGasClouds_elements; // Msol // Mass of elements in the cloud phase (in ColdGas)
  float mu_gas; // ? // fraction of molecular media




  struct DustRates DustColdGasRates; // ? // Rates of creation and destruction of dust

  struct elements DustColdGasDiff_elements; // Msol // Mass of elements in the diffused phase locked up in dust (in ColdGas)
  struct elements DustColdGasClouds_elements; // Msol // Mass of elements in the cloud phase locked up in dust (in ColdGas)
  float t_acc; // ? // dust growth rate
  float f_i[9]; // ? // fraction of dust in diffused media
  float f_c[9]; // ? // fraction of dust in molecular media
  float f_cmax[9]; // ? // max fraction of molecular media

};

// next only of interest to DB output, which generally requires complete tree

struct SFH_BIN {
 long long GalID; // ID of the galaxy
 short snapnum; // snapnum of the galaxy, repeated here for faster lookups of times etc
    short sfh_ibin; //Index of highest bin currently in use
//    float sfh_time; //time to present at the middle of bin in years.
//    float sfh_dt; //time width of bin in years.
  float sfh_DiskMass;
  float sfh_BulgeMass;
  float sfh_ICM;

  struct metals sfh_MetalsDiskMass; // Metals locked up in stars in disk.
  struct metals sfh_MetalsBulgeMass; // Metals locked up in stars in bulge.
  struct metals sfh_MetalsICM; // Metals locked up in stars in ICM.
# 273 "./code/h_galaxy_output.h"
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
  struct elements ColdGasDiff_elements;
  struct elements ColdGasClouds_elements;
  float mu_gas;

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

#pragma pack()
