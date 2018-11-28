struct LGalaxy {
   int Type;
   int HaloIndex;
   int SnapNum;
   float LookBackTimeToSnap;
   float CentralMvir;
   float CentralRvir;
   float DistanceToCentralGal[3];
   float Pos[3];
   float Vel[3];
   int Len;
   float Mvir;
   float Rvir;
   float Vvir;
   float Vmax;
   float GasSpin[3];
   float StellarSpin[3];
   float InfallVmax;
   float InfallVmaxPeak;
   int InfallSnap;
   float InfallHotGas;
   float HotRadius;
   float OriMergTime;
   float MergTime;
   float ColdGas;
   float StellarMass;
   float BulgeMass;
   float DiskMass;
   float HotGas;
   float EjectedMass;
   float BlackHoleMass;
   float ICM;
   struct metals MetalsColdGas;
   struct metals MetalsBulgeMass;
   struct metals MetalsDiskMass;
   struct metals MetalsHotGas;
   struct metals MetalsEjectedMass;
   struct metals MetalsICM;
   float PrimordialAccretionRate;
   float CoolingRadius;
   float CoolingRate;
   float CoolingRate_beforeAGN;
   float QuasarAccretionRate;
   float RadioAccretionRate;
   float Sfr;
   float SfrBulge;
   float XrayLum;
   float BulgeSize;
   float StellarDiskRadius;
   float GasDiskRadius;
   float CosInclination;
   int DisruptOn;
   int MergeOn;
   float MagDust[2];
   float Mag[2];
   float MagBulge[2];
   float MassWeightAge;
   struct elements DiskMass_elements;
   struct elements BulgeMass_elements;
   struct elements ColdGas_elements;
   struct elements HotGas_elements;
   struct elements ICM_elements;
   struct elements EjectedMass_elements;
   struct elements ColdGasDiff_elements;
  struct elements ColdGasClouds_elements;
  double mu_gas;
  struct DustRates DustISMRates;
  struct elements DustColdGasDiff_elements; 
  struct elements DustColdGasClouds_elements;
  double f_i[9];
  double f_c[9];
  double f_cmax[9];
  double t_acc;
 struct GALAXY_OUTPUT galaxy_output_hdf5[17][1000];
};
struct MoMaFGalaxy {
  long long GalID;
  short snapnum;
     short sfh_ibin;
   float sfh_DiskMass;
   float sfh_BulgeMass;
   float sfh_ICM;
   struct metals sfh_MetalsDiskMass;
   struct metals sfh_MetalsBulgeMass;
   struct metals sfh_MetalsICM;
   struct elements sfh_ElementsDiskMass;
   struct elements sfh_ElementsBulgeMass;
   struct elements sfh_ElementsICM;
   struct elements DiskMass_elements;
   struct elements BulgeMass_elements;
   struct elements ColdGas_elements;
   struct elements HotGas_elements;
   struct elements ICM_elements;
   struct elements EjectedMass_elements;
};
