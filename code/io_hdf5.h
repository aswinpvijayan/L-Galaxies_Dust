#include "allvars.h"
#include "proto.h"
#include "hdf5.h"
#include "hdf5_hl.h"
 
#define NRECORDS (hsize_t) 0
#define HDF5_STRING_SIZE 2048
 
/* File to set all the HDF5 table properties */ 
int ifield;
int rowsize=0;
hsize_t chunk_size=CHUNK_SIZE;
int * fill_data=NULL;
hid_t file_id;
size_t * output_offsets;
hid_t field_type;
hid_t * field_types;
size_t output_size;
size_t * output_sizes;
 
// Define the dimensions for the HDF5 table
// Should be as large as the number of possible array dimensions+1
hsize_t dims[6];
 
// The number of fields in the data
int nfields=74;
 
// The field names
const char * field_names[]={
"Type",
"SnapNum",
"HaloIndex",
"LookBackTimeToSnap",
"CentralMvir",
"CentralRvir",
"DistanceToCentralGal",
"Pos",
"Vel",
"Len",
"Mvir",
"Rvir",
"Vvir",
"Vmax",
"GasSpin",
"StellarSpin",
"InfallVmax",
"InfallVmaxPeak",
"InfallSnap",
"InfallHotGas",
"HotRadius",
"OriMergTime",
"MergTime",
"ColdGas",
"StellarMass",
"BulgeMass",
"DiskMass",
"HotGas",
"EjectedMass",
"BlackHoleMass",
"ICM",
"MetalsColdGas",
"MetalsBulgeMass",
"MetalsDiskMass",
"MetalsHotGas",
"MetalsEjectedMass",
"MetalsICM",
"PrimordialAccretionRate",
"CoolingRadius",
"CoolingRate",
"CoolingRate_beforeAGN",
"QuasarAccretionRate",
"RadioAccretionRate",
"Sfr",
"SfrBulge",
"XrayLum",
"BulgeSize",
"StellarDiskRadius",
"GasDiskRadius",
"StellarHalfMassRadius",
"StellarHalfLightRadius",
"CosInclination",
"DisruptOn",
"MergeOn",
"MagDust",
"Mag",
"MagBulge",
"MassWeightAge",
"DiskMass_elements",
"BulgeMass_elements",
"ColdGas_elements",
"HotGas_elements",
"ICM_elements",
"EjectedMass_elements",
"ColdGasDiff_elements",
"ColdGasClouds_elements",
"mu_gas",
"DustColdGasRates",
"DustColdGasDiff_elements",
"DustColdGasClouds_elements",
"t_acc",
"f_i",
"f_c",
"f_cmax",
};
 
// Information describing the datatypes and array flags
// that we will use to construct HDF5 field_types.
 
char types[]={
'i',
'i',
'i',
'f',
'f',
'f',
'f',
'f',
'f',
'i',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'i',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'm',
'm',
'm',
'm',
'm',
'm',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'f',
'i',
'i',
'f',
'f',
'f',
'f',
'e',
'e',
'e',
'e',
'e',
'e',
'e',
'e',
'f',
'd',
'e',
'e',
'f',
'f',
'f',
'f',
}; 
 
int flag3[]={
0,
0,
0,
0,
0,
0,
1,
1,
1,
0,
0,
0,
0,
0,
1,
1,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
}; 
 
int flagMag[]={
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
1,
1,
1,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
}; 
 
int flagRings[]={
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
}; 
 
int flagSFH[]={
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
}; 
 
int flagElements[]={
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
}; 