CREATE TABLE GALAXIES (
 type INTEGER NOT NULL 
,  haloIndex INTEGER NOT NULL 
,  snapNum INTEGER NOT NULL 
,  lookBackTimeToSnap REAL NOT NULL 
,  centralMvir REAL NOT NULL 
,  centralRvir REAL NOT NULL 
,  distanceToCentralGal_1 REAL NOT NULL 
,  distanceToCentralGal_2 REAL NOT NULL 
,  distanceToCentralGal_3 REAL NOT NULL 
,  pos_1 REAL NOT NULL 
,  pos_2 REAL NOT NULL 
,  pos_3 REAL NOT NULL 
,  vel_1 REAL NOT NULL 
,  vel_2 REAL NOT NULL 
,  vel_3 REAL NOT NULL 
,  len INTEGER NOT NULL 
,  mvir REAL NOT NULL 
,  rvir REAL NOT NULL 
,  vvir REAL NOT NULL 
,  vmax REAL NOT NULL 
,  gasSpin_1 REAL NOT NULL 
,  gasSpin_2 REAL NOT NULL 
,  gasSpin_3 REAL NOT NULL 
,  stellarSpin_1 REAL NOT NULL 
,  stellarSpin_2 REAL NOT NULL 
,  stellarSpin_3 REAL NOT NULL 
,  infallVmax REAL NOT NULL 
,  infallVmaxPeak REAL NOT NULL 
,  infallSnap INTEGER NOT NULL 
,  infallHotGas REAL NOT NULL 
,  hotRadius REAL NOT NULL 
,  oriMergTime REAL NOT NULL 
,  mergTime REAL NOT NULL 
,  coldGas REAL NOT NULL 
,  stellarMass REAL NOT NULL 
,  bulgeMass REAL NOT NULL 
,  diskMass REAL NOT NULL 
,  hotGas REAL NOT NULL 
,  ejectedMass REAL NOT NULL 
,  blackHoleMass REAL NOT NULL 
,  iCM REAL NOT NULL 
,  metals struct NOT NULL 
,  metals struct NOT NULL 
,  metals struct NOT NULL 
,  metals struct NOT NULL 
,  metals struct NOT NULL 
,  metals struct NOT NULL 
,  primordialAccretionRate REAL NOT NULL 
,  coolingRadius REAL NOT NULL 
,  coolingRate REAL NOT NULL 
,  coolingRate_beforeAGN REAL NOT NULL 
,  quasarAccretionRate REAL NOT NULL 
,  radioAccretionRate REAL NOT NULL 
,  sfr REAL NOT NULL 
,  sfrBulge REAL NOT NULL 
,  xrayLum REAL NOT NULL 
,  bulgeSize REAL NOT NULL 
,  stellarDiskRadius REAL NOT NULL 
,  gasDiskRadius REAL NOT NULL 
,  cosInclination REAL NOT NULL 
,  disruptOn INTEGER NOT NULL 
,  mergeOn INTEGER NOT NULL 
,  magDust_1 REAL NOT NULL 
,  magDust_2 REAL NOT NULL 
,  mag_1 REAL NOT NULL 
,  mag_2 REAL NOT NULL 
,  magBulge_1 REAL NOT NULL 
,  magBulge_2 REAL NOT NULL 
,  massWeightAge REAL NOT NULL 
,  elements struct NOT NULL 
,  elements struct NOT NULL 
,  elements struct NOT NULL 
,  elements struct NOT NULL 
,  elements struct NOT NULL 
,  elements struct NOT NULL 
,  dustRates struct NOT NULL 
,  elements struct NOT NULL 
,  gALAXY_OUTPUT_1 struct NOT NULL 
,  gALAXY_OUTPUT_2 struct NOT NULL 
,  gALAXY_OUTPUT_3 struct NOT NULL 
,  gALAXY_OUTPUT_4 struct NOT NULL 
,  gALAXY_OUTPUT_5 struct NOT NULL 
,  gALAXY_OUTPUT_6 struct NOT NULL 
,  gALAXY_OUTPUT_7 struct NOT NULL 
,  gALAXY_OUTPUT_8 struct NOT NULL 
,  gALAXY_OUTPUT_9 struct NOT NULL 
,  gALAXY_OUTPUT_10 struct NOT NULL 
,  gALAXY_OUTPUT_11 struct NOT NULL 
,  gALAXY_OUTPUT_12 struct NOT NULL 
,  gALAXY_OUTPUT_13 struct NOT NULL 
,  gALAXY_OUTPUT_14 struct NOT NULL 
,  gALAXY_OUTPUT_15 struct NOT NULL 
,  gALAXY_OUTPUT_16 struct NOT NULL 
,  gALAXY_OUTPUT_17 struct NOT NULL 
)
CREATE TABLE SFH_Times (
 snapnum INTEGER NOT NULL 
,  bin INTEGER NOT NULL 
,  lookbacktime double NOT NULL 
,  dt double NOT NULL 
,  nbins INTEGER NOT NULL 
 -- size = 20
)