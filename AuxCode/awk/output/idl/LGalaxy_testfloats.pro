;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION LGalaxy_testfloats, LGs, nstart 
; test whether floats are NaN or too small for SQLServer
; assumes the existence of a function testFloat accepting an array of floats
 badranges = []
 bad = 0
 sel = testFloat(LGs.LookBackTimeToSnap)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'LookBackTimeToSnap --- ', nstart+sel
     print, 'LookBackTimeToSnap --- ', LGs[sel].LookBackTimeToSnap
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CentralMvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CentralMvir --- ', nstart+sel
     print, 'CentralMvir --- ', LGs[sel].CentralMvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CentralRvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CentralRvir --- ', nstart+sel
     print, 'CentralRvir --- ', LGs[sel].CentralRvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DistanceToCentralGal(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DistanceToCentralGal[0] --- ', nstart+sel
     print, 'DistanceToCentralGal --- ', LGs[sel].DistanceToCentralGal(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DistanceToCentralGal(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DistanceToCentralGal[1] --- ', nstart+sel
     print, 'DistanceToCentralGal --- ', LGs[sel].DistanceToCentralGal(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DistanceToCentralGal(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DistanceToCentralGal[2] --- ', nstart+sel
     print, 'DistanceToCentralGal --- ', LGs[sel].DistanceToCentralGal(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Pos(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Pos[0] --- ', nstart+sel
     print, 'Pos --- ', LGs[sel].Pos(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Pos(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Pos[1] --- ', nstart+sel
     print, 'Pos --- ', LGs[sel].Pos(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Pos(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Pos[2] --- ', nstart+sel
     print, 'Pos --- ', LGs[sel].Pos(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vel(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vel[0] --- ', nstart+sel
     print, 'Vel --- ', LGs[sel].Vel(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vel(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vel[1] --- ', nstart+sel
     print, 'Vel --- ', LGs[sel].Vel(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vel(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vel[2] --- ', nstart+sel
     print, 'Vel --- ', LGs[sel].Vel(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mvir --- ', nstart+sel
     print, 'Mvir --- ', LGs[sel].Mvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Rvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Rvir --- ', nstart+sel
     print, 'Rvir --- ', LGs[sel].Rvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vvir)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vvir --- ', nstart+sel
     print, 'Vvir --- ', LGs[sel].Vvir
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Vmax)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Vmax --- ', nstart+sel
     print, 'Vmax --- ', LGs[sel].Vmax
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.GasSpin(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'GasSpin[0] --- ', nstart+sel
     print, 'GasSpin --- ', LGs[sel].GasSpin(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.GasSpin(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'GasSpin[1] --- ', nstart+sel
     print, 'GasSpin --- ', LGs[sel].GasSpin(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.GasSpin(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'GasSpin[2] --- ', nstart+sel
     print, 'GasSpin --- ', LGs[sel].GasSpin(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarSpin(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarSpin[0] --- ', nstart+sel
     print, 'StellarSpin --- ', LGs[sel].StellarSpin(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarSpin(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarSpin[1] --- ', nstart+sel
     print, 'StellarSpin --- ', LGs[sel].StellarSpin(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarSpin(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarSpin[2] --- ', nstart+sel
     print, 'StellarSpin --- ', LGs[sel].StellarSpin(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.InfallVmax)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'InfallVmax --- ', nstart+sel
     print, 'InfallVmax --- ', LGs[sel].InfallVmax
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.InfallVmaxPeak)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'InfallVmaxPeak --- ', nstart+sel
     print, 'InfallVmaxPeak --- ', LGs[sel].InfallVmaxPeak
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.InfallHotGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'InfallHotGas --- ', nstart+sel
     print, 'InfallHotGas --- ', LGs[sel].InfallHotGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HotRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HotRadius --- ', nstart+sel
     print, 'HotRadius --- ', LGs[sel].HotRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.OriMergTime)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'OriMergTime --- ', nstart+sel
     print, 'OriMergTime --- ', LGs[sel].OriMergTime
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MergTime)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MergTime --- ', nstart+sel
     print, 'MergTime --- ', LGs[sel].MergTime
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ColdGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ColdGas --- ', nstart+sel
     print, 'ColdGas --- ', LGs[sel].ColdGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarMass --- ', nstart+sel
     print, 'StellarMass --- ', LGs[sel].StellarMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.BulgeMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'BulgeMass --- ', nstart+sel
     print, 'BulgeMass --- ', LGs[sel].BulgeMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DiskMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DiskMass --- ', nstart+sel
     print, 'DiskMass --- ', LGs[sel].DiskMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HotGas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HotGas --- ', nstart+sel
     print, 'HotGas --- ', LGs[sel].HotGas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.EjectedMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'EjectedMass --- ', nstart+sel
     print, 'EjectedMass --- ', LGs[sel].EjectedMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.BlackHoleMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'BlackHoleMass --- ', nstart+sel
     print, 'BlackHoleMass --- ', LGs[sel].BlackHoleMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ICM)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ICM --- ', nstart+sel
     print, 'ICM --- ', LGs[sel].ICM
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.PrimordialAccretionRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'PrimordialAccretionRate --- ', nstart+sel
     print, 'PrimordialAccretionRate --- ', LGs[sel].PrimordialAccretionRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CoolingRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CoolingRadius --- ', nstart+sel
     print, 'CoolingRadius --- ', LGs[sel].CoolingRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CoolingRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CoolingRate --- ', nstart+sel
     print, 'CoolingRate --- ', LGs[sel].CoolingRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CoolingRate_beforeAGN)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CoolingRate_beforeAGN --- ', nstart+sel
     print, 'CoolingRate_beforeAGN --- ', LGs[sel].CoolingRate_beforeAGN
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.QuasarAccretionRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'QuasarAccretionRate --- ', nstart+sel
     print, 'QuasarAccretionRate --- ', LGs[sel].QuasarAccretionRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.RadioAccretionRate)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'RadioAccretionRate --- ', nstart+sel
     print, 'RadioAccretionRate --- ', LGs[sel].RadioAccretionRate
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Sfr)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Sfr --- ', nstart+sel
     print, 'Sfr --- ', LGs[sel].Sfr
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.SfrBulge)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'SfrBulge --- ', nstart+sel
     print, 'SfrBulge --- ', LGs[sel].SfrBulge
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.XrayLum)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'XrayLum --- ', nstart+sel
     print, 'XrayLum --- ', LGs[sel].XrayLum
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.BulgeSize)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'BulgeSize --- ', nstart+sel
     print, 'BulgeSize --- ', LGs[sel].BulgeSize
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarDiskRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarDiskRadius --- ', nstart+sel
     print, 'StellarDiskRadius --- ', LGs[sel].StellarDiskRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.GasDiskRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'GasDiskRadius --- ', nstart+sel
     print, 'GasDiskRadius --- ', LGs[sel].GasDiskRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarHalfMassRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarHalfMassRadius --- ', nstart+sel
     print, 'StellarHalfMassRadius --- ', LGs[sel].StellarHalfMassRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarHalfLightRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarHalfLightRadius --- ', nstart+sel
     print, 'StellarHalfLightRadius --- ', LGs[sel].StellarHalfLightRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.CosInclination)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'CosInclination --- ', nstart+sel
     print, 'CosInclination --- ', LGs[sel].CosInclination
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[0] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[1] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[0] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[1] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[0] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[1] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MassWeightAge)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MassWeightAge --- ', nstart+sel
     print, 'MassWeightAge --- ', LGs[sel].MassWeightAge
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.mu_gas)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'mu_gas --- ', nstart+sel
     print, 'mu_gas --- ', LGs[sel].mu_gas
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.t_acc)
 if(sel(0) gt -1) then begin
     bad=1
     print, 't_acc --- ', nstart+sel
     print, 't_acc --- ', LGs[sel].t_acc
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_i(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_i[0] --- ', nstart+sel
     print, 'f_i --- ', LGs[sel].f_i(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_i(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_i[1] --- ', nstart+sel
     print, 'f_i --- ', LGs[sel].f_i(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_i(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_i[2] --- ', nstart+sel
     print, 'f_i --- ', LGs[sel].f_i(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_i(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_i[3] --- ', nstart+sel
     print, 'f_i --- ', LGs[sel].f_i(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_i(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_i[4] --- ', nstart+sel
     print, 'f_i --- ', LGs[sel].f_i(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_i(5))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_i[5] --- ', nstart+sel
     print, 'f_i --- ', LGs[sel].f_i(5)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_i(6))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_i[6] --- ', nstart+sel
     print, 'f_i --- ', LGs[sel].f_i(6)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_i(7))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_i[7] --- ', nstart+sel
     print, 'f_i --- ', LGs[sel].f_i(7)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_i(8))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_i[8] --- ', nstart+sel
     print, 'f_i --- ', LGs[sel].f_i(8)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_c(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_c[0] --- ', nstart+sel
     print, 'f_c --- ', LGs[sel].f_c(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_c(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_c[1] --- ', nstart+sel
     print, 'f_c --- ', LGs[sel].f_c(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_c(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_c[2] --- ', nstart+sel
     print, 'f_c --- ', LGs[sel].f_c(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_c(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_c[3] --- ', nstart+sel
     print, 'f_c --- ', LGs[sel].f_c(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_c(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_c[4] --- ', nstart+sel
     print, 'f_c --- ', LGs[sel].f_c(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_c(5))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_c[5] --- ', nstart+sel
     print, 'f_c --- ', LGs[sel].f_c(5)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_c(6))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_c[6] --- ', nstart+sel
     print, 'f_c --- ', LGs[sel].f_c(6)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_c(7))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_c[7] --- ', nstart+sel
     print, 'f_c --- ', LGs[sel].f_c(7)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_c(8))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_c[8] --- ', nstart+sel
     print, 'f_c --- ', LGs[sel].f_c(8)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_cmax(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_cmax[0] --- ', nstart+sel
     print, 'f_cmax --- ', LGs[sel].f_cmax(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_cmax(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_cmax[1] --- ', nstart+sel
     print, 'f_cmax --- ', LGs[sel].f_cmax(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_cmax(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_cmax[2] --- ', nstart+sel
     print, 'f_cmax --- ', LGs[sel].f_cmax(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_cmax(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_cmax[3] --- ', nstart+sel
     print, 'f_cmax --- ', LGs[sel].f_cmax(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_cmax(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_cmax[4] --- ', nstart+sel
     print, 'f_cmax --- ', LGs[sel].f_cmax(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_cmax(5))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_cmax[5] --- ', nstart+sel
     print, 'f_cmax --- ', LGs[sel].f_cmax(5)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_cmax(6))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_cmax[6] --- ', nstart+sel
     print, 'f_cmax --- ', LGs[sel].f_cmax(6)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_cmax(7))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_cmax[7] --- ', nstart+sel
     print, 'f_cmax --- ', LGs[sel].f_cmax(7)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.f_cmax(8))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'f_cmax[8] --- ', nstart+sel
     print, 'f_cmax --- ', LGs[sel].f_cmax(8)
     badranges=[badranges,sel]
 endif
if(bad) then begin 
     print, 'badranges found: ',badranges
endif
return, badranges
end
