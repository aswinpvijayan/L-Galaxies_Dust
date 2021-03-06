;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO LGalaxy_zerofloats, LGs 
; test whether floats are NaN or too small for SQLServer
; if so, set offending values to 0
; assumes the existence of a function testFloat accepting an array of floats
 sel = testFloat(LGs.LookBackTimeToSnap)
 if(sel(0) gt -1) then begin
     LGs[sel].LookBackTimeToSnap = 0
 endif
 sel = testFloat(LGs.CentralMvir)
 if(sel(0) gt -1) then begin
     LGs[sel].CentralMvir = 0
 endif
 sel = testFloat(LGs.CentralRvir)
 if(sel(0) gt -1) then begin
     LGs[sel].CentralRvir = 0
 endif
 sel = testFloat(LGs.DistanceToCentralGal(0))
 if(sel(0) gt -1) then begin
     LGs[sel].DistanceToCentralGal(0) = 0
 endif
 sel = testFloat(LGs.DistanceToCentralGal(1))
 if(sel(0) gt -1) then begin
     LGs[sel].DistanceToCentralGal(1) = 0
 endif
 sel = testFloat(LGs.DistanceToCentralGal(2))
 if(sel(0) gt -1) then begin
     LGs[sel].DistanceToCentralGal(2) = 0
 endif
 sel = testFloat(LGs.Pos(0))
 if(sel(0) gt -1) then begin
     LGs[sel].Pos(0) = 0
 endif
 sel = testFloat(LGs.Pos(1))
 if(sel(0) gt -1) then begin
     LGs[sel].Pos(1) = 0
 endif
 sel = testFloat(LGs.Pos(2))
 if(sel(0) gt -1) then begin
     LGs[sel].Pos(2) = 0
 endif
 sel = testFloat(LGs.Vel(0))
 if(sel(0) gt -1) then begin
     LGs[sel].Vel(0) = 0
 endif
 sel = testFloat(LGs.Vel(1))
 if(sel(0) gt -1) then begin
     LGs[sel].Vel(1) = 0
 endif
 sel = testFloat(LGs.Vel(2))
 if(sel(0) gt -1) then begin
     LGs[sel].Vel(2) = 0
 endif
 sel = testFloat(LGs.Mvir)
 if(sel(0) gt -1) then begin
     LGs[sel].Mvir = 0
 endif
 sel = testFloat(LGs.Rvir)
 if(sel(0) gt -1) then begin
     LGs[sel].Rvir = 0
 endif
 sel = testFloat(LGs.Vvir)
 if(sel(0) gt -1) then begin
     LGs[sel].Vvir = 0
 endif
 sel = testFloat(LGs.Vmax)
 if(sel(0) gt -1) then begin
     LGs[sel].Vmax = 0
 endif
 sel = testFloat(LGs.GasSpin(0))
 if(sel(0) gt -1) then begin
     LGs[sel].GasSpin(0) = 0
 endif
 sel = testFloat(LGs.GasSpin(1))
 if(sel(0) gt -1) then begin
     LGs[sel].GasSpin(1) = 0
 endif
 sel = testFloat(LGs.GasSpin(2))
 if(sel(0) gt -1) then begin
     LGs[sel].GasSpin(2) = 0
 endif
 sel = testFloat(LGs.StellarSpin(0))
 if(sel(0) gt -1) then begin
     LGs[sel].StellarSpin(0) = 0
 endif
 sel = testFloat(LGs.StellarSpin(1))
 if(sel(0) gt -1) then begin
     LGs[sel].StellarSpin(1) = 0
 endif
 sel = testFloat(LGs.StellarSpin(2))
 if(sel(0) gt -1) then begin
     LGs[sel].StellarSpin(2) = 0
 endif
 sel = testFloat(LGs.InfallVmax)
 if(sel(0) gt -1) then begin
     LGs[sel].InfallVmax = 0
 endif
 sel = testFloat(LGs.InfallVmaxPeak)
 if(sel(0) gt -1) then begin
     LGs[sel].InfallVmaxPeak = 0
 endif
 sel = testFloat(LGs.InfallHotGas)
 if(sel(0) gt -1) then begin
     LGs[sel].InfallHotGas = 0
 endif
 sel = testFloat(LGs.HotRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].HotRadius = 0
 endif
 sel = testFloat(LGs.OriMergTime)
 if(sel(0) gt -1) then begin
     LGs[sel].OriMergTime = 0
 endif
 sel = testFloat(LGs.MergTime)
 if(sel(0) gt -1) then begin
     LGs[sel].MergTime = 0
 endif
 sel = testFloat(LGs.ColdGas)
 if(sel(0) gt -1) then begin
     LGs[sel].ColdGas = 0
 endif
 sel = testFloat(LGs.StellarMass)
 if(sel(0) gt -1) then begin
     LGs[sel].StellarMass = 0
 endif
 sel = testFloat(LGs.BulgeMass)
 if(sel(0) gt -1) then begin
     LGs[sel].BulgeMass = 0
 endif
 sel = testFloat(LGs.DiskMass)
 if(sel(0) gt -1) then begin
     LGs[sel].DiskMass = 0
 endif
 sel = testFloat(LGs.HotGas)
 if(sel(0) gt -1) then begin
     LGs[sel].HotGas = 0
 endif
 sel = testFloat(LGs.EjectedMass)
 if(sel(0) gt -1) then begin
     LGs[sel].EjectedMass = 0
 endif
 sel = testFloat(LGs.BlackHoleMass)
 if(sel(0) gt -1) then begin
     LGs[sel].BlackHoleMass = 0
 endif
 sel = testFloat(LGs.ICM)
 if(sel(0) gt -1) then begin
     LGs[sel].ICM = 0
 endif
 sel = testFloat(LGs.PrimordialAccretionRate)
 if(sel(0) gt -1) then begin
     LGs[sel].PrimordialAccretionRate = 0
 endif
 sel = testFloat(LGs.CoolingRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].CoolingRadius = 0
 endif
 sel = testFloat(LGs.CoolingRate)
 if(sel(0) gt -1) then begin
     LGs[sel].CoolingRate = 0
 endif
 sel = testFloat(LGs.CoolingRate_beforeAGN)
 if(sel(0) gt -1) then begin
     LGs[sel].CoolingRate_beforeAGN = 0
 endif
 sel = testFloat(LGs.QuasarAccretionRate)
 if(sel(0) gt -1) then begin
     LGs[sel].QuasarAccretionRate = 0
 endif
 sel = testFloat(LGs.RadioAccretionRate)
 if(sel(0) gt -1) then begin
     LGs[sel].RadioAccretionRate = 0
 endif
 sel = testFloat(LGs.Sfr)
 if(sel(0) gt -1) then begin
     LGs[sel].Sfr = 0
 endif
 sel = testFloat(LGs.SfrBulge)
 if(sel(0) gt -1) then begin
     LGs[sel].SfrBulge = 0
 endif
 sel = testFloat(LGs.XrayLum)
 if(sel(0) gt -1) then begin
     LGs[sel].XrayLum = 0
 endif
 sel = testFloat(LGs.BulgeSize)
 if(sel(0) gt -1) then begin
     LGs[sel].BulgeSize = 0
 endif
 sel = testFloat(LGs.StellarDiskRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].StellarDiskRadius = 0
 endif
 sel = testFloat(LGs.GasDiskRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].GasDiskRadius = 0
 endif
 sel = testFloat(LGs.StellarHalfMassRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].StellarHalfMassRadius = 0
 endif
 sel = testFloat(LGs.StellarHalfLightRadius)
 if(sel(0) gt -1) then begin
     LGs[sel].StellarHalfLightRadius = 0
 endif
 sel = testFloat(LGs.CosInclination)
 if(sel(0) gt -1) then begin
     LGs[sel].CosInclination = 0
 endif
 sel = testFloat(LGs.MagDust(0))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(0) = 0
 endif
 sel = testFloat(LGs.MagDust(1))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(1) = 0
 endif
 sel = testFloat(LGs.Mag(0))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(0) = 0
 endif
 sel = testFloat(LGs.Mag(1))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(1) = 0
 endif
 sel = testFloat(LGs.MagBulge(0))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(0) = 0
 endif
 sel = testFloat(LGs.MagBulge(1))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(1) = 0
 endif
 sel = testFloat(LGs.MassWeightAge)
 if(sel(0) gt -1) then begin
     LGs[sel].MassWeightAge = 0
 endif
 sel = testFloat(LGs.mu_gas)
 if(sel(0) gt -1) then begin
     LGs[sel].mu_gas = 0
 endif
 sel = testFloat(LGs.t_acc)
 if(sel(0) gt -1) then begin
     LGs[sel].t_acc = 0
 endif
 sel = testFloat(LGs.f_i(0))
 if(sel(0) gt -1) then begin
     LGs[sel].f_i(0) = 0
 endif
 sel = testFloat(LGs.f_i(1))
 if(sel(0) gt -1) then begin
     LGs[sel].f_i(1) = 0
 endif
 sel = testFloat(LGs.f_i(2))
 if(sel(0) gt -1) then begin
     LGs[sel].f_i(2) = 0
 endif
 sel = testFloat(LGs.f_i(3))
 if(sel(0) gt -1) then begin
     LGs[sel].f_i(3) = 0
 endif
 sel = testFloat(LGs.f_i(4))
 if(sel(0) gt -1) then begin
     LGs[sel].f_i(4) = 0
 endif
 sel = testFloat(LGs.f_i(5))
 if(sel(0) gt -1) then begin
     LGs[sel].f_i(5) = 0
 endif
 sel = testFloat(LGs.f_i(6))
 if(sel(0) gt -1) then begin
     LGs[sel].f_i(6) = 0
 endif
 sel = testFloat(LGs.f_i(7))
 if(sel(0) gt -1) then begin
     LGs[sel].f_i(7) = 0
 endif
 sel = testFloat(LGs.f_i(8))
 if(sel(0) gt -1) then begin
     LGs[sel].f_i(8) = 0
 endif
 sel = testFloat(LGs.f_c(0))
 if(sel(0) gt -1) then begin
     LGs[sel].f_c(0) = 0
 endif
 sel = testFloat(LGs.f_c(1))
 if(sel(0) gt -1) then begin
     LGs[sel].f_c(1) = 0
 endif
 sel = testFloat(LGs.f_c(2))
 if(sel(0) gt -1) then begin
     LGs[sel].f_c(2) = 0
 endif
 sel = testFloat(LGs.f_c(3))
 if(sel(0) gt -1) then begin
     LGs[sel].f_c(3) = 0
 endif
 sel = testFloat(LGs.f_c(4))
 if(sel(0) gt -1) then begin
     LGs[sel].f_c(4) = 0
 endif
 sel = testFloat(LGs.f_c(5))
 if(sel(0) gt -1) then begin
     LGs[sel].f_c(5) = 0
 endif
 sel = testFloat(LGs.f_c(6))
 if(sel(0) gt -1) then begin
     LGs[sel].f_c(6) = 0
 endif
 sel = testFloat(LGs.f_c(7))
 if(sel(0) gt -1) then begin
     LGs[sel].f_c(7) = 0
 endif
 sel = testFloat(LGs.f_c(8))
 if(sel(0) gt -1) then begin
     LGs[sel].f_c(8) = 0
 endif
 sel = testFloat(LGs.f_cmax(0))
 if(sel(0) gt -1) then begin
     LGs[sel].f_cmax(0) = 0
 endif
 sel = testFloat(LGs.f_cmax(1))
 if(sel(0) gt -1) then begin
     LGs[sel].f_cmax(1) = 0
 endif
 sel = testFloat(LGs.f_cmax(2))
 if(sel(0) gt -1) then begin
     LGs[sel].f_cmax(2) = 0
 endif
 sel = testFloat(LGs.f_cmax(3))
 if(sel(0) gt -1) then begin
     LGs[sel].f_cmax(3) = 0
 endif
 sel = testFloat(LGs.f_cmax(4))
 if(sel(0) gt -1) then begin
     LGs[sel].f_cmax(4) = 0
 endif
 sel = testFloat(LGs.f_cmax(5))
 if(sel(0) gt -1) then begin
     LGs[sel].f_cmax(5) = 0
 endif
 sel = testFloat(LGs.f_cmax(6))
 if(sel(0) gt -1) then begin
     LGs[sel].f_cmax(6) = 0
 endif
 sel = testFloat(LGs.f_cmax(7))
 if(sel(0) gt -1) then begin
     LGs[sel].f_cmax(7) = 0
 endif
 sel = testFloat(LGs.f_cmax(8))
 if(sel(0) gt -1) then begin
     LGs[sel].f_cmax(8) = 0
 endif
end
