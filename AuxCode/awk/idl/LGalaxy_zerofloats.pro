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
 sel = testFloat(LGs.MetalsColdGas)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsColdGas = 0
 endif
 sel = testFloat(LGs.MetalsStellarMass)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsStellarMass = 0
 endif
 sel = testFloat(LGs.MetalsBulgeMass)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsBulgeMass = 0
 endif
 sel = testFloat(LGs.MetalsDiskMass)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsDiskMass = 0
 endif
 sel = testFloat(LGs.MetalsHotGas)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsHotGas = 0
 endif
 sel = testFloat(LGs.MetalsEjectedMass)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsEjectedMass = 0
 endif
 sel = testFloat(LGs.MetalsICM)
 if(sel(0) gt -1) then begin
     LGs[sel].MetalsICM = 0
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
 sel = testFloat(LGs.MagDust(2))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(2) = 0
 endif
 sel = testFloat(LGs.MagDust(3))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(3) = 0
 endif
 sel = testFloat(LGs.MagDust(4))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(4) = 0
 endif
 sel = testFloat(LGs.MagDust(5))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(5) = 0
 endif
 sel = testFloat(LGs.MagDust(6))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(6) = 0
 endif
 sel = testFloat(LGs.MagDust(7))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(7) = 0
 endif
 sel = testFloat(LGs.MagDust(8))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(8) = 0
 endif
 sel = testFloat(LGs.MagDust(9))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(9) = 0
 endif
 sel = testFloat(LGs.MagDust(10))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(10) = 0
 endif
 sel = testFloat(LGs.MagDust(11))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(11) = 0
 endif
 sel = testFloat(LGs.MagDust(12))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(12) = 0
 endif
 sel = testFloat(LGs.MagDust(13))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(13) = 0
 endif
 sel = testFloat(LGs.MagDust(14))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(14) = 0
 endif
 sel = testFloat(LGs.MagDust(15))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(15) = 0
 endif
 sel = testFloat(LGs.MagDust(16))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(16) = 0
 endif
 sel = testFloat(LGs.MagDust(17))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(17) = 0
 endif
 sel = testFloat(LGs.MagDust(18))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(18) = 0
 endif
 sel = testFloat(LGs.MagDust(19))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(19) = 0
 endif
 sel = testFloat(LGs.MagDust(20))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(20) = 0
 endif
 sel = testFloat(LGs.MagDust(21))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(21) = 0
 endif
 sel = testFloat(LGs.MagDust(22))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(22) = 0
 endif
 sel = testFloat(LGs.MagDust(23))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(23) = 0
 endif
 sel = testFloat(LGs.MagDust(24))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(24) = 0
 endif
 sel = testFloat(LGs.MagDust(25))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(25) = 0
 endif
 sel = testFloat(LGs.MagDust(26))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(26) = 0
 endif
 sel = testFloat(LGs.MagDust(27))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(27) = 0
 endif
 sel = testFloat(LGs.MagDust(28))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(28) = 0
 endif
 sel = testFloat(LGs.MagDust(29))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(29) = 0
 endif
 sel = testFloat(LGs.MagDust(30))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(30) = 0
 endif
 sel = testFloat(LGs.MagDust(31))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(31) = 0
 endif
 sel = testFloat(LGs.MagDust(32))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(32) = 0
 endif
 sel = testFloat(LGs.MagDust(33))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(33) = 0
 endif
 sel = testFloat(LGs.MagDust(34))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(34) = 0
 endif
 sel = testFloat(LGs.MagDust(35))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(35) = 0
 endif
 sel = testFloat(LGs.MagDust(36))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(36) = 0
 endif
 sel = testFloat(LGs.MagDust(37))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(37) = 0
 endif
 sel = testFloat(LGs.MagDust(38))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(38) = 0
 endif
 sel = testFloat(LGs.MagDust(39))
 if(sel(0) gt -1) then begin
     LGs[sel].MagDust(39) = 0
 endif
 sel = testFloat(LGs.Mag(0))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(0) = 0
 endif
 sel = testFloat(LGs.Mag(1))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(1) = 0
 endif
 sel = testFloat(LGs.Mag(2))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(2) = 0
 endif
 sel = testFloat(LGs.Mag(3))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(3) = 0
 endif
 sel = testFloat(LGs.Mag(4))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(4) = 0
 endif
 sel = testFloat(LGs.Mag(5))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(5) = 0
 endif
 sel = testFloat(LGs.Mag(6))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(6) = 0
 endif
 sel = testFloat(LGs.Mag(7))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(7) = 0
 endif
 sel = testFloat(LGs.Mag(8))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(8) = 0
 endif
 sel = testFloat(LGs.Mag(9))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(9) = 0
 endif
 sel = testFloat(LGs.Mag(10))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(10) = 0
 endif
 sel = testFloat(LGs.Mag(11))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(11) = 0
 endif
 sel = testFloat(LGs.Mag(12))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(12) = 0
 endif
 sel = testFloat(LGs.Mag(13))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(13) = 0
 endif
 sel = testFloat(LGs.Mag(14))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(14) = 0
 endif
 sel = testFloat(LGs.Mag(15))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(15) = 0
 endif
 sel = testFloat(LGs.Mag(16))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(16) = 0
 endif
 sel = testFloat(LGs.Mag(17))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(17) = 0
 endif
 sel = testFloat(LGs.Mag(18))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(18) = 0
 endif
 sel = testFloat(LGs.Mag(19))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(19) = 0
 endif
 sel = testFloat(LGs.Mag(20))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(20) = 0
 endif
 sel = testFloat(LGs.Mag(21))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(21) = 0
 endif
 sel = testFloat(LGs.Mag(22))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(22) = 0
 endif
 sel = testFloat(LGs.Mag(23))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(23) = 0
 endif
 sel = testFloat(LGs.Mag(24))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(24) = 0
 endif
 sel = testFloat(LGs.Mag(25))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(25) = 0
 endif
 sel = testFloat(LGs.Mag(26))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(26) = 0
 endif
 sel = testFloat(LGs.Mag(27))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(27) = 0
 endif
 sel = testFloat(LGs.Mag(28))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(28) = 0
 endif
 sel = testFloat(LGs.Mag(29))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(29) = 0
 endif
 sel = testFloat(LGs.Mag(30))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(30) = 0
 endif
 sel = testFloat(LGs.Mag(31))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(31) = 0
 endif
 sel = testFloat(LGs.Mag(32))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(32) = 0
 endif
 sel = testFloat(LGs.Mag(33))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(33) = 0
 endif
 sel = testFloat(LGs.Mag(34))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(34) = 0
 endif
 sel = testFloat(LGs.Mag(35))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(35) = 0
 endif
 sel = testFloat(LGs.Mag(36))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(36) = 0
 endif
 sel = testFloat(LGs.Mag(37))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(37) = 0
 endif
 sel = testFloat(LGs.Mag(38))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(38) = 0
 endif
 sel = testFloat(LGs.Mag(39))
 if(sel(0) gt -1) then begin
     LGs[sel].Mag(39) = 0
 endif
 sel = testFloat(LGs.MagBulge(0))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(0) = 0
 endif
 sel = testFloat(LGs.MagBulge(1))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(1) = 0
 endif
 sel = testFloat(LGs.MagBulge(2))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(2) = 0
 endif
 sel = testFloat(LGs.MagBulge(3))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(3) = 0
 endif
 sel = testFloat(LGs.MagBulge(4))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(4) = 0
 endif
 sel = testFloat(LGs.MagBulge(5))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(5) = 0
 endif
 sel = testFloat(LGs.MagBulge(6))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(6) = 0
 endif
 sel = testFloat(LGs.MagBulge(7))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(7) = 0
 endif
 sel = testFloat(LGs.MagBulge(8))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(8) = 0
 endif
 sel = testFloat(LGs.MagBulge(9))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(9) = 0
 endif
 sel = testFloat(LGs.MagBulge(10))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(10) = 0
 endif
 sel = testFloat(LGs.MagBulge(11))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(11) = 0
 endif
 sel = testFloat(LGs.MagBulge(12))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(12) = 0
 endif
 sel = testFloat(LGs.MagBulge(13))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(13) = 0
 endif
 sel = testFloat(LGs.MagBulge(14))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(14) = 0
 endif
 sel = testFloat(LGs.MagBulge(15))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(15) = 0
 endif
 sel = testFloat(LGs.MagBulge(16))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(16) = 0
 endif
 sel = testFloat(LGs.MagBulge(17))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(17) = 0
 endif
 sel = testFloat(LGs.MagBulge(18))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(18) = 0
 endif
 sel = testFloat(LGs.MagBulge(19))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(19) = 0
 endif
 sel = testFloat(LGs.MagBulge(20))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(20) = 0
 endif
 sel = testFloat(LGs.MagBulge(21))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(21) = 0
 endif
 sel = testFloat(LGs.MagBulge(22))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(22) = 0
 endif
 sel = testFloat(LGs.MagBulge(23))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(23) = 0
 endif
 sel = testFloat(LGs.MagBulge(24))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(24) = 0
 endif
 sel = testFloat(LGs.MagBulge(25))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(25) = 0
 endif
 sel = testFloat(LGs.MagBulge(26))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(26) = 0
 endif
 sel = testFloat(LGs.MagBulge(27))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(27) = 0
 endif
 sel = testFloat(LGs.MagBulge(28))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(28) = 0
 endif
 sel = testFloat(LGs.MagBulge(29))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(29) = 0
 endif
 sel = testFloat(LGs.MagBulge(30))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(30) = 0
 endif
 sel = testFloat(LGs.MagBulge(31))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(31) = 0
 endif
 sel = testFloat(LGs.MagBulge(32))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(32) = 0
 endif
 sel = testFloat(LGs.MagBulge(33))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(33) = 0
 endif
 sel = testFloat(LGs.MagBulge(34))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(34) = 0
 endif
 sel = testFloat(LGs.MagBulge(35))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(35) = 0
 endif
 sel = testFloat(LGs.MagBulge(36))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(36) = 0
 endif
 sel = testFloat(LGs.MagBulge(37))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(37) = 0
 endif
 sel = testFloat(LGs.MagBulge(38))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(38) = 0
 endif
 sel = testFloat(LGs.MagBulge(39))
 if(sel(0) gt -1) then begin
     LGs[sel].MagBulge(39) = 0
 endif
 sel = testFloat(LGs.MassWeightAge)
 if(sel(0) gt -1) then begin
     LGs[sel].MassWeightAge = 0
 endif
 sel = testFloat(LGs.rbandWeightAge)
 if(sel(0) gt -1) then begin
     LGs[sel].rbandWeightAge = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(0) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(1) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(2) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(3))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(3) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(4))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(4) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(5))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(5) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(6))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(6) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(7))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(7) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(8))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(8) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(9))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(9) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(10))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(10) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(11))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(11) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(12))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(12) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(13))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(13) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(14))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(14) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(15))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(15) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(16))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(16) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(17))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(17) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(18))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(18) = 0
 endif
 sel = testFloat(LGs.sfh_DiskMass(19))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_DiskMass(19) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(0) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(1) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(2) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(3))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(3) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(4))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(4) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(5))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(5) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(6))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(6) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(7))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(7) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(8))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(8) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(9))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(9) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(10))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(10) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(11))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(11) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(12))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(12) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(13))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(13) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(14))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(14) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(15))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(15) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(16))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(16) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(17))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(17) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(18))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(18) = 0
 endif
 sel = testFloat(LGs.sfh_BulgeMass(19))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_BulgeMass(19) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(0) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(1) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(2) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(3))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(3) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(4))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(4) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(5))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(5) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(6))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(6) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(7))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(7) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(8))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(8) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(9))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(9) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(10))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(10) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(11))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(11) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(12))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(12) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(13))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(13) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(14))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(14) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(15))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(15) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(16))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(16) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(17))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(17) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(18))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(18) = 0
 endif
 sel = testFloat(LGs.sfh_ICM(19))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_ICM(19) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(0) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(1) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(2) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(3))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(3) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(4))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(4) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(5))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(5) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(6))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(6) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(7))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(7) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(8))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(8) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(9))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(9) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(10))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(10) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(11))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(11) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(12))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(12) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(13))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(13) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(14))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(14) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(15))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(15) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(16))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(16) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(17))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(17) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(18))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(18) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsDiskMass(19))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsDiskMass(19) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(0) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(1) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(2) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(3))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(3) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(4))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(4) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(5))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(5) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(6))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(6) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(7))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(7) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(8))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(8) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(9))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(9) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(10))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(10) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(11))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(11) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(12))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(12) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(13))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(13) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(14))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(14) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(15))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(15) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(16))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(16) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(17))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(17) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(18))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(18) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsBulgeMass(19))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsBulgeMass(19) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(0))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(0) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(1))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(1) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(2))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(2) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(3))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(3) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(4))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(4) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(5))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(5) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(6))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(6) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(7))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(7) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(8))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(8) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(9))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(9) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(10))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(10) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(11))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(11) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(12))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(12) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(13))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(13) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(14))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(14) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(15))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(15) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(16))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(16) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(17))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(17) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(18))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(18) = 0
 endif
 sel = testFloat(LGs.sfh_MetalsICM(19))
 if(sel(0) gt -1) then begin
     LGs[sel].sfh_MetalsICM(19) = 0
 endif
end
