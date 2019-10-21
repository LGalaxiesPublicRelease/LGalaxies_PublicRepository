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
 sel = testFloat(LGs.HaloSpin(0))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloSpin[0] --- ', nstart+sel
     print, 'HaloSpin --- ', LGs[sel].HaloSpin(0)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HaloSpin(1))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloSpin[1] --- ', nstart+sel
     print, 'HaloSpin --- ', LGs[sel].HaloSpin(1)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.HaloSpin(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'HaloSpin[2] --- ', nstart+sel
     print, 'HaloSpin --- ', LGs[sel].HaloSpin(2)
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
 sel = testFloat(LGs.H2fraction)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'H2fraction --- ', nstart+sel
     print, 'H2fraction --- ', LGs[sel].H2fraction
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.StellarMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'StellarMass --- ', nstart+sel
     print, 'StellarMass --- ', LGs[sel].StellarMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.DiskMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DiskMass --- ', nstart+sel
     print, 'DiskMass --- ', LGs[sel].DiskMass
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.BulgeMass)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'BulgeMass --- ', nstart+sel
     print, 'BulgeMass --- ', LGs[sel].BulgeMass
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
 sel = testFloat(LGs.MassFromInSitu)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MassFromInSitu --- ', nstart+sel
     print, 'MassFromInSitu --- ', LGs[sel].MassFromInSitu
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MassFromMergers)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MassFromMergers --- ', nstart+sel
     print, 'MassFromMergers --- ', LGs[sel].MassFromMergers
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MassFromBursts)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MassFromBursts --- ', nstart+sel
     print, 'MassFromBursts --- ', LGs[sel].MassFromBursts
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
 sel = testFloat(LGs.DiskRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'DiskRadius --- ', nstart+sel
     print, 'DiskRadius --- ', LGs[sel].DiskRadius
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.ColdGasRadius)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'ColdGasRadius --- ', nstart+sel
     print, 'ColdGasRadius --- ', LGs[sel].ColdGasRadius
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
 sel = testFloat(LGs.MagDust(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[2] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[3] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[4] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(5))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[5] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(5)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(6))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[6] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(6)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(7))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[7] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(7)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(8))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[8] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(8)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(9))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[9] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(9)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(10))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[10] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(10)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(11))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[11] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(11)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(12))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[12] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(12)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(13))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[13] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(13)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(14))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[14] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(14)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(15))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[15] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(15)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(16))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[16] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(16)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(17))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[17] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(17)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(18))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[18] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(18)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(19))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[19] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(19)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(20))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[20] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(20)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(21))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[21] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(21)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(22))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[22] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(22)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(23))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[23] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(23)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(24))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[24] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(24)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(25))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[25] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(25)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(26))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[26] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(26)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(27))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[27] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(27)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(28))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[28] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(28)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(29))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[29] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(29)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(30))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[30] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(30)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(31))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[31] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(31)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(32))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[32] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(32)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(33))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[33] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(33)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(34))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[34] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(34)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(35))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[35] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(35)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(36))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[36] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(36)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(37))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[37] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(37)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(38))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[38] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(38)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagDust(39))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagDust[39] --- ', nstart+sel
     print, 'MagDust --- ', LGs[sel].MagDust(39)
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
 sel = testFloat(LGs.Mag(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[2] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[3] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[4] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(5))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[5] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(5)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(6))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[6] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(6)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(7))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[7] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(7)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(8))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[8] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(8)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(9))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[9] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(9)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(10))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[10] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(10)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(11))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[11] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(11)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(12))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[12] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(12)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(13))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[13] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(13)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(14))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[14] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(14)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(15))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[15] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(15)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(16))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[16] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(16)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(17))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[17] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(17)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(18))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[18] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(18)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(19))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[19] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(19)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(20))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[20] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(20)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(21))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[21] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(21)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(22))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[22] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(22)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(23))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[23] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(23)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(24))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[24] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(24)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(25))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[25] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(25)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(26))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[26] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(26)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(27))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[27] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(27)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(28))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[28] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(28)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(29))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[29] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(29)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(30))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[30] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(30)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(31))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[31] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(31)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(32))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[32] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(32)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(33))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[33] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(33)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(34))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[34] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(34)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(35))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[35] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(35)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(36))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[36] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(36)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(37))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[37] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(37)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(38))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[38] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(38)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.Mag(39))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'Mag[39] --- ', nstart+sel
     print, 'Mag --- ', LGs[sel].Mag(39)
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
 sel = testFloat(LGs.MagBulge(2))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[2] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(2)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(3))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[3] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(3)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(4))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[4] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(4)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(5))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[5] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(5)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(6))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[6] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(6)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(7))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[7] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(7)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(8))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[8] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(8)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(9))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[9] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(9)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(10))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[10] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(10)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(11))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[11] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(11)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(12))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[12] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(12)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(13))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[13] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(13)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(14))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[14] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(14)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(15))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[15] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(15)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(16))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[16] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(16)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(17))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[17] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(17)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(18))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[18] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(18)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(19))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[19] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(19)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(20))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[20] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(20)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(21))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[21] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(21)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(22))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[22] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(22)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(23))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[23] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(23)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(24))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[24] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(24)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(25))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[25] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(25)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(26))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[26] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(26)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(27))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[27] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(27)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(28))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[28] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(28)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(29))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[29] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(29)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(30))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[30] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(30)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(31))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[31] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(31)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(32))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[32] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(32)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(33))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[33] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(33)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(34))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[34] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(34)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(35))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[35] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(35)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(36))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[36] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(36)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(37))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[37] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(37)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(38))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[38] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(38)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MagBulge(39))
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MagBulge[39] --- ', nstart+sel
     print, 'MagBulge --- ', LGs[sel].MagBulge(39)
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.MassWeightAge)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'MassWeightAge --- ', nstart+sel
     print, 'MassWeightAge --- ', LGs[sel].MassWeightAge
     badranges=[badranges,sel]
 endif
 sel = testFloat(LGs.rbandWeightAge)
 if(sel(0) gt -1) then begin
     bad=1
     print, 'rbandWeightAge --- ', nstart+sel
     print, 'rbandWeightAge --- ', LGs[sel].rbandWeightAge
     badranges=[badranges,sel]
 endif
if(bad) then begin 
     print, 'badranges found: ',badranges
endif
return, badranges
end
