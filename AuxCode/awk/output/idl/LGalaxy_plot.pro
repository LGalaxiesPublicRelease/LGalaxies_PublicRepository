;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO LGalaxy_hist, LGs 
; plot for each field in the input LGalaxy struct a histogram
 print, 'plotting  Type'
   plotHist, LGs.Type,'Type','I'
 print, 'plotting  HaloIndex'
   plotHist, LGs.HaloIndex,'HaloIndex','I'
 print, 'plotting  SnapNum'
   plotHist, LGs.SnapNum,'SnapNum','I'
 print, 'plotting  LookBackTimeToSnap'
   plotHist, LGs.LookBackTimeToSnap,'LookBackTimeToSnap','F'
 print, 'plotting  CentralMvir'
   plotHist, LGs.CentralMvir,'CentralMvir','F'
 print, 'plotting  CentralRvir'
   plotHist, LGs.CentralRvir,'CentralRvir','F'
 print, 'plotting  DistanceToCentralGal[0]'
   plotHist, LGs.DistanceToCentralGal(0),'DistanceToCentralGal[0]','F'
 print, 'plotting  DistanceToCentralGal[1]'
   plotHist, LGs.DistanceToCentralGal(1),'DistanceToCentralGal[1]','F'
 print, 'plotting  DistanceToCentralGal[2]'
   plotHist, LGs.DistanceToCentralGal(2),'DistanceToCentralGal[2]','F'
 print, 'plotting  Pos[0]'
   plotHist, LGs.Pos(0),'Pos[0]','F'
 print, 'plotting  Pos[1]'
   plotHist, LGs.Pos(1),'Pos[1]','F'
 print, 'plotting  Pos[2]'
   plotHist, LGs.Pos(2),'Pos[2]','F'
 print, 'plotting  Vel[0]'
   plotHist, LGs.Vel(0),'Vel[0]','F'
 print, 'plotting  Vel[1]'
   plotHist, LGs.Vel(1),'Vel[1]','F'
 print, 'plotting  Vel[2]'
   plotHist, LGs.Vel(2),'Vel[2]','F'
 print, 'plotting  Len'
   plotHist, LGs.Len,'Len','I'
 print, 'plotting  Mvir'
   plotHist, LGs.Mvir,'Mvir','F'
 print, 'plotting  Rvir'
   plotHist, LGs.Rvir,'Rvir','F'
 print, 'plotting  Vvir'
   plotHist, LGs.Vvir,'Vvir','F'
 print, 'plotting  Vmax'
   plotHist, LGs.Vmax,'Vmax','F'
 print, 'plotting  HaloSpin[0]'
   plotHist, LGs.HaloSpin(0),'HaloSpin[0]','F'
 print, 'plotting  HaloSpin[1]'
   plotHist, LGs.HaloSpin(1),'HaloSpin[1]','F'
 print, 'plotting  HaloSpin[2]'
   plotHist, LGs.HaloSpin(2),'HaloSpin[2]','F'
 print, 'plotting  InfallVmax'
   plotHist, LGs.InfallVmax,'InfallVmax','F'
 print, 'plotting  InfallVmaxPeak'
   plotHist, LGs.InfallVmaxPeak,'InfallVmaxPeak','F'
 print, 'plotting  InfallSnap'
   plotHist, LGs.InfallSnap,'InfallSnap','I'
 print, 'plotting  InfallHotGas'
   plotHist, LGs.InfallHotGas,'InfallHotGas','F'
 print, 'plotting  HotRadius'
   plotHist, LGs.HotRadius,'HotRadius','F'
 print, 'plotting  OriMergTime'
   plotHist, LGs.OriMergTime,'OriMergTime','F'
 print, 'plotting  MergTime'
   plotHist, LGs.MergTime,'MergTime','F'
 print, 'plotting  ColdGas'
   plotHist, LGs.ColdGas,'ColdGas','F'
 print, 'plotting  H2fraction'
   plotHist, LGs.H2fraction,'H2fraction','F'
 print, 'plotting  StellarMass'
   plotHist, LGs.StellarMass,'StellarMass','F'
 print, 'plotting  DiskMass'
   plotHist, LGs.DiskMass,'DiskMass','F'
 print, 'plotting  BulgeMass'
   plotHist, LGs.BulgeMass,'BulgeMass','F'
 print, 'plotting  HotGas'
   plotHist, LGs.HotGas,'HotGas','F'
 print, 'plotting  EjectedMass'
   plotHist, LGs.EjectedMass,'EjectedMass','F'
 print, 'plotting  BlackHoleMass'
   plotHist, LGs.BlackHoleMass,'BlackHoleMass','F'
 print, 'plotting  ICM'
   plotHist, LGs.ICM,'ICM','F'
 print, 'plotting  MassFromInSitu'
   plotHist, LGs.MassFromInSitu,'MassFromInSitu','F'
 print, 'plotting  MassFromMergers'
   plotHist, LGs.MassFromMergers,'MassFromMergers','F'
 print, 'plotting  MassFromBursts'
   plotHist, LGs.MassFromBursts,'MassFromBursts','F'
 print, 'plotting  PrimordialAccretionRate'
   plotHist, LGs.PrimordialAccretionRate,'PrimordialAccretionRate','F'
 print, 'plotting  CoolingRadius'
   plotHist, LGs.CoolingRadius,'CoolingRadius','F'
 print, 'plotting  CoolingRate'
   plotHist, LGs.CoolingRate,'CoolingRate','F'
 print, 'plotting  CoolingRate_beforeAGN'
   plotHist, LGs.CoolingRate_beforeAGN,'CoolingRate_beforeAGN','F'
 print, 'plotting  QuasarAccretionRate'
   plotHist, LGs.QuasarAccretionRate,'QuasarAccretionRate','F'
 print, 'plotting  RadioAccretionRate'
   plotHist, LGs.RadioAccretionRate,'RadioAccretionRate','F'
 print, 'plotting  Sfr'
   plotHist, LGs.Sfr,'Sfr','F'
 print, 'plotting  SfrBulge'
   plotHist, LGs.SfrBulge,'SfrBulge','F'
 print, 'plotting  XrayLum'
   plotHist, LGs.XrayLum,'XrayLum','F'
 print, 'plotting  BulgeSize'
   plotHist, LGs.BulgeSize,'BulgeSize','F'
 print, 'plotting  DiskRadius'
   plotHist, LGs.DiskRadius,'DiskRadius','F'
 print, 'plotting  ColdGasRadius'
   plotHist, LGs.ColdGasRadius,'ColdGasRadius','F'
 print, 'plotting  StellarHalfMassRadius'
   plotHist, LGs.StellarHalfMassRadius,'StellarHalfMassRadius','F'
 print, 'plotting  StellarHalfLightRadius'
   plotHist, LGs.StellarHalfLightRadius,'StellarHalfLightRadius','F'
 print, 'plotting  CosInclination'
   plotHist, LGs.CosInclination,'CosInclination','F'
 print, 'plotting  DisruptOn'
   plotHist, LGs.DisruptOn,'DisruptOn','I'
 print, 'plotting  MergeOn'
   plotHist, LGs.MergeOn,'MergeOn','I'
 print, 'plotting  MagDust[0]'
   plotHist, LGs.MagDust(0),'MagDust[0]','F'
 print, 'plotting  MagDust[1]'
   plotHist, LGs.MagDust(1),'MagDust[1]','F'
 print, 'plotting  MagDust[2]'
   plotHist, LGs.MagDust(2),'MagDust[2]','F'
 print, 'plotting  MagDust[3]'
   plotHist, LGs.MagDust(3),'MagDust[3]','F'
 print, 'plotting  MagDust[4]'
   plotHist, LGs.MagDust(4),'MagDust[4]','F'
 print, 'plotting  MagDust[5]'
   plotHist, LGs.MagDust(5),'MagDust[5]','F'
 print, 'plotting  MagDust[6]'
   plotHist, LGs.MagDust(6),'MagDust[6]','F'
 print, 'plotting  MagDust[7]'
   plotHist, LGs.MagDust(7),'MagDust[7]','F'
 print, 'plotting  MagDust[8]'
   plotHist, LGs.MagDust(8),'MagDust[8]','F'
 print, 'plotting  MagDust[9]'
   plotHist, LGs.MagDust(9),'MagDust[9]','F'
 print, 'plotting  MagDust[10]'
   plotHist, LGs.MagDust(10),'MagDust[10]','F'
 print, 'plotting  MagDust[11]'
   plotHist, LGs.MagDust(11),'MagDust[11]','F'
 print, 'plotting  MagDust[12]'
   plotHist, LGs.MagDust(12),'MagDust[12]','F'
 print, 'plotting  MagDust[13]'
   plotHist, LGs.MagDust(13),'MagDust[13]','F'
 print, 'plotting  MagDust[14]'
   plotHist, LGs.MagDust(14),'MagDust[14]','F'
 print, 'plotting  MagDust[15]'
   plotHist, LGs.MagDust(15),'MagDust[15]','F'
 print, 'plotting  MagDust[16]'
   plotHist, LGs.MagDust(16),'MagDust[16]','F'
 print, 'plotting  MagDust[17]'
   plotHist, LGs.MagDust(17),'MagDust[17]','F'
 print, 'plotting  MagDust[18]'
   plotHist, LGs.MagDust(18),'MagDust[18]','F'
 print, 'plotting  MagDust[19]'
   plotHist, LGs.MagDust(19),'MagDust[19]','F'
 print, 'plotting  MagDust[20]'
   plotHist, LGs.MagDust(20),'MagDust[20]','F'
 print, 'plotting  MagDust[21]'
   plotHist, LGs.MagDust(21),'MagDust[21]','F'
 print, 'plotting  MagDust[22]'
   plotHist, LGs.MagDust(22),'MagDust[22]','F'
 print, 'plotting  MagDust[23]'
   plotHist, LGs.MagDust(23),'MagDust[23]','F'
 print, 'plotting  MagDust[24]'
   plotHist, LGs.MagDust(24),'MagDust[24]','F'
 print, 'plotting  MagDust[25]'
   plotHist, LGs.MagDust(25),'MagDust[25]','F'
 print, 'plotting  MagDust[26]'
   plotHist, LGs.MagDust(26),'MagDust[26]','F'
 print, 'plotting  MagDust[27]'
   plotHist, LGs.MagDust(27),'MagDust[27]','F'
 print, 'plotting  MagDust[28]'
   plotHist, LGs.MagDust(28),'MagDust[28]','F'
 print, 'plotting  MagDust[29]'
   plotHist, LGs.MagDust(29),'MagDust[29]','F'
 print, 'plotting  MagDust[30]'
   plotHist, LGs.MagDust(30),'MagDust[30]','F'
 print, 'plotting  MagDust[31]'
   plotHist, LGs.MagDust(31),'MagDust[31]','F'
 print, 'plotting  MagDust[32]'
   plotHist, LGs.MagDust(32),'MagDust[32]','F'
 print, 'plotting  MagDust[33]'
   plotHist, LGs.MagDust(33),'MagDust[33]','F'
 print, 'plotting  MagDust[34]'
   plotHist, LGs.MagDust(34),'MagDust[34]','F'
 print, 'plotting  MagDust[35]'
   plotHist, LGs.MagDust(35),'MagDust[35]','F'
 print, 'plotting  MagDust[36]'
   plotHist, LGs.MagDust(36),'MagDust[36]','F'
 print, 'plotting  MagDust[37]'
   plotHist, LGs.MagDust(37),'MagDust[37]','F'
 print, 'plotting  MagDust[38]'
   plotHist, LGs.MagDust(38),'MagDust[38]','F'
 print, 'plotting  MagDust[39]'
   plotHist, LGs.MagDust(39),'MagDust[39]','F'
 print, 'plotting  Mag[0]'
   plotHist, LGs.Mag(0),'Mag[0]','F'
 print, 'plotting  Mag[1]'
   plotHist, LGs.Mag(1),'Mag[1]','F'
 print, 'plotting  Mag[2]'
   plotHist, LGs.Mag(2),'Mag[2]','F'
 print, 'plotting  Mag[3]'
   plotHist, LGs.Mag(3),'Mag[3]','F'
 print, 'plotting  Mag[4]'
   plotHist, LGs.Mag(4),'Mag[4]','F'
 print, 'plotting  Mag[5]'
   plotHist, LGs.Mag(5),'Mag[5]','F'
 print, 'plotting  Mag[6]'
   plotHist, LGs.Mag(6),'Mag[6]','F'
 print, 'plotting  Mag[7]'
   plotHist, LGs.Mag(7),'Mag[7]','F'
 print, 'plotting  Mag[8]'
   plotHist, LGs.Mag(8),'Mag[8]','F'
 print, 'plotting  Mag[9]'
   plotHist, LGs.Mag(9),'Mag[9]','F'
 print, 'plotting  Mag[10]'
   plotHist, LGs.Mag(10),'Mag[10]','F'
 print, 'plotting  Mag[11]'
   plotHist, LGs.Mag(11),'Mag[11]','F'
 print, 'plotting  Mag[12]'
   plotHist, LGs.Mag(12),'Mag[12]','F'
 print, 'plotting  Mag[13]'
   plotHist, LGs.Mag(13),'Mag[13]','F'
 print, 'plotting  Mag[14]'
   plotHist, LGs.Mag(14),'Mag[14]','F'
 print, 'plotting  Mag[15]'
   plotHist, LGs.Mag(15),'Mag[15]','F'
 print, 'plotting  Mag[16]'
   plotHist, LGs.Mag(16),'Mag[16]','F'
 print, 'plotting  Mag[17]'
   plotHist, LGs.Mag(17),'Mag[17]','F'
 print, 'plotting  Mag[18]'
   plotHist, LGs.Mag(18),'Mag[18]','F'
 print, 'plotting  Mag[19]'
   plotHist, LGs.Mag(19),'Mag[19]','F'
 print, 'plotting  Mag[20]'
   plotHist, LGs.Mag(20),'Mag[20]','F'
 print, 'plotting  Mag[21]'
   plotHist, LGs.Mag(21),'Mag[21]','F'
 print, 'plotting  Mag[22]'
   plotHist, LGs.Mag(22),'Mag[22]','F'
 print, 'plotting  Mag[23]'
   plotHist, LGs.Mag(23),'Mag[23]','F'
 print, 'plotting  Mag[24]'
   plotHist, LGs.Mag(24),'Mag[24]','F'
 print, 'plotting  Mag[25]'
   plotHist, LGs.Mag(25),'Mag[25]','F'
 print, 'plotting  Mag[26]'
   plotHist, LGs.Mag(26),'Mag[26]','F'
 print, 'plotting  Mag[27]'
   plotHist, LGs.Mag(27),'Mag[27]','F'
 print, 'plotting  Mag[28]'
   plotHist, LGs.Mag(28),'Mag[28]','F'
 print, 'plotting  Mag[29]'
   plotHist, LGs.Mag(29),'Mag[29]','F'
 print, 'plotting  Mag[30]'
   plotHist, LGs.Mag(30),'Mag[30]','F'
 print, 'plotting  Mag[31]'
   plotHist, LGs.Mag(31),'Mag[31]','F'
 print, 'plotting  Mag[32]'
   plotHist, LGs.Mag(32),'Mag[32]','F'
 print, 'plotting  Mag[33]'
   plotHist, LGs.Mag(33),'Mag[33]','F'
 print, 'plotting  Mag[34]'
   plotHist, LGs.Mag(34),'Mag[34]','F'
 print, 'plotting  Mag[35]'
   plotHist, LGs.Mag(35),'Mag[35]','F'
 print, 'plotting  Mag[36]'
   plotHist, LGs.Mag(36),'Mag[36]','F'
 print, 'plotting  Mag[37]'
   plotHist, LGs.Mag(37),'Mag[37]','F'
 print, 'plotting  Mag[38]'
   plotHist, LGs.Mag(38),'Mag[38]','F'
 print, 'plotting  Mag[39]'
   plotHist, LGs.Mag(39),'Mag[39]','F'
 print, 'plotting  MagBulge[0]'
   plotHist, LGs.MagBulge(0),'MagBulge[0]','F'
 print, 'plotting  MagBulge[1]'
   plotHist, LGs.MagBulge(1),'MagBulge[1]','F'
 print, 'plotting  MagBulge[2]'
   plotHist, LGs.MagBulge(2),'MagBulge[2]','F'
 print, 'plotting  MagBulge[3]'
   plotHist, LGs.MagBulge(3),'MagBulge[3]','F'
 print, 'plotting  MagBulge[4]'
   plotHist, LGs.MagBulge(4),'MagBulge[4]','F'
 print, 'plotting  MagBulge[5]'
   plotHist, LGs.MagBulge(5),'MagBulge[5]','F'
 print, 'plotting  MagBulge[6]'
   plotHist, LGs.MagBulge(6),'MagBulge[6]','F'
 print, 'plotting  MagBulge[7]'
   plotHist, LGs.MagBulge(7),'MagBulge[7]','F'
 print, 'plotting  MagBulge[8]'
   plotHist, LGs.MagBulge(8),'MagBulge[8]','F'
 print, 'plotting  MagBulge[9]'
   plotHist, LGs.MagBulge(9),'MagBulge[9]','F'
 print, 'plotting  MagBulge[10]'
   plotHist, LGs.MagBulge(10),'MagBulge[10]','F'
 print, 'plotting  MagBulge[11]'
   plotHist, LGs.MagBulge(11),'MagBulge[11]','F'
 print, 'plotting  MagBulge[12]'
   plotHist, LGs.MagBulge(12),'MagBulge[12]','F'
 print, 'plotting  MagBulge[13]'
   plotHist, LGs.MagBulge(13),'MagBulge[13]','F'
 print, 'plotting  MagBulge[14]'
   plotHist, LGs.MagBulge(14),'MagBulge[14]','F'
 print, 'plotting  MagBulge[15]'
   plotHist, LGs.MagBulge(15),'MagBulge[15]','F'
 print, 'plotting  MagBulge[16]'
   plotHist, LGs.MagBulge(16),'MagBulge[16]','F'
 print, 'plotting  MagBulge[17]'
   plotHist, LGs.MagBulge(17),'MagBulge[17]','F'
 print, 'plotting  MagBulge[18]'
   plotHist, LGs.MagBulge(18),'MagBulge[18]','F'
 print, 'plotting  MagBulge[19]'
   plotHist, LGs.MagBulge(19),'MagBulge[19]','F'
 print, 'plotting  MagBulge[20]'
   plotHist, LGs.MagBulge(20),'MagBulge[20]','F'
 print, 'plotting  MagBulge[21]'
   plotHist, LGs.MagBulge(21),'MagBulge[21]','F'
 print, 'plotting  MagBulge[22]'
   plotHist, LGs.MagBulge(22),'MagBulge[22]','F'
 print, 'plotting  MagBulge[23]'
   plotHist, LGs.MagBulge(23),'MagBulge[23]','F'
 print, 'plotting  MagBulge[24]'
   plotHist, LGs.MagBulge(24),'MagBulge[24]','F'
 print, 'plotting  MagBulge[25]'
   plotHist, LGs.MagBulge(25),'MagBulge[25]','F'
 print, 'plotting  MagBulge[26]'
   plotHist, LGs.MagBulge(26),'MagBulge[26]','F'
 print, 'plotting  MagBulge[27]'
   plotHist, LGs.MagBulge(27),'MagBulge[27]','F'
 print, 'plotting  MagBulge[28]'
   plotHist, LGs.MagBulge(28),'MagBulge[28]','F'
 print, 'plotting  MagBulge[29]'
   plotHist, LGs.MagBulge(29),'MagBulge[29]','F'
 print, 'plotting  MagBulge[30]'
   plotHist, LGs.MagBulge(30),'MagBulge[30]','F'
 print, 'plotting  MagBulge[31]'
   plotHist, LGs.MagBulge(31),'MagBulge[31]','F'
 print, 'plotting  MagBulge[32]'
   plotHist, LGs.MagBulge(32),'MagBulge[32]','F'
 print, 'plotting  MagBulge[33]'
   plotHist, LGs.MagBulge(33),'MagBulge[33]','F'
 print, 'plotting  MagBulge[34]'
   plotHist, LGs.MagBulge(34),'MagBulge[34]','F'
 print, 'plotting  MagBulge[35]'
   plotHist, LGs.MagBulge(35),'MagBulge[35]','F'
 print, 'plotting  MagBulge[36]'
   plotHist, LGs.MagBulge(36),'MagBulge[36]','F'
 print, 'plotting  MagBulge[37]'
   plotHist, LGs.MagBulge(37),'MagBulge[37]','F'
 print, 'plotting  MagBulge[38]'
   plotHist, LGs.MagBulge(38),'MagBulge[38]','F'
 print, 'plotting  MagBulge[39]'
   plotHist, LGs.MagBulge(39),'MagBulge[39]','F'
 print, 'plotting  MassWeightAge'
   plotHist, LGs.MassWeightAge,'MassWeightAge','F'
 print, 'plotting  rbandWeightAge'
   plotHist, LGs.rbandWeightAge,'rbandWeightAge','F'
 print, 'plotting  sfh_ibin'
   plotHist, LGs.sfh_ibin,'sfh_ibin','I'
 print, 'plotting  sfh_numbins'
   plotHist, LGs.sfh_numbins,'sfh_numbins','I'
end
