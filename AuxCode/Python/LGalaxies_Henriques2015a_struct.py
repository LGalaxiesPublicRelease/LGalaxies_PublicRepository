# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np

LGalaxiesStruct = np.dtype([
('Type',np.int32,1),
('HaloIndex',np.int32,1),
('SnapNum',np.int32,1),
('LookBackTimeToSnap',np.float32,1),
('CentralMvir',np.float32,1),
('CentralRvir',np.float32,1),
('DistanceToCentralGal',np.float32,3),
('Pos',np.float32,3),
('Vel',np.float32,3),
('Len',np.int32,1),
('Mvir',np.float32,1),
('Rvir',np.float32,1),
('Vvir',np.float32,1),
('Vmax',np.float32,1),
('GasSpin',np.float32,3),
('StellarSpin',np.float32,3),
('InfallVmax',np.float32,1),
('InfallVmaxPeak',np.float32,1),
('InfallSnap',np.int32,1),
('InfallHotGas',np.float32,1),
('HotRadius',np.float32,1),
('OriMergTime',np.float32,1),
('MergTime',np.float32,1),
('ColdGas',np.float32,1),
('StellarMass',np.float32,1),
('BulgeMass',np.float32,1),
('DiskMass',np.float32,1),
('HotGas',np.float32,1),
('EjectedMass',np.float32,1),
('BlackHoleMass',np.float32,1),
('ICM',np.float32,1),
('MetalsColdGas',np.float32,1),
('MetalsStellarMass',np.float32,1),
('MetalsBulgeMass',np.float32,1),
('MetalsDiskMass',np.float32,1),
('MetalsHotGas',np.float32,1),
('MetalsEjectedMass',np.float32,1),
('MetalsICM',np.float32,1),
('PrimordialAccretionRate',np.float32,1),
('CoolingRadius',np.float32,1),
('CoolingRate',np.float32,1),
('CoolingRate_beforeAGN',np.float32,1),
('QuasarAccretionRate',np.float32,1),
('RadioAccretionRate',np.float32,1),
('Sfr',np.float32,1),
('SfrBulge',np.float32,1),
('XrayLum',np.float32,1),
('BulgeSize',np.float32,1),
('StellarDiskRadius',np.float32,1),
('GasDiskRadius',np.float32,1),
('CosInclination',np.float32,1),
('DisruptOn',np.int32,1),
('MergeOn',np.int32,1),
('MagDust',np.float32,40),
('Mag',np.float32,40),
('MagBulge',np.float32,40),
('MassWeightAge',np.float32,1),
('rBandWeightAge',np.float32,1),
('sfh_ibin',np.int32,1),
('sfh_numbins',np.int32,1),
('sfh_DiskMass',np.float32,20),
('sfh_BulgeMass',np.float32,20),
('sfh_ICM',np.float32,20),
('sfh_MetalsDiskMass',np.float32,20),
('sfh_MetalsBulgeMass',np.float32,20),
('sfh_MetalsICM',np.float32,20)
])

PropertiesToRead = {}
for ii in LGalaxiesStruct.names:
	PropertiesToRead[ii] = False
           
PropertiesToRead['Type'] = True
#PropertiesToRead['HaloIndex'] = True
PropertiesToRead['SnapNum'] = True
#PropertiesToRead['LookBackTimeToSnap'] = True
#PropertiesToRead['CentralMvir'] = True
#PropertiesToRead['CentralRvir'] = True
PropertiesToRead['DistanceToCentralGal'] = True
PropertiesToRead['Pos'] = True
PropertiesToRead['Vel'] = True
#PropertiesToRead['Len'] = True
PropertiesToRead['Mvir'] = True
PropertiesToRead['Rvir'] = True
#PropertiesToRead['Vvir'] = True
#PropertiesToRead['Vmax'] = True
#PropertiesToRead['GasSpin'] = True
#PropertiesToRead['StellarSpin'] = True
#PropertiesToRead['InfallVmax'] = True
#PropertiesToRead['InfallVmaxPeak'] = True
#PropertiesToRead['InfallSnap'] = True
#PropertiesToRead['InfallHotGas'] = True
#PropertiesToRead['HotRadius'] = True
#PropertiesToRead['OriMergTime'] = True
#PropertiesToRead['MergTime'] = True
PropertiesToRead['ColdGas'] = True
PropertiesToRead['StellarMass'] = True
PropertiesToRead['BulgeMass'] = True
PropertiesToRead['DiskMass'] = True
PropertiesToRead['HotGas'] = True
#PropertiesToRead['EjectedMass'] = True
PropertiesToRead['BlackHoleMass'] = True
#PropertiesToRead['ICM'] = True
PropertiesToRead['MetalsColdGas'] = True
PropertiesToRead['MetalsStellarMass'] = True
#PropertiesToRead['MetalsBulgeMass'] = True
#PropertiesToRead['MetalsDiskMass'] = True
#PropertiesToRead['MetalsHotGas'] = True
#PropertiesToRead['MetalsEjectedMass'] = True
#PropertiesToRead['MetalsICM'] = True
#PropertiesToRead['PrimordialAccretionRate'] = True
#PropertiesToRead['CoolingRadius'] = True
#PropertiesToRead['CoolingRate'] = True
#PropertiesToRead['CoolingRate_beforeAGN'] = True
#PropertiesToRead['QuasarAccretionRate'] = True
#PropertiesToRead['RadioAccretionRate'] = True
PropertiesToRead['Sfr'] = True
#PropertiesToRead['SfrBulge'] = True
#PropertiesToRead['XrayLum'] = True
PropertiesToRead['BulgeSize'] = True
PropertiesToRead['StellarDiskRadius'] = True
#PropertiesToRead['GasDiskRadius'] = True
#PropertiesToRead['CosInclination'] = True
#PropertiesToRead['DisruptOn'] = True
#PropertiesToRead['MergeOn'] = True
PropertiesToRead['MagDust'] = True
#PropertiesToRead['Mag'] = True
#PropertiesToRead['MagBulge'] = True
PropertiesToRead['MassWeightAge'] = True
PropertiesToRead['rBandWeightAge'] = True
#PropertiesToRead['sfh_ibin'] = True
#PropertiesToRead['sfh_numbins'] = True
#PropertiesToRead['sfh_DiskMass'] = True
#PropertiesToRead['sfh_BulgeMass'] = True
#PropertiesToRead['sfh_ICM'] = True
#PropertiesToRead['sfh_MetalsDiskMass'] = True
#PropertiesToRead['sfh_MetalsBulgeMass'] = True
#PropertiesToRead['sfh_MetalsICM'] = True
        

# <codecell>


