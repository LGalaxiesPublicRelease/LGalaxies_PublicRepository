# numpy dtype for LGAL_GAL_STRUCT
import numpy
struct_dtype = numpy.dtype([
('Type',numpy.int32,1),
('HaloIndex',numpy.int32,1),
('SnapNum',numpy.int32,1),
('LookBackTimeToSnap',numpy.float32,1),
('CentralMvir',numpy.float32,1),
('CentralRvir',numpy.float32,1),
('DistanceToCentralGal',numpy.float32,3),
('Pos',numpy.float32,3),
('Vel',numpy.float32,3),
('Len',numpy.int32,1),
('Mvir',numpy.float32,1),
('Rvir',numpy.float32,1),
('Vvir',numpy.float32,1),
('Vmax',numpy.float32,1),
('ColdGasSpin',numpy.float32,3),
('DiskSpin',numpy.float32,3),
('InfallVmax',numpy.float32,1),
('InfallVmaxPeak',numpy.float32,1),
('InfallSnap',numpy.int32,1),
('InfallHotGas',numpy.float32,1),
('HotRadius',numpy.float32,1),
('OriMergTime',numpy.float32,1),
('MergTime',numpy.float32,1),
('ColdGas',numpy.float32,1),
('StellarMass',numpy.float32,1),
('DiskMass',numpy.float32,1),
('BulgeMass',numpy.float32,1),
('HotGas',numpy.float32,1),
('ReheatedGas',numpy.float32,1),
('EjectedMass',numpy.float32,1),
('BlackHoleMass',numpy.float32,1),
('ICM',numpy.float32,1),
('MassFromInSitu',numpy.float32,1),
('MassFromMergers',numpy.float32,1),
('MassFromBursts',numpy.float32,1),
('MetalsColdGas',numpy.float32,1),
('MetalsStellarMass',numpy.float32,1),
('MetalsDiskMass',numpy.float32,1),
('MetalsBulgeMass',numpy.float32,1),
('MetalsHotGas',numpy.float32,1),
('MetalsReheatedGas',numpy.float32,1),
('MetalsEjectedMass',numpy.float32,1),
('MetalsICM',numpy.float32,1),
('PrimordialAccretionRate',numpy.float32,1),
('CoolingRadius',numpy.float32,1),
('CoolingRate',numpy.float32,1),
('CoolingRate_beforeAGN',numpy.float32,1),
('QuasarAccretionRate',numpy.float32,1),
('RadioAccretionRate',numpy.float32,1),
('Sfr',numpy.float32,1),
('SfrBulge',numpy.float32,1),
('XrayLum',numpy.float32,1),
('BulgeSize',numpy.float32,1),
('DiskRadius',numpy.float32,1),
('ColdGasRadius',numpy.float32,1),
('StellarHalfMassRadius',numpy.float32,1),
('StellarHalfLightRadius',numpy.float32,1),
('CosInclination',numpy.float32,1),
('DisruptOn',numpy.int32,1),
('MergeOn',numpy.int32,1),
('MagDust',numpy.float32,40),
('Mag',numpy.float32,40),
('MagBulge',numpy.float32,40),
('MassWeightAge',numpy.float32,1),
('rbandWeightAge',numpy.float32,1),
('sfh_ibin',numpy.int32,1),
('sfh_numbins',numpy.int32,1),
('sfh_DiskMass',numpy.float32,SFH_NBIN),
('sfh_BulgeMass',numpy.float32,SFH_NBIN),
('sfh_ICM',numpy.float32,SFH_NBIN),
('sfh_MetalsDiskMass',numpy.float32,SFH_NBIN),
('sfh_MetalsBulgeMass',numpy.float32,SFH_NBIN),
('sfh_MetalsICM',numpy.float32,SFH_NBIN),
('ending','i4',0)
])
properties_used = {}
for el in struct_dtype.names:
	properties_used[el] = False
