;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO LGalaxy_allElements__define
tmp = {LGalaxy_allElements $
, Type : 0L $ 
, HaloIndex : 0L $ 
, SnapNum : 0L $ 
, CentralMvir : 0.0 $ 
, Pos : fltarr(3) $ 
, Vel : fltarr(3) $ 
, Len : 0L $ 
, Mvir : 0.0 $ 
, Rvir : 0.0 $ 
, Vvir : 0.0 $ 
, Vmax : 0.0 $ 
, GasSpin : fltarr(3) $ 
, StellarSpin : fltarr(3) $ 
, InfallVmax : 0.0 $ 
, InfallSnap : 0L $ 
, HotRadius : 0.0 $ 
, OriMergTime : 0.0 $ 
, MergTime : 0.0 $ 
, DistanceToCentralGal : fltarr(3) $ 
, ColdGas : 0.0 $ 
, BulgeMass : 0.0 $ 
, DiskMass : 0.0 $ 
, HotGas : 0.0 $ 
, EjectedMass : 0.0 $ 
, BlackHoleMass : 0.0 $ 
, BlackHoleGas : 0.0 $ 
, ICM : 0.0 $ 
, MetalsColdGas : fltarr(3) $ 
, MetalsBulgeMass : fltarr(3) $ 
, MetalsDiskMass : fltarr(3) $ 
, MetalsHotGas : fltarr(3) $ 
, MetalsEjectedMass : fltarr(3) $ 
, MetalsICM : fltarr(3) $ 
, Sfr : 0.0 $ 
, SfrBulge : 0.0 $ 
, XrayLum : 0.0 $ 
, BulgeSize : 0.0 $ 
, StellarDiskRadius : 0.0 $ 
, GasDiskRadius : 0.0 $ 
, CosInclination : 0.0 $
, DisruptOn : 0L $ 
, MergeOn : 0L $ 
, CoolingRadius : 0.0 $ 
, QuasarAccretionRate : 0.0 $ 
, RadioAccretionRate : 0.0 $ 
, Mag : fltarr(5) $ 
, MagBulge : fltarr(5) $ 
, MagDust : fltarr(5) $ 
, MassWeightAge : 0.0 $ 
, MagICL : fltarr(5) $ 
, sfh_ibin : 0L $ 
, sfh_time : fltarr(20) $ 
, sfh_dt : fltarr(20) $ 
, sfh_DiskMass : fltarr(20) $ 
, sfh_BulgeMass : fltarr(20) $ 
, sfh_ICM : fltarr(20) $ 
, sfh_MetalsDiskMass : fltarr(3,20) $ 
, sfh_MetalsBulgeMass : fltarr(3,20) $ 
, sfh_MetalsICM : fltarr(3,20) $
, sfh_ElementsDiskMass : fltarr(11,20) $
, sfh_ElementsBulgeMass : fltarr(11,20) $
, sfh_ElementsICM : fltarr(11,20) $
, DiskMass_elements : fltarr(11) $
, BulgeMass_elements : fltarr(11) $
, ColdGas_elements : fltarr(11) $
, HotGas_elements : fltarr(11) $
, ICM_elements : fltarr(11) $
, EjectedMass_elements : fltarr(11) $
}
end