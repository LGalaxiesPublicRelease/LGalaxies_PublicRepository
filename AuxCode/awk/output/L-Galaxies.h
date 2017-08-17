struct LGalaxy {
     int Type;
      int HaloIndex;
      int SnapNum;
      float LookBackTimeToSnap;
      float CentralMvir;
      float CentralRvir;
      float DistanceToCentralGal[3];
      float Pos[3];
      float Vel[3];
      int Len;
      float Mvir;
      float Rvir;
      float Vvir;
      float Vmax;
      float ColdGasSpin[3];
      float DiskSpin[3];
      float InfallVmax;
      float InfallVmaxPeak;
      int InfallSnap;
      float InfallHotGas;
      float HotRadius;
      float OriMergTime;
      float MergTime;
      float ColdGas;
      float StellarMass;
      float DiskMass;
      float BulgeMass;
      float HotGas;
      float ReheatedGas;
      float EjectedMass;
      float BlackHoleMass;
           float ICM;
      float MassFromInSitu;
      float MassFromMergers;
      float MassFromBursts;
      float MetalsColdGas;
      float MetalsStellarMass;
      float MetalsDiskMass;
      float MetalsBulgeMass;
      float MetalsHotGas;
      float MetalsReheatedGas;
      float MetalsEjectedMass;
      float MetalsICM;
      float PrimordialAccretionRate;
      float CoolingRadius;
           float CoolingRate;
      float CoolingRate_beforeAGN;
      float QuasarAccretionRate;
      float RadioAccretionRate;
      float Sfr;
      float SfrBulge;
      float XrayLum;
      float BulgeSize;
      float DiskRadius;
      float ColdGasRadius;
      float StellarHalfMassRadius;
      float StellarHalfLightRadius;
      float CosInclination;
      int DisruptOn;
      int MergeOn;
      float MagDust[40];
      float Mag[40];
      float MagBulge[40];
      float MassWeightAge;
      float rbandWeightAge;
      int sfh_ibin;
      int sfh_numbins;
                float sfh_DiskMass[SFH_NBIN];
      float sfh_BulgeMass[SFH_NBIN];
      float sfh_ICM[SFH_NBIN];
      float sfh_MetalsDiskMass[SFH_NBIN];
      float sfh_MetalsBulgeMass[SFH_NBIN];
      float sfh_MetalsICM[SFH_NBIN];
};
struct MoMaFGalaxy {
     long long GalID;
 // // ID of the galaxy     short snapnum;
 // // snapnum of the galaxy, repeated here for faster lookups of times etc     short sfh_ibin;
 // //Index of highest bin currently in use     //    float sfh_time;
 // yr // Lookback time at the middle of bin.     //    float sfh_dt;
 // yr // time width of bin.     float sfh_DiskMass;
 // 1e10 Msun/h // SFH of disk     float sfh_BulgeMass;
 // 1e10 Msun/h // SFH of bulge     float sfh_ICM;
 // 1e10 Msun/h // SFH of ICM     float sfh_MetalsDiskMass;
 // 1e10 Msun/h // Metals locked up in stars in disk.     float sfh_MetalsBulgeMass;
 // 1e10 Msun/h //Metals locked up in stars in bulge.     float sfh_MetalsICM;
};