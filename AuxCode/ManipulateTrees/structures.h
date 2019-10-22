// TODO add documentation, also for all fields
struct halo_data
{
	/* merger tree pointers */
	int Descendant;
	int FirstProgenitor;
	int NextProgenitor;
	int FirstHaloInFOFgroup;
	int NextHaloInFOFgroup;

  /* properties of halo */
	int Len;
	float M_Mean200, M_Crit200, M_TopHat;
	float Pos[3];
	float Vel[3];
	float VelDisp;
	float Vmax;
	float Spin[3];
	long long MostBoundID;

  /* original position in subfind output */
	int SnapNum;
	int FileNr;
	int SubhaloIndex;
	float SubHalfMass;
}
  *Halo, *Halo_Data;


struct  halo_ids_data
{
 long long HaloID;
 long long FileTreeNr;
 long long FirstProgenitor;
 long long LastProgenitor;
 long long NextProgenitor;
 long long Descendant;
 long long FirstHaloInFOFgroup;
 long long NextHaloInFOFgroup;
#ifdef MRII
 long long MainLeafID; 
#endif
 double    Redshift;
 int       PeanoKey;
 int       dummy;      /* need to use this padding for 64bit alignment */
} *HaloIDs, *HaloIDs_Data;



// TODO add documentation, also for all fields
struct halo_aux_data  /* auxiliary halo data */
{
	int DoneFlag;
	int HaloFlag;
	int NGalaxies;
	int FirstGalaxy;
	float M_Mean200_Unscaled;
	float M_Crit200_Unscaled;
	float Pos_Unscaled[3];
	float Vel_Unscaled[3];
	float Vmax_Unscaled;
	float Spin_Unscaled[3];
}
 *HaloAux;

