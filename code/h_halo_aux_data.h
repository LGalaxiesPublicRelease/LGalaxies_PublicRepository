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

