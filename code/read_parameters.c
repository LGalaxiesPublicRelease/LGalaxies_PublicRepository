#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#ifdef MCMC
#include "mcmc_vars.h"
#include "mcmc_proto.h"
#endif

#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

/**@file read_parameters.c reads all the parameters in input.par into global variables
 *       that can be used by the code. */


void read_parameter_file(char *fname)
{
  FILE *fd;
  char buf[400], buf1[400], buf2[400], buf3[400];
  int i, j, nt = 0;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int errorFlag = 0;

#ifdef PARALLEL
  if(ThisTask == 0)
    printf("\nreading parameter file:\n\n");
#else
  printf("\nreading parameter file:\n\n");
#endif

  strcpy(tag[nt], "OutputDir");
  addr[nt] = OutputDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileNameGalaxies");
  addr[nt] = FileNameGalaxies;
  id[nt++] = STRING;

  strcpy(tag[nt], "McFile");
  addr[nt] = McFile;
  id[nt++] = STRING;

  strcpy(tag[nt], "SimulationDir");
  addr[nt] = SimulationDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileWithOutputRedshifts");
  addr[nt] = FileWithOutputRedshifts;
  id[nt++] = STRING;

#ifdef SPECIFYFILENR
  strcpy(tag[nt], "FileNrDir");
  addr[nt] = FileNrDir;
  id[nt++] = STRING;
#endif

  strcpy(tag[nt], "SpecPhotDir");
  addr[nt] = SpecPhotDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "PhotPrefix");
  addr[nt] = PhotPrefix;
  id[nt++] = STRING;

  strcpy(tag[nt], "SpecPhotIMF");
  addr[nt] = SpecPhotIMF;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileWithFilterNames");
  addr[nt] = FileWithFilterNames;
  id[nt++] = STRING;

  strcpy(tag[nt], "CoolFunctionsDir");
  addr[nt] = CoolFunctionsDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "MaxMemSize");
  addr[nt] = &MaxMemSize;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Hashbits");
  addr[nt] = &Hashbits;
  id[nt++] = INT;

  /*strcpy(tag[nt], "NumWrittenInParallel");
  addr[nt] = &NumWrittenInParallel;
  id[nt++] = INT;*/

  strcpy(tag[nt], "MinGalOutputMass");
  addr[nt] = &MinGalOutputMass;
  id[nt++] = DOUBLE;

//Variables used in the MCMC
#ifdef MCMC
  strcpy(tag[nt], "MCMCStartingParFile");
  addr[nt] = MCMCStartingParFile;
  id[nt++] = STRING;

  strcpy(tag[nt], "MCMCParameterPriorsAndSwitches");
  addr[nt] = MCMCParameterPriorsAndSwitches;
  id[nt++] = STRING;

  strcpy(tag[nt], "MCMCObsConstraints");
  addr[nt] = MCMCObsConstraints;
  id[nt++] = STRING;

  strcpy(tag[nt], "MCMCWeightsObsConstraints");
  addr[nt] = MCMCWeightsObsConstraints;
  id[nt++] = STRING;

  strcpy(tag[nt], "CosmologyTablesDir");
  addr[nt] = CosmologyTablesDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "ObsConstraintsDir");
  addr[nt] = ObsConstraintsDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "MCMCSampleDir");
  addr[nt] = MCMCSampleDir;
  id[nt++] = STRING;

#ifdef MR_PLUS_MRII
  strcpy(tag[nt], "MCMCSampleFilePrefix_MR");
  addr[nt] = MCMCSampleFilePrefix_MR;
  id[nt++] = STRING;

  strcpy(tag[nt], "MCMCSampleFilePrefix_MRII");
  addr[nt] = MCMCSampleFilePrefix_MRII;
  id[nt++] = STRING;

  strcpy(tag[nt], "MCMCSampleFile_MR");
  addr[nt] = &MCMCSampleFile_MR;
  id[nt++] = INT;

  strcpy(tag[nt], "MCMCSampleFile_MRII");
  addr[nt] = &MCMCSampleFile_MRII;
  id[nt++] = INT;
#else
  strcpy(tag[nt], "MCMCSampleFilePrefix");
  addr[nt] = MCMCSampleFilePrefix;
  id[nt++] = STRING;

  strcpy(tag[nt], "MCMCSampleFile");
  addr[nt] = &MCMCSampleFile;
  id[nt++] = INT;
#endif

  strcpy(tag[nt], "MCMCTreeSampleFile");
  addr[nt] = &MCMCTreeSampleFile;
  id[nt++] = INT;

  strcpy(tag[nt], "MCMCHaloModelDir");
  addr[nt] = MCMCHaloModelDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "ChainLength");
  addr[nt] = &ChainLength;
  id[nt++] = INT;

  strcpy(tag[nt], "Sample_Physical_Parameters");
  addr[nt] = &Sample_Physical_Parameters;
  id[nt++] = INT;

  strcpy(tag[nt], "Time_Dependant_PhysPar");
  addr[nt] = &Time_Dependant_PhysPar;
  id[nt++] = INT;

  strcpy(tag[nt], "Sample_Cosmological_Parameters");
  addr[nt] = &Sample_Cosmological_Parameters;
  id[nt++] = INT;

  strcpy(tag[nt], "MCMCMode");
  addr[nt] = &MCMCMode;
  id[nt++] = INT;

  strcpy(tag[nt], "MCMC_LogStep_Size");
  addr[nt] = &MCMC_LogStep_Size;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "MCMC_Initial_Par_Displacement");
  addr[nt] = &MCMC_Initial_Par_Displacement;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "MCMC_Minimum_Obs_Error");
  addr[nt] = &MCMC_Minimum_Obs_Error;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "AddedErrOnMass");
  addr[nt] = &AddedErrOnMass;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "MachineTimeOut");
  addr[nt] = &MachineTimeOut;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "JobSubmitCommand");
  addr[nt] = JobSubmitCommand;
  id[nt++] = STRING;

  strcpy(tag[nt], "JobSubmitFile");
  addr[nt] = JobSubmitFile;
  id[nt++] = STRING;

  strcpy(tag[nt], "JobSubmitPipe");
  addr[nt] = JobSubmitPipe;
  id[nt++] = STRING;
#endif
 


   //Variables for the Scaling & Cosmological Parameters

  strcpy(tag[nt], "ScalePos");
  addr[nt] = &ScalePos;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ScaleMass");
  addr[nt] = &ScaleMass;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "BaryonFrac");
  addr[nt] = &BaryonFrac;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Sigma8");
  addr[nt] = &Sigma8;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Omega");
  addr[nt] = &Omega;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "OmegaLambda");
  addr[nt] = &OmegaLambda;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Hubble_h");
  addr[nt] = &Hubble_h;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Omega_OriginalCosm");
  addr[nt] = &Omega_OriginalCosm;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "OmegaLambda_OriginalCosm");
  addr[nt] = &OmegaLambda_OriginalCosm;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Hubble_h_OriginalCosm");
  addr[nt] = &Hubble_h_OriginalCosm;
  id[nt++] = DOUBLE;

#ifdef MR_PLUS_MRII  //OPTION for MCMC
  //MR
  strcpy(tag[nt], "FileWithZList_MR");
  addr[nt] = FileWithZList_MR;
  id[nt++] = STRING;

  strcpy(tag[nt], "PartMass_MR");
  addr[nt] = &PartMass_MR;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "BoxSize_MR");
  addr[nt] = &BoxSize_MR;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "FileWithZList_OriginalCosm_MR");
  addr[nt] = FileWithZList_OriginalCosm_MR;
  id[nt++] = STRING;

  strcpy(tag[nt], "PartMass_OriginalCosm_MR");
  addr[nt] = &PartMass_OriginalCosm_MR;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "BoxSize_OriginalCosm_MR");
  addr[nt] = &BoxSize_OriginalCosm_MR;
  id[nt++] = DOUBLE;

  //MRII
  strcpy(tag[nt], "FileWithZList_MRII");
  addr[nt] = FileWithZList_MRII;
  id[nt++] = STRING;

  strcpy(tag[nt], "PartMass_MRII");
  addr[nt] = &PartMass_MRII;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "BoxSize_MRII");
  addr[nt] = &BoxSize_MRII;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "FileWithZList_OriginalCosm_MRII");
  addr[nt] = FileWithZList_OriginalCosm_MRII;
  id[nt++] = STRING;

  strcpy(tag[nt], "PartMass_OriginalCosm_MRII");
  addr[nt] = &PartMass_OriginalCosm_MRII;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "BoxSize_OriginalCosm_MRII");
  addr[nt] = &BoxSize_OriginalCosm_MRII;
  id[nt++] = DOUBLE;
#else
  strcpy(tag[nt], "FileWithZList");
  addr[nt] = FileWithZList;
  id[nt++] = STRING;

  strcpy(tag[nt], "PartMass");
  addr[nt] = &PartMass;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "BoxSize");
  addr[nt] = &BoxSize;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "PartMass_OriginalCosm");
  addr[nt] = &PartMass_OriginalCosm;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "BoxSize_OriginalCosm");
  addr[nt] = &BoxSize_OriginalCosm;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "FileWithZList_OriginalCosm");
  addr[nt] = FileWithZList_OriginalCosm;
  id[nt++] = STRING;
#endif


#ifdef MR_PLUS_MRII  //OPTION for MCMC
  strcpy(tag[nt], "LastDarkMatterSnapShot_MR");
  addr[nt] = &LastDarkMatterSnapShot_MR;
  id[nt++] = INT;

  strcpy(tag[nt], "LastDarkMatterSnapShot_MRII");
  addr[nt] = &LastDarkMatterSnapShot_MRII;
  id[nt++] = INT;
#else
  strcpy(tag[nt], "LastDarkMatterSnapShot");
  addr[nt] = &LastDarkMatterSnapShot;
  id[nt++] = INT;
#endif

#ifndef MCMC
  strcpy(tag[nt], "FirstFile");
  addr[nt] = &FirstFile;
  id[nt++] = INT;

  strcpy(tag[nt], "LastFile");
  addr[nt] = &LastFile;
  id[nt++] = INT;
#endif



  //Physical Recipes

#ifdef EXCESS_MASS
  strcpy(tag[nt], "InfallModel");
  addr[nt] = &InfallModel;
  id[nt++] = INT;
#endif

  strcpy(tag[nt], "ReionizationModel");
  addr[nt] = &ReionizationModel;
  id[nt++] = INT;

  strcpy(tag[nt], "DiskRadiusModel");
  addr[nt] = &DiskRadiusModel;
  id[nt++] = INT;

  strcpy(tag[nt], "StarFormationModel");
  addr[nt] = &StarFormationModel;
  id[nt++] = INT;

  strcpy(tag[nt], "FeedbackReheatingModel");
  addr[nt] = &FeedbackReheatingModel;
  id[nt++] = INT;

  strcpy(tag[nt], "FeedbackEjectionModel");
  addr[nt] = &FeedbackEjectionModel;
  id[nt++] = INT;

  strcpy(tag[nt], "FeedbackEagleScaling");
  addr[nt] = &FeedbackEagleScaling;
  id[nt++] = INT;

  strcpy(tag[nt], "FeedbackReheatingDeansityScaling");
  addr[nt] = &FeedbackReheatingDeansityScaling;
  id[nt++] = INT;


  strcpy(tag[nt], "FateOfSatellitesGas");
  addr[nt] = &FateOfSatellitesGas;
  id[nt++] = INT;

  strcpy(tag[nt], "ReIncorporationModel");
  addr[nt] = &ReIncorporationModel;
  id[nt++] = INT;

  strcpy(tag[nt], "BlackHoleGrowth");
  addr[nt] = &BlackHoleGrowth;
  id[nt++] = INT;

  strcpy(tag[nt], "AGNRadioModeModel");
  addr[nt] = &AGNRadioModeModel;
  id[nt++] = INT;

  strcpy(tag[nt], "DiskInstabilityModel");
  addr[nt] = &DiskInstabilityModel;
  id[nt++] = INT;

  strcpy(tag[nt], "BHGrowthInDiskInstabilityModel");
  addr[nt] = &BHGrowthInDiskInstabilityModel;
  id[nt++] = INT;

  strcpy(tag[nt], "HotGasStripingModel");
  addr[nt] = &HotGasStripingModel;
  id[nt++] = INT;

  strcpy(tag[nt], "HotGasOnType2Galaxies");
  addr[nt] = &HotGasOnType2Galaxies;
  id[nt++] = INT;

  strcpy(tag[nt], "StarBurstModel");
  addr[nt] = &StarBurstModel;
  id[nt++] = INT;

  strcpy(tag[nt], "BulgeFormationInMinorMergersOn");
  addr[nt] = &BulgeFormationInMinorMergersOn;
  id[nt++] = INT;

  strcpy(tag[nt], "MetallicityOption");
  addr[nt] = &MetallicityOption;
  id[nt++] = INT;
  
#ifdef H2_AND_RINGS
  strcpy(tag[nt], "H2FractionRecipe");
  addr[nt] = &H2FractionRecipe;
  id[nt++] = INT;    
  
  strcpy(tag[nt], "SFRtdyn");
  addr[nt] = &SFRtdyn;
  id[nt++] = INT; 
#endif


  //Physical Parameters

  strcpy(tag[nt], "Reionization_z0");
  addr[nt] = &Reionization_z0;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Reionization_zr");
  addr[nt] = &Reionization_zr;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Yield");
  addr[nt] = &Yield;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "RecycleFraction");
  addr[nt] = &RecycleFraction;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ThreshMajorMerger");
  addr[nt] = &ThreshMajorMerger;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "MergerTimeMultiplier");
  addr[nt] = &MergerTimeMultiplier;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "RamPressureStrip_CutOffMass");
  addr[nt] = &RamPressureStrip_CutOffMass;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "RamPressureRadiusThreshold");
  addr[nt] = &RamPressureRadiusThreshold;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "SfrEfficiency");
  addr[nt] = &SfrEfficiency;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "SfrColdCrit");
  addr[nt] = &SfrColdCrit;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "SfrBurstEfficiency");
  addr[nt] = &SfrBurstEfficiency;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "SfrBurstSlope");
  addr[nt] = &SfrBurstSlope;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "AgnEfficiency");
  addr[nt] = &AgnEfficiency;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "BlackHoleGrowthRate");
  addr[nt] = &BlackHoleGrowthRate;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "BlackHoleDisruptGrowthRate");
  addr[nt] = &BlackHoleDisruptGrowthRate;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "BlackHoleSeedMass");
  addr[nt] = &BlackHoleSeedMass;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "BlackHoleAccretionRate");
  addr[nt] = &BlackHoleAccretionRate;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "BlackHoleCutoffVelocity");
  addr[nt] = &BlackHoleCutoffVelocity;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "FeedbackReheatingEpsilon");
  addr[nt] = &FeedbackReheatingEpsilon;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ReheatPreVelocity");
  addr[nt] = &ReheatPreVelocity;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ReheatSlope");
  addr[nt] = &ReheatSlope;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "FeedbackEjectionEfficiency");
  addr[nt] = &FeedbackEjectionEfficiency;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "EjectPreVelocity");
  addr[nt] = &EjectPreVelocity;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "EjectSlope");
  addr[nt] = &EjectSlope;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ReIncorporationFactor");
  addr[nt] = &ReIncorporationFactor;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ReincZpower");
  addr[nt] = &ReincZpower;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ReincVelocitypower");
  addr[nt] = &ReincVelocitypower;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "FracZSNIItoHot");
  addr[nt] = &FracZSNIItoHot;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "FracZSNIatoHot");
  addr[nt] = &FracZSNIatoHot;
  id[nt++] = DOUBLE;

//in the future ready in different energies for each type of SN
//#ifdef FEEDBACK_COUPLED_WITH_MASS_RETURN
//#else
  strcpy(tag[nt], "EnergySN");
  addr[nt] = &EnergySN;
  id[nt++] = DOUBLE;
//#endif

  strcpy(tag[nt], "EtaSN");
  addr[nt] = &EtaSN;
  id[nt++] = DOUBLE;

#ifdef H2_AND_RINGS
  strcpy(tag[nt], "Clumpingfactor");
  addr[nt] = &Clumpingfactor;
  id[nt++] = DOUBLE;  

  strcpy(tag[nt], "GasInflowVel");
  addr[nt] = &GasInflowVel;
  id[nt++] = DOUBLE;
#endif


  //UNITS

  strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
  addr[nt] = &UnitVelocity_in_cm_per_s;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "UnitLength_in_cm");
  addr[nt] = &UnitLength_in_cm;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "UnitMass_in_g");
  addr[nt] = &UnitMass_in_g;
  id[nt++] = DOUBLE;


  if((fd = fopen(fname, "r")))
    {
      while(!feof(fd))
        {
          *buf = 0;
          fgets(buf, 200, fd);
          if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
            continue;

          if(buf1[0] == '%')
            continue;

          for(i = 0, j = -1; i < nt; i++)
            if(strcmp(buf1, tag[i]) == 0)
              {
                j = i;
                tag[i][0] = 0;
                break;
              }

          if(j >= 0)
            {
#ifdef PARALLEL
              if(ThisTask == 0)
                printf("%35s\t%10s\n", buf1, buf2);
#else
              printf("%35s\t%10s\n", buf1, buf2);
#endif
              switch (id[j])
                {
                case DOUBLE:
                  *((double *) addr[j]) = atof(buf2);
                  break;
                case STRING:
                  strcpy(addr[j], buf2);
                  break;
                case INT:
                  *((int *) addr[j]) = atoi(buf2);
                  break;
                }
            }
          else
            {
              printf("Error in file %s:   Tag '%s' not allowed or multiple defined.\n", fname, buf1);
              errorFlag = 1;
            }
        }
      fclose(fd);

      i = strlen(OutputDir);
      if(i > 0)
        if(OutputDir[i - 1] != '/')
          strcat(OutputDir, "/");
    }
  else
    {
      printf("Parameter file %s not found.\n", fname);
      errorFlag = 1;
    }


  for(i = 0; i < nt; i++) {
    if(*tag[i]) {
      if (strcmp("MinGalOutputMass", tag[i]) == 0)
	MinGalOutputMass=0.;
#ifdef EXCESS_MASS
      else if (strcmp("InfallModel", tag[i]) == 0)
	InfallModel=0;
#endif
      else if (strcmp("FeedbackEagleScaling", tag[i]) == 0)
      	FeedbackEagleScaling=0;
      else if (strcmp("FeedbackReheatingDeansityScaling", tag[i]) == 0)
	FeedbackReheatingDeansityScaling=0;

      else {
	printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	errorFlag = 1;
      }
    }
  }

  if(errorFlag)
    terminate("parameterfile incorrect");

}
