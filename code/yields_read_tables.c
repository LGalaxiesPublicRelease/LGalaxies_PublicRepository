#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

void read_yield_tables(void)
{
	//------------------------------------------
	//READ LIFETIME MASS LIST:
	//------------------------------------------
	FILE *fd1;
	char buf1[100];
	int i1;
	float m1;
	static char *name1 = "stripped_interp_LifetimeMasses.txt";

	sprintf(buf1, "./YieldTables/%s", name1);

	if(!(fd1 = fopen(buf1, "r")))
        {
          printf("file `%s' not found.\n", buf1);
          exit(0);
        }

	for(i1=0; i1<LIFETIME_MASS_NUM; i1++)
        {
	  fscanf(fd1, "%f", &m1);
	  lifetimeMasses[i1] = m1;
        }
	fclose(fd1);
	//printf("Lifetime masses read.\n");

	//------------------------------------------
	//READ LIFETIME METALLICITY LIST:
	//------------------------------------------
	FILE *fd2;
	char buf2[100];
	int i2;
	float m2;
	static char *name2 = "stripped_LifetimeMetallicities.txt";

	sprintf(buf2, "./YieldTables/%s", name2);

	if(!(fd2 = fopen(buf2, "r")))
        {
          printf("file `%s' not found.\n", buf2);
          exit(0);
        }

	for(i2=0; i2<6; i2++)
        {
	  fscanf(fd2, "%f", &m2);
	  lifetimeMetallicities[i2] = m2;
        }
	fclose(fd2);
	//printf("Lifetime metallicities read.\n");

	//------------------------------------------
	//READ LIFETIME TABLE:
	//------------------------------------------
	FILE *fd3;
	char buf3[100];
	int i3,j3;
	float m3; //,mass;
	//float get_mass(float time, float Z0);
	static char *name3 = "stripped_interp_Lifetimes.txt";
	sprintf(buf3, "./YieldTables/%s", name3);

	if(!(fd3 = fopen(buf3, "r")))
        {
          printf("file `%s' not found.\n", buf3);
          exit(0);
        }

	for(i3=0; i3<LIFETIME_Z_NUM; i3++)
        {
          for(j3=0; j3<LIFETIME_MASS_NUM; j3++)
          {
        	  fscanf(fd3, "%f", &m3);
        	  lifetimes[i3][j3] = m3;
          }
        }
	fclose(fd3);
	//printf("Lifetimes read.\n");
	printf("Lifetime tables read.\n");

	//------------------------------------------
	//READ AGB MASS LIST:
	//------------------------------------------
	FILE *fd4;
	char buf4[100];
	int i4;
	float m4;
	static char *name4 = "stripped_interp_AGBMasses.txt";

	sprintf(buf4, "./YieldTables/%s", name4);

	if(!(fd4 = fopen(buf4, "r")))
        {
          printf("file `%s' not found.\n", buf4);
          exit(0);
        }

	for(i4=0; i4<AGB_MASS_NUM; i4++)
        {
	  fscanf(fd4, "%f", &m4);
	  AGBMasses[i4] = m4;
        }
	fclose(fd4);

	//------------------------------------------
	//READ AGB METALLICITY LIST:
	//------------------------------------------
	FILE *fd5;
	char buf5[100];
	int i5;
	float m5;
	static char *name5 = "stripped_AGBMetallicities.txt";

	sprintf(buf5, "./YieldTables/%s", name5);

	if(!(fd5 = fopen(buf5, "r")))
        {
          printf("file `%s' not found.\n", buf5);
          exit(0);
        }

	for(i5=0; i5<AGB_Z_NUM; i5++)
        {
	  fscanf(fd5, "%f", &m5);
	  AGBMetallicities[i5] = m5;
        }
	fclose(fd5);

	//------------------------------------------
	//READ AGB EJECTED MASS LISTS:
	//------------------------------------------
	FILE *fd6;
	char buf6[100];
	int i6,j6;
	float m6;
	static char *name6[] = {
  	"convolved_stripped_interp_AGB_Z004_EjectedMasses.txt",
  	"convolved_stripped_interp_AGB_Z008_EjectedMasses.txt",
  	"convolved_stripped_interp_AGB_Z019_EjectedMasses.txt"
	};

	for(i6 = 0; i6 < AGB_Z_NUM; i6++)
    {
		sprintf(buf6, "./YieldTables/%s", name6[i6]);

		if(!(fd6 = fopen(buf6, "r")))
        {
            printf("file `%s' not found.\n", buf6);
            exit(0);
        }

        for(j6=0; j6<AGB_MASS_NUM; j6++)
        {
        	  fscanf(fd6, "%f", &m6);
        	  AGBEjectedMasses[i6][j6] = m6 * Chabrier_IMF(AGBMasses[j6]);
        }
        fclose(fd6);
	 }

	//------------------------------------------
	//READ AGB TOTAL METALS LISTS: This is the sum of all the yields NOT including H and He:
	//------------------------------------------
	FILE *fd7;
	char buf7[100];
	int i7,j7;
	float m7;
	static char *name7[] = {
  	"convolved_stripped_interp_AGB_Z004_TotalMetals.txt",
  	"convolved_stripped_interp_AGB_Z008_TotalMetals.txt",
  	"convolved_stripped_interp_AGB_Z019_TotalMetals.txt"
	};

	for(i7 = 0; i7 < AGB_Z_NUM; i7++)
    	{
	  sprintf(buf7, "./YieldTables/%s", name7[i7]);

	  if(!(fd7 = fopen(buf7, "r")))
          {
            printf("file `%s' not found.\n", buf7);
            exit(0);
          }

            for(j7=0; j7<AGB_MASS_NUM; j7++)
	    {
	      fscanf(fd7, "%f", &m7);
	      AGBTotalMetals[i7][j7] = m7 * Chabrier_IMF(AGBMasses[j7]);
	    }
	  fclose(fd7);
	}

	//------------------------------------------
	//READ AGB YIELD TABLES:
	//------------------------------------------
	FILE *fd8;
	char buf8[100];
	int k8,i8,j8;
	float m8;
	static char *name8[] = {
  	"convolved_stripped_interp_AGB_Z004_Yields.txt",
  	"convolved_stripped_interp_AGB_Z008_Yields.txt",
  	"convolved_stripped_interp_AGB_Z019_Yields.txt"
	};

	for(k8 = 0; k8 < AGB_Z_NUM; k8++)
    {
	  sprintf(buf8, "./YieldTables/%s", name8[k8]);

	  if(!(fd8 = fopen(buf8, "r")))
          {
            printf("file `%s' not found.\n", buf8);
            exit(0);
          }

	  for(i8=0; i8<11; i8++) //Number of element species (inc H and He) = 11
          {
            for(j8=0; j8<AGB_MASS_NUM; j8++)
            {
            	fscanf(fd8, "%f", &m8);
            	AGBYields[k8][i8][j8] = m8 * Chabrier_IMF(AGBMasses[j8]);
            }
          }
    }

	printf("AGB yield tables read.\n");

	//------------------------------------------
	//READ SN-II MASS LIST:
	//------------------------------------------
	FILE *fd9;
	char buf9[100];
	int i9;
	float m9;
#ifdef PORTINARI
	static char *name9 = "stripped_interp_SNIIMasses.txt";
#endif
#ifdef CHIEFFI
	static char *name9 = "stripped_interp_CL_SNIIMasses.txt";
#endif

	sprintf(buf9, "./YieldTables/%s", name9);

	if(!(fd9 = fopen(buf9, "r")))
        {
          printf("file `%s' not found.\n", buf9);
          exit(0);
        }

	for(i9=0; i9<SNII_MASS_NUM; i9++)
        {
	  fscanf(fd9, "%f", &m9);
	  SNIIMasses[i9] = m9;
        }
	fclose(fd9);

	//------------------------------------------
	//READ SN-II METALLICITY LIST:
	//------------------------------------------
	FILE *fd10;
	char buf10[100];
	int i10;
	float m10;
#ifdef PORTINARI
	static char *name10 = "stripped_SNIIMetallicities.txt";
#endif
#ifdef CHIEFFI
	static char *name10 = "stripped_CL_SNIIMetallicities.txt";
#endif

	sprintf(buf10, "./YieldTables/%s", name10);

	if(!(fd10 = fopen(buf10, "r")))
        {
          printf("file `%s' not found.\n", buf10);
          exit(0);
        }

	for(i10=0; i10<SNII_Z_NUM; i10++)
        {
			fscanf(fd10, "%f", &m10);
			SNIIMetallicities[i10] = m10;
        }
	fclose(fd10);

	//------------------------------------------
	//READ SN-II EJECTED MASS LISTS:
	//------------------------------------------
	FILE *fd11;
	char buf11[100];
	int i11,j11;
	float m11;
#ifdef PORTINARI
	static char *name11[] = {
  	"convolved_stripped_interp_SNII_Z0004_EjectedMasses.txt",
	"convolved_stripped_interp_SNII_Z004_EjectedMasses.txt",
  	"convolved_stripped_interp_SNII_Z008_EjectedMasses.txt",
  	"convolved_stripped_interp_SNII_Z02_EjectedMasses.txt",
	"convolved_stripped_interp_SNII_Z05_EjectedMasses.txt"
	};
#endif
#ifdef CHIEFFI
	static char *name11[] = {
  	"convolved_stripped_interp_CL_SNII_Z0_EjectedMasses.txt",
	"convolved_stripped_interp_CL_SNII_Z000001_EjectedMasses.txt",
  	"convolved_stripped_interp_CL_SNII_Z0001_EjectedMasses.txt",
  	"convolved_stripped_interp_CL_SNII_Z001_EjectedMasses.txt",
	"convolved_stripped_interp_CL_SNII_Z006_EjectedMasses.txt",
  	"convolved_stripped_interp_CL_SNII_Z02_EjectedMasses.txt"
	};
#endif

	for(i11 = 0; i11 < SNII_Z_NUM; i11++)
    	{
	  sprintf(buf11, "./YieldTables/%s", name11[i11]);

	  if(!(fd11 = fopen(buf11, "r")))
          {
            printf("file `%s' not found.\n", buf11);
            exit(0);
          }
          for(j11=0; j11<SNII_MASS_NUM; j11++)
	  {
	    fscanf(fd11, "%f", &m11);
	    SNIIEjectedMasses[i11][j11] = m11 * Chabrier_IMF(SNIIMasses[j11]);
	  }
	  fclose(fd11);
	}

	//------------------------------------------
	//READ SN-II TOTAL METALS LISTS:
	//------------------------------------------
	FILE *fd12;
	char buf12[100];
	int i12,j12;
	float m12;
#ifdef PORTINARI
	static char *name12[] = {
  	"convolved_stripped_interp_SNII_Z0004_TotalMetals.txt",
	"convolved_stripped_interp_SNII_Z004_TotalMetals.txt",
  	"convolved_stripped_interp_SNII_Z008_TotalMetals.txt",
  	"convolved_stripped_interp_SNII_Z02_TotalMetals.txt",
	"convolved_stripped_interp_SNII_Z05_TotalMetals.txt"
	};
#endif
#ifdef CHIEFFI
	static char *name12[] = {
  	"convolved_stripped_interp_CL_SNII_Z0_TotalMetals.txt",
	"convolved_stripped_interp_CL_SNII_Z000001_TotalMetals.txt",
  	"convolved_stripped_interp_CL_SNII_Z0001_TotalMetals.txt",
  	"convolved_stripped_interp_CL_SNII_Z001_TotalMetals.txt",
	"convolved_stripped_interp_CL_SNII_Z006_TotalMetals.txt",
  	"convolved_stripped_interp_CL_SNII_Z02_TotalMetals.txt"
	};
#endif

	for(i12 = 0; i12 < SNII_Z_NUM; i12++)
    {
	  sprintf(buf12, "./YieldTables/%s", name12[i12]);

	  if(!(fd12 = fopen(buf12, "r")))
      {
            printf("file `%s' not found.\n", buf12);
            exit(0);
      }

      for(j12=0; j12<SNII_MASS_NUM; j12++)
	  {
	    fscanf(fd12, "%f", &m12);
	    SNIITotalMetals[i12][j12] = m12 * Chabrier_IMF(SNIIMasses[j12]);
	  }
	  fclose(fd12);
	}

	//------------------------------------------
	//READ SN-II YIELD TABLES:
	//------------------------------------------
	FILE *fd13;
	char buf13[100];
	int k13,i13,j13;
	float m13;
#ifdef PORTINARI
	static char *name13[] = {
  	"convolved_stripped_interp_SNII_Z0004_Yields.txt",
	"convolved_stripped_interp_SNII_Z004_Yields.txt",
  	"convolved_stripped_interp_SNII_Z008_Yields.txt",
  	"convolved_stripped_interp_SNII_Z02_Yields.txt",
	"convolved_stripped_interp_SNII_Z05_Yields.txt"
	};
#endif
#ifdef CHIEFFI
	static char *name13[] = {
  	"convolved_stripped_interp_CL_SNII_Z0_Yields.txt",
	"convolved_stripped_interp_CL_SNII_Z000001_Yields.txt",
  	"convolved_stripped_interp_CL_SNII_Z0001_Yields.txt",
  	"convolved_stripped_interp_CL_SNII_Z001_Yields.txt",
	"convolved_stripped_interp_CL_SNII_Z006_Yields.txt",
  	"convolved_stripped_interp_CL_SNII_Z02_Yields.txt"
	};
#endif

	for(k13 = 0; k13 < SNII_Z_NUM; k13++)
    	{
	  sprintf(buf13, "./YieldTables/%s", name13[k13]);

	  if(!(fd13 = fopen(buf13, "r")))
          {
            printf("file `%s' not found.\n", buf13);
            exit(0);
          }

	  for(i13=0; i13<11; i13++) //Number of element species (inc H and He) = 11
          {
            for(j13=0; j13<SNII_MASS_NUM; j13++) //Number of initial masses = 85 (for Portinari yields) (11 in initial yield tables)
	    {
	      fscanf(fd13, "%f", &m13);
	      SNIIYields[k13][i13][j13] = m13 * Chabrier_IMF(SNIIMasses[j13]);
	      /*if (i13 == 4 && k13 >= SNII_Z_NUM-2) {SNIIYields[k13][i13][j13] = 4.0 * m13 * Chabrier_IMF(SNIIMasses[j13]);} //Just a test!: Quadruple the oxygen production at high Z.
	      else {SNIIYields[k13][i13][j13] = m13 * Chabrier_IMF(SNIIMasses[j13]);}*/
	      /*if (i13 == 3 && k13 >= SNII_Z_NUM-2) {SNIIYields[k13][i13][j13] = 1.5 * m13 * Chabrier_IMF(SNIIMasses[j13]);} //Just a test!: 1.5 x nitrogen production at high Z.
	      else {SNIIYields[k13][i13][j13] = m13 * Chabrier_IMF(SNIIMasses[j13]);}*/
	      /*//Just a test!: Undo R. Wiersma's corrections to the P98 yields:
	      if (i13 == 2) {SNIIYields[k13][i13][j13] = (m13*2.0) * Chabrier_IMF(SNIIMasses[j13]);}
	      else if (i13 == 6) {SNIIYields[k13][i13][j13] = (m13*0.5) * Chabrier_IMF(SNIIMasses[j13]);}
	      else if (i13 == 10) {SNIIYields[k13][i13][j13] = (m13*2.0) * Chabrier_IMF(SNIIMasses[j13]);}
	      else {SNIIYields[k13][i13][j13] = m13 * Chabrier_IMF(SNIIMasses[j13]);}*/
	    }
	  }
	}

	printf("SN-II yield tables read.\n");

#ifndef DTD
	//------------------------------------------
	//READ SN-Ia MASS LIST:
	//------------------------------------------
	FILE *fd17;
	char buf17[100];
	int i17;
	float m17;
	static char *name17 = "stripped_interp_SNIaMasses.txt";

	sprintf(buf17, "./YieldTables/%s", name17);

	if(!(fd17 = fopen(buf17, "r")))
        {
          printf("file `%s' not found.\n", buf17);
          exit(0);
        }

	for(i17=0; i17<SNIA_MASS_NUM; i17++)
        {
	  fscanf(fd17, "%f", &m17);
	  SNIaMasses[i17] = m17;
        }
	fclose(fd17);

	//------------------------------------------
	//READ SN-Ia TOTAL METALS:
	//------------------------------------------
	FILE *fd15;
	char buf15[100];
	int i15;
	float m15;
	static char *name15 = "convolved_stripped_interp_SNIa_TotalMetals.txt";

	sprintf(buf15, "./YieldTables/%s", name15);

	 if(!(fd15 = fopen(buf15, "r")))
     {
        printf("file `%s' not found.\n", buf15);
        exit(0);
     }

	 for(i15 = 0; i15 < SNIA_MASS_NUM; i15++)
	 {
    	fscanf(fd15, "%f", &m15);
    	SNIaTotalMetals[i15] = m15 * Chabrier_IMF(SNIaMasses[i15]);
	 }

		//------------------------------------------
		//READ SN-Ia EJECTED MASSES:
		//------------------------------------------
		FILE *fd16;
		char buf16[100];
		int i16;
		float m16;
		static char *name16 = "convolved_stripped_interp_SNIa_EjectedMasses.txt";

		sprintf(buf16, "./YieldTables/%s", name16);

		 if(!(fd16 = fopen(buf16, "r")))
	     {
	        printf("file `%s' not found.\n", buf16);
	        exit(0);
	     }

		 for(i16 = 0; i16 < SNIA_MASS_NUM; i16++)
		 {
	    	fscanf(fd16, "%f", &m16);
	    	SNIaEjectedMasses[i16] = m16 * Chabrier_IMF(SNIaMasses[i16]);
		 }

	//------------------------------------------
	//READ SN-Ia YIELD TABLES:
	//------------------------------------------
		FILE *fd14;
		char buf14[100];
		int i14,j14;
		float m14;
		static char *name14 = "convolved_stripped_interp_SNIa_Yields.txt";

		sprintf(buf14, "./YieldTables/%s", name14);

		 if(!(fd14 = fopen(buf14, "r")))
	     {
	        printf("file `%s' not found.\n", buf14);
	        exit(0);
	     }

		 for(i14=0; i14<42; i14++) //Number of element species (inc H and He) = 42
	     {
	         for(j14=0; j14<SNIA_MASS_NUM; j14++) //Number of initial masses = 48 (30 in initial table)
		     {
		      fscanf(fd14, "%f", &m14);
		      SNIaYields[i14][j14] = m14 * Chabrier_IMF(SNIaMasses[j14]);
		    }
		  }

		printf("SN-Ia yield tables read.\n");
#else
		//------------------------------------------
		//READ SN-Ia YIELD TABLES:
		//------------------------------------------
			FILE *fd14;
			char buf14[100];
			int i14,j14;
			float m14;
			static char *name14 = "SNIaYields.txt";

			sprintf(buf14, "./YieldTables/%s", name14);

			 if(!(fd14 = fopen(buf14, "r")))
		     {
		        printf("file `%s' not found.\n", buf14);
		        exit(0);
		     }

			 for(i14=0; i14<42; i14++) //Number of element species (inc H and He) = 42
		     {
			      fscanf(fd14, "%f", &m14);
			      SNIaYields[i14] = m14;
			      //printf("%f, \n", SNIaYields[i14]);
			  }

			printf("SN-Ia yield tables read.\n");

#endif
}

double Chabrier_IMF(double M)
{
	double e,phi;

	//Coefficient values, normalising the mass-weighted IMF as a function of M for, giving a normalised mass fraction.
	//FOR x = 2.3:
	double A, B;
	//For an IMF normalised over 0.1 --> 120.0 Msun:
	if (SNII_MAX_MASS == 120.0)
	{
		A = 0.842984;
		B = 0.235480;
	}
	//For an IMF normalised over 0.1 --> 100.0 Msun:
	else if (SNII_MAX_MASS == 100.0)
	{
		A = 0.852023;
		B = 0.238004;
	}
	//For an IMF normalised over 0.1 --> 70.0 Msun:
	else if (SNII_MAX_MASS == 70.0)
	{
		A = 0.871761;
		B = 0.243518;
	}
	//For an IMF normalised over 0.1 --> 60.0 Msun:
	else if (SNII_MAX_MASS == 60.0)
	{
		A = 0.881259;
		B = 0.246171;
	}
	//For an IMF normalised over 0.1 --> 50.0 Msun:
	else if (SNII_MAX_MASS == 50.0)
	{
		A = 0.893355;
		B = 0.249550;
	}
	//For an IMF normalised over 0.1 --> 40.0 Msun:
	else if (SNII_MAX_MASS == 40.0)
	{
		A = 0.909581;
		B = 0.254083;
	}
	//For an IMF normalised over 0.1 --> 30.0 Msun:
	else if (SNII_MAX_MASS == 30.0)
	{
		A = 0.933161;
		B = 0.260669;
	}
	else {printf("Chabrier_IMF(): Normalization constants for IMF not known. Check upper mass limit (SNII_MAX_MASS)."); exit(1);}

	//FOR x = 2.0:
	//For an IMF normalised over 0.1 --> 120.0 Msun:
	/*const double A = 0.551390;
	const double B = 0.154026;
	A = 0.551390;
	B = 0.154026;*/

	const double x = 2.3; //Normal Chabrier IMF x = 2.3. Top-heavy IMF e.g. x = 2.0
	const double mc = 0.079;
	const double sigma = 0.69;

	if(M >= 1.0)
	{
		phi = B*M*pow(M,-x);
	}
	else
	{
		e = (1./M)*exp(-pow(log10(M)-log10(mc),2.)/(2.*pow(sigma,2.)));
		phi = A*M*e;
	}

	return phi/M;
}
