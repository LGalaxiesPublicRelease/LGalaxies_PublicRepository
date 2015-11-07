void generate_timebins()
{
  double previoustime, newtime, deltaT;
  int snap,step,sfh_ibin, i,j,sfh_dt[SFH_NBIN],sfh_t[SFH_NBIN];
  double timet,age_in_years;

  for (i=0;i<SFH_NBIN;i++)
  {
     sfh_dt[i]=0;
     sfh_t[i]=0;
   }
   sfh_ibin=0;
   sfh_dt[0]=1;
   sfh_t[0]=1;

  for(snap=0;snap<64;snap++)
  {
	  previoustime = NumToTime(snap);
	  newtime = NumToTime(snap+1);
	  deltaT = previoustime - newtime;


	  for(step=0;step<20;step++)
	  {
		  timet = previoustime - (step + 0.5) * (deltaT / STEPS);
		  age_in_years=(Age[0]-timet)*UnitTime_in_years/Hubble_h;

		  int ibin;
		  float t; // time in units of SFH_TIME_INTERVAL
		  int flag_merged_bins; // Boolean used to check whether have merged bins
		  int dt_merge; // Size of bins that we are checking for merging
		  int n_merge; // Number of bins of this size

		  t=age_in_years/SFH_TIME_INTERVAL;
		  ibin=sfh_ibin;
		  printf("Snap=%d Step=%d \n",snap,step);

		  while(sfh_t[ibin]<t)
		    {
		    // Add new bins, either to the correct time, or to the bin limit
		      for (i=ibin;sfh_t[i]<t && i<SFH_NBIN-1;i++)
		        {
		    	  if (sfh_dt[i+1]==0)
		    	    {
		    		  sfh_dt[i+1]=1;
		    		  sfh_t[i+1]=sfh_t[i]+1;
		    	    }
		        }
		    ibin=i;


		    // Now merge bins where we have three (or more) of the same size.
		    // Need to do this iteratively.
		    flag_merged_bins=1;
		    while(flag_merged_bins)
		      {
		    	flag_merged_bins=0;
		    	dt_merge=sfh_dt[0];
		    	i=0;
		    	// Will have checked all bins once dt_merge drops to zero
		    	while(!flag_merged_bins && dt_merge>0)
		    	  {
		    		// Count number of bins of this size
		    		n_merge=0;
		    		// The i=i below is to suppress a warning message
		    		for(i=i;sfh_dt[i]==dt_merge;i++)
		    			n_merge+=1;
		    		// If fewer than 3 bins then do nothing
		    		// (4 bins if dt_merge=1)
		    		// else exit loop and flag for merging
		    		if (n_merge<3 || (n_merge==3 && dt_merge==1))
		    		  {
		    			dt_merge/=2;
		    			n_merge=0;
		    		  }
		    		else
		    		  {
		    			flag_merged_bins=1;
		    			i=i-n_merge;
		    		  }
		    	  }
		    	// At this point, if flag_merged_bins is set then
		    	// we have to merge 3 bins into 2.
		    	if (flag_merged_bins)
		    	  {
		    		// Merge bins i and i+1
		    		sfh_dt[i]*=2;
		    		sfh_t[i]=sfh_t[i+1];
		    		// Relabel all the other bins
		    		for(i=i+1;i<ibin;i++)
		    		  {
		    			sfh_dt[i]=sfh_dt[i+1];
		    			sfh_t[i]=sfh_t[i+1];
		    		  }
		    		sfh_dt[i]=0;
		    		sfh_t[i]=0;
		    		ibin=i-1;
		    	  }
		      } // End loop over bin merging

			  if (sfh_t[ibin]<t && ibin==SFH_NBIN-1) {
				  printf("sfh_update_bins: too many bins required\n");
				  printf("...did you remember to call sfh_initialise?\n");
				  exit(1);
			  }
		  } // End while loop over time update

		  sfh_ibin=ibin;



		  printf("t=%f Snap=%d Step=%d Previous Time=%e Age=%e\n",
				  t, snap,step, previoustime*UnitTime_in_years/Hubble_h,age_in_years);
		  for(j=0;j<SFH_NBIN;j++)
			  printf("Bin[%d] t=%d t=%d\n",j,sfh_t[j],sfh_dt[j]);

	  }
  }

}
