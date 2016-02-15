IDL_LONG 
friendsoffriends(double x[], double y[], double z[], IDL_LONG nPoints,
								 double linkSep, IDL_LONG **nChunk, IDL_LONG ***chunkList,
								 IDL_LONG *nRa, IDL_LONG nDec, IDL_LONG *firstGroup,
								 IDL_LONG *multGroup, IDL_LONG *nextGroup,
								 IDL_LONG *inGroup, IDL_LONG *nGroups);
IDL_LONG 
chunkfriendsoffriends(double x[], double y[], double z[], IDL_LONG chunkList[],
											IDL_LONG nTargets, double linkSep, IDL_LONG firstGroup[],
											IDL_LONG multGroup[], IDL_LONG nextGroup[],
											IDL_LONG inGroup[], IDL_LONG *nGroups);
