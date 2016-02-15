typedef IDL_LONG CH_CODE;
#define CH_ERROR 0
#define CH_OK 1
#define CH_OUTOFRANGE 2

/* set up the chunks based on where the objects are; this code
 * is allowed to redefine the ra array at its convenience, to avoid
 * stuff landing on the 0./360. boundary; it defines the arrays 
 * decBounds and raBounds which define the chunks; it leaves an 
 * extra 1 cell width boundary, as well; use unsetchunks
 * to clean up the memory */
CH_CODE 
setchunks(double ra[], double dec[], IDL_LONG nPoints, double minSize,      
					double ***raBounds, double **decBounds, IDL_LONG **nRa,
					IDL_LONG *nDec, double *raOffset);
/* cleans up the memory allocated in setchunks */
CH_CODE
unsetchunks(double ***raBounds, double **decBounds, IDL_LONG **nRa,
						IDL_LONG *nDec);
/* take the objects and the chunks (already defined using setchunks)
 * and assign the objects to the appropriate chunks, with some leeway 
 * given by the parameter marginSize; use unassignchunks to get rid of 
 * the memory allocated here */
CH_CODE
assignchunks(double ra[], double dec[], IDL_LONG nPoints, double raOffset,  
						 double marginSize, double minSize, IDL_LONG ***nChunk,     
						 IDL_LONG ****chunkList, double **raBounds, double *decBounds,  
						 IDL_LONG *nRa, IDL_LONG nDec);
/* clean up memory allocated in assign_chunks */
CH_CODE 
unassignchunks(IDL_LONG ***nChunk, IDL_LONG ****chunkList, IDL_LONG *nRa, 
							 IDL_LONG nDec);
/* utility to find the set of chunks which a given point belongs to;
 * if raBounds wraps around 0/360, it allows -1 and nRa[decChunk] to 
 * be used as raChunkMin or raChunkMax */
CH_CODE
getchunkbounds(double ra, double dec, double marginSize,
							 IDL_LONG **raChunkMin, IDL_LONG **raChunkMax,
							 IDL_LONG *decChunkMin, IDL_LONG *decChunkMax,
							 double **raBounds, double *decBounds, IDL_LONG *nRa,
							 IDL_LONG nDec);
/* utility to find the chunk which a given point belongs to; the same as
 * getchunkbounds with marginSize==0., but a little faster */
CH_CODE
getchunk(double ra, double dec, IDL_LONG *raChunk, IDL_LONG *decChunk, 
				 double **raBounds, double *decBounds, IDL_LONG *nRa, IDL_LONG nDec);
/* gets minimum rarange for a sample, trying 6 different choices */
IDL_LONG 
rarange(double ra[], IDL_LONG nPoints, double minSize, double *raRangeMin,
				double *raOffset);
/* using raoffset determined by rarange, finds minimum and maximum ra's */
IDL_LONG
getraminmax(double ra[], double raOffset, IDL_LONG nPoints, double *raMin, 
						double *raMax);
