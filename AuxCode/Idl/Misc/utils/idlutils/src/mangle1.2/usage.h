/*------------------------------------------------------------------------------
  This is in-lined to ensure uniformity between all programs.
*/
    if (strchr(optstr, 'd')) printf("  -d\t\tadvise defaults and exit\n");

    if (strchr(optstr, 'q')) printf("  -q\t\texecute quietly\n");

    if (strchr(optstr, 'w')) printf("  -w<Wlmfile>\tname of file containing spherical harmonics\n");

    if (strchr(optstr, 'z')) printf("  -z<survey>\tname of survey, or of file containing list of weights\n");

    if (strchr(optstr, 'l')) printf("  -l<lmax>\tmaximum harmonic number\n");

    if (strchr(optstr, 'g')) printf("  -g<lsmooth>\tgaussian smoothing harmonic number (0. = no smooth)\n");

    if (strchr(optstr, 'a')) printf("  -a<angle>[u]\tangle within which to snap cap axes together\n");

    if (strchr(optstr, 'b')) printf("  -b<angle>[u]\tangle within which to snap cap latitudes together\n");

    if (strchr(optstr, 't')) printf("  -t<angle>[u]\tangle within which to snap edge to circle\n");

    if (strchr(optstr, 'y')) printf("  -y<ratio>\tsnap edge to circle only if closer than <ratio> x edge length\n");

    /* if (strchr(optstr, 'm')) printf("  -m<angle>[u]\tangle within which multiple intersections are coincident\n"); */

    if (strchr(optstr, 'c')) printf("  -c<seed>\tseed random number generator with integer <seed>\n");

    if (strchr(optstr, 'r')) printf("  -r<n>\t\tgenerate <n> random points\n");

    if (strchr(optstr, 'j')) printf("  -j[min][,max]\tkeep only polygons with weight in [min, max]\n");

    if (strchr(optstr, 'k')) printf("  -k[min][,max]\tkeep only polygons with area in [min, max] str\n");

    if (strchr(optstr, 'n')) printf("  -n\t\tintersect polygons of infile1 with those of same id in infile2\n");

    if (strchr(optstr, 's')) printf("  -s<n>\t\tskip first <n> characters of each line of infiles\n");

    if (strchr(optstr, 'e')) printf("  -e<n>\t\tread only to <n>'th character of each line (0 = no limit)\n");

    if (strchr(optstr, 'v')) printf("  -vo|-vn\tassign old (o) or new (n) polygon id numbers to output polygons\n");

    if (strchr(optstr, 'f')) {
	printf("  -f\t\tlist frames\n");
	printf("  -f<in>[,<ou>]\tinput, output angular frames\n");
	printf("  -f<azn>,<eln>,<azp>[u]\n\t\t<azn>,<eln> = azimuth, elevation of new pole wrt old frame\n\t\t<azp>,<eln> = azimuth, elevation of old pole wrt new frame\n");
    }

    if (strchr(optstr, 'u')) printf("  -u<in>[,<ou>]\tr radians, d degrees, m arcmin, s arcsec, h hms(RA) & dms(Dec)\n");

    if (strchr(optstr, 'p')) printf("  -p<precision>\tnumber of digits after the decimal place in output angles\n");

    if (strchr(optstr, 'i')) printf("  -i<f>[<n>][u]\tread infile in format <f>, with <n> objects per line\n");

    if (strchr(optstr, 'o')) printf("  -o<f>[u]\twrite outfile in format <f>\n");

    if (strchr(optstr, 'i') || strchr(optstr, 'o'))
	printf("  format <f>:\tc circle, e<i> edges, p polygon, r rectangle, R Region, v vertices\n");

    if (strchr(optstr, 'a') || strchr(optstr, 'b') || strchr(optstr, 't') || strchr(optstr, 'i') || strchr(optstr, 'o'))
	printf("  unit u:\tr radians, d degrees, m arcmin, s arcsec, h hms(RA) & dms(Dec)\n");

    printf("  - for infile\tmeans read from stdin\n");

    printf("  - for outfile\tmeans write to stdout\n");

    printf("mangle documentation is at http://casa.colorado.edu/~ajsh/mangle/\n");
