#include <string.h>

/*------------------------------------------------------------------------------
  Does string match (the initial characters of) any element of dictionary?

  Return value: index of matching element, or
		-1 if no match, or
		-2 if ambiguous.
*/
int strdict(str, dict)
char *str, *dict[];
{
    int i, match;

    /* initialize to no match */
    match = -1;
    /* go through each string in dictionary */
    for (i = 0; dict[i]; i++) {
	if (strncmp(str, dict[i], strlen(str)) == 0) {
	    /* exact match */
	    if (strcmp(str, dict[i]) == 0) {
		match = i;
		break;
	    /* match */
	    } else if (match == -1) {
		match = i;
	    /* ambiguous */
	    } else {
		match = -2;
		break;
	    }
	    
	}
    }
    return (match);
}
