#ifndef INPUTFILE_H
#define INPUTFILE_H

typedef struct {
    char *name;		/* filename */
    FILE *file;		/* file stream */
    char *line;		/* line buffer */
    size_t bufsize;	/* size of line buffer (will expand as necessary) */
    int line_number;	/* line number */
    int end;		/* maximum number of characters to read (0 = no limit) */
} inputfile;

#endif	/* INPUTFILE_H */
