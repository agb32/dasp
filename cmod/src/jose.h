#ifndef JOSESUBS_HEADER
#define JOSESUBS_HEADER
/*short int       *header;
 * removed, NAB 02/Apr/2013
 * reason: unused, conflicts upon multiple include */


/* functions */

int	**alloc2d_int(int, int);
float	**alloc2d_float(int, int);
int free2d_float(float **a,int m,int n);
FILE	*openfile(char *, char *);
FILE    *open_image(char *d);
#endif
