
short int       *header;


/* functions */

int	**alloc2d_int(int, int);
float	**alloc2d_float(int, int);
int free2d_float(float **a,int m,int n);
FILE	*openfile(char *, char *);
FILE    *open_image(char *d);
