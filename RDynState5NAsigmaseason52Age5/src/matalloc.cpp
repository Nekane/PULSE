#include <stdarg.h>
#include <stdlib.h>
#include "matalloc.h"
/*
   void *matalloc (size_t size, void *data, int ndim, ...)

   allocates memory and builds pointerstructure for
   ndim dimensional matrix. 

   size: size in bytes of matrix elements to be allocated
   data: if (void *)0: matalloc allocates memory
         else it is assumed data is already allocated
   ndim: number of dimensions of matrix to allocate
   ... : dimensions of the matrix to be allocated of type unsigned int

   return value: pointer to pointerstructure to access arrayelements
                 (void *)0 if something wrong with parameters or
                           if malloc doesn't work

   example:

    allocate 2-dimensional matrix, elements of type double,
    dimensions M and N:

      double **a;
      a=matalloc(sizeof(double),(void *)0,2,M,N);

      usage of a is then like: a[i][j] = 3;

    allocate 4-dimensional matrix, elements of type int,
    dimensions m,n,p,q:

      int ****a;
      m=8; n=123; p=9; q=11;
      a=matalloc(sizeof(int),(void *)0,4,m,n,p,q);
      
      a[i][j][k][l] = 12;

    3-dimensional matrix, elements of type double,
    dimensions l,m,n. Space for data is already allocated,
    so in fact, only the pointer structure is set up:

    double *a;
    double ***x;

    a=(double *) malloc(sizeof(double)*l*m*n);  
    
       a allocated outside matalloc 

    x=matalloc(sizeof(double),a,3,l,m,n);

       x[0][0][0] is the same element as a[0] 

  void matfree(void *p)
    
    Frees memory occupied by pointers and perhaps data allocated by
    matalloc. Only the space allocated by matalloc is freed.

    p:  returnvalue of matalloc

    example

    char ***a;
    a=matalloc(sizeof(char),(void *)0,3,10,20,15);
    
    matfree(a);
*/

void * matalloc_1(unsigned int i);
unsigned int *d;   /* array to hold dimensions */
int nd;   /* number of dimensions */
void **p;         /* points to p-array: array of pointers */
char *q;         /* points to q-array: array of characters */
void **pp;        /* points to first not allocated element of p-array */
unsigned int np;      /* number of elements in p-array */
unsigned int nq;      /* number of bytes in q-array */
unsigned int sz;      /* number of bytes per allocated matrix element */

void *matalloc(size_t size, void *data, int ndim, ...)
  {
    va_list argp;
    int i;
    int npad=64;
    size_t s;
    void *r;
    unsigned int pspace;

    if (ndim <=0)
	    return (void *)0;

    sz=size;
    nd=0;

    d=(unsigned int *)malloc(ndim*sizeof(int));
    if (d == 0)
      return (void *)0;

    va_start(argp, ndim);
    for (i=0; i<ndim; i++)
      {
        d[nd]  = va_arg(argp, unsigned int);
        nd++;
      }
    va_end(argp);

    if (nd >=2)
      {
        np = d[nd-2];
        for (i=nd-3; i>=0; i--)
	      np = d[i]*(1+np);
      }
    else
      np=0;

    nq=1;
    for (i=0; i<nd; i++)
	    nq *= d[i];

/* calculate length of p-array
   round it up to multiple of npad bytes if no data is given*/

    if (data)
      s=np*sizeof(void*);
    else
      s=((np * sizeof(void*) - 1)/npad + 1)*npad;

/* allocate in one step the space for p-array and q-array */
/* don't allocate space for q-array if data is given */

    pspace = s;
    if (!data)
      pspace = s+nq*size;
    p = pp = (void **)malloc(pspace);

   if ( p == 0)
    {
     free(d);
     return (void *)0;
    }

/* start of q-array after p-array if data not given */

    if (data)
      q = (char *) data;
    else
      q = (char *)p + s;

    r= matalloc_1(0);
    free(d);
    return r;
  }

void *matalloc_1(unsigned int ld)
  {
    int i;
    void **pr, **pl;
    char *ql;

    if (ld == nd - 1 )
	    /* simple: at this point we return a pointer to
	     * a space of d[ld]*size in the q-array
	     */
      {
        ql =  q;
	q += sz * d[ld];
	return (void *)ql;
      }
    /* This is more interesting: 
     * get d[ld] pointers from the p-array
     * and assign the values of matalloc-1(ld+1) to them */

    pr = pl = pp;
    pp  += d[ld];

    for (i=0; i < d[ld] ; i++)
       {
	 *pl = matalloc_1(ld+1);
	 pl++;
       }

    return pr;
  }    

void matfree(void *p)
  {
    free(p);
  }
